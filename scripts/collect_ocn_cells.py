"""Collect a per-cell table of ESA L2 OCN (owi) angles for IW / EW / SM products.

For each product the footprint and the ``orbitDirection`` attribute are read
from the CDSE daily Parquets (metadata only), and every ``step``-th ``owi``
cell is written with the ESA angles (incidence, elevation, heading), the
current s1-geometry estimates and the intermediate quantities the models need
(across-track position, footprint and orbit-model satellite distances, pass
direction from the footprint and from the catalogue).  The table feeds the
IW / EW / SM calibration and validation (``calibrate_iw_ew_sm.py``).

Usage::

    python scripts/collect_ocn_cells.py --out cells.parquet --mode IW --sat S1C \\
        --every 10 --per-day 6 --workers 12
"""

from __future__ import annotations

import argparse
import glob
import json
import math
import sys
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime
from functools import lru_cache
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from s1_geometry import (  # noqa: E402
    estimate_incidence_elevation,
    heading_from_lat,
    orbit_direction_from_ring,
    satellite_radius_km,
)
from s1_geometry.angles import _range_lines  # noqa: E402

LEGACY = Path("/legacy/project/cersat/public/cache/project/mpc-sentinel1/data/esa")
PARQUET = Path("/scale/project/lops-siam-sentinel1-workbench/data/ifremer/"
               "metadata_cdse_sentinel1/daily_parquet_files_v2")


@lru_cache(maxsize=8)
def _day_table(day: str) -> dict:
    import pyarrow.parquet as pq
    path = PARQUET / f"s1_{day[:4]}-{day[4:6]}-{day[6:8]}.parquet"
    t = pq.read_table(path, columns=["Name", "GeoFootprint", "Attributes"]).to_pylist()
    return {r["Name"]: r for r in t if "_OCN_" in r["Name"]}


def _near_far(ring, lines, direction: str):
    """(near, far) range-line midpoints for a given pass direction."""
    pts = np.asarray(json.loads(ring)["coordinates"][0] if isinstance(ring, str) else ring)[:4]
    lat0 = float(pts[:, 1].mean())
    heading = math.radians(heading_from_lat(lat0, direction))
    look = (math.cos(heading), -math.sin(heading))
    a, b = lines["near_mid"], lines["far_mid"]
    pa = a[0] * look[0] + a[1] * look[1]
    pb = b[0] * look[0] + b[1] * look[1]
    return (a, b) if pa <= pb else (b, a)


def collect(args) -> list[dict] | None:
    safe, step = args
    import netCDF4
    name = Path(safe).name
    parts = name.split("_")
    mode = parts[1]
    beam = mode if mode.startswith("S") and mode[1:].isdigit() else None
    lib_mode = "SM" if beam else mode
    day = parts[5][:8]
    try:
        row = _day_table(day)[name]
        attrs = {a["Name"]: a["Value"] for a in json.loads(row["Attributes"])}
        attr_dir = str(attrs.get("orbitDirection", "")).lower()
        ring, source = row["GeoFootprint"], "cdse"
    except (KeyError, OSError, FileNotFoundError):
        # not in the CDSE catalogue (e.g. every SM OCN product): the manifest
        # footprint, reversed to the CDSE vertex order, and its <s1:pass>
        try:
            text = Path(safe, "manifest.safe").read_text()
        except OSError:
            return None
        import re
        coords = re.findall(r"<gml:coordinates>([^<]+)</gml:coordinates>", text)
        pas = re.search(r"<s1:pass>(\w+)</s1:pass>", text)
        if not coords or not pas:
            return None
        pts = [[float(v) for v in q.split(",")][::-1] for q in coords[0].split()][::-1]
        ring, source = pts + [pts[0]], "manifest"
        attr_dir = pas.group(1).lower()
    try:
        lines = _range_lines(ring)
        ring_dir = orbit_direction_from_ring(ring)
    except ValueError:
        return None
    files = sorted(Path(safe, "measurement").glob("*-ocn-*.nc"))
    if not files:
        return None
    with netCDF4.Dataset(files[0]) as ds:
        get = {}
        for key, var in (("lon", "owiLon"), ("lat", "owiLat"), ("inc", "owiIncidenceAngle"),
                         ("el", "owiElevationAngle"), ("hd", "owiHeading")):
            if var in ds.variables:
                get[key] = np.ma.filled(ds[var][:].astype(float), np.nan)[::step, ::step].ravel()
        if "lon" not in get or "inc" not in get:
            return None
        n = get["lon"].size
        for key in ("el", "hd"):
            get.setdefault(key, np.full(n, np.nan))
    ok = np.isfinite(get["lon"]) & np.isfinite(get["lat"]) & np.isfinite(get["inc"])
    near_r, far_r = lines["near_mid"], lines["far_mid"]
    near_a, far_a = _near_far(ring, lines, attr_dir) if attr_dir in ("ascending", "descending") \
        else (near_r, far_r)
    out = []
    for i in np.flatnonzero(ok):
        lon, lat = float(get["lon"][i]), float(get["lat"][i])
        try:
            a = estimate_incidence_elevation(ring, lon, lat, lib_mode, beam)
            inc_mod, el_mod = a.incidence_deg, a.elevation_deg
        except ValueError:
            inc_mod = el_mod = np.nan

        def frac(nr, fr, lon=lon, lat=lat):
            ax = (fr[0] - nr[0], fr[1] - nr[1])
            return ((lon - nr[0]) * ax[0] + (lat - nr[1]) * ax[1]) / (ax[0] ** 2 + ax[1] ** 2)

        def frac_metric(nr, fr, lon=lon, lat=lat):
            # local tangent plane at the swath centre: longitude scaled by cos(lat)
            k = math.cos(math.radians(0.5 * (nr[1] + fr[1])))
            ax = ((fr[0] - nr[0]) * k, fr[1] - nr[1])
            v = ((lon - nr[0]) * k, lat - nr[1])
            return (v[0] * ax[0] + v[1] * ax[1]) / (ax[0] ** 2 + ax[1] ** 2)

        out.append({
            "product": name, "sat": name[:3], "mode": lib_mode, "beam": beam or "",
            "source": source,
            "day": day, "attr_dir": attr_dir, "ring_dir": ring_dir, "lon": lon, "lat": lat,
            "dfrac_ring": frac(near_r, far_r), "dfrac_attr": frac(near_a, far_a),
            "dfrac_attr_m": frac_metric(near_a, far_a),
            "near_lon": near_a[0], "near_lat": near_a[1], "far_lon": far_a[0], "far_lat": far_a[1],
            "width_km": lines["width_km"], "r_s_fp": lines["r_s_km"] or np.nan,
            "r_s_orb": float(satellite_radius_km(lat, attr_dir))
            if attr_dir in ("ascending", "descending") else np.nan,
            "inc_ocn": float(get["inc"][i]), "el_ocn": float(get["el"][i]),
            "hd_ocn": float(get["hd"][i]), "inc_mod": inc_mod, "el_mod": el_mod,
            "hd_mod_ring": heading_from_lat(lat, ring_dir),
            "hd_mod_attr": heading_from_lat(lat, attr_dir)
            if attr_dir in ("ascending", "descending") else np.nan,
        })
    return out


def _safes(mode: str, sat: str, every: int, per_day: int, years) -> list[str]:
    s = sat[-1].lower()
    if mode == "SM":
        root = LEGACY / f"sentinel-1{s}/L2/SM/{sat}_SM_OCN__2S"
    else:
        root = LEGACY / f"sentinel-1{s}/L2/{mode}/{sat}_{mode}_OCN__2S"
    days = []
    for y in sorted(p.name for p in root.iterdir()):
        if years and int(y) not in years:
            continue
        days += sorted(glob.glob(str(root / y / "*")))
    out = []
    for d in days[::every]:
        safes = sorted(glob.glob(d + "/*.SAFE"))
        if safes:
            idx = np.unique(np.linspace(0, len(safes) - 1, min(per_day, len(safes))).round().astype(int))
            out += [safes[k] for k in idx]
    return out


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--out", required=True)
    parser.add_argument("--mode", choices=("IW", "EW", "SM"), required=True)
    parser.add_argument("--sat", required=True, help="S1A, S1C or S1D")
    parser.add_argument("--every", type=int, default=10, help="take one day every N")
    parser.add_argument("--per-day", type=int, default=6)
    parser.add_argument("--step", type=int, default=6, help="owi grid decimation")
    parser.add_argument("--years", type=int, nargs="*", default=None)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args(argv)
    import pandas as pd
    safes = _safes(args.mode, args.sat, args.every, args.per_day, args.years)
    rows, n_ok = [], 0
    with ProcessPoolExecutor(args.workers) as pool:
        for res in pool.map(collect, [(s, args.step) for s in safes], chunksize=2):
            if res:
                rows += res
                n_ok += 1
    pd.DataFrame(rows).to_parquet(args.out)
    print(f"{args.sat} {args.mode}: {n_ok}/{len(safes)} products, {len(rows)} cells -> {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
