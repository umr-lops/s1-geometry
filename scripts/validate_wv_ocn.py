"""Validate the metadata-only WV angle chain against ESA L2 WV OCN products.

For every ``S1?_WV_OCN__2S*.SAFE`` given, the product footprint is read from
``manifest.safe`` (the same imagette rings, in the same time order, as the
CDSE OData ``Footprint``), the angles are estimated with
:func:`s1_geometry.estimate_wv_angles` at every ``owi`` cell of every
imagette -- ``r_s`` from the orbit model, i.e. metadata only -- and compared
with ``owiIncidenceAngle``, ``owiElevationAngle`` and ``owiHeading``.

Usage::

    python scripts/validate_wv_ocn.py --out stats.csv /path/to/*.SAFE
    python scripts/validate_wv_ocn.py --out stats.csv --sample 3 /path/S1C_WV_OCN__2S/2026/*/

(a directory that is not a SAFE is searched for SAFE products; ``--sample k``
keeps ``k`` products per directory, spread over the day).  Needs ``netCDF4``.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from datetime import datetime
from pathlib import Path

import numpy as np

from s1_geometry import estimate_wv_angles


def manifest_footprint(safe: Path) -> list[np.ndarray]:
    """Imagette rings (lon, lat) from ``manifest.safe`` (file order = time order)."""
    text = (safe / "manifest.safe").read_text()
    rings = []
    for block in re.findall(r"<gml:coordinates>([^<]+)</gml:coordinates>", text):
        latlon = np.array([p.split(",") for p in block.split()], dtype=float)
        ring = latlon[:, ::-1]
        rings.append(np.vstack([ring, ring[:1]]))
    return rings


def read_owi(safe: Path) -> dict[str, np.ndarray]:
    """Concatenated owi cells of all the imagettes of a product."""
    import netCDF4

    out: dict[str, list] = {k: [] for k in ("lon", "lat", "inc", "el", "hd", "sw")}
    for path in sorted((safe / "measurement").glob("*.nc")):
        sw = 1 if "-wv1-" in path.name else 2
        with netCDF4.Dataset(path) as ds:
            ds.set_auto_mask(True)
            arrays = [np.ma.filled(ds[v][:].astype(float), np.nan).ravel() for v in (
                "owiLon", "owiLat", "owiIncidenceAngle", "owiElevationAngle", "owiHeading")]
        for key, arr in zip(("lon", "lat", "inc", "el", "hd"), arrays, strict=True):
            out[key].append(arr)
        out["sw"].append(np.full(arrays[0].size, sw))
    data = {k: np.concatenate(v) for k, v in out.items()}
    ok = np.isfinite(data["lon"]) & np.isfinite(data["lat"]) & np.isfinite(data["inc"])
    return {k: v[ok] for k, v in data.items()}


def _stats(err: np.ndarray) -> tuple[float, float, float]:
    err = err[np.isfinite(err)]
    if err.size == 0:
        return (np.nan, np.nan, np.nan)
    return (float(np.sqrt(np.mean(err ** 2))), float(err.mean()), float(np.abs(err).max()))


def validate(safe: Path, window) -> dict:
    rings = manifest_footprint(safe)
    owi = read_owi(safe)
    res = estimate_wv_angles(rings, None, owi["lon"], owi["lat"], window=window,
                             on_outside="nan")
    name = safe.name
    t0, t1 = (datetime.strptime(s, "%Y%m%dT%H%M%S") for s in name.split("_")[5:7])
    row = {
        "product": name,
        "satellite": name[:3],
        "pass": re.search(r"<s1:pass>(\w+)", (safe / "manifest.safe").read_text()).group(1),
        "duration_min": round((t1 - t0).total_seconds() / 60.0, 2),
        "n_rings": len(rings),
        "window": str(window),
        "lat_min": float(owi["lat"].min()),
        "lat_max": float(owi["lat"].max()),
        "r_s_km": res.r_s_km,
        "n": int(np.isfinite(res.incidence_deg).sum()),
        "n_outside": int(np.isnan(res.incidence_deg).sum()),
        "subswath_agree": float(np.mean(res.subswath == owi["sw"])),
    }
    errors = {
        "inc": res.incidence_deg - owi["inc"],
        "el": res.elevation_deg - owi["el"],
        "hd": (res.heading_deg - owi["hd"] + 180.0) % 360.0 - 180.0,
    }
    for key, err in errors.items():
        row[f"{key}_rmse"], row[f"{key}_bias"], row[f"{key}_max"] = _stats(err)
        row[f"{key}_p95"] = float(np.nanpercentile(np.abs(err), 95))
        for sw in (1, 2):
            row[f"{key}_bias_wv{sw}"] = _stats(err[owi["sw"] == sw])[1]
    return row


def _expand(paths: list[str], sample: int | None) -> list[Path]:
    out = []
    for p in map(Path, paths):
        if p.name.endswith(".SAFE"):
            out.append(p)
            continue
        safes = sorted(p.glob("*WV_OCN__2S*.SAFE"))
        if sample and len(safes) > sample:
            safes = [safes[i] for i in np.linspace(0, len(safes) - 1, sample).round().astype(int)]
        out.extend(safes)
    return out


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("paths", nargs="+", help="SAFE products or directories of SAFEs")
    parser.add_argument("--out", required=True, help="output CSV (one row per product)")
    parser.add_argument("--sample", type=int, default=None, help="products per directory")
    parser.add_argument("--window", default="auto",
                        help="'auto' (default, as add_wv_angle_columns), 'none' or an integer")
    args = parser.parse_args(argv)
    window = {"auto": "auto", "none": None}.get(args.window.lower())
    if window is None and args.window.lower() != "none":
        window = int(args.window)

    rows = []
    for safe in _expand(args.paths, args.sample):
        try:
            rows.append(validate(safe, window))
        except (ValueError, OSError, KeyError) as error:
            print(f"SKIP {safe.name}: {error}", file=sys.stderr)
            continue
        r = rows[-1]
        print(f"{r['product'][:62]} {r['pass'][:4]} {r['n_rings']:3d} rings "
              f"lat {r['lat_min']:6.1f}..{r['lat_max']:6.1f}  inc {r['inc_rmse']:.3f} "
              f"el {r['el_rmse']:.3f} hd {r['hd_rmse']:.3f}  sw {r['subswath_agree']:.3f}")
    if rows:
        with open(args.out, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    return 0


if __name__ == "__main__":
    sys.exit(main())
