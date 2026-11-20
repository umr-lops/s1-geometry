"""Calibrate the WV look-angle law (``wv_theta_deg``) on ESA L2 WV OCN products.

For every imagette of the given products, the metadata-only chain
(:func:`s1_geometry.estimate_wv_angles`: footprint + orbit-model ``r_s``) is
run with the current law and with the law + 0.1 deg; one Newton step gives
the off-nadir angle ``theta_needed`` that cancels the mean incidence error
against ``owiIncidenceAngle`` (or ``oswIncidenceAngle``) over the imagette.
The law is then refitted as::

    theta(sw, r_s, lat) = base(sw, r_s) + P_sw(lat)

with ``base`` the documented altitude law and ``P_sw`` a polynomial in the
imagette latitude (one per subswath).  Usage::

    python scripts/calibrate_wv_look_angle.py --rings rings.csv --workers 12 \\
        /path/S1C_WV_OCN__2S/2026/0[0-9]0/ ...          # collect per-imagette rows
    python scripts/calibrate_wv_look_angle.py --rings rings.csv --fit 4   # fit only

Needs ``netCDF4``; ``validate_wv_ocn.py`` provides the readers.
"""

from __future__ import annotations

import argparse
import csv
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

from validate_wv_ocn import manifest_footprint, read_owi  # noqa: E402

import s1_geometry.wv as wv  # noqa: E402
from s1_geometry.orbit import satellite_radius_km  # noqa: E402

DTHETA = 0.1
FIELDS = ("product", "ring", "lat", "subswath", "direction", "r_s_km", "theta_used",
          "dinc", "dinc_dtheta", "theta_needed", "n")


def _base(sw: int, r_s: float) -> float:
    """The documented altitude law alone (no latitude term)."""
    c_min, c_max = wv.WV_THETA_DOC[sw]
    slope = (c_max - c_min) / (wv._ALT_MAX - wv._ALT_MIN)
    return 0.5 * (c_min + c_max) + slope * (r_s - wv.R_E_KM - wv._REF_ALT)


def collect(safe: Path) -> list[dict]:
    """Per-imagette Newton rows of one product (empty list if unusable)."""
    try:
        rings = manifest_footprint(safe)
        owi = read_owi(safe)
    except (OSError, KeyError, IndexError, ValueError):
        return []
    if len(rings) < 4 or owi["lon"].size == 0:
        return []
    step = max(1, owi["lon"].size // (len(rings) * 60))  # ~60 cells per imagette
    lon, lat, inc_ref = owi["lon"][::step], owi["lat"][::step], owi["inc"][::step]
    law = wv.wv_theta_deg
    try:
        res0 = wv.estimate_wv_angles(rings, None, lon, lat, on_outside="nan")
        wv.wv_theta_deg = lambda s, r, la: law(s, r, la) + DTHETA
        res1 = wv.estimate_wv_angles(rings, None, lon, lat, on_outside="nan")
    except ValueError:
        return []
    finally:
        wv.wv_theta_deg = law
    closed = [wv._closed(r) for r in wv.merge_antimeridian_splits(wv.as_multipolygon(rings))]
    centers = np.array([wv._to_3d(r).mean(axis=0) for r in closed])
    directions = wv._pass_direction_per_ring(centers)
    rows = []
    for i, ring in enumerate(closed):
        m = (res0.ring_index == i) & np.isfinite(res0.incidence_deg) & np.isfinite(inc_ref)
        if m.sum() < 5:
            continue
        ring_lat = float(ring[:, 1].mean())
        sw = int(res0.subswath[m][0])
        r_s = float(satellite_radius_km(ring_lat, directions[i]))
        dinc = float(np.mean(res0.incidence_deg[m] - inc_ref[m]))
        deriv = float(np.mean(res1.incidence_deg[m] - res0.incidence_deg[m])) / DTHETA
        theta = law(sw, r_s, ring_lat)
        rows.append({
            "product": safe.name, "ring": i, "lat": ring_lat, "subswath": sw,
            "direction": directions[i], "r_s_km": r_s, "theta_used": theta, "dinc": dinc,
            "dinc_dtheta": deriv, "theta_needed": theta - dinc / deriv, "n": int(m.sum()),
        })
    return rows


def fit(path: str, degree: int) -> None:
    data = np.genfromtxt(path, delimiter=",", names=True, dtype=None, encoding="utf-8")
    print(f"{data.size} imagettes, {np.unique(data['product']).size} products, "
          f"lat {data['lat'].min():.1f}..{data['lat'].max():.1f}")
    for sw in (1, 2):
        d = data[data["subswath"] == sw]
        base = np.array([_base(sw, r) for r in d["r_s_km"]])
        y = d["theta_needed"] - base
        coef = np.polynomial.polynomial.polyfit(d["lat"] / 90.0, y, degree)
        res = y - np.polynomial.polynomial.polyval(d["lat"] / 90.0, coef)
        coefs = tuple(round(float(c), 6) for c in coef)
        print(f"WV{sw}: n={d.size} P(x=lat/90) coefficients low->high {coefs}")
        print(f"     residual std {res.std():.3f} deg, p95 {np.percentile(np.abs(res), 95):.3f}")
        for direction in ("ascending", "descending"):
            m = d["direction"] == direction
            print(f"     {direction}: mean residual {res[m].mean():+.3f} (n={m.sum()})")
        bins = np.arange(-80, 81, 10)
        idx = np.digitize(d["lat"], bins)
        line = " ".join(f"{bins[k - 1]:+d}:{res[idx == k].mean():+.2f}"
                        for k in range(1, len(bins)) if (idx == k).sum() > 20)
        print(f"     residual by 10-deg lat bin: {line}")


def _expand(paths: list[str]) -> list[Path]:
    out = []
    for p in map(Path, paths):
        out.extend([p] if p.name.endswith(".SAFE") else sorted(p.glob("*WV_OCN__2S*.SAFE")))
    return out


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("paths", nargs="*", help="SAFE products or directories of SAFEs")
    parser.add_argument("--rings", required=True, help="per-imagette CSV (written / read)")
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--fit", type=int, default=None, help="polynomial degree to fit")
    args = parser.parse_args(argv)
    if args.paths:
        safes = _expand(args.paths)
        with open(args.rings, "w", newline="") as f, ProcessPoolExecutor(args.workers) as pool:
            writer = csv.DictWriter(f, fieldnames=FIELDS)
            writer.writeheader()
            for k, rows in enumerate(pool.map(collect, safes, chunksize=2)):
                writer.writerows(rows)
                if k % 50 == 0:
                    print(f"{k}/{len(safes)} products", flush=True)
    if args.fit is not None:
        fit(args.rings, args.fit)
    return 0


if __name__ == "__main__":
    sys.exit(main())
