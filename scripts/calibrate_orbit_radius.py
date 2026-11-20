"""Calibrate the satellite geocentric distance model ``r_s(lat, pass)``.

Fits the coefficients embedded in :mod:`s1_geometry.orbit` from Sentinel-1
precise orbit files (``AUX_POEORB``, ~4.6 MB per satellite-day).  They are
public on the ASF ``s1-orbits`` S3 bucket (no authentication)::

    python scripts/calibrate_orbit_radius.py --download S1A_OPER_AUX_POEORB_... ...
    python scripts/calibrate_orbit_radius.py /path/to/*.EOF

For every state vector the target point of a right-looking beam at the WV
mean off-nadir angle is computed, and ``r_s`` is fitted against the target
latitude (the latitude available in the metadata: the reference point), one
polynomial per pass direction (``VZ > 0`` ascending).  The script prints the
coefficients (to paste into ``orbit.py``) and the residual statistics per
file, and can write a subsampled CSV fixture for the unit tests.
"""

from __future__ import annotations

import argparse
import re
import sys
import urllib.request
from pathlib import Path

import numpy as np

S3_URL = "https://s1-orbits.s3.amazonaws.com/AUX_POEORB/{name}"
R_SPHERE_KM = 6371.0
#: Off-nadir angle (deg) of the fit target (midway between WV1 and WV2).
FIT_OFF_NADIR_DEG = 26.0
DEGREE = 10


def read_eof(path: Path) -> dict[str, np.ndarray]:
    """State vectors (ECEF, km and km/s) of a POEORB EOF file."""
    text = path.read_text()
    out = {}
    for key in ("X", "Y", "Z", "VX", "VY", "VZ"):
        unit = "m/s" if key.startswith("V") else "m"
        values = re.findall(rf'<{key} unit="{unit}">([-+\d.eE]+)</{key}>', text)
        out[key] = np.asarray(values, dtype=float) / 1e3
    return out


def target_latitude(sv: dict[str, np.ndarray], off_nadir_deg: float) -> np.ndarray:
    """Geocentric latitude (deg) of the right-looking beam target (sphere)."""
    pos = np.stack((sv["X"], sv["Y"], sv["Z"]), axis=1)
    vel = np.stack((sv["VX"], sv["VY"], sv["VZ"]), axis=1)
    r_s = np.linalg.norm(pos, axis=1)
    up = pos / r_s[:, None]
    right = np.cross(vel, up)
    right /= np.linalg.norm(right, axis=1)[:, None]
    sig = np.radians(off_nadir_deg)
    # central angle between nadir and the target on the sphere
    inc = np.arcsin(np.clip(r_s / R_SPHERE_KM * np.sin(sig), -1.0, 1.0))
    gamma = inc - sig
    target = np.cos(gamma)[:, None] * up + np.sin(gamma)[:, None] * right
    return np.degrees(np.arcsin(np.clip(target[:, 2], -1.0, 1.0)))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("files", nargs="+", help="EOF files (or names with --download)")
    parser.add_argument("--download", action="store_true",
                        help="download the named POEORB files from the s1-orbits bucket")
    parser.add_argument("--outdir", default=".", help="download directory")
    parser.add_argument("--fixture", help="write a subsampled CSV fixture here")
    args = parser.parse_args(argv)

    paths = []
    for name in args.files:
        path = Path(args.outdir) / Path(name).name
        if args.download and not path.exists():
            urllib.request.urlretrieve(S3_URL.format(name=Path(name).name), path)
        paths.append(path if args.download else Path(name))

    data = []
    for path in paths:
        sv = read_eof(path)
        r_s = np.sqrt(sv["X"] ** 2 + sv["Y"] ** 2 + sv["Z"] ** 2)
        data.append((path.name[:3], path.name, r_s,
                     target_latitude(sv, FIT_OFF_NADIR_DEG),
                     target_latitude(sv, 20.9), target_latitude(sv, 31.5),
                     sv["VZ"] > 0))

    for label, is_asc in (("ascending", True), ("descending", False)):
        lat = np.concatenate([d[3][d[6] == is_asc] for d in data])
        r_s = np.concatenate([d[2][d[6] == is_asc] for d in data])
        fit = np.polynomial.Polynomial.fit(lat / 90.0, r_s, DEGREE, domain=[-1, 1])
        coef = fit.convert(domain=[-1, 1]).coef
        print(f"{label}: lat range {lat.min():.2f} .. {lat.max():.2f}")
        print(f"    coefficients (x = lat/90, low to high): {tuple(float(c) for c in coef)}")
        for _sat, name, rs, _, lat1, lat2, asc in data:
            m = asc == is_asc
            for beam, lt in (("WV1", lat1), ("WV2", lat2)):
                x = np.clip(lt[m], lat.min(), lat.max()) / 90.0
                res = rs[m] - np.polynomial.polynomial.polyval(x, coef)
                print(f"    {name[:3]} {name[42:50]} {beam}: bias {res.mean():+.3f} "
                      f"std {res.std():.3f} max {np.abs(res).max():.3f} km")

    if args.fixture:
        rows = []
        for sat, name, rs, _, lat1, lat2, asc in data:
            sel = slice(0, None, 97)
            for a, b, c, d in zip(lat1[sel], lat2[sel], rs[sel], asc[sel], strict=True):
                rows.append(f"{sat},{name[42:50]},{a:.4f},{b:.4f},{c:.4f},{int(d)}")
        Path(args.fixture).write_text(
            "satellite,day,lat_wv1,lat_wv2,r_s_km,ascending\n" + "\n".join(rows) + "\n")
        print(f"fixture: {len(rows)} rows -> {args.fixture}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
