"""End-to-end: SAR incidence / elevation angles from a Parquet of footprints.

Usage:

    python examples/incidence_elevation_end2end.py [input.parquet] [output.parquet]

Without arguments, the bundled sample (``examples/data/
s1_sample_footprints.parquet``: 5 first IW + 5 first EW footprints of the
CDSE daily Parquet s1_2026-11-13, with random in-footprint reference
positions) is used, and the result is written next to it.

The input Parquet has one row per acquisition with the columns ``Name``,
``GeoFootprint`` (the CDSE GeoJSON Polygon string), ``ref_lon`` and
``ref_lat`` (the reference point, typically inside the footprint).  The
optional ``Mode`` column gives the imaging mode; when it is missing, the
mode is inferred from the SAFE Name (``S1[A-D]_IW_`` / ``S1[A-D]_EW_``).

The output Parquet adds the columns ``sar_incidence_angle``,
``sar_elevation_angle`` (deg) and the intermediate quantities ``mode``,
``dfrac``, ``swath_width_km`` and ``r_s_km``.
"""

import re
import sys
from pathlib import Path

import pandas as pd

from s1_geometry import estimate_incidence_elevation

_MODE_RE = re.compile(r"^S1[A-D]_(IW|EW)_")
EXAMPLE = Path(__file__).parent / "data" / "s1_sample_footprints.parquet"
REQUIRED_COLUMNS = ("Name", "GeoFootprint", "ref_lon", "ref_lat")


def _mode_from_name(name: str) -> str:
    """Imaging mode from the SAFE Name (S1[A-D]_IW_ / S1[A-D]_EW_)."""
    match = _MODE_RE.match(name)
    if not match:
        raise ValueError(f"cannot infer the imaging mode from the SAFE name: {name}")
    return match.group(1)


def add_angles(df: pd.DataFrame, mode_column: str = "Mode") -> pd.DataFrame:
    """Return the input with the angle columns appended."""
    missing = [column for column in REQUIRED_COLUMNS if column not in df.columns]
    if missing:
        raise ValueError(f"missing required column(s) in the input Parquet: {missing}")
    if mode_column in df.columns:
        modes = df[mode_column].astype(str).str.upper().tolist()
    else:
        modes = [_mode_from_name(name) for name in df["Name"]]
    rows = []
    for row, mode in zip(df.itertuples(index=False), modes):
        try:
            result = estimate_incidence_elevation(
                row.GeoFootprint, float(row.ref_lon), float(row.ref_lat), mode)
        except ValueError as error:
            raise RuntimeError(f"{row.Name}: {error}") from error
        rows.append({
            "mode": result.mode,
            "dfrac": result.dfrac,
            "swath_width_km": result.swath_width_km,
            "r_s_km": result.r_s_km,
            "sar_incidence_angle": result.incidence_deg,
            "sar_elevation_angle": result.elevation_deg,
        })
    return pd.concat([df, pd.DataFrame(rows)], axis=1)


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    in_path = Path(argv[0]) if argv else EXAMPLE
    out_path = (
        Path(argv[1])
        if len(argv) >= 2
        else in_path.with_name(f"{in_path.stem}_angles.parquet")
    )
    result = add_angles(pd.read_parquet(in_path))
    result.to_parquet(out_path, index=False)
    print(f"wrote {out_path} ({len(result)} rows)")
    print(result[
        ["Name", "mode", "ref_lon", "ref_lat", "swath_width_km", "r_s_km",
         "dfrac", "sar_incidence_angle", "sar_elevation_angle"]
    ].to_string(index=False, max_colwidth=52))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
