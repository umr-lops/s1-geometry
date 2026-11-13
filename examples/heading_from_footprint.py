"""Estimate the ground heading of a Sentinel-1 IW acquisition from its SAFE footprint.

Usage:

    python examples/heading_from_footprint.py [footprint.geojson] [ref_lat]

Without arguments, a real ascending IW SAFE footprint (CDSE GeoFootprint,
2026-11-09) is used. ``footprint.geojson`` must contain a GeoJSON Polygon
(e.g. the ``GeoFootprint`` column of the CDSE daily Parquets).
"""

import sys
from pathlib import Path

from s1_geometry import estimate_geometry

# Real S1D IW OCN SAFE footprint (CDSE GeoFootprint, 2026-11-09), ascending.
EXAMPLE_ASC_RING = [
    (-76.425255, 37.274273), (-76.00177, 35.512432),
    (-73.232109, 35.912037), (-73.591484, 37.671867),
    (-76.425255, 37.274273),
]


def main() -> int:
    if len(sys.argv) >= 2:
        polygon = Path(sys.argv[1]).read_text()
    else:
        polygon = EXAMPLE_ASC_RING
    ref_lat = float(sys.argv[2]) if len(sys.argv) >= 3 else None

    result = estimate_geometry(polygon, ref_lat=ref_lat)
    print(f"orbit_direction:        {result.orbit_direction}")
    print(f"ref_lat:                {result.ref_lat:.6f} deg")
    print(f"heading_deg:            {result.heading_deg:.3f} deg (clockwise from North)")
    print(f"heading_uncertainty_deg: {result.heading_uncertainty_deg:.3f} deg")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
