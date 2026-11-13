"""Command-line interface: ``s1-geometry``.

Estimate the ground heading of a Sentinel-1 IW acquisition from its SAFE
footprint polygon::

    s1-geometry heading footprint.geojson --lat 36.59
    s1-geometry heading '{"type": "Polygon", "coordinates": [[...]]}' --lat 36.59 --json

The polygon argument is a GeoJSON file path or an inline JSON string: a
GeoJSON Polygon object, a GeoJSON Feature, or a bare ring of [lon, lat]
vertices.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from s1_geometry import estimate_geometry


def _load_polygon(source: str) -> dict:
    """Read the polygon from a GeoJSON file path or an inline JSON string.

    Accepts a GeoJSON Polygon object, a GeoJSON Feature with a Polygon
    geometry, or a bare ring (a list of [lon, lat] vertices).
    """
    path = Path(source)
    text = path.read_text() if path.is_file() else source
    data = json.loads(text)
    if isinstance(data, dict) and data.get("type") == "Polygon":
        return data
    if isinstance(data, dict) and data.get("type") == "Feature":
        return data["geometry"]
    if isinstance(data, list) and data and isinstance(data[0], (list, tuple)):
        return {"type": "Polygon", "coordinates": [data]}
    raise ValueError(
        "expected a GeoJSON Polygon object, a GeoJSON Feature, or a ring of [lon, lat] vertices"
    )


def main(argv: list[str] | None = None) -> int:
    """Entry point of the ``s1-geometry`` console script."""
    parser = argparse.ArgumentParser(
        prog="s1-geometry",
        description="Estimate geometric quantities from a Sentinel-1 SAR footprint.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    heading = sub.add_parser(
        "heading", help="estimate the ground heading from a SAFE footprint polygon",
    )
    heading.add_argument(
        "polygon",
        help="GeoJSON Polygon file path or inline JSON: a Polygon object, a "
             "Feature, or a bare ring of [lon, lat] vertices",
    )
    heading.add_argument(
        "--lat", type=float, required=True,
        help="reference latitude in degrees (the heading is evaluated there)",
    )
    heading.add_argument(
        "--json", action="store_true", help="print the result as a JSON object",
    )

    args = parser.parse_args(argv)

    if args.command == "heading":
        try:
            polygon = _load_polygon(args.polygon)
        except (OSError, ValueError, json.JSONDecodeError) as error:
            parser.error(f"cannot read the polygon: {error}")
        result = estimate_geometry(polygon, ref_lat=args.lat)
        if args.json:
            print(json.dumps({
                "orbit_direction": result.orbit_direction,
                "ref_lat": result.ref_lat,
                "heading_deg": result.heading_deg,
                "heading_uncertainty_deg": result.heading_uncertainty_deg,
            }, indent=2))
        else:
            print(f"orbit_direction: {result.orbit_direction}")
            print(f"ref_lat: {result.ref_lat:.6f}")
            print(f"heading_deg: {result.heading_deg:.3f}")
            print(f"heading_uncertainty_deg: {result.heading_uncertainty_deg:.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
