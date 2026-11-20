"""Command-line interface: ``s1-geometry``.

Estimate the ground heading of a Sentinel-1 IW acquisition from its SAFE
footprint polygon::

    s1-geometry heading footprint.geojson --lat 36.59
    s1-geometry heading '{"type": "Polygon", "coordinates": [[...]]}' --lat 36.59 --json

or the incidence / elevation angles of an IW or EW acquisition at a
reference point::

    s1-geometry angles footprint.geojson --lon 10.5 --lat 36.59 --mode IW
    s1-geometry angles footprint.geojson --lon 10.5 --lat 36.59 --json

or the incidence / elevation / heading of a WV (Wave mode) acquisition at a
reference point, from the product footprint only (OData metadata)::

    s1-geometry wv "geography'SRID=4326;MULTIPOLYGON (((...)), ...)'" --lon -76.7 --lat -29.3

The polygon argument is a GeoJSON file path or an inline JSON string: a
GeoJSON Polygon object, a GeoJSON Feature, or a bare ring of [lon, lat]
vertices; for ``wv`` a GeoJSON MultiPolygon (object or Feature) or a WKT
``MULTIPOLYGON`` (the CDSE OData ``Footprint``).
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from s1_geometry import estimate_geometry, estimate_incidence_elevation
from s1_geometry.footprint import parse_wkt


def _load_polygon(source: str) -> dict:
    """Read the polygon from a GeoJSON file path or an inline JSON string.

    Accepts a GeoJSON Polygon object, a GeoJSON Feature with a Polygon
    geometry, or a bare ring (a list of [lon, lat] vertices).
    """
    try:
        is_file = Path(source).is_file()
    except OSError:  # an inline polygon longer than the OS path limit
        is_file = False
    text = Path(source).read_text() if is_file else source
    if text.lstrip()[:1] not in ("{", "["):
        return parse_wkt(text)
    data = json.loads(text)
    if isinstance(data, dict) and data.get("type") in ("Polygon", "MultiPolygon"):
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

    angles = sub.add_parser(
        "angles",
        help="estimate the incidence / elevation angles at a reference point",
    )
    angles.add_argument(
        "polygon",
        help="GeoJSON Polygon file path or inline JSON: a Polygon object, a "
             "Feature, or a bare ring of [lon, lat] vertices",
    )
    angles.add_argument(
        "--lon", type=float, required=True,
        help="reference longitude in degrees",
    )
    angles.add_argument(
        "--lat", type=float, required=True,
        help="reference latitude in degrees",
    )
    angles.add_argument(
        "--mode", choices=("IW", "EW", "SM"), default="IW",
        help="imaging mode (default: IW)",
    )
    angles.add_argument(
        "--beam", choices=tuple(f"S{i}" for i in range(1, 7)), default=None,
        help="SM beam S1..S6 (default: S3; used only with --mode SM)",
    )
    angles.add_argument(
        "--json", action="store_true", help="print the result as a JSON object",
    )

    wv = sub.add_parser(
        "wv",
        help="estimate the WV incidence / elevation / heading at a reference point",
    )
    wv.add_argument(
        "polygon",
        help="WV product footprint: GeoJSON MultiPolygon file path or inline JSON, "
             "or WKT MULTIPOLYGON (the CDSE OData Footprint)",
    )
    wv.add_argument("--lon", type=float, required=True, help="reference longitude in degrees")
    wv.add_argument("--lat", type=float, required=True, help="reference latitude in degrees")
    wv.add_argument(
        "--r-s", type=float, default=None, dest="r_s",
        help="satellite geocentric distance in km (default: metadata-only orbit model)",
    )
    wv.add_argument(
        "--json", action="store_true", help="print the result as a JSON object",
    )

    args = parser.parse_args(argv)

    if args.command == "heading":
        try:
            polygon = _load_polygon(args.polygon)
            result = estimate_geometry(polygon, ref_lat=args.lat)
        except (OSError, ValueError, json.JSONDecodeError) as error:
            parser.error(f"cannot estimate the heading: {error}")
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
    elif args.command == "angles":
        try:
            polygon = _load_polygon(args.polygon)
            result = estimate_incidence_elevation(
                polygon, args.lon, args.lat, args.mode, args.beam)
        except (OSError, ValueError, json.JSONDecodeError) as error:
            parser.error(f"cannot estimate the angles: {error}")
        if args.json:
            print(json.dumps({
                "mode": result.mode,
                "beam": result.beam,
                "ref_lon": result.ref_lon,
                "ref_lat": result.ref_lat,
                "dfrac": result.dfrac,
                "swath_width_km": result.swath_width_km,
                "r_s_km": result.r_s_km,
                "incidence_deg": result.incidence_deg,
                "elevation_deg": result.elevation_deg,
            }, indent=2))
        else:
            print(f"mode: {result.mode}")
            if result.beam:
                print(f"beam: {result.beam}")
            print(f"ref_lon: {result.ref_lon:.6f}")
            print(f"ref_lat: {result.ref_lat:.6f}")
            print(f"dfrac: {result.dfrac:.4f}")
            print(f"swath_width_km: {result.swath_width_km:.2f}")
            print(f"r_s_km: {result.r_s_km:.2f}")
            print(f"incidence_deg: {result.incidence_deg:.3f}")
            print(f"elevation_deg: {result.elevation_deg:.3f}")
    elif args.command == "wv":
        try:
            polygon = _load_polygon(args.polygon)
            result = estimate_geometry(polygon, ref_lat=args.lat, ref_lon=args.lon,
                                       mode="WV", r_s_km=args.r_s)
        except (OSError, ValueError, json.JSONDecodeError) as error:
            parser.error(f"cannot estimate the WV angles: {error}")
        fields = {
            "mode": result.mode,
            "ref_lon": args.lon,
            "ref_lat": result.ref_lat,
            "subswath": result.subswath,
            "orbit_direction": result.orbit_direction,
            "r_s_km": result.r_s_km,
            "incidence_deg": result.incidence_deg,
            "elevation_deg": result.elevation_deg,
            "heading_deg": result.heading_deg,
        }
        if args.json:
            print(json.dumps(fields, indent=2))
        else:
            for name, value in fields.items():
                print(f"{name}: {value:.3f}" if isinstance(value, float) else f"{name}: {value}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
