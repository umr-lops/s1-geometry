"""Sentinel-1 SAR footprint handling: SAFE polygon parsing and orbit direction.

The footprint of a Sentinel-1 IW acquisition is the 4-vertex SAFE polygon
(the ``GeoFootprint`` column of the CDSE daily Parquets, or the geometry of a
matchup STAC item). Its vertex order is fixed and independent of the orbit
direction; the orbit direction is derived from it here.

A ``ring`` / ``polygon`` can be passed anywhere in this package as a sequence
of (lon, lat) vertices, a GeoJSON Polygon *string* or a GeoJSON Polygon
*dict* (the STAC item representation).
"""

from __future__ import annotations

import json
from collections.abc import Iterable, Sequence
from pathlib import Path

import numpy as np

ASCENDING = "ascending"
DESCENDING = "descending"
TURNING = "turning"

#: \\|lat(v[1]) - lat(v[0])\\| below which a SAFE polygon is classified as
#: turning. Measured on 5012 IW OCN rings (2026-11-01..12): 11 products
#: (0.2 %, all 2SSH at 77.1-77.5 S, at the polar turn) have \\|dlat\\| < 0.3 deg,
#: the next smallest is 0.31 deg, and typical \\|dlat\\| is ~1.5 deg.
DLAT_TURN_THRESHOLD_DEG = 0.3


def _as_point_array(ring: Sequence[Sequence[float]] | str | dict) -> np.ndarray:
    """Validate a polygon ring and return it as an (N, 2) (lon, lat) array."""
    if isinstance(ring, str):
        try:
            ring = json.loads(ring)
        except json.JSONDecodeError as error:
            raise ValueError(f"ring is not a valid GeoJSON Polygon string: {error}") from error
    if isinstance(ring, dict):
        if ring.get("type") != "Polygon":
            raise ValueError(f"expected a GeoJSON Polygon, got type {ring.get('type')!r}")
        ring = ring["coordinates"][0]
    points = np.asarray(ring, dtype=float)
    if points.ndim != 2 or points.shape[1] != 2 or points.shape[0] < 4:
        raise ValueError(f"ring must have at least 4 (lon, lat) vertices, got {ring!r}")
    if not np.isfinite(points).all():
        raise ValueError(f"ring contains non-finite coordinates: {ring!r}")
    if np.abs(points[:, 1]).max() > 90.0:
        raise ValueError(f"ring latitude out of range: {ring!r}")
    return points


def _closed_ring_vertices(points: np.ndarray) -> list[tuple[float, float]]:
    """The first 4 ring vertices, closed, as (lon, lat) tuples."""
    vertices = [tuple(vertex) for vertex in points[:4]]
    if vertices[-1] != vertices[0]:
        vertices.append(vertices[0])
    return vertices


def orbit_direction_from_ring(
    ring: Sequence[Sequence[float]] | str | dict,
) -> str:
    """Classify the orbit direction from the SAR polygon (SAFE footprint).

    The ring is the 4-vertex SAFE footprint in the fixed CDSE vertex order
    ``[bottom-left, top-left, top-right, bottom-right]`` (first vertex
    repeated at the end); a GeoJSON Polygon string or dict is accepted.

    With ``dlat = lat(v[1]) - lat(v[0])`` (top-left minus bottom-left):

    * ``dlat < -0.3 deg`` -> :data:`ASCENDING`
    * ``dlat > +0.3 deg`` -> :data:`DESCENDING`
    * otherwise           -> :data:`TURNING`

    Validated on 3844 CDSE GeoFootprint rings with measured OWI headings
    (2026-11-01..12, S1C+S1D, both branches): 99.77 % agreement with the
    orbit relabel; all exceptions are the \\|dlat\\| < 0.3 deg turning products
    at -77.1..-77.5 deg S, correctly classified as turning.
    """
    points = _as_point_array(ring)
    dlat = points[1, 1] - points[0, 1]
    if dlat < -DLAT_TURN_THRESHOLD_DEG:
        return ASCENDING
    if dlat > DLAT_TURN_THRESHOLD_DEG:
        return DESCENDING
    return TURNING


def point_in_polygon(
    lon: float, lat: float, polygon: Sequence[Sequence[float]] | str | dict,
) -> bool:
    """Ray-casting point-in-polygon test against a SAFE footprint."""
    points = _as_point_array(polygon)
    if tuple(points[-1]) == tuple(points[0]):
        points = points[:-1]
    inside = False
    x, y = float(lon), float(lat)
    for index in range(len(points)):
        x1, y1 = points[index]
        x2, y2 = points[(index + 1) % len(points)]
        if (y1 > y) != (y2 > y):
            x_cross = (x2 - x1) * (y - y1) / (y2 - y1) + x1
            if x < x_cross:
                inside = not inside
    return inside


def polygon_index_for_point(
    lon: float, lat: float, polygons: Iterable[Sequence[Sequence[float]] | str | dict],
) -> int | None:
    """Index of the first SAFE polygon containing the point, else None."""
    for index, polygon in enumerate(polygons):
        if point_in_polygon(lon, lat, polygon):
            return index
    return None


def ring_to_geojson(ring: Sequence[Sequence[float]] | str | dict) -> str:
    """Serialize the 4-vertex SAFE ring as a closed compact GeoJSON Polygon."""
    closed = _closed_ring_vertices(_as_point_array(ring))
    return json.dumps(
        {"type": "Polygon", "coordinates": [[[float(lon), float(lat)] for lon, lat in closed]]},
        separators=(",", ":"),
    )


def load_sar_polygons(
    paths: Iterable[str | Path],
) -> dict[str, tuple[np.ndarray, str]]:
    """Index IW products of CDSE daily Parquets (v2) by SAFE basename.

    Returns ``SAFE basename -> (4-vertex (lon, lat) ring, compact GeoJSON
    string)`` for every IW product whose ``GeoFootprint`` parses to a valid
    closed ring. WV rows (and malformed rings) are skipped. This is the
    polygon source until the IW matchup STAC files are available.

    Requires the optional ``parquet`` extra (pandas + pyarrow).
    """
    import pyarrow.parquet as pq

    polygons: dict[str, tuple[np.ndarray, str]] = {}
    for path in paths:
        frame = pq.read_table(path, columns=["Name", "GeoFootprint"]).to_pandas()
        iw = frame["Name"].str.fullmatch(r"S1[A-D]_IW_.*\.SAFE", na=False)
        names = frame["Name"].to_numpy()[iw]
        geos = frame["GeoFootprint"].to_numpy()[iw]
        for name, geo in zip(names, geos):
            try:
                points = _as_point_array(geo)
                if tuple(points[-1]) != tuple(points[0]):
                    raise ValueError(f"ring not closed: {name}")
            except (ValueError, KeyError, TypeError, json.JSONDecodeError):
                continue
            polygons[str(name)] = (points, ring_to_geojson(points))
    return polygons
