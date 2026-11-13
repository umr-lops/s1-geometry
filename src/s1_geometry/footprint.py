"""Sentinel-1 SAR footprint handling: SAFE polygon parsing and orbit direction.

The footprint of a Sentinel-1 IW/EW acquisition is the 4-vertex SAFE polygon
(the ``GeoFootprint`` column of the CDSE daily Parquets, or the geometry of a
matchup STAC item).  The orbit (pass) direction is derived from the ring
here; see :func:`orbit_direction_from_ring` for the rule and the supported
vertex-order conventions.

A ``ring`` / ``polygon`` can be passed anywhere in this package as a sequence
of (lon, lat) vertices, a GeoJSON Polygon *string* or a GeoJSON Polygon
*dict* (the STAC item representation).
"""

from __future__ import annotations

import json
import math
from collections.abc import Iterable, Sequence
from pathlib import Path

import numpy as np

ASCENDING = "ascending"
DESCENDING = "descending"
TURNING = "turning"

#: Flight north component (in deg latitude) below which in magnitude a SAFE
#: polygon is classified as turning (equal to the historical
#: \\|lat(v[1]) - lat(v[0])\\| for the standard CDSE vertex order). Measured on
#: 5012 IW OCN rings (2026-11-01..12): 11 products (0.2 %, all 2SSH at
#: 77.1-77.5 S, at the polar turn) have a magnitude below 0.3 deg, the next
#: smallest is 0.31 deg, and the typical magnitude is ~1.5 deg.
DLAT_TURN_THRESHOLD_DEG = 0.3

#: Kilometres per degree of latitude on the mean-radius (6371 km) sphere that
#: the CDSE catalogue footprint is computed on.
_KM_PER_DEG = 2.0 * math.pi * 6371.0 / 360.0

#: Longest-edge length (km) above which a ring is treated as EW-like: IW
#: footprint edges are 150-274 km (burst / swath) whereas EW edges are
#: 412-486 km.
_EW_LIKE_MAX_EDGE_KM = 330.0

#: Relative length difference below which the two opposite-edge families are
#: considered indistinguishable (near-square rings, e.g. EW and RAW).
_FAMILY_AMBIGUITY_TOL = 0.10

#: An opposite-edge family shorter than this (km) is treated as degenerate.
_DEGEN_FAMILY_KM = 2.0


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


def _edge_vectors_km(points: np.ndarray) -> np.ndarray:
    """The 4 ring edge vectors as a ``(4, 2)`` array in (east, north) km.

    ``e[i] = v[i+1] - v[i]`` (indices modulo 4), with the longitude difference
    scaled by ``cos`` of the ring mean latitude to kilometres on the 6371 km
    mean-radius sphere.
    """
    lat_mid = float(points[:, 1].mean())
    east_scale = _KM_PER_DEG * math.cos(math.radians(lat_mid))
    edges = np.empty((4, 2))
    for i in range(4):
        delta = points[(i + 1) % 4] - points[i]
        edges[i] = (delta[0] * east_scale, delta[1] * _KM_PER_DEG)
    return edges


def orbit_direction_from_ring(
    ring: Sequence[Sequence[float]] | str | dict,
) -> str:
    """Classify the orbit (pass) direction from the SAR polygon (SAFE footprint).

    The ring is the 4-vertex SAFE footprint; a GeoJSON Polygon string or dict
    is accepted as well as a sequence of (lon, lat) vertices.

    The catalogue ``orbitDirection`` is the *pass* (flight) direction, not
    the orbit branch: a northbound pass is :data:`ASCENDING` and a southbound
    pass :data:`DESCENDING` (an EW product at 70 deg S on the northbound
    return leg is still ``ascending``).  The flight direction runs along the
    burst (along-track) edges of the footprint parallelogram.

    Rule
    ----
    1.  The edge vectors ``e0..e3`` are computed in (east, north) km at the
        ring mean latitude.  They form two opposite-edge families, ``{e0, e2}``
        and ``{e1, e3}``.  One family is the pair of **range lines** (the two
        opposite edges parallel to the flight, i.e. the near and far range
        lines); the other is the pair of **burst edges** (perpendicular to the
        flight).
    2.  For IW-like rings the range-line family is the **shorter** opposite-edge
        pair (the along-track burst extent, 150-180 km, is shorter than the
        across-track swath width, 248-274 km).  The flight runs parallel to the
        range lines: it is ``-e0`` when the range-line family is ``{e0, e2}``
        and ``+e1`` when it is ``{e1, e3}``.  The sign follows from the radar
        being left-looking (the near range line is on the left of the flight);
        for a ring in the standard CDSE vertex order this reduces exactly to
        the historical dlat rule.
    3.  EW-like rings (longest edge > 330 km) and rings whose two families are
        near-equal in length (relative difference < 10 %, e.g. near-square EW
        or RAW) are ambiguous: the standard CDSE vertex order is assumed and
        the flight is ``-e0`` (the historical dlat rule).
    4.  The pass is :data:`ASCENDING` when the flight's north component
        exceeds +:data:`DLAT_TURN_THRESHOLD_DEG` deg latitude,
        :data:`DESCENDING` when it is below the negative threshold, and
        :data:`TURNING` otherwise (polar turnaround).

    Vertex-order conventions
    ------------------------
    The flight direction is two-valued from the ring alone: the ring does not
    label the near/far or start/end corners, so the two opposite along-track
    edge directions are physically indistinguishable (a Z2 symmetry).  The
    rule resolves the ambiguity correctly for the two documented vertex
    orders:

    * the **standard CDSE** order (the ``GeoFootprint`` column of the daily
      Parquets, counter-clockwise), and
    * the **product-page** order (the clockwise ring copied from the CDSE
      product page).

    Validated against the ``orbitDirection`` attribute of 66243 IW/EW rings
    (13 daily v2 Parquets, 2026-11-01..13): IW 99.70 %, EW 92.05 % (both
    unchanged by this rule), ETA 97.82 % (the historical dlat rule gave
    55.66 %, because ETA rings use a rotated vertex order).  A ring in a
    vertex order other than the two documented ones (e.g. a fully reversed
    ring) may be misclassified, since the opposite flight direction is not
    distinguishable from the ring alone.
    """
    points = _as_point_array(ring)
    edges = _edge_vectors_km(points[:4])
    lengths = np.linalg.norm(edges, axis=1)
    len02 = float(lengths[0] + lengths[2])
    len13 = float(lengths[1] + lengths[3])
    if (
        min(len02, len13) < _DEGEN_FAMILY_KM
        or float(lengths.max()) > _EW_LIKE_MAX_EDGE_KM
        or abs(len02 - len13) < _FAMILY_AMBIGUITY_TOL * max(len02, len13)
    ):
        # Degenerate, EW-like or near-square ring: the along-track family is
        # not identifiable, assume the standard CDSE vertex order.
        flight = -edges[0]
    elif len02 < len13:
        # {e0, e2} is the range-line family (parallel to the flight).
        flight = -edges[0]
    else:
        # {e1, e3} is the range-line family (parallel to the flight).
        flight = edges[1]
    dlat = float(flight[1]) / _KM_PER_DEG
    if dlat > DLAT_TURN_THRESHOLD_DEG:
        return ASCENDING
    if dlat < -DLAT_TURN_THRESHOLD_DEG:
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
