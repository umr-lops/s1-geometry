"""Sentinel-1 SAR footprint handling: SAFE polygon parsing and orbit direction.

The footprint of a Sentinel-1 IW/EW acquisition is the 4-vertex SAFE polygon
(the ``GeoFootprint`` column of the CDSE daily Parquets, or the geometry of a
matchup STAC item).  A WV acquisition instead imaged two subswaths in a
sequence of bursts, so its ``GeoFootprint`` is a *MultiPolygon* of one small
ring per burst; see :func:`as_multipolygon` and :mod:`s1_geometry.wv`.  The
orbit (pass) direction of a ring is derived here; see
:func:`orbit_direction_from_ring` for the rule and the supported vertex-order
conventions.

A ``ring`` / ``polygon`` can be passed anywhere in this package as a sequence
of (lon, lat) vertices, a GeoJSON Polygon *string* or a GeoJSON Polygon
*dict* (the STAC item representation); a MultiPolygon is a GeoJSON
MultiPolygon string / dict or a sequence of such rings.
"""

from __future__ import annotations

import json
import math
import re
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

#: Across-track family length (km) band of the SM (stripmap) beams.  The
#: attitude-steering law keeps every SM beam's ground swath fixed at ~80 km
#: (independent of the orbit altitude), so an SM SAFE ring always carries one
#: opposite-edge family of ~80 km (the across-track burst edges); the other
#: family (the along-track range lines) is much longer (a few hundred km).
_SM_BAND_LO_KM = 72.0
_SM_BAND_HI_KM = 88.0

#: Longest-edge length (km) below which a ring is too small to be a SAFE
#: swath footprint (e.g. a ~20 km WV imagette): its orbit direction is not
#: identifiable from the ring alone.
_MIN_SWATH_EDGE_KM = 40.0

#: Long / short opposite-edge family ratio below which an SM-like ring is
#: treated as near-square instead (an SM burst shorter than ~15 s is nearly
#: square and falls back to the standard-CDSE-order assumption).
_SM_MIN_ASPECT_RATIO = 1.2


_WKT_PREFIX = re.compile(r"^\s*(?:geography\s*')?\s*(?:SRID=\d+\s*;)?\s*", re.IGNORECASE)
_WKT_PAIR = re.compile(r"([-+]?[\d.]+(?:[eE][-+]?\d+)?)\s+([-+]?[\d.]+(?:[eE][-+]?\d+)?)")


def parse_wkt(text: str) -> dict:
    """Parse a WKT / EWKT ``POLYGON`` or ``MULTIPOLYGON`` into a GeoJSON dict.

    Accepts the CDSE OData ``Footprint`` form
    ``geography'SRID=4326;MULTIPOLYGON (((lon lat, ...)), ...)'`` as well as
    plain (E)WKT.  Only 2D coordinates are supported.
    """
    body = _WKT_PREFIX.sub("", text).strip().rstrip("'").strip()
    match = re.match(r"^(MULTIPOLYGON|POLYGON)\s*(\(.*\))$", body, re.IGNORECASE | re.DOTALL)
    if match is None:
        raise ValueError(f"not a WKT POLYGON / MULTIPOLYGON: {text[:80]!r}")
    kind, coords = match.group(1).upper(), match.group(2)
    as_json = _WKT_PAIR.sub(r"[\1,\2]", coords).replace("(", "[").replace(")", "]")
    try:
        nested = json.loads(as_json)
    except json.JSONDecodeError as error:
        raise ValueError(f"malformed WKT coordinates: {error}") from error
    return {"type": "MultiPolygon" if kind == "MULTIPOLYGON" else "Polygon",
            "coordinates": nested}


def _loads_geometry(text: str):
    """GeoJSON or WKT string -> Python object."""
    stripped = text.lstrip()
    if stripped[:1] in ("{", "["):
        return json.loads(stripped)
    return parse_wkt(stripped)


def _as_point_array(ring: Sequence[Sequence[float]] | str | dict) -> np.ndarray:
    """Validate a polygon ring and return it as an (N, 2) (lon, lat) array."""
    if isinstance(ring, str):
        try:
            ring = _loads_geometry(ring)
        except json.JSONDecodeError as error:
            raise ValueError(f"ring is not a valid GeoJSON Polygon string: {error}") from error
    if isinstance(ring, dict):
        if ring.get("type") == "MultiPolygon":
            raise ValueError(
                "expected a Polygon, got a MultiPolygon: a WV footprint is handled by "
                "s1_geometry.wv.estimate_wv_angles (or estimate_geometry(mode='WV'))")
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


def as_multipolygon(
    geom: Sequence[Sequence[Sequence[float]]] | str | dict,
) -> list[np.ndarray]:
    """Validate a MultiPolygon and return its (lon, lat) sub-polygon rings.

    ``geom`` is a GeoJSON MultiPolygon *string* or *dict*, a WKT / EWKT
    ``MULTIPOLYGON`` string (e.g. the CDSE OData ``Footprint``), or a sequence of
    rings (each a sequence of (lon, lat) vertices).  Returns a list of
    ``(N_i, 2)`` (lon, lat) arrays, one per sub-polygon, in the given order.

    For a WV acquisition the CDSE ``GeoFootprint`` is a MultiPolygon of one
    small (~14 km) ring per burst; the ring order is the burst (time) order,
    alternating between the two subswaths (see :mod:`s1_geometry.wv`).
    """
    if isinstance(geom, str):
        try:
            geom = _loads_geometry(geom)
        except json.JSONDecodeError as error:
            raise ValueError(
                f"not a valid GeoJSON MultiPolygon string: {error}") from error
    if isinstance(geom, dict):
        if geom.get("type") != "MultiPolygon":
            raise ValueError(
                f"expected a GeoJSON MultiPolygon, got type {geom.get('type')!r}")
        coords = [poly[0] for poly in geom["coordinates"]]
    else:
        try:
            coords = list(geom)
        except TypeError as error:
            raise ValueError(f"cannot interpret {geom!r} as a MultiPolygon") from error
    if not coords:
        raise ValueError("MultiPolygon has no sub-polygons")
    rings: list[np.ndarray] = []
    for i, ring in enumerate(coords):
        points = np.asarray(ring, dtype=float)
        if points.ndim != 2 or points.shape[1] != 2 or points.shape[0] < 3:
            raise ValueError(
                f"sub-polygon {i} must have at least 3 (lon, lat) vertices, "
                f"got {ring!r}")
        if not np.isfinite(points).all():
            raise ValueError(f"sub-polygon {i} contains non-finite coordinates")
        if np.abs(points[:, 1]).max() > 90.0:
            raise ValueError(f"sub-polygon {i} latitude out of range")
        rings.append(points)
    return rings


def _on_antimeridian(ring: np.ndarray) -> np.ndarray:
    return np.abs(np.abs(ring[:, 0]) - 180.0) < 1e-6


def merge_antimeridian_splits(rings: list[np.ndarray]) -> list[np.ndarray]:
    """Re-join the WV imagette rings split at the antimeridian.

    The CDSE OData ``Footprint`` cuts an imagette crossing +/-180 deg into two
    consecutive polygons (one on each side), so the product has one ring more
    than imagettes and the WV1 / WV2 alternation is broken.  Two consecutive
    rings that both touch the antimeridian on opposite sides, with
    overlapping latitude ranges, are merged: the cut vertices (on +/-180 deg)
    are dropped and the remaining corners (first piece's, then second
    piece's, in their original order) are unwrapped next to the first piece's
    longitude.  The merged ring may have longitudes beyond +/-180 deg,
    which every computation of :mod:`s1_geometry.wv` (ECEF) handles.
    """
    out: list[np.ndarray] = []
    i = 0
    while i < len(rings):
        a = rings[i]
        if i + 1 < len(rings):
            b = rings[i + 1]
            on_a, on_b = _on_antimeridian(a), _on_antimeridian(b)
            if (
                on_a.any() and on_b.any()
                and np.sign(a[on_a, 0][0]) != np.sign(b[on_b, 0][0])
                and a[:, 1].min() <= b[:, 1].max() and b[:, 1].min() <= a[:, 1].max()
            ):
                side = np.sign(a[on_a, 0][0])  # the first piece is at side * 180
                pts = np.vstack([a[~on_a], b[~on_b]])
                pts[:, 0] = np.where(np.sign(pts[:, 0]) == side, pts[:, 0],
                                     pts[:, 0] + side * 360.0)
                _, first = np.unique(pts.round(9), axis=0, return_index=True)
                pts = pts[np.sort(first)]
                out.append(np.vstack([pts, pts[:1]]))
                i += 2
                continue
        out.append(a)
        i += 1
    return out


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
        being right-looking (the near range line is on the right of the
        flight); for a ring in the standard CDSE vertex order this reduces
        exactly to the historical dlat rule.
    3.  For SM (stripmap) rings the situation is reversed: the attitude-
        steering law keeps the beam's ground swath fixed at ~80 km, so the
        across-track burst edges are the **shorter** family and the range
        lines are the **longer** one (a few hundred km).  The SM ring is
        detected by one family in the 72-88 km band and a long / short ratio
        of at least 1.2; the flight is ``-e0`` when the range lines are
        ``{e0, e2}`` (standard CDSE order) and ``-e1`` when they are
        ``{e1, e3}`` (product-page order).
    4.  EW-like rings (longest edge > 330 km) and rings whose two families are
        near-equal in length (relative difference < 10 %, e.g. a near-square
        IW or EW ring) are ambiguous from the edge lengths alone: the
        range-line family is decided from the two candidate across-track
        widths (see :func:`_orbit_range_line_family` in :mod:`angles`).  When
        the widths are clearly separated (e.g. an EW ring) the width-inversion
        family (the same discriminator as :func:`swath_width_km`) is used; when
        they are a near-square knife-edge (a few km apart, e.g. a near-square
        IW ring) the flight-parallel family is the one with the larger
        north-south extent.  When the width does not invert to an orbit
        altitude (degenerate or partial-swath ring) the standard CDSE vertex
        order is assumed and the flight is ``-e0`` (the historical dlat rule).
    5.  The pass is :data:`ASCENDING` when the flight's north component
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
    short_km = 0.5 * min(len02, len13)
    long_km = 0.5 * max(len02, len13)
    if float(lengths.max()) < _MIN_SWATH_EDGE_KM:
        raise ValueError(
            f"footprint ring too small for a SAFE swath (longest edge "
            f"{float(lengths.max()):.1f} km < {_MIN_SWATH_EDGE_KM:.0f} km, e.g. a WV "
            "imagette): the orbit direction is not identifiable from it; use the "
            "orbitDirection metadata, or the WV product footprint")
    if (
        _SM_BAND_LO_KM <= short_km <= _SM_BAND_HI_KM
        and long_km >= _SM_MIN_ASPECT_RATIO * short_km
    ):
        # SM (stripmap) ring: the across-track burst edges are the ~80 km
        # family, so the range lines are the LONG opposite-edge pair (the
        # opposite of IW / EW, whose burst extent is shorter than the swath).
        # In the standard CDSE (counter-clockwise) order the range lines are
        # {e0, e2} (flight = -e0); in the product-page (clockwise) order they
        # are {e1, e3} (flight = -e1).  Note: a partial-swath IW ring whose
        # width falls in the SM band (e.g. an ~80 km IW30) in the product-page
        # order is misclassified (the flight direction is flipped).
        flight = -edges[0] if len02 > len13 else -edges[1]
    elif (
        min(len02, len13) < _DEGEN_FAMILY_KM
        or float(lengths.max()) > _EW_LIKE_MAX_EDGE_KM
        or abs(len02 - len13) < _FAMILY_AMBIGUITY_TOL * max(len02, len13)
    ):
        # Degenerate, EW-like or near-square ring: the along-track family is
        # not identifiable from the edge lengths alone.  Decide it by the
        # across-track widths (see :func:`_orbit_range_line_family`): width
        # inversion when the two family widths are clearly separated (e.g. EW),
        # else the flight-parallel (north-south) family for a near-square
        # knife-edge.  Fall back to the standard CDSE vertex order (flight =
        # -e0) when the width does not invert to an orbit altitude (degenerate
        # / partial-swath ring).  The import is lazy because angles.py imports
        # footprint.py at module level.
        try:
            from s1_geometry.angles import _orbit_range_line_family
            family = _orbit_range_line_family(ring)
        except Exception:  # keep the historical -e0 behaviour, never crash
            family = "standard"
        flight = -edges[0] if family == "standard" else edges[1]
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
        for name, geo in zip(names, geos, strict=True):
            try:
                points = _as_point_array(geo)
                if tuple(points[-1]) != tuple(points[0]):
                    raise ValueError(f"ring not closed: {name}")
            except (ValueError, KeyError, TypeError, json.JSONDecodeError):
                continue
            polygons[str(name)] = (points, ring_to_geojson(points))
    return polygons
