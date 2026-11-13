"""Generic estimation entry point for Sentinel-1 SAR footprint geometry.

New geometric quantities (incidence angle, elevation angle, ...) are added
to the :class:`S1Geometry` dataclass and computed by
:func:`estimate_geometry`; shared footprint handling stays in
:mod:`s1_geometry.footprint`.
"""

from __future__ import annotations

from dataclasses import dataclass

from s1_geometry.footprint import _as_point_array, orbit_direction_from_ring
from s1_geometry.heading import heading_from_lat, sar_ground_heading_uncertainty


@dataclass(frozen=True)
class S1Geometry:
    """Geometric quantities estimated from a Sentinel-1 IW SAFE footprint.

    Attributes
    ----------
    orbit_direction:
        ``"ascending" | "descending" | "turning"`` (from the footprint).
    ref_lat:
        Reference latitude (deg) at which the heading was evaluated.
    heading_deg:
        Ground heading (deg, clockwise from geographic North) at
        :attr:`ref_lat`.
    heading_uncertainty_deg:
        Latitude- and hemisphere-dependent heading uncertainty (deg); users
        decide instead of rejecting high-latitude acquisitions.
    """

    orbit_direction: str
    ref_lat: float
    heading_deg: float
    heading_uncertainty_deg: float


def estimate_geometry(sar_polygon, ref_lat: float | None = None) -> S1Geometry:
    """Estimate the geometric quantities of a Sentinel-1 IW acquisition.

    Parameters
    ----------
    sar_polygon:
        The SAFE footprint: a sequence of (lon, lat) vertices, a GeoJSON
        Polygon string (the CDSE ``GeoFootprint`` column) or a GeoJSON
        Polygon dict (the matchup STAC item representation).
    ref_lat:
        Reference latitude (deg) at which the heading is evaluated. Defaults
        to the mean latitude of the 4 footprint vertices.
    """
    if ref_lat is None:
        ref_lat = float(_as_point_array(sar_polygon)[:4, 1].mean())
    orbit = orbit_direction_from_ring(sar_polygon)
    return S1Geometry(
        orbit_direction=orbit,
        ref_lat=float(ref_lat),
        heading_deg=heading_from_lat(ref_lat, orbit),
        heading_uncertainty_deg=sar_ground_heading_uncertainty(ref_lat),
    )
