"""Generic estimation entry point for Sentinel-1 SAR footprint geometry.

New geometric quantities (incidence angle, elevation angle, ...) are added
to the :class:`S1Geometry` dataclass and computed by
:func:`estimate_geometry`; shared footprint handling stays in
:mod:`s1_geometry.footprint`, the incidence / elevation estimation in
:mod:`s1_geometry.angles`.
"""

from __future__ import annotations

from dataclasses import dataclass

from s1_geometry.angles import estimate_incidence_elevation
from s1_geometry.footprint import _as_point_array, orbit_direction_from_ring
from s1_geometry.heading import heading_from_lat, sar_ground_heading_uncertainty


@dataclass(frozen=True)
class S1Geometry:
    """Geometric quantities estimated from a Sentinel-1 SAFE footprint.

    Attributes
    ----------
    orbit_direction:
        ``"ascending" | "descending" | "turning"`` (from the footprint).
    ref_lat:
        Reference latitude (deg) at which the heading was evaluated.
    heading_deg:
        Ground heading (deg, clockwise from geographic North) at
        :attr:`ref_lat` (the IW model; a documented extrapolation for EW).
    heading_uncertainty_deg:
        Latitude- and hemisphere-dependent heading uncertainty (deg); users
        decide instead of rejecting high-latitude acquisitions.
    mode:
        Imaging mode (``"IW"`` or ``"EW"``) when the incidence / elevation
        angles were estimated, otherwise ``None``.
    incidence_deg:
        Incidence angle (deg, from the local vertical) at the reference
        point, or ``None`` when :attr:`mode` is ``None``.
    elevation_deg:
        Elevation angle (``90 - incidence``) at the reference point, or
        ``None`` when :attr:`mode` is ``None``.
    """

    orbit_direction: str
    ref_lat: float
    heading_deg: float
    heading_uncertainty_deg: float
    mode: str | None = None
    incidence_deg: float | None = None
    elevation_deg: float | None = None


def estimate_geometry(
    sar_polygon,
    ref_lat: float | None = None,
    mode: str | None = None,
    ref_lon: float | None = None,
) -> S1Geometry:
    """Estimate the geometric quantities of a Sentinel-1 acquisition.

    Parameters
    ----------
    sar_polygon:
        The SAFE footprint: a sequence of (lon, lat) vertices, a GeoJSON
        Polygon string (the CDSE ``GeoFootprint`` column) or a GeoJSON
        Polygon dict (the matchup STAC item representation).
    ref_lat:
        Reference latitude (deg) at which the heading is evaluated. Defaults
        to the mean latitude of the 4 footprint vertices.
    mode:
        Imaging mode, ``"IW"`` or ``"EW"``.  When given, the incidence and
        elevation angles are estimated at the reference point; the heading
        fields remain the IW model (a documented extrapolation for EW).
    ref_lon:
        Reference longitude (deg) of the point at which the incidence /
        elevation angles are evaluated.  Must be given together with
        ``ref_lat``; when both are omitted and ``mode`` is given, the
        footprint center (mean of the 4 vertices) is used.

    Raises
    ------
    ValueError
        If exactly one of ``ref_lon`` / ``ref_lat`` is given, the mode is
        unknown, or the footprint width does not invert to an orbit altitude
        in the documented band (see
        :func:`~s1_geometry.angles.estimate_incidence_elevation`).
    """
    corners = _as_point_array(sar_polygon)[:4]
    if mode is not None and (ref_lon is None) != (ref_lat is None):
        raise ValueError("ref_lon and ref_lat must be given together")
    if ref_lat is None:
        ref_lat = float(corners[:, 1].mean())
    if mode is not None and ref_lon is None:
        ref_lon = float(corners[:, 0].mean())
    orbit = orbit_direction_from_ring(sar_polygon)
    geometry = S1Geometry(
        orbit_direction=orbit,
        ref_lat=float(ref_lat),
        heading_deg=heading_from_lat(ref_lat, orbit),
        heading_uncertainty_deg=sar_ground_heading_uncertainty(ref_lat),
    )
    if mode is None:
        return geometry
    angles = estimate_incidence_elevation(sar_polygon, ref_lon, ref_lat, mode)
    return S1Geometry(
        orbit_direction=orbit,
        ref_lat=float(ref_lat),
        heading_deg=geometry.heading_deg,
        heading_uncertainty_deg=geometry.heading_uncertainty_deg,
        mode=angles.mode,
        incidence_deg=angles.incidence_deg,
        elevation_deg=angles.elevation_deg,
    )
