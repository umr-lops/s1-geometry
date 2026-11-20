"""Generic estimation entry point for Sentinel-1 SAR footprint geometry.

New geometric quantities (incidence angle, elevation angle, ...) are added
to the :class:`S1Geometry` dataclass and computed by
:func:`estimate_geometry`; shared footprint handling stays in
:mod:`s1_geometry.footprint`, the incidence / elevation estimation in
:mod:`s1_geometry.angles`.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from s1_geometry.angles import estimate_incidence_elevation
from s1_geometry.footprint import (
    ASCENDING,
    DESCENDING,
    _as_point_array,
    orbit_direction_from_ring,
)
from s1_geometry.heading import heading_from_lat, sar_ground_heading_uncertainty
from s1_geometry.wv import WV_HEADING_UNCERTAINTY_DEG, estimate_wv_angles


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
        Imaging mode (``"IW"``, ``"EW"``, ``"SM"`` or ``"WV"``) when the
        incidence / elevation angles were estimated, otherwise ``None``.  For
        ``"WV"`` the heading is the WV track model of :mod:`s1_geometry.wv`
        (not the IW latitude polynomial) and its uncertainty is
        :data:`~s1_geometry.wv.WV_HEADING_UNCERTAINTY_DEG` (0.25 deg, the
        p95 error of the worst product of the S1C validation).
    beam:
        SM beam (``"S1".."S6"``) when :attr:`mode` is ``"SM"``, otherwise
        ``None``.  The SM chain is spec-based and not validated against ESA
        L2 data (no SM L2 products in the archive).
    incidence_deg:
        Incidence angle (deg, from the local vertical) at the reference
        point, or ``None`` when :attr:`mode` is ``None``.
    elevation_deg:
        Radar elevation angle at the reference point (the ESA L2 OCN
        ``*ElevationAngle`` convention: the angle at the satellite between
        nadir and the line of sight, equal to the off-nadir angle), or
        ``None`` when :attr:`mode` is ``None``.
    subswath:
        WV subswath (1 or 2) of the reference point when :attr:`mode` is
        ``"WV"``, otherwise ``None``.
    r_s_km:
        Satellite geocentric distance (km) used for the angles (IW / EW: the
        footprint-width inversion; WV: the given or orbit-model value), or
        ``None`` when :attr:`mode` is ``None``.
    """

    orbit_direction: str
    ref_lat: float
    heading_deg: float
    heading_uncertainty_deg: float
    mode: str | None = None
    beam: str | None = None
    incidence_deg: float | None = None
    elevation_deg: float | None = None
    subswath: int | None = None
    r_s_km: float | None = None


def _estimate_wv_geometry(sar_polygon, ref_lon, ref_lat, r_s_km, orbit_direction):
    """The WV branch of :func:`estimate_geometry` (one reference point)."""
    if ref_lon is None or ref_lat is None:
        raise ValueError("mode='WV' needs ref_lon and ref_lat (a point inside an imagette)")
    wv = estimate_wv_angles(sar_polygon, r_s_km, ref_lon, ref_lat,
                            orbit_direction=orbit_direction)
    heading = float(wv.heading_deg[0])
    # pass direction of the reference point from its heading (northbound
    # flight = ascending); the WV track model already handles the turn
    orbit = ASCENDING if math.cos(math.radians(heading)) >= 0.0 else DESCENDING
    return S1Geometry(
        orbit_direction=orbit,
        ref_lat=float(ref_lat),
        heading_deg=heading,
        heading_uncertainty_deg=WV_HEADING_UNCERTAINTY_DEG,
        mode="WV",
        incidence_deg=float(wv.incidence_deg[0]),
        elevation_deg=float(wv.elevation_deg[0]),
        subswath=int(wv.subswath[0]),
        r_s_km=float(wv.r_s_km),
    )


def estimate_geometry(
    sar_polygon,
    ref_lat: float | None = None,
    mode: str | None = None,
    ref_lon: float | None = None,
    beam: str | None = None,
    *,
    r_s_km: float | None = None,
    orbit_direction: str | None = None,
) -> S1Geometry:
    """Estimate the geometric quantities of a Sentinel-1 acquisition.

    Parameters
    ----------
    sar_polygon:
        The SAFE footprint: a sequence of (lon, lat) vertices, a GeoJSON
        Polygon string (the CDSE ``GeoFootprint`` column), a GeoJSON
        Polygon dict (the matchup STAC item representation) or a WKT string
        (the CDSE OData ``Footprint``).  For ``mode="WV"`` it is the product
        MultiPolygon (one ring per imagette).
    ref_lat:
        Reference latitude (deg) at which the heading is evaluated. Defaults
        to the mean latitude of the 4 footprint vertices.
    mode:
        Imaging mode, ``"IW"``, ``"EW"``, ``"SM"`` or ``"WV"``.  When given,
        the incidence and elevation angles are estimated at the reference
        point; for IW / EW / SM the heading fields remain the IW model (a
        documented extrapolation for EW / SM), for WV they come from the WV
        track model (see :func:`~s1_geometry.wv.estimate_wv_angles`).
    ref_lon:
        Reference longitude (deg) of the point at which the incidence /
        elevation angles are evaluated.  Must be given together with
        ``ref_lat``; when both are omitted and ``mode`` is given, the
        footprint center (mean of the 4 vertices) is used.
    beam:
        SM beam (``"S1"`` to ``"S6"``, default ``"S3"``), used only when
        ``mode="SM"``.  The SM chain is spec-based (documented off-nadir
        bands) and not validated against ESA L2 data.
    r_s_km, orbit_direction:
        ``mode="WV"`` only: the satellite geocentric distance (km) and the
        pass direction; both default to the metadata-only orbit model (see
        :func:`~s1_geometry.orbit.satellite_radius_km`).

    Raises
    ------
    ValueError
        If exactly one of ``ref_lon`` / ``ref_lat`` is given, the mode is
        unknown, or the footprint width does not invert to an orbit altitude
        in the documented band (see
        :func:`~s1_geometry.angles.estimate_incidence_elevation`).
    """
    if mode is not None and str(mode).upper() == "WV":
        return _estimate_wv_geometry(sar_polygon, ref_lon, ref_lat, r_s_km, orbit_direction)
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
    angles = estimate_incidence_elevation(
        sar_polygon, ref_lon, ref_lat, mode, beam)
    return S1Geometry(
        orbit_direction=orbit,
        ref_lat=float(ref_lat),
        heading_deg=geometry.heading_deg,
        heading_uncertainty_deg=geometry.heading_uncertainty_deg,
        mode=angles.mode,
        beam=angles.beam,
        incidence_deg=angles.incidence_deg,
        elevation_deg=angles.elevation_deg,
        r_s_km=angles.r_s_km,
    )
