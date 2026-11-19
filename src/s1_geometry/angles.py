"""S1 **IW/EW** SAR incidence and elevation angles from the SAFE footprint.

The radar is right-looking and steered to a quasi-constant off-nadir angle
``sigma`` per sub-swath (the design goal is a quasi-constant incidence /
ground swath over the orbit).  Two consistent models are combined
(verified against the ESA / SentiWiki documented angles to <= 0.07 deg):

* **Constant slant range** (the CDSE catalogue footprint).  The across-track
  ground width ``W`` of the ``GeoFootprint`` is a monotonic function of the
  satellite geocentric distance ``R_s`` (the CDSE catalogue polygon is
  computed with the mean Earth radius 6371 km: the observed IW width 251 to
  274 km over the altitude cycle matches the 6371 km model, 248 to 271 km).
  Inverting ``W -> R_s`` (bisection) calibrates the satellite altitude from
  the footprint alone, no orbit file needed.
* **Off-nadir** (the actual radar steering).  With the satellite distance
  ``R_s``, the geocentric radius ``R_e`` of the reference point and the
  off-nadir angle ``sigma`` (the linear roll law between the documented
  min/max-altitude values per sub-swath), the geocentric incidence is the
  closed form ``theta_c = asin[(R_s/R_e) * sin(sigma)]``.  The local
  (ellipsoid-normal) incidence adds the deflection of the vertical
  (max ~0.19 deg near \\|lat\\| = 45 deg).  The **radar elevation angle**
  (the ESA L2 OCN ``*ElevationAngle`` convention: the angle at the satellite
  between nadir and the line of sight) is exactly the off-nadir angle
  ``sigma`` (equivalently ``incidence - gamma``, ``gamma`` being the
  sub-satellite to target central angle; the triangle angle sum gives
  ``sigma = incidence - gamma``).

The across-track width ``W`` is the spherical distance (6371 km sphere)
between the midpoints of the two range lines of the 4-vertex SAFE ring.  The
range lines are the two across-track edges of the parallelogram; for the
standard CDSE vertex order they are the ``v0-v1`` and ``v2-v3`` edges, and
which of the two is the near range line is resolved from the right-looking
geometry (the near range line is the one on the ground-track side).  Rings
with a rotated vertex order (e.g. ETA) fall back to the other edge pair when
the standard width does not invert to a valid orbit altitude.  Validated on
66243 CDSE IW/EW rings: 98.9 % (IW) and 90.8 % (EW) of the products invert an
altitude in the documented 697.6 to 725.8 km band;
the full chain (ring -> W -> R_s -> incidence) is simulated to <= 0.11 deg
(IW) and <= 0.05 deg (EW) maximum error (SOBA SAR-metadata retrieval,
2026-11-13).

The incidence depends only on ``(lat_ref, lon_ref, sub-swath position, W)``;
it does not depend on the heading (which sets the azimuth, not the
incidence).  Full-swath footprints are the supported case; a partial-swath
EW footprint whose width does not invert to an orbit altitude in the
documented band raises :class:`ValueError`.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from s1_geometry.footprint import (
    _as_point_array,
    orbit_direction_from_ring,
    point_in_polygon,
)
from s1_geometry.heading import heading_from_lat

#: Mean Earth radius (km).  The CDSE catalogue footprint is a simplified
#: constant-slant-range polygon computed with this radius; the width
#: calibration (W -> R_s) uses it too, so that the inverted altitude is
#: ``R_s - R_E_KM`` in the documented 697.6 to 725.8 km band.
R_E_KM = 6371.0

#: Documented S1 orbit altitude band (km) -> satellite geocentric distance.
R_S_MIN_KM = R_E_KM + 697.6
R_S_MAX_KM = R_E_KM + 725.8

#: Reference satellite geocentric distance (km) for the **SM** (stripmap)
#: chain.  The SM ground swath is (by design of the attitude-steering law)
#: nearly independent of the orbit altitude, so the footprint width does not
#: determine ``R_s``; the SM angles are evaluated at this mid-band reference
#: distance, which carries an intrinsic ~ +/-0.8 deg incidence uncertainty
#: over the documented altitude band.
R_S_REF_SM_KM = R_E_KM + 0.5 * (697.6 + 725.8)

#: Bisection search box for the W -> R_s inversion (km): the documented band
#: plus a margin, so that a width slightly outside the band still inverts
#: (the result is then rejected by the band check, which triggers the
#: rotated-ring fallback).
_R_S_BOX = (R_S_MIN_KM - 40.0, R_S_MAX_KM + 40.0)

#: Across-track width separation (km) below which the two candidate range-line
#: families are a "knife-edge" for the orbit-direction rule: the width
#: inversion can no longer tell them apart (both invert, or one is just outside
#: the band), so the flight-parallel (north-south) family is used instead.  For
#: near-square IW footprints the two widths differ by only a few km (~250 vs
#: ~253 km) whereas a real across-track / along-track split (e.g. EW) differs by
#: tens to hundreds of km, so a small absolute margin is safe.
_ORBIT_KNIFE_EDGE_KM = 5.0

#: WGS84 ellipsoid.
WGS84_A = 6378137.0
WGS84_F = 1.0 / 298.257223563
WGS84_E2 = WGS84_F * (2.0 - WGS84_F)

#: Documented off-nadir angles (deg) per mode,
#: ``(near-min, near-max, far-min, far-max)`` at min / max orbit altitude.
#: IW: near = IW1, far = IW3; EW: near = EW1, far = EW5.
MODE_SIGMA = {
    "IW": (27.53, 26.00, 40.40, 38.82),
    "EW": (17.94, 16.36, 41.20, 39.62),
}

#: Documented off-nadir angles (deg) for the **SM** (stripmap) beams, one entry
#: per beam: ``(near-min, near-max, far-min, far-max)`` at min / max orbit
#: altitude (ESA User Handbook, "Stripmap" table).  An SM image is acquired
#: with a *single* beam ``S1..S6``; the imaged swath runs from that beam's
#: near to far off-nadir edge.  Unlike IW / EW these bands are **not**
#: calibrated against ESA L2 OCN products: no SM L2 products exist in the CDSE
#: archive (2014-2026), so the SM chain is spec-based only and *not validated*
#: against L2 data.
SM_BEAM_SIGMA = {
    "S1": (17.93, 16.45, 23.53, 21.96),
    "S2": (21.00, 19.51, 26.33, 24.77),
    "S3": (26.18, 24.67, 30.99, 29.45),
    "S4": (30.87, 29.34, 35.15, 33.63),
    "S5": (35.07, 33.53, 38.85, 37.34),
    "S6": (37.53, 35.98, 41.01, 39.51),
}


def _sigma_constants(
    mode: str, beam: str | None = None,
) -> tuple[float, float, float, float]:
    """The ``(near-min, near-max, far-min, far-max)`` off-nadir band for a mode.

    ``beam`` selects the SM beam (``"S1".."S6"``, default ``"S3"``); it is
    ignored for IW / EW, whose near / far edges are fixed (IW1/IW3, EW1/EW5).
    """
    if mode in MODE_SIGMA:
        return MODE_SIGMA[mode]
    if mode == "SM":
        b = (beam or "S3").upper()
        if b not in SM_BEAM_SIGMA:
            raise ValueError(
                f"unknown SM beam {beam!r} (expected one of S1..S6)")
        return SM_BEAM_SIGMA[b]
    raise ValueError(f"unknown mode: {mode!r} (expected 'IW', 'EW' or 'SM')")


@dataclass(frozen=True)
class S1Angles:
    """Incidence / elevation angles estimated from an S1 SAFE footprint.

    Attributes
    ----------
    mode:
        Imaging mode, ``"IW"``, ``"EW"`` or ``"SM"``.
    beam:
        SM beam (``"S1".."S6"``) when :attr:`mode` is ``"SM"``, otherwise
        ``None``.  The SM chain is spec-based (documented off-nadir bands) and
        not validated against ESA L2 data (no SM L2 products in the archive).
    ref_lon, ref_lat:
        Reference point (deg) at which the angles are evaluated.
    dfrac:
        Across-track position of the reference point, 0 = near range line,
        1 = far range line (of the footprint).
    swath_width_km:
        Across-track ground width of the footprint (km).
    r_s_km:
        Inverted satellite geocentric distance (km); the altitude above the
        mean Earth is ``r_s_km - R_E_KM``.
    incidence_deg:
        Incidence angle at the reference point, measured from the local
        (ellipsoid) vertical, in degrees.
    off_nadir_deg:
        Off-nadir (radar steering) angle at the reference point, in degrees:
        the linear roll law between the documented near / far sub-swath values
        interpolated by :attr:`dfrac`.
    elevation_deg:
        Radar elevation angle, the ESA L2 OCN ``*ElevationAngle`` convention:
        the angle at the satellite between nadir and the line of sight, which
        equals the off-nadir angle :attr:`off_nadir_deg` (equivalently
        ``incidence - gamma``), in degrees.
    """

    mode: str
    ref_lon: float
    ref_lat: float
    dfrac: float
    swath_width_km: float
    r_s_km: float
    incidence_deg: float
    off_nadir_deg: float
    elevation_deg: float
    beam: str | None = None


# --------------------------------------------------------------------------- #
# WGS84 helpers
# --------------------------------------------------------------------------- #
def geocentric_radius_km(lon_deg: float, lat_deg: float) -> float:
    """WGS84 geocentric radius (km) of the point (lon_deg, lat_deg)."""
    lon = math.radians(lon_deg)
    lat = math.radians(lat_deg)
    n = WGS84_A / math.sqrt(1.0 - WGS84_E2 * math.sin(lat) ** 2)
    x = n * math.cos(lat) * math.cos(lon)
    y = n * math.cos(lat) * math.sin(lon)
    z = n * (1.0 - WGS84_E2) * math.sin(lat)
    return math.sqrt(x * x + y * y + z * z) / 1e3


def deflection_of_vertical_deg(lat_deg: float) -> float:
    """Angle (deg) between the geocentric direction and the ellipsoid normal.

    Maximum ~0.19 deg near \\|lat\\| = 45 deg; zero at the equator and poles.
    """
    lat = math.radians(lat_deg)
    lat_gc = math.atan((1.0 - WGS84_E2) * math.tan(lat))
    return math.degrees(lat - lat_gc)


# --------------------------------------------------------------------------- #
# Off-nadir model (the actual radar steering)
# --------------------------------------------------------------------------- #
def sigma_of_r_s(sigma_min: float, sigma_max: float, r_s_km: float) -> float:
    """Off-nadir angle (deg) at satellite distance ``r_s_km`` (km).

    Linear attitude-steering (roll) law between the documented min / max
    orbit-altitude values.
    """
    t = (r_s_km - R_S_MIN_KM) / (R_S_MAX_KM - R_S_MIN_KM)
    return sigma_min + t * (sigma_max - sigma_min)


def incidence_offnadir_deg(r_s_km: float, r_e_km: float, sigma_deg: float):
    """Geocentric incidence (deg) ``= asin[(R_s/R_e) * sin(sigma)]``.

    Returns ``None`` when the argument exceeds 1 (the beam does not reach
    the ground at that geometry, e.g. the EW near edge at high altitude).
    """
    x = (r_s_km / r_e_km) * math.sin(math.radians(sigma_deg))
    if x > 1.0:
        return None
    return math.degrees(math.asin(x))


# --------------------------------------------------------------------------- #
# Constant-slant-range model (the CDSE catalogue footprint)
# --------------------------------------------------------------------------- #
def slant_from_incidence(theta_deg: float, r_s_km: float, r_e_km: float) -> float:
    """Slant range (km) for a geocentric incidence ``theta`` (deg).

    From the law of cosines on the S-C-P triangle:
    ``R_s^2 = R_e^2 + R^2 + 2 R_e R cos(theta)``.
    """
    c = math.cos(math.radians(theta_deg))
    return -r_e_km * c + math.sqrt(r_s_km ** 2 - r_e_km ** 2 * (1.0 - c * c))


def mode_slant_ranges(mode: str, r_s_km: float = R_S_MIN_KM,
                      r_e_km: float = R_E_KM,
                      beam: str | None = None) -> tuple[float, float]:
    """Near / far sub-swath slant ranges (km) for a mode.

    Derived from the documented min-altitude off-nadir angles (near and far
    sub-swath edges); the slant ranges are properties of the pointing and
    of ``r_s_km`` only.  For ``mode="SM"`` the beam (``"S1".."S6"``, default
    ``"S3"``) selects the band.
    """
    sn_min, _, sf_min, _ = _sigma_constants(mode, beam)
    th_n = incidence_offnadir_deg(r_s_km, r_e_km, sn_min)
    th_f = incidence_offnadir_deg(r_s_km, r_e_km, sf_min)
    if th_n is None or th_f is None:
        raise ValueError(f"out-of-domain mode constants for mode {mode!r}")
    return (slant_from_incidence(th_n, r_s_km, r_e_km),
            slant_from_incidence(th_f, r_s_km, r_e_km))


def central_angle_rad(r_s_km: float, r_e_km: float, r_km: float) -> float:
    """Angle at the Earth center (rad) between the sub-satellite point and a
    ground point at slant range ``r_km`` (triangle S-C-P)."""
    c = (r_s_km ** 2 + r_e_km ** 2 - r_km ** 2) / (2.0 * r_s_km * r_e_km)
    return math.acos(min(max(c, -1.0), 1.0))


def width_const_slant_km(r_s_km: float, r_e_km: float, r_near_km: float,
                         r_far_km: float) -> float:
    """Across-track ground width (km) for constant slant ranges."""
    g_near = central_angle_rad(r_s_km, r_e_km, r_near_km) * r_e_km
    g_far = central_angle_rad(r_s_km, r_e_km, r_far_km) * r_e_km
    return g_far - g_near


def solve_rs_from_width(w_km: float, r_e_km: float, r_near_km: float,
                        r_far_km: float) -> float | None:
    """Invert the across-track ground width to the satellite geocentric
    distance ``R_s`` (km).

    ``W(R_s)`` is monotonic (constant-slant-range model), so a bisection on
    the documented altitude band plus a margin is robust.  Returns ``None``
    when the width is outside the reachable range.
    """
    def f(r_s: float) -> float:
        return width_const_slant_km(r_s, r_e_km, r_near_km, r_far_km) - w_km

    lo, hi = _R_S_BOX
    f_lo, f_hi = f(lo), f(hi)
    if f_lo * f_hi > 0:
        return None
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        f_mid = f(mid)
        if abs(f_mid) < 1e-9 or (hi - lo) < 1e-9:
            return mid
        if f_lo * f_mid <= 0:
            hi, f_hi = mid, f_mid
        else:
            lo, f_lo = mid, f_mid
    return 0.5 * (lo + hi)


# --------------------------------------------------------------------------- #
# Footprint width and across-track position
# --------------------------------------------------------------------------- #
def _ring_corners(ring) -> list[tuple[float, float]]:
    """The 4 (lon, lat) corners of the SAFE ring, in ring order."""
    points = _as_point_array(ring)
    return [tuple(vertex) for vertex in points[:4]]


def _spherical_km(p: tuple[float, float], q: tuple[float, float],
                  r_e_km: float = R_E_KM) -> float:
    """Spherical distance (km) between two (lon, lat) points on a sphere of
    radius ``r_e_km``."""
    lat1, lat2 = math.radians(p[1]), math.radians(q[1])
    lon1, lon2 = math.radians(p[0]), math.radians(q[0])
    c = (math.sin(lat1) * math.sin(lat2)
         + math.cos(lat1) * math.cos(lat2) * math.cos(lon1 - lon2))
    return math.acos(min(max(c, -1.0), 1.0)) * r_e_km


def _edge_midpoint(a: tuple[float, float],
                   b: tuple[float, float]) -> tuple[float, float]:
    return ((a[0] + b[0]) / 2.0, (a[1] + b[1]) / 2.0)


def swath_width_km(ring) -> float:
    """Across-track ground width (km) of the SAFE footprint.

    The spherical distance (6371 km sphere, the CDSE catalogue convention)
    between the midpoints of the two range lines.  In the standard CDSE
    vertex order the range lines are the ``v0-v1`` and ``v2-v3`` edges; a
    ring whose standard-order width does not invert to an orbit altitude in
    the documented band is assumed to use a rotated vertex order (e.g. ETA)
    and the other edge pair is used instead.  Which of the two is the near
    range line is resolved from the right-looking geometry (see
    :func:`_near_far_from_look`).
    """
    return _range_lines(ring)["width_km"]


def _range_line_family(ring) -> str:
    """Return the range-line family of the SAFE ring.

    ``"standard"`` when the range lines are the ``v0-v1`` / ``v2-v3`` edges
    (the standard CDSE vertex order) and ``"rotated"`` when they are the
    ``v1-v2`` / ``v3-v0`` edges (a rotated order such as ETA or the
    product-page order).

    The family is decided by width inversion: the range-line family is the
    one whose across-track width inverts to an orbit altitude in the
    documented band (or, for SM, to the SM width band).  This is the same
    discriminator :func:`swath_width_km` uses and it stays correct for
    near-square footprints, where the edge lengths of the two families are
    too close to tell them apart.  Unlike :func:`_range_lines` it does not
    resolve the near / far ordering, so it can be called from
    ``orbit_direction_from_ring`` without recursing back into it.

    Raises
    ------
    ValueError
        If neither edge pair inverts in band (partial-swath footprint or an
        unknown ring convention).
    """
    corners = _ring_corners(ring)
    v0, v1, v2, v3 = corners
    candidates = (
        ("standard", _edge_midpoint(v0, v1), _edge_midpoint(v2, v3)),
        ("rotated", _edge_midpoint(v1, v2), _edge_midpoint(v3, v0)),
    )
    for label, mid_a, mid_b in candidates:
        width = _spherical_km(mid_a, mid_b)
        r_s = _invert_width(width)
        if r_s is not None and R_S_MIN_KM <= r_s <= R_S_MAX_KM:
            return label
    raise ValueError(
        "the SAFE footprint width does not invert to an orbit altitude in "
        "the documented band: no range-line family identified "
        "(partial-swath footprint or unknown ring convention)"
    )


def _orbit_range_line_family(ring) -> str:
    """Decide the range-line family for the orbit-direction rule.

    For a near-square footprint the two edge families are too close in length
    to tell apart, so the family is decided from the two candidate across-track
    widths (the distance between the midpoints of each opposite-edge pair):

    * when the widths are clearly separated (>= :data:`_ORBIT_KNIFE_EDGE_KM`,
      e.g. an EW ring whose across-track swath is far from its along-track
      extent) the width-inversion family :func:`_range_line_family` is used;
    * when they are a knife-edge (< :data:`_ORBIT_KNIFE_EDGE_KM`, e.g. a
      near-square IW ring, ~250 vs ~253 km) the flight-parallel family is the
      one more aligned with the orbit, i.e. with the larger north-south extent
      (the Sentinel-1 orbit runs north-south away from the polar turn).  This
      is what resolves the near-square IW rings whose width is just below the
      inverting band, where the width inversion alone is a knife-edge.

    Raises
    ------
    ValueError
        Like :func:`_range_line_family`, if the widths are clearly separated
        but neither inverts in band (partial-swath footprint or unknown ring
        convention).
    """
    v0, v1, v2, v3 = _ring_corners(ring)
    std_w = _spherical_km(_edge_midpoint(v0, v1), _edge_midpoint(v2, v3))
    rot_w = _spherical_km(_edge_midpoint(v1, v2), _edge_midpoint(v3, v0))
    if abs(std_w - rot_w) >= _ORBIT_KNIFE_EDGE_KM:
        return _range_line_family(ring)
    # Knife-edge: the flight-parallel (range-line) family is the one with the
    # larger north-south extent (the orbit runs N-S away from the polar turn).
    ns_std = abs(v1[1] - v0[1]) + abs(v3[1] - v2[1])
    ns_rot = abs(v2[1] - v1[1]) + abs(v0[1] - v3[1])
    return "standard" if ns_std > ns_rot else "rotated"


def _range_lines(ring) -> dict:
    """Identify the range lines of the SAFE ring.

    Returns ``{"width_km", "near_mid", "far_mid", "r_s_km", "label"}`` with
    the across-track width, the (lon, lat) midpoints of the near / far range
    lines, the inverted orbit altitude, and the range-line family label.

    Raises
    ------
    ValueError
        If neither edge pair inverts to an orbit altitude in the documented
        band (partial-swath footprint, unknown ring convention, or a width
        outside the 697.6 to 725.8 km altitude reach).
    """
    corners = _ring_corners(ring)
    v0, v1, v2, v3 = corners
    candidates = {
        "standard": (_edge_midpoint(v0, v1), _edge_midpoint(v2, v3)),
        "rotated": (_edge_midpoint(v1, v2), _edge_midpoint(v3, v0)),
    }
    try:
        label = _range_line_family(ring)
    except ValueError:
        # Re-raise with the detailed per-family widths for diagnostics.
        raise ValueError(
            "the SAFE footprint width does not invert to an orbit altitude in "
            f"the documented 697.6 to 725.8 km band (standard width "
            f"{_spherical_km(_edge_midpoint(v0, v1), _edge_midpoint(v2, v3)):.1f} "
            f"km, rotated width "
            f"{_spherical_km(_edge_midpoint(v1, v2), _edge_midpoint(v3, v0)):.1f} "
            "km): partial-swath footprint or unknown ring convention"
        ) from None
    mid_a, mid_b = candidates[label]
    width = _spherical_km(mid_a, mid_b)
    r_s = _invert_width(width)
    # The near / far ordering is the same for every vertex order: the
    # right-looking radar puts the near range line on the ground-track
    # side and the far range line furthest along the look direction.
    near, far = _near_far_from_look(ring, mid_a, mid_b)
    return {"width_km": width, "near_mid": near, "far_mid": far,
            "r_s_km": r_s, "label": label}


def _sm_width_band_km() -> tuple[float, float]:
    """Across-track width (km) range of the SM beams over the altitude band.

    Each SM beam is imaged alone and the attitude-steering law keeps its
    ground swath fixed, so every beam inverts to a nearly altitude-invariant
    width (~80 km).  The band is the min / max over all beams at the min and
    max orbit altitude, with a small margin for model / product variation.
    """
    widths = []
    for beam in SM_BEAM_SIGMA:
        sn_min, sn_max, sf_min, sf_max = SM_BEAM_SIGMA[beam]
        for r_s in (R_S_MIN_KM, R_S_MAX_KM):
            rn = mode_slant_ranges("SM", r_s, R_E_KM, beam)[0]
            rf = mode_slant_ranges("SM", r_s, R_E_KM, beam)[1]
            widths.append(width_const_slant_km(r_s, R_E_KM, rn, rf))
    return min(widths) - 2.0, max(widths) + 2.0


def _invert_width(width_km: float) -> float | None:
    """Invert a candidate width to a satellite geocentric distance (km).

    The mode is not known from the ring alone.  The IW and EW widths invert
    to an orbit altitude by bisection (disjoint bands, [~248, ~271] km and
    [~412, ~486] km).  The SM beams carry a nearly altitude-invariant width
    (~80 km), so the SM width does not invert to an altitude; any width in
    the SM band is accepted and maps to the reference distance
    :data:`R_S_REF_SM_KM` (see :func:`_sm_width_band_km`).
    """
    for mode in MODE_SIGMA:
        r_near, r_far = mode_slant_ranges(mode)
        r_s = solve_rs_from_width(width_km, R_E_KM, r_near, r_far)
        if r_s is not None and R_S_MIN_KM <= r_s <= R_S_MAX_KM:
            return r_s
    sm_lo, sm_hi = _sm_width_band_km()
    if sm_lo <= width_km <= sm_hi:
        return R_S_REF_SM_KM
    return None


def _near_far_from_look(
    ring, mid_a: tuple[float, float], mid_b: tuple[float, float],
) -> tuple[tuple[float, float], tuple[float, float]]:
    """Order two range-line midpoints into (near, far).

    Sentinel-1 is a *right-looking* radar: the footprint extends from the
    ground track along the look direction (right of the flight heading).
    The near range line is the one closest to the track, i.e. with the
    smallest projection on the look direction (the far range line lies
    furthest along it).  Validated against the per-pixel incidence angles
    of ESA OCN products (the near edge carries the lowest incidence).
    """
    mean_lat = float(sum(vertex[1] for vertex in _ring_corners(ring)) / 4.0)
    orbit = orbit_direction_from_ring(ring)
    heading = math.radians(heading_from_lat(mean_lat, orbit))
    look = (math.cos(heading), -math.sin(heading))  # (east, north), right of flight
    proj = [mid[0] * look[0] + mid[1] * look[1] for mid in (mid_a, mid_b)]
    if proj[0] <= proj[1]:
        return mid_a, mid_b
    return mid_b, mid_a


def dfrac_of_point_from_lines(
    lines: dict, ref_lon: float, ref_lat: float,
) -> float:
    """Across-track fraction of a point, given precomputed range lines."""
    near = lines["near_mid"]
    axis = (lines["far_mid"][0] - near[0], lines["far_mid"][1] - near[1])
    norm = axis[0] ** 2 + axis[1] ** 2
    if norm == 0.0:
        raise ValueError("degenerate SAFE footprint (zero width)")
    return ((ref_lon - near[0]) * axis[0] + (ref_lat - near[1]) * axis[1]) / norm


def dfrac_of_point(ring, ref_lon: float, ref_lat: float) -> float:
    """Across-track position (0 = near range line, 1 = far range line) of
    the point ``(ref_lon, ref_lat)`` within the SAFE footprint."""
    return dfrac_of_point_from_lines(_range_lines(ring), ref_lon, ref_lat)


def _ref_point_in_footprint(ring, ref_lon: float, ref_lat: float) -> bool:
    """Whether the point lies inside (or on the boundary of) the SAFE
    footprint.

    The SAFE ring is a simple 4-vertex polygon.  It is only approximately a
    parallelogram in (lon, lat) space (the range lines are arcs of constant
    slant range on the sphere; the deviation from a parallelogram is ~0.1
    deg for IW and ~0.6 deg for an EW ring at 70 deg S), so a barycentric
    fraction test is not sufficient: a point exactly on the footprint
    boundary must be checked against the ring edges themselves.  Boundary
    points are *inside*: the near / far range lines lie on the ring edges
    and are imaged.  Away from the boundary the test is a plain ray cast,
    which is unambiguous.  The check is order-robust (any vertex order).
    """
    corners = _ring_corners(ring)
    x, y = float(ref_lon), float(ref_lat)
    for i in range(4):
        ax, ay = corners[i]
        bx, by = corners[(i + 1) % 4]
        dx, dy = bx - ax, by - ay
        norm2 = dx * dx + dy * dy
        if norm2 == 0.0:
            continue
        u = ((x - ax) * dx + (y - ay) * dy) / norm2
        if -1e-9 <= u <= 1.0 + 1e-9:
            dist = math.hypot(x - ax - u * dx, y - ay - u * dy)
            if dist <= 1e-5:  # ~1 m: the ring coordinates carry ~0.1 m
                return True
    return point_in_polygon(x, y, ring)


# --------------------------------------------------------------------------- #
# Full chain: footprint + reference point -> incidence & elevation
# --------------------------------------------------------------------------- #
def _resolve_mode(mode: str | None) -> str:
    if mode is None:
        return "IW"
    if mode not in ("IW", "EW", "SM"):
        raise ValueError(
            f"unknown mode: {mode!r} (expected 'IW', 'EW' or 'SM')")
    return mode


def estimate_incidence_elevation(
    sar_polygon, ref_lon: float, ref_lat: float, mode: str | None = None,
    beam: str | None = None,
) -> S1Angles:
    """Estimate the incidence / elevation angles of an S1 acquisition.

    Parameters
    ----------
    sar_polygon:
        The SAFE footprint: a sequence of (lon, lat) vertices, a GeoJSON
        Polygon string (the CDSE ``GeoFootprint`` column) or a GeoJSON
        Polygon dict (the STAC item representation).
    ref_lon, ref_lat:
        Reference point (deg) at which the angles are evaluated (typically
        a point inside the footprint).
    mode:
        Imaging mode, ``"IW"`` (default), ``"EW"`` or ``"SM"``.
    beam:
        SM beam (``"S1"`` to ``"S6"``, default ``"S3"``), used only when
        ``mode="SM"``: an SM image is acquired with a single beam and the
        across-track sigma runs from that beam's near to far off-nadir edge.
        Ignored for IW / EW.

        .. note::

           The SM off-nadir bands are the documented (spec) values; the SM
           chain is **not validated against ESA L2 OCN data** (no SM L2
           products exist in the CDSE archive).

    Returns
    -------
    S1Angles
        The incidence / elevation angles and the intermediate quantities
        (see :class:`S1Angles`).

    Raises
    ------
    ValueError
        If the mode or SM beam is unknown, the reference point is outside
        the SAFE footprint, the ring is malformed, or the footprint width
        does not invert to an orbit altitude in the documented band.
    """
    resolved_mode = _resolve_mode(mode)
    if not _ref_point_in_footprint(sar_polygon, float(ref_lon), float(ref_lat)):
        raise ValueError(
            f"the reference point ({float(ref_lon):.4f}, {float(ref_lat):.4f}) "
            f"is outside the SAFE footprint: the incidence / elevation angles "
            f"cannot be computed outside the imaged area"
        )
    lines = _range_lines(sar_polygon)
    dfrac = dfrac_of_point_from_lines(
        lines, float(ref_lon), float(ref_lat))
    # The band check in _range_lines is mode-agnostic; re-invert with the
    # requested mode, so that a width valid for one mode but not the
    # requested one (e.g. an EW ring passed with mode="IW") is an error.
    # SM is the exception: its ground swath is nearly altitude-invariant
    # (the attitude-steering law keeps the beam's ground band fixed), so
    # the width does not invert to an orbit altitude by bisection; the band
    # check in _range_lines already mapped the width to the reference
    # distance R_S_REF_SM_KM (intrinsic ~ +/-0.8 deg altitude ambiguity).
    if resolved_mode == "SM":
        sm_lo, sm_hi = _sm_width_band_km()
        if not (sm_lo <= lines["width_km"] <= sm_hi):
            raise ValueError(
                f"the footprint width {lines['width_km']:.1f} km is not in "
                f"the SM band ({sm_lo:.1f} to {sm_hi:.1f} km): an SM image "
                f"is acquired with a single beam (S1..S6), so its ground "
                f"swath is much narrower than an IW / EW swath"
            )
        r_s = lines["r_s_km"]
    else:
        r_near, r_far = mode_slant_ranges(resolved_mode)
        r_s = solve_rs_from_width(lines["width_km"], R_E_KM, r_near, r_far)
    if r_s is None or not (R_S_MIN_KM <= r_s <= R_S_MAX_KM):
        raise ValueError(
            f"the footprint width {lines['width_km']:.1f} km does not "
            f"invert to an orbit altitude in the 697.6 to 725.8 km band "
            f"for mode {resolved_mode!r}"
            + (f" beam {beam!r}" if resolved_mode == "SM" else "")
        )

    sn_min, sn_max, sf_min, sf_max = _sigma_constants(resolved_mode, beam)
    sigma = sigma_of_r_s(sn_min, sn_max, r_s) + dfrac * (
        sigma_of_r_s(sf_min, sf_max, r_s) - sigma_of_r_s(sn_min, sn_max, r_s)
    )
    r_e_ref = geocentric_radius_km(float(ref_lon), float(ref_lat))
    theta_gc = incidence_offnadir_deg(r_s, r_e_ref, sigma)
    if theta_gc is None:
        raise ValueError(
            f"out-of-domain geometry: the beam at off-nadir {sigma:.2f} deg "
            f"does not reach the ground at R_e = {r_e_ref:.1f} km"
        )
    incidence = theta_gc + deflection_of_vertical_deg(float(ref_lat))
    # The radar elevation angle (ESA OCN convention) is the off-nadir angle
    # sigma: the triangle angle sum gives sigma = incidence - gamma exactly.
    out_beam = (beam or "S3").upper() if resolved_mode == "SM" else None
    return S1Angles(
        mode=resolved_mode,
        ref_lon=float(ref_lon),
        ref_lat=float(ref_lat),
        dfrac=dfrac,
        swath_width_km=lines["width_km"],
        r_s_km=r_s,
        incidence_deg=incidence,
        off_nadir_deg=sigma,
        elevation_deg=sigma,
        beam=out_beam,
    )


def sar_incidence_angle(
    sar_polygon, ref_lon: float, ref_lat: float, mode: str | None = None,
    beam: str | None = None,
) -> float:
    """Incidence angle (deg, from the local vertical) of an S1 acquisition
    at the reference point.

    ``beam`` selects the SM beam when ``mode="SM"``.  This is the function
    behind the ``sar_incidence_angle`` column of the SAR catalogues.  See
    :func:`estimate_incidence_elevation`.
    """
    return estimate_incidence_elevation(
        sar_polygon, ref_lon, ref_lat, mode, beam).incidence_deg


def sar_elevation_angle(
    sar_polygon, ref_lon: float, ref_lat: float, mode: str | None = None,
    beam: str | None = None,
) -> float:
    """Radar elevation angle (deg) of an S1 acquisition at the reference
    point.

    ``beam`` selects the SM beam when ``mode="SM"``.  The ESA L2 OCN
    ``*ElevationAngle`` convention: the angle at the satellite between nadir
    and the line of sight, equal to the off-nadir (steering) angle.  This is
    the function behind the ``sar_elevation_angle`` column of the SAR
    catalogues.  See :func:`estimate_incidence_elevation`.
    """
    return estimate_incidence_elevation(
        sar_polygon, ref_lon, ref_lat, mode, beam).elevation_deg
