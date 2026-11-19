"""Tests for the incidence / elevation angle estimation (s1_geometry.angles)."""

import json
import math

import pytest

from s1_geometry import (
    R_E_KM,
    R_S_MAX_KM,
    R_S_MIN_KM,
    S1Angles,
    SM_BEAM_SIGMA,
    dfrac_of_point,
    estimate_incidence_elevation,
    sar_elevation_angle,
    sar_incidence_angle,
    swath_width_km,
)
from s1_geometry.angles import (
    MODE_SIGMA,
    R_S_REF_SM_KM,
    _range_line_family,
    _ref_point_in_footprint,
    central_angle_rad,
    deflection_of_vertical_deg,
    dfrac_of_point_from_lines,
    geocentric_radius_km,
    incidence_offnadir_deg,
    mode_slant_ranges,
    slant_from_incidence,
    solve_rs_from_width,
    sigma_of_r_s,
    width_const_slant_km,
)
from tests.rings import (
    ASC_RING,
    DESC_RING,
    EW_RING,
    F512_NEAR_SQUARE_RING,
    USER_DESC_RING_CW,
)


def _center(ring):
    """Mean of the 4 ring vertices (the footprint center)."""
    corners = ring[:4]
    return sum(p[0] for p in corners) / 4.0, sum(p[1] for p in corners) / 4.0


def _midpoint(ring, i, j):
    return (
        (ring[i][0] + ring[j][0]) / 2.0,
        (ring[i][1] + ring[j][1]) / 2.0,
    )


# --------------------------------------------------------------------------- #
# WGS84 helpers
# --------------------------------------------------------------------------- #
def test_geocentric_radius_km_wgs84():
    assert geocentric_radius_km(0.0, 0.0) == pytest.approx(6378.137, abs=1e-3)
    assert geocentric_radius_km(0.0, 90.0) == pytest.approx(6356.7523, abs=1e-3)
    assert geocentric_radius_km(0.0, 45.0) == pytest.approx(6367.4895, abs=1e-3)


def test_deflection_of_vertical_deg():
    assert deflection_of_vertical_deg(0.0) == pytest.approx(0.0, abs=1e-12)
    assert deflection_of_vertical_deg(45.0) == pytest.approx(0.1924, abs=1e-3)
    assert deflection_of_vertical_deg(-45.0) == pytest.approx(-0.1924, abs=1e-3)


# --------------------------------------------------------------------------- #
# Off-nadir model (the actual radar steering)
# --------------------------------------------------------------------------- #
def test_sigma_of_r_s_is_linear_between_band_edges():
    sn_min, sn_max, sf_min, sf_max = MODE_SIGMA["IW"]
    assert sigma_of_r_s(sn_min, sn_max, R_S_MIN_KM) == pytest.approx(sn_min)
    assert sigma_of_r_s(sn_min, sn_max, R_S_MAX_KM) == pytest.approx(sn_max)
    assert sigma_of_r_s(sf_min, sf_max, R_S_MIN_KM) == pytest.approx(sf_min)
    assert sigma_of_r_s(sf_min, sf_max, R_S_MAX_KM) == pytest.approx(sf_max)


def test_incidence_offnadir_flat_limit():
    # R_s = R_e: the off-nadir angle is the incidence itself.
    assert incidence_offnadir_deg(6371.0, 6371.0, 30.0) == pytest.approx(30.0, abs=1e-9)


def test_incidence_offnadir_out_of_domain_returns_none():
    # The beam at 80 deg off-nadir does not reach a ground point at 6371 km.
    assert incidence_offnadir_deg(8000.0, 6371.0, 80.0) is None


# --------------------------------------------------------------------------- #
# Constant-slant-range model (the CDSE catalogue footprint)
# --------------------------------------------------------------------------- #
def test_mode_slant_ranges_are_ordered():
    for mode in ("IW", "EW"):
        r_near, r_far = mode_slant_ranges(mode)
        assert r_near < r_far


def test_mode_slant_ranges_rejects_unknown_mode():
    with pytest.raises(ValueError, match="unknown mode"):
        mode_slant_ranges("XX")


def test_width_to_r_s_roundtrip():
    r_near, r_far = mode_slant_ranges("IW")
    width = width_const_slant_km(7085.0, R_E_KM, r_near, r_far)
    assert solve_rs_from_width(width, R_E_KM, r_near, r_far) == pytest.approx(
        7085.0, abs=1e-6)


def test_width_to_r_s_rejects_out_of_band_widths():
    r_near, r_far = mode_slant_ranges("IW")
    assert solve_rs_from_width(100.0, R_E_KM, r_near, r_far) is None
    assert solve_rs_from_width(1000.0, R_E_KM, r_near, r_far) is None


# --------------------------------------------------------------------------- #
# Footprint width and across-track position
# --------------------------------------------------------------------------- #
def test_swath_width_iw_rings():
    for ring in (ASC_RING, DESC_RING, USER_DESC_RING_CW):
        assert 240.0 <= swath_width_km(ring) <= 280.0


def test_swath_width_ew_ring():
    assert swath_width_km(EW_RING) == pytest.approx(414.2, abs=1.0)


def test_swath_width_rejects_out_of_band_ring():
    small_box = [(0.0, 0.0), (0.5, 0.0), (0.5, 0.5), (0.0, 0.5), (0.0, 0.0)]
    with pytest.raises(ValueError, match="does not invert"):
        swath_width_km(small_box)


def test_range_line_family():
    # standard CDSE order: the range lines are the {e0, e2} family
    assert _range_line_family(ASC_RING) == "standard"
    assert _range_line_family(DESC_RING) == "standard"
    # product-page (clockwise) order: the range lines are the {e1, e3} family
    assert _range_line_family(USER_DESC_RING_CW) == "rotated"
    # near-square IW ring: the edge-length heuristic assumes standard, but the
    # width-inversion family is rotated (the misclassification of issue #6)
    assert _range_line_family(F512_NEAR_SQUARE_RING) == "rotated"
    # near-square EW ring: both discriminators agree on standard
    assert _range_line_family(EW_RING) == "standard"


def test_range_line_family_rejects_out_of_band_ring():
    small_box = [(0.0, 0.0), (0.5, 0.0), (0.5, 0.5), (0.0, 0.5), (0.0, 0.0)]
    with pytest.raises(ValueError, match="does not invert"):
        _range_line_family(small_box)


def test_dfrac_of_point_from_lines_on_a_box():
    lines = {"near_mid": (0.5, 1.0), "far_mid": (0.5, 0.0)}
    assert dfrac_of_point_from_lines(lines, 0.5, 1.0) == pytest.approx(0.0, abs=1e-12)
    assert dfrac_of_point_from_lines(lines, 0.5, 0.5) == pytest.approx(0.5, abs=1e-12)
    assert dfrac_of_point_from_lines(lines, 0.5, 0.0) == pytest.approx(1.0, abs=1e-12)
    assert dfrac_of_point_from_lines(lines, 0.5, 0.1) == pytest.approx(0.9, abs=1e-12)


def test_dfrac_at_range_line_midpoints():
    # standard CDSE order, right-looking radar: v0-v1 near, v2-v3 far.
    mid_near = _midpoint(ASC_RING, 0, 1)
    mid_far = _midpoint(ASC_RING, 2, 3)
    assert dfrac_of_point(ASC_RING, *mid_near) == pytest.approx(0.0, abs=1e-9)
    assert dfrac_of_point(ASC_RING, *mid_far) == pytest.approx(1.0, abs=1e-9)


def test_ref_point_in_footprint_box():
    box = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0)]
    assert _ref_point_in_footprint(box, 0.3, 0.7)
    # boundary points are inside (the range lines lie on the ring edges)
    assert _ref_point_in_footprint(box, 0.0, 0.5)
    assert _ref_point_in_footprint(box, 1.0, 0.5)
    assert _ref_point_in_footprint(box, 0.5, 0.0)
    assert _ref_point_in_footprint(box, 0.5, 1.0)
    assert not _ref_point_in_footprint(box, 1.2, 0.5)
    assert not _ref_point_in_footprint(box, -0.1, 0.5)
    assert not _ref_point_in_footprint(box, 0.5, -0.2)
    assert not _ref_point_in_footprint(box, 0.5, 1.3)


# --------------------------------------------------------------------------- #
# Full chain: footprint + reference point -> incidence & elevation
# --------------------------------------------------------------------------- #
def test_estimate_iw_full_chain():
    clon, clat = _center(ASC_RING)
    result = estimate_incidence_elevation(ASC_RING, clon, clat, "IW")
    assert isinstance(result, S1Angles)
    assert result.mode == "IW"
    assert result.ref_lon == pytest.approx(clon)
    assert result.ref_lat == pytest.approx(clat)
    assert R_S_MIN_KM <= result.r_s_km <= R_S_MAX_KM
    assert result.swath_width_km == pytest.approx(254.0, abs=1.0)
    assert result.dfrac == pytest.approx(0.5, abs=1e-9)
    assert 25.0 <= result.incidence_deg <= 45.0
    # radar elevation = off-nadir (the ESA OCN convention), not 90 - incidence
    assert result.elevation_deg == pytest.approx(result.off_nadir_deg, abs=1e-12)
    assert result.elevation_deg != pytest.approx(90.0 - result.incidence_deg, abs=0.1)
    assert 28.0 <= result.elevation_deg <= 42.0


def test_estimate_default_mode_is_iw():
    clon, clat = _center(ASC_RING)
    with_mode = estimate_incidence_elevation(ASC_RING, clon, clat, "IW")
    default = estimate_incidence_elevation(ASC_RING, clon, clat)
    assert default.incidence_deg == pytest.approx(with_mode.incidence_deg, abs=1e-12)


def test_estimate_iw_incidence_increases_from_near_to_far():
    mid_near = _midpoint(ASC_RING, 0, 1)
    mid_far = _midpoint(ASC_RING, 2, 3)
    inc_near = sar_incidence_angle(ASC_RING, *mid_near, "IW")
    inc_far = sar_incidence_angle(ASC_RING, *mid_far, "IW")
    assert inc_near < inc_far
    # the far end is the far sub-swath (IW3, ~40 deg), the near end IW1 (~28 deg)
    assert inc_near < 32.0
    assert inc_far > 38.0


def test_estimate_ew_full_chain():
    clon, clat = _center(EW_RING)
    result = estimate_incidence_elevation(EW_RING, clon, clat, "EW")
    assert result.mode == "EW"
    assert R_S_MIN_KM <= result.r_s_km <= R_S_MAX_KM
    assert result.swath_width_km == pytest.approx(414.2, abs=1.0)
    assert result.dfrac == pytest.approx(0.5, abs=1e-9)
    assert 20.0 <= result.incidence_deg <= 45.0
    # radar elevation = off-nadir (the ESA OCN convention)
    assert result.elevation_deg == pytest.approx(result.off_nadir_deg, abs=1e-12)
    assert 20.0 <= result.elevation_deg <= 42.0
    # elevation agrees on the two entry points
    assert sar_elevation_angle(EW_RING, clon, clat, "EW") == pytest.approx(
        result.off_nadir_deg, abs=1e-12)


def test_estimate_descending_iw_ring():
    clon, clat = _center(DESC_RING)
    result = estimate_incidence_elevation(DESC_RING, clon, clat, "IW")
    assert 25.0 <= result.incidence_deg <= 45.0
    assert R_S_MIN_KM <= result.r_s_km <= R_S_MAX_KM


def test_estimate_product_page_cw_ring():
    # The CW ring copied from the CDSE product page (descending) must invert
    # to a valid orbit altitude and a mid-swath across-track position.
    clon, clat = _center(USER_DESC_RING_CW)
    result = estimate_incidence_elevation(USER_DESC_RING_CW, clon, clat, "IW")
    assert R_S_MIN_KM <= result.r_s_km <= R_S_MAX_KM
    assert result.dfrac == pytest.approx(0.5, abs=1e-9)
    assert 25.0 <= result.incidence_deg <= 45.0
    assert result.elevation_deg == pytest.approx(result.off_nadir_deg, abs=1e-12)


def test_estimate_accepts_a_geojson_string():
    clon, clat = _center(ASC_RING)
    geo = json.dumps({"type": "Polygon", "coordinates": [ASC_RING]})
    result = estimate_incidence_elevation(geo, clon, clat, "IW")
    assert 25.0 <= result.incidence_deg <= 45.0


# --------------------------------------------------------------------------- #
# Error cases
# --------------------------------------------------------------------------- #
def test_estimate_rejects_a_point_outside_the_footprint():
    clon, clat = _center(ASC_RING)
    # a point 2 deg north of the center is outside the ~1.8 deg footprint
    with pytest.raises(ValueError, match="outside the SAFE footprint"):
        estimate_incidence_elevation(ASC_RING, clon, clat + 2.0, "IW")
    with pytest.raises(ValueError, match="outside the SAFE footprint"):
        estimate_incidence_elevation(ASC_RING, clon - 3.0, clat, "IW")
    # the product-page CW ring is checked as well
    with pytest.raises(ValueError, match="outside the SAFE footprint"):
        estimate_incidence_elevation(
            USER_DESC_RING_CW, clon, clat, "IW")


def test_estimate_accepts_points_on_the_footprint_boundary():
    # the near / far range lines lie on the ring edges and are imaged:
    # boundary points must be accepted (a ray-casting test is inconsistent
    # for them).
    for i, j in ((0, 1), (1, 2), (2, 3), (3, 0)):
        result = estimate_incidence_elevation(ASC_RING, *_midpoint(ASC_RING, i, j), "IW")
        assert 0.0 <= result.dfrac <= 1.0


def test_estimate_rejects_unknown_mode():
    clon, clat = _center(ASC_RING)
    with pytest.raises(ValueError, match="unknown mode"):
        estimate_incidence_elevation(ASC_RING, clon, clat, "XX")


def test_estimate_rejects_an_invalid_polygon():
    with pytest.raises(ValueError):
        estimate_incidence_elevation("not a polygon", 0.0, 0.0, "IW")


def test_estimate_rejects_out_of_band_width():
    small_box = [(0.0, 0.0), (0.5, 0.0), (0.5, 0.5), (0.0, 0.5), (0.0, 0.0)]
    with pytest.raises(ValueError, match="does not invert"):
        estimate_incidence_elevation(small_box, 0.25, 0.25, "IW")


def test_estimate_rejects_an_ew_ring_passed_as_iw():
    clon, clat = _center(EW_RING)
    with pytest.raises(ValueError, match="does not invert"):
        estimate_incidence_elevation(EW_RING, clon, clat, "IW")


# --------------------------------------------------------------------------- #
# SM (stripmap) chain — spec-based, not validated against ESA L2 data
# --------------------------------------------------------------------------- #
def _sm_ring(beam: str, lon0: float, lat0: float, length_km: float = 300.0):
    """A physically consistent SM SAFE ring for ``beam``.

    Built from the documented (spec) off-nadir band of the beam with the
    constant-slant-range model at a mid-band orbit distance.  The ring is
    only consistent with the model up to the SM width, which is nearly
    altitude-invariant (see the ``R_S_REF_SM_KM`` docstring).
    """
    sn_min, sn_max, sf_min, sf_max = SM_BEAM_SIGMA[beam]
    sn = sigma_of_r_s(sn_min, sn_max, R_S_REF_SM_KM)
    sf = sigma_of_r_s(sf_min, sf_max, R_S_REF_SM_KM)
    rn = slant_from_incidence(
        incidence_offnadir_deg(R_S_REF_SM_KM, R_E_KM, sn), R_S_REF_SM_KM, R_E_KM)
    rf = slant_from_incidence(
        incidence_offnadir_deg(R_S_REF_SM_KM, R_E_KM, sf), R_S_REF_SM_KM, R_E_KM)
    g_near = central_angle_rad(R_S_REF_SM_KM, R_E_KM, rn) * R_E_KM
    g_far = central_angle_rad(R_S_REF_SM_KM, R_E_KM, rf) * R_E_KM
    km_per_deg = 111.19
    cosl = math.cos(math.radians(lat0))
    dl_n, dl_f = g_near / (km_per_deg * cosl), g_far / (km_per_deg * cosl)
    dd = (length_km / 2.0) / km_per_deg
    near_lon, far_lon = lon0 + dl_n, lon0 + dl_f
    # standard CDSE (counter-clockwise) vertex order: the first edge e0 runs
    # down the near range line for ascending (near = west) and up it for
    # descending (near = east), so e0 is the near range line in both cases.
    ring = [
        (near_lon, lat0 + dd),  # NW: north end of the near range line
        (near_lon, lat0 - dd),  # SW: south end of the near range line
        (far_lon, lat0 - dd),   # SE: south end of the far range line
        (far_lon, lat0 + dd),   # NE: north end of the far range line
    ]
    return ring + [ring[0]], g_near, g_far


def test_sm_beam_sigma_spec_values():
    # documented off-nadir bands (deg) per beam:
    # (near-min-alt, near-max-alt, far-min-alt, far-max-alt)
    assert SM_BEAM_SIGMA == {
        "S1": (17.93, 16.45, 23.53, 21.96),
        "S2": (21.00, 19.51, 26.33, 24.77),
        "S3": (26.18, 24.67, 30.99, 29.45),
        "S4": (30.87, 29.34, 35.15, 33.63),
        "S5": (35.07, 33.53, 38.85, 37.34),
        "S6": (37.53, 35.98, 41.01, 39.51),
    }
    for beam, (sn_min, sn_max, sf_min, sf_max) in SM_BEAM_SIGMA.items():
        # the near edge always looks closer than the far edge of the beam
        assert sn_min < sf_min and sn_max < sf_max
        # the whole SM incidence range is documented as 18.3 to 46.8 deg
        assert 17.0 <= sn_min and sf_max <= 48.0
        _ = beam


def test_sm_mode_slant_ranges_ordered():
    for beam in SM_BEAM_SIGMA:
        r_near, r_far = mode_slant_ranges("SM", beam=beam)
        assert r_near < r_far
    # S1 looks closer than S6: the whole S1 band sits inside near ranges
    s1n, s1f = mode_slant_ranges("SM", beam="S1")
    s6n, s6f = mode_slant_ranges("SM", beam="S6")
    assert s1f < s6n
    assert s6n < s6f


def test_sm_ring_width_is_accepted_by_the_band_check():
    ring, g_near, g_far = _sm_ring("S3", -10.0, 40.0)
    width = swath_width_km(ring)
    assert (g_far - g_near) - 2.0 <= width <= (g_far - g_near) + 2.0
    assert 75.0 <= width <= 85.0


def test_sm_full_chain():
    ring, _, _ = _sm_ring("S3", -10.0, 40.0)
    clon, clat = _center(ring)
    result = estimate_incidence_elevation(ring, clon, clat, "SM", beam="S3")
    assert result.mode == "SM"
    assert result.beam == "S3"
    # the SM width is altitude-invariant: r_s is pinned to the reference
    assert result.r_s_km == pytest.approx(R_S_REF_SM_KM, abs=1e-9)
    assert result.dfrac == pytest.approx(0.5, abs=1e-9)
    # S3 center at the reference distance: sigma ~27.8 deg, inc ~31.5 deg
    assert result.off_nadir_deg == pytest.approx(27.82, abs=0.1)
    assert result.incidence_deg == pytest.approx(31.45, abs=0.1)
    # radar elevation = off-nadir (the ESA OCN convention)
    assert result.elevation_deg == pytest.approx(result.off_nadir_deg, abs=1e-12)
    # the beam parameter is ignored (echoed as None) for IW / EW
    iw = estimate_incidence_elevation(ASC_RING, *_center(ASC_RING), "IW")
    assert iw.beam is None


def test_sm_incidence_increases_from_near_to_far():
    ring, g_near, g_far = _sm_ring("S3", -10.0, 40.0)
    cosl = math.cos(math.radians(40.0))
    lat0 = 40.0
    lon_n = -10.0 + g_near / (111.19 * cosl)
    lon_f = -10.0 + g_far / (111.19 * cosl)
    inc_near = sar_incidence_angle(ring, lon_n, lat0, "SM", beam="S3")
    inc_far = sar_incidence_angle(ring, lon_f, lat0, "SM", beam="S3")
    assert inc_near < inc_far
    # the near / far edges carry the documented S3 band (~29 / ~34 deg)
    assert 27.0 <= inc_near <= 31.0
    assert 32.5 <= inc_far <= 36.0


def test_sm_default_beam_is_s3():
    ring, _, _ = _sm_ring("S3", -10.0, 40.0)
    clon, clat = _center(ring)
    default = estimate_incidence_elevation(ring, clon, clat, "SM")
    explicit = estimate_incidence_elevation(ring, clon, clat, "SM", beam="S3")
    assert default.beam == "S3"
    assert default.incidence_deg == pytest.approx(
        explicit.incidence_deg, abs=1e-12)


def test_sm_rejects_an_unknown_beam():
    ring, _, _ = _sm_ring("S3", -10.0, 40.0)
    clon, clat = _center(ring)
    with pytest.raises(ValueError, match="unknown SM beam"):
        estimate_incidence_elevation(ring, clon, clat, "SM", beam="S9")
    with pytest.raises(ValueError, match="unknown SM beam"):
        mode_slant_ranges("SM", beam="S7")


def test_sm_rejects_an_iw_ring_passed_as_sm():
    # an IW swath is far wider than any SM beam's ground band
    clon, clat = _center(ASC_RING)
    with pytest.raises(ValueError, match="not in the SM band"):
        estimate_incidence_elevation(ASC_RING, clon, clat, "SM")


def test_sm_rejects_an_sm_ring_passed_as_iw():
    ring, _, _ = _sm_ring("S3", -10.0, 40.0)
    clon, clat = _center(ring)
    with pytest.raises(ValueError, match="does not invert"):
        estimate_incidence_elevation(ring, clon, clat, "IW")
