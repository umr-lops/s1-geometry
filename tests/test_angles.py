"""Tests for the incidence / elevation angle estimation (s1_geometry.angles)."""

import json

import pytest

from s1_geometry import (
    R_E_KM,
    R_S_MAX_KM,
    R_S_MIN_KM,
    S1Angles,
    dfrac_of_point,
    estimate_incidence_elevation,
    sar_elevation_angle,
    sar_incidence_angle,
    swath_width_km,
)
from s1_geometry.angles import (
    MODE_SIGMA,
    deflection_of_vertical_deg,
    dfrac_of_point_from_lines,
    geocentric_radius_km,
    incidence_offnadir_deg,
    mode_slant_ranges,
    solve_rs_from_width,
    sigma_of_r_s,
    width_const_slant_km,
)
from tests.rings import ASC_RING, DESC_RING, EW_RING


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
    for ring in (ASC_RING, DESC_RING):
        assert 240.0 <= swath_width_km(ring) <= 280.0


def test_swath_width_ew_ring():
    assert swath_width_km(EW_RING) == pytest.approx(414.2, abs=1.0)


def test_swath_width_rejects_out_of_band_ring():
    small_box = [(0.0, 0.0), (0.5, 0.0), (0.5, 0.5), (0.0, 0.5), (0.0, 0.0)]
    with pytest.raises(ValueError, match="does not invert"):
        swath_width_km(small_box)


def test_dfrac_of_point_from_lines_on_a_box():
    lines = {"near_mid": (0.5, 1.0), "far_mid": (0.5, 0.0)}
    assert dfrac_of_point_from_lines(lines, 0.5, 1.0) == pytest.approx(0.0, abs=1e-12)
    assert dfrac_of_point_from_lines(lines, 0.5, 0.5) == pytest.approx(0.5, abs=1e-12)
    assert dfrac_of_point_from_lines(lines, 0.5, 0.0) == pytest.approx(1.0, abs=1e-12)
    assert dfrac_of_point_from_lines(lines, 0.5, 0.1) == pytest.approx(0.9, abs=1e-12)


def test_dfrac_at_range_line_midpoints():
    # standard CDSE order: v0-v1 far, v2-v3 near.
    mid_near = _midpoint(ASC_RING, 2, 3)
    mid_far = _midpoint(ASC_RING, 0, 1)
    assert dfrac_of_point(ASC_RING, *mid_near) == pytest.approx(0.0, abs=1e-9)
    assert dfrac_of_point(ASC_RING, *mid_far) == pytest.approx(1.0, abs=1e-9)


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
    assert result.elevation_deg == pytest.approx(90.0 - result.incidence_deg, abs=1e-12)


def test_estimate_default_mode_is_iw():
    clon, clat = _center(ASC_RING)
    with_mode = estimate_incidence_elevation(ASC_RING, clon, clat, "IW")
    default = estimate_incidence_elevation(ASC_RING, clon, clat)
    assert default.incidence_deg == pytest.approx(with_mode.incidence_deg, abs=1e-12)


def test_estimate_iw_incidence_increases_from_near_to_far():
    mid_near = _midpoint(ASC_RING, 2, 3)
    mid_far = _midpoint(ASC_RING, 0, 1)
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
    assert result.elevation_deg == pytest.approx(90.0 - result.incidence_deg, abs=1e-12)
    # elevation / incidence agree on the two entry points
    assert sar_elevation_angle(EW_RING, clon, clat, "EW") == pytest.approx(
        90.0 - result.incidence_deg, abs=1e-12)


def test_estimate_descending_iw_ring():
    clon, clat = _center(DESC_RING)
    result = estimate_incidence_elevation(DESC_RING, clon, clat, "IW")
    assert 25.0 <= result.incidence_deg <= 45.0
    assert R_S_MIN_KM <= result.r_s_km <= R_S_MAX_KM


def test_estimate_accepts_a_geojson_string():
    clon, clat = _center(ASC_RING)
    geo = json.dumps({"type": "Polygon", "coordinates": [ASC_RING]})
    result = estimate_incidence_elevation(geo, clon, clat, "IW")
    assert 25.0 <= result.incidence_deg <= 45.0


# --------------------------------------------------------------------------- #
# Error cases
# --------------------------------------------------------------------------- #
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
