"""Tests for the generic core: estimate_geometry + the S1Geometry dataclass."""

import dataclasses
import json

import pytest

from s1_geometry import S1Geometry, estimate_geometry, estimate_incidence_elevation
from s1_geometry.footprint import ASCENDING, DESCENDING, TURNING
from s1_geometry.heading import heading_from_lat, sar_ground_heading_uncertainty
from tests.rings import ASC_RING, DESC_RING, EW_RING, TURN_RING


def test_estimate_geometry_explicit_reference_latitude():
    result = estimate_geometry(ASC_RING, ref_lat=36.592652)
    assert result.orbit_direction == ASCENDING
    assert result.ref_lat == pytest.approx(36.592652)
    assert result.heading_deg == pytest.approx(
        heading_from_lat(36.592652, ASCENDING), abs=1e-9)
    assert result.heading_uncertainty_deg == pytest.approx(
        sar_ground_heading_uncertainty(36.592652), abs=1e-9)


@pytest.mark.parametrize(
    ("ring", "orbit"),
    [
        (ASC_RING, ASCENDING),
        (DESC_RING, DESCENDING),
        (TURN_RING, TURNING),
    ],
)
def test_estimate_geometry_classifies_all_three_orbits(ring, orbit):
    result = estimate_geometry(ring)
    assert result.orbit_direction == orbit
    assert result.heading_deg == pytest.approx(heading_from_lat(result.ref_lat, orbit))


def test_estimate_geometry_default_reference_latitude_is_the_ring_mean():
    result = estimate_geometry(ASC_RING)
    expected = sum(lat for _, lat in ASC_RING[:4]) / 4
    assert result.ref_lat == pytest.approx(expected)


def test_estimate_geometry_accepts_a_geojson_string():
    geo = json.dumps({"type": "Polygon", "coordinates": [DESC_RING]})
    result = estimate_geometry(geo, ref_lat=36.4)
    assert result.orbit_direction == DESCENDING
    assert result.heading_deg == pytest.approx(heading_from_lat(36.4, DESCENDING), abs=1e-9)


def test_estimate_geometry_rejects_an_invalid_polygon():
    with pytest.raises(ValueError):
        estimate_geometry("not a polygon")


def test_s1_geometry_is_a_frozen_dataclass():
    result = estimate_geometry(ASC_RING)
    assert isinstance(result, S1Geometry)
    with pytest.raises(dataclasses.FrozenInstanceError):
        result.heading_deg = 0.0  # type: ignore[misc]


def test_estimate_geometry_without_mode_keeps_angles_none():
    result = estimate_geometry(ASC_RING)
    assert result.mode is None
    assert result.incidence_deg is None
    assert result.elevation_deg is None


def test_estimate_geometry_with_mode_fills_the_angles():
    result = estimate_geometry(ASC_RING, mode="IW")
    assert result.mode == "IW"
    assert result.incidence_deg is not None
    assert result.elevation_deg == pytest.approx(90.0 - result.incidence_deg, abs=1e-12)


def test_estimate_geometry_with_mode_uses_the_ring_center_by_default():
    corners = ASC_RING[:4]
    clon = sum(p[0] for p in corners) / 4.0
    clat = sum(p[1] for p in corners) / 4.0
    result = estimate_geometry(ASC_RING, mode="IW")
    direct = estimate_incidence_elevation(ASC_RING, clon, clat, "IW")
    assert result.ref_lat == pytest.approx(clat)
    assert result.incidence_deg == pytest.approx(direct.incidence_deg, abs=1e-12)


def test_estimate_geometry_with_mode_and_explicit_reference_point():
    result = estimate_geometry(ASC_RING, ref_lon=-74.8, ref_lat=36.6, mode="IW")
    direct = estimate_incidence_elevation(ASC_RING, -74.8, 36.6, "IW")
    assert result.incidence_deg == pytest.approx(direct.incidence_deg, abs=1e-12)
    assert result.elevation_deg == pytest.approx(direct.elevation_deg, abs=1e-12)


def test_estimate_geometry_with_mode_accepts_ew():
    result = estimate_geometry(EW_RING, mode="EW")
    assert result.mode == "EW"
    assert result.incidence_deg is not None


def test_estimate_geometry_requires_lon_and_lat_together():
    with pytest.raises(ValueError, match="together"):
        estimate_geometry(ASC_RING, mode="IW", ref_lon=-74.8)
    with pytest.raises(ValueError, match="together"):
        estimate_geometry(ASC_RING, mode="IW", ref_lat=36.6)
