"""Tests for the generic core: estimate_geometry + the S1Geometry dataclass."""

import dataclasses
import json

import pytest

from s1_geometry import S1Geometry, estimate_geometry
from s1_geometry.footprint import ASCENDING, DESCENDING, TURNING
from s1_geometry.heading import heading_from_lat, sar_ground_heading_uncertainty
from tests.rings import ASC_RING, DESC_RING, TURN_RING


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
