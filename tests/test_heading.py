"""Tests for the heading module (F(ref_lat, orbit direction) + uncertainty)."""

import pytest

from s1_geometry.heading import (
    heading_from_lat,
    sar_ground_heading,
    sar_ground_heading_uncertainty,
    sar_heading_and_uncertainty,
)
from tests.rings import ASC_RING, DESC_RING, TURN_RING

ASC = "ascending"
DESC = "descending"
TURN = "turning"


@pytest.mark.parametrize(
    ("lat", "orbit", "expected"),
    [
        (0.0, ASC, 348.0673),
        (20.0, ASC, 349.0534),
        (45.0, ASC, 350.0983),
        (68.0, ASC, 347.4889),
        (0.0, DESC, 191.9727),
        (20.0, DESC, 190.9730),
        (45.0, DESC, 189.9103),
        (68.0, DESC, 192.4448),
        (45.0, TURN, 270.0),
        (-77.4, TURN, 270.0),
    ],
)
def test_heading_from_lat_anchor_values(lat, orbit, expected):
    assert heading_from_lat(lat, orbit) == pytest.approx(expected, abs=1e-3)


def test_heading_from_lat_is_wrapped_and_extrapolates_beyond_70():
    assert 0.0 <= heading_from_lat(80.0, ASC) < 360.0
    assert 0.0 <= heading_from_lat(-85.0, DESC) < 360.0


def test_heading_from_lat_rejects_unknown_orbit():
    with pytest.raises(ValueError, match="orbit direction"):
        heading_from_lat(45.0, "sideways")


def test_sar_ground_heading_uses_the_polygon_orbit_direction():
    assert sar_ground_heading(36.592652, ASC_RING) == pytest.approx(
        heading_from_lat(36.592652, ASC), abs=1e-9)
    assert sar_ground_heading(36.592652, DESC_RING) == pytest.approx(
        heading_from_lat(36.592652, DESC), abs=1e-9)
    assert sar_ground_heading(-77.42, TURN_RING) == 270.0


@pytest.mark.parametrize(
    ("lat", "expected"),
    [
        (0.0, 0.185),
        (30.0, 1.0539),
        (60.0, 2.4761),
        (-30.0, 1.2038),
        (-60.0, 3.9187),
        (-68.0, 5.5336),
    ],
)
def test_sar_ground_heading_uncertainty_range1(lat, expected):
    assert sar_ground_heading_uncertainty(lat) == pytest.approx(expected, abs=1e-3)


def test_sar_ground_heading_uncertainty_grows_with_latitude():
    for lat in range(70):
        low = sar_ground_heading_uncertainty(float(lat))
        high = sar_ground_heading_uncertainty(float(lat + 1))
        assert low < high or low == pytest.approx(high, abs=1e-6)


def test_sar_ground_heading_uncertainty_high_latitude_and_clamps():
    # North: increasing through the turning latitude, then clamped at 84.5 N.
    assert sar_ground_heading_uncertainty(75.0) < sar_ground_heading_uncertainty(80.0)
    assert sar_ground_heading_uncertainty(84.5) == pytest.approx(35.6918, abs=1e-3)
    assert sar_ground_heading_uncertainty(85.0) == sar_ground_heading_uncertainty(84.5)
    # South: clamped at the turning latitude 77.5 S.
    assert sar_ground_heading_uncertainty(-77.5) == pytest.approx(21.5688, abs=1e-3)
    assert sar_ground_heading_uncertainty(-80.0) == sar_ground_heading_uncertainty(-77.5)


def test_sar_heading_and_uncertainty_uses_reference_latitude_and_polygon():
    heading, uncertainty = sar_heading_and_uncertainty(36.592652, ASC_RING)
    assert heading == pytest.approx(heading_from_lat(36.592652, ASC), abs=1e-9)
    assert uncertainty == pytest.approx(
        sar_ground_heading_uncertainty(36.592652), abs=1e-9)
