"""Tests for the metadata-only satellite distance model (s1_geometry.orbit)."""

from pathlib import Path

import numpy as np
import pytest

from s1_geometry import R_S_ORBIT_MAX_KM, R_S_ORBIT_MIN_KM, satellite_radius_km

# Subsample of POEORB state vectors (S1A 2026-06-29, S1C 2025-12-10, S1D
# 2026-05-11 and 2026-10-29) written by scripts/calibrate_orbit_radius.py:
# satellite, day, WV1 / WV2 target latitudes, r_s (km), ascending flag.
FIXTURE = Path(__file__).parent / "data" / "poeorb_r_s_sample.csv"


def _fixture():
    data = np.genfromtxt(FIXTURE, delimiter=",", names=True, dtype=None, encoding="utf-8")
    return data


@pytest.mark.parametrize("beam", ["lat_wv1", "lat_wv2"])
def test_model_reproduces_poeorb_state_vectors(beam):
    data = _fixture()
    asc = data["ascending"] == 1
    model = np.where(asc, satellite_radius_km(data[beam], "ascending"),
                     satellite_radius_km(data[beam], "descending"))
    err = model - data["r_s_km"]
    assert np.abs(err).max() < 0.35
    assert err.std() < 0.15
    for sat in np.unique(data["satellite"]):
        assert abs(err[data["satellite"] == sat].mean()) < 0.05


def test_scalar_and_array_inputs():
    value = satellite_radius_km(10.0, "ASCENDING")
    assert isinstance(value, float)
    arr = satellite_radius_km(np.array([-60.0, 0.0, 60.0]), "descending")
    assert arr.shape == (3,)
    assert (arr > R_S_ORBIT_MIN_KM).all() and (arr < R_S_ORBIT_MAX_KM).all()


def test_latitudes_beyond_the_orbit_are_clamped():
    assert satellite_radius_km(90.0, "ascending") == satellite_radius_km(84.91, "ascending")
    assert satellite_radius_km(-90.0, "descending") == satellite_radius_km(-78.65, "descending")


def test_unknown_direction():
    with pytest.raises(ValueError, match="unknown orbit direction"):
        satellite_radius_km(0.0, "turning")
