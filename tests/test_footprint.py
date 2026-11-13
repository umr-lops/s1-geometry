"""Tests for the footprint module (SAFE polygon parsing, orbit direction)."""

import json

import numpy as np
import pytest

from s1_geometry.footprint import (
    ASCENDING,
    DESCENDING,
    TURNING,
    load_sar_polygons,
    orbit_direction_from_ring,
    point_in_polygon,
    polygon_index_for_point,
    ring_to_geojson,
)
from tests.rings import ASC_RING, DESC_RING, TURN_RING


@pytest.mark.parametrize(
    ("ring", "expected"),
    [
        (ASC_RING, ASCENDING),
        (DESC_RING, DESCENDING),
        (TURN_RING, TURNING),
        (ASC_RING[:4], ASCENDING),  # open ring (no repeated first vertex)
    ],
)
def test_orbit_direction_from_ring(ring, expected):
    assert orbit_direction_from_ring(ring) == expected


def test_orbit_direction_from_geojson_string_and_dict():
    geo = json.dumps({"type": "Polygon", "coordinates": [ASC_RING]})
    assert orbit_direction_from_ring(geo) == ASCENDING
    assert orbit_direction_from_ring(json.loads(geo)) == ASCENDING  # STAC dict form


@pytest.mark.parametrize(
    "ring",
    [
        [(0.0, 0.0), (1.0, 1.0), (2.0, 2.0)],                 # < 4 vertices
        [(0.0, 0.0), (1.0, 1.0), (2.0, 2.0), (3.0, 91.0)],    # latitude > 90
        [(0.0, 0.0), (1.0, np.nan), (2.0, 2.0), (3.0, 3.0)],  # non-finite
        "not a polygon",
        {"type": "MultiPolygon", "coordinates": []},
    ],
)
def test_orbit_direction_rejects_invalid_ring(ring):
    with pytest.raises(ValueError):
        orbit_direction_from_ring(ring)


def test_point_in_polygon():
    inside = (-74.812654, 36.592652)  # centroid of the ascending ring
    assert point_in_polygon(*inside, ASC_RING)
    assert not point_in_polygon(-74.812654, 38.6, ASC_RING)  # ~1 deg north
    assert not point_in_polygon(-78.0, 36.592652, ASC_RING)  # ~3 deg west
    # open ring (no repeated first vertex) behaves the same
    assert point_in_polygon(*inside, ASC_RING[:4])
    # GeoJSON string and dict forms
    geo = json.dumps({"type": "Polygon", "coordinates": [ASC_RING]})
    assert point_in_polygon(*inside, geo)
    assert point_in_polygon(*inside, json.loads(geo))


def test_polygon_index_for_point():
    other = [
        (-10.0, 10.0), (-9.0, 11.0), (-7.0, 10.5), (-8.0, 9.5),
        (-10.0, 10.0),
    ]
    assert polygon_index_for_point(-74.812654, 36.592652, [other, ASC_RING]) == 1
    assert polygon_index_for_point(-8.5, 10.25, [other, ASC_RING]) == 0
    assert polygon_index_for_point(0.0, 0.0, [other, ASC_RING]) is None


def test_ring_to_geojson_closes_the_ring_compactly_once():
    geo = json.loads(ring_to_geojson(ASC_RING))
    ring = geo["coordinates"][0]
    assert geo["type"] == "Polygon"
    assert len(ring) == 5
    assert ring[0] == ring[-1]
    assert len(ring_to_geojson(ASC_RING[:4])) > 0  # open ring is closed too


def _footprint_parquet(path, rows):
    pyarrow = pytest.importorskip("pyarrow")
    parquet = pytest.importorskip("pyarrow.parquet")
    table = pyarrow.table(
        {
            "Name": pyarrow.array([name for name, _ in rows], type=pyarrow.string()),
            "GeoFootprint": pyarrow.array([geo for _, geo in rows], type=pyarrow.string()),
        }
    )
    parquet.write_table(table, path)


def test_load_sar_polygons_indexes_only_well_formed_iw_products(tmp_path):
    path = tmp_path / "s1_day.parquet"
    _footprint_parquet(path, [
        ("S1D_IW_OCN__2SDV_20261109T225729_005397_00A4B6_6EA4.SAFE",
         json.dumps({"type": "Polygon", "coordinates": [ASC_RING]})),
        ("S1C_IW_OCN__2SSH_20261109T043709_010207_0145F2_C0E0.SAFE",
         json.dumps({"type": "Polygon", "coordinates": [DESC_RING]})),
        # IW product with a malformed GeoFootprint: skipped
        ("S1B_IW_OCN__2SDV_20261109T010000_002000_002001_B000.SAFE", "not-a-polygon"),
        # IW product with a ring that is not closed: skipped
        ("S1A_IW_OCN__2SDV_20261109T000000_001000_001001_A000.SAFE",
         json.dumps({"type": "Polygon", "coordinates": [ASC_RING[:4]]})),
    ])

    polygons = load_sar_polygons([path])

    assert set(polygons) == {
        "S1D_IW_OCN__2SDV_20261109T225729_005397_00A4B6_6EA4.SAFE",
        "S1C_IW_OCN__2SSH_20261109T043709_010207_0145F2_C0E0.SAFE",
    }
    ring, geo = polygons["S1C_IW_OCN__2SSH_20261109T043709_010207_0145F2_C0E0.SAFE"]
    assert orbit_direction_from_ring(ring) == DESCENDING
    assert json.loads(geo)["type"] == "Polygon"


def test_load_sar_polygons_skips_wv_rows(tmp_path):
    path = tmp_path / "s1_day.parquet"
    _footprint_parquet(path, [
        ("S1D_WV_OCN__2SDV_20261109T225729_005397_00A4B6_6EA4.SAFE",
         json.dumps({"type": "Polygon", "coordinates": [DESC_RING]})),
    ])
    assert load_sar_polygons([path]) == {}
