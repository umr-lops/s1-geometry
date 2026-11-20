"""Metadata-only WV chain: estimate_geometry(mode="WV"), the table helper, the CLI."""

import json

import numpy as np
import pytest

from s1_geometry import estimate_geometry
from s1_geometry.cli import main
from tests.test_wv import _centres, _synthetic_wv


def _odata_wkt(rings):
    return "geography'SRID=4326;MULTIPOLYGON (" + ", ".join(
        "((" + ", ".join(f"{lo} {la}" for lo, la in r) + "))" for r in rings) + ")'"


@pytest.fixture(scope="module")
def product():
    rings, subswaths, truth = _synthetic_wv(r_s_km=7077.0, n_bursts=10, along_track_km=60.0)
    lons, lats = _centres(rings)
    return rings, subswaths, truth, lons, lats


def test_estimate_geometry_wv_from_odata_footprint(product):
    rings, subswaths, truth, lons, lats = product
    geo = estimate_geometry(_odata_wkt(rings), ref_lat=lats[3], ref_lon=lons[3], mode="WV")
    assert geo.mode == "WV"
    assert geo.subswath == subswaths[3]
    assert geo.orbit_direction == "ascending"  # the synthetic track flies north
    assert geo.heading_uncertainty_deg == pytest.approx(0.25)
    assert 7063.0 < geo.r_s_km < 7085.0
    assert abs(geo.incidence_deg - truth["inc"][3]) < 1.0
    assert abs((geo.heading_deg - truth["hd"][3] + 180.0) % 360.0 - 180.0) < 1.0
    assert geo.elevation_deg < geo.incidence_deg


def test_estimate_geometry_wv_needs_a_reference_point(product):
    rings = product[0]
    with pytest.raises(ValueError, match="ref_lon and ref_lat"):
        estimate_geometry(rings, mode="WV")


def test_add_wv_angle_columns(product):
    pd = pytest.importorskip("pandas")
    from s1_geometry.table import WV_COLUMNS, add_wv_angle_columns

    rings, subswaths, truth, lons, lats = product
    wkt = _odata_wkt(rings)
    df = pd.DataFrame({
        "Footprint": [wkt] * len(rings) + [wkt, "POLYGON ((0 0, 1 0, 1 1, 0 0))"],
        "ref_lon": list(lons) + [lons[0] + 30.0, 0.5],
        "ref_lat": list(lats) + [lats[0], 0.5],
    })
    with pytest.raises(ValueError):
        add_wv_angle_columns(df)
    out = add_wv_angle_columns(df, on_error="nan")
    assert set(WV_COLUMNS) <= set(out.columns)
    assert "sar_incidence_angle" not in df.columns  # input untouched
    n = len(rings)
    np.testing.assert_allclose(out["sar_incidence_angle"][:n], truth["inc"], atol=1.0)
    np.testing.assert_array_equal(out["sar_wv_subswath"][:n], subswaths)
    # outside point and unusable footprint -> NaN
    assert out.iloc[n:][list(WV_COLUMNS)].isna().all().all()


def test_cli_wv(product, capsys):
    rings, subswaths, _, lons, lats = product
    assert main(["wv", _odata_wkt(rings), "--lon", str(lons[2]), "--lat", str(lats[2]),
                 "--json"]) == 0
    result = json.loads(capsys.readouterr().out)
    assert result["mode"] == "WV"
    assert result["subswath"] == subswaths[2]
    assert 15.0 < result["incidence_deg"] < 45.0


def test_cli_wv_rejects_single_imagette(product, capsys):
    ring = product[0][0]
    polygon = json.dumps({"type": "Polygon", "coordinates": [ring]})
    with pytest.raises(SystemExit):
        main(["wv", polygon, "--lon", "0", "--lat", "0"])
    assert "MultiPolygon" in capsys.readouterr().err
