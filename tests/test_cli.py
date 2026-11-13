"""Tests for the s1-geometry command-line interface."""

import json

import pytest

from s1_geometry.cli import main
from tests.rings import ASC_RING, DESC_RING, EW_RING


def test_heading_inline_geojson(capsys):
    polygon = json.dumps({"type": "Polygon", "coordinates": [ASC_RING]})
    assert main(["heading", polygon, "--lat", "36.592652"]) == 0
    out = capsys.readouterr().out
    assert "orbit_direction: ascending" in out
    assert "heading_deg: 349." in out


def test_heading_inline_geojson_json_output(capsys):
    polygon = json.dumps({"type": "Polygon", "coordinates": [DESC_RING]})
    assert main(["heading", polygon, "--lat", "36.4", "--json"]) == 0
    result = json.loads(capsys.readouterr().out)
    assert result["orbit_direction"] == "descending"
    assert result["ref_lat"] == pytest.approx(36.4)
    assert set(result) == {
        "orbit_direction", "ref_lat", "heading_deg", "heading_uncertainty_deg",
    }


def test_heading_from_a_geojson_file(tmp_path, capsys):
    path = tmp_path / "footprint.geojson"
    path.write_text(json.dumps({"type": "Polygon", "coordinates": [ASC_RING]}))
    assert main(["heading", str(path), "--lat", "36.592652"]) == 0
    assert "orbit_direction: ascending" in capsys.readouterr().out


def test_heading_accepts_a_bare_ring(capsys):
    ring = json.dumps([[lon, lat] for lon, lat in DESC_RING])
    assert main(["heading", ring, "--lat", "36.4"]) == 0
    assert "orbit_direction: descending" in capsys.readouterr().out


def test_heading_rejects_an_unreadable_polygon(capsys):
    with pytest.raises(SystemExit):
        main(["heading", "not-a-polygon-and-not-a-file", "--lat", "36.592652"])


def test_angles_inline_geojson(capsys):
    polygon = json.dumps({"type": "Polygon", "coordinates": [ASC_RING]})
    assert main(["angles", polygon, "--lon", "-74.81", "--lat", "36.59"]) == 0
    out = capsys.readouterr().out
    assert "mode: IW" in out
    assert "incidence_deg:" in out
    assert "elevation_deg:" in out


def test_angles_inline_geojson_json_output(capsys):
    polygon = json.dumps({"type": "Polygon", "coordinates": [DESC_RING]})
    assert main(
        ["angles", polygon, "--lon", "-6.21", "--lat", "36.42", "--mode", "IW", "--json"],
    ) == 0
    result = json.loads(capsys.readouterr().out)
    assert result["mode"] == "IW"
    assert set(result) == {
        "mode", "ref_lon", "ref_lat", "dfrac",
        "swath_width_km", "r_s_km", "incidence_deg", "elevation_deg",
    }
    assert result["elevation_deg"] == pytest.approx(90.0 - result["incidence_deg"])


def test_angles_ew_mode(capsys):
    polygon = json.dumps({"type": "Polygon", "coordinates": [EW_RING]})
    assert main(
        ["angles", polygon, "--lon", "-123.10", "--lat", "-70.05", "--mode", "EW"],
    ) == 0
    out = capsys.readouterr().out
    assert "mode: EW" in out
    assert "swath_width_km: 414." in out


def test_angles_rejects_an_ew_ring_passed_as_iw(capsys):
    polygon = json.dumps({"type": "Polygon", "coordinates": [EW_RING]})
    with pytest.raises(SystemExit):
        main(["angles", polygon, "--lon", "-123.10", "--lat", "-70.05"])


def test_angles_rejects_an_unreadable_polygon(capsys):
    with pytest.raises(SystemExit):
        main(["angles", "not-a-polygon-and-not-a-file", "--lon", "0.0", "--lat", "0.0"])
