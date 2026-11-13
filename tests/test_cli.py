"""Tests for the s1-geometry command-line interface."""

import json

import pytest

from s1_geometry.cli import main
from tests.rings import ASC_RING, DESC_RING


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
