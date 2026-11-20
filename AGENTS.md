# s1-geometry

## Project
Library to estimate geometric quantities (ground heading, incidence and
elevation angles, for IW / EW / SM / WV) from the Sentinel-1 SAR footprint
metadata only (the SAFE polygon, e.g. the CDSE OData `Footprint` /
`GeoFootprint`): no SAR data, no orbit file.

## Conventions
- Python >= 3.11, hatchling + `src/s1_geometry` layout, MIT (LOPS).
- Ruff line-length 100 (`select = E, F, W, I, UP, B, SIM, NPY`); pre-commit before commits.
- Version = git tag (`hatch-vcs`, calendar `YYYY.MM.DD`); update `CHANGELOG.md`.
- Tests: `python -m pytest -q` (dev-agent env on the HPC sandbox).
- Public API is re-exported from `s1_geometry/__init__.py`; the catalogue
  column name `sar_ground_heading` must stay a public function name.
- Embedded coefficients (`heading.py`, `angles.py`, `wv.py`, `orbit.py`)
  are validated deliverables; do not change them without refitting and
  revalidating (see README "The model" / "Validation" and
  `scripts/calibrate_orbit_radius.py` for the orbit model).
- Estimators plug into `estimate_geometry` in `geometry.py`; shared footprint
  handling (GeoJSON / WKT parsing) belongs in `footprint.py`.

## Environment
- Any Python >= 3.11 environment with `pip install -e '.[parquet,test,docs]'`.
