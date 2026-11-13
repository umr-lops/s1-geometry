# s1-geometry

## Project
Library to estimate geometric quantities (ground heading today; incidence
and elevation angles planned) from a Sentinel-1 SAR footprint (the SAFE
polygon, e.g. the CDSE `GeoFootprint` column).

## Conventions
- Python >= 3.11, hatchling + `src/s1_geometry` layout, MIT (LOPS).
- Ruff line-length 100 (`select = E, F, W`); pre-commit before commits.
- Tests: `python -m pytest -q` (dev-agent env on the HPC sandbox).
- Public API is re-exported from `s1_geometry/__init__.py`; the catalogue
  column name `sar_ground_heading` must stay a public function name.
- Coefficients live in `src/s1_geometry/heading.py` (embedded, validated on
  11107 IW OCN products); do not change them without refitting and
  revalidating (see README "The model" and "Validation").
- Future estimators (incidence, elevation) plug into `estimate_geometry`
  in `geometry.py`; shared footprint handling belongs in `footprint.py`.

## Environment
- Use the dev-agent interpreter:
  `/scratch/agrouaze/opencode_sandbox/micromamba/envs/dev-agent/bin/python`.
