# s1-geometry

Tool to estimate geometric quantities (such as incidence, elevation or
ground heading angle) from a Sentinel-1 SAR footprint.

The current release (v0.1) provides the **IW ground heading** as a function
of the reference point and the SAFE footprint polygon, with a
latitude-dependent uncertainty. Incidence angle and elevation angle
estimators are planned on the same footprint-based architecture.

## Install

```bash
# from PyPI
python -m pip install s1-geometry

# from a checkout (editable)
python -m pip install -e .

# optional: CDSE daily-Parquet footprint indexer (load_sar_polygons)
python -m pip install -e '.[parquet]'
```

Python >= 3.11 is required; the only hard dependency is `numpy`.

## Usage

### Python

The generic entry point estimates every available geometric quantity from
the SAFE footprint:

```python
from s1_geometry import estimate_geometry

# sar_polygon: the CDSE GeoFootprint — a (lon, lat) vertex sequence,
# a GeoJSON Polygon string, or a GeoJSON Polygon dict (STAC geometry)
result = estimate_geometry(sar_polygon, ref_lat=36.592652)

result.orbit_direction          # "ascending" | "descending" | "turning"
result.ref_lat                  # reference latitude (deg)
result.heading_deg              # ground heading (deg, clockwise from North)
result.heading_uncertainty_deg  # latitude-dependent uncertainty (deg)
```

When `ref_lat` is omitted it defaults to the mean latitude of the footprint
vertices. Catalogue-oriented named functions are also exported (this is the
function behind the SOBA `sar_ground_heading` column):

```python
from s1_geometry import (
    sar_ground_heading,                 # heading (deg) at ref_lat
    sar_ground_heading_uncertainty,     # U(lat) (deg), per hemisphere
    sar_heading_and_uncertainty,        # (heading, U) at once
    orbit_direction_from_ring,          # the dlat rule
)
```

### CLI

```bash
s1-geometry heading footprint.geojson --lat 36.59
s1-geometry heading '{"type": "Polygon", "coordinates": [[...]]}' --lat 36.59 --json
```

The polygon argument is a GeoJSON file path or an inline JSON string: a
GeoJSON Polygon object, a GeoJSON Feature, or a bare ring of [lon, lat]
vertices; `--json` prints a machine-readable result.

A runnable example lives in [`examples/heading_from_footprint.py`](examples/heading_from_footprint.py).

## The model

The heading of an IW acquisition depends almost only on the latitude and the
orbit direction; the four satellites agree within a fraction of a degree, so
two pooled degree-4 polynomials (ascending / descending) describe the whole
constellation. Beyond |lat| = 70 deg the polynomials extrapolate and the
uncertainty model carries the extra error.

The **orbit direction** is derived from the footprint: the CDSE
`GeoFootprint` ring is a 4-vertex polygon in a fixed vertex order
`[bottom-left, top-left, top-right, bottom-right]` (first vertex repeated).
With `dlat = lat(v[1]) - lat(v[0])`:

* `dlat < -0.3 deg` -> ascending (slice oriented NW-SE)
* `dlat > +0.3 deg` -> descending (slice oriented NE-SW)
* otherwise -> **turning** (polar turnaround): near-due-west heading of
  270 deg (the satellite keeps its westward ground drift) with the large
  turning-latitude uncertainty.

The **uncertainty** is per hemisphere: the p95 of the fit error binned by
|lat| (degree-4 for |lat| <= 70, piecewise-linear beyond, clamped at the
measured turning latitudes 84.5 N / 77.5 S). It grows with latitude so that
users can decide for themselves instead of rejecting high-latitude
acquisitions.

## Validation

Fitted on 13 consecutive days of S1 IW OCN OWI data from all four satellites
(11107 products, 44428 grid corners; fit window |lat| <= 70 deg, headings
wrapped to the per-orbit median): ascending n = 4739 residual std 0.30 deg
(p95 0.64), descending n = 5104 residual std 0.32 deg (p95 0.65).

Against the measured OWI heading of the 11107 products: mean |error|
0.749 deg, median 0.139 deg, p95 2.90 deg; **98.88 % of the headings fall
within the returned U(lat)** and 99.96 % within 2 U(lat). In the fitted
window |lat| <= 70 deg the mean error is 0.14-0.25 deg per latitude band and
the p95 stays below 1.19 deg.

The orbit-direction rule was validated on 3844 CDSE GeoFootprint rings with
measured OWI headings (2026-11-01..12, S1C+S1D, both branches): 99.77 %
agreement with the orbit relabel; all exceptions are the |dlat| < 0.3 deg
turning products at 77.1-77.5 deg S, correctly classified as turning.

## Roadmap

- [x] IW ground heading `sar_ground_heading = F(ref_lat, SAR_polygon)` + U(lat)
- [ ] Incidence angle estimator
- [ ] Elevation angle estimator

## Development

```bash
python -m pip install -e '.[parquet]' pytest pre-commit
pre-commit run --all-files   # ruff + file hygiene
python -m pytest -q          # test suite
```

GitHub workflows: `CI` (pre-commit + pytest), `Build` (wheel + CLI smoke
test), `Docs` (sphinx), `Publish` (PyPI on release). Container image:
`apptainer build s1-geometry.sif apptainer.def`.
