# s1-geometry

[![CI](https://github.com/umr-lops/s1-geometry/actions/workflows/ci.yml/badge.svg)](https://github.com/umr-lops/s1-geometry/actions/workflows/ci.yml)
[![Docs](https://readthedocs.org/projects/s1-geometry/badge/?version=latest)](https://s1-geometry.readthedocs.io)
[![PyPI](https://img.shields.io/pypi/v/s1-geometry.svg)](https://pypi.org/project/s1-geometry/)
[![Python](https://img.shields.io/pypi/pyversions/s1-geometry.svg)](https://pypi.org/project/s1-geometry/)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)

Tool to estimate geometric quantities (such as incidence, elevation or
ground heading angle) from a Sentinel-1 SAR footprint, **from the catalogue
metadata only** (e.g. the CDSE OData `Footprint`): no SAR data, no orbit
file.

It provides:

- the **IW ground heading** as a function of the reference point and the
  SAFE footprint polygon, with a latitude-dependent uncertainty;
- the **incidence / elevation angles** of IW, EW (validated) and SM
  (spec-based) acquisitions at a reference point;
- the **WV (Wave mode) incidence / elevation / heading** at reference points,
  from the product footprint (one ring per imagette) and a metadata-only model
  of the satellite distance.

## Install

```bash
# from PyPI
python -m pip install s1-geometry

# from a checkout (editable)
python -m pip install -e .

# optional: CDSE daily-Parquet footprint indexer (load_sar_polygons)
python -m pip install -e '.[parquet]'
```

Python >= 3.11 is required; the only hard dependency is `numpy` (`pandas` for the
table helper and the Parquet indexer, via the `parquet` extra).

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
    orbit_direction_from_ring,          # pass direction (flight rule)
)
```

The incidence / elevation angles of an IW or EW acquisition are estimated at
a reference point from the SAFE footprint and the imaging mode:

```python
from s1_geometry import estimate_incidence_elevation

# sar_polygon as above; mode "IW" (default) or "EW"
angles = estimate_incidence_elevation(sar_polygon, ref_lon=10.5, ref_lat=36.59, mode="IW")

angles.mode             # "IW" | "EW"
angles.dfrac            # 0 = near range line, 1 = far range line
angles.swath_width_km   # across-track ground width of the footprint (km)
angles.r_s_km           # inverted satellite geocentric distance (km)
angles.incidence_deg    # incidence angle from the local vertical (deg)
angles.elevation_deg    # radar elevation angle = off-nadir angle (deg, ESA OCN convention)
```

The same quantities are available on the generic entry point
(`estimate_geometry(sar_polygon, ref_lon=10.5, ref_lat=36.59, mode="IW")`)
and through the catalogue-oriented functions `sar_incidence_angle(...)` /
`sar_elevation_angle(...)` (the functions behind the `sar_incidence_angle` /
`sar_elevation_angle` columns).

#### WV (Wave mode), from the OData metadata only

The WV product footprint (the CDSE OData `Footprint` WKT, or the
`GeoFootprint` GeoJSON) is a MultiPolygon with one small ring per imagette,
in time order. The satellite geocentric distance comes from a metadata-only
orbit model (`satellite_radius_km(lat, orbit_direction)`, 0.1 km accuracy
against the POEORB state vectors of S1A/S1C/S1D), so nothing else is needed:

```python
from s1_geometry import estimate_geometry, estimate_wv_angles

geo = estimate_geometry(odata_footprint, ref_lon=-76.7, ref_lat=-29.3, mode="WV")
geo.incidence_deg, geo.elevation_deg, geo.heading_deg, geo.subswath, geo.r_s_km

wv = estimate_wv_angles(odata_footprint, None, lons, lats)  # vectorised, per point
```

Add the columns to a table of reference points (one row per point, the
product footprint in a column):

```python
import pandas as pd
from s1_geometry.table import add_wv_angle_columns

df = add_wv_angle_columns(pd.read_parquet("wv.parquet"),
                          footprint="Footprint", lon="ref_lon", lat="ref_lat")
# + sar_incidence_angle, sar_elevation_angle, sar_ground_heading, sar_wv_subswath
```

A **single imagette polygon is not enough** (the two subswath chains of the
product locate the track plane): pass the product footprint (at least 4
imagettes), with the rings **in time order** — the CDSE `Footprint` /
`GeoFootprint`, not the re-ordered `geometry` column of the daily Parquets
(it is detected and rejected). Imagettes split at the antimeridian are
re-joined.

The reference point must lie **inside the SAFE footprint** (the 4-vertex
parallelogram, boundary included).  If it is outside,
`estimate_incidence_elevation` raises a `ValueError` ("outside the SAFE
footprint") instead of returning an extrapolated, meaningless value.

### CLI

```bash
s1-geometry heading footprint.geojson --lat 36.59
s1-geometry heading '{"type": "Polygon", "coordinates": [[...]]}' --lat 36.59 --json
s1-geometry angles footprint.geojson --lon 10.5 --lat 36.59 --mode IW
s1-geometry angles footprint.geojson --lon 10.5 --lat 36.59 --json
s1-geometry wv "geography'SRID=4326;MULTIPOLYGON (((...)), ...)'" --lon -76.7 --lat -29.3
```

The polygon argument is a GeoJSON file path or an inline JSON string: a
GeoJSON Polygon object, a GeoJSON Feature, or a bare ring of [lon, lat]
vertices; `--json` prints a machine-readable result.

Runnable examples: [`examples/heading_from_footprint.py`](examples/heading_from_footprint.py)
and [`examples/incidence_elevation_end2end.py`](examples/incidence_elevation_end2end.py)
(the end-to-end script turns a Parquet of footprints + reference positions
into a Parquet with the `sar_incidence_angle` / `sar_elevation_angle` columns;
a bundled 10-row sample lives in `examples/data/`).

## The model

The heading of an IW acquisition depends almost only on the latitude and the
orbit direction; the four satellites agree within a fraction of a degree, so
two pooled degree-4 polynomials (ascending / descending) describe the whole
constellation. Beyond |lat| = 70 deg the polynomials extrapolate and the
uncertainty model carries the extra error.

The **orbit direction** is the *pass* (flight) direction derived from the
footprint, not the orbit branch: a northbound pass is `ascending`, a
southbound pass `descending`. The flight runs along the burst (along-track)
edges of the 4-vertex SAFE parallelogram, which are the *shorter* opposite-
edge family for IW. The rule is order-robust: it handles both the standard
CDSE `GeoFootprint` vertex order (counter-clockwise) and the clockwise ring
copied from the CDSE product page, and reduces to the historical
`dlat = lat(v[1]) - lat(v[0])` rule for the standard order. The pass is
**turning** (polar turnaround, near-due-west heading of 270 deg) when the
flight's north component is below 0.3 deg in magnitude.

The **uncertainty** is per hemisphere: the p95 of the fit error binned by
|lat| (degree-4 for |lat| <= 70, piecewise-linear beyond, clamped at the
measured turning latitudes 84.5 N / 77.5 S). It grows with latitude so that
users can decide for themselves instead of rejecting high-latitude
acquisitions.

The **incidence / elevation angles** combine two consistent models (verified
against the ESA / SentiWiki documented angles to <= 0.07 deg):

* **Constant slant range** (the CDSE catalogue footprint): the across-track
  ground width W of the 4-vertex ring is a monotonic function of the
  satellite geocentric distance R_s; the catalogue polygon is computed with
  the mean Earth radius (6371 km), so inverting W -> R_s (bisection)
  calibrates the satellite altitude from the footprint alone, with no orbit
  file. In the standard CDSE vertex order (GRDH/GRDM/SLC/OCN) the range lines
  are the v0-v1 (far) and v2-v3 (near) edges, ascending and descending alike;
  rings with a rotated vertex order (e.g. ETA) fall back to the other edge
  pair when the standard width does not invert to a valid orbit altitude.
* **Off-nadir** (the actual radar steering): with R_s, the WGS84 geocentric
  radius R_e of the reference point and the off-nadir angle sigma (the linear
  roll law between the documented min / max-altitude values per sub-swath),
  the geocentric incidence is the closed form
  `theta_c = asin[(R_s/R_e) * sin(sigma)]`; the local (ellipsoid-normal)
  incidence adds the deflection of the vertical (max ~0.19 deg near
  |lat| = 45 deg). The returned elevation is the radar elevation angle of
  the ESA L2 OCN `*ElevationAngle` convention (the angle at the satellite
  between nadir and the line of sight, i.e. the off-nadir angle sigma), not
  `90 deg - incidence`.

The reference point selects sigma through its across-track position
`dfrac` (0 = near range line, 1 = far range line). Full-swath footprints are
the supported case; a partial-swath EW footprint whose width does not invert
to an orbit altitude in the documented 697.6-725.8 km band raises
`ValueError`.

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

The orbit-direction rule was validated against the `orbitDirection`
attribute of 66243 CDSE IW/EW rings (13 daily v2 Parquets,
2026-11-01..13): **IW 99.70 % and EW 92.05 %** (both unchanged by the
order-robust rule), **ETA 97.82 %** (the historical dlat rule gave 55.66 %,
because ETA rings use a rotated vertex order). The remaining EW errors are
the pre-existing high-latitude ones. The flight direction is two-valued from
the ring alone (the near/far and start/end corners are unlabeled, a Z2
symmetry): the rule resolves it for the two documented vertex orders
(standard CDSE and product-page); a fully reversed ring may be
misclassified.

The angle model was validated on 66243 CDSE IW/EW rings: **98.9 % (IW) and
90.8 % (EW) of the products invert an orbit altitude in the documented
697.6-725.8 km band**; the off-nadir model reproduces the documented
per-sub-swath angles to <= 0.07 deg, and the full chain (ring -> W -> R_s ->
incidence) is simulated to <= 0.11 deg (IW) and <= 0.05 deg (EW) maximum
error.

The **WV satellite distance** model (`s1_geometry.orbit`) is fitted on 37444
POEORB state vectors (S1A 2026-06-29, S1C 2025-12-10, S1D 2026-05-11 and
2026-10-29; `scripts/calibrate_orbit_radius.py`, public `s1-orbits` bucket):
residual std 0.10 km, max 0.29 km at the WV1 / WV2 target latitudes, bias
per satellite below 0.04 km (~0.02 deg of incidence). The WV angle chain is
validated, from the footprint metadata only, against the ESA L2 OCN grid of
**164 products of Sentinel-1A (2015-2026, not used in the calibration), -1C
and -1D** (latitudes -80 to +52 deg, up to 143 imagettes, 2.97 million
cells): **RMSE incidence 0.036 deg, elevation 0.041 deg, heading 0.10 deg**
(S1A alone: 0.021 / 0.033 / 0.10 deg), subswath labels 100 % correct
(`scripts/validate_wv_ocn.py`). The look-angle law latitude term is fitted
on 62588 imagettes of other S1C / S1D days (`scripts/calibrate_wv_look_angle.py`).
The heading uncertainty returned for WV is 0.4 deg (p95 of the worst product,
through the polar turn). Long products need the sliding-window track fit
(`window="auto"`, default). The orbit model is unit- and slot-independent
(checked across the 2026 S1C relocation). Not covered: products with fewer
than 4 imagettes, the commissioning phase of a satellite.

## Roadmap

- [x] IW ground heading `sar_ground_heading = F(ref_lat, SAR_polygon)` + U(lat)
- [x] Incidence angle estimator `sar_incidence_angle` (IW + EW, SM spec-based)
- [x] Elevation angle estimator `sar_elevation_angle` (IW + EW, SM spec-based)
- [x] WV incidence / elevation / heading from the OData metadata only
- [x] WV validation on S1A / S1C / S1D (164 products, -80 to +52 deg)
- [x] WV look-angle law recalibrated (northern-hemisphere bias removed)
- [x] WV heading uncertainty (0.4 deg)

## Development

```bash
python -m pip install -e '.[parquet,test]' pre-commit
pre-commit run --all-files   # ruff + file hygiene
python -m pytest -q          # test suite
```

GitHub workflows: `CI` (pre-commit, pytest + coverage on Python 3.11-3.13,
examples), `Build` (wheel + CLI smoke test), `Docs` (sphinx -W), `Publish`
(PyPI on release; the version is the git tag, via `hatch-vcs`). See
[CONTRIBUTING.md](CONTRIBUTING.md) and [CHANGELOG.md](CHANGELOG.md). Container image:
`apptainer build s1-geometry.sif apptainer.def`.
