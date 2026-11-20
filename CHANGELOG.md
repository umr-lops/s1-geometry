# Changelog

All notable changes to this project are documented here. Versions are the
git tags (calendar versioning, `YYYY.MM.DD`).

## Unreleased

### Added
- `s1_geometry.orbit.satellite_radius_km(lat, orbit_direction)`: the satellite
  geocentric distance from metadata only (fitted on POEORB state vectors of
  S1A / S1C / S1D, 0.10 km std), and `scripts/calibrate_orbit_radius.py`.
- WV angles from the CDSE OData metadata only: `estimate_wv_angles(..., r_s_km=None)`
  uses the orbit model; per-ring pass direction (products crossing the polar
  turn); optional sliding-window track fit for long products (`window=`);
  `on_outside="raise" | "nan"`.
- `estimate_geometry(..., mode="WV")` and the `s1-geometry wv` CLI subcommand.
- `s1_geometry.table.add_wv_angle_columns`: add the WV angle columns to a
  pandas table of reference points.
- `parse_wkt`: WKT / EWKT (`geography'SRID=4326;MULTIPOLYGON (...)'`, the OData
  `Footprint`) accepted wherever a polygon is.

### Fixed
- WV: the accepted satellite-distance band was the WGS84 *altitude* band
  added to 6371 km (7068.6-7096.8 km) whereas the real geocentric distance is
  7065.3-7082.4 km: 22 % of the orbit (half of |lat| > 50 deg) was rejected.
- WV: target points outside the footprint are rejected (they silently gave
  meaningless angles, e.g. 95.8 deg incidence).
- `orbit_direction_from_ring` (and `sar_ground_heading`) raise on rings too
  small to be a swath footprint (a WV imagette was classified `turning`,
  heading 270 deg).
- A MultiPolygon passed to the polygon functions now points to the WV API.
- CLI: long inline polygons crashed with "File name too long".
- Docs: the elevation is the radar elevation angle (off-nadir, ESA OCN
  convention), not `90 - incidence`; WV is the Wave mode, not "wide swath".
- Single source of truth for the version (git tag via `hatch-vcs`; the
  package said 0.3.0, the wheel 0.2.0).

### Changed
- WV per-point computation is vectorised.
- CI: Python 3.11-3.13 matrix, coverage, examples run, concurrency;
  Dependabot; ruff rules `I, UP, B, SIM, NPY`.

## 2026.11.18
- WV and SM angle estimation, ESA L2 OCN certification, near-square orbit
  direction fix.

## 2026.11.16
- IW ground heading, IW / EW incidence and elevation angles.
