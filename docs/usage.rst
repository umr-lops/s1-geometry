Install and run
===============

Install the package in a Python 3.11+ environment::

   python -m pip install s1-geometry

or, from a checkout, in editable mode::

   python -m pip install -e .

The CDSE daily-Parquet footprint indexer (``load_sar_polygons``) needs the
optional ``parquet`` extra::

   python -m pip install s1-geometry[parquet]

Python API
----------

The generic entry point estimates every available geometric quantity from the
SAFE footprint::

   from s1_geometry import estimate_geometry

   result = estimate_geometry(sar_polygon, ref_lat=36.592652)
   # result.orbit_direction      "ascending" | "descending" | "turning"
   # result.ref_lat              reference latitude (deg)
   # result.heading_deg           ground heading (deg, clockwise from North)
   # result.heading_uncertainty_deg  latitude-dependent uncertainty (deg)

``sar_polygon`` is the 4-vertex SAFE footprint (the CDSE ``GeoFootprint``
column): a sequence of (lon, lat) vertices, a GeoJSON Polygon string, or a
GeoJSON Polygon dict (STAC item geometry). When ``ref_lat`` is omitted it
defaults to the mean latitude of the footprint vertices.

Catalogue-oriented named functions are also exported:

* ``sar_ground_heading(ref_lat, sar_polygon)`` -- heading in degrees; this is
  the function behind the SOBA ``sar_ground_heading`` column.
* ``sar_ground_heading_uncertainty(ref_lat)`` -- latitude- and
  hemisphere-dependent uncertainty in degrees.
* ``sar_heading_and_uncertainty(ref_lat, sar_polygon)`` -- both at once.
* ``orbit_direction_from_ring(sar_polygon)`` -- the orbit-direction rule.

CLI
---

::

   s1-geometry heading footprint.geojson --lat 36.59
   s1-geometry heading '{"type": "Polygon", "coordinates": [[...]]}' --lat 36.59 --json

The polygon argument is a GeoJSON file path or an inline JSON string: a
GeoJSON Polygon object, a GeoJSON Feature, or a bare ring of [lon, lat]
vertices; ``--json`` prints a machine-readable result.

The model
---------

The heading of an IW acquisition depends almost only on the latitude and the
orbit direction; the four satellites agree within a fraction of a degree, so
two pooled degree-4 polynomials (ascending / descending) describe the whole
constellation (fit window |lat| <= 70 deg; beyond it they extrapolate and the
uncertainty model carries the extra error). The orbit direction is derived
from the footprint: with ``dlat = lat(v[1]) - lat(v[0])`` on the 4-vertex
ring (fixed CDSE order), ``dlat < -0.3`` is ascending, ``dlat > +0.3`` is
descending, otherwise the product is at the polar turnaround (**turning**)
and receives a near-due-west heading of 270 deg.

Validation
----------

Against the measured OWI heading of the 11107 products used for the fit:
mean |error| 0.749 deg, median 0.139 deg, p95 2.90 deg; 98.88 % of the
headings fall within the returned U(lat) and 99.96 % within 2 U(lat). In the
fitted window |lat| <= 70 deg the mean error is 0.14-0.25 deg per latitude
band and the p95 stays below 1.19 deg.
