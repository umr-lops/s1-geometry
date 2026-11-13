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

The incidence / elevation angles of an IW or EW acquisition are estimated at
a reference point from the SAFE footprint and the imaging mode::

   from s1_geometry import estimate_incidence_elevation

   result = estimate_incidence_elevation(sar_polygon, 10.5, 36.59, "IW")
   # result.mode             "IW" or "EW"
   # result.ref_lon          reference longitude (deg)
   # result.ref_lat          reference latitude (deg)
   # result.dfrac            0 = near range line, 1 = far range line
   # result.swath_width_km   across-track ground width of the footprint (km)
   # result.r_s_km           inverted satellite geocentric distance (km)
   # result.incidence_deg    incidence angle from the local vertical (deg)
   # result.elevation_deg    90 - incidence (deg, exact)

``mode`` defaults to ``"IW"``. A dedicated generic entry point is also
available: ``estimate_geometry(sar_polygon, ref_lon=10.5, ref_lat=36.59,
mode="IW")`` fills ``incidence_deg`` / ``elevation_deg`` on the
:class:`~s1_geometry.geometry.S1Geometry` result (the heading fields remain
the IW model). Catalogue-oriented named functions are exported as well:

* ``sar_incidence_angle(sar_polygon, ref_lon, ref_lat, mode)`` -- incidence
  in degrees; the function behind the ``sar_incidence_angle`` column.
* ``sar_elevation_angle(sar_polygon, ref_lon, ref_lat, mode)`` -- elevation
  in degrees; the function behind the ``sar_elevation_angle`` column.
* ``swath_width_km(sar_polygon)`` -- the across-track ground width of the
  footprint (km).
* ``dfrac_of_point(sar_polygon, ref_lon, ref_lat)`` -- the across-track
  position of the point (0 = near range line, 1 = far range line).

The reference point must lie **inside the SAFE footprint** (the 4-vertex
parallelogram, boundary included).  ``estimate_incidence_elevation`` (and
therefore ``estimate_geometry`` and the catalogue-oriented functions) raises
a ``ValueError`` when the point is outside, instead of returning an
extrapolated value.

CLI
---

::

   s1-geometry heading footprint.geojson --lat 36.59
   s1-geometry heading '{"type": "Polygon", "coordinates": [[...]]}' --lat 36.59 --json
   s1-geometry angles footprint.geojson --lon 10.5 --lat 36.59 --mode IW
   s1-geometry angles footprint.geojson --lon 10.5 --lat 36.59 --json

The polygon argument is a GeoJSON file path or an inline JSON string: a
GeoJSON Polygon object, a GeoJSON Feature, or a bare ring of [lon, lat]
vertices; ``--json`` prints a machine-readable result.

The model
---------

The heading of an IW acquisition depends almost only on the latitude and the
orbit direction; the four satellites agree within a fraction of a degree, so
two pooled degree-4 polynomials (ascending / descending) describe the whole
constellation (fit window \|lat\| <= 70 deg; beyond it they extrapolate and the
uncertainty model carries the extra error). The orbit direction is the
*pass* (flight) direction derived from the footprint: the flight runs along
the burst (along-track) edges of the 4-vertex SAFE parallelogram, the
shorter opposite-edge family for IW. The rule is order-robust -- it handles
both the standard CDSE ``GeoFootprint`` vertex order (counter-clockwise) and
the clockwise ring copied from the CDSE product page, and reduces to the
historical ``dlat = lat(v[1]) - lat(v[0])`` rule for the standard order.
A northbound pass is ascending, a southbound pass descending; a product at
the polar turnaround (**turning**, flight north component below 0.3 deg in
magnitude) receives a near-due-west heading of 270 deg.

Validation
----------

Against the measured OWI heading of the 11107 products used for the fit:
mean \|error\| 0.749 deg, median 0.139 deg, p95 2.90 deg; 98.88 % of the
headings fall within the returned U(lat) and 99.96 % within 2 U(lat). In the
fitted window \|lat\| <= 70 deg the mean error is 0.14-0.25 deg per latitude
band and the p95 stays below 1.19 deg.

The angle model
---------------

The radar is left-looking and steered to a quasi-constant **off-nadir angle**
``sigma`` per sub-swath. Two consistent models are combined (verified against
the ESA / SentiWiki documented angles to <= 0.07 deg):

* **Constant slant range** -- the CDSE catalogue footprint. The across-track
  ground width ``W`` of the ``GeoFootprint`` is a monotonic function of the
  satellite geocentric distance ``R_s``; the catalogue polygon is computed
  with the mean Earth radius (6371 km), so inverting ``W -> R_s`` (bisection)
  calibrates the satellite altitude from the footprint alone (``R_s - 6371``
  km) with no orbit file needed.
* **Off-nadir** -- the actual radar steering. With the satellite distance
  ``R_s``, the WGS84 geocentric radius ``R_e`` of the reference point and the
  off-nadir angle ``sigma`` (the linear roll law between the documented
  min / max-altitude values per sub-swath), the geocentric incidence is the
  closed form ``theta_c = asin[(R_s/R_e) * sin(sigma)]``. The local
  (ellipsoid-normal) incidence adds the deflection of the vertical (max
  ~0.19 deg near \|lat\| = 45 deg) and the elevation is exactly
  ``90 deg - incidence``.

The across-track width ``W`` is the spherical distance (6371 km sphere)
between the midpoints of the two range lines of the 4-vertex SAFE ring. In
the standard CDSE vertex order (GRDH/GRDM/SLC/OCN) the range lines are the
``v0-v1`` (far) and ``v2-v3`` (near) edges, for ascending and descending
alike; rings with a rotated vertex order (e.g. ETA) fall back to the other
edge pair when the standard width does not invert to a valid orbit altitude.
The across-track position ``dfrac`` of the reference point (0 = near, 1 = far)
selects ``sigma`` by a linear roll law between the near and far sub-swath
edges.

Angle validation
----------------

Validated on 66243 CDSE IW/EW rings: 98.9 % (IW) and 90.8 % (EW) of the
products invert an orbit altitude in the documented 697.6 to 725.8 km band;
the full chain (ring -> W -> R_s -> incidence) is simulated to <= 0.11 deg
(IW) and <= 0.05 deg (EW) maximum error. The off-nadir model reproduces the
documented per-sub-swath angles to <= 0.07 deg. Full-swath footprints are the
supported case; a partial-swath EW footprint whose width does not invert to
an orbit altitude in the documented band raises :class:`ValueError`.
