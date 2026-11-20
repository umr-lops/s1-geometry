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
   # result.elevation_deg    radar elevation angle = off-nadir angle sigma (deg,
   #                         ESA L2 OCN *ElevationAngle convention)

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

WV (Wave mode) from the OData metadata only
-------------------------------------------

A WV product footprint (the CDSE OData ``Footprint``, a WKT
``MULTIPOLYGON``; or the ``GeoFootprint`` GeoJSON) holds one small ring per
imagette, in time order, alternating between the WV1 and WV2 subswaths.  The
angles are returned per reference point; the satellite geocentric distance
comes from the metadata-only orbit model
(:func:`~s1_geometry.orbit.satellite_radius_km`), so nothing else is needed::

   from s1_geometry import estimate_geometry, estimate_wv_angles

   geo = estimate_geometry(odata_footprint, ref_lon=-76.7, ref_lat=-29.3, mode="WV")
   # geo.incidence_deg, geo.elevation_deg, geo.heading_deg, geo.subswath, geo.r_s_km

   wv = estimate_wv_angles(odata_footprint, None, lons, lats)  # vectorised

To add the columns to a table (one row per reference point, the product
footprint in a column; needs the ``parquet`` extra)::

   from s1_geometry.table import add_wv_angle_columns

   df = add_wv_angle_columns(df, footprint="Footprint", lon="ref_lon", lat="ref_lat")
   # new columns: sar_incidence_angle, sar_elevation_angle, sar_ground_heading,
   #              sar_wv_subswath

When the product footprint is not available **in time order** but the
imagettes are -- the OData ``Imagettes`` expansion (each with its ``index``
and ``GeoFootprint``) or a set of STAC imagette items -- rebuild it first::

   from s1_geometry import estimate_wv_angles, wv_footprint_from_imagettes

   rings = wv_footprint_from_imagettes(product["Imagettes"])   # sorted by "index"
   wv = estimate_wv_angles(rings, None, lons, lats, on_outside="nan")
   # wv.ring_index == imagette index - 1, wv.subswath == 1 / 2 (WV1 / WV2)

A **single imagette polygon is not enough**: the two subswath chains of the
product are needed to locate the track plane, so a product needs at least 4
imagettes.  The rings must be **in time order** (the order of the CDSE
``Footprint`` / ``GeoFootprint``, where ring *i* is imagette ``index`` *i + 1*
of the ``WVImagettes`` list); a re-ordered geometry (e.g. the ``geometry``
column of the daily Parquets, normalised by a GIS library) is detected and
rejected.  Imagettes split in two at the antimeridian in the OData
``Footprint`` are re-joined.  Validated on 164 Sentinel-1A / -1C / -1D products
(2.97 million ESA L2 OCN cells): RMSE 0.036 / 0.041 / 0.10 deg for the
incidence / elevation / heading (see :doc:`certification`).  Points outside every imagette
raise ``ValueError`` (or give NaN with ``on_outside="nan"``).

The orbit model reproduces the POEORB state vectors of Sentinel-1A, -1C and
-1D to 0.10 km (std), 0.29 km (max), i.e. ~0.02 deg of incidence.  Note
that the satellite geocentric distance (7065 to 7082 km) is not
``6371 km +`` the documented 697.6-725.8 km altitude (a height above the
WGS84 ellipsoid).

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
   s1-geometry angles footprint.geojson --lon 10.5 --lat 36.59 --mode SM --beam S3 --json
   s1-geometry wv "geography'SRID=4326;MULTIPOLYGON (((...)), ...)'" --lon -76.7 --lat -29.3

The polygon argument is a file path or an inline string: a GeoJSON Polygon
or MultiPolygon object, a GeoJSON Feature, a bare ring of [lon, lat]
vertices, or WKT (the CDSE OData ``Footprint``).  ``--json`` prints a
machine-readable result; ``wv`` also accepts ``--r-s`` to force the satellite
distance (km).

How the angles are computed, how the models were calibrated and how accurate
they are is explained in :doc:`methodology` and :doc:`certification`.
