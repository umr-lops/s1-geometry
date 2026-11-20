s1-geometry
===========

Incidence angle, radar elevation angle and ground heading of a Sentinel-1
SAR acquisition, computed **from the catalogue metadata only**: the
footprint polygon of the CDSE catalogue (OData ``Footprint`` /
``GeoFootprint``) and a reference point.  No SAR data, no orbit file.

* **IW / EW / SM**: incidence and elevation at a reference point from the
  4-corner footprint, and the IW ground heading model with a
  latitude-dependent uncertainty (SM is specification-based, not validated);
* **WV (Wave mode)**: incidence, elevation and heading at any point of the
  product, from the imagette footprints and a metadata-only orbit model --
  0.04° incidence RMSE against 164 ESA L2 OCN products of S1A, S1C and S1D.

::

   from s1_geometry import estimate_geometry

   iw = estimate_geometry(footprint, ref_lon=10.5, ref_lat=36.59, mode="IW")
   wv = estimate_geometry(odata_footprint, ref_lon=-76.7, ref_lat=-29.3, mode="WV")
   print(wv.incidence_deg, wv.elevation_deg, wv.heading_deg)

.. code-block:: bash

   pip install s1-geometry
   s1-geometry angles footprint.geojson --lon 10.5 --lat 36.59 --mode IW
   s1-geometry wv "<OData Footprint WKT>" --lon -76.7 --lat -29.3

.. image:: _static/figures/method_geometry.png
   :class: only-light
   :width: 100%
   :alt: Viewing geometry: satellite, nadir, imaged point and the angles

.. image:: _static/figures/method_geometry_dark.png
   :class: only-dark
   :width: 100%
   :alt: Viewing geometry: satellite, nadir, imaged point and the angles

Start with :doc:`methodology` (what goes in, how it works, how it was
calibrated), then :doc:`usage`; :doc:`certification` gives the comparison
with the ESA products.  Cite the software with its Zenodo DOI (see the
README).

.. toctree::
   :maxdepth: 2

   methodology
   usage
   certification
   api/s1_geometry
