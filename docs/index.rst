s1-geometry
===========

Tool to estimate geometric quantities (such as incidence, elevation or ground
heading angle) from a Sentinel-1 SAR footprint.

The current release (v0.2) provides the **IW ground heading** as a function
of the reference point and the SAFE footprint polygon, together with a
latitude-dependent uncertainty, and the **incidence / elevation angles** of
IW and EW acquisitions at a reference point::

   from s1_geometry import estimate_geometry, estimate_incidence_elevation

   geometry = estimate_geometry(footprint_polygon, ref_lat=36.592652)
   print(geometry.heading_deg, geometry.heading_uncertainty_deg)

   angles = estimate_incidence_elevation(footprint_polygon, 10.5, 36.59, "IW")
   print(angles.incidence_deg, angles.elevation_deg)

.. code-block:: bash

   s1-geometry heading footprint.geojson --lat 36.59
   s1-geometry angles footprint.geojson --lon 10.5 --lat 36.59 --mode IW

The model is a per-orbit-direction degree-4 polynomial in latitude, fitted on
11107 IW OCN products; 98.88 % of the measured headings fall within the
returned uncertainty. See :doc:`usage`.

Test suite: ``python -m pytest -q``.

.. toctree::
   :maxdepth: 2

   usage
   api/s1_geometry
