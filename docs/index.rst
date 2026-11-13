s1-geometry
===========

Tool to estimate geometric quantities (such as incidence, elevation or ground
heading angle) from a Sentinel-1 SAR footprint.

The current release (v0.1) provides the **IW ground heading** as a function
of the reference point and the SAFE footprint polygon, together with a
latitude-dependent uncertainty::

   from s1_geometry import estimate_geometry

   geometry = estimate_geometry(footprint_polygon, ref_lat=36.592652)
   print(geometry.heading_deg, geometry.heading_uncertainty_deg)

Incidence angle and elevation angle estimators are planned.

.. code-block:: bash

   s1-geometry heading footprint.geojson --lat 36.59

The model is a per-orbit-direction degree-4 polynomial in latitude, fitted on
11107 IW OCN products; 98.88 % of the measured headings fall within the
returned uncertainty. See :doc:`usage`.

Test suite: ``python -m pytest -q``.

.. toctree::
   :maxdepth: 2

   usage
   api/s1_geometry
