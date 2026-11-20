Validation against ESA L2 OCN
=============================

The library is checked against the ESA Sentinel-1 Level-2 OCN products,
which carry, on their ``owi`` grid, the incidence angle
(``owiIncidenceAngle``), the radar elevation angle (``owiElevationAngle``)
and the ground heading (``owiHeading``) computed by ESA from the full orbit
and attitude data.  s1-geometry is run on the same grid cells **from the
catalogue footprint only** and the difference (model − ESA) is analysed.
How the models work is explained in :doc:`methodology`.

All the figures and tables of this page are produced by
``scripts/make_figures.py`` (examples, IW / EW sample) and
``scripts/validate_wv_ocn.py`` (WV); the per-example statistics are in
``_static/certification/stats.json``.

.. include:: _generated/example_products.inc


Examples
--------

One ascending and one descending product per mode.  For IW and EW, each row
shows the ESA field, the error map and the error distribution; for WV, the
mean of each imagette (ESA: circles, model: lines) and the error
distribution over all the cells.

IW
^^

Ascending: |iw_asc_product|

.. image:: _static/figures/example_iw_asc.png
   :class: only-light
   :width: 100%
   :align: center
   :alt: IW ascending example: ESA fields, error maps and error distributions

.. image:: _static/figures/example_iw_asc_dark.png
   :class: only-dark
   :width: 100%
   :align: center
   :alt: IW ascending example: ESA fields, error maps and error distributions

Descending: |iw_desc_product|

.. image:: _static/figures/example_iw_desc.png
   :class: only-light
   :width: 100%
   :align: center
   :alt: IW descending example: ESA fields, error maps and error distributions

.. image:: _static/figures/example_iw_desc_dark.png
   :class: only-dark
   :width: 100%
   :align: center
   :alt: IW descending example: ESA fields, error maps and error distributions

EW
^^

Ascending: |ew_asc_product|

.. image:: _static/figures/example_ew_asc.png
   :class: only-light
   :width: 100%
   :align: center
   :alt: EW ascending example: ESA fields, error maps and error distributions

.. image:: _static/figures/example_ew_asc_dark.png
   :class: only-dark
   :width: 100%
   :align: center
   :alt: EW ascending example: ESA fields, error maps and error distributions

Descending: |ew_desc_product|

.. image:: _static/figures/example_ew_desc.png
   :class: only-light
   :width: 100%
   :align: center
   :alt: EW descending example: ESA fields, error maps and error distributions

.. image:: _static/figures/example_ew_desc_dark.png
   :class: only-dark
   :width: 100%
   :align: center
   :alt: EW descending example: ESA fields, error maps and error distributions

The IW / EW error grows towards the edges of the swath and of each
sub-swath: the off-nadir angle is interpolated as a straight line from the
near to the far range line, whereas IW and EW are made of 3 and 5 sub-swaths,
each with its own beam.  EW, twice as wide, carries the larger error.

WV
^^

Ascending: |wv_asc_product|

.. image:: _static/figures/example_wv_asc.png
   :class: only-light
   :width: 100%
   :align: center
   :alt: WV ascending example: imagette means and error distributions

.. image:: _static/figures/example_wv_asc_dark.png
   :class: only-dark
   :width: 100%
   :align: center
   :alt: WV ascending example: imagette means and error distributions

Descending: |wv_desc_product|

.. image:: _static/figures/example_wv_desc.png
   :class: only-light
   :width: 100%
   :align: center
   :alt: WV descending example: imagette means and error distributions

.. image:: _static/figures/example_wv_desc_dark.png
   :class: only-dark
   :width: 100%
   :align: center
   :alt: WV descending example: imagette means and error distributions

Statistics of the examples
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. include:: _generated/example_stats.inc


WV: Sentinel-1A, -1C and -1D
----------------------------

``scripts/validate_wv_ocn.py`` runs the metadata-only chain -- the product
footprint (the ``manifest.safe`` imagette rings, identical to the CDSE OData
``Footprint``), ``r_s`` from the orbit model of :mod:`s1_geometry.orbit`,
``window="auto"`` -- on every ``owi`` cell of every imagette (the ``osw``
imagette centre for the products before mid-2020, whose WV ``owi`` grid is
empty).  Sample: one day per month (three products per day spread over the
day) for S1C (2025-2026) and S1D (2026), S1A 2025-2026 plus one day per year
2015-2024; 164 products, latitudes -80 to +52°, up to 143 imagettes,
2,974,963 cells.  Per product statistics:
``_static/certification/wv_validation_s1acd.csv``.

Look-angle law recalibration.  The first law (a linear latitude term fitted
on 64 mostly southern-hemisphere products) had an incidence bias growing in
the northern hemisphere (+0.6° at 30 N, +0.9° at 45-50 N).  The latitude
term is now a°ree-5 polynomial, fitted by
``scripts/calibrate_wv_look_angle.py`` (one Newton step per imagette on the
incidence error) on 62,588 imagettes of 1,205 S1C / S1D products of *other*
days than the validation sample (latitudes -79.5 to +86°); S1A is not used
in the fit.

.. list-table::
   :header-rows: 1

   * - satellite
     - products
     - incidence RMSE
     - elevation RMSE
     - heading RMSE
     - incidence max error
   * - S1A
     - 75
     - 0.021°
     - 0.033°
     - 0.097°
     - 0.159°
   * - S1C
     - 57
     - 0.019°
     - 0.032°
     - 0.096°
     - 0.163°
   * - S1D
     - 32
     - 0.069°
     - 0.065°
     - 0.123°
     - 0.475°
   * - all
     - 164
     - 0.036°
     - 0.041°
     - 0.103°
     - 0.475°
   * - first law
     - 163
     - 0.233°
     - 0.224°
     - 0.098°
     - 1.634°

Subswath (WV1 / WV2) labels are 100 % correct and no cell falls outside the
footprint.  The incidence bias is below 0.08° in every latitude band.  The
heading p95 error of the worst product is 0.36° (S1D products through the
polar turn), hence the WV heading uncertainty of 0.4°.  Long products need
the sliding-window track fit (``window="auto"``): with a single track geometry
their heading RMSE is 0.36 instead of 0.09°.

Excluded from the statistics:

* the S1C commissioning phase (before 2025-03-08, no S1C POEORB, orbit not
  yet the reference one);
* end-of-datatake products with 2 imagettes (at least 4 are needed);
* reference headings known to be wrong: ``owiHeading`` of the IPF 003.31
  products (2021: 357° on a descending pass), imagettes whose
  ``owiHeading`` is constant at 90.0004° (an ESA default value, seen near
  the antimeridian), and ``oswHeading`` (before mid-2020; a different
  quantity, 1-6° away from the ground-track heading).

The centre of an imagette ring is the mean of its *closed* ring vertices (the
first vertex counted twice), the convention of the CDSE footprints the
look-angle law is calibrated with.

.. image:: _static/figures/validation_wv_latitude.png
   :class: only-light
   :width: 100%
   :align: center
   :alt: WV incidence bias and heading RMSE per product against latitude

.. image:: _static/figures/validation_wv_latitude_dark.png
   :class: only-dark
   :width: 100%
   :align: center
   :alt: WV incidence bias and heading RMSE per product against latitude


IW and EW: sample of 2026 products
----------------------------------

Two Sentinel-1D products every fourth day of 2026 for IW and for EW (every
fourth grid cell).  Products whose footprint-derived pass direction differs
from the catalogue ``orbitDirection``, or whose ESA grid is not covered by
the model (partial-swath footprints), are left out; the per-product values
are in ``_static/certification/iw_ew_validation_2026.csv``.

Kept: 41 of 46 IW and 19 of 29 EW products.  Median per-product RMSE below
60° of latitude: IW incidence 0.54°, heading 0.50°; EW incidence 1.44°,
heading 0.50°.  Beyond 60°: IW 2.5° / 2.5°, EW 2.1° / 5.5° (heading
polynomials extrapolated, polar turn).  Left out: 5 IW products (72-77° S)
and 8 EW products (mostly beyond 60° N) whose footprint-derived pass
direction is wrong, and 2 EW products only partly covered.

.. image:: _static/figures/validation_iw_ew_latitude.png
   :class: only-light
   :width: 100%
   :align: center
   :alt: IW and EW incidence bias and heading RMSE per product against latitude

.. image:: _static/figures/validation_iw_ew_latitude_dark.png
   :class: only-dark
   :width: 100%
   :align: center
   :alt: IW and EW incidence bias and heading RMSE per product against latitude


Orbit model stability
---------------------



The figure is on the :doc:`methodology` page.


SM (stripmap)
-------------

No SM Level-2 product exists in the CDSE archive (2014-2026), so the SM
chain is **specification-based and not validated** against ESA data (see
:doc:`methodology`).
