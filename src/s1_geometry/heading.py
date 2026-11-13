"""S1 **IW** SAR ground heading as a function of the reference point and the orbit.

The heading of an IW acquisition depends almost only on the latitude and the
orbit direction; the four satellites agree within a fraction of a degree, so
two pooled degree-4 polynomials (ascending / descending) describe the whole
constellation. Fitted on 13 consecutive days of S1 IW OCN OWI data from all
four satellites (11107 products, 44428 grid corners; fit window
|lat| <= 70 deg, headings wrapped to the per-orbit median):

* ascending  : n = 4739, fit residual std 0.30 deg (p95 0.64)
* descending : n = 5104, fit residual std 0.32 deg (p95 0.65)

Beyond |lat| = 70 the polynomials extrapolate and the uncertainty model
carries the extra error. Rare **turning** products at the polar turnaround
receive a near-due-west heading of 270 deg (the satellite keeps its westward
ground drift) with the large turning-latitude uncertainty.

The coefficients are embedded below; they are the validated deliverable of
the SOBA SAR-metadata retrieval work (2026-11-12). Do not change them
without refitting and revalidating.
"""

from __future__ import annotations

import numpy as np

from s1_geometry.footprint import (
    ASCENDING,
    DESCENDING,
    TURNING,
    orbit_direction_from_ring,
)

#: Heading (deg, clockwise from North) of turning products at the polar
#: turnaround: the satellite always keeps its westward ground drift.
TURN_HEADING_DEG = 270.0

#: Latitude (deg, absolute value) beyond which the heading polynomials
#: extrapolate and the uncertainty switches from its degree-4 fit to the
#: piecewise-linear high-latitude table.
HEADING_FIT_CUTOFF_ABS_LAT = 70.0

#: Measured turning latitudes (deg, absolute value), strongly asymmetric:
#: 77.5 S (products reach the turn) / 84.5 N (products stop just before).
TURNING_LATITUDES = {"N": 84.5, "S": 77.5}

# (reference_heading_deg, degree-4 coefficients high-to-constant) fitted on
# the per-orbit median-wrapped headings of the 11107 IW OCN products.
_HEADING_FIT: dict[str, tuple[float, tuple[float, ...]]] = {
    ASCENDING: (
        349.1359,
        (-5.386606409525314e-07, 2.6534165330268556e-05,
         -1.0060918442871438e-04, 4.501394674817113e-02,
         -1.068609606919818),
    ),
    DESCENDING: (
        191.1954,
        (5.33163332650265e-07, -2.6570027177813648e-05,
         1.2037546134848548e-04, -4.6026655901221264e-02,
         7.772727186048888e-01),
    ),
}

# Per-hemisphere heading uncertainty (deg): (range-1 degree-4 coefficients
# high-to-constant in |lat|, (|lat|, p95) table for range 2, clamp latitude).
_UNCERTAINTY: dict[str, tuple[tuple[float, ...], list[list[float]], float]] = {
    "N": (
        (1.5252441589102288e-07, 3.158429132186365e-06,
         -9.378091381906103e-04, 5.0136914275963895e-02,
         1.8501203797241456e-01),
        [[70.5, 4.054490372763409], [71.5, 4.348602435788118],
         [72.5, 4.722937510378863], [73.5, 5.394590532889538],
         [74.5, 5.895073821413549], [75.5, 6.639534042218128],
         [76.5, 7.54677108362001], [77.5, 8.919053883363013],
         [78.5, 10.391993024788878], [79.5, 12.742941355980616],
         [80.5, 15.337671630804026], [81.5, 19.76108050355947],
         [82.5, 24.618360683146246], [83.5, 31.763574779138647],
         [84.5, 35.69177237473567]],
        84.5,
    ),
    "S": (
        (-2.0092368424904635e-07, 6.0944483375833336e-05,
         -3.273997278825826e-03, 8.257959296137626e-02,
         1.9030338660016521e-01),
        [[70.5, 7.640712701689401], [71.5, 9.529721533554614],
         [72.5, 12.068154604405782], [73.5, 15.623203751600194],
         [74.5, 19.96623526419088], [75.5, 27.50049254160709],
         [76.5, 44.29388010065139], [77.5, 21.568830322320565]],
        77.5,
    ),
}


def heading_from_lat(lat: float, orbit_direction: str) -> float:
    """Ground heading (deg, clockwise from North) of an IW acquisition.

    ``orbit_direction`` is one of :data:`~s1_geometry.footprint.ASCENDING`,
    :data:`~s1_geometry.footprint.DESCENDING` or
    :data:`~s1_geometry.footprint.TURNING`; turning products receive
    :data:`TURN_HEADING_DEG`.
    """
    if not isinstance(orbit_direction, str) or orbit_direction not in (
        ASCENDING, DESCENDING, TURNING,
    ):
        raise ValueError(f"unknown orbit direction: {orbit_direction!r}")
    if orbit_direction == TURNING:
        return TURN_HEADING_DEG
    reference, coef = _HEADING_FIT[orbit_direction]
    value = reference + np.polyval(coef, float(lat))
    return float(np.remainder(value, 360.0))


def sar_ground_heading(
    ref_lat: float, sar_polygon,
) -> float:
    """Ground heading (deg) of the IW acquisition at the reference point.

    The orbit direction is derived from ``sar_polygon`` (the SAFE footprint);
    the heading itself depends only on ``ref_lat`` (the reference point
    latitude) and the orbit direction. This is the function behind the
    ``sar_ground_heading`` column of the SOBA catalogue.
    """
    return heading_from_lat(ref_lat, orbit_direction_from_ring(sar_polygon))


def sar_ground_heading_uncertainty(ref_lat: float) -> float:
    """Latitude- and hemisphere-dependent heading uncertainty (deg).

    Per hemisphere, the p95 of the |corner error| versus the degree-4 fit,
    binned by |lat|: |lat| <= 70 is a degree-4 polynomial in |lat|;
    |lat| > 70 is a piecewise-linear table on 1-deg bins, clamped at the
    measured turning latitude (84.5 N / 77.5 S). Users can decide for
    themselves whether the heading of a given (high-latitude) acquisition is
    reliable instead of rejecting it.
    """
    hemisphere = "N" if ref_lat >= 0 else "S"
    abs_lat = abs(float(ref_lat))
    coef, table, _clamp = _UNCERTAINTY[hemisphere]
    if abs_lat <= HEADING_FIT_CUTOFF_ABS_LAT:
        return float(np.polyval(coef, abs_lat))
    table = np.asarray(table, dtype=float)
    return float(np.interp(abs_lat, table[:, 0], table[:, 1],
                           left=table[0, 1], right=table[-1, 1]))


def sar_heading_and_uncertainty(
    ref_lat: float, sar_polygon,
) -> tuple[float, float]:
    """``(heading, heading_uncertainty)`` of the IW acquisition at ref_lat."""
    orbit = orbit_direction_from_ring(sar_polygon)
    return (heading_from_lat(ref_lat, orbit), sar_ground_heading_uncertainty(ref_lat))
