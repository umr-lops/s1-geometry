"""Sentinel-1 satellite geocentric distance from metadata only.

The CDSE OData metadata of a Sentinel-1 product carries no orbit state
vectors, but Sentinel-1 flies a frozen, tightly maintained orbit: the
satellite geocentric distance ``r_s`` is, to a fraction of a kilometre, a
function of the latitude and the pass direction alone, identical for
Sentinel-1A, -1C and -1D.  :func:`satellite_radius_km` evaluates that
relation, so the WV angle chain (which needs ``r_s``) runs from the
``Footprint`` and ``orbitDirection`` metadata, without orbit files.

The latitude argument is the latitude of the **imaged point** (the reference
point), not of the sub-satellite point: the model is fitted against the
target of a right-looking beam at 26 deg off-nadir (midway between WV1 and
WV2).

Calibration (``scripts/calibrate_orbit_radius.py``): 37444 POEORB state
vectors, Sentinel-1A (2026-06-29), -1C (2025-12-10) and -1D (2026-05-11,
2026-10-29), evaluated at the WV1 (20.9 deg) and WV2 (31.5 deg) targets:
residual std 0.10-0.11 km, maximum 0.29 km, bias per satellite / day below
0.04 km.  An ``r_s`` error of 0.3 km moves the WV incidence by ~0.02 deg.

Note that ``r_s`` (7065.3 to 7082.4 km) is a geocentric distance; the
documented 697.6 to 725.8 km Sentinel-1 altitude band is the height above
the WGS84 ellipsoid, *not* ``r_s - 6371 km``.
"""

from __future__ import annotations

import numpy as np

from s1_geometry.footprint import ASCENDING, DESCENDING

__all__ = [
    "R_S_ORBIT_MAX_KM",
    "R_S_ORBIT_MIN_KM",
    "satellite_radius_km",
]

#: Range (km) of the Sentinel-1 satellite geocentric distance over the orbit
#: (POEORB, all satellites): 7065.3 to 7082.4 km, padded by 2 km.
R_S_ORBIT_MIN_KM = 7063.0
R_S_ORBIT_MAX_KM = 7085.0

# Polynomial coefficients (low to high degree) in x = lat / 90 of r_s (km),
# per pass direction, fitted by scripts/calibrate_orbit_radius.py.  Do not
# change them without refitting and revalidating.
_R_S_FIT: dict[str, tuple[float, ...]] = {
    ASCENDING: (
        7076.9765835159615, -13.186957368688967, -7.795022836716421,
        5.354418438834789, 6.550925803174891, -0.5939857765482541,
        -2.4460358440142844, 0.021552186640470632, 0.6867636251797369,
        -0.007846343710712007, -0.1418717232867006,
    ),
    DESCENDING: (
        7076.974486095424, -13.191514914827016, -7.764834343485362,
        5.343435049395258, 6.444060405591145, -0.5498051456494308,
        -2.2709750567085734, -0.0009186917569631949, 0.5463103859044868,
        -0.022176652969448114, -0.09302015717889646,
    ),
}

#: Latitude range (deg) of the fit; latitudes beyond are clamped.
_FIT_LAT_RANGE = (-78.65, 84.91)


def satellite_radius_km(lat_deg, orbit_direction: str):
    """Satellite geocentric distance ``r_s`` (km) when imaging ``lat_deg``.

    Parameters
    ----------
    lat_deg:
        Latitude (deg) of the imaged point, scalar or array.  Values outside
        the orbit's reach (-78.65 to 84.91 deg) are clamped.
    orbit_direction:
        ``"ascending"`` or ``"descending"`` (case-insensitive, e.g. the CDSE
        ``orbitDirection`` attribute ``"ASCENDING"``).

    Returns
    -------
    float or numpy.ndarray
        ``r_s`` in km (a float for a scalar input).
    """
    key = str(orbit_direction).lower()
    if key not in _R_S_FIT:
        raise ValueError(
            f"unknown orbit direction {orbit_direction!r} "
            "(expected 'ascending' or 'descending')")
    lat = np.clip(np.asarray(lat_deg, dtype=float), *_FIT_LAT_RANGE)
    value = np.polynomial.polynomial.polyval(lat / 90.0, _R_S_FIT[key])
    return float(value) if np.ndim(value) == 0 else value
