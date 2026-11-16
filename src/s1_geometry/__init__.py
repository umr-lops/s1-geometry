"""s1_geometry: geometric quantities from a Sentinel-1 SAR footprint.

Current release (v0.3): the IW **ground heading** as a function of the
reference point and the SAFE footprint (see :mod:`s1_geometry.heading`),
the **incidence / elevation angles** of IW, EW and SM acquisitions from the
SAFE footprint and a reference point (see :mod:`s1_geometry.angles`; the SM
stripmap chain is spec-based and not validated against ESA L2 data), and the
**per-point incidence / elevation / heading** of WV (wide-swath) acquisitions
from the burst-ring MultiPolygon and the satellite distance (see
:mod:`s1_geometry.wv`).
"""

from s1_geometry.angles import (
    MODE_SIGMA,
    R_E_KM,
    R_S_MAX_KM,
    R_S_MIN_KM,
    S1Angles,
    SM_BEAM_SIGMA,
    dfrac_of_point,
    estimate_incidence_elevation,
    sar_elevation_angle,
    sar_incidence_angle,
    swath_width_km,
)
from s1_geometry.footprint import (
    ASCENDING,
    DESCENDING,
    DLAT_TURN_THRESHOLD_DEG,
    TURNING,
    as_multipolygon,
    load_sar_polygons,
    orbit_direction_from_ring,
    point_in_polygon,
    polygon_index_for_point,
    ring_to_geojson,
)
from s1_geometry.geometry import S1Geometry, estimate_geometry
from s1_geometry.wv import (
    S1WVAngles,
    WV_THETA_DOC,
    estimate_wv_angles,
    wv_theta_deg,
)
from s1_geometry.heading import (
    HEADING_FIT_CUTOFF_ABS_LAT,
    TURN_HEADING_DEG,
    TURNING_LATITUDES,
    heading_from_lat,
    sar_ground_heading,
    sar_ground_heading_uncertainty,
    sar_heading_and_uncertainty,
)

__version__ = "0.3.0"

__all__ = [
    "ASCENDING",
    "DESCENDING",
    "DLAT_TURN_THRESHOLD_DEG",
    "HEADING_FIT_CUTOFF_ABS_LAT",
    "MODE_SIGMA",
    "R_E_KM",
    "R_S_MAX_KM",
    "R_S_MIN_KM",
    "S1Angles",
    "S1Geometry",
    "S1WVAngles",
    "SM_BEAM_SIGMA",
    "TURNING",
    "TURN_HEADING_DEG",
    "TURNING_LATITUDES",
    "WV_THETA_DOC",
    "__version__",
    "as_multipolygon",
    "dfrac_of_point",
    "estimate_geometry",
    "estimate_incidence_elevation",
    "estimate_wv_angles",
    "heading_from_lat",
    "load_sar_polygons",
    "orbit_direction_from_ring",
    "point_in_polygon",
    "polygon_index_for_point",
    "ring_to_geojson",
    "sar_elevation_angle",
    "sar_ground_heading",
    "sar_ground_heading_uncertainty",
    "sar_heading_and_uncertainty",
    "sar_incidence_angle",
    "swath_width_km",
    "wv_theta_deg",
]
