"""s1_geometry: geometric quantities from a Sentinel-1 SAR footprint.

Current release (v0.2): the IW **ground heading** as a function of the
reference point and the SAFE footprint (see :mod:`s1_geometry.heading`),
and the **incidence / elevation angles** of IW and EW acquisitions from the
SAFE footprint and a reference point (see :mod:`s1_geometry.angles`).
"""

from s1_geometry.angles import (
    MODE_SIGMA,
    R_E_KM,
    R_S_MAX_KM,
    R_S_MIN_KM,
    S1Angles,
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
    load_sar_polygons,
    orbit_direction_from_ring,
    point_in_polygon,
    polygon_index_for_point,
    ring_to_geojson,
)
from s1_geometry.geometry import S1Geometry, estimate_geometry
from s1_geometry.heading import (
    HEADING_FIT_CUTOFF_ABS_LAT,
    TURN_HEADING_DEG,
    TURNING_LATITUDES,
    heading_from_lat,
    sar_ground_heading,
    sar_ground_heading_uncertainty,
    sar_heading_and_uncertainty,
)

__version__ = "0.2.0"

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
    "TURNING",
    "TURN_HEADING_DEG",
    "TURNING_LATITUDES",
    "__version__",
    "dfrac_of_point",
    "estimate_geometry",
    "estimate_incidence_elevation",
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
]
