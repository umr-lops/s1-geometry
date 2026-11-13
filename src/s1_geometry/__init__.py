"""s1_geometry: geometric quantities from a Sentinel-1 SAR footprint.

Current release (v0.1): the IW **ground heading** as a function of the
reference point and the SAFE footprint (see :mod:`s1_geometry.heading`).
Incidence angle and elevation angle estimators are planned.
"""

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

__version__ = "0.1.0"

__all__ = [
    "ASCENDING",
    "DESCENDING",
    "DLAT_TURN_THRESHOLD_DEG",
    "HEADING_FIT_CUTOFF_ABS_LAT",
    "S1Geometry",
    "TURNING",
    "TURN_HEADING_DEG",
    "TURNING_LATITUDES",
    "__version__",
    "estimate_geometry",
    "heading_from_lat",
    "load_sar_polygons",
    "orbit_direction_from_ring",
    "point_in_polygon",
    "polygon_index_for_point",
    "ring_to_geojson",
    "sar_ground_heading",
    "sar_ground_heading_uncertainty",
    "sar_heading_and_uncertainty",
]
