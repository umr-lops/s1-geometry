"""s1_geometry: geometric quantities from a Sentinel-1 SAR footprint.

The IW **ground heading** as a function of the
reference point and the SAFE footprint (see :mod:`s1_geometry.heading`),
the **incidence / elevation angles** of IW, EW and SM acquisitions from the
SAFE footprint and a reference point (see :mod:`s1_geometry.angles`; the SM
stripmap chain is spec-based and not validated against ESA L2 data), and the
**per-point incidence / elevation / heading** of WV (Wave mode) acquisitions
from the imagette MultiPolygon (see :mod:`s1_geometry.wv`), with the satellite
distance from the metadata-only orbit model (see :mod:`s1_geometry.orbit`):
everything runs from the CDSE OData metadata, no orbit file or SAR data.
"""

from importlib.metadata import PackageNotFoundError, version

from s1_geometry.angles import (
    MODE_SIGMA,
    R_E_KM,
    R_S_MAX_KM,
    R_S_MIN_KM,
    SM_BEAM_SIGMA,
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
    as_multipolygon,
    load_sar_polygons,
    orbit_direction_from_ring,
    parse_wkt,
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
from s1_geometry.orbit import R_S_ORBIT_MAX_KM, R_S_ORBIT_MIN_KM, satellite_radius_km
from s1_geometry.wv import (
    WV_THETA_DOC,
    S1WVAngles,
    estimate_wv_angles,
    wv_theta_deg,
)

try:
    __version__ = version("s1-geometry")
except PackageNotFoundError:  # running from a source tree without install
    __version__ = "0+unknown"

__all__ = [
    "ASCENDING",
    "DESCENDING",
    "DLAT_TURN_THRESHOLD_DEG",
    "HEADING_FIT_CUTOFF_ABS_LAT",
    "MODE_SIGMA",
    "R_E_KM",
    "R_S_MAX_KM",
    "R_S_MIN_KM",
    "R_S_ORBIT_MAX_KM",
    "R_S_ORBIT_MIN_KM",
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
    "parse_wkt",
    "point_in_polygon",
    "polygon_index_for_point",
    "ring_to_geojson",
    "sar_elevation_angle",
    "sar_ground_heading",
    "sar_ground_heading_uncertainty",
    "sar_heading_and_uncertainty",
    "sar_incidence_angle",
    "satellite_radius_km",
    "swath_width_km",
    "wv_theta_deg",
]
