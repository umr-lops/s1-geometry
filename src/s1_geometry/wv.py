"""S1 **WV** SAR incidence, elevation and heading angles from the SAFE footprint.

The wide-swath (WV) acquisition imaged by the C-band radar (Sentinel-1D) uses
two fixed off-nadir subswaths (WV1 ~ 21 deg, WV2 ~ 33 deg) sampled in a
sequence of short bursts that alternate between the two subswaths.  The CDSE
``GeoFootprint`` is a *MultiPolygon*: one small (~14 km) ring per burst, the
rings forming two interleaved chains (the subswaths) along the ground track.

Unlike IW / EW, the WV footprint does not encode a single constant-slant-range
swath, so the satellite geocentric distance ``r_s_km`` is an *input* (from the
orbit or the metadata), and the angles are returned **per ground point** (the
look geometry differs across each burst's imaged cell).  The model:

* fit the track-plane normal ``n`` from the paired chain midpoints and build
  the in-plane along-track basis ``u`` (the max eigenvector of the ring centres
  *projected onto the track plane* -- the raw eigenvector is tilted out of
  plane by the two-chain zig-zag, up to ~1.3 deg, which would corrupt the
  heading), oriented by the ring time progression, and ``rperp = n x u``;
* break the **mirror ambiguity**: the two subswath assignments are both
  self-consistent (the sub-satellite point can be placed on either side of the
  track plane).  The radar is right-looking, so the look direction is
  ``Right = u x Up`` (``Up`` the scene radial); the normal is oriented opposite
  to the look (``n . Right < 0``), which fixes the physical placement.  In that
  convention the near subswath WV1 is the *larger*-offset chain;
* label the two subswath chains from the gap in the ring centres' track-plane
  offset (WV1 = the larger-offset chain, closer to nadir);
* for each burst, the ground-range offset ``d_k = h tan(theta_k)`` with the
  calibrated look-angle law :func:`wv_theta_deg`; the sub-satellite point is
  placed beyond the WV1 chain (``delta0 = d_m + d_km``);
* for each target point: the incidence (angle at the ground, from the local
  normal), the radar elevation angle (``incidence - gamma``, the ESA L2 OCN
  ``*ElevationAngle`` convention) and the heading (the in-plane tangent
  azimuth, clockwise from North).

Validated point-by-point against the ESA L2 OCN ``owi`` grid (incidence,
elevation, heading) on four products (two ascending, two descending): RMSE
0.19-0.44 deg (incidence), 0.18-0.41 deg (elevation), 0.09-0.32 deg (heading).
See :doc:`/certification`.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from s1_geometry.angles import R_E_KM, R_S_MAX_KM, R_S_MIN_KM
from s1_geometry.footprint import as_multipolygon

__all__ = [
    "S1WVAngles",
    "WV_THETA_DOC",
    "estimate_wv_angles",
    "wv_theta_deg",
]

#: Documented WV off-nadir band centers (deg) at min / max orbit altitude
#: (SentiWiki / S1 Product Definition; reference altitude 711.7 km).
WV_THETA_DOC = {
    1: (21.715, 20.11),   # WV1: (min-altitude, max-altitude) band centers
    2: (33.09, 31.49),    # WV2
}

_ALT_MIN, _ALT_MAX, _REF_ALT = 697.6, 725.8, 711.7

#: Calibrated latitude term ``(deg/deg, deg)`` of the look-angle law, fitted on
#: 64 WV OCN products (2026-01..06); residual std 0.086 deg (WV1) / 0.065 deg
#: (WV2).  Do not change without refitting and revalidating.
_WV_THETA_LAT = {1: (0.012076, 0.2450), 2: (0.009995, 0.6560)}


def wv_theta_deg(subswath: int, r_s_km: float, lat_deg: float) -> float:
    """WV look (off-nadir) angle (deg) for subswath ``subswath`` (1 or 2).

    The documented band-center at the reference altitude (711.7 km), linear in
    the satellite geocentric distance ``r_s_km``, plus a calibrated latitude
    term.  ``r_s_km`` is the satellite geocentric distance (km); the altitude
    is ``r_s_km - R_E_KM``.
    """
    if subswath not in WV_THETA_DOC:
        raise ValueError(f"unknown WV subswath {subswath!r} (expected 1 or 2)")
    c_min, c_max = WV_THETA_DOC[subswath]
    slope = (c_max - c_min) / (_ALT_MAX - _ALT_MIN)
    base = 0.5 * (c_min + c_max) + slope * (float(r_s_km) - R_E_KM - _REF_ALT)
    m, b = _WV_THETA_LAT[subswath]
    return base + m * float(lat_deg) + b


def _to_3d(lonlat: np.ndarray, r_e_km: float = R_E_KM) -> np.ndarray:
    """(N, 2) (lon, lat) deg -> (N, 3) ECEF (km) on a sphere of radius
    ``r_e_km`` (the CDSE catalogue convention)."""
    lon = np.radians(lonlat[:, 0])
    lat = np.radians(lonlat[:, 1])
    return r_e_km * np.stack(
        (np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)),
        axis=1,
    )


def _enu_basis(slon_deg: float, slat_deg: float):
    """East / North / Up unit vectors (ECEF) at (slon_deg, slat_deg)."""
    lonr, latr = math.radians(slon_deg), math.radians(slat_deg)
    e = np.array([-math.sin(lonr), math.cos(lonr), 0.0])
    n = np.array([-math.sin(latr) * math.cos(lonr),
                  -math.sin(latr) * math.sin(lonr), math.cos(latr)])
    u = np.array([math.cos(latr) * math.cos(lonr),
                  math.cos(latr) * math.sin(lonr), math.sin(latr)])
    return e, n, u


def _fit_plane_min(p: np.ndarray):
    """Normal (min-eigenvector) and centroid of the point cloud ``p``."""
    mu = p.mean(axis=0)
    s = p - mu
    _, v = np.linalg.eigh(s.T @ s)
    n = v[:, 0]
    n /= np.linalg.norm(n)
    return n, mu


@dataclass(frozen=True)
class S1WVAngles:
    """Per-point WV angles and the intermediate quantities.

    Attributes
    ----------
    mode:
        ``"WV"``.
    r_s_km:
        Satellite geocentric distance (km); the altitude is ``r_s_km - R_E_KM``.
    mean_lat_deg:
        Mean latitude (deg) of the footprint ring centres.
    n_rings:
        Number of burst rings in the footprint.
    lon_deg, lat_deg:
        (N,) target point longitudes / latitudes (deg), as passed in.
    incidence_deg:
        (N,) incidence angle at each point, from the local (geocentric)
        normal, in degrees.
    elevation_deg:
        (N,) radar elevation angle at each point (``incidence - gamma``, the
        ESA L2 OCN ``*ElevationAngle`` convention), in degrees.
    heading_deg:
        (N,) heading (deg, clockwise from North) of the burst each point maps
        to (the in-plane tangent azimuth at the sub-satellite point).
    ring_index:
        (N,) index of the burst ring each point maps to (nearest ring centre).
    subswath:
        (N,) 1 or 2, the subswath each point belongs to.
    """

    mode: str
    r_s_km: float
    mean_lat_deg: float
    n_rings: int
    lon_deg: np.ndarray
    lat_deg: np.ndarray
    incidence_deg: np.ndarray
    elevation_deg: np.ndarray
    heading_deg: np.ndarray
    ring_index: np.ndarray
    subswath: np.ndarray


def estimate_wv_angles(
    sar_multipolygon, r_s_km: float, lon, lat,
) -> S1WVAngles:
    """Estimate the WV incidence / elevation / heading angles at target points.

    Parameters
    ----------
    sar_multipolygon:
        The WV SAFE footprint: a GeoJSON MultiPolygon string / dict, or a
        sequence of (lon, lat) rings (one small ring per burst, in burst /
        time order, alternating between the two subswaths).
    r_s_km:
        Satellite geocentric distance (km).  WV cannot invert ``R_s`` from the
        footprint width (there is no single constant-slant-range swath), so it
        is an input; it must be in the documented 697.6 to 725.8 km altitude
        band.
    lon, lat:
        Target point longitudes / latitudes (deg): a scalar or an array (e.g.
        the centres of the ESA OCN ``owi`` grid cells).  The angles are
        returned per point.

    Returns
    -------
    S1WVAngles
        Per-point incidence / elevation / heading and the intermediate
        quantities (see :class:`S1WVAngles`).

    Raises
    ------
    ValueError
        If ``r_s_km`` is outside the documented altitude band, the footprint
        is malformed, or it does not show two distinct subswath chains.
    """
    r_s = float(r_s_km)
    if not (R_S_MIN_KM <= r_s <= R_S_MAX_KM):
        raise ValueError(
            f"r_s_km {r_s:.2f} is outside the documented 697.6 to 725.8 km "
            f"altitude band (satellite geocentric distance {R_S_MIN_KM:.1f} "
            f"to {R_S_MAX_KM:.1f} km)")
    rings = as_multipolygon(sar_multipolygon)
    if len(rings) < 4:
        raise ValueError(
            f"WV footprint needs at least 4 burst rings (two per subswath), "
            f"got {len(rings)}")

    centers = np.array([_to_3d(r).mean(axis=0) for r in rings])
    n_ring = centers.shape[0]
    mean_center = centers.mean(axis=0)
    mean_lat = float(np.mean([r[:, 1].mean() for r in rings]))

    # Track-plane normal from the paired chain midpoints (consecutive rings
    # are opposite subswaths).
    mids = (centers[:-1] + centers[1:]) / 2.0
    n_fp, _ = _fit_plane_min(mids)

    # In-plane along-track direction: the max eigenvector of the ring centres
    # (the cloud is elongated along track), oriented by the ring time
    # progression (the CDSE ring order is the burst order), and projected onto
    # the track plane (the raw eigenvector of the two-chain zig-zag is tilted
    # out of plane by up to ~1.3 deg, which would corrupt the heading).
    sc = centers - mean_center
    _, v = np.linalg.eigh(sc.T @ sc)
    u = v[:, -1]
    if u @ (centers[-1] - centers[0]) < 0:
        u = -u
    u = u - (u @ n_fp) * n_fp
    u /= np.linalg.norm(u)

    # The radar is right-looking: the look direction is Right = u x Up (Up
    # the scene radial).  The two self-consistent sub-satellite placements
    # are mirror images about the track plane (both subswath assignments give
    # the same nadir-to-chain distances); orienting the normal opposite to
    # the look (n_fp . Right < 0) picks the physical one, in which
    # convention the near subswath WV1 is the larger-offset chain and the
    # sub-satellite point lies beyond it (delta0 = d_m + d_km).
    up = mean_center / np.linalg.norm(mean_center)
    right = np.cross(u, up)
    right /= np.linalg.norm(right)
    if n_fp @ right > 0:
        n_fp = -n_fp
    u = u - (u @ n_fp) * n_fp
    u /= np.linalg.norm(u)
    rperp = np.cross(n_fp, u)
    rperp /= np.linalg.norm(rperp)

    # Label the two subswath chains from the gap in the track-plane offset.
    offset = centers @ n_fp
    srt = np.sort(offset)
    gaps = np.diff(srt)
    k = int(np.argmax(gaps))
    if gaps[k] < 10.0:
        raise ValueError(
            "WV footprint does not show two distinct subswath chains "
            f"(largest track-plane offset gap {gaps[k]:.1f} km)")
    thr = 0.5 * (srt[k] + srt[k + 1])
    is_wv1 = offset >= thr  # WV1 = larger offset (closer to the nadir)
    if not is_wv1.any() or is_wv1.all():
        raise ValueError("WV subswath chain labeling failed (single chain)")
    d1m = float(offset[is_wv1].mean())
    d2m = float(offset[~is_wv1].mean())

    a = centers @ u
    x = centers @ rperp
    phi = np.arctan2(x, a)
    dphi = np.array([
        (phi[i + 2] - phi[i] + np.pi) % (2.0 * np.pi) - np.pi
        for i in range(0, n_ring - 2, 2)
    ])
    sense = 1.0 if dphi.mean() > 0 else -1.0

    # Per-burst sub-satellite point Q_i and heading.
    h = r_s - R_E_KM
    q = np.empty((n_ring, 3))
    burst_heading = np.empty(n_ring)
    for i in range(n_ring):
        sw = 1 if is_wv1[i] else 2
        th = math.radians(wv_theta_deg(sw, r_s, mean_lat))
        d_km = h * math.tan(th)
        d_m = d1m if sw == 1 else d2m
        delta0 = d_m + d_km
        if delta0 >= R_E_KM:
            raise ValueError(
                f"degenerate WV geometry at burst {i} "
                f"(delta0 = {delta0:.1f} km >= R_E = {R_E_KM:.1f} km)")
        rho_t = math.sqrt(R_E_KM ** 2 - delta0 ** 2)
        qi = (delta0 * n_fp
              + rho_t * (math.cos(phi[i]) * u + math.sin(phi[i]) * rperp))
        q[i] = R_E_KM * qi / np.linalg.norm(qi)
        t = sense * (-math.sin(phi[i]) * u + math.cos(phi[i]) * rperp)
        slon = math.degrees(math.atan2(q[i, 1], q[i, 0]))
        slat = math.degrees(math.asin(np.clip(q[i, 2] / np.linalg.norm(q[i]),
                                             -1.0, 1.0)))
        e, nn, uu = _enu_basis(slon, slat)
        th2 = t - (t @ uu) * uu
        th2 /= np.linalg.norm(th2)
        burst_heading[i] = math.degrees(math.atan2(th2 @ e, th2 @ nn)) % 360.0

    # Target points -> nearest burst ring, then per-point angles.
    lon_arr = np.atleast_1d(np.asarray(lon, dtype=float))
    lat_arr = np.atleast_1d(np.asarray(lat, dtype=float))
    if lon_arr.shape != lat_arr.shape:
        raise ValueError(f"lon and lat must have the same shape, got "
                         f"{lon_arr.shape} and {lat_arr.shape}")
    p = _to_3d(np.column_stack([lon_arr.ravel(), lat_arr.ravel()]))
    p = R_E_KM * p / np.linalg.norm(p, axis=1, keepdims=True)
    d = np.linalg.norm(p[:, None, :] - centers[None, :, :], axis=2)
    ring_index = np.argmin(d, axis=1)

    inc = np.empty(p.shape[0])
    el = np.empty(p.shape[0])
    for j in range(p.shape[0]):
        qi = q[ring_index[j]]
        gamma = math.acos(np.clip((qi @ p[j]) / R_E_KM ** 2, -1.0, 1.0))
        rk = math.sqrt(r_s ** 2 + R_E_KM ** 2
                       - 2.0 * r_s * R_E_KM * math.cos(gamma))
        ci = (r_s ** 2 - rk ** 2 - R_E_KM ** 2) / (2.0 * rk * R_E_KM)
        inc[j] = math.degrees(math.acos(np.clip(ci, -1.0, 1.0)))
        el[j] = inc[j] - math.degrees(gamma)
    heading = burst_heading[ring_index]
    subswath = np.where(is_wv1[ring_index], 1, 2)

    return S1WVAngles(
        mode="WV",
        r_s_km=r_s,
        mean_lat_deg=mean_lat,
        n_rings=n_ring,
        lon_deg=lon_arr,
        lat_deg=lat_arr,
        incidence_deg=inc,
        elevation_deg=el,
        heading_deg=heading,
        ring_index=ring_index,
        subswath=subswath,
    )
