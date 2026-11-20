"""S1 **WV** SAR incidence, elevation and heading angles from the SAFE footprint.

The Wave mode (WV) of the Sentinel-1 C-band radar acquires small (~20 km)
imagettes on two fixed off-nadir subswaths (WV1 ~ 21 deg, WV2 ~ 32 deg
off-nadir), alternating between the two subswaths ("leapfrog") along the
ground track.  The product footprint (the CDSE OData ``Footprint`` /
``GeoFootprint``) is a *MultiPolygon*: one small ring per imagette, in time
order, the rings forming two interleaved chains (the subswaths).  A single
imagette polygon does not carry enough information: the chains are needed to
locate the track plane.

Unlike IW / EW, the WV footprint does not encode a single constant-slant-range
swath, so the satellite geocentric distance ``r_s`` cannot be inverted from
it.  It is either given, or (by default) taken from the metadata-only orbit
model :func:`~s1_geometry.orbit.satellite_radius_km` (latitude and pass
direction, ~0.1 km accuracy), so the chain runs from the OData metadata
alone.  The angles are returned **per ground point**.  The model:

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
elevation, heading), from the footprint metadata only (``r_s`` from the
orbit model): 63 Sentinel-1C products (2025-01 to 2026-11, latitudes -74 to
+52 deg, up to 118 imagettes / 29 min), 1.43 million cells: RMSE 0.24 deg
(incidence), 0.23 deg (elevation), 0.09 deg (heading), subswath labels 100 %
correct; and two Sentinel-1D products.  See :doc:`/certification`.

Limits: a product needs at least 4 imagettes (end-of-datatake products with
2 imagettes are rejected); the rings must be in time order; the orbit model
assumes the operational reference orbit (not the commissioning phase).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from s1_geometry.angles import R_E_KM
from s1_geometry.footprint import (
    ASCENDING,
    DESCENDING,
    as_multipolygon,
    merge_antimeridian_splits,
)
from s1_geometry.orbit import R_S_ORBIT_MAX_KM, R_S_ORBIT_MIN_KM, satellite_radius_km

__all__ = [
    "WV_AUTO_WINDOW",
    "WV_AUTO_WINDOW_MIN_RINGS",
    "WV_HEADING_UNCERTAINTY_DEG",
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

#: ``window="auto"``: products with more imagettes than this are fitted with a
#: sliding window of +/- ``WV_AUTO_WINDOW`` imagettes (S1C validation).
WV_AUTO_WINDOW_MIN_RINGS = 40
WV_AUTO_WINDOW = 10

#: WV heading uncertainty (deg): the 95th percentile of the heading error of
#: the worst product of the S1C validation (0.234 deg), rounded up.
WV_HEADING_UNCERTAINTY_DEG = 0.25

#: Calibrated latitude term ``(deg/deg, deg)`` of the look-angle law, fitted on
#: 64 WV OCN products (2026-01..06); residual std 0.086 deg (WV1) / 0.065 deg
#: (WV2).  Do not change without refitting and revalidating.
_WV_THETA_LAT = {1: (0.012076, 0.2450), 2: (0.009995, 0.6560)}


def wv_theta_deg(subswath: int, r_s_km: float, lat_deg: float) -> float:
    """WV look (off-nadir) angle (deg) for subswath ``subswath`` (1 or 2).

    The documented band-center at the reference altitude (711.7 km), linear in
    the satellite geocentric distance ``r_s_km``, plus a calibrated latitude
    term.  ``r_s_km`` is the satellite geocentric distance (km).

    Convention: the law uses ``r_s_km - R_E_KM`` (6371 km sphere) as the
    "altitude", whereas the documented band (697.6 to 725.8 km) is the height
    above the WGS84 ellipsoid; the two differ by up to ~14 km at high
    latitude.  The calibrated latitude term was fitted with this convention,
    so keep it unchanged unless the law is refitted.
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


def _closed(ring: np.ndarray) -> np.ndarray:
    """The ring with its first vertex repeated at the end.

    The ring centre is the mean of the closed ring's vertices (the first
    vertex counted twice): this is the convention of the CDSE footprints the
    look-angle law was calibrated with (the unweighted vertex mean gives a
    0.33 deg instead of 0.24 deg incidence RMSE on the S1C validation), so it
    is enforced whatever the input form.
    """
    if np.allclose(ring[0], ring[-1]):
        return ring
    return np.vstack([ring, ring[:1]])


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
        Mean satellite geocentric distance (km) over the imagette rings.
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


def _pass_direction_per_ring(centers: np.ndarray) -> np.ndarray:
    """Pass direction of each ring from the ring time progression.

    The latitude change between the same-subswath neighbours ``i - 2`` and
    ``i + 2`` (the ring order is the burst / time order) gives the flight's
    north component; a product crossing the polar turn gets the right
    direction on each side.
    """
    lat = np.degrees(np.arcsin(np.clip(
        centers[:, 2] / np.linalg.norm(centers, axis=1), -1.0, 1.0)))
    n = len(lat)
    idx = np.arange(n)
    dlat = lat[np.minimum(idx + 2, n - 1)] - lat[np.maximum(idx - 2, 0)]
    return np.where(dlat >= 0.0, ASCENDING, DESCENDING)


def _solve_window(centers: np.ndarray, r_s: np.ndarray, mean_lat: float):
    """Sub-satellite point, heading and subswath of each ring of a window.

    ``centers`` are the (n, 3) ring centres (km, ECEF, time order) and
    ``r_s`` the (n,) satellite geocentric distances (km).  Returns
    ``(q, heading_deg, is_wv1)``: the (n, 3) sub-satellite points, the (n,)
    in-plane tangent headings and the (n,) WV1 flags.
    """
    n_ring = centers.shape[0]
    mean_center = centers.mean(axis=0)

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
            f"(largest track-plane offset gap {gaps[k]:.1f} km); are the rings in "
            "time order (the CDSE Footprint / GeoFootprint order)?")
    thr = 0.5 * (srt[k] + srt[k + 1])
    is_wv1 = offset >= thr  # WV1 = larger offset (closer to the nadir)
    if not is_wv1.any() or is_wv1.all():
        raise ValueError("WV subswath chain labeling failed (single chain)")
    # In time order the imagettes alternate between WV1 and WV2 (checked on
    # every CDSE WV product); a footprint whose rings are not in time order
    # (e.g. re-ordered by a GIS union / normalisation) cannot be used.
    if n_ring >= 4 and np.mean(np.diff(is_wv1.astype(int)) != 0) < 0.9:
        raise ValueError(
            "WV footprint rings are not in time order (consecutive imagettes do not "
            "alternate between WV1 and WV2): use the CDSE Footprint / GeoFootprint, "
            "not a re-ordered geometry")
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
    q = np.empty((n_ring, 3))
    burst_heading = np.empty(n_ring)
    for i in range(n_ring):
        sw = 1 if is_wv1[i] else 2
        h = r_s[i] - R_E_KM
        th = math.radians(wv_theta_deg(sw, r_s[i], mean_lat))
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
    return q, burst_heading, is_wv1


def estimate_wv_angles(
    sar_multipolygon,
    r_s_km: float | None,
    lon,
    lat,
    *,
    orbit_direction: str | None = None,
    window: int | str | None = "auto",
    on_outside: str = "raise",
) -> S1WVAngles:
    """Estimate the WV incidence / elevation / heading angles at target points.

    Parameters
    ----------
    sar_multipolygon:
        The WV SAFE footprint: a GeoJSON MultiPolygon string / dict, a WKT
        ``MULTIPOLYGON`` (e.g. the CDSE OData ``Footprint``), or a sequence of
        (lon, lat) rings (one small ring per imagette, **in time order**,
        alternating between the two subswaths).  Imagettes split at the
        antimeridian (two polygons in the OData ``Footprint``) are re-joined,
        so ``ring_index`` is the imagette index (0-based, time order).
    r_s_km:
        Satellite geocentric distance (km), or ``None`` to derive it per
        ring from the metadata-only orbit model
        :func:`~s1_geometry.orbit.satellite_radius_km` (latitude + pass
        direction).  When given, it must lie in the Sentinel-1 orbit range
        :data:`~s1_geometry.orbit.R_S_ORBIT_MIN_KM` to
        :data:`~s1_geometry.orbit.R_S_ORBIT_MAX_KM`.
    lon, lat:
        Target point longitudes / latitudes (deg): a scalar or an array (e.g.
        the centres of the ESA OCN ``owi`` grid cells, or the imagette
        reference points).  The angles are returned per point.
    orbit_direction:
        Pass direction used by the orbit model when ``r_s_km`` is ``None``
        (``"ascending"`` / ``"descending"``, case-insensitive, e.g. the CDSE
        ``orbitDirection``).  By default it is derived per ring from the ring
        time progression, which also handles products crossing the polar
        turn (the catalogue attribute only describes the product start).
    window:
        ``"auto"`` (default): one track geometry for products of up to
        :data:`WV_AUTO_WINDOW_MIN_RINGS` imagettes, else a sliding window of
        +/- :data:`WV_AUTO_WINDOW` imagettes.  ``None``: one track geometry
        for the whole product.  An integer ``w``: the geometry of each ring
        is fitted on the ``2 w + 1`` rings around it.  Long products (up to
        ~30 min, 120 imagettes) have a ground track that departs from a
        single great circle (Earth rotation): on the S1C validation the
        sliding window halves their error (incidence RMSE 0.22 instead of
        0.45 deg, heading 0.09 instead of 0.36 deg).
    on_outside:
        ``"raise"`` (default) raises ``ValueError`` when a target point lies
        outside every imagette ring; ``"nan"`` returns NaN angles for it.

    Returns
    -------
    S1WVAngles
        Per-point incidence / elevation / heading and the intermediate
        quantities (see :class:`S1WVAngles`).

    Raises
    ------
    ValueError
        If ``r_s_km`` is outside the Sentinel-1 orbit range, the footprint is
        malformed or does not show two distinct subswath chains, or (with
        ``on_outside="raise"``) a target point is outside the footprint.
    """
    if on_outside not in ("raise", "nan"):
        raise ValueError(f"on_outside must be 'raise' or 'nan', got {on_outside!r}")
    rings = [_closed(r) for r in merge_antimeridian_splits(as_multipolygon(sar_multipolygon))]
    if len(rings) < 4:
        raise ValueError(
            f"WV footprint needs at least 4 imagette rings (two per subswath), "
            f"got {len(rings)}; a single imagette polygon is not enough, pass "
            "the product footprint (the CDSE OData Footprint MultiPolygon)")
    centers = np.array([_to_3d(r).mean(axis=0) for r in rings])
    n_ring = centers.shape[0]
    ring_lat = np.array([float(r[:, 1].mean()) for r in rings])
    mean_lat = float(ring_lat.mean())

    if r_s_km is None:
        if orbit_direction is None:
            directions = _pass_direction_per_ring(centers)
        else:
            directions = np.full(n_ring, str(orbit_direction).lower())
        r_s = np.array([satellite_radius_km(lt, d)
                        for lt, d in zip(ring_lat, directions, strict=True)])
    else:
        r_s = np.full(n_ring, float(r_s_km))
        if not (R_S_ORBIT_MIN_KM <= r_s[0] <= R_S_ORBIT_MAX_KM):
            raise ValueError(
                f"r_s_km {r_s[0]:.2f} is outside the Sentinel-1 orbit range "
                f"{R_S_ORBIT_MIN_KM:.1f} to {R_S_ORBIT_MAX_KM:.1f} km (satellite "
                "geocentric distance, not the altitude above the ellipsoid)")

    if window == "auto":
        window = WV_AUTO_WINDOW if n_ring > WV_AUTO_WINDOW_MIN_RINGS else None
    if window is None or n_ring <= 2 * int(window) + 1:
        q, burst_heading, is_wv1 = _solve_window(centers, r_s, mean_lat)
    else:
        w = int(window)
        if w < 2:
            raise ValueError(f"window must be >= 2 rings, got {window!r}")
        size = 2 * w + 1
        q = np.empty((n_ring, 3))
        burst_heading = np.empty(n_ring)
        is_wv1 = np.empty(n_ring, dtype=bool)
        for i in range(n_ring):
            lo = min(max(i - w, 0), n_ring - size)
            sub = slice(lo, lo + size)
            qw, hw, sw = _solve_window(centers[sub], r_s[sub],
                                       float(ring_lat[sub].mean()))
            q[i], burst_heading[i], is_wv1[i] = qw[i - lo], hw[i - lo], sw[i - lo]

    # Target points -> nearest imagette ring, then per-point angles.
    lon_arr = np.atleast_1d(np.asarray(lon, dtype=float))
    lat_arr = np.atleast_1d(np.asarray(lat, dtype=float))
    if lon_arr.shape != lat_arr.shape:
        raise ValueError(f"lon and lat must have the same shape, got "
                         f"{lon_arr.shape} and {lat_arr.shape}")
    p = _to_3d(np.column_stack([lon_arr.ravel(), lat_arr.ravel()]))
    d = np.linalg.norm(p[:, None, :] - centers[None, :, :], axis=2)
    ring_index = np.argmin(d, axis=1)

    # A point is inside the footprint when it is within its nearest ring's
    # circumscribed radius (centre to farthest vertex), with a 10 % margin.
    ring_radius = np.array([
        np.linalg.norm(_to_3d(r) - c, axis=1).max() for r, c in zip(rings, centers, strict=True)
    ])
    outside = d[np.arange(len(p)), ring_index] > 1.1 * ring_radius[ring_index]
    if outside.any() and on_outside == "raise":
        j = int(np.flatnonzero(outside)[0])
        raise ValueError(
            f"{int(outside.sum())} target point(s) outside the WV footprint, e.g. "
            f"({lon_arr.ravel()[j]:.4f}, {lat_arr.ravel()[j]:.4f}), "
            f"{d[j, ring_index[j]]:.1f} km from the nearest imagette centre")

    qi = q[ring_index]
    rs = r_s[ring_index]
    gamma = np.arccos(np.clip(np.einsum("ij,ij->i", qi, p) / R_E_KM ** 2, -1.0, 1.0))
    rk = np.sqrt(rs ** 2 + R_E_KM ** 2 - 2.0 * rs * R_E_KM * np.cos(gamma))
    ci = (rs ** 2 - rk ** 2 - R_E_KM ** 2) / (2.0 * rk * R_E_KM)
    inc = np.degrees(np.arccos(np.clip(ci, -1.0, 1.0)))
    el = inc - np.degrees(gamma)
    heading = burst_heading[ring_index].copy()
    subswath = np.where(is_wv1[ring_index], 1, 2)
    if outside.any():
        inc[outside] = el[outside] = heading[outside] = np.nan

    return S1WVAngles(
        mode="WV",
        r_s_km=float(r_s.mean()),
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
