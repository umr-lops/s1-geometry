"""Tests for the WV (Wave mode) angle estimation (s1_geometry.wv).

The synthetic footprint generator places a satellite on a circular orbit at
geocentric distance ``r_s`` looking right at two fixed off-nadir angles and
records the ground-range point of each burst as a small ring.  The estimator
must invert that footprint (footprint + ``r_s`` only) back to the per-point
incidence / elevation / heading, which it is checked against the physical
ground truth.
"""

import json

import numpy as np
import pytest

from s1_geometry import (
    R_E_KM,
    R_S_MAX_KM,
    R_S_MIN_KM,
    S1WVAngles,
    as_multipolygon,
    estimate_wv_angles,
    merge_antimeridian_splits,
    wv_theta_deg,
)

_REF_ALT_KM = 711.7


# --------------------------------------------------------------------------- #
# Synthetic WV footprint (physical ground truth)
# --------------------------------------------------------------------------- #
def _ecef_to_lonlat(p):
    lon = np.degrees(np.arctan2(p[1], p[0]))
    lat = np.degrees(np.arcsin(np.clip(p[2] / np.linalg.norm(p), -1.0, 1.0)))
    return lon, lat


def _ground_range(r_s, beta, theta_deg):
    """Ground-range point (ECEF, on the R_E sphere) for a satellite at central
    angle ``beta`` (along track from the scene centre) looking right at
    off-nadir ``theta_deg``."""
    sat = r_s * np.array([np.cos(beta), 0.0, np.sin(beta)])
    up = sat / r_s
    fwd = np.array([-np.sin(beta), 0.0, np.cos(beta)])
    right = np.cross(fwd, up)
    right /= np.linalg.norm(right)
    beam = (np.cos(np.radians(theta_deg)) * (-up)
            + np.sin(np.radians(theta_deg)) * right)
    # slant range: |sat + t beam| = R_E.  The beam points toward Earth, so the
    # imaged point is the NEAR-side intersection (the smaller root).
    t = -(sat @ beam) - np.sqrt((sat @ beam) ** 2 - (r_s ** 2 - R_E_KM ** 2))
    g = sat + t * beam
    return R_E_KM * g / np.linalg.norm(g), sat


def _azimuth_clockwise_from_north(d, p):
    """Azimuth (deg, clockwise from North) of the ECEF vector ``d`` at point
    ``p`` (ECEF, on the sphere)."""
    up = p / R_E_KM
    east = np.cross(np.array([0.0, 0.0, 1.0]), up)
    if np.linalg.norm(east) < 1e-9:
        east = np.array([1.0, 0.0, 0.0])
    east /= np.linalg.norm(east)
    north = np.cross(up, east)
    return np.degrees(np.arctan2(d @ east, d @ north)) % 360.0


def _synthetic_wv(r_s_km=7085.0, n_bursts=20, theta1_deg=21.0, theta2_deg=33.0,
                 along_track_km=150.0, ring_radius_km=7.0, n_ring_vertices=8):
    """Build a physically-consistent WV MultiPolygon.

    Returns ``(multipolygon, subswaths, truth)`` where ``multipolygon`` is a
    list of (lon, lat) rings (burst / time order, alternating subswath),
    ``subswaths`` is an (n_bursts,) array of 1 / 2, and ``truth`` is a dict of
    per-burst ground-truth arrays (incidence, elevation, heading, in deg).
    """
    r_s = float(r_s_km)
    rings = []
    subswaths = []
    inc_true = []
    el_true = []
    hd_true = []
    for i in range(n_bursts):
        beta = (i * along_track_km) / R_E_KM          # central angle along track
        theta = theta1_deg if i % 2 == 0 else theta2_deg
        sw = 1 if i % 2 == 0 else 2
        g, sat = _ground_range(r_s, beta, theta)
        # ground truth angles at the ground-range point
        up_g = g / R_E_KM
        to_sat = sat - g
        inc = np.degrees(np.arccos(np.clip(to_sat @ up_g / np.linalg.norm(to_sat),
                                          -1.0, 1.0)))
        gamma = np.degrees(np.arccos(np.clip((sat / r_s) @ up_g, -1.0, 1.0)))
        inc_true.append(inc)
        el_true.append(inc - gamma)
        # heading: instantaneous tangent of the ground-range motion (a small
        # forward step along the track)
        g_fwd, _ = _ground_range(r_s, beta + 0.005, theta)
        hd_true.append(_azimuth_clockwise_from_north(g_fwd - g, g))
        # small ring around the ground-range point (local tangent circle)
        east = np.cross(np.array([0.0, 0.0, 1.0]), up_g)
        if np.linalg.norm(east) < 1e-9:
            east = np.array([1.0, 0.0, 0.0])
        east /= np.linalg.norm(east)
        north = np.cross(up_g, east)
        verts = []
        for k in range(n_ring_vertices):
            a = 2.0 * np.pi * k / n_ring_vertices
            p = g + ring_radius_km * (np.cos(a) * east + np.sin(a) * north)
            p = R_E_KM * p / np.linalg.norm(p)
            verts.append(_ecef_to_lonlat(p))
        verts.append(verts[0])
        rings.append(verts)
        subswaths.append(sw)
    return rings, np.array(subswaths), dict(
        inc=np.array(inc_true), el=np.array(el_true), hd=np.array(hd_true))


def _lonlat_to_ecef(v):
    lon, lat = np.radians(v[0]), np.radians(v[1])
    return R_E_KM * np.array([np.cos(lat) * np.cos(lon),
                              np.cos(lat) * np.sin(lon),
                              np.sin(lat)])


def _daz(a, b):
    return (a - b + 180.0) % 360.0 - 180.0


# --------------------------------------------------------------------------- #
# Look-angle law
# --------------------------------------------------------------------------- #
def test_wv_theta_deg_documented_band_centers():
    # at the reference altitude the law returns the documented off-nadir band
    # centers (WV1 ~ 21 deg, WV2 ~ 33 deg), plus the calibrated latitude term
    r_ref = R_E_KM + _REF_ALT_KM
    th1_ref = wv_theta_deg(1, r_ref, 0.0)
    th2_ref = wv_theta_deg(2, r_ref, 0.0)
    assert 20.0 < th1_ref < 23.0
    assert 31.0 < th2_ref < 35.0
    # WV2 always looks farther from nadir than WV1
    assert th2_ref > th1_ref


def test_wv_theta_deg_decreases_with_altitude():
    # higher orbit (larger R_s) -> smaller off-nadir angle (documented law)
    assert wv_theta_deg(1, R_S_MAX_KM, 0.0) < wv_theta_deg(1, R_S_MIN_KM, 0.0)
    assert wv_theta_deg(2, R_S_MAX_KM, 0.0) < wv_theta_deg(2, R_S_MIN_KM, 0.0)


def test_wv_theta_deg_rejects_unknown_subswath():
    with pytest.raises(ValueError, match="unknown WV subswath"):
        wv_theta_deg(3, R_E_KM + _REF_ALT_KM, 0.0)


# --------------------------------------------------------------------------- #
# MultiPolygon parsing
# --------------------------------------------------------------------------- #
def test_as_multipolygon_from_list_of_rings():
    rings, _, _ = _synthetic_wv(n_bursts=8)
    out = as_multipolygon(rings)
    assert len(out) == 8
    assert out[0].shape == (9, 2)


def test_as_multipolygon_from_geojson_string():
    rings, _, _ = _synthetic_wv(n_bursts=8)
    geo = json.dumps({"type": "MultiPolygon",
                      "coordinates": [[r] for r in rings]})
    out = as_multipolygon(geo)
    assert len(out) == 8


def test_as_multipolygon_rejects_single_polygon():
    box = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0)]
    with pytest.raises(ValueError, match="MultiPolygon"):
        as_multipolygon({"type": "Polygon", "coordinates": [box]})


def test_as_multipolygon_rejects_empty():
    with pytest.raises(ValueError, match="no sub-polygons"):
        as_multipolygon([])


# --------------------------------------------------------------------------- #
# Full chain: footprint + r_s -> per-point angles
# --------------------------------------------------------------------------- #
def test_estimate_wv_angles_inverts_synthetic_footprint():
    rings, subswaths, truth = _synthetic_wv(r_s_km=7085.0, n_bursts=20)
    # target points = the ring centres (one per burst)
    lons = np.array([np.mean([v[0] for v in r]) for r in rings])
    lats = np.array([np.mean([v[1] for v in r]) for r in rings])
    result = estimate_wv_angles(rings, 7085.0, lons, lats)
    assert isinstance(result, S1WVAngles)
    assert result.mode == "WV"
    assert result.n_rings == 20
    assert result.incidence_deg.shape == (20,)
    assert result.elevation_deg.shape == (20,)
    assert result.heading_deg.shape == (20,)
    assert not (np.isnan(result.incidence_deg).any()
                or np.isnan(result.elevation_deg).any())
    # the estimator recovers the physical ground truth (footprint + r_s only)
    assert np.sqrt((result.incidence_deg - truth["inc"]) ** 2).mean() < 0.5
    assert np.sqrt((result.elevation_deg - truth["el"]) ** 2).mean() < 0.5
    assert np.sqrt(_daz(result.heading_deg, truth["hd"]) ** 2).mean() < 1.0
    # the two subswaths are correctly separated, WV1 nearer nadir than WV2
    inc1 = result.incidence_deg[result.subswath == 1]
    inc2 = result.incidence_deg[result.subswath == 2]
    assert inc1.mean() < inc2.mean()
    assert 15.0 < inc1.mean() < 30.0
    assert 25.0 < inc2.mean() < 45.0


def test_estimate_wv_angles_subswath_labels_match_generator():
    rings, subswaths, _ = _synthetic_wv(r_s_km=7085.0, n_bursts=20)
    lons = np.array([np.mean([v[0] for v in r]) for r in rings])
    lats = np.array([np.mean([v[1] for v in r]) for r in rings])
    result = estimate_wv_angles(rings, 7085.0, lons, lats)
    # the footprint-only chain labeling must agree with the generator's
    # subswath assignment (up to the WV1 / WV2 naming, which is fixed by the
    # look direction, so it should match exactly here)
    assert np.array_equal(result.subswath, subswaths)


def test_estimate_wv_angles_accepts_geojson_string():
    rings, _, _ = _synthetic_wv(r_s_km=7085.0, n_bursts=12)
    geo = json.dumps({"type": "MultiPolygon",
                      "coordinates": [[r] for r in rings]})
    lons = np.array([np.mean([v[0] for v in r]) for r in rings])
    lats = np.array([np.mean([v[1] for v in r]) for r in rings])
    result = estimate_wv_angles(geo, 7085.0, lons, lats)
    assert result.n_rings == 12


# --------------------------------------------------------------------------- #
# Error cases
# --------------------------------------------------------------------------- #
def test_estimate_wv_angles_rejects_out_of_band_r_s():
    rings, _, _ = _synthetic_wv(n_bursts=8)
    lons = np.array([np.mean([v[0] for v in r]) for r in rings])
    lats = np.array([np.mean([v[1] for v in r]) for r in rings])
    with pytest.raises(ValueError, match="orbit range"):
        estimate_wv_angles(rings, 5000.0, lons, lats)


def test_estimate_wv_angles_rejects_too_few_rings():
    rings, _, _ = _synthetic_wv(n_bursts=3)
    lons = np.array([np.mean([v[0] for v in r]) for r in rings])
    lats = np.array([np.mean([v[1] for v in r]) for r in rings])
    with pytest.raises(ValueError, match="at least 4 imagette rings"):
        estimate_wv_angles(rings, 7085.0, lons, lats)


def _centres(rings):
    lons = np.array([np.mean([v[0] for v in r[:-1]]) for r in rings])
    lats = np.array([np.mean([v[1] for v in r[:-1]]) for r in rings])
    return lons, lats


def test_estimate_wv_angles_accepts_high_latitude_r_s():
    # 7066 km is a real Sentinel-1 geocentric distance (polar part of the
    # orbit); the former band (6371 + 697.6 km) wrongly rejected it
    rings, _, _ = _synthetic_wv(r_s_km=7066.0, n_bursts=8)
    lons, lats = _centres(rings)
    result = estimate_wv_angles(rings, 7066.0, lons, lats)
    assert np.isfinite(result.incidence_deg).all()


def test_estimate_wv_angles_r_s_from_orbit_model():
    # metadata-only chain: no r_s given, the orbit model provides it
    # short product (8 imagettes, 400 km): r_s barely varies along it
    rings, _, truth = _synthetic_wv(r_s_km=7077.0, n_bursts=8, along_track_km=50.0)
    lons, lats = _centres(rings)
    result = estimate_wv_angles(rings, None, lons, lats)
    assert R_S_MIN_KM - 10.0 < result.r_s_km < R_S_MAX_KM
    explicit = estimate_wv_angles(rings, result.r_s_km, lons, lats)
    np.testing.assert_allclose(result.incidence_deg, explicit.incidence_deg, atol=0.05)
    # same order of accuracy as with the true r_s
    assert np.abs(result.incidence_deg - truth["inc"]).mean() < 1.0
    # an explicit orbit direction is accepted too
    asc = estimate_wv_angles(rings, None, lons, lats, orbit_direction="ASCENDING")
    assert np.isfinite(asc.incidence_deg).all()


def test_estimate_wv_angles_window_matches_full_fit_on_great_circle():
    # on a great-circle track the sliding-window fit equals the global fit
    rings, subswaths, _ = _synthetic_wv(r_s_km=7077.0, n_bursts=30)
    lons, lats = _centres(rings)
    full = estimate_wv_angles(rings, 7077.0, lons, lats)
    win = estimate_wv_angles(rings, 7077.0, lons, lats, window=4)
    # same track geometry, and the look-angle law uses each imagette's own
    # latitude in both cases
    np.testing.assert_allclose(win.incidence_deg, full.incidence_deg, atol=0.02)
    assert np.abs(_daz(win.heading_deg, full.heading_deg)).max() < 0.1
    assert np.array_equal(win.subswath, subswaths)
    with pytest.raises(ValueError, match="window"):
        estimate_wv_angles(rings, 7077.0, lons, lats, window=1)


def test_estimate_wv_angles_outside_point():
    rings, _, _ = _synthetic_wv(r_s_km=7077.0, n_bursts=8)
    lons, lats = _centres(rings)
    far_lon, far_lat = np.append(lons, lons[0] + 30.0), np.append(lats, lats[0])
    with pytest.raises(ValueError, match="outside the WV footprint"):
        estimate_wv_angles(rings, 7077.0, far_lon, far_lat)
    result = estimate_wv_angles(rings, 7077.0, far_lon, far_lat, on_outside="nan")
    assert np.isnan(result.incidence_deg[-1]) and np.isnan(result.heading_deg[-1])
    assert np.isfinite(result.incidence_deg[:-1]).all()
    with pytest.raises(ValueError, match="on_outside"):
        estimate_wv_angles(rings, 7077.0, lons, lats, on_outside="skip")


def test_estimate_wv_angles_accepts_odata_wkt():
    rings, _, _ = _synthetic_wv(r_s_km=7077.0, n_bursts=8)
    lons, lats = _centres(rings)
    wkt = "geography'SRID=4326;MULTIPOLYGON (" + ", ".join(
        "((" + ", ".join(f"{lo} {la}" for lo, la in r) + "))" for r in rings) + ")'"
    from_wkt = estimate_wv_angles(wkt, 7077.0, lons, lats)
    from_list = estimate_wv_angles(rings, 7077.0, lons, lats)
    np.testing.assert_allclose(from_wkt.incidence_deg, from_list.incidence_deg)


def test_single_imagette_is_rejected_with_a_hint():
    rings, _, _ = _synthetic_wv(n_bursts=8)
    with pytest.raises(ValueError, match="single imagette"):
        estimate_wv_angles(rings[:1], 7077.0, 0.0, 0.0)


def _split_at_antimeridian(ring):
    """Cut a ring (unwrapped lon around 180) into the two OData-style pieces."""
    pts = np.asarray(ring)[:-1]
    pieces = []
    for side in (1.0, -1.0):  # +180 side keeps lon <= 180, -180 side lon >= 180
        keep = (lambda p: p[0] <= 180.0) if side > 0 else (lambda p: p[0] >= 180.0)
        out = []
        for k in range(len(pts)):
            p, q = pts[k], pts[(k + 1) % len(pts)]
            if keep(p):
                out.append(p)
            if keep(p) != keep(q):
                t = (180.0 - p[0]) / (q[0] - p[0])
                out.append(np.array([180.0, p[1] + t * (q[1] - p[1])]))
        out = np.array(out)
        out[:, 0] = np.where(out[:, 0] > 180.0 - 1e-12, out[:, 0] - 360.0, out[:, 0]) \
            if side < 0 else out[:, 0]
        pieces.append(np.vstack([out, out[:1]]))
    # time order of the pieces does not matter; -180 piece first like the CDSE
    return [pieces[1], pieces[0]]


def test_antimeridian_split_imagette_is_rejoined():
    rings, subswaths, _ = _synthetic_wv(r_s_km=7077.0, n_bursts=10, along_track_km=60.0)
    lons, lats = _centres(rings)
    k = 5
    shift = 180.013 - lons[k]  # ring k straddles the antimeridian (no corner on it)
    shifted = [[(lo + shift, la) for lo, la in r] for r in rings]
    wrap = [[((lo + 180.0) % 360.0 - 180.0, la) for lo, la in r] for r in shifted]
    pieces = _split_at_antimeridian(shifted[k])
    split = wrap[:k] + pieces + wrap[k + 1:]
    assert len(split) == len(rings) + 1
    merged = merge_antimeridian_splits(as_multipolygon(split))
    assert len(merged) == len(rings)
    # the merged ring holds exactly the original corners (unwrapped)
    orig = np.asarray(shifted[k])[:-1] - [360.0, 0.0]
    got = merged[k][:-1]
    assert len(got) == len(orig)
    assert np.abs(np.sort(got, axis=0) - np.sort(orig, axis=0)).max() < 1e-6
    # end to end: same result as the unsplit footprint with that merged ring
    plon = (lons + shift + 180.0) % 360.0 - 180.0
    res = estimate_wv_angles(split, 7077.0, plon, lats)
    ref = estimate_wv_angles(wrap[:k] + [merged[k]] + wrap[k + 1:], 7077.0, plon, lats)
    assert res.n_rings == len(rings)
    np.testing.assert_array_equal(res.subswath, subswaths)
    np.testing.assert_array_equal(res.ring_index, np.arange(len(rings)))
    np.testing.assert_allclose(res.incidence_deg, ref.incidence_deg, atol=1e-9)


def test_open_and_closed_rings_give_the_same_result():
    rings, _, _ = _synthetic_wv(r_s_km=7077.0, n_bursts=8, along_track_km=60.0)
    lons, lats = _centres(rings)
    closed = estimate_wv_angles(rings, 7077.0, lons, lats)
    opened = estimate_wv_angles([r[:-1] for r in rings], 7077.0, lons, lats)
    np.testing.assert_allclose(opened.incidence_deg, closed.incidence_deg)


def test_rings_not_in_time_order_are_rejected():
    rings, _, _ = _synthetic_wv(r_s_km=7077.0, n_bursts=12, along_track_km=60.0)
    lons, lats = _centres(rings)
    order = list(range(0, 12, 2)) + list(range(1, 12, 2))  # WV1 chain then WV2 chain
    with pytest.raises(ValueError, match="time order"):
        estimate_wv_angles([rings[i] for i in order], 7077.0, lons, lats)


def test_wv_footprint_from_imagettes_orders_and_keeps_dateline_pieces():
    from s1_geometry import wv_footprint_from_imagettes

    rings, subswaths, _ = _synthetic_wv(r_s_km=7077.0, n_bursts=8, along_track_km=60.0)
    lons, lats = _centres(rings)
    imagettes = [{"index": k + 1, "GeoFootprint": {"type": "Polygon", "coordinates": [r]}}
                 for k, r in enumerate(rings)]
    shuffled = imagettes[::-1]  # e.g. a list re-ordered by a GIS tool
    built = wv_footprint_from_imagettes(shuffled)
    assert len(built) == len(rings)
    res = estimate_wv_angles(built, 7077.0, lons, lats)
    np.testing.assert_array_equal(res.ring_index, np.arange(len(rings)))
    np.testing.assert_array_equal(res.subswath, subswaths)
    # bare geometries (sequence order = time order), WKT and a dateline MultiPolygon
    split = {"type": "MultiPolygon", "coordinates": [
        [[[179.9, 0.0], [180.0, 0.0], [180.0, 0.2], [179.9, 0.2], [179.9, 0.0]]],
        [[[-180.0, 0.0], [-179.9, 0.0], [-179.9, 0.2], [-180.0, 0.2], [-180.0, 0.0]]]]}
    out = wv_footprint_from_imagettes([split, "POLYGON ((0 0, 1 0, 1 1, 0 0))"])
    assert len(out) == 3 and out[0][0, 0] == 179.9 and out[2].shape == (4, 2)
    with pytest.raises(ValueError, match="no geometry"):
        wv_footprint_from_imagettes([{"index": 1, "Name": "x"}])
