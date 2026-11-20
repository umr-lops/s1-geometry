"""Regenerate every figure and statistics table of the documentation.

Each figure is written twice into ``docs/_static/figures/``: ``<name>.png``
(light theme) and ``<name>_dark.png`` (dark theme, shown by the Furo theme
in dark mode).  The ESA L2 OCN comparisons also write
``docs/_static/certification/stats.json`` and the reStructuredText tables
included by ``docs/certification.rst`` (``docs/_generated/*.inc``), so the
numbers in the documentation always come from the same run as the figures.

Usage::

    python scripts/make_figures.py                 # everything (needs the Ifremer data)
    python scripts/make_figures.py --only method   # schematic / model figures only
    python scripts/make_figures.py --only examples --only iw_ew_lat

Groups: ``method`` (viewing geometry, WV footprint, orbit model, heading
model, WV look-angle law, WV validation, SM profiles: only the repository
data, plus one WV footprint from the CDSE daily Parquets), ``examples``
(IW / EW / WV examples against the ESA L2 OCN products) and ``iw_ew_lat``
(IW / EW accuracy vs latitude over a sample of 2026 products).  Needs
``matplotlib``, ``pandas``, ``pyarrow`` and ``netCDF4``.
"""

from __future__ import annotations

import argparse
import glob
import json
import math
import sys
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
from matplotlib.patches import Arc, Circle, FancyArrowPatch, Polygon  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

import s1_geometry.wv as wv  # noqa: E402
from s1_geometry import (  # noqa: E402
    HEADING_FIT_CUTOFF_ABS_LAT,
    SM_BEAM_SIGMA,
    as_multipolygon,
    estimate_incidence_elevation,
    estimate_wv_angles,
    heading_from_lat,
    orbit_direction_from_ring,
    sar_ground_heading_uncertainty,
    satellite_radius_km,
)
from s1_geometry.angles import (  # noqa: E402
    R_E_KM,
    R_S_REF_SM_KM,
    incidence_offnadir_deg,
    sigma_of_r_s,
)

FIG = ROOT / "docs" / "_static" / "figures"
CERT = ROOT / "docs" / "_static" / "certification"
GEN = ROOT / "docs" / "_generated"

LEGACY = Path("/legacy/project/cersat/public/cache/project/mpc-sentinel1/data/esa")
PARQUET = Path("/scale/project/lops-siam-sentinel1-workbench/data/ifremer/"
               "metadata_cdse_sentinel1/daily_parquet_files_v2")

EXAMPLES = {
    "iw_asc": ("IW", "S1D_IW_OCN__2SDV_20260110T042252"),
    "iw_desc": ("IW", "S1D_IW_OCN__2SDV_20260110T063235"),
    "ew_asc": ("EW", "S1D_EW_OCN__2SDV_20260221T230012"),
    "ew_desc": ("EW", "S1D_EW_OCN__2SDV_20260416T081150"),
    "wv_asc": ("WV", "S1D_WV_OCN__2SSV_20260828T233712"),
    "wv_desc": ("WV", "S1D_WV_OCN__2SSV_20260826T015607"),
}


# --------------------------------------------------------------------------- #
# Themes (reference data-viz palette: first three categorical slots, blue ramp,
# blue-grey-red diverging pair; light and dark steps selected separately)
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class Theme:
    name: str
    surface: str
    ink: str
    ink2: str
    muted: str
    grid: str
    series: tuple[str, str, str]
    neutral: str
    seq: tuple[str, ...]
    div: tuple[str, ...]

    @property
    def cmap_seq(self):
        return LinearSegmentedColormap.from_list(f"seq_{self.name}", self.seq)

    @property
    def cmap_div(self):
        return LinearSegmentedColormap.from_list(f"div_{self.name}", self.div)


_BLUE = ("#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b")
LIGHT = Theme("light", "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9",
              ("#2a78d6", "#eb6834", "#1baf7a"), "#f0efec", _BLUE,
              ("#184f95", "#2a78d6", "#86b6ef", "#f0efec", "#f2a3a2", "#e34948", "#a8302f"))
DARK = Theme("dark", "#1a1a19", "#f0efec", "#c3c2b7", "#898781", "#2c2c2a",
             ("#3987e5", "#d95926", "#199e70"), "#383835", _BLUE[::-1],
             ("#86b6ef", "#3987e5", "#1c5cab", "#383835", "#a8302f", "#e66767", "#f2a3a2"))
THEMES = (LIGHT, DARK)
SIX = {"light": ("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"),
       "dark": ("#3987e5", "#d95926", "#199e70", "#c98500", "#d55181", "#0ca30c")}


def _apply(t: Theme) -> None:
    plt.rcdefaults()
    plt.rcParams.update({
        "figure.facecolor": t.surface, "axes.facecolor": t.surface,
        "savefig.facecolor": t.surface, "axes.edgecolor": t.muted,
        "axes.labelcolor": t.ink2, "axes.titlecolor": t.ink, "text.color": t.ink,
        "xtick.color": t.muted, "ytick.color": t.muted, "xtick.labelcolor": t.ink2,
        "ytick.labelcolor": t.ink2, "grid.color": t.grid, "grid.linewidth": 0.6,
        "axes.grid": True, "axes.axisbelow": True, "axes.spines.top": False,
        "axes.spines.right": False, "font.size": 8.5, "axes.titlesize": 9.5,
        "axes.titleweight": "bold", "axes.titlelocation": "left",
        "axes.labelsize": 8.5, "legend.frameon": False, "legend.fontsize": 8,
        "lines.linewidth": 2.0, "lines.markersize": 4, "figure.dpi": 160,
    })


def _save(fig, name: str, t: Theme) -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    suffix = "" if t.name == "light" else "_dark"
    fig.savefig(FIG / f"{name}{suffix}.png", bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)


def themed(name: str):
    """Decorator: draw ``fn(theme)`` once per theme and save both PNGs."""
    def deco(fn):
        def run(*args, **kwargs):
            for t in THEMES:
                _apply(t)
                fig = fn(t, *args, **kwargs)
                _save(fig, name, t)
            print(f"figure {name}")
        return run
    return deco


def _stats(err: np.ndarray) -> dict:
    err = err[np.isfinite(err)]
    if err.size == 0:
        return {"n": 0}
    return {"n": int(err.size), "rmse_deg": float(np.sqrt(np.mean(err ** 2))),
            "bias_deg": float(err.mean()), "p95_abs_deg": float(np.percentile(np.abs(err), 95)),
            "max_abs_deg": float(np.abs(err).max())}


def _wrap(d):
    return (np.asarray(d) + 180.0) % 360.0 - 180.0


# --------------------------------------------------------------------------- #
# Data access (Ifremer HPC paths)
# --------------------------------------------------------------------------- #
def _safe(prefix: str) -> Path:
    sat, mode = prefix[:3], prefix[4:6]
    day = prefix.split("_")[5][:8]
    from datetime import datetime
    doy = datetime.strptime(day, "%Y%m%d").strftime("%Y/%j")
    hits = sorted(glob.glob(str(LEGACY / f"sentinel-1{sat[2].lower()}" / "L2" / mode /
                                f"{sat}_{mode}_OCN__2S" / doy / f"{prefix}*.SAFE")))
    if not hits:
        raise FileNotFoundError(prefix)
    return Path(hits[0])


def _odata_row(name: str) -> dict:
    """The CDSE daily-Parquet row of a product (footprints, attributes)."""
    import pyarrow.parquet as pq
    day = name.split("_")[5][:8]
    path = PARQUET / f"s1_{day[:4]}-{day[4:6]}-{day[6:8]}.parquet"
    table = pq.read_table(path, columns=["Name", "Footprint", "GeoFootprint", "Attributes"],
                          filters=[("Name", "==", name)])
    if table.num_rows == 0:
        raise KeyError(name)
    row = table.to_pylist()[0]
    row["attrs"] = {a["Name"]: a["Value"] for a in json.loads(row["Attributes"])}
    return row


def _read_ocn_grid(safe: Path):
    """IW / EW owi grid: lon, lat, incidence, elevation, heading (2D)."""
    import netCDF4
    f = sorted((safe / "measurement").glob("*-ocn-*.nc"))[0]
    with netCDF4.Dataset(f) as ds:
        out = [np.ma.filled(ds[v][:].astype(float), np.nan) for v in (
            "owiLon", "owiLat", "owiIncidenceAngle", "owiElevationAngle", "owiHeading")]
    return out


def _model_iw_ew(ring, mode, lon, lat, step=1):
    inc = np.full(lon.shape, np.nan)
    el = np.full(lon.shape, np.nan)
    orbit = orbit_direction_from_ring(ring)
    hd = np.full(lon.shape, np.nan)
    for i in range(0, lon.shape[0], step):
        for j in range(0, lon.shape[1], step):
            if not (np.isfinite(lon[i, j]) and np.isfinite(lat[i, j])):
                continue
            try:
                a = estimate_incidence_elevation(ring, float(lon[i, j]), float(lat[i, j]), mode)
            except ValueError:
                continue
            inc[i, j], el[i, j] = a.incidence_deg, a.elevation_deg
            hd[i, j] = heading_from_lat(float(lat[i, j]), orbit)
    return inc, el, hd, orbit


# --------------------------------------------------------------------------- #
# Methodology figures
# --------------------------------------------------------------------------- #
@themed("method_geometry")
def fig_geometry(t: Theme):
    fig, ax = plt.subplots(figsize=(7.4, 3.9))
    ax.set_aspect("equal")
    ax.axis("off")
    r_e, r_s, sigma = 1.0, 2.15, 19.0
    s = np.array([0.0, r_s])
    beam = np.array([math.sin(math.radians(sigma)), -math.cos(math.radians(sigma))])
    b = s @ beam
    dist = -b - math.sqrt(b * b - (r_s ** 2 - r_e ** 2))
    p = s + dist * beam
    gamma = math.degrees(math.atan2(p[0], p[1]))
    ax.add_patch(Circle((0, 0), r_e, facecolor=t.neutral, edgecolor=t.muted, lw=1.2))
    ax.plot([0, 0], [0, r_s], color=t.muted, lw=1.0, ls=(0, (4, 3)))
    ax.plot([0, p[0] * 1.55], [0, p[1] * 1.55], color=t.muted, lw=1.0, ls=(0, (4, 3)))
    ang = np.radians(np.linspace(90, 90 - gamma, 40))
    ax.plot(r_e * np.cos(ang), r_e * np.sin(ang), color=t.series[1], lw=3.2,
            solid_capstyle="round")
    ax.plot([s[0], p[0]], [s[1], p[1]], color=t.series[0], lw=2.2)
    ax.add_patch(Circle(s, 0.05, color=t.ink, zorder=5))
    ax.add_patch(Circle(p, 0.04, color=t.series[1], zorder=5))
    ax.add_patch(Circle((0, r_e), 0.03, color=t.ink2, zorder=5))
    ax.add_patch(Arc(s, 0.9, 0.9, theta1=270, theta2=270 + sigma, color=t.ink, lw=1.2))
    ax.text(s[0] + 0.07, s[1] - 0.62, "σ", fontsize=13, color=t.ink)
    ax.add_patch(Arc((0, 0), 0.5, 0.5, theta1=90 - gamma, theta2=90, color=t.ink, lw=1.2))
    ax.text(0.05, 0.30, "γ", fontsize=13, color=t.ink)
    up = math.degrees(math.atan2(p[1], p[0]))
    to_s = math.degrees(math.atan2(s[1] - p[1], s[0] - p[0]))
    ax.add_patch(Arc(p, 0.5, 0.5, theta1=up, theta2=to_s, color=t.ink, lw=1.2))
    ax.text(p[0] + 0.13, p[1] + 0.26, "θ", fontsize=13, color=t.ink)
    ax.text(-0.12, r_s + 0.02, "satellite", fontsize=8.5, color=t.ink2, ha="right")
    ax.text(-0.08, r_e + 0.06, "nadir", fontsize=8.5, color=t.ink2, ha="right")
    ax.text(p[0] + 0.09, p[1] - 0.14, "imaged point", fontsize=8.5, color=t.ink2)
    ax.text(p[0] * 1.55 + 0.06, p[1] * 1.55, "local vertical", fontsize=8, color=t.muted,
            ha="left", va="center")
    ax.text(0.0, -0.18, "Earth centre", fontsize=8, color=t.muted, ha="center")
    ax.add_patch(Circle((0, 0), 0.02, color=t.muted))
    ax.text(0.30, 1.88, "line of sight", fontsize=8, color=t.series[0])
    ax.text(1.65, 1.95,
            "Inputs from the metadata\n"
            "  footprint  →  where the swath is\n"
            "  r_s  →  satellite distance from Earth centre\n\n"
            "Angles at one point\n"
            "  σ  off-nadir angle at the satellite\n"
            "      (= ESA \"radar elevation angle\")\n"
            "  θ  incidence angle, from the local vertical\n"
            "  γ  Earth-centre angle, nadir → point\n"
            "  orange arc: ground range (nadir → point)\n\n"
            "  θ = σ + γ        sin θ = (r_s / R) · sin σ",
            fontsize=8.5, color=t.ink, va="top")
    ax.set_xlim(-1.15, 4.6)
    ax.set_ylim(-1.05, 2.35)
    ax.set_title("Viewing geometry in the plane across the track (not to scale)")
    return fig


def _wv_example_geometry(name="S1D_WV_OCN__2SSV_20260828T233712"):
    row = _odata_row(next(n for n in _odata_names(name)))
    rings = [wv._closed(r) for r in wv.merge_antimeridian_splits(as_multipolygon(row["Footprint"]))]
    centers = np.array([wv._to_3d(r).mean(axis=0) for r in rings])
    ring_lat = np.array([r[:, 1].mean() for r in rings])
    dirs = wv._pass_direction_per_ring(centers)
    r_s = np.array([satellite_radius_km(la, d) for la, d in zip(ring_lat, dirs, strict=True)])
    q, heading, is_wv1 = wv._solve_window(centers, r_s, ring_lat)
    q_lon = np.degrees(np.arctan2(q[:, 1], q[:, 0]))
    q_lat = np.degrees(np.arcsin(q[:, 2] / np.linalg.norm(q, axis=1)))
    return row, rings, is_wv1, q_lon, q_lat, heading


def _odata_names(prefix):
    import pyarrow.parquet as pq
    day = prefix.split("_")[5][:8]
    path = PARQUET / f"s1_{day[:4]}-{day[4:6]}-{day[6:8]}.parquet"
    names = pq.read_table(path, columns=["Name"]).column(0).to_pylist()
    return [n for n in names if n.startswith(prefix)]


@themed("method_wv_footprint")
def fig_wv_footprint(t: Theme, data):
    row, rings, is_wv1, q_lon, q_lat, _ = data
    lat0 = float(np.mean([r[:, 1].mean() for r in rings]))
    fig, ax = plt.subplots(figsize=(4.8, 5.6))
    for k, (ring, w1) in enumerate(zip(rings, is_wv1, strict=True)):
        c = t.series[0] if w1 else t.series[1]
        ax.add_patch(Polygon(ring, closed=True, facecolor=c, edgecolor=c, alpha=0.85, lw=0.8))
        if k in (0, 1, len(rings) - 2, len(rings) - 1):
            cx, cy = ring[:-1].mean(axis=0)
            ax.text(cx + (0.35 if w1 else 0.35), cy, f"#{k + 1}", fontsize=7, color=t.ink2,
                    va="center")
    ax.plot(q_lon, q_lat, color=t.ink2, lw=1.2, ls=(0, (4, 3)))
    ax.plot([], [], color=t.series[0], lw=6, label="WV1 imagettes (near, ~23° incidence)")
    ax.plot([], [], color=t.series[1], lw=6, label="WV2 imagettes (far, ~36° incidence)")
    ax.plot([], [], color=t.ink2, lw=1.2, ls=(0, (4, 3)),
            label="estimated sub-satellite track")
    i0, i1 = len(q_lon) // 3, len(q_lon) // 3 + 3
    ax.add_patch(FancyArrowPatch((q_lon[i0], q_lat[i0]), (q_lon[i1], q_lat[i1]),
                                 arrowstyle="-|>", mutation_scale=14, color=t.ink, lw=1.5))
    ax.text(q_lon[i1] - 0.9, q_lat[i1], "flight", fontsize=8, color=t.ink, ha="right")
    mid = len(rings) // 2
    tgt = rings[mid + 1][:-1].mean(axis=0)
    ax.add_patch(FancyArrowPatch((q_lon[mid], q_lat[mid]), tuple(tgt), arrowstyle="-|>",
                                 mutation_scale=12, color=t.muted, lw=1.0, ls="--"))
    ax.text((q_lon[mid] + tgt[0]) / 2, (q_lat[mid] + tgt[1]) / 2 + 0.25,
            "radar looks right", fontsize=8, color=t.muted, ha="center")
    ax.set_aspect(1.0 / math.cos(math.radians(lat0)))
    ax.autoscale_view()
    ax.margins(0.08)
    ax.set_xlabel("longitude (°)")
    ax.set_ylabel("latitude (°)")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=1)
    ax.set_title(f"WV footprint: {len(rings)} imagettes in time order\n{row['Name'][:32]}…",
                 fontsize=9)
    return fig


@themed("method_iw_footprint")
def fig_iw_footprint(t: Theme, ring_text: str, name: str):
    from s1_geometry.angles import _range_lines
    from s1_geometry.footprint import _as_point_array
    pts = _as_point_array(ring_text)[:4]
    lines = _range_lines(ring_text)
    near, far = np.array(lines["near_mid"]), np.array(lines["far_mid"])
    lat0 = float(pts[:, 1].mean())
    fig, ax = plt.subplots(figsize=(6.0, 4.6))
    ax.add_patch(Polygon(pts, closed=True, facecolor=t.neutral, edgecolor=t.muted, lw=1.0))
    edges = [(pts[i], pts[(i + 1) % 4]) for i in range(4)]
    for mid, c, label in ((near, t.series[0], "near range line (lowest incidence)"),
                          (far, t.series[1], "far range line (highest incidence)")):
        a, b = min(edges, key=lambda e, m=mid: np.linalg.norm((e[0] + e[1]) / 2 - m))
        ax.plot([a[0], b[0]], [a[1], b[1]], color=c, lw=3.2, label=label,
                solid_capstyle="round")
    ax.add_patch(FancyArrowPatch(tuple(near), tuple(far), arrowstyle="<|-|>",
                                 mutation_scale=12, color=t.ink, lw=1.2))
    mid = (near + far) / 2
    ax.text(mid[0], mid[1] + 0.12, f"swath width W = {lines['width_km']:.0f} km",
            fontsize=8, color=t.ink, ha="center")
    ref = near + 0.3 * (far - near) + np.array([0.0, -0.45])
    ax.plot(*ref, "o", ms=7, color=t.ink, mfc=t.surface, mew=1.6)
    ax.text(ref[0] + 0.06, ref[1] - 0.18, "reference point\nacross-track position 0.3",
            fontsize=8, color=t.ink2)
    track = (pts[1] - pts[0]) if np.linalg.norm(pts[1] - pts[0]) < np.linalg.norm(
        pts[2] - pts[1]) else (pts[2] - pts[1])
    if (orbit_direction_from_ring(ring_text) == "ascending") != (track[1] > 0):
        track = -track
    start = pts.min(axis=0) + np.array([0.0, 0.2])
    ax.add_patch(FancyArrowPatch(tuple(start), tuple(start + 0.45 * track), arrowstyle="-|>",
                                 mutation_scale=14, color=t.ink, lw=1.5))
    ax.text(*(start + 0.5 * track + np.array([0.05, 0.0])), "flight", fontsize=8, color=t.ink)
    ax.set_aspect(1.0 / math.cos(math.radians(lat0)))
    ax.autoscale_view()
    ax.margins(0.12)
    ax.set_xlabel("longitude (°)")
    ax.set_ylabel("latitude (°)")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=1)
    ax.set_title(f"IW footprint (4 corners) of {name[:32]}…", fontsize=9)
    return fig


@themed("method_orbit_model")
def fig_orbit_model(t: Theme):
    import pandas as pd
    d = pd.read_csv(ROOT / "tests" / "data" / "poeorb_r_s_sample.csv")
    lat = np.linspace(-78.6, 84.9, 400)
    fig, axes = plt.subplots(1, 2, figsize=(7.6, 3.1), gridspec_kw={"wspace": 0.3})
    for asc, c, label in ((1, t.series[0], "ascending"), (0, t.series[1], "descending")):
        m = d.ascending == asc
        axes[0].plot(d.lat_wv1[m], d.r_s_km[m], "o", ms=3, color=c, alpha=0.55, mew=0)
        axes[0].plot(lat, satellite_radius_km(lat, label), color=c, lw=2, label=label)
        res = d.r_s_km[m] - satellite_radius_km(d.lat_wv1[m].to_numpy(), label)
        axes[1].plot(d.lat_wv1[m], res, "o", ms=3, color=c, alpha=0.7, mew=0, label=label)
    axes[0].set_xlabel("latitude of the imaged point (°)")
    axes[0].set_ylabel("satellite distance from Earth centre $r_s$ (km)")
    axes[0].set_title("Orbit model vs precise orbits")
    axes[0].legend(loc="lower center")
    for v in (-0.29, 0.29):
        axes[1].axhline(v, color=t.muted, lw=1, ls=(0, (4, 3)))
    axes[1].axhline(0, color=t.ink2, lw=0.8)
    axes[1].set_ylim(-0.6, 0.6)
    axes[1].set_xlabel("latitude of the imaged point (°)")
    axes[1].set_ylabel("precise orbit − model (km)")
    axes[1].set_title("Residual: 0.10 km std, 0.29 km max")
    axes[1].text(0.02, 0.04, "S1A, S1C, S1D precise orbit files (POEORB)",
                 transform=axes[1].transAxes, fontsize=7.5, color=t.muted)
    return fig


@themed("method_orbit_stability")
def fig_orbit_stability(t: Theme):
    import pandas as pd
    d = pd.read_csv(CERT / "orbit_r_s_drift_2025_2026.csv", parse_dates=["day"])
    fig, ax = plt.subplots(figsize=(7.6, 2.8))
    ax.axvspan(pd.Timestamp("2026-06-09"), pd.Timestamp("2026-06-23"), color=t.grid, lw=0)
    ax.text(pd.Timestamp("2026-06-16"), 1.55, "S1C moves to\nthe S1A slot\n(no acquisitions)",
            fontsize=7.5, color=t.ink2, ha="center", va="top")
    for k, sat in enumerate(("S1A", "S1C", "S1D")):
        g = d[d.sat == sat]
        ax.plot(g.day, g.bias_wv1, "o-", color=t.series[k], ms=4, lw=1.2, label=sat)
    ax.axhline(0, color=t.ink2, lw=0.8)
    ax.set_ylim(-2.2, 2.2)
    ax.set_ylabel("daily mean of\nprecise orbit − model (km)")
    ax.set_title("The orbit model does not depend on the satellite or its orbital slot")
    ax.legend(loc="lower right", ncol=3)
    ax.annotate("S1C commissioning", (pd.Timestamp("2025-03-15"), -1.93),
                xytext=(pd.Timestamp("2025-05-20"), -1.75), fontsize=7.5, color=t.ink2,
                arrowprops={"arrowstyle": "-", "color": t.muted, "lw": 0.8})
    return fig


@themed("method_heading_model")
def fig_heading_model(t: Theme):
    lat = np.linspace(-78, 82, 400)
    fig, axes = plt.subplots(1, 3, figsize=(7.8, 2.9), gridspec_kw={"wspace": 0.45})
    for ax, d, c in ((axes[0], "ascending", t.series[0]), (axes[1], "descending", t.series[1])):
        h = np.array([heading_from_lat(x, d) for x in lat])
        u = np.array([sar_ground_heading_uncertainty(x) for x in lat])
        ax.fill_between(lat, h - u, h + u, color=c, alpha=0.18, lw=0)
        ax.plot(lat, h, color=c)
        for v in (-HEADING_FIT_CUTOFF_ABS_LAT, HEADING_FIT_CUTOFF_ABS_LAT):
            ax.axvline(v, color=t.muted, lw=0.8, ls=(0, (4, 3)))
        ax.set_xlabel("latitude (°)")
        ax.set_title(f"{d} pass")
    axes[0].set_ylabel("ground heading (° from North)")
    for ax, d in ((axes[0], "ascending"), (axes[1], "descending")):
        core = np.array([heading_from_lat(x, d) for x in np.linspace(-75, 80, 200)])
        ax.set_ylim(core.min() - 6, core.max() + 6)
    ab = np.linspace(0, 84.5, 300)
    for k, (hemi, sign) in enumerate((("north", 1), ("south", -1))):
        lim = 84.5 if sign > 0 else 77.5
        x = ab[ab <= lim]
        axes[2].plot(x, [sar_ground_heading_uncertainty(sign * v) for v in x],
                     color=t.series[k], label=hemi)
    axes[2].axvline(HEADING_FIT_CUTOFF_ABS_LAT, color=t.muted, lw=0.8, ls=(0, (4, 3)))
    axes[2].set_yscale("log")
    axes[2].set_xlabel("|latitude| (°)")
    axes[2].set_ylabel("uncertainty U (°)")
    axes[2].set_title("uncertainty")
    axes[2].legend(loc="upper left")
    fig.suptitle("IW heading model: one polynomial per pass direction, band = ±U(lat); "
                 "dashed: end of the fitted range (70°)", x=0.01, ha="left", fontsize=8,
                 color=t.ink2, y=1.02)
    return fig


@themed("method_wv_law")
def fig_wv_law(t: Theme):
    import pandas as pd
    b = pd.read_csv(CERT / "wv_look_angle_training_binned.csv")
    lat = np.linspace(-79.5, 86, 300)
    first = {1: (0.012076, 0.2450), 2: (0.009995, 0.6560)}
    fig, axes = plt.subplots(1, 2, figsize=(7.6, 3.0), sharey=True,
                             gridspec_kw={"wspace": 0.08})
    for ax, sw in zip(axes, (1, 2), strict=True):
        g = b[b.subswath == sw]
        c = t.series[sw - 1]
        ax.fill_between(g.lat_bin, g["mean"] - g["std"], g["mean"] + g["std"], color=c,
                        alpha=0.18, lw=0)
        ax.plot(g.lat_bin, g["mean"], "o", ms=3, color=c, mew=0,
                label="needed by the ESA data (1° bins, ±1 std)")
        ax.plot(lat, np.polynomial.polynomial.polyval(lat / 90, wv._WV_THETA_LAT[sw]),
                color=t.ink, lw=1.6, label="calibrated law (degree-5 polynomial)")
        m, k = first[sw]
        ax.plot(lat, m * lat + k, color=t.muted, lw=1.4, ls=(0, (4, 3)), label="first law (linear)")
        ax.set_title(f"WV{sw}")
        ax.set_xlabel("latitude of the imagette (°)")
    axes[0].set_ylabel("off-nadir correction (°)\nadded to the documented beam angle")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.1))
    return fig


@themed("validation_wv_latitude")
def fig_wv_validation(t: Theme):
    import pandas as pd
    new = pd.read_csv(CERT / "wv_validation_s1acd.csv")
    old = pd.read_csv(CERT / "wv_validation_first_law.csv")
    for d in (new, old):
        d["latc"] = (d.lat_min + d.lat_max) / 2
    fig, axes = plt.subplots(1, 2, figsize=(7.6, 3.0), gridspec_kw={"wspace": 0.28})
    axes[0].plot(old.latc, old.inc_bias, "o", ms=4, mew=0, color=t.muted, alpha=0.6,
                 label="first law")
    for k, sat in enumerate(("S1A", "S1C", "S1D")):
        g = new[new.satellite == sat]
        axes[0].plot(g.latc, g.inc_bias, "o", ms=4, mew=0, color=t.series[k], label=sat)
        axes[1].plot(g.latc, g.hd_rmse, "o", ms=4, mew=0, color=t.series[k], label=sat)
    axes[0].axhline(0, color=t.ink2, lw=0.8)
    axes[0].set_ylabel("incidence bias per product (°)")
    axes[0].set_xlabel("product mid latitude (°)")
    axes[0].set_title("Incidence: calibrated vs first law")
    axes[0].legend(loc="upper left", ncol=2)
    axes[1].set_ylabel("heading RMSE per product (°)")
    axes[1].set_xlabel("product mid latitude (°)")
    axes[1].set_ylim(0, None)
    axes[1].set_title("Heading")
    return fig


@themed("method_sm_profiles")
def fig_sm(t: Theme):
    fig, ax = plt.subplots(figsize=(6.4, 3.0))
    x = np.linspace(0, 1, 50)
    cols = SIX[t.name]
    for k, (beam, (nmin, nmax, fmin, fmax)) in enumerate(SM_BEAM_SIGMA.items()):
        sn = sigma_of_r_s(nmin, nmax, R_S_REF_SM_KM)
        sf = sigma_of_r_s(fmin, fmax, R_S_REF_SM_KM)
        inc = [incidence_offnadir_deg(R_S_REF_SM_KM, R_E_KM, sn + v * (sf - sn)) for v in x]
        ax.plot(x, inc, color=cols[k])
        ax.text(1.02, inc[-1], beam, color=t.ink2, fontsize=8, va="center")
    ax.set_xlim(0, 1.1)
    ax.set_xlabel("across-track position in the swath (0 = near edge, 1 = far edge)")
    ax.set_ylabel("incidence (°)")
    ax.set_title("SM beams S1-S6: incidence across the swath (from the ESA specification)")
    return fig


# --------------------------------------------------------------------------- #
# ESA L2 OCN examples
# --------------------------------------------------------------------------- #
def _example_iw_ew(key, mode, prefix):
    safe = _safe(prefix)
    row = _odata_row(safe.name)
    ring = row["GeoFootprint"]
    lon, lat, inc, el, hd = _read_ocn_grid(safe)
    m_inc, m_el, m_hd, orbit = _model_iw_ew(ring, mode, lon, lat)
    err = {"incidence": m_inc - inc, "elevation": m_el - el, "heading": _wrap(m_hd - hd)}
    ref = {"incidence": inc, "elevation": el, "heading": hd}
    return {"key": key, "mode": mode, "product": safe.name, "pass": orbit, "lon": lon,
            "lat": lat, "ref": ref, "err": err,
            "lat_mean": float(np.nanmean(lat))}


def _example_wv(key, prefix):
    from validate_wv_ocn import read_owi
    safe = _safe(prefix)
    row = _odata_row(safe.name)
    owi = read_owi(safe)
    res = estimate_wv_angles(row["Footprint"], None, owi["lon"], owi["lat"], on_outside="nan")
    err = {"incidence": res.incidence_deg - owi["inc"], "elevation": res.elevation_deg - owi["el"],
           "heading": _wrap(res.heading_deg - owi["hd"])}
    ref = {"incidence": owi["inc"], "elevation": owi["el"], "heading": owi["hd"]}
    return {"key": key, "mode": "WV", "product": safe.name,
            "pass": row["attrs"].get("orbitDirection", "").lower(), "lon": owi["lon"],
            "lat": owi["lat"], "ref": ref, "err": err, "ring": res.ring_index,
            "sw": res.subswath, "model": {"incidence": res.incidence_deg,
                                          "elevation": res.elevation_deg,
                                          "heading": res.heading_deg},
            "lat_mean": float(np.nanmean(owi["lat"])), "r_s_km": res.r_s_km}


def _hist(ax, t, err, stats):
    e = err[np.isfinite(err)]
    if e.size == 0:
        ax.text(0.5, 0.5, "no reference data", transform=ax.transAxes, ha="center",
                color=t.muted)
        ax.set_yticks([])
        return
    lim = max(0.05, float(np.percentile(np.abs(e), 99.5)))
    ax.hist(e, bins=50, range=(-lim, lim), color=t.series[0], alpha=0.9, lw=0)
    ax.axvline(0, color=t.ink2, lw=0.8)
    ax.grid(axis="y", visible=False)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_xlabel(f"RMSE {stats['rmse_deg']:.3f}°  bias {stats['bias_deg']:+.3f}°  "
                  f"max {stats['max_abs_deg']:.3f}°", fontsize=7.5, color=t.ink)


def fig_example(ex):
    names = ("incidence", "elevation", "heading")
    stats = {n: _stats(ex["err"][n]) for n in names}

    @themed(f"example_{ex['key']}")
    def draw(t: Theme):
        if ex["mode"] == "WV":
            fig, axes = plt.subplots(3, 2, figsize=(7.6, 6.8),
                                     gridspec_kw={"width_ratios": [2.3, 1], "hspace": 0.7,
                                                  "wspace": 0.18})
            ring = ex["ring"]
            for r, n in enumerate(names):
                ax = axes[r, 0]
                for sw, c in ((1, t.series[0]), (2, t.series[1])):
                    m = ex["sw"] == sw
                    idx = np.unique(ring[m])
                    ref = [np.nanmean(ex["ref"][n][m & (ring == i)]) for i in idx]
                    mod = [np.nanmean(ex["model"][n][m & (ring == i)]) for i in idx]
                    ax.plot(idx + 1, mod, color=c, lw=1.6, label=f"WV{sw} model")
                    ax.plot(idx + 1, ref, "o", ms=4.5, color=c, mfc=t.surface, mew=1.3,
                            label=f"WV{sw} ESA OCN")
                ax.set_ylabel(f"{n} (°)")
                ax.set_title(f"{n}: imagette means")
                if r == 0:
                    ax.legend(ncol=2, loc="center right", fontsize=7)
                _hist(axes[r, 1], t, ex["err"][n], stats[n])
                axes[r, 1].set_title("model − ESA, all cells", fontsize=8.5)
            axes[2, 0].set_xlabel("imagette number (time order)")
        else:
            lon, lat = ex["lon"], ex["lat"]
            fig, axes = plt.subplots(3, 3, figsize=(8.4, 7.8),
                                     gridspec_kw={"width_ratios": [1, 1, 0.9], "hspace": 0.62,
                                                  "wspace": 0.62})
            aspect = 1.0 / math.cos(math.radians(ex["lat_mean"]))
            for r, n in enumerate(names):
                ref, err = ex["ref"][n], ex["err"][n]
                im = axes[r, 0].pcolormesh(lon, lat, ref, cmap=t.cmap_seq, shading="auto",
                                           rasterized=True)
                cb = fig.colorbar(im, ax=axes[r, 0], fraction=0.05, pad=0.03)
                cb.outline.set_visible(False)
                lim = max(0.05, float(np.nanpercentile(np.abs(err), 99)))
                im2 = axes[r, 1].pcolormesh(lon, lat, err, cmap=t.cmap_div, vmin=-lim,
                                            vmax=lim, shading="auto", rasterized=True)
                cb2 = fig.colorbar(im2, ax=axes[r, 1], fraction=0.05, pad=0.03)
                cb2.outline.set_visible(False)
                for ax in axes[r, :2]:
                    ax.set_aspect(aspect)
                    ax.grid(False)
                    ax.tick_params(labelsize=7)
                axes[r, 0].set_title(f"ESA {n} (°)", fontsize=9)
                axes[r, 1].set_title("model − ESA (°)", fontsize=9)
                _hist(axes[r, 2], t, err, stats[n])
                axes[r, 2].set_title("error distribution", fontsize=8.5)
        fig.suptitle(f"{ex['mode']} {ex['pass']} — {ex['product'][:-5]}", x=0.01, ha="left",
                     fontsize=8.5, color=t.ink2, y=0.995)
        return fig

    draw()
    return stats


# --------------------------------------------------------------------------- #
# IW / EW accuracy vs latitude (sample of 2026 S1D products)
# --------------------------------------------------------------------------- #
def _iw_ew_product(args):
    mode, safe = args
    try:
        row = _odata_row(Path(safe).name)
        lon, lat, inc, el, hd = _read_ocn_grid(Path(safe))
        m_inc, m_el, m_hd, orbit = _model_iw_ew(row["GeoFootprint"], mode, lon, lat, step=4)
    except (ValueError, KeyError, OSError, IndexError):
        return None
    attr = row["attrs"].get("orbitDirection", "").lower()
    e_inc = (m_inc - inc)[::4, ::4]
    e_hd = _wrap(m_hd - hd)[::4, ::4]
    ok = np.isfinite(e_inc)
    if ok.sum() < 50:
        return None
    return {"mode": mode, "product": Path(safe).name, "pass_attr": attr, "pass_ring": orbit,
            "lat_mid": float(np.nanmean(lat)), "n": int(ok.sum()),
            "inc_rmse": float(np.sqrt(np.nanmean(e_inc ** 2))),
            "inc_bias": float(np.nanmean(e_inc)),
            "hd_rmse": float(np.sqrt(np.nanmean(e_hd ** 2))),
            "covered": float(np.isfinite(m_inc[::4, ::4]).sum()
                             / max(np.isfinite(inc[::4, ::4]).sum(), 1))}


def run_iw_ew_lat(workers: int) -> None:
    import pandas as pd
    jobs = []
    for mode in ("IW", "EW"):
        days = sorted(glob.glob(str(LEGACY / f"sentinel-1d/L2/{mode}/S1D_{mode}_OCN__2S/2026/*")))
        for day in days[::4]:
            safes = sorted(glob.glob(day + "/*.SAFE"))
            jobs += [(mode, s) for s in safes[:: max(1, len(safes) // 2)][:2]]
    with ProcessPoolExecutor(workers) as pool:
        rows = [r for r in pool.map(_iw_ew_product, jobs, chunksize=1) if r]
    d = pd.DataFrame(rows)
    d.to_csv(CERT / "iw_ew_validation_2026.csv", index=False)
    print(f"iw_ew_lat: {len(d)} products of {len(jobs)}")
    fig_iw_ew_lat()


@themed("validation_iw_ew_latitude")
def fig_iw_ew_lat(t: Theme):
    import pandas as pd
    d = pd.read_csv(CERT / "iw_ew_validation_2026.csv")
    d = d[(d.pass_attr == d.pass_ring) & (d.covered > 0.9)]
    fig, axes = plt.subplots(1, 2, figsize=(7.6, 3.0), gridspec_kw={"wspace": 0.28})
    for k, mode in enumerate(("IW", "EW")):
        g = d[d["mode"] == mode]
        axes[0].plot(g.lat_mid, g.inc_bias, "o", ms=4, mew=0, color=t.series[k],
                     label=f"{mode} ({len(g)} products)")
        axes[1].plot(g.lat_mid, g.hd_rmse, "o", ms=4, mew=0, color=t.series[k], label=mode)
    axes[0].axhline(0, color=t.ink2, lw=0.8)
    axes[0].set_ylabel("incidence bias per product (°)")
    axes[1].set_ylabel("heading RMSE per product (°)")
    for ax in axes:
        ax.set_xlabel("product mid latitude (°)")
    axes[0].set_title("IW / EW incidence")
    axes[1].set_title("IW / EW heading")
    axes[1].set_yscale("log")
    axes[0].legend(loc="lower left")
    return fig


# --------------------------------------------------------------------------- #
# Statistics tables
# --------------------------------------------------------------------------- #
def write_tables(all_stats: dict) -> None:
    CERT.mkdir(parents=True, exist_ok=True)
    (CERT / "stats.json").write_text(json.dumps(all_stats, indent=2) + "\n")
    GEN.mkdir(parents=True, exist_ok=True)
    lines = [".. Generated by scripts/make_figures.py -- do not edit.", "",
             "+-------+-------------+-----------+---------+-----------+-----------+-----------+",
             "| mode  | pass        | angle     | cells   | RMSE (°)  | bias (°)  | max (°)   |",
             "+=======+=============+===========+=========+===========+===========+===========+"]
    for s in all_stats.values():
        for angle in ("incidence", "elevation", "heading"):
            a = s[angle]
            lines.append(f"| {s['mode']:<5} | {s['pass']:<11} | {angle:<9} | {a['n']:>7,} | "
                         f"{a['rmse_deg']:>9.3f} | {a['bias_deg']:>+9.3f} | "
                         f"{a['max_abs_deg']:>9.3f} |")
            lines.append("+-------+-------------+-----------+---------+-----------+"
                         "-----------+-----------+")
    (GEN / "example_stats.inc").write_text("\n".join(lines) + "\n")
    products = [".. Generated by scripts/make_figures.py -- do not edit.", ""]
    for key, s in all_stats.items():
        products.append(f".. |{key}_product| replace:: ``{s['product']}``")
    (GEN / "example_products.inc").write_text("\n".join(products) + "\n")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--only", action="append",
                        choices=("method", "examples", "iw_ew_lat"), default=None)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args(argv)
    groups = set(args.only or ("method", "examples", "iw_ew_lat"))
    if "method" in groups:
        fig_geometry()
        fig_wv_footprint(_wv_example_geometry())
        iw_name = _odata_names(EXAMPLES["iw_asc"][1])[0]
        fig_iw_footprint(_odata_row(iw_name)["GeoFootprint"], iw_name)
        fig_orbit_model()
        fig_orbit_stability()
        fig_heading_model()
        fig_wv_law()
        fig_wv_validation()
        fig_sm()
    if "examples" in groups:
        all_stats = {}
        for key, (mode, prefix) in EXAMPLES.items():
            ex = _example_wv(key, prefix) if mode == "WV" else _example_iw_ew(key, mode, prefix)
            stats = fig_example(ex)
            all_stats[key] = {"mode": mode, "pass": ex["pass"], "product": ex["product"],
                              "lat_mean": ex["lat_mean"], **stats}
        write_tables(all_stats)
    if "iw_ew_lat" in groups:
        run_iw_ew_lat(args.workers)
    return 0


if __name__ == "__main__":
    sys.exit(main())
