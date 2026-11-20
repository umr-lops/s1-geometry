"""Add the WV angle columns to a table of reference points (pandas).

Typical use: a Parquet of WV imagettes with the product footprint (the CDSE
OData ``Footprint``, one MultiPolygon per product) and a reference point per
row (``ref_lon`` / ``ref_lat``, e.g. the imagette centre)::

    import pandas as pd
    from s1_geometry.table import add_wv_angle_columns

    df = pd.read_parquet("wv_imagettes.parquet")
    df = add_wv_angle_columns(df, footprint="Footprint", lon="ref_lon", lat="ref_lat")
    df.to_parquet("wv_imagettes_with_angles.parquet")

Only metadata is needed: the satellite distance comes from the orbit model of
:mod:`s1_geometry.orbit`.  Requires the optional ``parquet`` extra (pandas).
"""

from __future__ import annotations

import json

import numpy as np

from s1_geometry.wv import estimate_wv_angles

__all__ = ["WV_COLUMNS", "add_wv_angle_columns"]

#: Names of the columns written by :func:`add_wv_angle_columns`.
WV_COLUMNS = (
    "sar_incidence_angle",
    "sar_elevation_angle",
    "sar_ground_heading",
    "sar_wv_subswath",
)

def _key(value) -> str:
    return value if isinstance(value, str) else json.dumps(value, sort_keys=True, default=list)


def add_wv_angle_columns(
    df,
    footprint: str = "Footprint",
    lon: str = "ref_lon",
    lat: str = "ref_lat",
    *,
    r_s_km: float | None = None,
    window: int | str | None = "auto",
    on_outside: str = "nan",
    on_error: str = "raise",
):
    """Return a copy of ``df`` with the WV incidence / elevation / heading columns.

    Parameters
    ----------
    df:
        A pandas DataFrame, one row per reference point.
    footprint:
        Column holding the WV **product** footprint (WKT, GeoJSON string /
        dict or a list of rings).  Rows sharing a footprint are computed
        together.  A single-imagette polygon is not enough (see
        :func:`~s1_geometry.wv.estimate_wv_angles`).
    lon, lat:
        Columns holding the reference point longitude / latitude (deg).
    r_s_km:
        Satellite geocentric distance (km); ``None`` (default) uses the
        metadata-only orbit model.
    window:
        Passed to :func:`~s1_geometry.wv.estimate_wv_angles` (default
        ``"auto"``).
    on_outside:
        ``"nan"`` (default) or ``"raise"`` for reference points outside every
        imagette of their footprint.
    on_error:
        ``"raise"`` (default) or ``"nan"``: what to do when a footprint cannot
        be processed (malformed, fewer than 4 imagettes, ...).

    Returns
    -------
    pandas.DataFrame
        A copy of ``df`` with the :data:`WV_COLUMNS` added (the elevation is
        the radar elevation angle of the ESA L2 OCN convention, the heading is
        clockwise from North).
    """
    if on_error not in ("raise", "nan"):
        raise ValueError(f"on_error must be 'raise' or 'nan', got {on_error!r}")
    out = df.copy()
    n = len(out)
    values = {name: np.full(n, np.nan) for name in WV_COLUMNS}
    keys = np.array([_key(v) for v in out[footprint].to_numpy()], dtype=object)
    lons = out[lon].to_numpy(dtype=float)
    lats = out[lat].to_numpy(dtype=float)
    import pandas as pd

    groups = pd.Series(keys).groupby(keys, sort=False).indices
    for rows in groups.values():
        geom = out[footprint].iloc[rows[0]]
        try:
            res = estimate_wv_angles(geom, r_s_km, lons[rows], lats[rows],
                                     window=window, on_outside=on_outside)
        except ValueError:
            if on_error == "raise":
                raise
            continue
        values["sar_incidence_angle"][rows] = res.incidence_deg
        values["sar_elevation_angle"][rows] = res.elevation_deg
        values["sar_ground_heading"][rows] = res.heading_deg
        sub = res.subswath.astype(float)
        sub[np.isnan(res.incidence_deg)] = np.nan
        values["sar_wv_subswath"][rows] = sub
    for name, column in values.items():
        out[name] = column
    return out
