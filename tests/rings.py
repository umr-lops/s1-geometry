"""Real S1 IW/EW SAFE footprint rings, 4 vertices (first vertex repeated at
the end).

Standard CDSE ``GeoFootprint`` order (counter-clockwise):

* ascending : S1D_IW_OCN__2SDV_20261109T225729_005397_00A4B6_6EA4.SAFE
* descending: S1D_IW_OCN__2SDV_20261101T062649_005270_00A057_FEC8.SAFE
* turning   : synthetic, flat along-track (dlat -0.05) like the 2SSH rings
  seen at the southern polar turn (77.1-77.5 S).
* EW        : S1C_EW_GRDM_1SSH_20261113T043704_20261113T043809_010313_
  014992_F824_COG.SAFE (ascending, ~70 S, 65 s burst).

Product-page order (clockwise, as copied from the CDSE product page):

* USER_DESC_RING_CW: S1D_IW_GRDH_1SDV_20261113T062354_20261113T062419_005445_
  00A65A_AC02.SAFE (descending, ~47 N, south Brittany).
* ASC_RING_CW      : the ascending ring above, re-expressed in the
  product-page vertex order (same product, expected direction unchanged).
"""

ASC_RING = [
    (-76.425255, 37.274273), (-76.00177, 35.512432),
    (-73.232109, 35.912037), (-73.591484, 37.671867),
    (-76.425255, 37.274273),
]
DESC_RING = [
    (-4.947158, 35.468987), (-4.588374, 36.973103),
    (-7.510267, 37.383934), (-7.81269, 35.881718),
    (-4.947158, 35.468987),
]
TURN_RING = [
    (-174.2, -77.40), (-173.9, -77.45), (-171.6, -77.48), (-171.9, -77.43),
    (-174.2, -77.40),
]
EW_RING = [
    (-130.805969, -69.542892), (-124.031662, -72.842934),
    (-115.130974, -70.386246), (-122.45948, -67.414169),
    (-130.805969, -69.542892),
]
USER_DESC_RING_CW = [
    (-5.260309, 47.881615), (-1.791353, 47.47398),
    (-2.225839, 45.976742), (-5.596841, 46.383617),
    (-5.260309, 47.881615),
]
# ASC_RING re-expressed in the product-page (clockwise) vertex order.
ASC_RING_CW = [
    (-73.232109, 35.912037), (-76.00177, 35.512432),
    (-76.425255, 37.274273), (-73.591484, 37.671867),
    (-73.232109, 35.912037),
]
