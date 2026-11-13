"""Real S1 IW/EW SAFE footprint rings from the CDSE GeoFootprint, 4 vertices
in the fixed [BL, TL, TR, BR] order (first vertex repeated at the end).

* ascending : S1D_IW_OCN__2SDV_20261109T225729_005397_00A4B6_6EA4.SAFE
* descending: S1D_IW_OCN__2SDV_20261101T062649_005270_00A057_FEC8.SAFE
* turning   : synthetic, flat along-track (dlat -0.05) like the 2SSH rings
  seen at the southern polar turn (77.1-77.5 S).
* EW        : S1C_EW_GRDM_1SSH_20261113T043704_20261113T043809_010313_
  014992_F824_COG.SAFE (ascending, ~70 S, 65 s burst).
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
