# Adversarial verification verdict — bounded-spread-exact-floor (LRC14 lemma B(k))

## ENGINE: validated rigorous (gap=2/7 breakpoints added). Matches all anchors
mu({0,1,2,3})=19/21, mu(consecutive-13)=829/4620, L1 scale-invariance, sandwich within 1e-6.

## CLAIMED CONSTANTS — ALL REPRODUCED EXACTLY
- cap-14 mu_min(k), k=3..13: 1,19/21,9/14,4/7,13/35,71/220,164/735,468/2695,409/2548,5367/35035,5367/35035  [OK]
- minimizers non-consecutive (perforated) for k>=7  [OK]
- F(k) iid ceiling k=4..13 incl F(13)=3132376013/13841287201  [OK]
- refuting witnesses: dilation(round 3/2) = 6547/49980; spread-18 = 7037/59976  [OK, both < cap-14]
- cap-stability flat for k=4,5,6,7; k=7 large-spread descent cannot beat 13/35  [OK]
- rho*=0 examples EXACT: P={1,2,3},E=(0..8,10): mu=121/490, meas(GP)=29/42, rho*=0 (disjoint) [OK]
  P={1,2,3,4},E=(0,2..8,10): mu=27/98, meas(GP)=13/21, rho*=0 [OK]

## DECISIVE NEW FINDING (beyond the claim): mu(k=13) plunges FAR below claimed values
Heuristic multi-restart descent + exact confirmation (2e6-sample sandwich, diff ~1e-6):
  mu = 3303/52780 ~ 0.0626  at spread 40
  mu = 23059/412335 ~ 0.0559 at spread 45
  mu = 1314101/28198716 ~ 0.0466 at spread 62
  mu = 3736/85785 ~ 0.04355 at spread 60   <-- LOWEST FOUND
ALL of these are < 1/14 ~ 0.07143.  Claim's headline "mu_min(13) <= 7037/59976 ~ 0.117,
in a moderate-spread dip strictly below F(13) then RISES BACK" is FALSE: mu keeps dropping,
~0.044, no sign of a positive floor near F(13). Equidistributing shapes (Sidon, Beatty,
squares) give HIGH mu ~0.18-0.24 (near ceiling); the minimizers are CLUSTERED shapes.

## CONSEQUENCE
The lemma "B(k): inf_E mu(E) >= c0 > 0" is NOT established by this work and the claimed
floor 5367/35035 (and the spread-18 0.117) are NOT lower bounds — mu goes to <= 0.0436.
inf_E mu(E) is empirically << 1/14; whether it is exactly 0 is open (descent is heuristic,
not exhaustive), but it is certainly far below every "floor" the angle proposed.
The G_P-intersected target rho* hitting 0 (verified exact) is also genuine.
NET: the angle's NEGATIVE conclusions are SOUND and well-supported; its POSITIVE
constants are correct only as bounded-spread minima (not infima); it does NOT close B(k).
