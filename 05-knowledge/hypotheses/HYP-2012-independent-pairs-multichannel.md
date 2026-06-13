---
id: HYP-2012
status: PARTIALLY-TRUE
source: oracle-2026-06-01-S532
related:
  - HYP-2011
  - HYP-2008
  - HYP-2003
---

# HYP-2012: the LRC multi-channel metric is the independent-pair count floor(n/2) (a perfect matching)

**Independent pair** = two disjoint arcs = orthogonal roots in A_{n-1} (90-degree
'independent' arcs). Max set = perfect matching, size floor(n/2).

**VERIFIED (`independent_pairs_channels_s532.py`):**
1. n=4 (2 independent pairs): with a suitable fixed frame (the 4 cross arcs), the
   2 matching-arc bits index ALL 4 iso classes A000568(4) — 8 of 16 frames give the
   bijection. The 4-tournament iso class = state of its 2 independent pairs (the
   user's claim). 
2. independent pairs = floor(n/2): 1,2,2,3,3,4,7 for n=3,4,5,6,7,8,14. **n=14 -> 7 =
   opus-S524's 7 CRT classes** (6 pairs {i,i+7} + singleton {7} = speed n/2 = the
   diameter pair, S530). The CRT channels ARE the independent pairs.
3. The n=4 parity law (S531) is the one-channel degeneration: a+b+c = parity(pair1)
   XOR parity(pair2), so the 2 pairs fuse to a single Z/2 bit -> half of n=4 proved.
4. n=6 is irreducibly multi-channel: of 120 primitive 5-sets, EVERY one has an active
   full-support resonance -> no single congruence kills the inside debt; the vanishing
   is a JOINT state of the 3 independent pairs.

**CLAIM (the multi-channel generalization):** LRC(n) is a condition on the joint
STATE of the floor(n/2) independent pairs; the inside debt vanishes / is bounded for
joint states satisfying a congruence mod (n or n/2). n=4: 1 bit (parity). n=14: 7
channels = opus CRT; "musical chairs" = rotation of which pair is the last blocker.

**OPEN / next:** compute the n=6 inside debt as a function of the 3-pair joint state;
find the 3-channel analogue of "a+b+c odd" (a condition mod 6). Characterize which
of the 8/16 n=4 frames give the iso bijection (likely the strongly-connected 4-cycle
frames). Relate the special pair (speed n/2) to the source-sink/diameter (HYP-2008).

**Files:** `04-computation/independent_pairs_channels_s532.py` (+.out),
`04-computation/lrc_independent_pairs_s532.py` (+.out). Reflection:
`07-reflections/independent-pairs-are-the-channels-s532.md`.
