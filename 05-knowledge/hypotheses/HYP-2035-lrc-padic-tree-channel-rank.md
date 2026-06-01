---
id: HYP-2035
status: PARTIALLY-TRUE
source: oracle-2026-06-01-S542
related:
  - HYP-2033
  - HYP-2012
  - HYP-2003
---

# HYP-2035: the p-adic tree unifies LRC's channels/apex/moat/center; channel rank = omega(n/2); n=2p^k are the rank-one cases

**Setup.** For each prime p|n the LRC speeds/times live on the SL_2(Q_p) Bruhat-Tits
(p-adic) tree; the system sits on the PRODUCT of these trees (S399/S400/S410).

**VERIFIED (`lrc_padic_tree_channels_s542.py`):**
- The CHANNELS (S532) = residues mod n/2 = the vertices at depth v_p(n/2) of the
  product p-adic tree (CRT); count = n/2.
- The APEX / source-sink (S530), speed n/2, = the residue-0 branch (the singleton
  channel; toward the root/gate).
- n=14 is the clean case: 7-adic tree level-1 = 7 channels, residue-0 = [7] = apex,
  residues 1..6 = the 6 live pairs {i,i+7}; 6/7 nodes split at level 2 (the
  coupling = the tree depth). Moat (S410) at {2:1, 7:1}.
- The CHANNEL RANK = omega(n/2) (distinct primes of n/2). Rank-one cases are
  n = 2p^k: n=4(2),6(3),8(4),10(5),14(7),16(8),18(9),22(11); n=14 is the simplest
  non-trivial (one prime, depth one). Composite n/2 (n=12,20,30) = product of trees
  = inter-tower coupling.

**CLAIMS (open):**
- (A) The difficulty is set by the prime-power SHAPE of n/2 (its tree position),
  not the size of n. Rank-one n=2p^k may be provable by a single p-adic tower
  (moat-descent) argument; n=14 (rank-one, depth-one) is the first test.
- (B) The center=shift (S541) makes LRC a BOUNDARY condition on the tree (basepoint-
  free); the moat (S410) is the obstruction at each boundary layer, exported to
  child vertices by a gate (protector speed). LRC = the recursive moat descent
  always leaves a clear boundary ray.
- (C) The S526 covering inside-debt is a sum over the omega(n/2) prime-towers,
  vanishing exactly when the towers decouple.

**Files:** `04-computation/lrc_padic_tree_channels_s542.py` (+.out). Reflection:
`07-reflections/lrc-padic-tree-channels-are-branches-s542.md`. Unifies S410, S532,
S530, S541, S526, S525.
