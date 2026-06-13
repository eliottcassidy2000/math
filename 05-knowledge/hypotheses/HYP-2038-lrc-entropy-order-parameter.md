---
id: HYP-2038
status: PARTIALLY-TRUE
source: oracle-2026-06-01-S543
related:
  - HYP-2036
  - HYP-2035
  - HYP-2033
  - HYP-2003
  - THM-391
---

# HYP-2038: LRC is a phase transition — the loneliness box-dimension is the order parameter; the tight AP is the critical point

**Setup.** On the base-p tree (depth-d node = interval [j/p^d,(j+1)/p^d); child map
x->p, entropy log p), the CENTER=shift (S541, the +1 odometer) is equicontinuous =
ZERO entropy, so LRC's entropy is transverse to the center (the x->p descent).

**VERIFIED (`lrc_tree_entropy_s543.py`, p=2,3):**
- Box-entropies of S={t: min_i ||v_i t|| >= 1/n}: h_full = lim (1/d)log_p(#full
  nodes), h_touch = lim (1/d)log_p(#nodes meeting S).
- Generic systems: h_full = h_touch = 1 (lonely set has box-dim 1, positive
  measure), base-p independent.
- Tight AP (n=4,5,6): |S|=0, h_full=-inf (no safe interval), h_touch=0 (N_touch
  BOUNDED -> a finite set of wall points). The critical/zero-entropy line.
- PHASE TRANSITION: D(theta) = dim of {min>=theta} steps 1->0 at theta_c = max
  collar = the loneliness radius (S541). AP: theta_c = 1/n EXACTLY; generic:
  theta_c > 1/n.

**CLAIM:** LRC(n) <=> theta_c(v) >= 1/n for every speed system <=> the dimension-1
(positive-entropy) phase of the lonely set reaches the threshold 1/n; equivalently
the order parameter h_touch(1/n) >= 0. The extremal AP is the system whose critical
point lands exactly on 1/n (theta_c=1/n) -- LRC is tight at the phase boundary.

**OPEN / vigorous next:**
- (A) Thermodynamic formalism: the pressure P(beta) for the potential -beta*log(collar)
  on the x->p system; its phase transition at beta_c <-> theta_c=1/n should be
  singular exactly at the AP. Compute P(beta) and its singularity.
- (B) The critical set at theta_c degenerates to the wall points a/n = the S410
  p-adic MOAT (depth v_p(n)); its multifractal spectrum is a p-adic measure (p|n)
  encoding how tightly the AP fails to clear 1/n. Compute the moat dimension.
- (C) Tie to HYP-2035: the rank omega(n/2) of the channel tree vs the critical-set
  dimension.
- (D) Tie to THM-391: entropy/box-dimension sees global criticality, but a bare
  q-grid zero branch has empty local endpoint core; the multifractal object must
  remember descendant critical-wall labels, not only branch occupancy.

**Files:** `04-computation/lrc_tree_entropy_s543.py` (+.out). Reflection:
`07-reflections/lrc-tree-entropy-the-order-parameter-s543.md`.
