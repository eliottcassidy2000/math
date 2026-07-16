---
id: THM-907  # renumbered from 906 (death-star Bernoulli-B4 first-pushed)
title: THE ANOVA ASSEMBLY FOR RESIDUE SIX — exact channel decomposition of codex's 84-weight certificate: q(a,b,c) = β₀ + Σ D(pairs) + T with β₀ = 12653/34300 ≈ 0.3689 (margin 0.1011 to the 47/100 target), singleton channels integrating to zero for EVERY speed, all large pair channels NEGATIVE (D(1:3) = −0.1057, D(1:2) = −0.0620; max |D| decays: 0.0132 beyond p+q > 20), and the zero-marginal triple channel T carrying the entire risk — CALIBRATED: T is large exactly on small-relation triples (+0.1007 at (1,2,7), relation (1,3,−1); +0.0781 at (1,4,7), relation (3,1,−1)); the universal bound (3) reduces to ONE lemma: |T(a,b,c)| ≤ C_rel/|k₁k₂k₃|_min (BV-Fourier relation tail), whose inputs are now all pinned
status: channel decomposition EXACT (rationals; q(1,4,7) = 81/175 reproduced); pair table exact to p,q ≤ 40 with visible 1/pq decay; T calibrated on 12 triples; the relation-tail lemma is REDUCED-AND-NAMED with constants, not yet proved — the honest residual (one lemma, one page of BV-Fourier)
source: mac-mini-2026-07-16-S118 (owner: build the reflection-even cubic certificate and close negative residue 6)
depends_on: [codex THM-904 (the 84-weight pointwise reduction + target (3)), THM-905, THM-903 (reflection frame: the 84 states carry the s ↦ 6−s involution with 4 fixed multisets — 44 even weights suffice)]
script: 04-computation/residue6_closure_assembly_macmini_S118.py -> 05-knowledge/results/residue6_closure_assembly_macmini_S118.out
---

# THM-907 — the ANOVA assembly

Hoeffding-decompose β on (ℤ/7)³: β = β₀ + Σ singletons + Σ pairs + triple. Then for
distinct speeds, q(a,b,c) = β₀ + Σ_{pairs} D(primitive ratio) + T(a,b,c), because sector
occupation is exactly uniform for every speed (singletons die), pair terms see only the
joint law at the primitive ratio (exact breakpoint sweep), and T is the zero-marginal
remainder.

**Exact inputs.** β₀ = 12653/34300 ≈ 0.36889; target margin 47/100 − β₀ ≈ 0.10111.
Pair table (p, q ≤ 40, exact): the large channels are NEGATIVE (D(1:3) = −0.10570,
D(1:2) = −0.06195, D(3:4) = −0.04437); largest POSITIVE observed +0.01585 (D(1:4));
max |D| beyond p + q > 20 is 0.01321 — the 1/pq decay is live. Triple samples:
T(1,2,7) = +0.10071 (relation (1,3,−1), |k₁k₂k₃| = 3), T(2,4,7) = +0.09791,
T(1,4,7) = +0.07811 (relation (3,1,−1)), generic-looking triples |T| < 0.04 — **T is a
relation-lattice functional, exactly as THM-903 predicted; the mass sits on the small
|k₁k₂k₃| strata.**

**The one remaining lemma (named, constants pinned).** |T(a,b,c)| ≤ C_rel /
min{|k₁k₂k₃| : k·(a,b,c) = 0, k ≠ 0 componentwise} with C_rel from the BV-Fourier
coefficients of the triple channel (pointwise |T-channel| ≤ 2.017; coordinate slices are
7-step functions of total variation ≤ 2·Σ|jumps|). Given the lemma with C_rel ≤ 0.31-ish
(the ζ(3)/π³ heuristic), every triple with min-relation product > C_rel/(0.101 − ΣD⁺)
satisfies (3); the complement is an EXPLICIT finite relation box, and codex's exact scan
(≤ 60, 28,876 triples, max 81/175) already covers its speed-realizations. Then (3) closes,
hence −F₆(E) ≤ 47/490 < 0.097 — negative residue 6 done. The reflection-even reduction
(THM-903): the involution s ↦ 6−s on the 84 states has 4 fixed multisets, so the
certificate search space is 44 even weights — available if the tail constant needs
tightening.
