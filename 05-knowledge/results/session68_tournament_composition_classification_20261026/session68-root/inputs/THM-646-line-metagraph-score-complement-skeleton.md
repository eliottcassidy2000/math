---
id: THM-646
title: THE LINE METAGRAPH'S GENERATING PATTERN — (i) the SCORE-COMPLEMENT LAW s + s' = c = (n−2, n−1, …, n−1, n) holds pointwise for every line (proved + exhaustive n≤7); (ii) the line metagraph is the class fibration pulled over the DETERMINISTIC labeled-score involution s ↦ c − s; (iii) the even-graph shadow is translation by K* with K* = K_n (odd n) or K_n minus the consecutive perfect matching (even n); (iv) the kernel is CONSTRAINT-DOMINATED, not quasi-random (klein-S161 C5 refuted at n=7)
status: PROVED ((i) with proof below, verified on all 2^m tilings n=4..7; (ii) immediate from (i); (iii) proved by the parity count below, verified n=4..7). REFUTATION with data ((iv): corr(N, fiber-product) = 0.50/0.47/0.29 at n=5/6/7, max deviation 1309×). DATA: transitive-partner H = 5, 9, 17, 31 — the 2^{n−2}+1 pattern BREAKS at n=7.
source: mac-mini-2026-07-07-S47 (HYP-5017; owner: the line metagraph + the meta-abstraction move)
depends_on:
  - THM-643   # the parity layer (mac-mini-S46)
  - THM-644   # the fiber law g|Aut| = H_anti (opus-S139)
related:
  - HYP-4851  # klein-S161 (the concurrent line-side theorems; C5 refined here)
  - MISTAKE-033/065
external: Landau's theorem (score sequences); the cycle space of K_n.
---

# THM-646 — The line metagraph's generating pattern

The line metagraph L_n: nodes = iso classes, one edge per line (τ-orbit {t, flip(t)}),
colored blue/black, with multiplicity kernel N(C,D) = #lines between fibers.

## (i) The score-complement law (the deterministic skeleton)

> For EVERY tiling t, with labeled score vectors s of t and s' of flip(t):
> **s(v) + s'(v) = c(v) for every vertex v, where c = (n−2, n−1, n−1, …, n−1, n)**
> (base-path labeling; c(1) = n−2, c(n) = n, all middle entries n−1).

*Proof.* Flipping all m tiles reverses exactly the non-base-path arcs at each vertex,
so s(v) + s'(v) = (base-path wins of v, twice) + (tile degree of v)
= 2·[v ≥ 2] + (v−2)₊ + (n−1−v)₊, which evaluates to n−2 at v=1, n at v=n, n−1
otherwise. ∎ (Verified on all 2^m tilings, n = 4..7: zero exceptions.)

Sanity: Σc = n(n−1) = 2·C(n,2) ✓. Corollaries: N(C,D) = 0 unless some labeled
representatives satisfy the c-complement relation (a hard support constraint — at n=7
only 12342 of 104196 class-pairs are realized); the transitive class's partner stratum
is score multiset {1, 2, …, n−3, n−2, n−2} (verified n=4..7 exactly).

## (ii) The meta-abstraction: where the chaos comes from

On LABELED score vectors the line map is the EXACT involution s ↦ c − s — fully
deterministic, zero entropy. All apparent randomness of L_n is the geometry of the
CLASS-OVER-SCORES fibration (which iso classes realize which labeled score vectors
against the base path — Landau's polytope refined by the path). Projecting the
involution to score MULTISETS already loses determinism (unique partner multiset for
only 17/59 strata at n=7, degree tail up to 19) because classes embed their scores
against the path in many ways. **The generating pattern: L_n = (Landau score geometry
over the base path) ⋉ (the involution s ↦ c−s), with the THM-643/644 parity laws as
the mod-2 skeleton.** The "supposed chaos" is exactly the fibration's geometry; the
matching itself carries no additional information.

## (iii) The even-graph shadow: translation by K*

Every line projects, in cycle-space addresses, to the translation
addr(flip t) = addr(t) ⊕ K*, where K* = Σ_tiles (fundamental cycle). The coefficient
of a non-path edge is 1; of the path edge (k, k+1) is k(n−k) − 1 mod 2. Since
k + (n−k) = n: for ODD n the product k(n−k) is always even, so every path edge is
present: **K* = K_n** (which is an even graph exactly when n is odd). For EVEN n the
path edges with k odd drop out: **K* = K_n minus the consecutive perfect matching
{(1,2), (3,4), …, (n−1,n)}**. (Verified n=4..7: edge counts 4/10/12/21 of 6/10/15/21.)
The flip layer of the metagraph is thus, on the even-graph side, a single group
translation — maximal at odd n, matching-deficient at even n — the even/odd duality in
its sharpest form yet, and the natural mechanism site for THM-643's C2 (blue
self-loops only at even n).

## (iv) The kernel is constraint-dominated, NOT quasi-random

Against the fiber-product model N ≈ fiber(C)·fiber(D)/2^m on the realized support:
corr = 0.50 / 0.47 / 0.29 and support-mass ratio 1.42 / 2.25 / 4.75 at n = 5/6/7, with
per-pair deviations up to 1309×. klein-S161's C5 ("near-fiber-proportional with
positive assortativity") holds in its assortativity half but FAILS as proportionality:
the score-complement support constraint plus within-stratum concentration dominate.
The correct model is (i)+(ii): hard support from the score law, weights from the
score-fibration geometry — not uniform mixing.

## Data point (small-case-pattern caution)

The transitive tiling's line partner has H = 5, 9, 17, 31 at n = 4,5,6,7: equal to
2^{n−2}+1 for n = 4,5,6, then 31 ≠ 33 at n = 7 (while the partner's score stratum
follows (i) exactly at every n). Another three-case mirage caught by the fourth case
(MISTAKE-115/119 family); the CLAUDE.md principal-line neighbor (a flip-metagraph
object) and the line partner (a τ-object) genuinely diverge at n = 7.

## Files
`04-computation/gn_line_metagraph_macmini_S47.py` (+ `.out`): the law verification,
kernel census, strata quotient, K* verification, transitive-partner identification.
