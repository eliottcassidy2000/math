# HYP-7046 — The dilate rate law (exact) + the diagonal-polynomial n=5 metagraph test

**Status:** RESOLVED (death-star-2026-07-16-S23). Part 1 → **THM-898 PROVED + referee 8/8
exact** (two bases × c = 2..5; parallel-pairs = n₀ exactly; Hamiltonian-base corollary
rate ≈ 1 − 1/c explains the census's 0.801 = 1009/1260 on the nose; the doublet subcode =
a decidable dilate detector). Part 2: Lemma 2.3's bounds HOLD but are LOOSE at n=5 scale
(tight only on singletons; cannot see zero-silent-flip classes; c₃-parity sets: bound 8 vs
actual 3.6/2.0) — NOT yet a theory of level widths; the tool's regime is larger n. THE
UNEXPECTED FIND: d(W) (diagonal complexity, min delta-realizing degree) is a NEW GEOMETRIC
class invariant — it separates equal-size classes (two |W|=5, c₃=2 classes need d=2 while
the |W|=5, c₃=4 class needs d=1; the counting bound is tight for 10/12 classes, the
near-regular pair exceptional). Named follow-ups: n=6 (m=10, 56 classes) is the real
Lemma-2.3 test; the near-dilate V(c) asymptote; d(W) column for the n=6 atlas.

## Part 1 — the dilate rate law (paper derivation, to referee)

The movie of cE₀ is the base movie traversed c times (dilate covariance on states + the
(1/(7c))-grid collision structure): merged events n = c·n₀, states V = V₀, hence
> **rate(cE₀) = 1 − (V₀ − 1)/(c·n₀), k = (c−1)n₀ + k₀, d = 2 for c ≥ 2**,
and structurally **palette code(cE₀) = (c−1)·n₀ parallel-edge doublets ⊕ lift(base code)**.
(Census anchor already fits: 5·[2,3,5,7,11,13]: 1260 = 5·252, V = 252, k = 1009 = 4·252+1.)
Referee: c-ladders on two bases + doublet-count check + near-dilate (consecutive) V(c)
scaling measured against the law.

## Part 2 — the diagonal-polynomial n=5 test (2607.14068 Lemma 2.3 applied to the metagraph)

For W ⊆ F₂⁶ (tilings of an iso class / an invariant level set), the minimal degree d(W)
with full evaluation rank (delta functions realizable by deg ≤ d polynomials on W) gives,
by Lemma 2.3 on Q₆[W]: **avg internal wiggly degree of W ≤ 2·d(W)** — a candidate THEORY
for silent-mutation density and level widths (none exists). Compute at n = 5 (64 tilings,
12 classes): d(W) for all classes, H-levels, and c₃-parity sets (note: c₃ mod 2 is a
degree-≤3 polynomial in tile bits — 3-cycle indicators are cubic); compare 2d(W) vs exact
internal degrees; report tightness/informativeness honestly.

-> HYP-7036, 2607.14068 reflection, T1540, THM-584; death-star-S23.
