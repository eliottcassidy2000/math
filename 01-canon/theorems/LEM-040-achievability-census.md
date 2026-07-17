---
id: LEM-040
title: THE (φ, o) ACHIEVABILITY CENSUS (LEM-036's named open resolved at the real-cluster strata, with two discoveries). Exhaustive enumeration of all 7⁶ offset vectors on four moving strata (P6 permutation (1..6); C1 one-constant (1,2,3,4,5,0); R1 repeated slope (1,1,2,3,4,5); C2 two-constant (1,2,3,4,0,0)) + the all-constant stratum: (1) THE ACHIEVABILITY GAP: |S_r| ∈ {4, 5} is IMPOSSIBLE in every moving stratum tested — the survivor-count spectrum is {0,1,2,3,6} (P6), {0,1,2} (C1, R1 — one constant or one repeat already kills triples), {0,1,2,3} (C2); full columns (|S| = 6) occur ONLY in P6, and the FULL-COLUMN THEOREM is re-proved by exhaustion (the 7 configs with |S| = 6 are exactly o = tφ). |S| = 7 only via the all-constant stratum (2520 covering offsets — no moving coordinate, as LEM-036(B) requires). (2) THE QR-TRIPLE LAW: the P6-achievable triples are EXACTLY the 14 affine images of the quadratic-residue triple {1, 2, 4} (= the orbit with point-stabilizer of order 3; the non-residue coset {3,5,6} generates the same family) — the sporadic triples of the permutation stratum are the QR-coset triples. C2 achieves all 35 triples. (3) Across moving strata: 71 distinct achievable S = 5 affine orbits (singletons, pairs, QR-triples, all-triples via C2, six-sets); all singletons and all pairs achievable everywhere
status: PROVED BY EXHAUSTION on the five strata (4 × 117,649 + all-constant; each a finite decide) + the QR identification verified exactly; the gap for ARBITRARY moving φ-multisets is CONJECTURED (finite check: 924 φ-multisets × 7⁶ — decide-shaped, named)
source: boxeph-2026-07-17-S67 (finishing sweep; LEM-036's named open)
depends_on: [LEM-036 (orbit reduction + full-column theorem)]
script: 04-computation/lrc14_achievability_census_boxeph_S67.py -> 05-knowledge/results/lrc14_achievability_census_boxeph_S67.out
---

# LEM-040 — the achievability census

The abstract landscape of survivor sets: a GAP at sizes 4–5, full columns
only for shifted multiplication permutations, and the achievable triples of
the permutation stratum carrying quadratic-residue arithmetic — the mod-7
character theory surfacing inside the pure combinatorics of the survivor
criterion. Sporadics are: singletons, pairs, QR-triples (P6), arbitrary
triples (two constants), nothing else.

## Evidence log
- [x] four moving strata exhaustive + all-constant; gap + spectrum per stratum
- [x] full-column theorem by exhaustion; QR-triple identification exact
- [ ] named (decide-shaped): the gap over all 924 φ-multisets
