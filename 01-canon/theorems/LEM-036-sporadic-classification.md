---
id: LEM-036
title: THE SPORADIC CLASSIFICATION FOR ARBITRARY CLUSTERS (LEM-035's named open RESOLVED). (A) ORBIT REDUCTION: at class-0 crossings x = k/(7m) of a 7-full owner in ANY cluster, k = mj + r, the left-section vector is v_f(j,r) = φ_f·j + o_f(r) (mod 7) with φ_f = f mod 7, o_f(r) = ⌊fr/m⌋ − [m|fr] — the survivor set S_r = {j : 0 ∉ v, {1..5} ⊆ v} depends ONLY on (φ, o(r)) ∈ Z₇⁶ × Z₇⁶; runners with 7 | f are CONSTANT coordinates. (B) FULL-COLUMN THEOREM: |S_r| = 6 ⟺ φ is a permutation of {1..6} AND o(r) = t·φ (the shifted multiplication permutation; S_r = Z₇∖{−t}); |S_r| = 7 impossible whenever a moving coordinate exists. Proof: the six zero-hits j_f = −o_f/φ_f must coincide (forcing o = tφ), and non-injective φ has value-set ≤ 5 whose (j+t)-scalings equal {1..5} for at most one j. LEM-035's clean/rescued columns are exactly the t = 0 instances; ALL 12 real full columns measured have t = 0 and occur only on family clusters (elsewhere φ is not a permutation). (C) CONSTANT-COORDINATE LAWS: o_f = 0 at 7|f kills the column; otherwise the constant fills one required section — why two-large (φ = (1,2,3,4,5,0)) has only partial columns yet several (its r = 1 columns get S = {1, j*} with 6j* ≡ const). (D) TWO RECURRING SPORADIC SHAPES decoded: the ALL-BOUNDARY column r = 0 (o = (6,…,6), alive at j = 6 exactly when moving support = {1..5} + a nonzero constant — near-AP k=12, two-large k=48) and the NEAR-TOP column r = m−2 with carry (0,1,2,3,4,4) (the LEM-035(D) mechanism; recurs at two-large-84 r=10 with the constant coordinate playing the f=6 role). CENSUS: 126 columns over 11 clusters = 12 full + 15 sporadic + 99 dead; count law Σ_r|S_r| = #survivors exact on every owner
status: PROVED ((A) 2 lines; (B) both directions as above; (C)/(D) short) + REFEREED EXACT (11 clusters, every 7-full owner, every column vs the census-verified S64 formula; full-column theorem checked both directions on all 126 columns; |S_r| = 7 never)
source: boxeph-2026-07-17-S65 (owner directive: sporadic classification for arbitrary clusters)
depends_on: [LEM-035 (the named open; the formula), LEM-034 (lattice)]
script: 04-computation/lrc14_sporadic_classification_boxeph_S65.py -> 05-knowledge/results/lrc14_sporadic_classification_boxeph_S65.out
---

# LEM-036 — the sporadic classification

Columns are graded by |S_r| ∈ {0,…,6}: **6 = shifted multiplication
permutations (o = tφ, φ a permutation of {1..6}) — the only full columns,
and only family-shaped clusters can have them; 7 is impossible; everything
in between is a sporadic, computable from (φ, o) by a 7-step mod-7 check.**
The 15 real sporadic columns across all referee geometries reduce to the
all-boundary shape (r = 0), the near-top carry shape (r = m−2), and
constant-coordinate fills — each mechanism identified in the tables.

## Evidence log
- [x] orbit reduction exact everywhere; count law exact on every owner
- [x] full-column theorem both directions, 126 columns; all real t = 0
- [x] sporadic tables with mechanisms (15 columns, 18 survivors)
- [x] named open RESOLVED at the real-cluster strata (LEM-040, S67): achievability gap |S| ∉ {4,5}; QR-triple law; full-column theorem by exhaustion; residual: the gap over all 924 φ-multisets (decide-shaped)
