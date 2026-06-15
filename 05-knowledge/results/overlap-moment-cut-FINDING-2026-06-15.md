# Can an OVERLAP-density moment matrix cut a baby-Hodge hole? — FINDING

**Session:** mac-mini-2026-06-15 | **Verdict:** NO. The hole is an integrality gap, not a moment-positivity gap.

## The c3=8 fiber at n=6 (the (2,2,2,3,3,3) near-regular score class)

Fresh enumeration via `gentourng` (n=6, 56 tournaments). The c3=8 fiber is EXACTLY
the 5 tournaments of score (2,2,2,3,3,3) — forced by Kendall–Moran
c3 = C(n,3) − Σ C(s_i,2) (verified on all 56; only this score gives c3=8).

| idx | c5 | p33 (overlap) | alpha_2 (disjoint) | H |
|-----|----|---------------|--------------------|---|
| 27  | 6  | 24 | 4 | 45 |
| 41  | 8  | 26 | 2 | 41 |
| 29  | 11 | 27 | 1 | 43 |
| 42  | 11 | 27 | 1 | 43 |
| 36  | 12 | 27 | 1 | 45 |

(p33 = intersecting cyclic-triangle pairs; alpha_2 = vertex-disjoint cyclic-triangle pairs; H = OCF/Rédei–Berge weighted-independent-set count over odd cycles.)

**Structural collapse:** p33 + alpha_2 = C(8,2) = 28 identically in the fiber. So the
two "overlap carriers" are ONE degree of freedom, not two: alpha_2 = 28 − p33.
The overlap moment matrix can see at most the single coordinate (p33), equivalently c5.

## Realized c5 vs holes

- c3=8 fiber: c5 realized = {6, 8, 11, 12}; **holes = {7, 9, 10}**.
- All n=6 (any c3): c5 realized = {0..9, 11, 12}; **c5=10 is a GLOBAL hole** at n=6.
  (c5=7 and c5=9 are fiber-holes but ARE realized at other c3.)

## Does an overlap PSD/Cauchy–Schwarz moment matrix cut c5=10? NO.

1. **(8,10) is the midpoint of (8,8) and (8,12), both realized.** Its
   convex combination (50/50 random-tournament blend of T8 and T12) is a genuine
   tournament-limit object with triangle density = c3=8's and 5-cycle density = c5=10's.
2. The continuous moment-feasible region for any finite family of densities is EXACTLY
   the closed convex hull of the finite density vectors (tournament-limit theory). c5=10
   lies strictly inside that hull (p33 ∈ [26.0, 26.667] at c5=10).
3. An explicit 4×4 OVERLAP triangle-type flag Gram matrix (root types: T_src, T_mid,
   T_sink, cyclic) is PSD at EVERY point of the c5 axis including the holes
   (min eigenvalue ≥ 0 at c5 = 6,7,8,9,10,11,12). Because it is literally an
   accumulated outer product Σ v⊗v of real density data, it is PSD by construction
   at every genuine limit point — so it can never exclude an interior point.

**Therefore no finite PSD relaxation built from overlap densities (or skew spectra,
or any continuous moment) can cut c5=10. The hole is an INTEGRALITY GAP**: c5=10 is
forbidden only because no INTEGER tournament realizes it, while its continuous
density relaxation is perfectly feasible (it is a midpoint of two real tournaments).

### Honest caveat
A naive 2×2 "arc co-occurrence" matrix I built came out non-PSD (off-diag 0.3 >
diag 0.2) — but that is a CONSTRUCTION ARTIFACT (raw co-occurrence counts are not a
true Gram matrix E[f_i f_j]); it is not a valid certificate and is discarded. The
genuinely Gram-structured 4×4 overlap matrix is PSD everywhere, as it must be.

### Why this is the key finding
The skew Hankel/spectral moments are PSD always and never bind (given). We now know
the OVERLAP carriers don't bind either — for the same fundamental reason: ANY moment
matrix on genuine density data is automatically PSD at every limit point, and the
holes are interior limit points. The baby-Hodge holes at n=6 are integrality gaps;
a moment/SOS hierarchy can only cut them by ADDING integrality (Boolean/rank-1
constraints), not by any finite PSD relaxation of continuous densities.

### Scripts
- `04-computation/overlap_moment_matrix_n6.py` — fiber table + c5 holes
- `04-computation/overlap_moment_psd_test.py` — convex-hull / moment-region test
- `04-computation/overlap_gram_explicit_psd.py` — explicit 4×4 overlap Gram eigenvalues
- Outputs in `05-knowledge/results/overlap_*.out`
