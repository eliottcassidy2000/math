---
id: THM-868
title: THE GEOMETRIC LADDER OF THE SUPPORT FILTRATION — the run-truncation towers' row and diagonal GFs are partial sums of geometric series in v = x²/(1−x)² (rows) and u = x³/((1−x)²(1+x)) (diagonals); the ladders telescope exactly to 1/(1−2x) and 1/(1−x−x²), and the deficits ("missing Moser regions" / "missing Fibonacci mass") have exact product GFs with leading exponents 2j+1 / 3j+2
status: PROVED (three-line algebra + telescopes; every identity machine-exact to order 40, j ≤ 6) + independently confirmed by OEIS (the stored GFs of A002624/A060099 are literally the s = 2,3 layer terms)
source: mac-mini-2026-07-15-S110 (owner: "more sequences like this, especially new ones; keep the rank-truncation perspective abstractly")
depends_on:
  - S109 figurate note (T_j tower, exact-run law C(m,s)C(k−1,s−1))
related: [THM-854(II) kps rank-2 law, opus-S317 Vandermonde truncation law, klein-S313 draft, A000127/A002664/A002620/A002624/A060099/A038505]
script: 04-computation/truncation_ladder_thm868_macmini_S110.py -> 05-knowledge/results/truncation_ladder_thm868_macmini_S110.out
---

# THM-868 — the geometric ladder

Let T_j be the run-truncation triangle (entries = #k-subsets of [r] with ≤ j maximal runs;
T_1 ⊂ T_2 = polygonal ⊂ … → Pascal), R_j(x) = Σ_r (row sum) x^r and G_j(x) = Σ_n
(shallow-diagonal sum) x^n.

**Theorem.** With v = x²/(1−x)² and u = x³/((1−x)²(1+x)):

> R_j(x) = 1/(1−x) + (1/(x(1−x)))·(v + v² + … + v^j)
> G_j(x) = 1/(1−x) + (1/(x(1−x)))·(u + u² + … + u^j)

i.e. **the support filtration is a geometric series in one substituted variable per
direction** — the s-th run-layer contributes the s-th geometric term, explicitly
x^{2s−1}/(1−x)^{2s+1} (rows) and x^{3s−1}/((1−x)^{2s+1}(1+x)^s) (diagonals). The ladders
telescope by the two one-line polynomial identities

> (1−x)² − x² = 1 − 2x   and   (1−x)²(1+x) − x³ = 1 − x − x²,

so the full sums are **exactly** 1/(1−2x) (all subsets) and 1/(1−x−x²) (Fibonacci): the
doubling law and the Fibonacci recurrence are the j → ∞ limits of one geometric ladder.
The deficit (tail) GFs are exact products:

> 2^r − R_j :  x^{2j+1} / ((1−x)^{2j+1}(1−2x))
> F − G_j  :  x^{3j+2} / ((1−x)^{2j+1}(1+x)^j(1−x−x²))

whose leading exponents re-prove the S109 first-miss law (rows r = 2j+1, diagonals
n = 3j+2, deficit 1) — Moser's missing 32nd region is the x⁵ coefficient of the j = 2 row
tail; the region deficit thereafter grows like 2^r (times a degree-2j polynomial), and the
diagonal deficit like Fibonacci (times a degree-2j quasi-polynomial).

*Proof.* Layer s of the triangle has entries C(m,s)C(k−1,s−1). Rows: Σ_r x^r Σ_k = the
convolution of the two binomial columns = x^{2s−1}/(1−x)^{2s+1} (standard). Diagonals:
m = n−2k+1 gives Σ_n x^n Σ_k C(n−2k+1, s)C(k−1, s−1) = x^{3s−1}/((1−x)^{2s+1}(1+x)^s)
(substitute y = x² in the second column's GF). Both are geometric with ratios v, u; sum
the geometric series and clear denominators with the two displayed identities. ∎

## The sequence table (OEIS sweep 2026-07-15)

| object | head | status |
|---|---|---|
| row tail j=2 (2^r − Moser) | 1, 7, 29, 93, 256, 638, 1486 | = **A002664** exactly (= Σ_{k≥5}C(r,k) = Σ_{i≥3}C(r+1,2i)) |
| diagonal tail j=2 (F − G) | 1, 4, 13, 33, 76, 159, 315, 594 | **NEW** (opus-S317's "deviations", now with exact GF) |
| layer-2 diagonals | 1, 3, 8, 16, 30, 50, 80, 120 | = **A002624**, whose stored GF is our layer formula (s=2) |
| layer-3 diagonals | 1, 4, 13, 32, 71, 140, 259 | = **A060099**, ditto (s=3) — independent confirmation |
| layer-1 diagonals | quarter-squares | = A002620 |
| d=3 Moser rows (1+C(r+2,3)+C(r+2,5)) | 1, 2, 5, 12, 27, 57, 113 | = **A362193** (Grassmannian pattern-avoiders — bijection wanted) |
| d=4 Moser rows (1+C(r+3,4)+C(r+3,6)) | 1, 2, 6, 17, 43, 99, 211 | **NEW** |
| d=3 diagonals | 1, 1, 2, 4, 8, 15, 27, 47, 79 | **NEW** (3D Fibonacci–Moser analog) |
| centered-polygonal rows (1+r+C(r+2,4)) | 1, 2, 4, 9, 20, 41, 77 | = A101338 |
| centered-polygonal diagonals | 1, 1, 2, 3, 6, 11, 20, 34, 55, 85 | **NEW** |
| odd-support rows ((2^r − Re(1+i)^{r+1})/2) | 1, 3, 6, 10, 16, 28, 56, 120 | = **A038505** (every 4th binomial — the period-8 Gaussian-integer law) |
| odd-support diagonals | 1, 2, 4, 6, 9, 12, 17, 24, 38, 62 | **NEW** (Fibonacci split) |
| G_3, G_4 | …143,228… / …144,233,377,609… | **NEW** (G_2 = G was new in S109) |
| level-multiplicity maxima (THM-866 poset) | 1, 2, 3, 4, 8, 14, 34, 72, 175 | **NEW** (tournament-native; L(n, top−1) = n−2 proved) |

Also proved en route: the **d-dimensional Moser law** — row sums of the d-dimensional
polygonal triangle are 1 + C(r+d−1, d) + C(r+d−1, d+2) (double hockey stick, exact
d ≤ 5, all r); and the **majorization refinement** — the dominance order on score
sequences is NOT graded by x/8 (covers with Δx = 16 exist from n = 6: the level lattice
of THM-866 is strictly finer than the majorization order; the walk, not the poset,
carries the completeness).

## S109 erratum (caught by this machinery)

S109's printed "Fibonacci defect" 0,1,1,2,3,5,8,13,22,… used F(n+2); the true deficit
F(n+1)-indexed is 0,…,0,1,4,13,33,76,159 (first nonzero at n = 8), matching the exact
tail GF above and opus-S317's independently computed deviations. Logged in MISTAKES.
