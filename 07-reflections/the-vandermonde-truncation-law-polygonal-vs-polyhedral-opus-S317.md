# The Vandermonde truncation law — polygonal vs polyhedral, Moser vs Fibonacci, and the self-similar difference triangle

*opus-2026-07-15-S317; owner: generalize the triangulars along both axes
(polygonal, polyhedral), the two triangles' row/skip sums (2ⁿ & Fibonacci vs
A000127 & 1,1,2,3,5,8,13,21,33,51,76,111,157,218), find the diagonal
difference patterns (1,3,6,10,15 then 4,13,28,50,80, and deeper), understand
the fundamental objects.*

## The two triangles, precisely

- **Polyhedral (simplex) triangle** = Pascal: column j holds C(i+j−1, j),
  i = 1, 2, … — ones, naturals, triangulars, tetrahedrals, pentatopes, …
- **Polygonal triangle**: column j holds the (j+1)-gonal numbers
  P_{j+1}(i) — ones, naturals, triangulars, squares, pentagonals, … Both
  columns start 1, j+1 ("one, then a natural index"), and the first three
  columns coincide — exactly as the owner observed.

Verified exactly: polyhedral row sums = 2ⁿ, shallow-diagonal (skip) sums =
Fibonacci; polygonal row sums = **A000127** (1,2,4,8,16,31,57,99,163,256 —
Moser's circle sequence), polygonal skip sums = **1,1,2,3,5,8,13,21,33,51,
76,111,157,218** (first Fibonacci deviation at the 9th term, mirroring
A000127's first 2ⁿ deviation at the 6th).

## The one identity underneath (the fundamental object)

Vandermonde:  **C(i+j−1, j) = Σ_{k≥0} C(i, k+1)·C(j−1, k)**,
and the polygonal closed form  **P_{j+1}(i) = i + (j−1)·C(i, 2)**
is EXACTLY the k ≤ 1 truncation of that sum. Therefore:

> **Polygonal = the first two Vandermonde layers of polyhedral.** The whole
> difference triangle is the tail:
> **Δ(i, j) = Σ_{k≥2} C(i, k+1)·C(j−1, k)**  (verified i = 3..7, all j).

Consequences, each matching an owner observation:
- i = 3 diagonal: Δ = C(j−1, 2) — the triangular numbers 1, 3, 6, 10, 15. ✓
- i = 4 diagonal: Δ = C(j−1, 3) + 4·C(j−1, 2) = (j−1)(j−2)(j+9)/6 —
  the sequence 4, 13, 28, 50, 80. ✓
- i = 5: Δ = C(j−1,4) + 5C(j−1,3) + 10C(j−1,2) → 10, 35, 81, 155, 265. ✓
- **General: the coefficient row of the i-th diagonal is Pascal's own row**
  (C(i,3), C(i,4), …, C(i,i)) — the difference between the two triangles is
  a self-similar copy of Pascal, shifted two layers deep. The deeper you
  compare, the more of Pascal you find inside the discrepancy.

## The truncation cascade (Moser : 2ⁿ :: this : Fibonacci)

The same k ≤ 1-vs-tail split explains all four sums at once:
- Row direction: Σ_j C(n−j+1+j−1, j) = 2ⁿ; truncating each entry to its
  polygonal (k ≤ 1) part gives Σ ≡ A000127 — which is the classical
  Σ_{k≤4} C(n, k): Moser's sequence IS truncated-2ⁿ, and the polygonal
  triangle implements the truncation entry-wise.
- Diagonal direction: Σ_j C(n−j, j) = F(n+1); the entry-wise k ≤ 1
  truncation gives the owner's sequence:
  **skip(n) = Σ_j [(n−2j+1) + (j−1)·C(n−2j+1, 2)]**, i.e.
  skip(n) = F(n+1) − Σ_j Σ_{k≥2} C(n−2j+1, k+1)C(j−1, k), with deviations
  0,…,0,1,4,13,33,76,159 beginning at n = 8 — **the exact Fibonacci analogue
  of Moser's relation to the powers of 2**, as predicted. (The deviations
  are themselves tail-sums of the same shape — the cascade recurses.)

## Why this belongs in THIS repo

The staircase is δ_{n−2} with m = T_{n−2} tiles; n² = T_{n−1} + T_n is the
two-staircase gluing (Mode A/B); the transitive ceiling (n³−n)/3 is a
polyhedral-column value; and S316's level law is the squares-vs-triangulars
parity split of the axis. The polygonal/polyhedral divergence — flat
(quadratic, k ≤ 1) growth vs simplex (full-binomial) growth — is the same
divergence the metagraph shows between per-level structure (polynomial
data: scores, x) and full-class structure (exponential data: fibers, 2^m):
the Vandermonde layer index k is the "depth of interaction," and every
truncation-at-k object in this message (polygonal numbers, A000127, the
skip sequence) is a depth-≤-k shadow of a depth-∞ object (simplex numbers,
2ⁿ, Fibonacci). Same grammar as S315's finite-depth resonance (Farey-14 vs
golden) and S316's 2-adic OCF tower (each digit needs one more layer):
**the repo's recurring meta-pattern is exact truncation laws — finite-depth
shadows with Pascal-structured discrepancies.**

Cross-refs: THM-865 (the level walk this session), HYP-6935 (S316 parity
probes), the-farey-14-row (S315), everything-is-the-triangle (the
foundation), eureka-zeckendorf-simplex-cuboid.
