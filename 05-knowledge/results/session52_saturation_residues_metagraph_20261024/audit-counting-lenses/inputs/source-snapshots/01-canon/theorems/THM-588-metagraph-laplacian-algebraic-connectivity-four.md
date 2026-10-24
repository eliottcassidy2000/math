---
id: THM-588
title: The dominance-reversal tournament metagraph Laplacian has spectrum {2k : mult(k)>0} (mult = THM-587 signed cycle index); since mult(1)=0 (proved) and the cyclicity 3-cycle count is a nonzero level-2 invariant (mult(2)>=1), its ALGEBRAIC CONNECTIVITY is exactly 4 for every n>=3, the arc-flip random walk P=A/d has spectral gap exactly 4/C(n,2), and the Fiedler (slowest) mode is the 3-cycle count
status: PROVED (mult(1)=0 by the transposition argument; mult(2)>=1 since cyclicity is a degree-2 relabeling invariant with nonzero level-2 part; alg.conn = 2*min{k>=1:mult(k)>0} = 4); VERIFIED n=3..9 (mult(1)=0, mult(2)=1, alg.conn=4, gap=4/d) and the Fiedler mode = cyclicity (|corr|=1.000 n=4,5; 0.998 n=6)
source: klein-2026-06-29-S3
depends_on:
  - THM-587   # mult(k) = coeff of the signed cycle index; eigenvalues of A are d-2k
  - THM-584   # complement = antipodal; metagraph = S_n-quotient of Q_d
related:
  - HYP-3552  # the LRC level-2 / pairwise-second-moment transfer this motivates
  - HYP-3543  # one R, three spectra
  - eigenvalues-of-the-merged-metagraph.md  # opus-S268: the DIFFERENT merged-wiggly metagraph (H=2nd eigenvector, gap ~c/n conjectural); this sharpens the canonical arc-flip normalization
results:
  - 04-computation/metagraph-laplacian-spectral-invariants.py
  - 05-knowledge/results/metagraph-laplacian-spectral-invariants.out
---

# THM-588 — Algebraic connectivity of the dominance-reversal metagraph is 4

## Setup

The **arc-flip (dominance-reversal) metagraph** is the canonical `d`-regular weighted graph on the
`A000568(n)` iso classes (`d = C(n,2)`): one edge per single arc reversal (= one hypercube `Q_d` edge),
`A[i][j] = #{arcs whose flip sends a rep of class i into class j}`, row sums `d`. By THM-584/587 its
adjacency eigenvalues are `d-2k` with multiplicity `mult(k) = [x^k] P_n(x)` (the signed cycle index).

## Statement

1. **Laplacian spectrum.** `L = dI - A` has spectrum `{2k : k=0..d}` with multiplicity `mult(k)`.
2. **`mult(1) = 0`** for all `n>=2`: there is no `S_n`-invariant at hypercube level 1.
3. **`mult(2) >= 1`** for all `n>=3`: the cyclicity (3-cycle count) is a degree-2 relabeling invariant
   with nonzero level-2 component. (Computation: `mult(2)=1` exactly for n=3..9 — the Fiedler space is
   simple.)
4. **Algebraic connectivity = 4** for all `n>=3`: the smallest nonzero Laplacian eigenvalue is
   `2 * min{k>=1 : mult(k)>0} = 2*2 = 4`.
5. **Arc-flip walk gap.** The simple random walk `P = A/d` has spectral gap `1 - lambda_2 = 4/d
   = 4/C(n,2)` exactly (relaxation `~ C(n,2)/4`, mixing `O(n^2 log n)`).
6. **Fiedler mode = cyclicity.** The eigenvector at `A`-eigenvalue `d-4` (level 2) is the 3-cycle count
   `c(T) = C(n,3) - sum_i C(s_i,2)` (verified `|corr| = 1.000` n=4,5; `0.998` n=6). It is the
   slowest-equilibrating observable of the dominance-reversal walk.

## Proof

(1) is `L = dI - A` applied to `A`'s spectrum `d-2k` (THM-587).

(2) `mult(1) = dim of S_n-invariants in span{chi_a : a an arc}`. For any arc `a={i,j}`, the transposition
`(i j) in S_n` fixes the pair but reverses its orientation, so it acts as `-1` on `chi_a` and sends no
other arc to `a`. Hence in any invariant `v = sum_a c_a chi_a`, the `chi_{ij}`-component of `(i j)v` is
`-c_{ij}`; invariance forces `c_{ij} = -c_{ij}`, so `c_{ij}=0` for every arc, i.e. `v=0`. Thus `mult(1)=0`.

(3) The cyclicity `c(T)` is invariant under relabeling and is a degree-2 polynomial in the arc variables
(`sum_i C(s_i,2)` is quadratic in the scores, which are linear in the arcs). Its Boolean-Fourier expansion
therefore has a level-2 component, and that component is nonzero (it genuinely depends on pairs of arcs).
Being `S_n`-invariant, its level-2 projection is a nonzero `S_n`-invariant at level 2, so `mult(2) >= 1`.

(4) The metagraph is connected (arc flips connect all tournaments), so `0` is a simple Laplacian
eigenvalue (`mult(0)=1`, the constant). The next eigenvalue is `2 * min{k>=1 : mult(k)>0}`; by (2),(3)
that minimum is `2`, giving algebraic connectivity `4`. (5) is `(min nonzero L-eigenvalue)/d`. (6) the
`d-4` eigenspace is the level-2 invariant space, spanned (for n=3..9) by the cyclicity. ∎

## Significance and relation to prior work

- **A clean, exact spectral invariant.** Unlike the merged-wiggly metagraph of `eigenvalues-of-the-merged-metagraph.md`
  (opus-S268), whose Markov gap was only conjectured `~ c/n` and whose spectrum is semicircular/messy, the
  canonical **arc-flip weighted** metagraph has the EXACT gap `4/d` and integer Laplacian spectrum `{2k}`.
  These are two different normalizations of "the metagraph": S268's slow mode is the H-gradient (`H ≈ 2nd
  eigenvector`, a degree-`(n-1)` observable); the arc-flip metagraph's slow mode is the degree-2 cyclicity.
  Both are real; THM-588 fixes the canonical one and removes the conjecture.
- **No-linear / single-quadratic structure.** `mult(1)=0` (no first-order invariant) and `mult(2)=1` (one
  quadratic invariant) say the lowest nontrivial content of the metagraph is a SINGLE pairwise/quadratic
  mode — the cyclicity. This is the structural template behind the LRC transfer HYP-3552: the binding
  obstruction is purely the second moment, with the first moment forced to vanish.
- **Engineering.** All of `mult(k)` (hence the full Laplacian spectrum, complexity `kappa(n)`, neutral-flip
  trace `tr(A) = sum mult(k)(d-2k) = 2,6,16,60,328,3160,54928` for n=3..9, heat trace, spectral zeta) come
  from `n!` permutations via THM-587 — a closed-form spectral generator for the metagraph at any `n`.

Caveat (MISTAKE-033): `A` is the full arc-flip graph; "arc flip" = reverse one dominance (one `Q_d` edge),
not a tiling-model tile flip.
