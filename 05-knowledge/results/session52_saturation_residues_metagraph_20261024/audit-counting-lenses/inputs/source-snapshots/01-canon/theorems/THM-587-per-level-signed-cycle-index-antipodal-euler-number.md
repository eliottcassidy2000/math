---
id: THM-587
title: The per-level metagraph multiplicities are the PER-LEVEL SIGNED CYCLE INDEX P_n(x) = (1/n!) sum_sigma prod_cycles (1 + s_c x^{ell_c}) of the vertex-induced signed S_n action on the C(n,2) arcs; its two antipodal evaluations are P_n(1)=A000568(n) (the all-odd-cycle tournament Burnside) and P_n(-1)=SC(n) (self-converse count = the ANTIPODAL EULER/LEFSCHETZ number); it computes the full metagraph spectrum from n! permutations, past the 2^{C(n,2)} enumeration wall
status: PROVED (Molien/Burnside for the level-graded character + Lefschetz fixed-point count of the complement involution); VERIFIED exactly -- reproduces the n<=6 block spectra (cross-check MATCH), P_n(1)=A000568 and P_n(-1)=SC=2,2,8,12,88,176 for n=3..8 (SC cross-checked vs CLAUDE.md NS-merged 0,1,2,22,184)
source: klein-2026-06-29-S2
depends_on:
  - THM-584   # complement = antipodal map; metagraph = S_n-quotient of Q_d; mult(k)=dim invariants at level k
related:
  - HYP-3540  # this CLOSES it (the level-multiplicity sequence = this signed cycle index)
  - HYP-3543  # mac-mini: metagraph spectrum = LRC cap spectrum (one R); identified the signed-Burnside nature
  - HYP-3538  # the R +-1 organizing principle
  - HYP-2080  # self-complementary / Burnside counts
  - HYP-3544  # the Kaczynski/computational-homology refinement of P_n(-1) into Betti numbers
results:
  - 04-computation/signed-cycle-index-metagraph-spectrum.py
  - 05-knowledge/results/signed-cycle-index-metagraph-spectrum.out
---

# THM-587 — The per-level signed cycle index, and its two antipodal evaluations

## Statement

Let `d = C(n,2)`. `S_n` acts on the `d` unordered vertex-pairs as a **signed** permutation group
(the vertex-induced subgroup of the hyperoctahedral `B_d`): a swap that reverses a pair flips its
orientation bit. For `sigma in S_n`, decompose its action on the pairs into signed cycles; a cycle of
length `ell_c` has orientation sign `s_c in {+1,-1}` (the product of the per-step order-reversals around
the cycle). Then the metagraph eigenvalue multiplicities of THM-584 (`mult(k) = dim of S_n-invariants at
hypercube level k = multiplicity of eigenvalue d-2k`) are the coefficients of the **per-level signed
cycle index**

```
P_n(x) = sum_{k=0}^{d} mult(k) x^k  =  (1/n!) * sum_{sigma in S_n}  prod_{cycles c of sigma}  (1 + s_c x^{ell_c}).
```

Its two evaluations at the antipodal fixed signs are:

1. **`P_n(1) = A000568(n)`** — the number of iso classes. At `x=1` a cycle factor is `2` if `s_c=+1`
   and `0` if `s_c=-1`, so only `sigma` whose pair-cycles are all orientation-preserving survive,
   contributing `2^{#cycles}` — exactly the classical all-odd-cycle tournament Burnside (CLAUDE.md:
   "Fix=0 for even cycles, 2^{orbit-pairs} for all-odd").
2. **`P_n(-1) = SC(n)`** — the number of **self-converse** (self-complementary) tournament classes.
   This is the **Lefschetz fixed-point count** of the complement involution `R` on the iso-class space:
   `P_n(-1) = sum_k (-1)^k mult(k) = dim V_+ - dim V_- = trace(R) = #{classes fixed by R} = SC(n)`.
   It is the **equivariant Euler characteristic** of the level-graded invariant complex (the antipodal
   Euler number), since `R` acts as `(-1)^k` on level `k` (THM-584).

Consequently `dim V_+ = (P_n(1)+P_n(-1))/2 = (A000568+SC)/2 = V_merged` (the merged metagraph
`G_n/Z_2`) and `dim V_- = (P_n(1)-P_n(-1))/2 = (A000568-SC)/2 = #NS pairs`.

## Proof

By THM-584, `mult(k)` is the dimension of `S_n`-invariants in the level-`k` Boolean-Fourier subspace
`W_k = span{chi_S : |S| = k}` of functions on `Q_d`. By Burnside/Molien,
`sum_k mult(k) x^k = (1/n!) sum_sigma (sum_k chi_k(sigma) x^k)`, where `chi_k(sigma) = trace of sigma on
W_k`. A `sigma`-invariant `S` (with `sigma(S)=S`) is a union of whole pair-cycles, and `sigma chi_S =
eps(sigma,S) chi_S` with `eps(sigma,S) = prod_{cycles c subset S} s_c` (each included cycle contributes
its orientation sign). Hence `sum_k chi_k(sigma) x^k = prod_c (1 + s_c x^{ell_c})` (each cycle is either
excluded, factor `1`, or included, factor `s_c x^{ell_c}`). Averaging over `S_n` gives `P_n(x)`.
(NB: this is the SYMMETRIC-subset trace; it differs from the exterior `det(I + x M_sigma)` by reordering
signs — only this version matches the verified spectrum.)
`P_n(1)`: Molien at `x=1` is `(1/n!) sum_sigma dim(W^sigma) = #orbits on Q_d = A000568(n)`.
`P_n(-1) = sum_k (-1)^k mult(k)`; since `R = antipodal` acts as `(-1)^k` on level `k`, this is
`trace(R)` on the class space, whose value for an involution is the number of fixed points `= SC(n)`. ∎

## Verification (exact)

| n | P_n(1)=A000568 | P_n(-1)=SC(n) | metagraph multiplicities mult(k), k=0.. |
|---|----|----|----|
| 4 | 4 | 2 | 1,0,1,1,1,0,0  (MATCHES verified n=4 block spectrum) |
| 5 | 12 | 8 | 1,0,1,1,4,1,3,0,1,0,0  (MATCHES) |
| 6 | 56 | 12 | 1,0,1,1,5,5,10,8,12,6,4,2,1,0,0,0  (MATCHES) |
| 7 | 456 | 88 | 1,0,1,1,5,6,22,23,55,56,85,61,69,31,28,6,5,0,1,0,0,0  (NEW) |
| 8 | 6880 | 176 | 1,0,1,1,5,6,26,44,115,197,390,562,825,941,1027,917,755,504,310,154,66,24,6,2,1,0,0,0,0 (NEW) |

`SC(n) = 2,2,8,12,88,176` (n=3..8) cross-checks against CLAUDE.md's independently-computed NS-merged
sequence `0,1,2,22,184` via `SC = A000568 - 2*NS` (matches at n=3..7). SC is the classical self-converse
tournament count (Eplett); n=7,8 values extend the project's table.

## Significance

- **CLOSES HYP-3540.** The level-multiplicity sequence is exactly this signed cycle index — not a
  graph-by-edges OEIS row (those sum to A000088; the bit-flip makes the action signed, summing to
  A000568). mac-mini (HYP-3543) identified the signed-Burnside nature; this gives the closed generating
  function and proves both evaluations.
- **Breaks the enumeration wall.** `P_n` needs only `n!` permutations, not `2^{C(n,2)}` tournaments, so
  it delivers the FULL metagraph spectrum (eigenvalues `d-2k` with multiplicities) at `n=7` (`2^21`),
  `n=8` (`2^28`), and beyond — where direct metagraph construction is infeasible. Engineering: a
  closed-form spectral generator for the tournament dominance-reversal metagraph.
- **`P_n(-1) = SC(n)` is the antipodal Euler number.** Complement `R` is free on the `2^d` labeled
  tournaments (`x != x XOR 1`), so `(Q_d, R)` is a free `Z_2`-space (the Borsuk-Ulam setting); the
  self-converse census is its equivariant Euler characteristic. This is the entry point for the
  Ky-Fan / computational-homology refinement (HYP-3544): refine `chi = SC(n)` into the `Z_2`-equivariant
  Betti numbers of the level-graded complex; the R-odd Betti numbers are the Borsuk-Ulam obstruction.
- **The two evaluations unify two Burnsides.** `P_n(1)` is the tournament Burnside (all-odd-cycle), and
  `P_n(-1)` is the self-converse Burnside, as the `x = +1 / -1` antipodal specializations of one signed
  cycle index. The whole metagraph spectrum interpolates between them, level by level.

See the reflection `the-per-level-signed-cycle-index-borsuk-ulam-ky-fan.md`. Caveat (MISTAKE-033):
`R` = reverse ALL arcs = antipodal of the full `Q_d`, not the tiling complement (flip tiles, fix base path).
