---
id: THM-3023
title: "The Newton-ratio transform as a dynamical system: Phi_6 neutrality and the metallic discriminant eigenvalue"
status: PROVED + VERIFIED-EXACT / AWAITING INDEPENDENT HOSTILE AUDIT
source: klein-S428
depends_on:
  - THM-3010-ballot-column-newton-ratios-and-metallic-alternation
related:
  - THM-3000-fixed-edge-cumulant-curvature-universality-and-bounded-jet-transfer
  - THM-3004-circuit-sign-change-cluster-law-and-classifier-refutation
  - THM-2997-first-gap-wall-stripped-all-width-second-edge-circuit-positivity
script: 04-computation/newton_ratio_transform_dynamics_thm3023.py
output: 05-knowledge/results/newton_ratio_transform_dynamics_thm3023.out
script_sha256: 82e111c082c0809ec80ac83a4c83e2d46dcf43cd6332b59ae82a459e3705b5f0
output_sha256: 9e5c6220d0593efe636c4f8da5596a39fd38ea5e71f10a63b3a655a1fc4c0eee
hash_basis: LF-normalized bytes
---

# THM-3023 -- the Newton-ratio transform as a dynamical system

**PROVED + VERIFIED-EXACT.**

A new object.  On positive sequences define the **Newton-ratio transform**

    T(h)_k = h_k^2 / (h_(k-1) h_(k+1)).                            (1)

Its output is the Newton ratio `R_k` of THM-2997 (1), so the entire circuit lane
(THM-3000 through THM-3015) studies **one application** of `T`.  Iterating it is
new, and the iteration has clean structure.

## 1. `T` is minus the second difference in log coordinates

`log T(h)_k = -Delta^2(log h)_(k-1)`, hence

    T^m(h)_k = prod_(i=0)^(2m) h_(k-m+i)^((-1)^(m+i) binom(2m,i)),   (2)

multiplicative iterated binomial differencing.  Verified exactly against direct
iteration for `m=1..4` on Catalan.  (The circuit of THM-3004 is one more
difference: `c_k = Delta(log T(h))`.)

## 2. Fixed points are exactly the Lyness period-6 sequences

`T(h)=h` iff `h_(k-1)h_(k+1) = h_k`, i.e.

    h_(k+1) = h_k / h_(k-1),                                        (3)

whose orbits are **6-periodic**: `a, b, b/a, 1/a, 1/b, a/b, a, ...`  Verified on
three seeds.  This is the Lyness/cluster-periodicity recursion.

## 3. Eigensequences, and `Phi_6`

With `log h_k = A x^k`, `T(h) = h^mu` where

    mu(x) = 2 - x - 1/x,    eigen-condition   x^2 + (mu-2)x + 1 = 0. (4)

The eigen-condition is a **reciprocal** quadratic (constant term `+1`, root
product `+1`) -- exactly complementary to THM-3010's metallic family
`x^2 - nx - 1` (root product `-1`).  At `mu = 1` it is

    x^2 - x + 1 = Phi_6(x),                                         (5)

the **sixth cyclotomic polynomial**, whose roots are the primitive 6th roots of
unity -- which is section 2's period 6, arrived at independently.

## 4. Spectrum: the discrete Laplacian

On `|x| = 1`, `x = e^(i theta)`,

    mu = 2 - 2 cos theta = 4 sin^2(theta/2)  in  [0,4],             (6)

the discrete-Laplacian dispersion relation.  The **neutral frequency is
`theta = pi/3`**: log-modes slower than that contract under `T`, faster ones
expand, and `pi/3` is fixed.

## 5. The metallic Simson eigenvalue IS the discriminant

For the metallic recurrence `a_(k+1) = n a_k + a_(k-1)`, with `lambda` the root
of `x^2 - nx - 1`, the Simson/alternating mode is `x = -lambda^(-2)`.  Since
`lambda^2 + lambda^(-2) = n^2 + 2`,

    mu = 2 + lambda^2 + lambda^(-2) = n^2 + 4 = disc(x^2 - nx - 1).  (7)

    golden  5,   silver  8,   BRONZE 13,   copper 20.

The `13` is the same `13` as the `sqrt(13)` in THM-3010's bronze log-concavity
threshold `(3+sqrt13)/2`, and `5 = disc(phi)`.  Numerically the growth ratio of
`|log T^m|` is `5.0000` for golden and `8.0000` for silver.

> **This restates THM-3010 section 3 as an eigenvalue statement.**  "Metallic
> recurrences attain maximal circuit alternation" becomes: the alternation is the
> **unstable eigendirection** of `T`, with growth rate exactly the discriminant.

## 6. Orbits

Hypergeometric sequences have smooth (low-frequency) logs and **contract** to the
fixed point `1`: measured `max|log T^m(h)|` falls monotonically for Catalan,
central binomial and factorial.  Fibonacci **expands**, at rate `5`.

## 7. Boundaries, and one flag

- (2), (3), (4), (5), (7) are identities, proved.  (6) is the standard symbol
  computation.  Section 6 is FINITE-EXACT on the displayed window: contraction is
  observed for `m<=4`, not proved to persist.
- The eigensequences `log h_k = A x^k` are doubly exponential and are **not**
  classical combinatorial sequences; they are the spectral coordinates, not
  examples.  Classical sequences are superpositions dominated by the `x=1`
  (`mu=0`) mode.
- `T` is not injective (any geometric factor is killed), so "orbit" means forward
  orbit only.
- **Flag, not a claim:** `Phi_6` is the neutral locus here, and the LRC lane
  records `Phi_6` resonance as universal in a completely different setting.  I
  have not checked whether these are the same `Phi_6` in any structural sense and
  am not asserting a connection -- only that the coincidence is now on the record
  and cheap for someone to test.
