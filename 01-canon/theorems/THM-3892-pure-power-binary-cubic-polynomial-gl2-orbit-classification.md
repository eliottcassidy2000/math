---
id: THM-3892
title: "Pure-power binary-cubic carriers are polynomial GL2 orbits"
status: >
  RESERVED / PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT; awaiting
  independent hostile audit.  Over an algebraically closed characteristic-zero
  field, every homogeneous degree-D binary-cubic coefficient row with
  discriminant kappa C^(4D) is, after dehomogenizing C=1, exactly a constant
  squarefree cubic acted on by a polynomial GL2 matrix of constant
  determinant.  Conversely every such orbit has constant discriminant and
  homogenizes to a pure-power carrier.  Since k[t] is Euclidean, the matrix
  has a finite elementary/Cohn shear word.  Thus pure single-support leading
  carriers are intrinsically elementary-word objects; genuinely new leading
  geometry must use several discriminant directions or leave the stated
  field/characteristic hypotheses.
source: jc_sparse_direct_search / all-degree successor to THM-3891, 2026-08-23
audit: >
  SELF-AUDITED proof candidate.  The exact companion verifies binary-cubic
  GL2 covariance with determinant exponent six, scalar exponent four, the
  pair-determinant coordinate reconstruction, the one-branch
  Riemann--Hurwitz inequalities, central-scalar/contact accounting, a
  nontrivial triangular degree-two control, and an alternating three-shear
  degree-eight control.  Normal and -O replays are byte-identical to the
  frozen 146-check transcript.  Independent audit must recheck finite-etale
  trivialization, polynomial-section factorization, and the elementary-word
  corollary's exact source-versus-target scope.
related:
  - THM-3709-cohn-alternating-two-by-two-elementary-decoration-nonentry
  - THM-3721-automorphic-cohn-one-right-shear-nonentry
  - THM-3736-automorphic-cohn-complete-constant-sl2-polynomial-exposure-classification
  - THM-3889-maximally-confluent-quadratic-binary-cubic-two-place-obstruction
  - THM-3891-split-quadratic-c8-pencil-branch-value-two-place-obstruction
script: 04-computation/jc2_pure_power_binary_cubic_polynomial_gl2_orbit_thm3892.py
output: 05-knowledge/results/jc2_pure_power_binary_cubic_polynomial_gl2_orbit_thm3892.out
script_sha256: 298d0f56806e9b7d7d7e750587036c61e57bf06cd8f1f37d8e43cb70aeb4d61e
output_sha256: e9051f3bd5da9ad6db88cbf739830dbbac053509968a075984cb46c3e401ab48
semantic_sha256: f83b74a00542c1c943e5677c993038ea46106a9ae41031dd653fac3eddc1c6b2
hash_basis: raw LF bytes
---

# THM-3892 -- pure-power carriers are exactly polynomial matrix orbits

**RESERVED / PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT; awaiting
independent hostile audit.**  Work over an algebraically closed field `k` of
characteristic zero.  Fix `D>=0`, target coordinates `A,C`, and binary-form
coordinates `U,V`.  Let

```text
Phi_D in k[A,C]_D tensor k[U,V]_3                              (1)
```

be a homogeneous binary cubic whose discriminant is

```text
Disc(Phi_D)=kappa C^(4D),                         kappa!=0.  (2)
```

Dehomogenize in the complementary chart:

```text
phi(t;U,V)=Phi_D(t,1;U,V),                         t=A/C.   (3)
```

Then there are a squarefree constant binary cubic `F_0`, a scalar `c in
k^*`, and a matrix

```text
M(t) in GL_2(k[t]),                    det M in k^*,          (4)
```

such that, for the column `X=(U,V)^T`,

```text
phi(t;U,V)=c F_0(M(t)X).                                  (5)
```

Conversely, every expression `(5)` whose coefficient degree is at most `D`
homogenizes by

```text
Phi_D(A,C)=C^D phi(A/C)                                   (6)
```

and satisfies `(2)`, with

```text
kappa=c^4 det(M)^6 Disc(F_0).                              (7)
```

Thus `(2)` is equivalent to an isotrivial polynomial `GL_2` orbit.  After a
constant rescaling, `M` may be taken in `SL_2(k[t])`; Euclidean row reduction
then writes it as a finite alternating product of upper and lower elementary
shears and a constant matrix.

## 1. The root divisor is a trivial etale cover

Equation `(2)` gives

```text
Disc(phi)=kappa.                                           (8)
```

The zero divisor `Z(phi)` in `P^1_[U:V] times A^1_t` is therefore finite
etale of degree three over `A^1_t`.  Indeed, every fibre has three distinct
points; the divisor is projective and quasi-finite over the base, hence
finite, and the invertible discriminant is the relative etale certificate.

The triviality needed here has a short characteristic-zero proof.  Let a
connected component have degree `r` and normalize its projective completion
over `P^1_t`.  All ramification lies over `t=infinity`.  Riemann--Hurwitz
requires total ramification at least

```text
2r-2,                                                       (9)
```

whereas one fibre contributes at most `r-1`.  Thus `r=1`.  Consequently

```text
Z(phi)=Gamma_0 disjoint_union Gamma_1 disjoint_union Gamma_2, (10)
```

where every `Gamma_i` is the graph of a polynomial section
`A^1 -> P^1`.

## 2. Pair determinants force one common polynomial matrix

Since `Pic(A^1)=0`, write the three sections by primitive polynomial rows

```text
ell_i(t)=(p_i(t),q_i(t)),                 gcd(p_i,q_i)=1.   (11)
```

The two binary cubics with the same relative zero divisor differ by a unit of
`k[t]`, hence by a scalar:

```text
phi=c product_{i=0}^2 (p_iU+q_iV),                c in k^*. (12)
```

The product-of-determinants formula and `(8)` give

```text
kappa=c^4 product_{i<j} det(ell_i,ell_j)^2.                (13)
```

Every determinant in `(13)` is a nonzero polynomial and their product is a
unit.  Therefore

```text
det(ell_i,ell_j) in k^*                       for i!=j.     (14)
```

Take `M(t)` to be the matrix with rows `ell_0,ell_1`.  It lies in
`GL_2(k[t])`.  Write

```text
ell_2=a(t)ell_0+b(t)ell_1.                                (15)
```

The determinant with `ell_0` is `b det M`, and the determinant with `ell_1`
is `-a det M`.  Equation `(14)` forces

```text
a,b in k^*.                                                (16)
```

Thus `(12)` becomes `(5)` with

```text
F_0(W_0,W_1)=W_0 W_1(aW_0+bW_1).                          (17)
```

Both `a,b` are nonzero, so `F_0` is squarefree.  This proves the forward
classification without choosing roots anew after the first etale split.

## 3. Converse and covariance exponents

For a binary cubic `F` and a constant or polynomial matrix `M`, the universal
covariance identity is

```text
Disc(F(MX))=det(M)^6 Disc(F).                              (18)
```

Multiplying the cubic by `c` multiplies its discriminant by `c^4`.  Since
`det M` is a scalar unit, `(5)` has constant discriminant `(7)`.  Formula
`(6)` multiplies each of the four coefficients by the appropriate total
homogenizing power, so its discriminant becomes `C^(4D)` times `(7)`.  This
proves the converse.

Let `E=deg_t phi<=D`.  The exact information lost by dehomogenization is the
central scalar layer

```text
Phi_D=C^(D-E) Psi_E,                                      (19)
```

where `Psi_E` is the primitive homogenization of the orbit.  Its discriminant
contributes `C^(4E)`, while the central scalar contributes the remaining
`C^(4(D-E))`.  This sidecar separates genuine section motion from a common
coefficient vanishing and prevents the pure-power exponent from being
mistaken for root contact alone.

## 4. Elementary-word corollary and its scope

Scale a constant row of `M` so that `det M=1`.  The ring `k[t]` is Euclidean.
Applying the Euclidean algorithm to the first column and clearing the second
by determinant one performs left multiplication by matrices

```text
E_+(v)=[[1,v],[0,1]],                 E_-(u)=[[1,0],[u,1]], (20)
```

with `u,v in k[t]`, until a constant matrix remains.  Reversing the row
operations expresses `M` as a finite elementary word.  Hence every pure-power
carrier is obtained from one constant cubic by a finite polynomial Cohn word.

This is a **source-binary** matrix statement.  The matrix acts on `(U,V)` and
depends on the dehomogenized target parameter `t`; it is not automatically a
polynomial symplectic target transformation of a planar Keller pair.  On
rehomogenizing it may carry a pole at `C=0`, and the central layer `(19)` is
lost.  The theorem identifies the leading carrier grammar, not a JC(2)
equivalence or counterexample.

For `D=2`, THM-3891 analyzes the two possible low-degree word profiles and
proves that every irreducible completed discriminant still has at least two
normalization places at infinity.  In higher degree, the matrix word is now
the exact search object: its alternating depth, degree profile, central
layer, and the lower-order perturbation determine whether the two-place
obstruction persists.

## 5. Boundaries

Algebraic closedness and characteristic zero are essential to the stated
split form.  Over a nonclosed field the three sections can carry finite-etale
Galois descent, and in positive characteristic Artin--Schreier covers make
`A^1` nontrivially etale-covered.  A leading discriminant with two or more
projective support lines also lies outside `(2)`.  Higher-order completion,
the polynomiality of a full Keller atlas, nonproperness, and JC(2) remain
**OPEN**.

Reproduce the exact packet with

```bash
python3 04-computation/jc2_pure_power_binary_cubic_polynomial_gl2_orbit_thm3892.py
python3 -O 04-computation/jc2_pure_power_binary_cubic_polynomial_gl2_orbit_thm3892.py
```

Both streams must byte-match
`05-knowledge/results/jc2_pure_power_binary_cubic_polynomial_gl2_orbit_thm3892.out`.
