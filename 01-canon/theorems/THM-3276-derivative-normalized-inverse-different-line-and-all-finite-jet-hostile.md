---
id: THM-3276
title: "Derivative-normalized inverse-different line and all-finite-jet hostile"
status: >
  PROVED + VERIFIED-EXACT.  For a chosen monic primitive quartic coordinate
  xi with derivative star D=f_T(xi), the Keller cofactors form exactly the
  one-dimensional K-line L_(f,xi)=K D^(-1) inside the derivative-typed
  inverse different.  At a 1+3 factorization A=K x L, a cofactor lies on
  this line iff its physical Jacobian cD is diagonal, iff the pointed ratio
  rho is one, iff Norm_(L/K)(rho-1)=0.  The norm-one gauge
  c_lambda=(lambda^(-3),lambda)D^(-1) meets the line exactly when
  lambda^4=1.  More sharply, for every N the choice lambda=1+t^N preserves
  the total norm, all valuations and leading units, and the complete
  inverse-spectral numerator and cofactor jets through order N-1, yet has
  nonzero Keller defect of order v(Delta)=3N.  Thus no bounded finite-jet or
  associated-graded refinement of the THM-3275 conductor clock can recognize
  Keller incidence.  The exact diagonal line is the smallest K-linear
  condition, but attaching the true affine chain-rule cofactor to it remains
  the JC gate.  No polynomial-map hostile realization, C3 exclusion, or
  JC(2) theorem follows.
source: jc-inverse-different-line-2026-08-02
depends_on:
  - THM-3275-off-resonant-centered-packet-covariant-order-index-clock
  - THM-3274-graph-quartic-centered-packet-fixed-decoder-and-cofactor-independence
  - THM-3064-pointed-cubic-norm-keller-decoder-and-inverse-different-boundary
  - THM-2621-planar-degree-four-inverse-spectral-keller-congruence-and-sheet-defect-pole-ledger
related:
  - THM-3038-split-monogenic-order-cross-resultant-conductor-and-affine-owner-boundary
  - THM-3066-k4-initial-face-product-quotient-blind-to-keller-sheetwise-cofactor
script: 04-computation/jc_inverse_different_diagonal_line_finite_jet_hostile_thm3276.py
output: 05-knowledge/results/jc_inverse_different_diagonal_line_finite_jet_hostile_thm3276.out
script_sha256: 9e12a51eb119ec0ce0a6d7ad9c6ff07fca993953343862b50b82c0cc2630c4b3
output_sha256: d0082d4cd0d7a315a17b8fa51404d7edef20abffbf3b12db015720366e4643ac
hash_basis: LF-normalized bytes
---

# THM-3276 -- Keller incidence is one exact inverse-different line

**PROVED + VERIFIED-EXACT.**

## 1. Type the line by the chosen primitive coordinate

Let `K` be a characteristic-zero field and let

```text
f(T) in K[T] be monic, separable, and quartic,
A=K[T]/(f),                         xi=T mod f.          (1)
```

The chosen primitive coordinate supplies its derivative star

```text
D=f_T(xi) in A*.                                        (2)
```

Over an integral monogenic order, `D^(-1)A` is the codifferent.  Over the
fraction algebra it remains the correctly typed derivative-trivialized
inverse-different module.  Define the **diagonal inverse-different line**

```text
L_(f,xi)=D^(-1)(K*1_A) subset D^(-1)A.                 (3)
```

This notation is intentionally not `L_f` alone.  The line depends on the
chosen primitive graph coordinate through `(2)`.  If that coordinate is
changed, the derivative star and primitive-element cofactor must transform
reciprocally; only their physical product is coordinate-independent.

Let `c in A*` be a primitive-element cofactor and put

```text
J=cD in A*.                                              (4)
```

Then directly from `(3)`,

```text
c in L_(f,xi)  iff  J in K*1_A.                         (5)
```

Thus `(3)` is not a decorative projectivization.  It is exactly the set of
cofactors whose physical Jacobian packet is constant on all four sheets.
It is also minimal in a literal linear sense: every possible nonzero Keller
constant `kappa` gives `c=kappa D^(-1)`, so any `K`-linear subspace containing
all Keller cofactors contains `(3)`; conversely `(3)` consists only of those
cofactors.

## 2. The pointed cubic equation of the line

Assume now a local `1+3` factorization

```text
A=K x L,                    [L:K]=3,
D=(D_0,D_C),                c=(c_0,c_C),                (6)
```

with `L` a field.  The physical packet and pointed ratio are

```text
J_0=c_0D_0 in K*,
J_C=c_CD_C in L*,
rho=J_C/J_0 in L*.                                     (7)
```

Equations `(3)--(7)` give

```text
c in L_(f,xi)
 iff J_C=J_0
 iff rho=1.                                              (8)
```

Because `L` is a field, THM-3064 supplies the one-scalar descent

```text
Delta_J=Norm_(L/K)(rho-1),
c in L_(f,xi) iff Delta_J=0.                            (9)
```

The field hypothesis is load bearing.  For a split cubic product, a zero
norm detects one vanishing component, not diagonal equality; there the
element `rho-1` or its three coordinates must be retained.

If `q=c^(-1)` is the inverse-spectral numerator, `(3)` is equivalently

```text
q in K*D,               [q]=[D] in P(A).               (10)
```

For the true affine inverse pair of THM-2621,

```text
q(T)=f_v(T)b_u(T)-f_u(T)b_v(T) mod f,                  (11)
```

and `(10)` is exactly the chain-rule Keller congruence

```text
q == kappa^(-1)f_T mod f.                              (12)
```

Consequently `(3)`, `(9)`, and `(10)` are three forms of the same object:
an inverse-different line, a pointed cubic norm equation, and an affine
inverse-spectral incidence.

## 3. The determinant-one gauge and its exact intersection

Let `e` be the fixed-sheet idempotent decoded in THM-3274 and begin from the
unit-Jacobian cofactor `D^(-1)`.  For `lambda in K*`, define

```text
u_lambda=lambda^(-3)e+lambda(1-e),
c_lambda=u_lambda D^(-1).                              (13)
```

Its norm is one:

```text
Norm_(A/K)(u_lambda)=lambda^(-3)lambda^3=1.             (14)
```

Thus `(13)` preserves the total cofactor norm and, for a valuation-unit
`lambda`, every cofactor valuation.  It also leaves the root quartic,
centered packet, fixed projector, collision covariant, and THM-3275 order
clock untouched.  But

```text
J_lambda=c_lambda D=(lambda^(-3),lambda),
rho_lambda=lambda^4.                                   (15)
```

By `(5)` or `(8)`,

```text
c_lambda in L_(f,xi)
 iff lambda^(-3)=lambda
 iff lambda^4=1.                                        (16)
```

The residual fourth roots in `(16)` are not failures of the line test.  If
`lambda^4=1`, both physical components equal `lambda`, so the packet is
genuinely Keller with constant `lambda`.  Projectively, the norm-one gauge
orbit meets the single point `[D^(-1)]` exactly at these equivalent scalar
representatives.

## 4. An exact integral polynomial model

To test more than valuations, let `k` be a characteristic-zero field
containing `zeta_3`, and set

```text
K=k((t)),                L=K[y]/(y^3-t),
g(T)=T^3-t,              f(T)=(T-1)g(T),
d=g(1)=1-t.                                               (17)
```

The fixed idempotent and derivative star are

```text
e_0=(T^3-t)/(1-t) mod f,
D=f_T(T) mod f.                                         (18)
```

Both are integral in the split order because `1-t` is a unit.  For an
arbitrary base unit `lambda`, define the degree-at-most-three polynomial

```text
q_lambda(T)
 =lambda^(-1)f_T(T)
  +(lambda^3-lambda^(-1))d e_0(T).                      (19)
```

Its fixed and cubic values are exactly

```text
q_lambda(1)=lambda^3D_0,
q_lambda(y)=lambda^(-1)D_C.                             (20)
```

Therefore `c_lambda=q_lambda^(-1)` is precisely the gauge `(13)`.  In
particular, the hostile is already represented by an honest reduced
quartic-algebra polynomial; what is not asserted is that `(19)` arises from
an affine companion `b` through `(11)`.

The THM-3064 resultant computes the line equation without choosing `y`:

```text
 Res_T(g(T), q_lambda(1)(T-1)g'(T)-d q_lambda(T))
 ---------------------------------------------------
 d^3 Res_T(g(T),q_lambda(T))
   =(lambda^4-1)^3.                                    (21)
```

Thus the exact scalar cutting the diagonal line agrees with `(15)--(16)`.

## 5. Every bounded finite jet fails

Fix any integer `N>=1` and specialize

```text
lambda_N=1+t^N.                                         (22)
```

All four branch scales are congruent to one modulo `t^N`:

```text
lambda_N^3 ==lambda_N^(-1)
 ==lambda_N^(-3)==lambda_N ==1 mod t^N.                (23)
```

Equations `(19)--(23)` imply the module-typed congruences

```text
q_(lambda_N) ==D mod t^N A,
c_(lambda_N) ==D^(-1) mod t^N D^(-1)A.                 (24)
```

Hence the twisted packet agrees with the unit-Keller packet in

- every valuation;
- every leading residue and normalized leading unit;
- the full reduced inverse-spectral polynomial through order `N-1`;
- the full cofactor jet through order `N-1`;
- total norm, root packet, collision covariant, and order-index clock.

Nevertheless characteristic zero and residue characteristic different from
two give

```text
rho_N-1=(1+t^N)^4-1=4t^N+O(t^(2N)),
v(rho_N-1)=N,                                            (25)

Delta_N=Norm_(L/K)(rho_N-1)=(rho_N-1)^3 !=0,
v(Delta_N)=3N.                                          (26)
```

Given any proposed finite jet depth `R`, choose `N>R`.  Then `(24)` makes
the hostile indistinguishable from the Keller line through depth `R`, while
`(26)` proves exact nonincidence.  This refutes, uniformly, every criterion
which uses only a bounded associated graded, any finite list of leading-unit
corrections, or a fixed truncation of the cofactor Laurent series.

The hostile is stronger than merely saying that valuations fail: the
Keller defect may be postponed arbitrarily far without changing any of the
root-side mathematics.

## 6. What the affine gate now is

The local information hierarchy is exact:

```text
root packet + C:
  fixed component, first nonbase depth, and order/different valuations;

all bounded cofactor jets:
  arbitrarily accurate approximation to the diagonal line;

exact affine datum:
  the true q from (11) and exact incidence [q]=[f_T].    (27)
```

Thus the smallest exact **linear** datum is the typed line `(3)`, or its
single pointed equation `(9)` after the cubic factor is known to be a field.
But root and order invariants do not identify the true cofactor produced by
the same affine polynomial graph.  One must still prove that the element
from `(11)` has, or cannot have, incidence `(10)` while satisfying the
global polynomial and boundary constraints.

That recognition/realization step is the remaining JC gate.  For an actual
Keller map, the chain rule puts `q` on `(10)` automatically; progress must
therefore make that exact line incompatible with nonproper graph geometry,
not try to distinguish it using another bounded local jet.  The family
`(22)` proves that no such bounded local refinement can suffice.

## 7. Exact reproduction

Run

```bash
python3 04-computation/jc_inverse_different_diagonal_line_finite_jet_hostile_thm3276.py
python3 -O 04-computation/jc_inverse_different_diagonal_line_finite_jet_hostile_thm3276.py
```

Both modes must reproduce
`05-knowledge/results/jc_inverse_different_diagonal_line_finite_jet_hostile_thm3276.out`
byte for byte.  The companion verifies the two branch evaluations of
`q_lambda`, determinant-one norm, fourth-power line intersection, exact
THM-3064 resultant, sixteen increasingly deep finite-jet hostiles, all jet
and defect valuations, and a fifth-order explicit control.
