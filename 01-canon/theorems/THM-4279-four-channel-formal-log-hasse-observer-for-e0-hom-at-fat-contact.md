---
id: THM-4279
title: "Four-channel formal-log Hasse observer for E0-Hom at the fat contact"
status: >
  PROVED RELATIVE TO THM-4272 + FORMAL-LOCAL. After translating a morphism
  C_0->E_0 by its value at Q_epsilon and applying the elliptic formal
  logarithm, the four b-adic coefficients in degrees 1,2,4,7 recover its
  complete Hom class. Equivalently, with the target-value sidecar, these four
  normal-Hasse channels have the same fibres as the full length-twelve
  restriction on global curve morphisms. They form a rank-four linear
  observer on the chosen-embedding
  Hom space, and four scalar linear channels are minimal there. The observer
  is not faithful on arbitrary maps from the fat contact: exp_E0(b^3) is a
  sharp sidecar hostile. No new incidence, raw Keller descent, exact-M=12
  entry, JC(2), or DC(2) is proved.
source: root/cross-frontier-bridge/2026-08-27
depends_on:
  - THM-4272-lambda-zero-a23-contact-and-e0-infinity-jet-obstruction
related:
  - THM-4274-confined-confluent-observer-fibre-and-density-transport
  - THM-4275-opposite-parity-attachment-observers-and-confluent-sample-matroids
audit: >
  Direct corollary audit. THM-4272 supplies the complete chosen-embedding
  pullback-differential space and its four distinct b-orders 0,1,3,6.
  Successive leading-coefficient elimination makes the corresponding four
  coefficient functionals triangular and injective. The characteristic-zero
  elliptic formal logarithm converts them to map channels 1,2,4,7 and turns
  target addition into subtraction of observer vectors. Dimension gives the
  linear four-channel lower bound; exp_E0(b^3) checks the indispensable
  global-Hom sidecar.
---

# THM-4279 -- four-channel formal-log Hasse observer for `E_0`-Hom at the fat contact

**PROVED RELATIVE TO THM-4272 + FORMAL-LOCAL. NO NEW INCIDENCE IS EXCLUDED;
`JC(2)` AND `DC(2)` REMAIN OPEN.**

## 1. Statement

Work over `C` with THM-4272's curve, target, point, and toric uniformizer

```text
C_0:x^6+y^4=1,                 E_0:Y^2=X^3+1,
Q=Q_epsilon,                   b=1/S at Q.             (1)
```

Put `omega_E=dX/Y`. For a curve morphism `m:C_0->E_0`, translate its target
value to the origin,

```text
mbar=m-m(Q),                   mbar(Q)=O.               (2)
```

Let `log_E` be the characteristic-zero formal logarithm at `O`, normalized
by

```text
d log_E=omega_E,
```

and write

```text
L_m(b)=log_E(mbar(b)) in b C[[b]].                     (3)
```

Define the four-channel normal-Hasse observer

```text
H(m)=([b]L_m,[b^2]L_m,[b^4]L_m,[b^7]L_m) in C^4.     (4)
```

> **Theorem.** For all morphisms `m,n:C_0->E_0`,
>
> ```text
> m=n
> iff
> m(Q)=n(Q) and H(m)=H(n).                             (5)
> ```
>
> Thus, after retaining the target-value sidecar, the four channels in `(4)`
> have exactly the same fibres as the complete morphism and as its full
> scheme-theoretic restriction to `12Q`. In fact they factor through the
> restriction to `8Q`.

Let `O_E=Z[omega]` be the Eisenstein endomorphism ring and choose the complex
embedding through which `O_E` acts on `omega_E`. If

```text
M=Hom(J(C_0),E_0),             M_sigma=M tensor_(O_E,sigma) C,    (6)
```

then `(4)` extends to a complex-linear isomorphism from the four-dimensional
space `M_sigma` onto `C^4`. Consequently fewer than four complex-linear
scalar channels cannot be a complete observer on `M_sigma`.

The chosen-embedding tensor product in `(6)` is essential. This theorem does
not identify it with `M tensor_Z C`, which has the wrong dimension and is the
type error already excluded in THM-4272.

## 2. Filtered pullback-differential lemma

We use the following elementary filtered-space fact. Let `V` be a subspace of
`C[[b]] db` with a basis `eta_0,...,eta_(s-1)` satisfying

```text
ord_b eta_0=e_0 < e_1 < ... < e_(s-1).                (7)
```

Then

```text
J:V->C^s,
eta |-> ([b^(e_0)]eta/db,...,[b^(e_(s-1))]eta/db)     (8)
```

is injective. Indeed, the first coordinate in `(8)` kills the coefficient of
`eta_0`; after that, the second kills the coefficient of `eta_1`, and so on.
This is triangular leading-term elimination, not a generic-point argument.

THM-4272 proves that the complete pullback-differential space for maps to
`E_0`, through the chosen CM embedding, is

```text
V=span_C {dx/y^3, y dx/y^3, x y dx/y^3, y^2 dx/y^3}, (9)
```

and that these four displayed vectors have respective `b`-orders

```text
6,3,1,0.                                               (10)
```

Reordering `(9)` by increasing order, `(7)` is therefore

```text
(e_0,e_1,e_2,e_3)=(0,1,3,6).                         (11)
```

The four functionals in `(8)` are injective on `V`. Postcomposition by
`alpha in O_E` multiplies the pullback differential by `sigma(alpha)`, so the
pullback map extends over exactly the tensor product `(6)`. Since THM-4272's
full `O_E`-basis has rank four, it identifies `M_sigma` with `V`: injectivity
also follows directly because a characteristic-zero homomorphism with zero
differential is zero. This proves the rank-four and linear-minimality
assertions.

## 3. Formal logarithm and fibre equality

For every based morphism `(2)`, translation invariance of `omega_E` gives

```text
dL_m=mbar^*(omega_E)=m^*(omega_E).                    (12)
```

For `e in {0,1,3,6}` one has

```text
[b^e](dL_m/db)=(e+1)[b^(e+1)]L_m.                    (13)
```

The factors `e+1` are nonzero over `C`. Combining `(8)`, `(11)`, and `(13)`
proves that `H(m)=0` forces `m^*(omega_E)=0`. A morphism of smooth complex
curves with zero pullback of a nonzero invariant differential is constant;
the based condition makes `mbar` the zero map, so `m` itself is constant.

The formal logarithm is a homomorphism for the elliptic formal group. Hence

```text
overline{m-n}=(m-n)-(m(Q)-n(Q))=mbar-nbar,
L_(m-n)=L_m-L_n.                                      (14)
```

If the right side of `(5)` holds, the based difference `m-n` has zero
observer and is therefore zero. The reverse implication is immediate. This
proves `(5)`.

All exponents in `(4)` are less than eight, so the observer depends only on
the restriction to

```text
8Q=Spec C[b]/(b^8).                                   (15)
```

In particular it factors through THM-4272's `12Q` contact. Equality on the
full fat contact implies equality of the four channels, but the proof uses
only the four labelled coefficients, not twelve generic point values and not
a transverse `W`-jet.

## 4. Sharp boundary: the Hom sidecar cannot be dropped

The four-channel compression is faithful because global `C_0->E_0` maps
land in the four-dimensional subspace `(9)`. It is false for arbitrary local
maps from the fat contact. Let `exp_E` be the inverse of `log_E`, the formal
exponential, and
consider over `C[b]/(b^8)` the based local map

```text
s(b)=exp_E(b^3).                                       (16)
```

It is nonconstant, but

```text
log_E(s(b))=b^3,
H(s)=(0,0,0,0).                                       (17)
```

Thus `(4)` is not a complete observer on the full contact algebra or on all
local `8Q->E_0` maps. THM-4275's full normal-Hasse warning remains correct:
the sparse observer becomes lossless only after retaining membership in the
global Hom subspace. The precise connection ledger is

```text
source:              global C_0->E_0 morphisms modulo translation;
target:              four formal-log Hasse coefficients;
map:                 m |-> degrees 1,2,4,7 of log_E(m-m(Q));
preserved predicate: equality of the complete Hom class;
destroyed data:      target value and arbitrary local-contact directions;
restoring sidecars:  m(Q), global-Hom membership, b, and omega_E;
cheapest hostile:    exp_E(b^3).                       (18)
```

## 5. Scope

This theorem compresses an already excluded special-fibre incidence; it does
not exclude a new degree-`34/42` class. It supplies no equations for the
relative Hom scheme, no effective equation or radius for THM-4272's collar,
no lift of `(4)` away from `W=0`, and no descent of a resolved Keller map to
the raw `A_23` contact. Exact-`M=12` entry, `JC(2)`, and `DC(2)` remain open.
**QED.**
