---
id: THM-3962
title: "Constant-q product curves have a universal ramification-unit obstruction"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. If q=q(P) is
  independent of t and
  T^3-3PT-q(P) is irreducible, the maximal normalization is the cylinder
  over the normalized coefficient curve. A planar Keller open would force
  that curve to be rational. Riemann--Hurwitz then gives at least two
  infinity-or-ramification punctures, so the maximal etale cylinder carries
  a nonconstant unit, impossible on A2. This closes all constant-q members,
  including every nonnormal conductor-debt order. Explicit normalizations
  type the repeated-factor P^2+P and double-P cP^2 controls.
source: jc-degree6-one-place / post-THM-3961 conductor-debt closure, 2026-08-24
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (jc-zero-debt-lift, 2026-08-24). The audit
  independently reconstructed normalization of the cylinder as the
  cylinder over the normalized coefficient curve; checked the dense-plane
  projection and positive-genus exclusion; rederived the degree-three
  Riemann--Hurwitz support bound and its punctured-rational-curve unit
  contradiction; and verified both explicit normalizations, singular
  conductor addresses, ramification packets, unit ranks, and class ledgers.
  Normal and optimized runs byte-match the frozen 46-gate output; raw and
  semantic hashes and documentation checks pass. No repair was required.
depends_on: []
related:
  - THM-3961-arbitrary-q-hidden-repetition-normality-and-conductor-debt
  - THM-3922-affine-plane-open-finite-normal-completion-class-group
script: 04-computation/jc2_constant_q_product_curve_ramification_unit_thm3962.py
output: 05-knowledge/results/jc2_constant_q_product_curve_ramification_unit_thm3962.out
script_sha256: 1babd10dbf3eea74afae74c10c63caa16b0b0865a883b4bc7b18297ed260c1cf
output_sha256: a3265d75a6178d4607a35692d18a7eb5db0fd5f3aa8491e61142bb289502402c
semantic_sha256: 4afd29f2c4a6c0e1c4e97467abce7839474e2785e589269816d6120604e451b6
hash_basis: raw LF bytes
---

# THM-3962 -- every coefficient-constant cubic cylinder is closed

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.** Work over an
algebraically closed field `k` of characteristic zero. Here “constant `q`”
means constant in the cylinder coordinate `t`, not necessarily constant in
`P`. Let

```text
f(P,T)=T^3-3PT-q(P),                    q in k[P],        (1)
C=Spec k[P,T]/(f),
X0=C x A1_t=Spec k[P,t,T]/(f).
```

Assume `f` is irreducible over `k(P)`. Then no same-function-field
nontrivial planar Keller chart has `X0` as its finite cubic order, even when
`X0` is nonnormal. Equivalently, every constant-`q(P)` member of the
arbitrary-`q` grammar of THM-3961 is closed; its two conductor-debt loci do
not survive the maximal normalization.

The proof is an exact product-curve obstruction that applies more generally
to every separable finite curve cover of degree greater than one.

## 1. Normalization commutes with the cylinder

Let

```text
B0=k[P,T]/(f),
Bbar=normalization of B0 in k(C),
Ctilde=Spec Bbar.                                         (2)
```

Since `B0` is an affine curve algebra over a characteristic-zero field, its
normalization is finite. The ring `Bbar[t]` is normal, integral and finite
over `B0[t]`, and has the same fraction field. Conversely, an element of
`Frac(Bbar[t])` integral over `B0[t]` is also integral over `Bbar[t]`,
because its monic equation has coefficients in the smaller ring and hence
in the larger one. Normality of `Bbar[t]` then contains it. Therefore the
full finite normalization is exactly

```text
X=Ctilde x A1_t=Spec Bbar[t].                            (3)
```

The target function `P` induces a finite separable morphism

```text
g:Ctilde -> A1_P,                   degree(g)=3,          (4)
```

and the normalized surface map is `pi=g x id_(A1_t)`.

## 2. A Keller plane would force the normalized curve to be rational

Suppose a same-function-field Keller source existed. The normal source
`A2` factors through `(3)`. Because its map to the target is etale and
therefore quasi-finite, normalization-form Zariski Main identifies it with
a dense open immersion

```text
j:U isomorphic to A2 -> Ctilde x A1.                    (5)
```

Projection of the dense open `U` to `Ctilde` is dominant. After completing
`Ctilde` to a smooth projective curve `Cbar`, this gives a dominant rational
map `P2 --> Cbar`. Restriction to a general projective line is a nonconstant
rational map `P1 --> Cbar`. In characteristic zero, Riemann--Hurwitz forbids
such a map when `g(Cbar)>0`. Hence

```text
Cbar is isomorphic to P1.                                (6)
```

This step closes every positive-genus normalization before any class-group
or conductor calculation is needed.

## 3. The rational case has an unavoidable nonconstant unit

The finite morphism `(4)` extends uniquely to a degree-three morphism

```text
gbar:P1=Cbar -> P1_P.                                    (7)
```

Let `S_inf=Cbar\Ctilde`, and let `R_aff` be the support of the ramification
divisor of `g` on `Ctilde`. The maximal affine etale curve locus is

```text
V=Ctilde\R_aff=P1\S,                  S=S_inf union R_aff. (8)
```

Riemann--Hurwitz gives

```text
sum_(x in P1) (e_x-1)=2*3-2=4.                           (9)
```

At one source point the contribution is at most `3-1=2`; thus the
ramification support of `gbar` contains at least two points. Every such
point is either in `R_aff` or in `S_inf`, so

```text
|S| >= 2.                                                (10)
```

The punctured rational curve `P1\S` consequently has a nonconstant unit:
after choosing two points of `S`, a ratio of their linear forms has divisor
supported on `S`. Let `u` be such a unit. Then

```text
u in O(V)^*\k^*,                 pr_1^*u in O(V x A1)^*. (11)
```

The Keller open `(5)` lies in `V x A1`, because `pi` is etale on it.
Therefore `pr_1^*u|_U` is a unit on `A2`, hence a scalar. Since `U` is
dense in `X`, that equality would make `u` scalar in `k(C)`, contradicting
`(11)`. This proves

```text
no irreducible constant-q(P) cubic admits a same-field Keller A2 chart. (12)
```

The argument sees the maximal normalization, so conductor modifications of
the original monogenic order cannot evade it. It also avoids assuming that
the Keller open is a product open.

## 4. Repeated-adjusted-factor control: `q=P^2+P`

The first hostile from THM-3961 has

```text
f=T^3-3PT-P^2-P,
K=h^2(h-1)^2.                                            (13)
```

Its curve is irreducible: as a quadratic in `P`, `-f` has discriminant
`(T+1)^2(4T+1)`, which is nonsquare. Its full normalization is `k[v]` under

```text
P=v^3,                         T=v(v+1),
v=(P+T)/(T+1) in k(C).                                  (14)
```

Indeed `(14)` annihilates `f`; `v` is integral by `v^3-P=0`, the inverse
formula gives equality of function fields, and `k[v]` is normal.

The Jacobian equations first give `P=T^2` and then
`(T+1)(2T+1)=0`; the second candidate is not on `f=0`. Thus the unique
singular point `(P,T)=(1,-1)` has two normalization addresses

```text
v=omega, omega^2,                  omega^2+omega+1=0.     (15)
```

With `x=T+1,y=P-1`, the local equation is
`x^3-3x^2-3xy-y^2`; its quadratic tangent cone has discriminant `-3`, so
the two addresses are two distinct nodal branches. Thus the repeated hidden
factor is a conductor gluing, not genuine
ramification after normalization. The normalized target map is
`P=v^3`; its projective ramification points are `v=0,infinity`, both of
index three. The affine etale locus is `Gm_v`, with unit `v`. Hence the
universal proof reduces here to the immediate forbidden-`v` unit. Moreover

```text
Cl(k[v,t])=0,                       k[v,t]^*=k^*,          (16)
```

The ramification cylinder `(v=0)` is principal and primitive; deleting it
makes `v` a unit. Thus there is no hidden class-group boundary choice.

## 5. Double-`P` control: `q=cP^2`

For every `c in k^*`, the second hostile has

```text
f=T^3-3PT-cP^2,
K=h^3(ch-2).                                             (17)
```

The quadratic-in-`P` discriminant is `T^2(9+4cT)`, so the curve is again
irreducible. Its full normalization is `k[v]` under

```text
P=v^2(3+cv),                    T=v(3+cv),
v=P/T,                          cv^2+3v-T=0.             (18)
```

The last monic equation after division by the unit `c` proves integrality;
the rational inverse proves birationality. After `P=T^2`, the remaining
Jacobian row is `T(2cT+3)=0`, and the nonzero candidate is not on `f=0`.
Thus the singular origin is unique and has two conductor addresses

```text
v=0, -3/c.                                               (19)
```

Its tangent cone is `-P(3T+cP)`, so these are again two distinct nodal
branches.

The normalized polynomial map `P(v)` has

```text
dP/dv=3v(2+cv).                                          (20)
```

Its finite ramification points are `v=0,-2/c`, both of index two, and
`v=infinity` has index three. The affine etale curve is

```text
P1\{0,-2/c,infinity},
O(V)^*/k^* = Z*class(v) direct-sum Z*class(v+2/c).       (21)
```

Again the normalization cylinder is factorial with scalar global units,
while the principal ramification cylinders `(v=0)` and `(cv+2=0)` have
trivial class and deleting them creates the two forbidden units predicted
by `(11)`.

The constant `c=0` is intentionally excluded: then
`f=T(T^2-3P)` is reducible and fails the theorem's domain gate.

## 6. Scope and next escape

The theorem closes every `q=q(P)` specialization, not merely the two
controls `(13)` and `(17)`. It includes arbitrary degree, arbitrary hidden
factorization, arbitrary conductor gluing, and both normal and nonnormal
monogenic orders, provided the cubic function field is integral.

It does **not** prove JC(2). The mechanism depends critically on the product
decomposition `(3)`. A surviving conductor-debt design must make `q` vary
genuinely with `t`; zeros and poles of its moving coefficients can then add
vertical normalization components, classes, and boundary incidences. Those
vertical seams are the next exact target.

## Reproduction

```bash
python3 04-computation/jc2_constant_q_product_curve_ramification_unit_thm3962.py
python3 -O 04-computation/jc2_constant_q_product_curve_ramification_unit_thm3962.py
sha256sum 04-computation/jc2_constant_q_product_curve_ramification_unit_thm3962.py \
  05-knowledge/results/jc2_constant_q_product_curve_ramification_unit_thm3962.out
python3 agents/check_docs.py
```
