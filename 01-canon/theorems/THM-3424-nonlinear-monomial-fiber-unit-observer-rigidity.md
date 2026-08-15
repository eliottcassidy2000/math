---
id: THM-3424
title: "Nonlinear monomial-fiber unit-observer rigidity"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  Let K have characteristic
  zero, d>=2, a!=0, and
  P=ax+b+g(x)z^d with nonconstant g.  The generic Hamiltonian unit class
  vanishes exactly when g=c(x-alpha)^(1+qd) for some q>=1.  Its exact
  integral K[P]-annihilator is then ((P-(a alpha+b))^q), and is zero for
  every other nonconstant g.  No full multiroot module decomposition,
  polynomial mate, new Keller case, or JC(2) conclusion is claimed.
source: root-2608-jc-unit-observer-rigidity-2026-08-15
audit: independent weight-sector/valuation/two-infinity proof reconstruction; THM-3422 annihilator transport and boundary audit; normal/optimized/stored-output replay; hash and documentation audit clean
depends_on:
  - THM-3418-one-monomial-nonlinear-fiber-keller-classification
  - THM-3419-generic-kummer-response-regular-sector-rank
  - THM-3422-one-root-nonlinear-integral-hamiltonian-response
related:
  - THM-3348-linear-z-generic-puncture-response-and-one-root-valuation
  - THM-3354-inequivalent-h1-carriers-and-typed-obstruction-cospan
script: 04-computation/jc_nonlinear_multiroot_unit_rigidity_thm3424.py
output: 05-knowledge/results/jc_nonlinear_multiroot_unit_rigidity_thm3424.out
script_sha256: 22c29c3e0086dbe1a592d39027b41e1626484d032e3f7ddc9499629d4f294e38
output_sha256: 2d14bebf31cc736d158a74ce135d36f57c42663668630799e5d89d163f545921
hash_basis: LF-normalized bytes
---

# THM-3424 -- nonlinear monomial-fiber unit-observer rigidity

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**  The generic rigidity
proof below is self-contained.  The exact positive integral annihilator
imports the proved and independently audited one-root module calculation in
THM-3422.

## 1. Statement and inheritance

Let `K` be a field of characteristic zero, let `d>=2`, and put

```text
P=ax+b+g(x)z^d,                   a in K*,
D=D_P=Jac(P,-),                   C_P=K[x,z]/D(K[x,z]),
theta=[1] in C_P.                                          (1)
```

Assume `g` is nonconstant.  Let

```text
theta_gen=theta tensor_(K[P]) K(P).                        (2)
```

Then

```text
theta_gen=0
iff g(x)=c(x-alpha)^(1+qd)
    for c in K*, alpha in K, and an integer q>=1.          (3)
```

Combining `(3)` with the proved one-root integral classification in
[THM-3422](THM-3422-one-root-nonlinear-integral-hamiltonian-response.md)
gives the exact annihilator

```text
Ann_(K[P])(theta)=
  ((P-(a alpha+b))^q)       in the resonant one-root case (3),
  (0)                       for every other nonconstant g. (4)
```

The class in the first line of `(4)` is nonzero vertical torsion, not a
polynomial mate.  THM-3418 already proves that `[1]` is nonzero integrally
for every nonconstant `g` in this family.  The present theorem identifies
exactly when localization kills it and, with THM-3422, its complete
annihilator.

THM-3419 gives generic rank `N=deg(rad(g))` in every fiber sector.  Equal
rank does not locate the unit vector inside that packet.  The missing
coordinate is the weight-one generic-fiber equation below.

| item | exact content |
|---|---|
| source | the integral unit class `theta=[1]` |
| target | its generic Kummer class and its `K[P]` annihilator |
| map | localization and projection to the weight-one input sector |
| preserved | a nonzero generic class forces zero integral annihilator |
| destroyed | resonant one-root vertical torsion |
| required sidecars | every root valuation, `x=infinity` degree, and `t=infinity` order |
| cheapest hostile | `d=2`, `g=x(x-1)^3` passes the leading equation but fails the full infinity gate |

This is an observer classification only.  No local-root direct sum or full
multiroot integral module is asserted.

## 2. Normalization and the weight-one equation

Replace `P` by `(P-b)/a` and `g` by `g/a`.  This changes the Hamiltonian
differential by a nonzero scalar and replaces the polynomial generator by an
affine one, so both generic vanishing and annihilator exponents are unchanged.
It is enough to prove `(3)` for

```text
P=x+g(x)z^d.                                               (5)
```

Write `t` for the abstract generic value of `P` and `F=K(t)`.  The generic
fiber from THM-3419 is

```text
A=F tensor_(K[P]) K[x,z]
 ~=F[x,g^(-1),z]/(z^d-(t-x)/g).                           (6)
```

It is free over `F[x,g^(-1)]` with basis `1,z,...,z^(d-1)`.  The Hamiltonian
differential lowers fiber weight by one.  Hence, if `[1]=0` in `A/D(A)`, the
weight-one component of any primitive is already a primitive.  It is uniquely
of the form

```text
Q=q(x,t)z,                    q in F[x,g^(-1)].            (7)
```

Because `D(t)=0`, differentiation at fixed `t` and the identity
`gz^d=t-x` give

```text
D(qz)
 =q(1+g'z^d)-dgz^d q'
 =q+(t-x)((g'/g)q-dq').                                  (8)
```

Therefore

```text
theta_gen=0
iff q+(t-x)((g'/g)q-dq')=1
    for some q in F[x,g^(-1)].                            (9)
```

Equation `(9)` is both necessary and sufficient; no de Rham dimension count
is being substituted for exactness of this particular class.

## 3. Root valuations remove every finite denominator

Extend the coefficient field faithfully so that

```text
g=c product_(i=1)^N (x-alpha_i)^(e_i),                    (10)
```

with distinct roots.  Generic vanishing and exactness commute with this field
extension.  Put `y=x-alpha_i`.  Then `t-x` is a unit in the local Laurent
field and

```text
g'/g=e_i/y+(regular).                                     (11)
```

If `q=c_0 y^(-p)+...` has a pole of order `p>0`, the most singular term in
the parentheses in `(9)` is

```text
(e_i+dp)c_0 y^(-p-1).                                    (12)
```

It is nonzero in characteristic zero and cannot be cancelled by the first
`q` in `(9)`.  Thus `q` is regular at every root.  A nonzero constant term
would create the same uncancelled simple pole, so write

```text
q=c_0 y^n+...,                    n>=1.                    (13)
```

The first possible term contributed by the parentheses is

```text
(t-alpha_i)(e_i-dn)c_0 y^(n-1).                           (14)
```

Only `n=1` can equal the constant right side of `(9)`.  When `e_i=d`, its
coefficient vanishes and all remaining terms have positive order, so that
case is impossible as well.  Consequently every root is a simple zero of
`q`; since the only allowed finite poles were roots of `g`,

```text
q in F[x],                    rad(g) divides q.            (15)
```

## 4. The first infinity lock

Set

```text
r=deg(g),                    m=deg_x(q).                  (16)
```

Equation `(15)` gives `m>=N>=1`.  If `c_q x^m` is the leading term of `q`,
then

```text
(g'/g)q-dq'=(r-dm)c_q x^(m-1)+lower terms.                (17)
```

In `(9)`, the coefficient of `x^m` is therefore

```text
(1-r+dm)c_q.                                               (18)
```

The term multiplied by `t` has smaller `x`-degree and cannot cancel it.
Hence

```text
r=dm+1.                                                    (19)
```

In particular `r=1 mod d`.

## 5. The second infinity lock forces one root

Define the `t`-independent operator

```text
B(h)=(g'/g)h-dh'.                                         (20)
```

On polynomials divisible by `rad(g)`, `B` has no nonzero kernel under
`(19)`.  Indeed, `B(h)=0` implies

```text
(h^d/g)'=0,
```

so `h^d=Cg`; comparison of degrees would give `d|r`, contrary to `(19)`.

View the coefficients of `q` as Laurent series at `t=infinity`, and write
the leading term as

```text
q=t^L q_0(x)+lower t-orders,             q_0!=0.           (21)
```

Every Laurent coefficient remains divisible by `rad(g)`.  Equation `(9)` is

```text
q+(t-x)B(q)=1.                                             (22)
```

Since `B(q_0)!=0`, a leading order `L>=0` would leave an uncancelled positive
power `t^(L+1)`, while `L<=-2` could not produce the constant right side.
Thus

```text
L=-1,                         B(q_0)=1.                   (23)
```

Let `n=deg(q_0)>=N`.  Multiplying the second equation in `(23)` by `g`
gives

```text
g'q_0-dgq_0'=g.                                            (24)
```

If `n>1`, the coefficient in degree `r+n-1` forces

```text
r=dn,                                                      (25)
```

again contradicting `(19)`.  Therefore `n=1`.  Since `rad(g)|q_0` and `g`
is nonconstant,

```text
N=1.                                                       (26)
```

The unique geometric root is Galois fixed and hence belongs to `K`.  Write
`g=c(x-alpha)^e`.  Now `(19)` and `m>=1` give

```text
e=1+qd,                       q=m>=1.                     (27)
```

This proves the necessary direction of `(3)`.

## 6. Converse and the integral exponent

Suppose `(27)` holds and put

```text
y=x-alpha,                    tau=t-alpha.                (28)
```

The following coefficient belongs to `K(t)[x]`:

```text
q_gen=tau^(-q) sum_(j=0)^(q-1)
      [binom(q-1,j)/(1+jd)] y^(q-j)(tau-y)^j.             (29)
```

Direct substitution in `(9)` gives one.  Equivalently, `(29)` is the
THM-3422 endpoint primitive divided by the invertible generic scalar
`tau^q`.  This proves the sufficient direction of `(3)` independently of
the integral module splitting.

If the generic class is nonzero, no nonzero polynomial in `P` can annihilate
the integral class, so the second line of `(4)` follows immediately.  In the
resonant one-root case, THM-3422's orbit calculation places `[1]` exactly `q`
arrows before the unique zero arrow and proves

```text
Ann_(K[P])(theta)=((P-(a alpha+b))^q).                    (30)
```

This last step is precisely the integral information imported from THM-3422;
the generic rigidity argument does not silently infer an integral exponent
from localization.

## 7. Hostiles, computation, and boundaries

The two infinity locks are separately necessary.  For

```text
d=2,                         g=x(x-1)^3,                  (31)
```

the leading equation `(24)` has the exact solution `q_0=x(x-1)`, but
`deg(g)=4` violates `(19)`.  Thus a local or leading-asymptotic solution does
not imply generic exactness.  Direct multiroot systems that pass the initial
degree/root-count screen and still fail include multiplicity profiles

```text
(d;e_1,e_2)=(2;3,4), (3;2,5),
(d;e_1,e_2,e_3)=(4;3,5,5).                               (32)
```

The exact companion checks `775` ordered profiles, solves `107` exact linear
systems over `Q(t)`, verifies `144` translated one-root primitives, and
performs `367692` local/infinity coefficient checks.  Its declared finite
universe is `2<=d<=6`, one to three distinct roots at the canonical positions
`-1,1,3`, and multiplicities one through five, supplemented by
`(31)--(32)`.  It is evidence for the proof, not a replacement for its
unbounded quantifiers.

The boundaries are sharp:

- constant nonzero `g` has `C_P=0` integrally, so `[1]=0`; it is outside the
  nonconstant statement;
- `g=0` is a separate non-Kummer zero-response boundary;
- at `d=1`, the conclusion matches THM-3348's classification: the generic
  unit vanishes exactly for a single repeated root.  The present weight-sector
  proof is stated only for `d>=2`;
- characteristic zero is load bearing in `(12)--(14)`;
- even in `(3)`, generic vanishing is not a polynomial mate.  MISTAKE-374
  requires retention of the nonzero vertical torsion in `(4)`.

No full multiroot `K[P]`-module decomposition, polynomial mate, new Keller
stratum, or conclusion about the remaining cases of `JC(2)` follows.  QED.
