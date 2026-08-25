---
id: THM-4124
title: "Planar Keller integral-degree-ratio all-vertex target shear"
status: >
  PROVED FROM CITED LANG NEWTON SIMILARITY + elementary leading-form and
  convexity arguments + VERIFIED-EXACT controls. If a nonautomorphic planar
  Keller pair has deg Q=r*deg P with integral r, one scalar target shear
  Q-cP^r lowers deg Q and cancels every nonzero vertex of the scaled Newton
  polygon, after which Lang similarity gives strict radial contraction. A
  degree-sum-minimal representative under triangular target automorphisms
  therefore has neither degree dividing the other, so both entries of its
  reduced degree pair exceed one. The tempting 2:3 expression Q^2-cP^3 is
  not a target automorphism and does not preserve the Keller condition.
  JC(2) remains OPEN.
source: planar-jacobian-squeeze / 2026-08-25
audit: >
  PASS. Two independent agents rederived the common leading scalar, the
  all-vertex coefficient synchronization, strict radial containment, and
  response-fibre minimality consequence. They checked the r=1 boundary,
  rejected any extension to nonintegral 2:3 resonance, and supplied both a
  triangular-automorphism boundary control and a scaled-polygon non-Keller
  hostile. Normal and optimized exact-certificate runs byte-match the frozen
  output.
depends_on: []
related:
  - THM-2230-planar-jacobian-response-fiber-and-exact-target-shear-quotient
  - THM-3544-planar-keller-target-pencil-total-degree-six-floor
  - THM-4122-planar-keller-asymptotic-width-and-resonant-shear-contraction
external:
  - "J. Lang, Newton polygons of Jacobian pairs, J. Pure Appl. Algebra 72 (1991), 39--51, DOI:10.1016/0022-4049(91)90128-O."
  - "Nguyen Van Chau, Non-proper value set and the Jacobian condition, Ann. Polon. Math. 84 (2004), Theorem 1, arXiv:math/0305088."
script: 04-computation/jc2_asymptotic_width_shear_controls_thm4122_4124.py
output: 05-knowledge/results/jc2_asymptotic_width_shear_controls_thm4122_4124.out
script_sha256: 7fd2f9a2e32b199489506c22bf0f34c849ab5b0f4750c20e78d40038bd491f5a
output_sha256: d1bf27ab87eb78b79667854c3c57181942e248d83145a729e3460e641c8e2ca4
hash_basis: raw LF bytes
---

# THM-4124 -- an integral degree ratio admits one global vertex shear

**PROVED FROM CITED LANG SIMILARITY + elementary arguments + VERIFIED-EXACT
controls; JC(2) OPEN.** Let

```text
F=(P,Q):C^2 -> C^2,                  Jac(P,Q) in C*,      (1)
m=deg P>1,                 n=deg Q=r m>1,      r in N.    (2)
```

Assume `F` is not a polynomial automorphism. For a polynomial `R`, write

```text
N_0(R)=conv({(0,0)} union supp(R)).                       (3)
```

There is a single `c in C*` such that for every nonzero vertex `v` of
`N_0(P)`,

```text
[x^(r v)]Q=c([x^v]P)^r.                                  (4)
```

For that scalar,

```text
Q_1=Q-cP^r,
deg Q_1<n,
N_0(Q_1)=(deg Q_1/m)N_0(P) subsetneq rN_0(P).            (5)
```

Thus one genuine triangular target automorphism cancels every outer vertex,
not merely the ordinary top homogeneous face.

## 1. The leading face chooses the scalar

The highest homogeneous part of `Jac(P,Q)` vanishes. The standard homogeneous
dependence lemma, together with `n=rm`, gives a unique `c in C*` for which

```text
Q_n=c(P_m)^r.                                             (6)
```

Consequently the triangular target change

```text
(P,Q) -> (P,Q-cP^r)                                      (7)
```

preserves the constant Jacobian and strictly lowers the second total degree.
The target shear is invertible, so the new pair remains nonautomorphic.
Degree zero for `Q_1` contradicts its nonzero Jacobian with `P`. If
`deg Q_1=1`, an affine source change makes `Q_1=y`; then
`Jac(P,y)=P_x` is constant, whence `P=ax+h(y)` and the pair is an
automorphism. Thus `deg Q_1>1`, including when `r=1`.

## 2. Lang similarity propagates one cancellation to all vertices

Lang's theorem applied first to `(P,Q)` and then to `(P,Q_1)` gives

```text
N_0(Q)=rN_0(P),
N_0(Q_1)=lambda N_0(P),          lambda=deg Q_1/m<r.     (8)
```

If `v!=0` is a vertex of `N_0(P)`, then `rv` does not belong to
`lambda N_0(P)`: otherwise `rv=lambda z` for some `z` in `N_0(P)`, and

```text
v=(lambda/r)z+(1-lambda/r)0
```

would be a nontrivial convex combination. Hence the coefficient of
`x^(rv)` in `Q_1` is zero.

Choose a linear functional exposing `v`. In the product `P^r`, a sum of `r`
support exponents can equal `rv` only when all of them equal `v`. Therefore

```text
[x^(rv)]P^r=([x^v]P)^r.                                  (9)
```

The vanishing coefficient of `Q_1` and `(9)` prove `(4)`. This coefficient
sidecar is stronger than polygon similarity by itself.

## 3. A target-orbit normal form removes divisibility

Fix `P` and minimize the degree of the second component in the exact response
fibre

```text
Q+C[P].                                                   (10)
```

If `m` divided `deg Q`, equations `(5)` and `(7)` would produce a member of
`(10)` of smaller degree. Thus a degree-minimal response cannot have
`m | deg Q`.

More symmetrically, let `O(F)` be the orbit generated by the two triangular
target automorphisms

```text
(P,Q)->(P,Q+h(P)),              (P,Q)->(P+h(Q),Q).       (11)
```

The attained degree sums form a nonempty subset of the positive integers, so
`O(F)` has a global degree-sum minimizer. Every orbit member remains a
nonautomorphic Keller pair and has both degrees greater than one. At a global
minimizer, either divisibility invokes `(5)` or its swapped version and
strictly lowers the sum, a contradiction. Thus neither displayed degree
divides the other. If

```text
(deg P,deg Q)=G(d,e),                   gcd(d,e)=1,       (12)
```

it follows that

```text
d>=2,                              e>=2.                 (13)
```

This is an orbit-normal-form statement. It does not say that an arbitrary
presentation of the same pair has nondividing degrees.

## 4. Asymptotic synchronization and the 2:3 firewall

In Nguyen's monic source chart, if `deg Q=r deg P` and the leading `y`
coefficients are `A,B`, his Theorem 1 parametrizes each nonproperness
component with leading target terms

```text
(A xi^M, B xi^(rM)).                                      (14)
```

The scalar in `(6)` is `c=B/A^r`, so `(7)` cancels the first resonant target
term on every component. It does not synchronize the intrinsic multiplier,
the next Puiseux coefficient, or the amount of degree drop. The separate
coefficient identity printed later in that paper is not used here.

For the reduced `2:3` pair, the analogous expression

```text
R=Q^2-cP^3                                               (15)
```

is only a polynomial observable. The chain rule gives

```text
Jac(P,R)=2Q Jac(P,Q),                                    (16)
```

which is nonconstant. Thus `(15)` is not a target automorphism and no degree
descent follows. In particular this theorem does not contract the live
`(72,108)` cell.

## 5. Controls and boundary

- For `P=x+y^2`, `Q=y+P^r`, one has `Jac(P,Q)=1`; the polygons scale, all
  vertex ratios in `(4)` equal one, and the shear leaves `Q_1=y`. This is a
  triangular-automorphism boundary control, not a nonautomorphic example.
- For `P=x+y^2`, `Q=y+x^2+2y^4`, one still has
  `N_0(Q)=2N_0(P)`, but the two outer vertex ratios are `1,2` and the
  Jacobian is `1-4xy+8y^3`. Polygon similarity without the Keller equation
  does not synchronize coefficients.

The theorem excludes integral reduced degree ratios only after target-orbit
minimization. It does not handle rational nonintegral resonance, prove a
uniform degree drop, make the nonproperness set empty, or prove `JC(2)`.
