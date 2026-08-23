---
id: THM-3876
title: "Higher monomial marked-root profiles fail branch-ring descent"
status: >
  RESERVED / PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT; awaiting
  independent hostile audit.  For the entire monomial normalization tower
  (A,C)=(t^m,6t(1+t^(m+1))), the cubic cusp identity uniquely forces the
  marked value B=2t^2(3+4t^(m+1)).  This value descends to the plane branch
  ring for m=1 and m=2, but for every m>=3 an explicit primitive-root pair
  of distinct normalization addresses has the same (A,C) and different B.
  Hence no polynomial carrier profile b(A,C) can contain that marked branch.
source: jc_sparse_direct_search / post-THM-3873 monomial normalization tower, 2026-08-23
audit: >
  SELF-AUDITED proof candidate.  The assertion-free exact companion checks
  the universal cusp square and forced marked value, the primitive-root
  collision and unequal-value formula, the m=1 forward-graph and m=2
  triangular-parabola boundaries, and exact cyclotomic hostile controls for
  m=3,...,20.  Normal and optimized runs byte-match the frozen 82-gate
  transcript.  The all-m argument itself is symbolic and does not depend on
  the finite replay range.  Independent hostile audit remains required.
related:
  - THM-3866-all-polynomial-graph-branches-force-projective-companion
  - THM-3873-first-nongraph-triangular-parabola-companion
script: 04-computation/jc2_monomial_marked_root_nondescent_thm3876.py
output: 05-knowledge/results/jc2_monomial_marked_root_nondescent_thm3876.out
script_sha256: 3d4fc0cd92152e84dd8d95127350c2372848a6f997bf4bc7027406d98d5a3e08
output_sha256: aff11d64e97cbd2be4b9c35742c01f483ff5346311ffe3ef7c5f8eae29952124
semantic_sha256: b42121616115392fe501a51994d7e489a1431d62acc20ac5373dbc43b42962df
hash_basis: raw LF bytes
---

# THM-3876 -- the monomial marked-root tower stops at the first non-graph branch

**RESERVED / PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT; awaiting
independent hostile audit.**  Work over an algebraically closed field `k` of
characteristic zero.  Put

```text
Delta_b=-27A^2b^2+8AC^3-54ACb+9C^2-54b.                       (1)
```

For an integer `m>=1`, let `Gamma_m` be the irreducible plane curve with
normalization

```text
nu_m: A1_t -> A2_(A,C),
A=t^m,                 C=6t(1+t^(m+1)).                        (2)
```

If `b in k[A,C]` and `Delta_b` vanishes on `Gamma_m`, then its restriction
to the normalization is forced to be

```text
b(A(t),C(t))=B_m(t):=2t^2(3+4t^(m+1)).                        (3)
```

The normalization map `(2)` is finite and birational for every `m`.  The
forced value `(3)` descends to the branch ring exactly at the first two
levels of this tower:

```text
m=1: B_1=2A^2(3+4A^2),
     C=6A(1+A^2);                                      forward graph,

m=2: T=C/6-A^2=t,
     B_2=6A+8A^2T;                              triangular parabola,

m>=3: B_m is not an element of k[t^m,t+t^(m+2)].                (4)
```

Consequently, for every `m>=3` there is **no** polynomial `b(A,C)` whose
discriminant packet `V(Delta_b)` contains `Gamma_m`.  The obstruction occurs
before any choice of transverse quotient: the uniquely forced function on
the normalization does not descend to the singular plane branch.

This theorem is confined to the monomial family `(2)`.  It does not close
arbitrary singular plane branches having affine-line normalization, nor
non-graph branches outside this tower.

## 1. The branch ring and its normalization

The coordinate ring of the image of `(2)` is

```text
R_m=k[t^m,t+t^(m+2)] subset k[t].                              (5)
```

The harmless factor `6` in `C` has been removed in `(5)`.  The element `t`
is integral over `R_m`, since it is a root of the monic polynomial
`X^m-A in R_m[X]`, where `A=t^m`.  Thus `k[t]` is finite over `R_m`.

It remains to check birationality rather than assume that an affine-line
parametrization is an embedding.  Away from `t=0`, two points `t` and `t'`
with equal `A` have

```text
t'=zeta t,                     zeta^m=1.                       (6)
```

Equality of their `C`-coordinates is equivalent to

```text
(zeta-1)+(zeta^2-1)t^(m+1)=0.                                 (7)
```

For `zeta=1` this is the diagonal.  For `zeta=-1`, equation `(7)` is
impossible.  Every other `m`th root of unity gives at most one value of
`t^(m+1)`.  Hence all off-diagonal coincidences form a finite set and the
generic fibre of `nu_m` has one point.  The finite map is therefore
birational, so

```text
Frac(R_m)=k(t).                                                 (8)
```

Since `k[t]` is integrally closed and integral over `R_m`, it is exactly the
integral closure of `R_m` in `(8)`: any element of `k(t)` integral over
`R_m` is also integral over `k[t]`, hence lies in `k[t]`.

## 2. The cusp identity forces one marked value

The universal cubic identity behind `(1)` is

```text
A^2 Delta_b=27(P^3-u^2),
P=1+(2/3)AC,                 u=1+AC+A^2b.                      (9)
```

Set

```text
U=t^(m+1).                                                     (10)
```

Along `(2)`, one has

```text
P=1+4U(1+U)=(1+2U)^2.                                        (11)
```

If `Delta_b` vanishes on `Gamma_m`, then in the domain `k[t]`

```text
u^2=(1+2U)^6,
```

so `u=(1+2U)^3` or `u=-(1+2U)^3`.  At `t=0`, however, `A=C=0` and therefore `u(0)=1`.
Characteristic zero excludes the negative sign.  Solving for the
restriction `B(t)=b(A(t),C(t))` gives

```text
t^(2m) B(t)
 =(1+2U)^3-1-6U(1+U)
 =2U^2(3+4U),                                                  (12)

B(t)=2t^2(3+4t^(m+1))=B_m(t).                                 (13)
```

Thus `(3)` is necessary for every polynomial carrier profile; there is no
choice of sign or transverse term left on the marked component.

## 3. The two descent boundaries

At `m=1`, the normalization coordinate is already `A=t`, and `(2)--(3)`
become

```text
C=6A(1+A^2),                 B_1=2A^2(3+4A^2).                 (14)
```

This is the forward-graph boundary covered by THM-3866.

At `m=2`, the triangular polynomial coordinate

```text
T=C/6-A^2                                                         (15)
```

restricts to `t+t^4-t^4=t`.  Consequently

```text
R_2=k[T],                       A=T^2,
B_2=6A+8A^2T.                                                   (16)
```

This is exactly the first fixed-coordinate non-graph triangular parabola
of THM-3873.  In particular, the distinction in `(4)` is a branch-ring
distinction: normalization by `A1` alone does not decide descent.

## 4. A uniform double address for every m>=3

Assume `m>=3` and choose a primitive `m`th root of unity `zeta`.  Then

```text
zeta != 1,                 zeta != -1,                 zeta+1 !=0.   (17)
```

The last two assertions follow because a primitive root has order `m>=3`.
Choose `t in k^*` satisfying

```text
t^(m+1)=-1/(zeta+1).                                          (18)
```

Such a nonzero `t` exists because `k` is algebraically closed.  The two
normalization addresses `t` and `zeta t` are distinct.  They have equal
`A`-coordinates because `zeta^m=1`.  Their `C/6` difference, divided by
`t`, is

```text
(zeta-1)+(zeta^2-1)t^(m+1)=0                                  (19)
```

by `(18)`.  Thus

```text
nu_m(t)=nu_m(zeta t).                                          (20)
```

The forced marked values do not agree.  Using `zeta^(m+1)=zeta`, one finds

```text
B_m(zeta t)-B_m(t)
 =2t^2[3(zeta^2-1)+4(zeta^3-1)t^(m+1)]
 =-2t^2 (zeta-1)^3/(zeta+1) !=0.                              (21)
```

Every factor in the final expression is nonzero by `(17)--(18)`.

## 5. Non-descent and the polynomial-carrier contradiction

Every polynomial `b(A,C)` restricts along `(2)` to an element of the branch
ring `R_m`.  Equivalently, its pullback must take the same value at any two
normalization addresses with the same `(A,C)`.  Equations `(20)--(21)` show
that `B_m` fails this necessary condition for every `m>=3`; hence

```text
B_m notin R_m.                                                  (22)
```

But Section 2 proves that a polynomial carrier with `Delta_b|_(Gamma_m)=0`
would have pullback exactly `B_m`.  This contradiction proves the claimed
all-degree nonentry.  Notice that no squarefreeness, reducedness, or
transverse-factor assumption is used: the obstruction is already present
in the restriction of `b` to the reduced image branch.

## 6. Exact replay and scope

The exact companion

```text
python3 04-computation/jc2_monomial_marked_root_nondescent_thm3876.py
python3 -O 04-computation/jc2_monomial_marked_root_nondescent_thm3876.py
```

checks the universal identities `(11)--(13)` and `(19)--(21)`, the two
descent boundaries `(14)--(16)`, and exact cyclotomic exclusions for every
`3<=m<=20`.  Those bounded cyclotomic rows are hostile controls only.  The
proof for arbitrary `m>=3` is the order argument `(17)` followed by the
single universal calculation `(18)--(21)`.

The successful mechanism is **normalization-address separation**: a cusp
identity may force a perfectly polynomial function of the normalization
parameter while the plane branch ring forgets enough address information
to prevent that function from descending.  THM-3873 is the last monomial
level where a polynomial triangular sidecar restores the lost parameter;
from `m=3` onward the explicit double address kills the carrier before any
projective-companion or Jacobian-mate analysis begins.
