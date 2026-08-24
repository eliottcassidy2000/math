---
id: THM-3956
title: "Split hidden cubic integrality and repeated-root trichotomy"
status: >
  RESERVED / PROVISIONAL PROOF CANDIDATE / NOT CANON. If the one-parameter
  hidden cubic E+C h^2-2h^3 with C,E in k[t] splits over k(t), then every
  root is already polynomial: its monic normalization makes each root
  integral over the normal ring k[t]. Three distinct roots are exactly
  THM-3953. If a nonzero root is repeated, the missing h-row forces the
  third root to be minus one half of it and the natural depressed cubic is
  reducible. The only irreducible repeated-root exception has hidden roots
  0,0,r. Its normal cubic surface is polynomially isomorphic to
  u v=(u-2r)^3; the divisor P=0 is three times its unique reduced prime and
  is total ramification. Any Keller affine-plane open must delete that
  prime, making P a forbidden nonconstant unit. An independent Nagata
  calculation gives scalar units and class group
  Z^(s-1) plus Z/(3 gcd(m_j)) when r is nonconstant. Thus the complete
  k(t)-split hidden-cubic lane has no same-field nontrivial planar Keller
  chart. This candidate is not usable until independent hostile audit and
  status promotion.
source: jc-degree6-one-place + jc-extra-debt-local / post-THM-3953 split-field closure, 2026-08-24
depends_on:
  - THM-3953-rationally-split-hidden-cubic-ramification-triangle-nonentry
  - THM-3922-affine-plane-open-boundary-basis-class-group-obstruction
related:
  - THM-3954-extra-common-debt-creates-A-singularities-and-nonunibranch-residual-boundary
script: 04-computation/jc2_split_hidden_cubic_integrality_repeated_trichotomy_thm3956.py
output: 05-knowledge/results/jc2_split_hidden_cubic_integrality_repeated_trichotomy_thm3956.out
script_sha256: 1665a199083d668a0036089d3a80a9e82fef16a810c5b37608f5ac0b86c6a03e
output_sha256: a7afeba9c3af5c78d26423688cdc3ed276123f63e77c3bea3b4dfa02addfd1a2
semantic_sha256: d91f3715e9c4e6044650c878342402ae9d1a2ebf7ce244f8b3e4be4e61ccf0a7
hash_basis: raw LF bytes
---

# THM-3956 -- the split hidden-cubic lane is closed

**RESERVED / PROVISIONAL PROOF CANDIDATE / NOT CANON.** Work over an
algebraically closed field `k` of characteristic zero. Let

```text
C,E in k[t],
G(h,t)=E+C h^2-2h^3,                                      (1)
F(T,P,t)=T^3-3PT-(E+CP),
X=Spec k[P,t,T]/(F),             pi:X -> A2_(P,t).         (2)
```

Assume that `G` splits completely over `k(t)`. This theorem removes the two
formal escape clauses left by THM-3953: rational denominators and duplicate
roots. Its conclusion is narrowly about the natural finite cubic `(2)` and a
same-function-field planar Keller chart. It is not a classification of an
arbitrary cubic extension or of JC(2).

## 1. Rational denominators cannot occur

Dividing `(1)` by its leading unit gives the monic polynomial

```text
h^3-(C/2)h^2-E/2 in k[t][h].                              (3)
```

Every root in `k(t)` is integral over `k[t]`. Since the UFD `k[t]` is
integrally closed in `k(t)`, every such root lies in `k[t]`. Therefore

```text
G split over k(t)  implies  all three hidden roots lie in k[t].   (4)
```

In particular, clearing rational denominators cannot create a new vertical
boundary packet while retaining the monic polynomial model `(1)`. If the
three polynomial roots are distinct, THM-3953 applies verbatim and excludes
a same-field Keller affine-plane atlas. It remains only to classify repeated
roots.

## 2. The repeated-root trichotomy

Let the roots be `x,x,y in k[t]`. The missing linear `h`-row is the Vieta
identity

```text
x^2+2xy=x(x+2y)=0.                                       (5)
```

Because `k[t]` is a domain, there are exactly two branches.

### 2.1. The nonzero repeated-root branch is reducible

If `x` is not the zero polynomial, `(5)` gives `y=-x/2`. Consequently

```text
C=3x,                      E=-x^3,                       (6)
F=(T+x)(T^2-xT+x^2-3P).                                 (7)
```

Thus `(2)` is reducible and cannot be the integral finite normalization of a
connected cubic function field. More precisely, the normalization of the
reduced surface is the disjoint union of two affine planes:

```text
T=-x(t)                                      gives k[P,t],
P=(T^2-xT+x^2)/3                            gives k[T,t].   (8)
```

On the quadratic component, the coordinate `w=2T-x` gives

```text
P=x^2/4+w^2/12.                                           (9)
```

The two components are glued along `T=-x,P=x^2`; that squared conductor
factor is not a third cubic sheet. The specialization `x=0` in `(7)` is the
triple-zero case

```text
F=T(T^2-3P),                                               (10)
```

and is reducible for the same reason.

### 2.2. The exceptional repeated-zero branch is integral and normal

The other branch of `(5)` is

```text
x=0,                  roots 0,0,r(t),                    (11)
```

where `r=y`. The case `r=0` is already `(10)`, so assume `r!=0`. Then

```text
C=2r,                    E=0,
F=T^3-P(3T+2r).                                           (12)
```

The polynomial `(12)` is irreducible. Viewed as a primitive polynomial
linear in `P`, its two coefficients are `T^3` and `3T+2r`; their gcd in
`k[t,T]` is one because `r` is a nonzero polynomial independent of `T`.
Gauss's lemma proves that `X` is integral.

There is an exact polynomial change of variables

```text
u=3T+2r(t),                 v=27P,                        (13)
T=(u-2r(t))/3,              P=v/27,                       (14)
```

under which `(12)` becomes

```text
X isomorphic to Spec A_r,
A_r=k[t,u,v]/(u v-(u-2r(t))^3).                           (15)
```

This integral surface is normal. For

```text
H=u v-(u-2r)^3,                                           (16)
```

the singular equations are

```text
H_v=u,
H_u=v-3(u-2r)^2,
H_t=6r'(u-2r)^2.                                         (17)
```

Together with `H=0`, they force

```text
u=v=r(t)=0.                                               (18)
```

Hence the singular locus is finite. The hypersurface domain is `S2`, and a
finite singular locus on a surface has codimension two, so it is `R1`.
Serre's criterion proves that `A_r` is normal. Multiple roots of `r` may
worsen the isolated singularity but do not create a height-one defect.

## 3. Total ramification makes `P` a forbidden unit

On the repeated-zero surface `(12)`, the scheme-theoretic base divisor is

```text
X x_(A2) V(P)=Spec k[t,T]/(T^3).                          (19)
```

It therefore has exactly one reduced height-one prime

```text
D=(P,T),                      div_X(P)=3D.                (20)
```

At the generic point of `D`, `r` is a unit and

```text
P=T^3/(3T+2r),                                            (21)
```

so `(20)` is genuine total ramification of index three for `pi`. The other
ramification prime is obtained from

```text
F|_(P=T^2)=-2T^2(T+r),                                   (22)
```

but it is not needed for the obstruction.

Suppose a same-function-field planar Keller map existed. Since `X` is the
normal finite model of its function field over `k[P,t]`, normalization-form
Zariski Main would identify the source affine plane with an open

```text
j:U isomorphic to A2 -> X                               (23)
```

on which `pi` is etale. The generic point of the ramified prime `D` cannot
lie in `U`, so `D` is a boundary divisor. But `(20)` says that `P` has no
zero on `U`; hence

```text
P|_U in Gamma(U,O_U)^*.                                  (24)
```

The only units on `A2` are scalars, whereas dominance makes the target
coordinate `P` nonconstant. This contradiction excludes the Keller chart.
Notice that the argument allows arbitrary additional deleted divisors: once
`D` is deleted, `(20)` already makes `P` a unit.

## 4. Independent unit and class-group anatomy

The model `(15)` also exposes the exact vertical debt. Localizing at `u`
gives

```text
A_r[u^(-1)]=k[t,u,u^(-1)].                               (25)
```

If `r` is a nonzero scalar, `u` is already a unit on all of `X`: expanding
`(15)` gives

```text
8r^3=u(u^2-6ru+12r^2-v).                                 (26)
```

Thus `A_r=k[t,u,u^(-1)]`, its class group is zero, and the nonconstant unit
`u` independently forbids a dense affine-plane open.

Now suppose

```text
r=kappa product_(j=1)^s (t-alpha_j)^(m_j),
kappa in k*,             alpha_i distinct,               (27)
```

with `s>=1` and `m_j>=1`. The height-one primes above `u=0` are exactly

```text
Q_j=(u,t-alpha_j).                                        (28)
```

At the generic point of `Q_j`, `v` is a unit. If `nu_j` is the normalized
valuation, `(15)` first forces

```text
nu_j(u)>m_j nu_j(t-alpha_j),
nu_j(u)=3m_j nu_j(t-alpha_j).                             (29)
```

Since the maximal ideal is generated by `u,t-alpha_j`, the latter element
is a uniformizer. Therefore

```text
nu_j(t-alpha_j)=1,             nu_j(u)=3m_j.              (30)
```

Nagata localization applied to `(25)` now gives

```text
A_r^*=k*,
Cl(A_r)=Z^s / <(3m_1,...,3m_s)>
       isomorphic to Z^(s-1) plus Z/(3 gcd(m_1,...,m_s)). (31)
```

The total-ramification prime `D` from `(20)` has the complementary divisor
ledger

```text
div(u-2r)=D+sum_j m_j Q_j,
div(v)=3D.                                                  (32)
```

Thus `[D]` is the class of `-sum_j m_jQ_j` in `(31)` and has exact order
three: its first and second multiples cannot lie in the one relation
generated by `(3m_j)`, while its third multiple is principal. The nonzero
`3`-torsion in `(31)` is therefore a second, independent THM-3922
obstruction to a dense affine-plane open. It also records exactly why
isolated zeros of `r` do not repair the total-ramification failure: they add
vertical primes with a common factor three in their valuation relation.

## 5. Exhaustion and scope

The complete split-field decision tree is therefore

```text
G splits over k(t)
  -> every root is polynomial;
  -> three distinct roots: THM-3953 boundary triangle or class/unit gate;
  -> repeated nonzero root: reducible natural cubic;
  -> repeated zero, third root nonzero: normal surface, forbidden-P-unit;
  -> triple zero: reducible natural cubic.                            (33)
```

There is no rational-denominator survivor and no duplicate-root survivor in
this natural one-parameter model. The proof does **not** say that every
hidden cubic must split, nor does it exclude changes of target or of finite
model that leave the same-function-field natural-cubic scope. Irreducible
hidden cubics, non-`A1` primary branches, and general planar Keller maps
remain outside. **QED (candidate pending independent audit).**
