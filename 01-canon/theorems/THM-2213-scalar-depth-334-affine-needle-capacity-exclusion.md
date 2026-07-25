---
id: THM-2213
title: "Scalar depth-(3,3,4) affine-needle capacity exclusion"
status: >
  PROVED + VERIFIED-EXACT. In the scalar five-unit/three-blocker branch,
  the actual blocker valuation profile (3,3,4) is empty. On the primitive
  13^5 guard-safe annulus, depth-three blockers are constant on the 156
  oriented residue stalks modulo 13^2. An exact integer carrier records
  the guard-safe dangerous-root capacity of every one of 171,366 unit sign
  classes on every stalk. All 3,081 unordered depth-three blocker pairs,
  including 78 diagonals, have positive conditional top-five deficit. The
  unique minimum margin is 24,022 at blocker labels (1,84). Direct
  full-torsion replay verifies the hostile residual and five capacities.
  This closes one depth-four profile only; it is not a proof of LRC(14).
source: klein-2026-07-24-scalar-depth334-affine-needle
depends_on:
  - THM-2192-scalar-five-plus-three-root-sheet-chord-invoice
  - THM-2198-scalar-five-plus-three-image-pump-and-first-depth-exclusion
related:
  - THM-2201-cyclic-root-fibre-hasse-jet-transition-carrier
  - THM-2203-fixed-dyadic-coordinate-section-and-covector-intersection
  - THM-2204-scalar-depth-223-thirteen-lift-capacity-law
  - THM-2207-scalar-depth-123-labelled-guard-hole-exclusion
  - THM-2216-residual-capacity-hinge-gram-law
  - THM-2218-labelled-guard-hole-fourier-and-signed-lift-energy
script: 04-computation/lrc14_scalar_depth334_cyclic_capacity_certificate_thm2213.py
output: 05-knowledge/results/lrc14_scalar_depth334_cyclic_capacity_certificate_thm2213.out
script_sha256: 21d8fcfd9f0ff983394183edb1e61b15372638e00a1ade3b4e594e4b4a2dae92
output_sha256: 661cd15fab3d72d6e4dd250151595e26ad8c61d7cd2ee6f9699c1269574db090
hash_basis: working-tree bytes (LF)
---

# THM-2213 -- scalar depth-`(3,3,4)` exclusion

Put

```text
D_a={t in R/Z:||at||<1/14},
C_H={t in R/Z:||Ht||>1/7}.                           (1)
```

In the scalar `5+3` branch of THM-2192 and THM-2198, suppose

```text
C_H subset union_(i=1)^5 D_(q_i)
             union union_(j=1)^3 D_(c_j)             (2)
```

almost everywhere.  The guard `H` and five coefficients `q_i` are
positive thirteen-units; the three actual blockers have positive
thirteen-valuation.  After relabelling,

```text
1<=lambda_1<=lambda_2<lambda_3,
lambda_j=nu_13(c_j).                                 (3)
```

This theorem excludes

```text
(lambda_1,lambda_2,lambda_3)=(3,3,4).                (4)
```

## 1. Primitive depth-five carrier

Assume (4) and put

```text
N=13^5=371293,        Q=13^4=28561,
R=13^2=169.                                           (5)
```

Multiplication by `H^(-1) modulo N` normalizes the guard to one and
preserves all three valuations.  Every primitive residue modulo `N` has a
unique root-fibre presentation

```text
z=r+kQ,
r in {1,...,Q-1}, 13 does not divide r,
k in F_13.                                          (6)
```

Retain precisely the guard-safe sheets

```text
U_N={z:13 does not divide z and 7||z||_N>N}.         (7)
```

There are

```text
26364 primitive image phases r,
|U_N|=244810,
18830 nine-sheet fibres and 7534 ten-sheet fibres.  (8)
```

A depth-four blocker is harmless on this layer.  Its normalized
coefficient is `13^4w`, so at a primitive `z`

```text
13^4wz/N=wz/13
```

is a nonzero thirteenth root and has norm at least `1/13>1/14`.
The script checks all six sign classes directly.

## 2. The affine-needle stalks

Write the image phases by their oriented residue

```text
s=r mod R,        s in (Z/RZ)^*.                     (9)
```

There are 156 such residues, and each stalk

```text
P_s={r:same residue s mod R}                         (10)
```

contains `Q/R=169` image phases.  A depth-three blocker has normalized
coefficient `13^3u`; its action on (6) is determined by `uz/R`.
Consequently it is constant on `P_s`, with active-stalk indicator

```text
A_u(s)=1_{14||us||_R<R}.                             (11)
```

Up to sign there are 78 possible depth-three unit parts `u`.

For a unit terminal label `q modulo N`, define its exact stalk capacity

```text
H(q,s)
 =#{z in U_N:
      z lies over an image phase in P_s,
      14||qz||_N<N}.                                (12)
```

The root-sheet law gives at most two terminal sheets per image phase, so

```text
0<=H(q,s)<=2*169=338.                                (13)
```

The 171,366 unit sign classes modulo `N` split into the thirteen raw lifts
of the 13,182 base sign classes modulo `Q`.  The script constructs (12)
family by family and verifies that the canonical signs form exactly the
full unit universe.  This is the finite affine-needle carrier: `s` is the
stalk, a unit label is the needle, and `H(q,s)` is their weighted
incidence.  No Euclidean Kakeya theorem is invoked.

## 3. Exact two-blocker capacities

For depth-three blocker labels `u,v`, put

```text
P_(u,v)={s:A_u(s)=A_v(s)=0},                         (14)

W_(u,v)=sum_(s in P_(u,v))
          #{guard-safe sheets over P_s},             (15)

C_(u,v)(q)=sum_(s in P_(u,v))H(q,s).                (16)
```

If the five unit masks in (2) covered the residual, the union bound would
force

```text
W_(u,v)<=Top_5(C_(u,v)),                             (17)
```

where the right side is the sum of the five largest unit-label
capacities.

For fast exact evaluation, write

```text
F(q)=sum_s H(q,s),
X_u(q)=sum_(s:A_u(s)=1)H(q,s).                      (18)
```

Then

```text
C_(u,v)(q)
 =F(q)-X_u(q)-X_v(q)
   +sum_(s:A_u(s)=A_v(s)=1)H(q,s).                  (19)
```

Equations (13) and (19) give the pointwise upper envelope

```text
C_(u,v)(q)<=U_(u,v)(q)

U_(u,v)(q)=min(
  F(q)-X_u(q),
  F(q)-X_v(q),
  F(q)-X_u(q)-X_v(q)
    +338|A_u intersection A_v|
).                                                   (20)
```

The first two terms are the residual capacities after just one blocker;
the third is inclusion--exclusion with the sharp stalkwise height bound.

## 4. Complete prefix certificate

For a diagonal pair, (16) is a singleton residual row and its top five are
selected directly.  For `u!=v`, take the 64 largest values of (20), sort
them by the upper envelope, and evaluate the exact expression (19) until
the next upper bound is strictly below the current fifth exact value.
Every unevaluated label then has smaller exact capacity, so the retained
five are the global exact top five.  This is branch-and-bound, not
sampling.

The prefix never expands.  Across all

```text
C(78+1,2)=3081
```

unordered pairs with repetition, only 35,974 exact candidate evaluations
are needed.  Symmetry in `(u,v)` covers all `78^2` ordered pairs.

Every pair satisfies

```text
W_(u,v)-Top_5(C_(u,v))>0.                            (21)
```

The unique minimum row is

```text
(u,v)=(1,84),
W_(1,84)=188290,

five largest (capacity,label):
(33928,2198),
(33918,2196),
(32518,1098),
(32496,1099),
(31408,6),

minimum margin=24022.                                (22)
```

The complete pair-table digest is

```text
97ae53fb2a59bf4bd8a34e18dfb919063ef502eed938fda9fb30c0c901fdea68.
```

The branch-trace digest is

```text
3f92815e67a77af35d97d44cbf50cd5622b444b0d1e13917e8fd44bc549aac3f.
```

## 5. Independent controls

The certificate includes:

1. definition-level agreement between the grouped indicator (11) and all
   78 depth-three blockers on every root sheet;
2. direct full-torsion reconstruction of the hostile residual in (22) and
   all five reported capacities, without the grouped representation;
3. direct emptiness of all six depth-four blocker sign classes on `U_N`;
4. the exact thirteen-lift family-sum identity on the hostile residual for
   all 13,182 base families;
5. explicit fixed-width accumulator bounds and load-bearing `require`
   checks that remain active under `python -O`.

Normal and optimized runs are byte-identical.  The computation uses exact
integer arithmetic, with no floating point.

## 6. Exclusion and scope

For the actual pair of depth-three blockers, (21) contradicts the
necessary cover inequality (17).  Hence a primitive guard-safe residue
avoids the two depth-three blockers and all five unit masks.  It also
avoids the depth-four blocker by Section 1.  Since powers of thirteen are
coprime to seven and fourteen, no equality endpoint occurs; the residue
thickens to an uncovered open interval, contradicting (2).  Therefore
the profile `(3,3,4)` is empty.

This theorem settles no other depth-four profile.  In particular it does
not infer `(2,3,4)` from monotonicity, and it gives no uniform statement
for deepest valuation at least five.  LRC(14) remains open.  QED.
