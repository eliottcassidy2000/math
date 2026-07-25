---
id: THM-2215
title: "Scalar depth-(2,3,4) affine-needle capacity exclusion"
status: >
  PROVED + VERIFIED-EXACT. In the scalar five-unit/three-blocker branch,
  the actual blocker valuation profile (2,3,4) is empty. On the primitive
  13^5 guard-safe annulus, the exact unit-capacity carrier is grouped into
  2,028 oriented residue stalks modulo 13^3. Every one of the 1,014
  depth-two blocker classes and 78 depth-three blocker classes agrees
  definitionally with its grouped mask. All 79,092 typed pairs have
  positive conditional top-five deficit. The unique minimum margin is
  27,440 at blocker labels (844,1). Direct full-torsion replay verifies the
  hostile residual and five capacities. This closes one depth-four profile
  only; it is not a proof of LRC(14).
source: klein-2026-07-24-scalar-depth234-affine-needle
depends_on:
  - THM-2192-scalar-five-plus-three-root-sheet-chord-invoice
  - THM-2198-scalar-five-plus-three-image-pump-and-first-depth-exclusion
related:
  - THM-2201-cyclic-root-fibre-hasse-jet-transition-carrier
  - THM-2203-fixed-dyadic-coordinate-section-and-covector-intersection
  - THM-2204-scalar-depth-223-thirteen-lift-capacity-law
  - THM-2207-scalar-depth-123-labelled-guard-hole-exclusion
  - THM-2213-scalar-depth-334-affine-needle-capacity-exclusion
  - THM-2216-residual-capacity-hinge-gram-law
  - THM-2218-labelled-guard-hole-fourier-and-signed-lift-energy
script: 04-computation/lrc14_scalar_depth234_cyclic_capacity_certificate_thm2215.py
output: 05-knowledge/results/lrc14_scalar_depth234_cyclic_capacity_certificate_thm2215.out
script_sha256: a3636453f4819be375e3b673ec8b801bacf40235eaf9cc657a1a75d4efc529ab
output_sha256: 809d235eb59e8cd888780224d310f526ac34c9835713ea3eba4734e00f909cb5
hash_basis: working-tree bytes (LF)
---

# THM-2215 -- scalar depth-`(2,3,4)` exclusion

Use

```text
D_a={t in R/Z:||at||<1/14},
C_H={t in R/Z:||Ht||>1/7}.                           (1)
```

In the scalar `5+3` branch of THM-2192 and THM-2198, suppose

```text
C_H subset union_(i=1)^5 D_(q_i)
             union union_(j=1)^3 D_(c_j)             (2)
```

almost everywhere.  The guard and five `q_i` are positive
thirteen-units, while

```text
1<=lambda_1<=lambda_2<lambda_3,
lambda_j=nu_13(c_j).                                 (3)
```

This theorem excludes

```text
(lambda_1,lambda_2,lambda_3)=(2,3,4).                (4)
```

## 1. Primitive carrier and stalk scale

Assume (4) and normalize the guard to one modulo

```text
N=13^5=371293,       Q=13^4=28561.                  (5)
```

Every primitive residue has a unique presentation

```text
z=r+kQ,
r in {1,...,Q-1}, 13 does not divide r,
k in F_13.                                          (6)
```

The guard-safe carrier

```text
U_N={z:13 does not divide z and 7||z||_N>N}          (7)
```

has 244,810 sheets over 26,364 primitive image phases: 18,830
nine-sheet fibres and 7,534 ten-sheet fibres.

The depth-four blocker is harmless throughout `U_N`, since its normalized
value is a nonzero thirteenth root.  The script checks all six depth-four
sign classes directly.

For the mixed shallower pair use the oriented residue

```text
s=r mod 13^3,
s in (Z/13^3Z)^*.                                   (8)
```

There are 2,028 such stalks, each containing thirteen image phases.  A
depth-two blocker `13^2u` is constant on a stalk because its value depends
only on `uz/13^3`.  A depth-three blocker `13^3v` is also constant there,
with its mask pulled back from `s modulo 13^2`.

Up to sign there are

```text
1014 depth-two labels u modulo 13^3,
  78 depth-three labels v modulo 13^2.              (9)
```

## 2. Exact unit incidence matrix

For a unit terminal label `q modulo N`, define

```text
H(q,s)
 =#{z in U_N:
      z lies over an image phase in stalk s,
      14||qz||_N<N}.                                (10)
```

Each stalk has thirteen image phases and every terminal mask has at most
two active sheets per phase, so

```text
0<=H(q,s)<=26.                                      (11)
```

The 171,366 unit sign classes modulo `N` split into thirteen raw lifts of
each of the 13,182 base sign classes modulo `Q`.  The script constructs
all rows of (10), verifies that their canonical signs partition the full
unit universe, and checks every grouped blocker mask against its original
definition on every root:

```text
1014/1014 depth-two classes PASS,
  78/78   depth-three classes PASS.                 (12)
```

This is again an affine-needle incidence model, now at the finer stalk
scale `13^3`.  It invokes no Euclidean Kakeya theorem.

## 3. Mixed residual capacities

Let `A_u` and `B_v` be the active stalk sets of the depth-two and
depth-three blockers.  Put

```text
P_(u,v)=(Z/13^3Z)^* \(A_u union B_v),                (13)

W_(u,v)=sum_(s in P_(u,v))
          #{guard-safe sheets over s},              (14)

C_(u,v)(q)=sum_(s in P_(u,v))H(q,s).                (15)
```

If five unit masks covered the residual, then

```text
W_(u,v)<=Top_5(C_(u,v)).                            (16)
```

Write

```text
F(q)=sum_s H(q,s),
X_u(q)=sum_(s in A_u)H(q,s),
Y_v(q)=sum_(s in B_v)H(q,s).                        (17)
```

Exact inclusion--exclusion gives

```text
C_(u,v)(q)
 =F(q)-X_u(q)-Y_v(q)
  +sum_(s in A_u intersection B_v)H(q,s).           (18)
```

By (11),

```text
C_(u,v)(q)<=U_(u,v)(q)

U_(u,v)(q)=min(
  F(q)-X_u(q),
  F(q)-Y_v(q),
  F(q)-X_u(q)-Y_v(q)
    +26|A_u intersection B_v|
).                                                   (19)
```

Every removal total is a sub-sum of a full unit capacity below `2^16`.
Thus the two large `uint16` matrix products in the artifact are literal
integer accumulations without wraparound.  Four hostile columns of each
type are independently recomputed as `uint32` sums.

## 4. Complete typed-pair certificate

For each of the

```text
1014*78=79092                                        (20)
```

typed pairs `(u,v)`, select the 64 largest values of (19), evaluate the
exact expression (18) in decreasing upper-bound order, and stop only when
the next bound is strictly below the current fifth exact value.  The
unevaluated labels then cannot enter the global exact top five.

No prefix expands.  The full scan uses 878,784 exact candidate
evaluations, and every pair satisfies

```text
W_(u,v)-Top_5(C_(u,v))>0.                            (21)
```

The unique minimum row is

```text
(u,v)=(844,1),
W_(844,1)=175248,

five largest (capacity,label):
(29664,142637),
(29646,2198),
(29646,142635),
(29628,2196),
(29224,6),

minimum margin=27440.                                (22)
```

The complete pair-table digest is

```text
f141ef307efc773fb90e6a20c534f4a2d754b31bd1f63658591d34551161af3f.
```

The branch-trace digest is

```text
60c2bdcfbe8b91e2b2bdead424d2b3bcde29ff69af62aaae3f58356d1c74580e.
```

## 5. Independent controls

The artifact also:

1. directly reconstructs the hostile residual in (22) on the full
   primitive torsion layer and recomputes all five reported capacities;
2. verifies emptiness of all six depth-four blocker classes on `U_N`;
3. checks the exact thirteen-lift family-sum identity on the hostile
   residual for all 13,182 base families;
4. freezes the pair and branch traces, candidate histogram, universe
   sizes, endpoint gates, and fixed-width bounds using `require`, not
   `assert`.

All decisions use integers; there is no floating point.  Normal and
optimized runs are byte-identical.

## 6. Exclusion and scope

For the actual mixed pair, (21) contradicts the necessary cover inequality
(16).  Therefore a primitive guard-safe residue avoids both shallower
blockers and all five unit masks, and Section 1 makes it safe from the
depth-four blocker as well.  Since powers of thirteen are coprime to seven
and fourteen, it lies off every equality endpoint and thickens to an
uncovered interval.  This contradicts (2), so `(2,3,4)` is empty.

This theorem does not follow from THM-2213, nor does it settle the four
remaining scalar depth-four profiles `(1,1,4)`, `(1,2,4)`, `(1,3,4)`, and
`(2,2,4)`.  It gives no uniform statement for larger deepest valuation.
LRC(14) remains open.  QED.
