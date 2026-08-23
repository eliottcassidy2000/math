---
id: THM-3700
title: "Minimal equal-step two-weight adjunction Pluecker gate"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Enlarge the
  seven-function collision-ring span of THM-3698 by A^2B and C^4, the minimal
  target monomials of weights -5 and 4.  The resulting nine functions occupy
  the equal-step weights {-5,-2,1,4}.  Their bracket-to-one bivector fibre is
  a seven-plane in a 36-dimensional exterior square, but it misses Gr(2,9):
  two explicit Pluecker minors already force 0=-7/24.  Hence no two linear
  combinations in this span form a Darboux pair.  This is an exact finite
  compression model, not an all-degree closure of the equal-step 3x4 cell or
  a counterexample to JC(2).
source: root / jc-sparse-direct-search audit, 2026-08-22
audit: >
  PASS.  The independent audit reproduced all nine functions and weights,
  linear independence, pair indexing, 48 coefficient rows, ranks 29/29, the
  seven free wedge coordinates, all 126 Pluecker coordinates, both decisive
  minors, their -7/24 contradiction, the collision, and bracket orientation.
depends_on:
  - THM-3686-y0-collision-normalization-and-bracket-anatomy
  - THM-3698-y0-collision-seven-function-pluecker-compression-no-go
related:
  - THM-3695-y0-collision-ring-danielewski-embedding-and-seven-piece-floor
  - THM-3699-y0-consecutive-four-weight-three-by-four-nonentry
script: 04-computation/jacobian_y0_equal_step_two_weight_adjunction_thm3700.py
output: 05-knowledge/results/jacobian_y0_equal_step_two_weight_adjunction_thm3700.out
script_sha256: 822a0e8db94286e1e5a358f258738b78b843004657002025340c1bbb2d2d118c
output_sha256: b22792a799ee8c3376b388b87d0f1aaf35402e404243fb9d6288a9fe520c8d82
hash_basis: LF-normalized bytes
---

# THM-3700 -- the minimal equal-step enlargement still misses the Grassmannian

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  THM-3698 showed
that one new homogeneous weight can never suffice.  This theorem performs the
first nonvacuous two-new-weight compression test and extracts a two-minor hand
certificate from the complete Pluecker calculation.

On the source plane use

```text
A=3z(2-x^2z),       B=2xz(2-x^2z),       C=x(1-x^2z). (1)
```

Order the nine collision-ring functions as

```text
F=(C,BC^2,AC^3, A,B^2,A^2C^2,ABC, A^2B,C^4).          (2)
```

They are linearly independent in `C[x,z]`.  Their weights are

```text
(1,1,1, -2,-2,-2,-2, -5,4),                           (3)
```

so the distinct weight set is the step-three progression

```text
{-5,-2,1,4}.                                           (4)
```

The last two functions are the minimal target-degree monomials in their new
weight directions.  All nine still retain the source collision

```text
(x,z)=(1,0),(-1,2).                                    (5)
```

## 1. The affine bivector fibre

As in THM-3698, define

```text
L:Lambda^2 span(F) -> C[x,z],
L(sum_(i<j) w_ij F_i wedge F_j)=sum_(i<j)w_ij{F_i,F_j}. (6)
```

There are `binom(9,2)=36` wedge variables.  Expanding `L(w)=1` in source
monomials gives

```text
48 coefficient rows,
rank(coefficient matrix)=rank(augmented matrix)=29.     (7)
```

Thus the bracket-to-one locus is a nonempty affine seven-plane.  In the
lexicographic pair order its free coordinates are

```text
{w_7,w_15,w_18,w_25,w_26,w_30,w_35},                  (8)
```

where

```text
w_18=w_(2,6),                    w_35=w_(7,8).          (9)
```

The rational three-bracket identity from THM-3686 is an explicit point of
this fibre, so nonentry below is purely decomposability debt.

## 2. Two Pluecker minors are already incompatible

A bivector in `Lambda^2 C^9` is one wedge `P wedge Q` exactly when all
`binom(9,4)=126` Pluecker quadrics vanish.  Substitution of the affine
seven-plane into those quadrics leaves 106 distinct nonzero equations, whose
complete Groebner basis is `[1]`.

The mechanism is much smaller.  Two substituted quadrics are

```text
Delta_(0,1,3,4)=-17(w_18-2w_35)/56,                   (10)

Delta_(0,1,3,6)=
 -(2304w_18^2-9216w_18w_35+9216w_35^2+343)/1176.      (11)
```

Decomposability would make `(10)` zero, hence

```text
w_18=2w_35.                                            (12)
```

After `(12)`, every quadratic term in `(11)` cancels and

```text
Delta_(0,1,3,6)=-343/1176=-7/24 !=0.                  (13)
```

This contradiction proves that the affine fibre does not meet the cone over
`Gr(2,9)`.  Equivalently,

```text
P,Q in span_C(F)  ==>  {P,Q}!=1.                       (14)
```

Since every `F_i` identifies the two points in `(5)`, a positive solution
would have been a genuine planar Jacobian counterexample.  None exists in
this span.

## 3. What this changes, and what it does not

The calculation is an independent exact view of a low-target-filtration
sector.  THM-3695's cited height floor also excludes any counterexample whose
two outputs both have such low target filtration, but it does not expose the
affine bivector fibre or the rank-two obstruction `(10)--(13)`.  Those data
are useful for choosing and auditing larger compression spaces.

This theorem closes only the particular nine-dimensional span `(2)`.  It
does not close arbitrary coefficient modules at the four weights `(4)`, an
all-degree `3 x 4` support word, or any higher-filtration enlargement.  No
counterexample and no proof of `JC(2)` is claimed.

## Reproduction

```bash
python3 -B 04-computation/jacobian_y0_equal_step_two_weight_adjunction_thm3700.py
python3 -B -O 04-computation/jacobian_y0_equal_step_two_weight_adjunction_thm3700.py
```

Both commands must byte-match the frozen transcript.  **QED.**
