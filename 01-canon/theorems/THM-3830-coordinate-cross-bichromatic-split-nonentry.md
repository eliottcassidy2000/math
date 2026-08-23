---
id: THM-3830
title: "A transverse coordinate cross cannot realize the cubic-pseudoplane bichromatic split"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For
  every root alpha of the five-slope polynomial in THM-3827, every nonzero
  scalar c, and arbitrary-degree q,d, the normalized unimodular-row ansatz
  k=c+xyq, h=alpha k+xy cannot satisfy the sidecar factor equation while d
  contains exactly one of x,y.  Equivalently, every solution has x|d iff
  y|d.  This closes only the elementary transverse-coordinate model of a
  bichromatic fixed fibre; disjoint and non-coordinate factorizations, the
  second row, and the Keller equation remain open.
source: jc_sparse_direct_search / THM-3827 bichromatic-passport constructive lane, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (root / jc-cohn-boundary and
  thm3830_hostile_audit, 2026-08-23).  The audits check the equal-origin
  fibre-product unit argument, determinant-one
  completion, both axis restrictions, reflection of x/y divisibility under
  finite scalar extension, and the exact reduced-cross scope.  Its custom
  Sylvester/Bareiss arithmetic reproduces all three characteristic-zero
  invariants; an exhaustive F_11 boundary universe has 39,930 candidates and
  exactly 60 solutions, all zero or the forced nonzero scalar.  Line-unit,
  trivial-idempotent, and CRT-connectedness controls agree with THM-3827.
  Primary and independent companions contain no effective assertion gap;
  normal and optimized raw LF streams match both frozen transcripts.
depends_on:
  - THM-3827-generic-fibre-genus-floor-for-nonlinear-cubic-plane-atlases
related:
  - THM-3822-nonlinear-cubic-plane-atlas-sl2-and-punctured-arm-gate
script: 04-computation/jc2_coordinate_cross_bichromatic_split_thm3830.py
output: 05-knowledge/results/jc2_coordinate_cross_bichromatic_split_thm3830.out
script_sha256: 84bce55eb0344aad1dea319adc5fd39c60fbbff175c88cf9f574bd569ae36913
output_sha256: 6d062ed9bfcd13005fd55d790fa5bf8d64005337b29572b5a502f94e3b50b529
semantic_sha256: 3ce7d10d25f49bbebe1ce8e6d7530783fe7a1a5e9045ef9b37a3678d3213d2c5
independent_script: 04-computation/jc2_coordinate_cross_bichromatic_split_thm3830_independent_audit.py
independent_output: 05-knowledge/results/jc2_coordinate_cross_bichromatic_split_thm3830_independent_audit.out
independent_script_sha256: 3a5abd63b8b60a2c7814bbe185f494fccc9cfcb829fab6f0b94dfa0305eafe4b
independent_output_sha256: e56f9c8d45d3f41289ba56de1de7ce67e9651c590360456d50a19ea5072a696a
independent_semantic_sha256: 7144dd1572ee7003d17d1faff15cc675a5116febfc9c6387d5b54167633c90b5
hash_basis: raw LF bytes
---

# THM-3830 -- a coordinate cross cannot carry opposite sidecar signs

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  Work over a
field `K` of characteristic zero.  After a finite scalar
extension, fix a root `alpha` of

```text
a(z)=(7z^2+3)(3z^3+7z^2+1).                             (1)
```

For `c in K*` and arbitrary `q,d in K[x,y]`, put

```text
k=c+xyq,                 h=alpha k+xy,                  (2)
A_5=(7h^2+3k^2)(3h^3+7h^2k+k^3),
B_3=(h+k)(2h+k)(3h-k).                                  (3)
```

If

```text
A_5=d(kB_3+h^2d),                                       (4)
```

then

```text
x divides d    if and only if    y divides d.            (5)
```

In particular, `(4)` cannot make `d` select exactly one irreducible
component of the transverse fixed-slope fibre

```text
h-alpha k=xy=0.                                         (6)
```

The quantifier over `q,d` is all-degree; no bounded-support hypothesis is
present.

## 1. The five slope constants are harmless

Set

```text
b(z)=(z+1)(2z+1)(3z-1).                                 (7)
```

The exact one-variable certificates are

```text
disc(a)=353831803500,
Res_z(a,b)=-31298700,
Res_z(a,z b)=93896100.                                  (8)
```

Thus every root `alpha` of `a` is simple and satisfies

```text
alpha!=0,                    b(alpha)!=0.                (9)
```

This is the only characteristic-zero arithmetic used below.  Notice also
that on the fibre `h=alpha k`,

```text
kB_3=k^4 b(alpha).                                      (10)
```

## 2. Why `(2)` is the whole normalized coordinate-cross grammar

Suppose more generally that `(k,h)` is a unimodular row and
`h-alpha k=xy`.  Reducing a Bezout identity for the row modulo `(xy)` shows
that the class of `k` is a unit in

```text
K[x,y]/(xy).                                            (11)
```

The normalization of `(11)` is `K[x] times K[y]`, and the image of `(11)`
in this normalization consists exactly of pairs whose values at the two
origins agree.  A unit therefore restricts to a unit on both affine lines,
hence to the same nonzero scalar `c`.  Consequently

```text
k=c+xyq                                                  (12)
```

for a polynomial `q`.  Conversely, every row `(2)` has the explicit
determinant-one completion

```text
m=q/c,                  C=(1+alpha q)/c,
Ck-mh=1.                                                   (13)
```

Thus the theorem is not discarding a hidden polynomial choice inside the
exact normalized fibre `(6)`.

## 3. Restriction to the unselected component

Assume for contradiction that

```text
x|d,                         y does not divide d.         (14)
```

Write a bar for restriction to `y=0`.  Then

```text
bar(k)=c,             bar(h)=alpha c,
bar(A_5)=c^5 a(alpha)=0,
bar(kB_3)=c^4 b(alpha).                                  (15)
```

Equation `(4)` becomes the identity in the domain `K[x]`

```text
0=bar(d)(c^4 b(alpha)+alpha^2 c^2 bar(d)).                (16)
```

The condition `y` does not divide `d` says exactly that `bar(d)!=0`.
Equations `(9)` and `(16)` force

```text
bar(d)=-c^2 b(alpha)/alpha^2 in K*.                       (17)
```

But `x|d` implies `x|bar(d)`, contradicting `(17)`.  Interchanging `x,y`
proves the reverse implication and hence `(5)`.

The mechanism is geometric as well as algebraic.  Along any fixed-slope
fibre, unimodularity makes `k` a unit, and `(10)` makes the two square-root
values `+kB_3` and `-kB_3` distinct everywhere.  Components carrying
opposite signs therefore cannot meet.  The two axes in `(6)` do meet, so
they cannot be the bichromatic pair required by THM-3827.

## 4. Exact scope and the surviving construction lane

This theorem rules out only the elementary ansatz in which the reducible
spectral member is the transverse coordinate cross.  It does **not** exclude

```text
* d containing both x and y, or neither;
* a fixed-slope fibre with disjoint components;
* a non-coordinate or higher-component factorization;
* solving the remaining intrinsic D/second-row equations; or
* satisfying a nonzero constant planar Jacobian.                 (18)
```

The nearest counterexample-positive replacement is therefore a disconnected
fixed fibre, for example a product of comaximal factors, rather than another
intersecting cross.  Equivalently, `(11)` has only trivial idempotents, so its
connected coordinate cross cannot realize THM-3827's canonical nontrivial CRT
split.  Thickened and merely set-theoretic crosses remain outside the literal
reduced-fibre scope.  No Jacobian counterexample is claimed.  **QED.**
