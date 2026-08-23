---
id: THM-3697
title: "Two-or-three-pure-x-Q normalized two-by-three support closure"
status: >
  PROVED + FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  In the normalized planar chart P=x plus exactly two distinct nonlinear
  monomials and Q=y plus exactly three distinct nonlinear monomials, suppose
  both P supports have positive x exponent and at least two Q supports are
  pure x.  Then Jac(P,Q) is not constant.  In the exactly-two branch no
  hypothesis is made on the six cross determinants, and every one of the 68
  exact principal-open saturation charts is the unit ideal.  The all-pure-x
  branch follows from an exact coordinate shear.  Pure-y-P branches and JC(2)
  remain open; no counterexample is claimed.
source: jc-sparse-direct-search / 2026-08-22
audit: >
  PASS -- root independently checked all eleven labels and their linked
  activities, partition necessity, every affine exit, the complex superset,
  the exact two-chart distinctness cover, 38,000=190*20*10 census, and both
  controls.  Root also supplied and checked the all-pure-x coordinate-shear
  argument.  No typed dual is inferred.
depends_on: []
related:
  - THM-3418-one-monomial-nonlinear-fiber-keller-classification
  - THM-3551-one-ray-planar-jacobian-mate-no-go
  - THM-3689-fully-transverse-two-by-two-sparse-support-gate
  - THM-3690-complete-normalized-two-by-two-sparse-support-closure
  - THM-3692-ordinary-two-by-three-shift-root-peeling
  - THM-3694-one-pure-x-q-two-by-three-support-closure
script: 04-computation/jc2_two_pure_x_q_two_by_three_support_closure_thm3697.py
output: 05-knowledge/results/jc2_two_pure_x_q_two_by_three_support_closure_thm3697.out
script_sha256: e3584d45bec94f69a6db9e35b4e715e34b16e31f394187d21c6143c096a73863
output_sha256: d481825c1708ec6192f80a9835a4a68a97f5a0bcdbc0ea509a50da5e834d2add
semantic_sha256: 5c6b4e2cb727dc0621c30faa2f8e5f6595d3b92f56247712b161ac75a29aecf2
hash_basis: raw LF bytes for files; affine-exit ledger, saturated Groebner bases, hostile debt, and bounded census for semantic hash
---

# THM-3697 -- two or three pure-`x` terms leave no two-by-three Keller packet

**PROVED + FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**
This theorem closes the next typed axis cell after THM-3694.  Unlike the
half-integral residue in that theorem, the present cell has no point over the
complex numbers after all activity and distinctness conditions are imposed.
The terminal all-pure-`x` cell then collapses under a coordinate shear.

Work over a characteristic-zero field `K`.  Write `X^(r,s)=x^r y^s` and
`|(r,s)|=r+s`.

## 1. Statement

Let

```text
A,C,B in Z_{>=0}^2,                  k,l in Z_{>=2},
A != C,                              B,(k,0),(l,0) pairwise distinct,
|A|,|C|,|B| >=2,                     k != l,
A_1 C_1 B_2 !=0.                                         (1)
```

For arbitrary nonzero coefficients `a,c,b,d,e in K*`, put

```text
P=x+a X^A+c X^C,
Q=y+b X^B+d x^k+e x^l.                                  (2)
```

Then

```text
Jac(P,Q) is not constant.                                 (3)
```

The same conclusion holds if all three nonlinear supports of `Q` are
distinct pure-`x` supports:

```text
Q=y+b x^k+d x^l+e x^m,              k,l,m distinct and >=2.             (3a)
```

No cross-transversality is assumed.  No member of

```text
det(A,B), det(A,(k,0)), det(A,(l,0)),
det(C,B), det(C,(k,0)), det(C,(l,0))                       (4)
```

is required to be nonzero; every vanishing pattern compatible with the shared
pure-`x` directions is included.  The theorem does not cover a pure-`y`
support in `P`, nor does it infer a statement obtained by interchanging the
two components.  Sections 2--4 prove the exactly-two branch `(2)`; Section 5
proves `(3a)`.

## 2. Eleven labels and exact activity

The eleven possible terms of `Jac(P,Q)-1` have exponent buckets

```text
A-e_1, C-e_1, B-e_2, (k,-1), (l,-1),
A+B-e_1-e_2, A+(k,0)-e_1-e_2, A+(l,0)-e_1-e_2,
C+B-e_1-e_2, C+(k,0)-e_1-e_2, C+(l,0)-e_1-e_2.            (5)
```

The first three divergence labels are active because their coefficients are
`a A_1`, `c C_1`, and `b B_2`.  The fourth and fifth are inactive because the
two pure-`x` monomials have zero `y` derivative.  Each of the six cross labels
is active if and only if its determinant in `(4)` is nonzero: its coefficient
is the product of that determinant with two nonzero source coefficients.
Thus activity has no hidden coefficient exception in characteristic zero.

Every active term in `(5)` has positive total degree.  Consequently a constant
Jacobian would be exactly `1`, and every active bucket would have multiplicity
at least two.  For each of the `64` declared cross-activity masks, enumerate
all set partitions of the three fixed active divergence labels and the
declared-active cross labels in which no block is a singleton.  This produces
exactly

```text
11,155 activity/partition systems.                         (6)
```

The enumeration is a necessary support condition only; it grants arbitrary
coefficient cancellation inside every repeated bucket.  Ruling out even this
relaxed universe therefore rules out the full coefficient equations.

## 3. Affine exits

For every partition, impose the two axis equations for the pure-`x` supports
and the two coordinate equalities within each bucket block.  Exact rational
affine-linear solving removes `11,121` systems with the ordered ledger

```text
inconsistent                                      9,504
A=C                                                1,233
two Q supports equal                                 294
|B| forced below 2                                    34
|A| forced below 2                                    24
|C| forced below 2                                    24
declared-active determinant forced to zero              8.              (7)
```

The ordering assigns only the first diagnostic to a system that may have
several defects.  Each listed exit contradicts `(1)` or the declared activity
mask.  Exactly `34` affine residual systems remain.

## 4. Lossless principal-open cover

On each residual, impose every declared-inactive determinant from `(4)` as an
equation.  Let `S_0` be the product of

```text
- A_1, C_1, and B_2;
- all declared-active determinants in (4);
- |U|(|U|-1) for U in {A,C,B,(k,0),(l,0)};
- k-l.                                                     (8)
```

For nonnegative integer exponent vectors, the total-degree factors in `(8)`
are nonzero exactly when the corresponding support is nonlinear.  The factor
`k-l` enforces distinctness of the two pure-`x` terms, while `B_2 !=0`
automatically distinguishes `B` from both.  Hence this replaces the integer
lattice by a larger complex algebraic universe without discarding an
admissible integer packet.

The remaining condition `A!=C` is the union of the two principal opens

```text
A_1-C_1 !=0                         or A_2-C_2 !=0.         (9)
```

For each witness `w` in `(9)`, introduce a fresh inverse variable `s` and add

```text
s S_0 w - 1 = 0.                                           (10)
```

This is saturation by a checked product, not division by a possibly zero
factor.  The two charts cover all of `A!=C`.  Their exact ledger is

```text
34 residual systems x 2 witness charts = 68 saturations,
0 identically empty witness products,
68 unit ideals,
0 nonunit saturated ideals.                               (11)
```

Thus no complex exponent packet satisfies even the relaxed support
conditions.  In particular no packet in `(Z_{>=0}^2)^5` satisfies `(1)` and
the no-singleton condition forced by a constant Jacobian.  This proves `(3)`.

## 5. The all-pure-`x` branch is a shear equation

It remains to prove `(3a)`.  Write

```text
Q=y+g(x),                           q=Q,
P_tilde(x,q)=P(x,q-g(x)).                                 (12)
```

The change `(x,y) -> (x,q)` is a polynomial automorphism with Jacobian one.
By the chain rule,

```text
P_x=(P_tilde)_x+g'(x)(P_tilde)_q,      P_y=(P_tilde)_q,
Q_x=g'(x),                             Q_y=1,

Jac(P,Q)=(P_tilde)_x.                                    (13)
```

All nonconstant terms in the normalized Jacobian have positive total degree,
so a constant Jacobian would equal one.  Equation `(13)` then gives

```text
P_tilde=x+H(q),                     P=x+H(y+g(x))          (14)
```

for a univariate polynomial `H`.  If `deg(H)=r>=2`, the leading term of
`H(y+g(x))` contributes a nonzero pure-`y` monomial of degree `r`; no lower
power of `y+g(x)` can cancel it.  This contradicts the hypothesis that both
nonlinear supports of `P` have positive `x` exponent.  If `r=1`, `(14)` has
a nonzero linear `y` term, contradicting the normalized linear part `P=x+...`.
Thus `H` is constant, and the absence of a constant term in `(1)--(2)` forces
that constant to be zero.  Both alleged nonlinear `P` supports disappear,
again a contradiction.  This proves `(3a)` without a degree bound.

## 6. Independent census and boundary controls

An independent direct census through total degree six constructs the active
buckets directly from exponent vectors.  It uses neither activity partitions,
affine solving, nor Groebner saturation.  The exact typed universe is

```text
25 nonlinear monomials,
38,000 choices of {A,C}, B, and {k,l},
0 no-singleton survivors.                                 (15)
```

A close hostile packet explains why a product-like axis branch almost works:

```text
P=x+a x^2+c x^3,
Q=y+bxy+d x^2+e x^3.                                     (16)
```

Here the two pure-`x` terms in `Q` are bracket-invisible because `P` is pure
`x`, and

```text
Jac(P,Q)-1
 =(2a+b)x +(3c+2ab)x^2 +3bc x^3.                          (17)
```

Choosing `b=-2a` and `c=4a^2/3` cancels the first two layers but forces the
terminal debt `-8a^3 x^3`, nonzero under the chart hypotheses.  This is a
near-survivor, not a counterexample.

The support hypotheses are essential.  The tame pairs

```text
P=x+alpha y^2,                 Q=y+delta P^2,
P=x,                           Q=y+d x^2+e x^3             (18)
```

both have Jacobian one; the first drops one `P` support, and the second drops
both.  They are positive controls that the computation does not confuse
support sparsity itself with an obstruction.

## 7. Reproduction and scope

Run

```bash
python3 -B 04-computation/jc2_two_pure_x_q_two_by_three_support_closure_thm3697.py
python3 -O -B 04-computation/jc2_two_pure_x_q_two_by_three_support_closure_thm3697.py
```

Both streams must agree byte-for-byte with the stored output.  The companion
uses exact rational polynomial arithmetic, checks every finite-exact count in
`(6)--(11)` and `(15)`, verifies the hostile and tame controls, and rejects
inactive Python `assert` statements.

Together with THM-3692 and THM-3694, this closes the entire normalized
two-by-three chart in which both `P` divergence labels are active: Sections
2--4 handle exactly two inactive `Q` divergence labels, Section 5 handles all
three, THM-3694 handles exactly one, and THM-3692 handles none.  This does
**not** prove any cell with a pure-`y` `P` support or a typed dual.  No Keller
counterexample is constructed here.  **QED.**
