---
id: THM-3698
title: "Seven-function Pluecker compression gate and two-new-weight adjunction floor"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For the seven collision-ring functions underlying
  the additive bracket identity of THM-3686, the affine fibre of bivectors
  mapping to 1 is an exact two-plane in a 21-dimensional exterior square.
  Its intersection with the decomposable Gr(2,7) cone is empty: after
  substitution, the 35 Pluecker quadrics generate the unit ideal.  Thus no
  two linear combinations of these seven functions form a Darboux pair,
  although two brackets do sum to 1.  The functions occupy only weights 1
  and -2; THM-3695 therefore also excludes every enlargement by functions of
  only one additional homogeneous weight.  At least two new distinct weights
  are needed before compression can possibly evade the inherited support
  floor.  Arbitrary larger spans and JC(2) remain OPEN.
source: root / 2026-08-22
audit: >
  PASS.  The independent audit replayed the source formulas, bracket
  orientation, retained collision, seven-dimensional span, 19/19 affine
  ranks, all 35 Pluecker equations, and the three bipartite rank-one charts.
  It also extracted the two-minor hand contradiction recorded in Section 2.
depends_on:
  - THM-3686-y0-collision-normalization-and-bracket-anatomy
  - THM-3695-y0-collision-ring-danielewski-embedding-and-seven-piece-floor
related:
  - THM-3696-y0-collision-ring-three-branch-conductor-and-graded-modules
script: 04-computation/jacobian_y0_seven_function_pluecker_compression_thm3698.py
output: 05-knowledge/results/jacobian_y0_seven_function_pluecker_compression_thm3698.out
script_sha256: 62da571a7b603b35d67ac301d75632d2b5a1da97bc44e6e7a9634e60a65a45b7
output_sha256: 5479b6d2c80af72023a2c29a81e7a550e02ca122216bbb1a318f775300f000a8
hash_basis: LF-normalized bytes
---

# THM-3698 -- the additive identity is an affine plane missing the Grassmannian

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  This theorem turns one-bracket compression into
a linear fibre followed by a decomposability test.  It gives a geometric
explanation of why the smallest known additive identity does not itself
produce a counterexample, and it identifies the minimum kind of enlargement
that is worth testing.

Use the collision-ring generators on the source plane from THM-3686:

```text
A=3z(2-x^2z),       B=2xz(2-x^2z),       C=x(1-x^2z). (1)
```

Put

```text
U=(C,BC^2,AC^3),
V=(A,B^2,A^2C^2,ABC),
E=span_C(U,V).                                           (2)
```

The seven displayed pullbacks are linearly independent, so `dim E=7`.  Every
one belongs to the collision ring and hence takes the same value at the two
source points `(1,0)` and `(-1,2)`.

## 1. Linearize first, impose decomposability second

Fix the ordered basis `F_0,...,F_6=(U,V)`.  The bracket defines a linear map

```text
L: Lambda^2 E -> C[x,z],
L(sum_(i<j) w_ij F_i wedge F_j)
  =sum_(i<j) w_ij {F_i,F_j}.                            (3)
```

Thus the bracket-to-one locus

```text
mathcal A={w in Lambda^2 E : L(w)=1}                   (4)
```

is affine-linear.  By contrast, `w` comes from one pair `P,Q in E` exactly
when it is decomposable:

```text
w=P wedge Q.                                           (5)
```

In coordinates, `(5)` is the affine cone over `Gr(2,7)`, cut out by the 35
quadrics

```text
w_ij w_kl-w_ik w_jl+w_il w_jk=0,
                         0<=i<j<k<l<=6.                (6)
```

This separates the linear bracket equation from the genuinely nonlinear
part of the counterexample search.

## 2. Exact affine fibre and Pluecker gate

Expanding `(3)--(4)` in source monomials gives

```text
21 wedge variables,
29 coefficient rows,
rank(L-system)=rank(augmented system)=19.               (7)
```

Hence `mathcal A` is a nonempty affine two-plane.  Exact rational elimination
expresses all 21 wedge coordinates in two free parameters.  Substitute that
parametrization into every quadric in `(6)`.  Their Groebner basis over `Q` is

```text
[1].                                                    (8)
```

There is also a transparent certificate behind `(8)`.  Every same-side
`U wedge U` and `V wedge V` coordinate vanishes on the whole affine fibre.
Writing its `U wedge V` block in parameters `a,b` gives

```text
M(a,b)=
[ 1/6,   -1/4, (1-9a)/6, -12b/7 ]
[ 8b/7,  3b/28, -15b/14,     7/4 ]
[ a,     -13/8,        0,       b ].                   (8a)
```

The associated skew matrix has rank twice `rank(M)`.  A nonzero decomposable
bivector would therefore force `rank(M)=1`.  But two of its minors are

```text
Delta_(01|01)=17b/56,
Delta_(01|03)=(2304b^2+343)/1176.                      (8b)
```

The first forces `b=0`, where the second is `7/24`, a contradiction.  Thus
the hand certificate and the complete Pluecker Groebner certificate agree.

Therefore

```text
mathcal A intersect Cone(Gr(2,7)) is empty.             (9)
```

Because `(8)` is a unit-ideal certificate over `Q`, the conclusion holds
after extension to `C`.  Equations `(5)` and `(9)` prove

```text
P,Q in E  ==>  {P,Q} != 1.                             (10)
```

The companion also performs an independent bipartite replay.  If
`P in span(U)` and `Q in span(V)`, the `3 x 4` coefficient-matrix fibre has
11 rows, rank 10, and affine dimension two.  A nonzero rank-one matrix has a
nonzero `U` coordinate, so three normalized pivot charts cover its rank-one
cone.  Every chart has Groebner basis `[1]`.

## 3. The fibre is positive; its minimum skew rank is four

The emptiness in `(9)` is not a disguised failure of the linear equation.
The rational identity

```text
1={C,-B^2/4+A/6+A^2C^2/6}
  +(7/4){BC^2,ABC}-(13/8){AC^3,B^2}                    (11)
```

is a point of `mathcal A`; its alternating matrix has rank six.  More
sharply, THM-3686 gives two brackets with every entry in `E` whose sum is one,
so `mathcal A` contains a bivector of rank at most four.  Equation `(9)` rules
out rank two, while a nonzero skew rank is even.  Thus

```text
min{rank(w):w in mathcal A}=4.                         (12)
```

Equivalently, the restricted additive bracket length of one in the span `E`
is exactly two.  This does not determine the additive bracket length in the
whole collision ring.

## 4. A uniform adjunction floor from the grading

The seven functions in `(2)` have an unexpectedly rigid weight profile:

```text
wt(C)=wt(BC^2)=wt(AC^3)=1,
wt(A)=wt(B^2)=wt(A^2C^2)=wt(ABC)=-2.                  (13)
```

Now enlarge `E` by any number of homogeneous collision-ring functions all of
one additional weight `r`.  Every element of the enlarged span has at most
the three active weights

```text
{-2,1,r}.                                              (14)
```

Two alleged Darboux outputs would therefore have at most six active pieces
in total.  THM-3695 excludes every collision-ring Darboux pair through total
support six.  Consequently no such one-new-weight enlargement can contain a
Darboux pair.

This prevents a large but uninformative computation: adjoining one target
monomial to `(2)` can never work, regardless of its degree or coefficients.
The first meaningful compression test must add at least two distinct new
grading weights, or one genuinely nonhomogeneous function carrying at least
two new weights.  By THM-3695's repaired coordinatewise floor, the first live
support shape is `3 x 4`; this theorem neither closes it nor constructs a
counterexample.

## Reproduction

```bash
python3 -B 04-computation/jacobian_y0_seven_function_pluecker_compression_thm3698.py
python3 -B -O 04-computation/jacobian_y0_seven_function_pluecker_compression_thm3698.py
```

Both commands must agree byte for byte with the frozen transcript.  **QED.**
