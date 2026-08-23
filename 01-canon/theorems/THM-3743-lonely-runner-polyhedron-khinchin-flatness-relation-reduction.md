---
id: THM-3743
title: "Lonely-runner quotient zonotope and Khinchin-flatness relation reduction"
status: >
  PROVED FROM CITED INPUTS + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  A primitive hypothetical LRC(14) counterexample makes the projection of
  [1/14,13/14]^13 along its speed line a full-dimensional lattice-free
  twelve-zonotope.  Its dual lattice is exactly the integer relation lattice,
  and the width in relation direction a is (6/7)||a||_1.  Khinchin flatness
  therefore forces a primitive nonzero speed relation with
  ||a||_1<=(7/6)Flt(12); the classical explicit Flt(d)<=d^(5/2) bound gives
  ||a||_1<=581.  The numerical cap is weaker than THM-2144's existing 367,
  but the proof exposes a restricted centrally symmetric, thirteen-generator,
  two-torsion-centred zonotope class and a recursive slice interface.  No
  LRC(14) conclusion follows.
source: root + lrc_khinchin_archaeology / 2026-08-23
audit: >
  PASS.  An independent audit checked the primal quotient versus dual
  annihilator normalization, integral rather than Euclidean-unit width,
  full-dimensionality, the closed-boundary convention, primitive-line
  saturation, the 581 floor, and comparison with THM-2051/2144/2164/2190.
  Normal and optimized exact companions are byte-identical.
depends_on:
  - THM-2144-anisotropic-selberg-kraft-relation-box
related:
  - THM-2051-fejer-bv-small-relation-alternative-for-lrc14
  - THM-2164-relative-packet-rank-harvesting
  - THM-2190-basis-safe-floor-and-height-500-rank-six-harvest
  - THM-3718-lrc-complete-atom-orbit-defect-saturation-and-semantic-boundary
script: 04-computation/lrc14_khinchin_flatness_relation_thm3743.py
output: 05-knowledge/results/lrc14_khinchin_flatness_relation_thm3743.out
script_sha256: 13aec7024ae297783f3cfc43d4d290a6ff515550dad6f03fc8bb3401661cbdaa
output_sha256: 5d13bc9940ed8050a1f68e2e37c9dc271cbd422b094de27548b5b0786491712b
hash_basis: raw LF bytes
---

# THM-3743 -- an LRC(14) counterexample would be flat in an exact relation direction

**PROVED FROM CITED INPUTS + VERIFIED-EXACT + INDEPENDENTLY
HOSTILE-AUDITED.**  The cited inputs are the standard flatness theorem and
its explicit classical bound.  The quotient-lattice identification, width
formula, restricted body class, and comparison with the current LRC relation
canon are proved here.

This is **Khinchin's geometry-of-numbers flatness theorem**.  It is not the
metric continued-fraction theorem, Khinchin's constant, Khintchine
transference, the Khintchine inequality, or Wiener--Khinchin autocorrelation.

## 1. Inheritance and conventions

Let

```text
v=(v_1,...,v_13) in Z_(>0)^13                         (1)
```

be a vector of distinct speeds with `gcd(v_1,...,v_13)=1`.  Primitive
normalization does not change its speed line or lonely-runner maximum, but it
makes the quotient lattice visibly saturated.

The closest geometric mechanism is the lonely-runner polyhedron formulation
of Beck--Hoşten--Schymura.  The closest proved in-repo relation mechanisms are:

- THM-2051, which forces a support-`3..5` relation of height at most `2^20`;
- THM-2144, which forces an anisotropic relation with total cap `367`;
- THM-2164/2190, which force several independent bounded relations.

The canonical hostiles are the arithmetic progression `{1,...,13}` and its
far lift `{1,...,12,5460}`: both have the tiny relation
`(2,-1,0,...)`, yet only the first has the AP geometry.  Relation existence by
itself cannot encode loneliness.  The least-used sidecars here are the
supporting facets and translate of the flat slice.

Use the **closed** target box

```text
B=[1/14,13/14]^13 subset R^13.                         (2)
```

A boundary point represents a valid witness with minimum distance exactly
`1/14`.  Replacing `(2)` by an open box would silently change the conjecture.

## 2. Exact quotient-zonotope formulation

Let

```text
Q=R^13/Rv,                  pi:R^13 -> Q,
Lambda=pi(Z^13),            K=pi(B).                   (3)
```

Since `v` is primitive, `Lambda` is a full lattice in the twelve-dimensional
space `Q`.

### Lemma 2.1 (closed LRC witness iff quotient lattice point)

There is a time `t in R` satisfying

```text
||t v_i||_(R/Z)>=1/14        for every i              (4)
```

if and only if

```text
K intersection Lambda != empty.                       (5)
```

*Proof.*  Condition `(4)` is equivalent to the existence of
`z in Z^13` with

```text
t v-z in B.                                            (6)
```

Equation `(6)` holds exactly when some `x in B` and `z in Z^13` have
`pi(x)=pi(z)`, which is `(5)`.  QED

Thus a hypothetical counterexample has

```text
K intersection Lambda=empty.                          (7)
```

The linear map `pi` is surjective, so it maps the interior of `B` to a
nonempty open subset of `Q`.  Hence `K` is full-dimensional.  It is also the
centrally symmetric quotient zonotope

```text
K=pi((1/2,...,1/2))
  +sum_(i=1)^13 [-3/7,3/7] pi(e_i).                    (8)
```

It has at most thirteen generators in dimension twelve.  Its centre class is
two-torsion modulo `Lambda`, because

```text
2 pi((1/2,...,1/2))=pi((1,...,1)) in Lambda.           (9)
```

These three restrictions—central symmetry, thirteen generators, and a
two-torsion centre—are extra structure absent from the general flatness
constant.

## 3. The dual lattice is exactly the speed-relation lattice

A linear functional on `Q` pulls back to a linear functional
`x|->a.x` on `R^13` which annihilates `v`.  It is integral on `Lambda` if and
only if it is integral on every `pi(e_i)`, equivalently `a in Z^13`.
Therefore

```text
Lambda*={a in Z^13:a.v=0}.                             (10)
```

This is the exact normalization.  The primal lattice is `pi(Z^13)`; the
intersection `Z^13 intersection v^perp` is its dual, not another expression
for the primal lattice.  Interchanging them introduces false determinant
factors.

### Lemma 3.1 (exact width formula)

For every `a in Lambda*`,

```text
width_a(K)=(6/7)||a||_1.                               (11)
```

*Proof.*  The value of `a` on `pi(x)` is well-defined because `a.v=0`.
Consequently

```text
width_a(K)
 =max_(x in B) a.x-min_(x in B) a.x
 =(13/14-1/14) sum_i |a_i|
 =(6/7)||a||_1.                                        (12)
```

Each positive coefficient selects the high endpoint for the maximum and the
low endpoint for the minimum; each negative coefficient does the reverse.
QED

Lattice width minimizes `(11)` over **integral covectors**
`a in Lambda*\{0}`.  There is no division by the Euclidean norm `||a||_2`.
That alternate normalization would make `(11)` false.

## 4. Khinchin flatness forces a bounded relation

Define `Flt(d)` as the supremal lattice width of a full-dimensional hollow
convex body in dimension `d`, under the integral-dual normalization used in
Section 3.  A body avoiding all lattice points, as in `(7)`, is in particular
hollow, so the flatness theorem applies to `K` and `Lambda`.

There is therefore a nonzero `a in Lambda*` satisfying

```text
width_a(K)<=Flt(12).                                   (13)
```

Divide `a` by the gcd of its coordinates.  The resulting vector is still an
integer relation and cannot have greater width.  Combining `(10)--(13)` gives
the main conclusion:

```text
There is primitive 0!=a in Z^13 with
a.v=0,                  ||a||_1<=(7/6)Flt(12).         (14)
```

The classical explicit estimate

```text
Flt(d)<=d^(5/2)                                          (15)
```

specializes to

```text
||a||_1<=(7/6)12^(5/2)=336 sqrt(3)<582.
```

Since the left side is integral,

```text
||a||_1<=581.                                          (16)
```

This numerical specialization deliberately uses the simple published
constant `(15)`, not an unspecified asymptotic improvement.

## 5. What is new, and what current canon already does better

THM-2144 applies to every zero-measure safe set, hence to every hypothetical
counterexample.  Its cap profile with ten coordinates bounded by `28` and
three bounded by `29` supplies a relation with

```text
||a||_1<=10*28+3*29=367.                               (17)
```

Thus `(16)` is numerically weaker than the existing explicit total-cap
certificate.  It must not be advertised as a new coefficient record.

The mechanisms are nevertheless orthogonal:

- THM-2144 is a signed Selberg--Kraft Fourier certificate in an anisotropic
  coefficient box.
- THM-2051 gives the much stronger support type `3..5`, but with height
  `2^20`.
- The flatness proof selects a supporting direction of the actual quotient
  zonotope and therefore comes with a geometric slice and active-facet
  question, although `(14)` alone does not record those sidecars.

This suggests the restricted constant

```text
Flt_zon,13,2tor(12)                                    (18)
```

obtained by taking the supremum only over centrally symmetric
twelve-zonotopes with at most thirteen generators and centre of order at most
two modulo the lattice.  Tautologically

```text
Flt_zon,13,2tor(12)<=Flt(12).                          (19)
```

A useful improvement would be an explicit bound for `(18)` together with a
description of its minimizing direction and slice—not merely a smaller
generic number.

## 6. Hostile controls and missing sidecars

Both

```text
v_AP=(1,2,...,13),
v_far=(1,2,...,12,5460)                                (20)
```

are primitive, distinct-speed rows and obey

```text
(2,-1,0,...,0).v_AP=(2,-1,0,...,0).v_far=0,
||(2,-1,0,...)||_1=3.                                  (21)
```

Their very different lonely-runner behaviour proves that neither `(14)`,
`(16)`, nor relation size alone is the target predicate.

For a recursive flatness argument, a selected direction must be upgraded to
the packet

```text
(positive and negative active facet sets,
 slice intercept modulo lattice hyperplanes,
 induced slice lattice and remaining generators,
 affine translate / two-torsion centre,
 physical root, owner, phase, ties, and temporal word). (22)
```

The first five entries are geometric; the remaining entries are exactly the
semantic debt left by THM-3718.  Dropping the intercept is especially
dangerous: parallel slices with the same normal can have different lattice
occupancy.

The best scoped next test is therefore not another relation census.  It is:

1. compute or bound the restricted constant `(18)`;
2. classify the possible active-facet sign patterns at a minimizing direction;
3. recurse on the induced eleven-dimensional lattice slice while retaining
   its translate; and
4. test whether that packet determines any THM-3718 owner/root/word class.

Failure of step 4 would give a precise stopping reason even if steps 1--3
sharpen the coefficient cap.

## 7. Exact companion and scope

Reproduce the arithmetic audit with

```bash
python3 -B 04-computation/lrc14_khinchin_flatness_relation_thm3743.py
python3 -B -O 04-computation/lrc14_khinchin_flatness_relation_thm3743.py
```

The companion checks the side length, exact corner widths on hostile signed
directions, the integer floor `floor(336 sqrt(3))=581`, the AP/far-AP control
relation, THM-2144's `367` comparison, and the two-torsion centre.  It does not
claim to computationally prove the cited flatness theorem.

The result is a lawful geometry-of-numbers reduction and a new restricted
flatness target.  It proves neither a recursive slice theorem nor LRC(14).
