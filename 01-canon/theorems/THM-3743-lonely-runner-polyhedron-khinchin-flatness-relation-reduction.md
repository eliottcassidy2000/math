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
  ||a||_1<=(7/6)Flt(12).  The explicit general estimate recorded by
  Averkov--Hofscheier--Nill gives ||a||_1<=356, improving THM-2144's total
  coefficient cap 367 by 11.  An l1-minimal choice is a Graver relation; its
  support-two branch is a finite atlas of 19,314 reduced ratios, while its
  higher-support branch is genuinely multiway.  No LRC(14) conclusion follows.
source: root + lrc_khinchin_archaeology / 2026-08-23
audit: >
  PASS.  An independent audit checked the primal quotient versus dual
  annihilator normalization, integral rather than Euclidean-unit width,
  full-dimensionality, the closed-boundary convention, primitive-line
  saturation, the exact 356 floor, the Graver split, and comparison with
  THM-2051/2052/2144/2169.  Two independently written exact companions agree;
  normal and optimized output is byte-identical.
depends_on:
  - THM-2144-anisotropic-selberg-kraft-relation-box
related:
  - THM-2051-fejer-bv-small-relation-alternative-for-lrc14
  - THM-2052-finite-height-forces-high-rank-bounded-relation-code
  - THM-2164-relative-packet-rank-harvesting
  - THM-2169-bounded-relation-on-every-lrc-deletion
  - THM-2190-basis-safe-floor-and-height-500-rank-six-harvest
  - THM-3718-lrc-complete-atom-orbit-defect-saturation-and-semantic-boundary
script: 04-computation/lrc14_khinchin_flatness_relation_thm3743.py
output: 05-knowledge/results/lrc14_khinchin_flatness_relation_thm3743.out
script_sha256: 391cfe3b920c82bf5698d8f2e191111097b31c181c35c43aca4274826cce02e5
output_sha256: b6902e00563091e2ec1de0bb5514d71e2a9f315b416dd5d8b77e0957fe409041
secondary_script: 04-computation/lrc14_khinchin_flatness_relation_audit_20260823.py
secondary_output: 05-knowledge/results/lrc14_khinchin_flatness_relation_audit_20260823.out
secondary_script_sha256: 16358a45d3ee6fc6c4ad6a6fb5780e1cf37a245e17ea78a051794abea1a68397
secondary_output_sha256: 7e657f35b943704a14b3557120482f8637fcac0b4e6d192db337bde9c33ba8bb
hash_basis: raw LF bytes
---

# THM-3743 -- an LRC(14) counterexample would be flat in an exact relation direction

**PROVED FROM CITED INPUTS + VERIFIED-EXACT + INDEPENDENTLY
HOSTILE-AUDITED.**  The cited inputs are the standard flatness theorem and the
explicit general bound recorded by Averkov--Hofscheier--Nill.  The
quotient-lattice identification, width formula, Graver split, restricted body
class, and comparison with the current LRC relation canon are proved here.

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

The cited width convention, flatness theorem, and explicit estimate are pinned
under
[`CORE-PAPERS.md`](../../05-knowledge/reference/CORE-PAPERS.md#khinchin-flatness-and-an-explicit-general-bound);
the quotient and coefficient algebra below are proved in-repo.

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

The explicit general estimate recorded by Averkov--Hofscheier--Nill is

```text
Flt(d)<=sqrt((d+1)(2d+1)/6) d^(3/2).                  (15)
```

At `d=12`, this gives

```text
Flt(12)^2<=((13)(25)/6)12^3=93600,
Flt(12)<=60 sqrt(26),
||a||_1<=70 sqrt(26)<357.
```

Since the left side is integral,

```text
||a||_1<=356.                                          (16)
```

This is an exact numerical consequence of the displayed published inequality,
not an appeal to an unspecified asymptotic improvement.

## 5. What is new, and what current canon already does better

THM-2144 applies to every zero-measure safe set, hence to every hypothetical
counterexample.  Its cap profile with ten coordinates bounded by `28` and
three bounded by `29` supplies a relation with

```text
||a||_1<=10*28+3*29=367.                               (17)
```

Thus `(16)` improves the existing explicit total-coefficient cap by `11`.
It does not improve THM-2051's sparse-support conclusion: the two statements
control different grades of the relation code.

The mechanisms are nevertheless orthogonal:

- THM-2144 is a signed Selberg--Kraft Fourier certificate in an anisotropic
  coefficient box.
- THM-2051 gives the much stronger support type `3..5`, but with height
  `2^20`.
- The flatness proof selects a supporting direction of the actual quotient
  zonotope and therefore comes with a geometric slice and active-facet
  question, although `(14)` alone does not record those sidecars.

### 5.1 The minimizing direction is a Graver relation

Choose a nonzero relation `a` with minimum `l1` norm.  If it decomposed
conformally as `a=b+c` with nonzero integer kernel vectors `b,c`, then

```text
||a||_1=||b||_1+||c||_1,
```

so both summands would be strictly shorter relations, a contradiction.  Thus
`a` is conformally indecomposable: it is a Graver element of the one-row speed
matrix.

If its support is two, on speeds `v_i,v_j`, primitiveness forces, up to sign,

```text
(a_i,a_j)=(v_j/g,-v_i/g),      g=gcd(v_i,v_j),
(v_i+v_j)/g<=356.                                  (18)
```

There are exactly `19,314` unordered distinct coprime positive ratios whose
reduced numerator plus denominator is at most `356`.  This finite branch can
be routed through THM-778's ordered continued-fraction word and its sidecars.

The support-at-least-three branch is genuine rather than disguised pair
arithmetic.  On speeds `(3,4,5)`, the relation `(1,-2,1)` has `l1=4`, smaller
than every primitive pair relation on those speeds.  This branch is a bounded
multiway partition identity and belongs with the relation-code/Fourier canon.

### 5.2 Exact joins with the existing relation code

The zero-safe premise of THM-2169 holds for a hypothetical counterexample.
Choose an index `i` in the support of `a`.  The THM-2169 relation on the
`i`-deletion has `i`th coordinate zero, whereas `a_i` is nonzero, so the two
relations are independent.  Flatness therefore upgrades the deletion theorem
to a global rank-two packet, but no further rank follows automatically.

There is also a precise conditional join with THM-2052.  In its rank-eleven
branch, let `W` be the rational span of the support-at-most-three,
height-`91^6` relations.  Then

```text
a notin W  => rank(W+Qa)=12 and the cofactor terminal applies;
a in W     => W itself contains a nonzero vector of l1 norm <=356.       (19)
```

In the first branch, Hadamard gives the explicit finite speed cap
`floor(356*3^(11/2)*(91^6)^11)`.  The second branch is the real structural
survivor: classify rank-eleven star relation spaces containing a short dense
Graver vector.  This paragraph inherits every conditional dependency already
declared by THM-2052; it is not an unconditional rank-eleven theorem.

This suggests the restricted constant

```text
Flt_zon,13,2tor(12)                                    (20)
```

obtained by taking the supremum only over centrally symmetric
twelve-zonotopes with at most thirteen generators and centre of order at most
two modulo the lattice.  Tautologically

```text
Flt_zon,13,2tor(12)<=Flt(12).                          (21)
```

A useful improvement would be an explicit bound for `(20)` together with a
description of its minimizing direction and slice—not merely a smaller
generic number.

## 6. Hostile controls and missing sidecars

Both

```text
v_AP=(1,2,...,13),
v_far=(1,2,...,12,5460)                                (22)
```

are primitive, distinct-speed rows and obey

```text
(2,-1,0,...,0).v_AP=(2,-1,0,...,0).v_far=0,
||(2,-1,0,...)||_1=3.                                  (23)
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
 physical root, owner, phase, ties, and temporal word). (24)
```

The first five entries are geometric; the remaining entries are exactly the
semantic debt left by THM-3718.  Dropping the intercept is especially
dangerous: parallel slices with the same normal can have different lattice
occupancy.

The best scoped next test is therefore not another relation census.  It is:

1. compute or bound the restricted constant `(20)`;
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

The main companion checks the side length, exact corner widths on hostile
signed directions, the exact floor `floor(70 sqrt(26))=356`, the AP/far-AP
control, the `367` comparison, the `19,314` pair atlas, a genuine triple
hostile, and the two-torsion centre.  The independently written secondary
companion exhausts projected vertices in small dimensions and checks the
rank-eleven/cofactor join.  Neither script claims to prove the cited flatness
theorem.

The result is a lawful geometry-of-numbers reduction and a new restricted
flatness target.  It proves neither a recursive slice theorem nor LRC(14).
