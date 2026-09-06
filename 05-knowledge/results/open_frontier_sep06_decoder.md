# Cross-component divisor cancellation supplies nonunit decoder closure

**Status: PROVED ELEMENTARILY ON THE STATED DOMAIN + INDEPENDENTLY AUDITED.
FINITE-EXACT controls.** This is a sufficient closure criterion for actual
two-component equality entries. It neither decides all entries nor proves
LRC(14). The 2+11 specialization is inherited, rather than claimed anew.

Root independently audited the full native and cross-divisor proofs and
source: PASS. The audit checked the physical-box pigeonhole, cleared
coefficient, exact radius, strict grid interior, hereditary weak/strict
split, and connected-tree ceilings. The separate Laurent-lane referee
also independently accepted Sections 2--3.

## 1. Inheritance and connection contract

The full native object is **THM-3818**,
`01-canon/theorems/THM-3818-scaled-inert-cubeclass-support-two-pair-packet.md`,
Sections 6.4--6.5. The exact two-generator box is inherited from
[the twelfth signed-box theorem](overnight12_20260906_lrc_gcd_semigroup.md),
and complete support orientations from
[the fourteenth general decoder](overnight14_20260906_lrc_general_decoder.md).
The closest closure mechanism is
[the fifteenth larger-unit theorem](overnight15_20260906_lrc_larger_unit.md).
It covers smaller size at most five when the larger primitive component
contains 1. The present result instead uses a divisor shared by a pair in
the larger component and one speed in the smaller component.

The inherited subset-gcd ceilings are
[the joint-shadow theorem](lrc14_joint_shadow_empty_core_next_sep06.md):
in a primitive strict counterexample, subset sizes 12,11,10,9,8,7 have gcd
at most 1,2,4,9,30,90 respectively. The exact endpoint-walk sidecar remains
[the fourteenth walk theorem](overnight14_20260906_lrc_endpoint_walk.md):
collective gcd does not determine endpoint gcd, as 6 -> 2 -> 3 shows.

The concept board is:

* full actual decoder equality and its physical box;
* the exact signed radius, with both support orientations retained;
* the cross-component divisor `gcd(D,v)`;
* the larger component's whole maximum and safe-arc length;
* hereditary scale ceilings and connected-tree height bounds;
* endpoint-walk cancellation, whose missing divisibility is still unpaid.

The source-to-target map sends a selected larger pair and one smaller label
to a bounded three-coordinate relation. It preserves the distinguished
coefficient and physical coordinate. The retained common divisor lowers
that coefficient and raises the scale forced by failure of the box test.
Cyclic gluing then converts the scale lower bound into a common safe time.
Replacing the whole larger maximum by the selected pair maximum would lose
the Lipschitz sidecar. Replacing the selected pair gcd by a collective gcd
would lose exactly the walk cancellation factor. Neither replacement is made.

The 2+11 version is already a specialization of the twelfth native test.
For 3+10, 4+9, and 5+8, the nonunit larger-component classes below are not
stated in the cited entry or larger-unit notes. The result is an extension
of those mechanisms; there is no external priority claim.

## 2. Domain and exact native criterion

Fix `Q=91^6=567869252041`. Let thirteen distinct positive integers n be
primitive and satisfy `sum(n)<=Q^2`. Build the actual all-scale inert-cube
decoder graph of THM-3818. Its primitive edge ratios have admissible sum at
most 356, and primitive edge coefficients at most 355. Suppose the graph
has exactly two connected components, written in primitive form as

    n = tV union gU,   |V|=a <= b=|U|,   a+b=13,
    gcd(V)=gcd(U)=gcd(t,g)=1.

Let W be the rational span of **all** zero integer rows supported on at most
three labels, with coefficient height at most Q. Assume `W=V_dec`, the span
of actual decoder edge rows. A proposed partition or a selected list of
negative support tests is insufficient for this assumption.

Write `K=max U`. Select distinct `u<w` in U and any `v` in V, and set

    D=gcd(u,w),  A=u/D, B=w/D,
    R=Q(A+B)-(A-1)(B-1),
    delta=gcd(gD,tv),  c=gD/delta,  x=tv/delta.       (1)

**Native criterion.** The conditions

    c<=Q,      a delta R >= 7(b+1) K v                (2)

imply a common physical time with every clearance strictly greater than
1/14. In particular, u,w need not include a unit or the maximum.

**Proof.** Every internal primitive pair height is at most Q. Indeed, if
one exceeded Q, adjoin a label in the other component. The `(Q+1)^3`
nonnegative coefficient triples have at most `Q^3+1` weighted sums because
the physical sum is at most `Q^2`. A collision gives a nonzero bounded
relation. Its distinguished coefficient cannot vanish, since that would
be an internal pair relation of height at most Q. The mixed relation lies
outside the direct sum of the two weighted internal kernels, contradicting
entry. Thus `B<=Q` in (1).

The exact signed-box theorem gives every integer in `[-R,R]` as
`rA+sB`, with `|r|,|s|<=Q`. If `c<=Q` and `x<=R`, then

    c(tv) - r(gu) - s(gw) = 0                        (3)

is an excluded bounded mixed relation. Its coefficient at the smaller label
is positive. Consequently entry and the first part of (2) force

    x>R,     t>delta R/v.                            (4)

The inherited lower-runner supplier gives a U-safe phase of clearance at
least `1/(b+1)`. Its closed 1/14-safe arc has length

    ell = a/[7(b+1)K].                              (5)

Here the entire U maximum is necessary. The t physical lifts of a V-safe
phase, whose V-clearance is at least `1/(a+1)`, form a complete translated
t-grid in the U clock because `gcd(t,g)=1`. By (2)--(4), `t>1/ell`, so a
grid point lies in the interior of the safe arc. Both components then have
strict clearance greater than 1/14. This proves the native criterion.

The equality boundary in (2) is allowed because (4) is strict. A negative
box membership result outside the central interval need not give a scale
bound stronger than (4); the exact radius is not replaced by the full
support radius `Q(A+B)`.

## 3. A nonunit endpoint criterion and automatic divisibility classes

Now take `w=K`, put `d=gcd(D,v)`, and let M_b be the inherited gcd ceiling
for b labels:

| b | 12 | 11 | 10 | 9 | 8 | 7 |
|---|---:|---:|---:|---:|---:|---:|
| M_b | 1 | 2 | 4 | 9 | 30 | 90 |

**Cross-divisor criterion.** In the domain above, the two conditions

    M_b D/d <= Q,
    lcm(D,v) <= aQ/[7(b+1)]                          (6)

imply weak LRC(14) safety. If also `g<=M_b`, the conclusion is strict safety.

To prove this, a strict counterexample must have `g<=M_b`, since the gcd
of its b physical larger labels equals g. The divisor d divides delta, so

    c = gD/delta <= M_b D/d <= Q.

Moreover `R>=QB=QK/D`. Indeed,
`R-QB=QA-(A-1)(B-1)>0` for `B<=Q`. Therefore

    delta R/v >= dQK/(Dv) = QK/lcm(D,v)
                   >= 7(b+1)K/a.

The native criterion contradicts strict failure. If g is above M_b, the
inherited gcd theorem already supplies weak safety; strictness is not
silently inherited from that separate branch.

**Divisibility corollary.** If `a<=5` and an endpoint-pair gcd D in the
larger primitive component divides any speed v in the smaller primitive
component, the full row is weakly safe. The larger component may be unitless.

Here `d=D`, so the coefficient gate in (6) is simply `M_b<=Q`. The connected
tree bound gives every primitive smaller speed

    v <= max V <= 355^(a-1).

For each `a=1,...,5` this is at most `aQ/[7(b+1)]`, and `lcm(D,v)=v`.
For a singleton V the statement has its literal meaning: v=1, so D=1.
The divisor need not be one for the other splits. In particular, at 5+8 it
may be as large as `355^4=15882300625`. This bound is conditional on the
stated divisibility, not a claim that every endpoint gcd satisfies it.

A second convenient corollary, without shared divisibility, is an endpoint
pair with D no larger than the following uniform cutoffs:

| split | max V tree bound | sufficient D ceiling |
|---|---:|---:|
| 1+12 | 1 | 6240321451 |
| 2+11 | 355 | 38086468 |
| 3+10 | 126025 | 175558 |
| 4+9 | 44738875 | 725 |
| 5+8 | 15882300625 | 2 |
| 6+7 | 5638216721875 | no positive automatic ceiling |

The formula is `min(floor(Q/M_b), floor(aQ/[7(b+1)355^(a-1)]))`.
For 2+11 this is weaker than the already inherited uniform endpoint result
after relabeling the two-component; no new claim is attached to that row.
For 6+7, (2) and (6) remain available. In particular, `D|v` and
`v<=60843134147` suffice; the connected tree bound alone does not force it.

For a non-endpoint selected pair, the same proof replaces the second line
of (6) by `lcm(D,v) K/w <= aQ/[7(b+1)]`. Its exact radius version (2) can
be stronger. This retains the physical whole-maximum factor K/w.

## 4. Actual entries, hostile boundaries, and computation

The source verifies four small nonunit examples of splits 3+10, 4+9, 5+8,
6+7. Each has larger maximum 30 and selected endpoint pair (3,30), so D=3;
the smaller shape includes 3. No larger shape contains 1. The 5+8 larger
shape `(3,5,6,9,10,12,15,30)` has no endpoint pair of gcd below 3, so even
the displayed uniform D<=2 corollary does not account for that control.
All four small controls pass every inherited necessary joint-shadow profile
for every body subset of sizes 7 through 12, not just the maximum-clock
ceilings. Their gcd maxima in descending subset size are `(1,1,1,3,3,3)`.
They are safe positive controls; passing necessary profiles is not an
assertion of failure or of novelty for those individual rows.

A larger nonunit 5+8 control is built from the incoming fifteenth star:

    V=(1013861907,1036929995,1060507875,1084612635,1096868145),
    v=1013861907,
    U=V union {3v,4v,9v},   g=1,   t=22Q+1.

The larger endpoint pair `(v,9v)` has D=v, which divides the smaller label
v and exceeds the old 11+2 uniform constant 76388115. The present split is
5+8; comparing these numbers does not assert a contradiction to the old
theorem. The actual full row has physical sum

    66123361394683029422040 < Q^2=322475487413604782665681.

Both graph components are verified from the actual atlas. The original
five-star edges have sum 356, and the three new links from v have sums
4,5,10. All U coordinates are coprime to t. Hence a two-small/one-large
crossing would have distinguished coefficient at least t>Q. A crossing in
the opposite orientation is impossible because `t min(V)>2Q max(U)`.
Independently, all 220 mixed supports pass the exact negative box test.
This supplies actual equality without inferring it from a selected subset.
This larger control is **already closed by inherited subset-gcd bounds**:
its gcd maxima for subset sizes 12 down to 7 are
`(1,183,183,33123,33123,5929017)`. It fails 1558 inherited profile checks.
It demonstrates large-divisor actual entry, rather than a residual shape
unpaid by the old gcd theorem. This distinction is checked in the source.

The explicit full-row phase

    12493123545929 / 24986247089806

has minimum clearance `3131122695665/24986247089806 > 1/14`.
It is obtained by gluing the simple component phases `eta=1/2` and
`zeta=1/2+1/(24v)`. The source verifies the phase directly against all
thirteen physical coordinates.

Hostiles retain the coefficient gate (`Q=3`, physical pair 8,12 and
distinguished 1); the missing boundary point 14 in the `(2,3),Q=3` box;
the failed replacement of the whole maximum by a selected maximum
(`{1,...,10,85}`); a noncoprime 4/2 clock; and the 6->2->3 collective-gcd
loss. An exact toy comparison `lcm(10,10)=10 <= 5*200/63 <100` confirms
that retaining d gives a strictly stronger arithmetic condition than Dv.
None of these is called an unsafe LRC row.

The finite universe consists of 67 complete signed boxes with budgets 2..8,
29484 coprime scale/divisor controls, and five named actual entries. Every
mixed support in both orientations is checked in the entries, for a total
of 1034 supports. All inherited profile checks are replayed using the frozen
joint-shadow JSON (SHA256
`935f3f687b6d7c89cc099e536f536238fd753bcc4c1747906d213cef387ca93f`).
The computation is not a census of the physical box.

[Source](../../04-computation/open_frontier_sep06_decoder.py),
[normal output](open_frontier_sep06_decoder.out):

```
python3 -B 04-computation/open_frontier_sep06_decoder.py
python3 -B -O 04-computation/open_frontier_sep06_decoder.py
```

Normal and optimized outputs are byte-identical; 125030 explicit gates pass.
Semantic digest:
`1dd73a034cfd521ba0a7c68c3f35388e27ea8a47b153d7e7e81c8d04980c3cb9`.
Raw SHA256 pins:

* source: `0e9ff0799542c55a88584a0e2bac3d66fb031035a144a9fa4ac4219a153bf17d`;
* output: `991d8161ab38777db0ab3788bcbf9ef3f455fe7c4d0fcc9610612e4d79ed36b9`.

## 5. Exact remaining obligation

Failure of the sufficient conditions is not failure of LRC. The residual
structural question is whether actual component connectivity, the complete
hereditary profiles, and full equality entry force an endpoint pair and
opposite label satisfying (6), or a non-endpoint pair satisfying (2).
The walk factor J remains the missing information in converting collective
divisor control into this selected-pair predicate. The cheapest next test
is therefore a profile-admissible connected shape with every cross-divisor
score in (6) failing, while retaining all actual equality tests and the
physical box. A new height census or another entry algorithm would not
answer that structural question.
