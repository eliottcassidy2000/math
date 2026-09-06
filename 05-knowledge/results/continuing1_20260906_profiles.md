# Ranked component capacities for all-height inert tails

**Status: PROVED relative to the inherited envelope; FINITE-EXACT;
INDEPENDENTLY AUDITED.** The [independent review](continuing1_20260906_lrc_profiles_audit.md)
rebuilds the complete finite head and literal odd-scale carriers. The large
body scout remains separately producer FINITE-EXACT.
No theorem ID or external priority claim. LRC(14) remains **OPEN**.

**Concurrent recovery.** Incoming commit8e560f214 contains the same ranked
envelope and packing consumer in
[open_frontier_sep06_components.md](open_frontier_sep06_components.md),
plus its sharp square-length consequence. The shared envelope here is an
independent concurrent proof, not a separate new theorem. The paired-body
rank-four simplification and the larger bounded scout are retained as
additional scoped checks; neither supplies a new physical body closure.

## Inheritance and precise connection

The closest proved mechanism is the five-profile mass/width envelope in
[creative_20260906_inert_pareto.md](creative_20260906_inert_pareto.md).
The recovered sidecar is the full open-component multiset in
[THM-4052 / affine-component width escape](../../01-canon/theorems/THM-4052-lrc14-affine-component-width-escape-cones.md).
The canonical hostile is forgetting component locations: even the complete
length multiset does not determine insertion, as the two physical bases in
[THM-799 / transverse-tooth cap](../../01-canon/theorems/THM-799-good-set-state-transverse-tooth-cap-lemma.md)
show. The corrected near miss is merging touching strict danger intervals;
their shared safe point must survive. No correction to those suppliers is claimed.

The map takes a compact body safe set and an open tail-spoilage set to their
decreasing lists of positive component lengths. It preserves necessary
containment inequalities at every rank, but loses positions, gap sizes and
endpoint owners. The target theorem is a sufficient escape certificate,
not an equivalence. The required arithmetic sidecar is the common odd tail
scale. The cheapest test is exact partial-sum comparison of the complete
46-pair head, followed by a physical body search.

## 1. The ranked containment obstruction

For a finite nonnegative decreasing list v define T_k(v)=sum of its first
k entries, padding with zeros. Let G be a nonempty compact subset of the
circle with finitely many connected components, and C a proper open subset
with finitely many connected components. Let u and v be their positive
component-length lists. If G is contained in C, then

    T_k(u) < T_k(v) whenever T_k(u)>0.                 (1)

To prove this for k no larger than the number of positive G-components,
select its k longest components. Each lies in a single C-component; at
most k C-components are used. Within each used open component, the selected
closed pieces leave positive length at its ends. Their total length is
strictly smaller than that component's length. Summing gives (1). For
larger k use the last positive index and monotonicity of T_k(v). Isolated
points are retained geometrically but contribute no positive length.

Thus any positive partial sum reaching or exceeding the corresponding
open-carrier partial sum certifies escape. This includes equality. It
contains both the largest-width and total-mass certificates.

Write u <=_w v when T_k(u)<=T_k(v) at every positive integer rank. For an
integer g>=1 let D_g(v) repeat each v_i/g exactly g times. Then

    T_(gr+s)(D_g(v))=T_r(v)+(s/g)*v_(r+1), 0<=s<g.   (2)

Consequently weak majorization is preserved by D_g, and D_g(v)<=_w v.
The latter also follows because dividing one interval into equal pieces
cannot increase any sum of a fixed number of the largest lengths.

## 2. The same five profiles control every rank

Let p<q be positive coprime odd 3-units, with every prime in p+q congruent
to 2 mod3 and of exponent at most two. Let C_(p,q) be the literal open set
of y for which the two half-lifts y/2 and (y+1)/2 are both spoiled by p,q
at strict danger radius 1/14. Its entire positive length list is weakly
majorized by one of precisely these five maximal lists:

| Pair | Complete decreasing component lengths |
|---|---|
| (7,13) | 1/49 repeated twice |
| (19,25) | 37/3325 twice; 23/3325 twice; 9/3325 twice |
| (5,41) | 2/287 repeated six times |
| (5,53) | 2/371 repeated six times; 9/1855 twice |
| (1,67) | 2/469 repeated ten times |

The inherited sharp mass bound is 20/469, and every component has length
at most 2/(7q). For q>67 all lengths are below 2/469 and their sum is below
20/469. Any such list is weakly majorized by ten copies of 2/469: use the
entry bound for k<=10 and the mass bound thereafter. Thus only q<=67 remains.
There are exactly 46 admissible pairs in that head. Exact opposite-parity
interval intersections and every partial sum give precisely the five
lists displayed. The finite head is exhaustive under the stated filters.
They are mutually incomparable because their masses increase and their
largest components decrease. Hence the theorem has no height restriction.

## 3. A sufficient physical-body certificate

For any nonempty finite set H of positive integer speeds, put

    G_H={y mod1: ||hy||>=1/14 for every h in H},

and let u be its positive component list. Let g be positive and odd.
If for each of the five lists v there is a rank k such that

    0 < T_k(u) >= T_k(D_g(v)),                         (3)

then every row 2H union {gp,gq}, for an admissible pair p,q as above, has
a phase of clearance at least 1/14. Its number of distinct speeds is the
actual size of this union; the LRC14 application uses an eleven-body.

Indeed strict failure would put the compact G_H in the inverse image of
C_(p,q) under y -> gy. Since g is odd, multiplication by g permutes the two
physical half-lifts, so this is the correct spoiled carrier. Its lengths
are D_g(v_actual). Section 2 and (2) place these below D_g(v) for one of
the five lists. That list's inequality (3) contradicts (1).

A certificate at g=1 works at every positive odd g, by D_g(v)<=_w v.
This is an all-height sufficient criterion. It is not a universal lower
bound on G_H, and it does not turn a chosen physical split into a decoder
component partition.

## 4. What is and is not a gain

The ranked test can exceed the mass/width test. For example the abstract
reflection-symmetric positive length list

    u=(1/100,1/100,1/100,1/100,1/2000,1/2000)

has mass 41/1000 and maximum 1/100. It fails the (19,25) mass/width
disjunction, but its first four lengths sum to 1/25>120/3325. Mass pays
the (7,13) profile and width pays each of the last three profiles. These
lengths can be realized by disjoint closed arcs on the circle. They are
**not asserted to be G_H for an eleven-speed body**.

If H contains an even speed, 1/2 is dangerous, while 0 is always dangerous.
Reflection therefore pairs every positive component with a different one.
The positive lengths occur in equal pairs. In this case (3) at g=1 is
exactly the old five mass-or-width disjunctions, except that the (19,25)
row gains the third alternative

    T_4(u)>=120/3325.                                (4)

For each target equal-value block, its partial sums are linear. A sorted
paired source is linear between consecutive even ranks. On each interval
between target even breakpoints, the bound by the source maximum and mass
already controls uniform target blocks; for (19,25) the only extra even
breakpoint is four. More directly: target (19,25) has breakpoints 2,4,6;
T_2(u)<=2 max(u), T_6(u)<=sum(u); target (5,53) has breakpoints 6,8,
with T_6(u)<=6 max(u), T_8(u)<=sum(u). The remaining targets are uniform.
Checking breakpoints and interpolation proves the assertion, including
strict/equality boundaries. The full criterion (3) applies without the
even-speed condition and at every odd scale.

## 5. Bounded physical search and stopping reason

The [body scout](../../04-computation/continuing1_20260906_profiles_body_scout.py)
examines every one of the 2,496,144 eleven-subsets of {1,...,24}, both as H
and as 2H. It rebuilds their closed safe intervals on the exact integer grid
of denominator 14*lcm(1,...,24)=74,959,204,320. No body in this universe
gains a certificate from (4) while missing the five old mass/width gates.
For the undilated bodies, exactly four pass the necessary mass and width
filters; their largest top-four sum is 387/14896, below 120/3325. The
[frozen scout output](continuing1_20260906_profiles_body_scout.out) gives
all four sets and their exact profiles. This is a bounded negative result
about the rank-four improvement, not an all-height redundancy theorem.

The main [source](../../04-computation/continuing1_20260906_profiles.py)
and [output](continuing1_20260906_profiles.out) check the whole finite
46-pair head, every partial sum, dilations 1,3,5, an additional tail bank,
462 paired-list controls on a declared seven-value alphabet, and strict
packing/equality/isolated-point hostiles. All 9,724 gates pass. Both programs
have byte-identical normal and optimized outputs:

    python -B 04-computation/continuing1_20260906_profiles.py
    python -B -O 04-computation/continuing1_20260906_profiles.py
    python -B 04-computation/continuing1_20260906_profiles_body_scout.py
    python -B -O 04-computation/continuing1_20260906_profiles_body_scout.py

The missing coordinate for further progress remains actual body geometry.
No new physical eleven-body family closure is claimed. This is a recorded
stopping reason for broader blind body enumeration; a future continuation
should recover a structural body family or retain phase incidence.
