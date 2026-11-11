# Two births: a boundary-event proof of the triangle-chain matching law

**PROVED IN NOTE, independently reviewed within this session.** The complete
finite coefficient checks below are **FINITE-EXACT**. The root reviewer
independently audited the coefficient grouping, expansion signs, endpoint
three-event walks, all four-event walks, double-event block cases, and the
nu=0 ESES pigeonhole; the root also reported the Berggren lane's independent
review before status promotion. Unselected endpoint labels are summed
within each group, not fixed before invoking the unmarked-triangle rotation.
The full arbitrary-k matching formula was subsequently proved and reviewed
in [TRIANGLE_ALL_K_PROOF.md](TRIANGLE_ALL_K_PROOF.md); this note retains the
earlier independent two-birth proof and its complete finite certificates.

Let T be an ordinal chain of r directed triangles, with arbitrary singleton
gaps g_0,...,g_r. Let nu be the maximum matching size in the path on these
triangles whose edge i is present precisely when g_i=0. Let C_2(T) be the
gcd of the Hamilton path counts of all tournaments formed by adjoining two
vertices while preserving every old edge.

The claim is

    C_2(T) = 3^max(r-2-min(2,nu),0).                    (1)

The constructor upper bound is already proved in
[TRIANGLE_FILTRATION.md](../fibres/TRIANGLE_FILTRATION.md): a new vertex
extends an adjacent triangle pair to an H41 graph, or a single triangle to
H5; use disjoint available pairs, then singles, then pad with sinks. The
all-sink two-birth extension has count3^r. We prove the matching lower bound.

## 1. The exact Boolean polynomial

Name the new vertices x and y. For each old vertex v, let X_v,Y_v indicate
the arcs x->v,y->v. Let Z indicate x->y. The Hamilton count is a multilinear
integer polynomial P in these 2n+1 orientation variables, of degree at most4:
a Hamilton path uses at most four edges incident to the two new vertices.
An edge in the opposite direction contributes1 minus the corresponding
variable. No path uses an edge twice.

The gcd of P's values on the Boolean cube equals the gcd of its integer
coefficients. One direction follows by evaluation; the other is Boolean
Möbius inversion

    [monomial S]P = sum_(A subset S) (-1)^(|S|-|A|) P(1_A).

It therefore suffices to prove divisibility of every coefficient.

Fix a coefficient and one of the two orders of x,y in the Hamilton path.
Deleting x,y leaves three ordered old paths P0,P1,P2, possibly empty.
Each old path visits the old strong components in their ordinal order.
For x before y, if all three paths are nonempty the four old-new factors are

    (1-X_end(P0)) X_start(P1)
      (1-Y_end(P1)) Y_start(P2).                      (2)

Empty P0 or P2 removes the corresponding factor. Empty P1 replaces its
two factors by Z. Reversing the new-vertex order swaps X,Y and replaces
that Z by1-Z. This is an exact partition of all Hamilton paths.

Group terms contributing to the fixed coefficient by the locations of the
boundary events (which old triangle, or which singleton), the identity of
the event, and the old vertex wherever that endpoint is selected by the
monomial. Unselected endpoint labels remain summed. Every positive factor
in (2) must be selected; an unselected endpoint can only come from a
negative factor's constant term. The expansion sign is constant within
each such group. Inside each old block, assignments to the old paths and
their internal path orders are independent of the assignments in the other
blocks. Thus a group factors into local integer counts, with one common
sign. A group with incompatible boundary conditions is zero.

Call a triangle marked if the coefficient monomial names any of its old
vertices. A triangle that is unmarked supplies a factor3: cyclically rotate
its three labels in the local covers, including any unselected endpoints.
There is no selected orientation variable in that triangle to distinguish
the rotation, and the action is free. This argument does not require the
triangle's vertices to occur consecutively in a full Hamilton path.

If the coefficient marks t triangles, it is therefore divisible by
`3^(r-t)`. Since t<=4, this recovers the general weak bound. Only t=3 and
t=4 require improvement.

## 2. Local boundary counts and availability walks

At a block, an old path is available if its prescribed start is at or before
that block and its prescribed end is at or after it. P0 is available from
the left until its end; P2 from its start to the right; P1 between its start
and end. A start event raises the number q of available paths by one; an
end event lowers it by one. Availability is only a permission to use a path
in a block; its local portion may be empty.

A triangle with exactly one boundary event, at a specified vertex, and
exactly two paths available at that event contributes **3**. The path with
the fixed start or end may contain respectively1,2,3 vertices of the
triangle. Each case has exactly one permitted cyclic order; the other path
contains the remaining vertices. These three cases exhaust the local covers.
For one or three available paths the corresponding count is1 or7, but only
the two-path value is needed.

An interval with q=0 cannot contain any old vertex. In particular, a zero
availability interval between distinct triangle blocks is possible only
when those triangles are consecutive and their intervening singleton gap is
zero. Such an interval certifies an edge in the zero-gap path defining nu.

When all four old-new boundary events occur in distinct blocks, the
availability walk begins and ends at1. Write S for a start and E for an end.
The possible sequences, together with the number available at each event,
are contained in the following table. Additional named-path constraints can
delete rows; they cannot add rows.

| Event order | Available counts at the events | Single-event two-path factors |
|:---|:---|:---|
|SSEE|2,3,3,2|first and fourth|
|SESE|2,2,2,2|all four|
|SEES|2,2,1,1|first and second|
|ESSE|1,1,2,2|third and fourth|
|ESES|1,1,1,1|none; requires two zero-availability intervals|

EESS makes availability negative and is impossible. A zero interval between
an end and a start in one of the last three rows must be empty of old
vertices; it cannot be ignored merely because the endpoints are named.

## 3. Four marked triangles

If t=4, all four old-new factors are selected, their endpoints lie in four
distinct triangles, and all three old paths are nonempty. In every row
other than ESES, at least two marked triangles have the local factor3 from
Section2. The ESES row requires two disjoint zero-gap adjacencies between
the four marked triangles, so it is impossible when nu<=1.

Consequently, for nu<=1, every such group is divisible by

    3^(r-4) * 3^2 = 3^(r-2).                            (3)

This is stronger than needed when nu=1. It is a termwise local factor,
not an asserted cancellation between unrelated paths.

## 4. Three marked triangles when nu=0

There are at least three selected old-new incidences. If x,y are consecutive
in the Hamilton path there are at most two such incidences, so that case
contributes zero here. If one new vertex is at an endpoint of the full path,
there are three boundary events, all selected and in the three marked
triangles. The availability walk begins at0 and ends at1, or vice versa.
With every internal triangle gap positive, the admissible three-event walk
has a marked event with two available paths: the alternative walk would
have a zero-availability interval between two distinct marked triangles.

It remains to treat four boundary events. There is at most one unselected
event, since three different marked triangles require at least three
selected incidences.

If all four events occur in distinct old blocks, every row in the preceding
table except ESES contains at least two single-event two-path positions.
At most one can be unmarked, so at least one is a marked triangle and
contributes3. In ESES, the two zero-availability intervals pair the four
event locations disjointly. At least three locations are marked triangles,
so one pair consists of two marked triangles. It would supply a zero-gap
edge, contradicting nu=0. This argument allows the fourth event to be a
singleton or an unmarked triangle.

Otherwise the four events occupy exactly three blocks, namely the three
marked triangles. One block contains two events and the other two contain
one each. Availability in every interval between these distinct triangles
must stay positive. If the double-event block is first, its net change is
either+2 or0 (a negative net change is impossible); in either case at least
one subsequent single event has exactly two available paths. The case of
the last block is symmetric. If the double-event block is in the middle,
the first single event must be S and the last E, and both have two available
paths. Thus in every case there is a marked single-event triangle with
local factor3. No ordering of two events inside a cyclic triangle is assumed.

Together with the unmarked triangles, every t=3 coefficient for nu=0 is
therefore divisible by `3^(r-3)*3=3^(r-2)`.

## 5. Assemble the lower bounds

If nu=0, coefficients with t<=2 have factor `3^(r-2)` by free rotations;
Sections3–4 give that same factor for t=3,4. If nu=1, coefficients with
t<=3 have factor `3^(r-3)`, and Section3 gives an even stronger factor for
t=4. If nu>=2, the general weak bound `3^max(r-4,0)` suffices.
Negative exponents are replaced by0. These are exactly the three cases
of the exponent in (1). The inherited disjoint constructors give the reverse
divisibility, completing the proof.

## 6. Complete finite polynomial certificates and hostile probes

`TRIANGLE_exact.cpp` independently constructs P by enumerating all ordered
three-path covers of the old source, multiplying the two new-vertex order
factors, and collecting coefficients. Each C3 has33 ordered three-path
covers:9 with all vertices in one bin,18 with occupancy2+1, and6 with
occupancy1+1+1. Each singleton has3 choices. Thus the complete cover universe
has size `33^r*3^sum(g_i)`. Endpoint-identical covers are aggregated without
discarding their integer multiplicity.

`TRIANGLE_check.py` computes complete C2 values for sixteen named sources:
all four binary internal-gap patterns for r=3, all eight for r=4, both
patterns for r=2, and two extra positive/exterior-gap controls. For r=4,
the internal-gap patterns000,001,010,011,100,101,110,111 give respectively

    C2 = 1,3,1,3,3,3,3,9.

The all-separated r4 source has15 old vertices, so this is a complete
coefficient certificate for all `2^31` two-birth orientation assignments.
It has14767 nonzero coefficients. The all-separated r3 source has4557
nonzero coefficients and C2=3, covering all `2^23` assignments. These are
complete finite-universe statements, not sampled gcds.

Six seeded polynomial evaluations per source (96 in total) are independently
checked by the literal subset Hamilton counter. Separate hostile k3 probes
give sample gcd3 on64 extensions of four separated triangles and sample
gcd9 on12 extensions of five separated triangles. These latter observations
are only finite samples and do not prove the k3 matching formula.

Reproduce with `python -B TRIANGLE_check.py --write`, then
`python -O -B TRIANGLE_check.py`. Coefficient hashes, minimal-valuation
coefficients, marked-triangle coefficient gcds, all literal evaluation
controls, and the actual k3 graphs are retained in TRIANGLE_receipt.json.
