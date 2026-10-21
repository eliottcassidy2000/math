# Two unit gaps still have an exact owner; the first three-record collision

**PROVED, pending independent audit:** for k>=3 the unified mass formula
in the frozen `zero_heavy.md`, section 6, remains valid throughout

    0<=delta<=2k-1.

This extends, rather than contradicts, that section's conservative
0<=delta<=2k-3 range. Its stated open two-unit-cell boundary is discharged
by the retained-port analysis below. The first actual cross-family owner
collision occurs at delta=2k, with an explicit fibre of size three.
No unrestricted class-count recurrence is claimed.

## Incoming inheritance and nonduplication

During this continuation, incoming `zero_heavy_next_20261020h/README.md`
(at commit e2654366789, read on main through 426f28797523) independently
proved the complete deficit-k formula and corrected the shared-endpoint
onset. Its k=2 proof already supplies the full-matching two-unit-gap
four/five-cell mechanism used in section 3 below. Those statements are
inherited and checked here, not claimed as newly discovered. The new
continuation is their use throughout the deficit corridor, the mixed
missing/unit four-cell port distinction, and the first cross-family
owner/correction at delta=2k. The incoming packet explicitly makes no
formula claim at k+1 or beyond. Its hashes are pinned in the addendum
manifest. The older heavy-deficiency-two owner remains separately
attributed, with its distinct parameter.

## 1. The unchanged grammar and formula

Keep the total-weight parameter and all definitions of the frozen packet.
Since delta<2k, the sharpened identity

    delta=(k+1)d+kh+U

allows only (d,h)=(0,0),(0,1),(1,0). These are full matching gaps, one
omitted matching gap, and a full-gap P3-plus-matching respectively. No
second omitted gap or extra backward-edge collision enters yet.
The binomial count is consequently unchanged. Set

    general: m=2k+1, E=F_k, J=F_k, P=(k-1)F_k;
    strong:  m=2k-1, E=I_k, J=(2k-3)I_{k-1}, P=P_k.

The class count is

    E [binom(m+delta-1,delta)
       -m binom(m+delta-k-1,delta-k)]
    +(m E-J) binom(m+delta-k-2,delta-k)
    +P binom(m+delta-k-3,delta-k-1).                     (1)

A negative lower index means the entire term is absent. Every nonabsent
binomial has its usual nonnegative-integer meaning. The earlier counting
argument is valid for all delta<2k. What needs extending is the exact
presentation-to-class quotient, not stars-and-bars.

## 2. Intrinsic cells are small but their ports must be retained

All weight>1 vertices still avoid backward arcs in every good order, so
their transitive order is intrinsic. Their incidence vectors partition
unit vertices into fixed ordered cells. Every exterior endpoint outside
an exceptional cell lies in a singleton cell, except for other exceptional
cells explicitly addressed below. Isomorphisms must retain all incidences
to these outside vertices.

For a full matching, at most two gaps are unit: three would cost
3(k-1)>2k-1 for k>=3. For the one-missing-gap family, the remaining deficit
is at most k-1, so at most one occupied gap is unit. For the P3 family the
remaining deficit is at most k-2, so every gap still has weight at least
two. Its previous rigid quotient proof therefore requires no change.

## 3. Full matching with two unit gaps

If the two unit gaps do not touch, treat their separate two- or three-
vertex cells. A two-vertex cell is ordered by its arc. A three-vertex
cell with no internal matched pair is transitive and has a unique order.
If its two endpoints are matched, it is C3, all its vertices have identical
outside incidences, and its rotations preserve the very same record.
This works even when matching edges join the two exceptional cells:
in that case neither relevant cell has an internal matched pair and
both local transitive orders fix their endpoint roles first. No circular
external-partner choice is used.

Two consecutive internal unit gaps give a cell

    e1, u1, e2, u2, e3,

where e's are endpoints and u's are gap vertices. There can be at most one
internal matching edge. With none, this is a transitive five-vertex cell,
whose order fixes all roles. With one, the unmatched endpoint has a fixed
external partner and is thereby a distinguished port of the cell. If the
internal pair is (e1,e2), (e1,e3), or (e2,e3), the marked endpoint's outdegree
inside this cell is respectively 0,2,4. Thus these three records cannot
collide. In the first and third cases the remaining possible C3 rotation
is internal and changes no record. In the middle case the marked endpoint
has degree two; its incidence on the other four vertices fixes their roles
(the unique middle mark in the ordered extreme-edge reversal).

If the consecutive unit gaps include an exterior gap, the cell has four
vertices, with two endpoints and two gaps. If the endpoints have outside
partners, they are distinguished and the two gap roles follow. If the
endpoints are matched, the cell is a source joined to C3, or C3 joined to
a sink, according to the fixed end. Its C3 rotations again preserve the
same record. These exhaust two-unit-gap cells, including strong-end cells,
where the first or last two interior unit gaps give the five-vertex case.

## 4. One missing gap and one unit gap

If these special gaps are separate, the missing gap yields a two-endpoint
cell. A backward edge between them is the rejected homogeneous-pair case;
otherwise its forward arc fixes both endpoint roles. The unit-gap cell
is then recovered as in section 3. An endpoint of the latter may match an
endpoint of the former, but the latter pair was fixed first by its own
arc, so the port recovery is not circular.

If the gaps are adjacent in the interior, there is a four-vertex cell of
one of the two forms

    X: e1, e2, u, e3;       Y: e1, u, e2, e3.

With no internal matching edge, all three endpoints have outside partners;
their complete external incidences and the local transitive order identify
the gap and the missing position. With one internal matching edge, its
third endpoint is a distinguished external port. There are exactly four
legal possibilities, because matching the consecutive endpoint pair at
the empty gap would create clones:

    X, pair(e1,e3): marked e2 has internal outdegree 2;
    X, pair(e2,e3): marked e1 has internal outdegree 3;
    Y, pair(e1,e2): marked e3 has internal outdegree 0;
    Y, pair(e1,e3): marked e2 has internal outdegree 1.

These marked degrees distinguish the four records. Importantly, the two
extreme-edge-reversal cells without their marked port are isomorphic;
erasing the external matching port would manufacture a false collision.
The extreme marked-degree cases leave only an internal C3 rotation. In
the other two cases the fixed marked vertex distinguishes the remaining
roles in the four-vertex tournament.

At an exterior end, the overlap cell has three vertices. On the left its
forms are e1,u,e2 (empty exterior gap, then unit gap) and u,e1,e2 (unit
exterior gap, then empty gap); right-end forms are their reverses. If the
endpoints match internally, the first is the allowed C3 and the second
is the rejected clone case. Otherwise their outside partners distinguish
the endpoint roles, and the local transitive order distinguishes the two
forms. The strong case has no exterior gaps and is already covered by
the four-vertex analysis at either end.

This exhausts the missing/unit interactions. The full matching, missing-
gap matching and P3 families have distinct numbers of vertices, so they
cannot collide with each other below delta=2k. Pair freeness and the proper-
cut strongness criterion are exactly the ones already proved. Thus every
record is recovered, up to the explicitly harmless rotations, proving (1)
through delta=2k-1.

## 5. The actual obstruction at delta=2k

Here two mechanisms acquire the same vertex count: a matching with two
empty gaps, and a P3-plus-matching with one unit gap. Their coincidence is
not an artefact of dropping ports. The following family retains an actual
external port and has precisely three good records in this grammar.

For every k>=2, use 2k ordered matching endpoints, with edges

    (0,2), (1,2k-1), (3,4), (5,6), ..., (2k-3,2k-2).

Give all gaps weight k except gaps 1 and 2, which are empty. The first
three endpoints form C3. Its endpoint 1 is distinguished by its external
matching partner 2k-1. All other local vertices see the outside identically.
The long edge (1,2k-1), together with (0,2), crosses every endpoint cut.
Thus omitting both exterior gaps gives a strong version as well.

The same weighted class has two P3 descriptions on 2k-1 endpoints:

    first:  (0,1), (0,2k-2), (2,3), (4,5), ...;
    second: (0,1), (1,2k-2), (2,3), (4,5), ... .

Every gap has weight k except gap 1, which has weight one. The internal
C3 now consists of the two first endpoints and their unit gap. Its shared
P3 endpoint is the same distinguished external port as before. In the
first description it is the first local endpoint; in the second it is
the last. Rotating the ordered presentation of the pointed C3 gives the
explicit isomorphisms. Actual vertices and external arcs are preserved.

All three records are pair free, have displayed feedback cost k, and have
mass M-2k (or S-2k after omitting exterior gaps). Their optimal feedback
cost is k by the same deficit identity. This is not a clone normalization.

Why the fibre is exactly three: large-weight vertices fix all positions
outside this one pointed C3. In the two-empty-gap matching record the
external port must be its middle endpoint, leaving exactly one cyclic
order. In a P3 record the port must be one of its two endpoint positions,
leaving two choices; the other roles are then forced. Full-gap matching
and one-missing-gap matching records have different vertex numbers. A
P3 record with no unit gap has one more large vertex and cannot represent
this class. The local pointed C3 is thus one class with one matching and
two P3 records, rather than three classes.

This recovers the mechanism of the older heavy-deficiency-two owner rule
in `defect_transitions_20261020/strata.md`, but at its distinct total-mass
threshold. The next task is a complete delta=2k correction including
multiple omissions and these owner fibres, not another naive continuation
of (1).

## 6. Exact scope of the companion checker

The checker exhausts the enlarged layers (k,delta)=(3,4),(3,5),(4,6),
both general and strong, canonicalizing every unit-cell permutation and
retaining clone rejection. It also tests the three-record collision at
k=2,3,4,5, with both exterior conventions. Larger all-k statements rest
on the preceding cell proof, not on an enlarged finite census.
