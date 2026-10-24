# An atlas-free score-mass closure of the regular-eleven correction

24 October 2026. Status: PROVED counting identities and input reduction;
FINITE-EXACT evaluation by two integer recurrences. The complete nonrigid
regular-eleven classification is an explicit source-pinned inherited input.
No almost-regular-ten or regular-eleven atlas is needed by this packet.

## 1. Results and scope

The exact labelled regular-eleven tournament count is

    L11=48,251,508,480.

The inherited classification gives eighteen nonrigid ordinary classes, with
full automorphism orders and orbit counts:

    three of Aut11, q=1; one of Aut55, q=1;
    seven of Aut5, q=3;
    one of Aut9, q=7; four of Aut3, q=7;
    two of Aut3, q=9.

Their labelled orbit mass is

    M=3*(11!/11)+(11!/55)+7*(11!/5)+(11!/9)+6*(11!/3)
     =151,764,480.

Every remaining regular class is rigid and has exactly11! distinct labelled
versions. Therefore

    rigid regular11 classes=(L11−M)/11!=1205;
    all regular11 classes=1205+18=1223;
    almost-regular10 classes=11*1205+78=13333;
    E_regular11=110*1205+660=133210.

The regular-eleven defect distribution, including error zero, is

    4 +7 t^22+5 t^66+2 t^88+1205 t^110.

The independently pinned nonregular/SCC constants then give

    E_s,11=133210+9166=142376,
    E_11=133210+19480=152690,
    a_11=(U11−152690)/11.

The unrestricted upper trace U11 and unrestricted order-eleven count a11 are
not computed. This is a labelled score-count and automorphism-mass sidecar,
not a claim that a bare metagraph or scalar sequence alone predicts a11.
The numerical regular11 and almost-regular10 class counts are outputs, not
inputs, calibration targets, or expected atlas lengths.

## 2. Fixed labelled score counts

For an integer vector d=(d_1,...,d_n), let N(d) be the number of tournaments
on a fixed labelled vertex set with those specified outdegrees. It is zero
unless each degree lies in0,...,n−1 and the sum is n(n−1)/2. The empty score
vector and the one-vertex zero vector each have count one.

Relabeling vertices gives a bijection between the counted sets for two
permutations of d. Hence N depends only on the multiset of score entries.
Sorting this vector is a legitimate memoization key; it does NOT turn the
objects being counted into ordinary isomorphism classes. No factorial is
divided out or multiplied in when sorting.

Distinguish one particular vertex v with required degree a, and let the
remaining required degree r occur m_r times. Choose x_r of these m_r vertices
as outneighbors of v. The constraints are0<=x_r<=m_r and sum_r x_r=a.
There are product_r binom(m_r,x_r) such choices of actual named neighbors.
A vertex beaten by v keeps required residual degree r within the remaining
n−1 vertices. A vertex beating v has residual degree r−1. Thus

    N(d)=sum_(sum x_r=a) [ product_r binom(m_r,x_r) ]
          *N( multiset union_r { r repeated x_r,
                                  r−1 repeated m_r−x_r } ).

This is exact: every counted tournament has one unique outgoing-neighbor set
at the distinguished vertex, and deleting that vertex leaves precisely one
residual tournament. Conversely every such neighbor choice and residual
tournament extend uniquely. Sorting only identifies equal-sized subproblems
by a relabeling bijection; the binomial factors retain their labelled
multiplicities. There is no extra factor for choosing v, which was already a
particular fixed labelled vertex.

The grouped implementation distinguishes the first, smallest-score vertex.
Its bounds additionally skip an incoming neighbor of required degree zero,
and an outgoing neighbor of required degree n−1, since these would create an
invalid residual degree. These are immediate impossibility filters, not
heuristic pruning or a future-multiple assumption.

## 3. Independent literal-subset recurrence

The second implementation distinguishes the last, largest-score vertex and
visits every actual subset X of its remaining named neighbors of size a.
It uses unit weight for each subset and residual degrees

    d'_j=d_j−1_(j not in X).

It does not group equal degrees or use binomial weights. Both recurrences use
exact Python integers and memoize only the degree multiset. They agree on
all regular cases through eleven and on the almost-regular ten profile.
Their agreement is a separate implementation check; the proof in Section2
is what establishes correctness and exhaustive labelled coverage.

The computed labelled regular counts for orders1,3,5,7,9,11 are

    1, 2, 24, 2640, 3230080, 48251508480.

The first four have the independently certified tiny regular controls as an
additional check. Neither implementation is supplied a numerical target for
order eleven or the eventual class count1223.

## 4. Independent fixed-root deletion check

Deleting one particular labelled vertex of a regular-eleven tournament leaves
a tournament on the other ten labels with exactly five degrees4 and five
 degrees5. Conversely every such labelled tournament has a unique regular
extension at that particular root: degree-four vertices beat the root, and
the root beats degree-five vertices.

There are binom(10,5)=252 possible choices of the labels with degree4, and
all profiles have the same count by relabeling. Both recurrences give

    N(4,4,4,4,4,5,5,5,5,5)=191474240,
    252*191474240=48251508480=L11.

This provides a distinct structural check on the score counting, without
constructing the almost-regular ordinary-class atlas or doing a reroot trace.
The factor252 counts disjoint labelled degree profiles, not ordinary classes.

## 5. Why subtracting the nonrigid mass is exhaustive

The pinned nonrigid theorem proves, by full automorphism orbit partitions,
that there are no other nonrigid regular-eleven classes. Its candidate groups
are full ordinary automorphism groups, not a displayed subgroup or merely a
stabilizer of a chosen equitable presentation. This completeness and full-group
information are load-bearing inputs to the present calculation.

For an ordinary class T, the number of distinct labelled tournaments on the
fixed eleven labels is11!/|Aut(T)|. Classes are disjoint and all represented
by labelled tournaments. Subtracting M therefore removes exactly every
nonrigid labelled object from L11. Every remaining object lies in a rigid
class of labelled orbit size11!, proving the exact quotient in Section1.
A missing nonrigid type or use of a proper subgroup would invalidate this
subtraction; no mere agreement with a target class count could repair it.

The eighteen known classes contribute78 pointed types, one for each full
vertex orbit, and total refinement error660. Each rigid regular class has
11 pointed types and error110. Those facts prove all remaining arithmetic
in Section1. Equivalently E_regular11=11(P10−a_regular11)=10P10−120, providing
an independent consistency identity from the pinned structural theorem.

## 6. Validation and limits

The only data files read are three source-pinned records: the nonrigid
classification receipt, its proof, and the independently audited nonregular
/SCC proof. The score recurrences have no graph atlas input. Every check uses
an explicit exception and survives `python3 -O`.

Normal and optimized outputs are byte-identical. Six physical-copy hostile
mutations test lost binomial multiplicity, reversed residual-degree ownership,
a missing252 profile factor, duplicate hidden Aut55 control, replacing group
orders by vertex-orbit counts in labelled mass, and confusing rigid pointed
orbits with one class. Each must fail semantically and against the manifest.
The final replay copies only manifest files into a fresh directory and repeats
both Python modes and the hostile suite from an unrelated working directory.

The input reduction is substantial but scoped: a complete nonrigid regular11
classification plus labelled score counting replaces an almost-regular10
atlas for this correction. No such nonrigid completeness theorem or compact
counting closure is asserted for arbitrary larger orders.
