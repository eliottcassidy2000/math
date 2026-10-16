# Component flips and the six colour palette matching test

2026-10-16. **PROVED:** the component palette extension lemma and the
complete nine-blocker criterion for its unrestricted palette banks.
**FINITE-EXACT:** all 32,768 pair banks, matching/cover witnesses, all 570
geometric boundary patterns, and an abstract four-chromatic control.
The actual nine-point incidence bank permits at most one blocked pair,
so two admissible pairs suffice even when they overlap.
**OPEN:** an authenticated unit-distance
forcing context in `Q(zeta55)`, Hadwiger--Nelson, uniform Collatz payment,
and planar Jacobian. No new geometric embedding is asserted.
Self-audited; no independent referee is claimed.

## 1. Recover the stronger boundary mechanism before searching again

Keep exactly the eight pentagon/rhombus ports of
[THM-4614, three-colourable-contexts-cannot-change-the-six-colour-port-relation-in-q-zeta55](../../01-canon/theorems/THM-4614-three-colourable-contexts-cannot-change-the-six-colour-port-relation-in-q-zeta55.md).
Every exterior point has at most two port neighbours, and the whole field
has just nine possible exterior double-port points. A three-colourable
exterior preserves all 570 boundary patterns.

The [interface profile lemma](six_colour_interface_profile_20261015f.md)
extends this to certain four- and five-colourable exteriors. The
[safe singleton lemma H2](compression_contexts_20261015g.md)
handles a four-colouring with one port-free class and a globally missing
boundary colour on another. Those results are sufficient extension tests,
not sufficient conditions for forcing when they fail.

The new sidecar is a connected component of the subgraph induced by **two**
exterior colour classes. Its two-colouring can be swapped independently
of other components. Replacing all components by one global class-level
palette loses exactly those independent choices.

## 2. Component lists and a bipartite palette graph

Let H be the exterior, with a supplied proper four-colouring
`C1,C2,C3,C4`. Fix a proper six-colouring beta of the marked ports.
For each connected component K of `H[C3 union C4]`, let

\[
 S_K=\{\beta(p):p\text{ has a neighbour in }K\cap C3\},\qquad
 T_K=\{\beta(p):p\text{ has a neighbour in }K\cap C4\}.
\]

An unordered colour pair `{a,b}` is **admissible** if, for every K,

\[
 (a\notin S_K\ \&\ b\notin T_K)
 \quad\text{or}\quad
 (b\notin S_K\ \&\ a\notin T_K).                    \tag{1}
\]

Call the set of such pairs Q. The orientation may depend on K. This
criterion is exact for using the supplied bipartition with just the two
colours `{a,b}` in each component, including isolated vertices.

Build a bipartite graph `P_Q`. Both shores are labelled by the fifteen
unordered pairs from the six colours. Join the left pair P1 to the right
pair P2 precisely when they are disjoint and

\[
             [6]\setminus(P1\cup P2)\in Q.            \tag{2}
\]

An edge is a candidate palette allocation: P1 for C1, P2 for C2, and
the complementary pair for the two-class components. A vertex in C1
with two differently coloured port neighbours forbids just one left
pair; a vertex in C2 similarly forbids just one right pair. Vertices
with zero or one prescribed neighbour colour never exhaust a pair.

**H3 (PROVED).** Delete the actual forbidden left/right palette vertices.
If one edge of `P_Q` remains, beta extends to the entire graph. In
particular, if `d12` counts the distinct forbidden class/pair types in
C1 and C2 and `nu(P_Q)` is its matching number, then

\[
                     d12<\nu(P_Q)                    \tag{3}
\]

is sufficient.

**Proof.** Use the surviving edge's palettes. In C1 and C2 choose an
allowed colour at each vertex. In each K choose one of the orientations
certified by (1). All edges between different palette blocks are proper
because the palettes are disjoint; all edges inside K are proper by its
bipartition; different components have no edge between them. Port edges
are proper by construction. A matching consists of options with distinct
left and right pairs, so each forbidden vertex destroys at most one of
them. Fewer than `nu(P_Q)` forbidden vertices cannot destroy every option.
QED.

This is a literal lift of the prescribed port colouring. Merely counting
nonempty local lists would not suffice: in a path `u-w-v` with u,v in C3,
if u forbids a and v forbids b, every vertex has a nonempty `{a,b}` list
but the whole component cannot use that pair. Conversely, two separate
components may require opposite orientations without any conflict.

## 3. Exactly which pair banks withstand nine arbitrary blockers

**H4 (PROVED).** A pair bank Q survives **every** set of at most nine
class/pair blockers in the palette construction if and only if either

1. Q contains two disjoint pairs; or
2. Q contains all five pairs incident with one colour.

The criterion concerns unrestricted syntactic blocker banks. A bank that
fails it need not be realizable by the nine actual geometric double-port
points, and failure of this palette construction is not failure of a
six-colour extension by other means.

**Proof of sufficiency.** Relabel two disjoint pairs as `01,23`. The
following ten edges form a matching in `P_Q`:

| left | right |
|---|---|
| 01 | 45 |
| 04 | 15 |
| 05 | 14 |
| 14 | 05 |
| 15 | 04 |
| 24 | 35 |
| 25 | 34 |
| 34 | 25 |
| 35 | 24 |
| 45 | 01 |

Every complement is 01 or 23; no left or right pair repeats. Nine
blockers cannot cover these ten independent options.

For the full five-edge star, remove its central colour. The ten pairs
on the remaining five colours form each shore of a three-regular
bipartite graph, where disjoint pairs are adjacent. For any left set A,
counting its incident edges gives `3|A|<=3|N(A)|`. The usual augmenting-path
argument therefore supplies a perfect matching of size ten. Equivalently,
the output file includes an explicit ten-edge matching. Again nine
blockers cannot cover it.

**Proof of necessity within the palette construction.** If Q contains
no disjoint pairs, its edges are pairwise intersecting. Such a family is
contained in a star or a triangle: after fixing edges `ab,ac`, any edge
not containing a must be bc; in that case every edge lies in that triangle.
If the full five-edge star is absent, Q is contained in a four-edge star
or a triangle.

For the star `01,02,03,04`, forbidding the four pairs `15,25,35,45` on
each shore covers every edge of `P_Q`, using eight blockers. Indeed a
complement in this star leaves colour5 in exactly one of P1,P2, and
colour0 in neither. For the triangle `01,02,12`, forbidding `34,35,45`
on each shore covers every edge, using six blockers: the four remaining
colours contain 3,4,5, so one of the two disjoint pairs lies entirely in
that triple. Subfamilies inherit the covers. This proves necessity. QED.

Every bank with at least five pairs passes H4. The exhaustive control
finds 32,582 safe banks and 186 unsafe banks among all `2^15` possibilities.
These counts are for an abstract finite palette space, not a distribution
of geometric contexts. Matching sizes and equal-sized vertex-cover
certificates are checked for every bank; no integer-programming oracle
or floating-point acceptance is needed.

## 4. The stronger screen at the unchanged geometric ports

For the inherited `Q(zeta55)` ports, `d12<=9` automatically. Consequently,
**two disjoint component-admissible pairs are enough to extend beta**.
This can work when all four exterior classes meet the ports and when
each class, in aggregate, sees all six boundary colours.

The old H2 is recovered exactly: if C4 is port-free and colour c is
missing globally from C3's port neighbours, then all five pairs `{c,a}`
are component-admissible. The four-class/two-interface case has every
pair admissible. The three-colourable theorem is also recovered by
taking C4 empty: each C3 component is a single vertex, and only a
double-port blocker can exclude its one forbidden pair.

For a prescribed boundary pattern to be removed by a four-colourable
exterior, **every** supplied proper four-colouring and every choice of
two classes to merge must fail the extension test. In particular, its Q
must be contained in a four-edge star or a triangle; the actual remaining
blockers must additionally cover the surviving palette graph. This is a
necessary search screen. It supplies neither a four-chromatic unit graph
in the field nor a forcing gadget.

The concrete control is the abstract Mycielskian of C5, with six isolated
single-port neighbours added to each of its four colour classes. The
eleven-vertex core is checked to be triangle-free and not three-colourable.
In the supplied four-colouring every class sees all six boundary colours,
but each two-class component sees at most a single colour on its constrained
side, so all fifteen pairs are admissible. The script constructs and checks
a complete six-colour lift at the eight-port graph. **No unit-distance
realization of that exterior in the field is asserted.**

## 5. The actual nine-point bank is stronger than a count of nine

The inherited double-port pairs, in the fixed port order, are

    02, 03, 13, 14, 15, 17, 24, 45, 47.

For a fixed prescribed boundary pattern beta, let `c_P` be the number
of these nine points whose two port colours are exactly the distinct pair P.
For a candidate component pair q, its four-colour complement has three
partitions `P union P'`. A **necessary and sufficient** condition for
some assignment of the nine points to C1/C2 to block every palette option
with complementary pair q is

\[
                 c_P+c_{P'}\ge2
       \quad\hbox{for each of the three partitions}.   \tag{4}
\]

This is an unrestricted assignment test, not a supplied exterior colouring.
For necessity, the two oriented options `(P,P')` and `(P',P)` need two
different blockers. For sufficiency, two copies of one type can be placed
on opposite shores; if both types occur, put one copy of each on the
same shore. The three groups use different pair types, so these choices
do not compete for points.

**G1 (FINITE-EXACT complete boundary lemma).** For every one of the 570
proper six-colour patterns, **at most one** pair q satisfies (4). The
exact exhaustive distribution is:

| colours used on the ports | no exceptional pair | one exceptional pair |
|---:|---:|---:|
| 3 | 30 | 0 |
| 4 | 144 | 36 |
| 5 | 223 | 24 |
| 6 | 111 | 2 |

The capacity test (4) and a separate exhaustive assignment of the actual
nine point types agree for every pattern: 96,688 C1/C2 assignments of
the active distinct-colour blockers, and 54,392 distinct pattern/bank
combinations. Same-colour points cannot block a pair and are omitted from
that enumeration. Every exceptional pattern has a saved
blocking assignment. Proper patterns are enumerated by restricted-growth
words and reproduce the inherited `30/180/247/113` counts; all fifteen
candidate pairs, including unused palette colours, are retained.

This is a finite exhaustive lemma on the canonically specified port graph,
conditional only on the inherited exact nine-point geometric census. It
does not enumerate exterior graphs. Allowing every point to lie in C1
or C2 is an overestimate of the actual blocker bank: missing points and
points in C3/C4 can only remove blockers.

**Consequence of H3 and G1:** if Q has **any two distinct pairs**, the
boundary pattern extends. For 508 patterns, any nonempty Q suffices.
For the other 62, a nonempty Q can fail this screen only when it is the
singleton containing that pattern's explicitly listed exceptional pair.
The actual class assignment must still realize its blocking cover.
Failure of the screen remains weaker than geometric forcing.

The geometry explains part of this strengthening: one colour class in
the port graph is an independent set. There are 46 such sets, and none
is incident with more than six distinct-colour double-port points.
The arbitrary-blocker covers from H4 need not obey that incidence bound.
The complete capacity test retains still more information than this
one-colour statistic.

### All boundary patterns using at most four colours extend across a four-colourable exterior

**H5 (PROVED, elementary).** More generally suppose an exterior has a
proper four-colouring, each exterior vertex has at most two port neighbours,
at most eleven exterior vertices have two differently coloured port
neighbours, and the prescribed boundary uses at most four of six colours.
Then the boundary colouring extends.

**Proof.** Two colours a,b are unused on the boundary. Among the four
exterior classes choose C1,C2 with the fewest double-port blockers.
Their total is at most five. Assign a,b to the other two classes, which
is safe on every port edge. There are six ways to split the remaining
four colours into ordered pairs for C1,C2. Each of their blocking vertices
destroys at most one option, so an option survives. Choose available
colours inside its independent classes. QED.

For the inherited nine-point bank this preserves all **210 patterns using
three or four colours** across every four-colourable exterior, without a
port-free-class hypothesis. It complements the earlier H2: that result
preserves 457 patterns under a stronger interface assumption. No optimality
of the numerical bound eleven is asserted.

The boundary-colour condition cannot simply be discarded in arbitrary
graphs, even when there are no double-port vertices. Take a complete
four-partite exterior with five vertices per part; in each part one vertex
forbids each of five used boundary colours via a single port neighbour.
Only the sixth colour can colour a whole part alone. At most one part
can use it as a singleton, while the other three need two colours each;
seven colours would be needed. This exterior is four-colourable, but its
given five-colour boundary does not extend with six. It has triangles
and is **not** a unit-distance graph in `Q(zeta55)`; it is a hostile to
dropping the component/geometry hypotheses, not a field counterexample.

## 6. What transfers to Collatz and inverse clocks

This session's [refinement theorem](collatz_refinement_debt_20261016a.md)
has a specific shared-label analogue. Physical source averaging commutes
with a common native map before norms are taken. The positive Gram
remainder survives that square and has a uniform budget; independently
maximizing each source row loses compatibility and creates sign-mixture
debt. Here component orientations replace an overly rigid global palette,
but each component's full allowed relation is retained before gluing.

The useful research move is to pack **mutually noncompeting witnesses**
before counting obstructions. H4 obtains ten palette options with distinct
blocking coordinates; the geometric census supplies at most nine blockers.
For Collatz, a genuine counterpart would be disjoint opposite-sign source
matches in actual terminal fibres, with their original weights and literal
clocks. The missing weighted lower bound cannot be imported from the
unweighted number ten or the palette theorem.

For planar Jacobian, retain the already proved inverse-clock scope from
[the Keller note](keller_inverse_clock_20261015f.md) and the
[compression invoice](compression_contexts_20261015g.md): on the actual
cokernel the source operators satisfy `[E,t]=I`, while a finite projection
has a missing-excursion operator of trace equal to the retained dimension.
Small rank does not make that cost small. A factorially small inverse orbit
still owes the return operation `E^m` before it represents the same source.
The present Collatz refinement term contracts because its actual native
maps are signed one-target maps. No corresponding contraction of the
Keller return operators has been established.

| source | representation/map | preserved predicate | necessary sidecar |
|---|---|---|---|
| annular source rows | common native transport and physical averaging | the same linear forms and Gram remainder | original weights, shared terminal and integer clock |
| exterior four-colouring | two-class components to pair-admissibility graph | a lift of the fixed six-colour port pattern | per-component orientation and actual blocker vertices |
| Keller cokernel | finite projection of source operations | the commutator only with its excursion term | the omitted excursions and return operator |

These are operation-level comparisons, not equivalences of the three
open problems. The next geometric question is whether an actual field
context can make Q empty, or the one pattern-specific exceptional singleton,
for every choice of two classes in every proper four-colouring, while
realizing the required geometric blocker assignments. H5 restricts any
removed boundary pattern in this scope to five or six colours.
The next Collatz question remains
an independent lower bound on the weighted capacity of authenticated
opposite-sign source matches.

## Reproduction

    python3 04-computation/experiments/six_colour_component_pairs_20261016a.py
    python3 -O 04-computation/experiments/six_colour_component_pairs_20261016a.py

[Output](six_colour_component_pairs_20261016a.json) includes all histogram
counts, explicit matching/cover examples, the abstract graph's supplied
colouring, its literal lift, and the unchanged eight-port colouring.
The output also pins the inherited geometric-input hash and saves all 62
exceptional patterns with their class/pair blocker witnesses.
H5 also has 1,365 load-composition and 1,586 blocker-bank controls; its
five-colour abstract hostile rejects all 4,096 disjoint palette allocations.
SHA256: `686b12dec8e9c82cd5bcd544882c9fd258fe725bf00e9edcbc7ee3a46c37ae76`.
