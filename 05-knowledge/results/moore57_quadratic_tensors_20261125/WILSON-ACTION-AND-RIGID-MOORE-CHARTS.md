# Yang–Mills and QCD lens: Wilson collision action for Moore charts

2026-11-26. **PROVED** finite action and rigidity statements; **FINITE-EXACT**
small controls and marked-board optimization. No claim about a continuum
mass gap, particle spectrum, confinement, or completion of the degree-57
Moore graph.

## What comes from the linked group and from the repo

The [Yang–Mills group's index](https://github.com/yang-mills-group/submissions)
separates the received manuscripts, checked mathematical statements, and
the interpretation of those statements. Its [October pure Yang–Mills
assessment](https://github.com/yang-mills-group/YM) credits a Wightman
construction by Wilson-lattice scaling, not uniqueness or a mass gap.
The [25 November QCD assessment](https://github.com/yang-mills-group/QCD)
describes an SU(3) construction with six massive fundamental Dirac flavors
and Wilson fermion determinants. It distinguishes continuum construction
from identifying operators, the fitted masses, and the coupling with
conventionally normalized QCD. Its summary does not establish a mass gap,
confinement, or a particle spectrum. These are reports of the assessments,
not independent checks of the continuum proofs.

The repo's [simplicial Wilson work](../session34_gauge_representations_20261023/gauge/SIMPLICIAL_WILSON_FANO.md)
already separates a gauge-invariant loop from a Hamiltonian count on a
marked slice. The [spectral and flux work](../session272_yang_mills_bridges_20261111/spectral_loops/PROOFS.md)
shows why local observables need not recover long holonomies, and why
energy estimates are not particle identification. The [SU(3) color
work](../session77_ternary_color_20261027/COLOR_FUSION.md) supplies a different
literal connection to Petersen through a signed six-port Gram matrix.
Here we transport the Wilson method using **permutation-valued sheet
matchings**, not SU(3) representations or physical quarks.

## 1. A gauge-invariant action with exactly the Moore charts at zero

Root a putative diameter-two Moore graph of valency k. Its k branches
have d=k−1 leaves each, and every two leaf branches carry a perfect
matching: a leaf in branch i has a unique two-edge path to the head of
branch j, through a leaf in j. Conversely start from a full system of
such matchings, indexed by all unordered pairs i,j, with no assumption
of the Moore property. Let P_ij be the permutation matrix of the link,
P_ji=P_ij^T. Add the rooted tree if recovering the full graph.

For i≠j define the path-count matrix and the defect action by

    R_ij = P_ij + Σ_{l not i,j} P_il P_lj,
    E(P) = Σ_{i<j} ||R_ij - J_d||_F² .

Each row of R_ij counts d routes; its entries sum to d. The gauge is
independent relabeling of the d sheets in each branch, so R_ij and its
defect transform by the same row and column permutations. E(P) is
exactly invariant. The full chart is Moore iff E(P)=0: all cross-branch
leaf pairs are then reached exactly once by a direct link or a leaf
intermediate. Same-branch leaves have the common head and no additional
common neighbor by bijectivity; the remaining rooted-tree cases are
already unique.

For a fixed pair i,j, let m_xy be the number of routes from x to y.
Since Σ m_xy=d², the squared defect is 2Σ binom(m_xy,2). This is a
literal Wilson count: for two route permutations U,V, the number of
coincident routes is tr(U V^T), i.e. the fixed-sheet count of the
relative holonomy. Pairing a direct route with a length-two route
counts a triangle three times, and pairing distinct length-two routes
counts a square twice, once for each opposite pair. A square in the
leaf graph cannot have two opposite vertices in one branch because a
leaf has only one neighbor in that branch. Thus, for the leaf graph,

    E(P) = 6 * (# triangles) + 4 * (# squares).

These loop weights are derived from Moore route multiplicities, not
the usual Wilson lattice action. Its zero seeks **no fixed sheets
on the comparison loops**, not identity holonomy. This sign distinction
is essential. For a finite determinant parallel, write X for the
leaf adjacency with N=kd and d-regular rows; then this action is

    tr(X³) + (tr(X⁴) − N d(2d−1))/2 .

These traces are the cubic and quartic terms in the formal closed-walk
expansion of −log det(I−zX), with the quartic backtracking contribution
subtracted as above. This is an exact finite determinant identity,
not the fermion determinant or a claim about renormalization.

### Exact response to changing one link

From a zero-action chart, replace one permutation by another moving
m of its d sheets. The change Δ has squared Frobenius norm 2m. There
is one affected direct defect and 2(k−2) affected length-two defects;
each is Δ times unchanged permutation matrices (or transposes), so
has the same norm. Thus one changed link gives

    E = 2m(2k−3),       triangles = m,       squares = m(k−3).

Every added edge closes precisely one triangle: the old Moore graph
had a unique length-two route between that formerly nonadjacent pair,
through a third branch, and neither leg is changed. These are all the
new triangles since a triangle cannot contain two links of the same
branch pair. The number of squares then follows from the action.
A transposition at degree 57 would cost **444** action units (two
triangles and 108 squares), conditional on starting at an exact chart.
The verifier checks all 318 single-link transpositions across
Petersen and Hoffman–Singleton by two independent cycle counts.

## 2. Frozen branch stars force a collective change

Take two distinct zero-action charts on the same labelled sheets.
Say a branch is untouched if *all* its incident matchings agree.
Let there be m untouched branches. At any leaf x where some link
changes, **at least m+1** branch neighbors of x must change.

Proof: choose a changed neighbor branch j of x, and let u be x's old
neighbor in j. For each untouched branch k, follow u to its neighbor
w_k in k. In the new chart, the original route x-u-w_k is gone.
The direct route and other unchanged routes do not become a second
route to w_k, and the new neighbor in j cannot hit w_k by bijectivity.
So the unique new route must go through the new neighbor y_l in some
other changed branch l of x. Both u and y_l were adjacent to w_k in
the old chart (links to k are fixed). The same l can serve at most
one k: otherwise u and y_l would have two common neighbors in the old
Moore graph. Thus x needs at least m other changed neighbors in
addition to the one in j. This also gives a concrete additional
constraint when equality holds: all changed exits other than j must
be assigned to distinct untouched receiver branches.

In particular the graph of changed branch links has minimum degree
m+1 on its k−m nonisolated vertices. Thus

    m <= floor((k−2)/2),     changed links >= ceil((k−m)(m+1)/2).

For Petersen, one frozen branch star already determines the chart;
for Hoffman–Singleton, three do. In degree 57, **two distinct
completions of the same three fully fixed branch stars would have to
differ on at least 108 whole matchings**: minimize the bound over
3<=m<=27. Each differing permutation removes at least two original
leaf edges, so at least 216 edges must be removed and 216 added.
These statements compare charts in common labels. They do not
establish that there are two completions, or even one.

## 3. The gauge group cannot be flattened to an abelian sheet action

Set all transports out of one arbitrarily fixed branch to identity
using the gauge. Let H be the group generated by the resulting
off-star transports, equivalently the based loop holonomy group.
At zero action H acts transitively, since any route list in a chart
row visits every sheet once. For d>2, **the commutator subgroup H'
must also be transitive**.

If not, take its sheet orbits, each of size s, and the induced abelian
transitive group on those orbits. It is regular on that quotient.
Write its translations additively, a_ij for the projected transport,
and σ for the sum of the elements of the quotient group. Each
ordered route list on the orbits is a list with each group element
occurring s times, so its sum is sσ. For every branch i other than
the gauge root, the sum of nonroot outgoing link labels equals sσ
by the route list from i to the root. Also a_ji=−a_ij. Thus the sum
of route labels between nonroot branches i,j is −a_ij: the row sums
cancel, and the intermediate path through the root contributes zero.
It follows every nonroot link has the same projected label −sσ.
But each projected route list to the root has one identity and
d−1 copies of that one label, so a nontrivial uniform quotient is
possible only if d=2. This proves the claim without assuming an
automorphism of the Moore graph.

In particular a regular H would have to be perfect. For d=56 this
is impossible: the number of Sylow 7-subgroups is 1 or 8. In the
first case there is a nontrivial quotient of order eight, not
perfect. In the second case the 48 nonidentity elements of order
seven leave just eight elements for a Sylow 2-subgroup, making it
unique and giving a quotient of order seven. In particular **no
regular 56-sheet permutation group can carry all the transports**.
This extends the repo's [cyclic fixed-star obstruction](CYCLIC-FOUR.md)
to a statement about the full holonomy action, without imposing
equivariance by a cyclic sheet symmetry. The small controls have
holonomy Z/2 for Petersen and S6 for Hoffman–Singleton.

## 4. The actual marked board has a positive small-edit barrier

For the [physical marked board](HADWIGER-NELSON-PHYSICAL-SLOTS.md),
pick one individually legal contact in each remaining branch.
Define its local squared defect by Σ over physical center slots
of (occupancy−1)^2. With the marker fixed all group totals agree.
Half this local energy is the number of coincident paths, or pairs
of new edges making a physical square with the old graph. This
objective accepts conflicting contacts while scoring them, rather
than pretending that they are valid extensions.

Retain the 50-reference patch after the marker and fill the two
vacant branches. At zero old-contact edits the minimum is three
squares, achieved by the new contacts (13,16) and (33,19).
Exact bounded optimization finds the same minimum through **six
reference changes**, so not even a one-square or a two-square
full selection is that close. The search enumerates a subset of
changed reference branches, removes their old contacts, and inserts
replacement choices before filling the two holes. Adding a contact
cannot decrease the current collision count, so pruning at the
incumbent is sound. The independent replay inserts no invalid claim
of girth: it counts the exact three new squares by actual common
neighbors in the 3,250-vertex graph. Zero on the full domain remains
open; the edit barrier is not a continuum mass gap.

## Replay

From this directory:

```sh
python nested_core_check.py --board > board.txt
c++ -std=c++17 -O3 marked_wilson_min.cpp -o /tmp/marked-wilson-min
/tmp/marked-wilson-min 6 < board.txt > action.txt
python wilson_chart_check.py board.txt action.txt
```
