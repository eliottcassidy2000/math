# Independent referee: unmerged incidence blocks and exact line-capacity states

**Status: PASS — independent analytic proof audit and FINITE-EXACT controls.**
The original-cell recurrence, arbitrary three-entry degree-cost compiler,
ordinary-matching call bound and connected geometric separation are accepted.
The finite referee imports no producer. It replaces blossom by exhaustive
ordinary matching and tests the block states using independent biconnected
components and literal local retained-subset enumeration.

## 1. Frozen targets and retained scope

The target is `continuing14_20260908_no3_line_states`, frozen before the
parent's status-only acceptance update:

    report 13ecfaf619853a4c049a1d8011ae00e3f424027a2c78257e91b5d4eac810e1af
    source 6fb957fc9be403b6242dac995e976184f9d52e3f103a732524b160521898d1ed
    output 84e5c9e9bea94e2e3b91aed7cc0d36cbebd556953b3bf56a1e9444d1472f2f1c
    cert   3f4aac7e677fc3fea5b1eca06f69b14d8bcf6ecb4d5c3886e587c69565c538c6

The inherited source is the literal originally overfull-line/original-cell
incidence graph, as in
[continuing13 cell separators](continuing13_20260908_no3_cell_separators.md)
and [its separator audit](continuing13_20260908_no3_separators_audit.md).
The earlier matching consumer is
`continuing12_20260908_no3_matching_completion.md`. Originally nonoverfull
lines remain nonoverfull under deletion, and cells outside the overfull-line
incidence can all be retained in an optimum. Thus this is exactly the
original-cell-loss problem for the chosen finite direction set.

The old construction merged through articulation lines. The new map keeps
ordinary biconnected blocks separate while retaining a line's exact local
contribution count. Binary cell state and three-valued line state carry
different information. Neither is replaced by a scalar repair cost.
The ordinary triangle still requires an integral graph consumer; the
concurrent triple and theta retain genuine shared decisions. None of the
claims concerns random-board frequency, an extremal no-three-in-line
constant, or polynomial time for arbitrary inseparable high-incidence cores.

## 2. Independent proof of both separator recurrences

Every incidence edge belongs to one ordinary block; bridges count as
two-vertex blocks. Distinct blocks can meet only at articulation vertices,
and their block-cut incidence graph is a forest. Root each component at a
block. A nonarticulation cell is charged in its unique block. An articulation
cell has one binary retained state shared by every incident block and incurs
its deletion cost once, either at the articulation node as in the producer
or in its nearest-to-root owning block as in this referee.

For a parent cell state s, each child block excludes the cost of that parent
cell. Its contribution is conditional on the same s. Hence the full cell
message is 1-s plus the sum of its child-block costs. This pays one original
cell once, even if it has arbitrarily many incidences or block occurrences.

For an articulation line, let t be the number of retained incidences in
the parent block. A child-block message specifies an exact contribution
r in {0,1,2}, not an upper bound. The child choices are independent once
their counts are fixed, and their sum must be at most 2-t. A min-plus
convolution with states 0,1,2 therefore computes the exact line message.
An impossible block count has infinite cost. Every line continuation has
finite cost for each parent t: deleting every child cell realizes zero
child contribution.

Inside a block, unary cell costs include all descendant cell continuations;
line reward profiles include all descendant line continuations. A parent
cell is fixed and excluded from the local deletion charge. A parent line
requires its exact prescribed contribution. Every nonparent line has its
ordinary capacity two. These statements account for all edges, cells and
constraints exactly once. Compatible block states therefore glue to one
feasible retained original set, and any feasible retained set restricts to
compatible states with the same cost. Forest induction proves the minimum
deletion recurrence in both directions.

The referee's implementation uses a different cost convention: it charges
a shared cell in the nearest-to-root block, leaving its separator message
as a sum without a unary term. It enumerates each block's complete retained
cell subsets directly, without a graph compiler. All tested roots agree
with whole original-subset optimization. Thus this is an independent test
of the charging convention as well as of the count states.

## 3. Arbitrary line rewards have an exact ordinary-matching gadget

After fixing all local high-incidence cells and any parent cell, remaining
decision cells have one or two incidences. Reject fixed choices that already
exceed capacity. If a remaining decision has costs c(0),c(1), start with
c(0) and assign retention reward c(0)-c(1). Shift each line profile by the
already retained fixed-cell count. This leaves arbitrary finite integer
reward profiles R_k for residual degrees k up to a capacity at most two.

For residual capacity two put

    L=min(R1-R2,R0-R1),
    a=R1-R2-L, b=R0-R2-2L.

Then a>=0 and b-a=R0-R1-L>=0. Join two line slots to one shared private
leaf with weight a on each edge, and join the slots to each other with
weight b. With zero, one or two externally occupied slots, the maximum
available internal reward is respectively b,a,0. The shared leaf prevents
two leaf edges being used simultaneously. Add constant R2+2L and charge
-L on each externally occupied slot. The resulting rewards are exactly
R0,R1,R2. No convexity, linearity or fractional relaxation is used.

The one-slot case uses constant R0 and incidence charge R1-R0; the
zero-slot case contributes only R0. These also preserve every profile
entry, including negative rewards.

For each original decision cell use two private ports joined internally.
For a two-incidence cell its two ports connect to the respective line
slots. For a one-incidence cell one port connects to its line and the
other to a private free vertex. Its retention reward is charged exactly
once, on the first port's external edges. Requiring both ports to be
covered leaves precisely two possibilities: the internal edge deletes
the cell, or two external edges retain it with all its incidences. A
half-retained original cell is excluded.

To enforce coverage in an ordinary unconstrained weighted matching, give
each covered port an additional reward M. An internal port edge gets 2M;
an external one gets M. If S is the sum of absolute noncoverage edge
weights, the difference between any two base matching weights is at most
S. Since all ports can be covered by internal edges, M>S forces every
maximum matching to cover every port. This works with arbitrary signed
weights; it is a dominance argument, not a maximum-cardinality convention.

If the parent line requires exact degree r, truncate its capacity to r,
then give each selected incidence there a bonus N exceeding the total
absolute base weight. Among fully covered matchings the degree is
maximized first. If that maximum is below r, the branch is infeasible;
otherwise it attains exactly r and then maximizes the true objective.
Only after adding N choose M larger than all resulting noncoverage
weights. This ordering makes port coverage the first priority, prescribed
degree the second, and the original objective the third. Subtracting the
artificial rewards recovers the original objective exactly.

The integer weights have polynomial bit length in the finite input data.
The compiler requires an integral ordinary weighted-matching optimum.
The independent finite oracle enumerates all ordinary matchings by vertex
mask recursion, so it does not depend on the producer's blossom call.

Two explicit boundaries are retained. Without mandatory coverage, a
two-incidence cell of reward ten can use the one available line slot while
its other line has capacity zero; the false half-cell gets ten although the
true retention optimum is zero. The coverage bonus repairs it. Separately,
for R=(0,1,0), the formula gives a=b=2 and L=-1. Replacing the shared
leaf by two separate leaves gives free-two reward four instead of two,
changing R0. The shared leaf is mathematically essential.

## 4. The exact branching parameter and complexity scope

Let h_B count cells with at least three incidences inside one ordinary
block B. Enumerating those cells leaves a one/two-incidence graph kernel.
A parent cell is already fixed and need not be branched. Each parent cell
has two possible states; each parent line has three; a root block has one.
Thus the number of ordinary matching calls is at most

    3 * sum_B 2^(h_B).

All finite-state convolutions have polynomial overhead, and each enumerated
high-cell assignment requires polynomial additional work. Equivalently,
total nonoracle work is polynomial per branch plus the polynomial separator
convolutions, in addition to the displayed exponential branching. This
does not claim that total work is polynomial in the incidence size alone
when h_B is unbounded.
With a polynomial-time ordinary weighted-matching consumer, the proved
total runtime is a polynomial factor in the finite input size times
sum_B 2^(h_B). The parent may clarify the corresponding sentence in the
filed primary during status promotion, recording that prose-only byte
transition; the preacceptance primary pins above remain the audited target.

In a nonbridge ordinary block every vertex has degree at least two. If
beta_B=|E_B|-|V_B|+1, summing degree minus two gives 2*beta_B-2. Every
counted high-incidence cell consumes at least one unit of this sum, so
h_B<=2*beta_B-2. A bridge has h_B=0. A uniform bound beta_B<=b on each
nonbridge block therefore gives the stated bound
3*number_of_blocks*2^max(0,2*b-2), regardless of the total cycle rank.
This is a sufficient condition. Large cycle rank carried by degree-two
cells does not itself create a hyperedge branching obligation.

No hidden line merge remains in this argument. High-incidence counts are
local to ordinary blocks, while articulation lines communicate their
capacity through exact counts. Genuine high-incidence cells inside a
single block still require the stated joint enumeration.

## 5. The actual geometric family and its exact optimum

An independent reconstruction of the eleven-point theta gives these six
overfull triples, indexed in the producer's displayed point order:

    slope1, label6:       {0,2,8}
    slope-1, label18:     {0,3,5}
    slope2, label0:       {0,4,10}
    slope1, label0:       {1,4,7}
    slope-1, label24:     {1,2,6}
    slope2, label-12:     {1,3,9}.

The exceptional cells are precisely 0 and 1. Three paths join them through
cells 2,3,4 and their line vertices, giving one theta block of cycle rank
two and h=2. The six private cells give six bridge blocks with h=0. The
common slope-one line has cells {1,4,7} and contains exactly one exceptional
centre.

Translate each copy by (100,100) and scale its coordinates by 1000^j.
The common slope-one label remains zero. Every other slope-one label is
a bounded nonzero integer times 1000^j; slope-minus-one labels lie in
[200,256]*1000^j; slope-two labels lie in a fixed negative interval
inside [-200,-1]*1000^j. These ranges do not overlap at distinct scales.
The occupied row and column ranges likewise separate. Consequently no
new cross-copy line coincidence occurs except on the declared common
line. There are 11k cells, 5k+1 overfull lines and 7k ordinary blocks,
with k values h=2 and 6k values h=0. The new bound is exactly 30k.

All these blocks merge into one bag under the inherited line-articulation
contraction; there is no separating articulation cell. Its 2k exceptional
cells have at most two retained common-line centres, while the other k
centres can be retained independently. Every other overfull line contains
at most one such centre. Thus the inherited capacity-filtered branch count
is exactly 2^k*(1+k+binom(k,2)). This is an exponential separation in k
between the stated oracle-call bounds, with actual integer coordinates.

For one copy let r be its retained contribution to the common line. The
exact minimum deletion costs for r=0,1,2 are 4,3,2. Upper bounds follow
by deleting the two exceptional cells and, as needed, one or two further
common-line cells. For r=0, the three compulsory common-line deletions
do not clear both remaining lines through the other centre, so a fourth
deletion is necessary. For r=1, two common-line deletions cannot clear all
three lines through that other centre, requiring a third deletion. For
r=2, one common-line deletion cannot clear the other five lines, requiring
two deletions. The referee also exhausts all 2,048 retained subsets of
this one literal copy to verify the complete table.

Summing the per-copy cost 4-r with total retained common-line contribution
at most two gives tau=4k-2, attained by allocating two retained common-line
cells. This proves the all-k formula without extrapolating the finite
family checks. Adding the isolated-copy scalar optimum two would instead
give 2k and lose the original common-line capacity.

These boards have at most two points on each occupied row and column;
their bounding squares contain empty rows and columns. They are sparse
finite geometric boards, not saturated squares or a random frequency model.

## 6. The nonconvex continuation and finite audit universe

For one five-point branch of the ten-point geometric control, its two
private overfull triples have point indices {0,1,2} and {2,3,4}, and its
two common-line cells are {0,3}. Literal enumeration gives exact
contribution costs (2,2,1) for r=0,1,2. Allowing at most 2-t retained
child contributions gives C(t)=(1,2,2). The first difference falls from
one to zero, so convexity is false on an actual geometric continuation.
The complete ten-point board has original repair cost three. This pays
the need for the arbitrary-profile compiler with a physical example.

The standalone referee checks:

* All 27 reward triples in {-7,2,11}^3, two original-cell patterns and
  four exact-degree modes: 216 comparisons of literal ordinary matching
  against exhaustive retained-subset optimization. It also checks an
  impossible prescribed degree, the half-cell hostile, and the shared-leaf
  hostile.
* All 97 simple two-regular boards in dimensions 2,3,4, generated by a
  full Cartesian row-choice enumeration and exact column-degree filtering.
  Their original retained-subset optima agree with independent block DP.
* All 109 linear triple-incidence systems on zero through four lines,
  with counts 1,1,2,9,96. Every root choice is checked. They are typed
  abstract incidences, not asserted to be geometric boards.
* The actual theta family for k=1 through 6, including line incidences,
  all ordinary block parameters, the exact old branch count and tau=4k-2.
  The complete one-copy cost table and the actual nonconvex ten-point
  continuation are separately enumerated.

NetworkX is used here for biconnected-component decomposition, a different
path from the producer's hand-written Tarjan decomposition. It is not used
as the matching oracle. All mathematical gates raise explicit exceptions.
The all-board and all-k claims rest on Sections 2--5; the finite controls
are an independent check of their implementation predicates and hostiles.

The actual adaptive-loss consequence is unchanged: a final feasible set
retains a feasible subset of the original cells, so its original-cell loss
is at least the original repair optimum. It does not turn this sparse
family into an extremal-density or random-board theorem.

## 7. Frozen reproduction

    python 04-computation/continuing14_20260908_no3_line_states_audit.py
    python -O 04-computation/continuing14_20260908_no3_line_states_audit.py

The source writes the same-stem certificate beside itself outside the
repository, or in `05-knowledge/results` after relocation. It passes
4,458 always-active exact gates in both normal and optimized Python.
Raw LF stdout and regenerated certificate bytes are identical.

    source 648cf2ec1671a85ff896adf3c84601bcceb8814d2d3299dd38baeb94a6083061
    output 4cdc75257b1770a74cba655b4457b4df2a1eac36fddddca308d38128c895324b
    cert   90f01990c256b185531c04ca87c4549a0c836991ba60c2ea44871a7968622041

The target is accepted with the explicit runtime wording in Section 4.
No mathematical repair to its recurrences, gadget, call bound or geometric
family is requested.
