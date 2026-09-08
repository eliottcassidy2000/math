# Independent audit: finite-direction repair as matching with retained concurrency

**Verdict: ACCEPTED after two scope clarifications in the primary prose.**
The [matching completion](continuing12_20260908_no3_matching_completion.md)
proves exact unbounded combinatorial reductions, with separately declared
finite geometric controls. This audit imports no producer. It makes no
claim that the matching identities themselves are new graph theory.

## 1. Original cells, arbitrary blocks and the exact graph correspondence

For a finite board B and finite directions S, let V be its originally
overfull lines. A cell with no incidence in V is always safe to retain.
If every other cell meets at most two such lines, send two-incidence
cells to graph edges on V and one-incidence cells to edges with distinct
private leaves. Vertex capacities are two on V and one on each private
leaf, and each edge has capacity one.

This is an exact correspondence for feasible retained sets CONTAINING all
zero-incidence cells; maximum retained sets may always be chosen this way.
An arbitrary feasible set can omit a safe cell, so the unrestricted
bijection wording in the candidate was too broad. The optimum formula
is unchanged. Every original overfull-line capacity appears literally,
and originally safe lines cannot become overfull in a retained subset.
This proves the exact capacity-two graph b-matching formula for all
original line sizes under the at-most-two-incidences hypothesis.

The safe-set block inequalities are also exact lower bounds. Undoing
deletions on a cell that touches no originally overfull old line cannot
spoil the old repair. These undoable deletions form a separate budget for
new directions on that restricted original set. Iterating uses nested
ORIGINAL safe sets; the successive strata are disjoint. Recomputing a
safe set after a chosen repair would not justify the same independent
charge. For a selected old dual-line collection Q, the separate regions
N(Q) and B minus N(Q) similarly pay |N(Q)|-2|Q| and the new restricted
repair cost. None of these bounds presumes probabilistic independence or
claims that the optimum is independent of block order.

## 2. All-triple matching, Boolean dual and fractional completion

When each original overfull line is a triple, a deletion must cover every
line vertex at least once. Shared cells cover an edge's two endpoints;
private cells cover one endpoint. A maximum graph matching covers2nu
vertices with nu cells; each unmatched vertex requires at most one further
cell, yielding v-nu. Conversely an inclusion-minimal cell cover consists
of private singletons and shared-edge stars. Choosing one edge per
nontrivial star gives a matching of size v minus the cell-cover cost.
This proves tau=v-nu even when the graph is not bipartite.

For selected line vertices Q, the original union has3|Q|-|E(Q)| cells,
since there is no three-way incidence. Thus the Boolean value is
|Q|-|E(Q)|, whose maximum is the graph independence number alpha. Removing
at most one endpoint for every induced edge proves the upper bound;
independent sets attain it. The original cell labels pay this identity.

In the fractional deletion LP the upper bounds x_p<=1 may be dropped
without changing the optimum in this TRIPLE regime: capping any larger
variable at1 preserves each line's required deletion sum>=1. This
justifies the simpler displayed LP dual with vertex variables y>=0,
y_a+y_b<=1 on a shared edge, and y_a<=1 at a private cell. A nonisolated
vertex gets that same upper bound from an incident edge and nonnegativity;
isolated vertices have private cells. Setting z=1-y yields the fractional
vertex-cover dual, and finite LP duality gives tau_LP=v-nu_f. The
candidate omitted the upper-bound justification; it is now explicit.

The bipartite double-cover compiler is sound: copying each fractional
edge value in both orientations gives a feasible flow of twice the
matching weight; symmetrizing an integral double-cover matching gives
a fractional matching of half its size. Integral bipartite flow makes
the bounds equal. This restores the precise LP value, while ordinary
nonbipartite matching restores the integer value. A third direction does
not rule out every graph method; it changes which graph theorem applies.

## 3. Genuine hyperedges, higher demands and the complete branching reduction

For an arbitrary finite direction system let H be cells meeting at least
three originally overfull lines, and C those meeting one or two. Fix the
entire retained subset R of H; reject it if it already puts three cells
on a line. Every remaining line has residual capacity2-|R intersect line|.
The cells in C form the exact residual graph b-matching described above.
Free cells may again all be retained. Every feasible retained set fixes
one valid R and one graph solution, and conversely. Maximizing proves

    tau=|H|+|C|-max_valid_R (|R|+mu_b(R)).

There are at most2^|H| integral graph calls. This is an exact reduction,
not a bound on |H|, a polynomial-time assertion for all boards, or a
substitute of fractional graph solutions for the required integer ones.
One concurrent cell really can repair three lines at once; branching
retains that action. A pairwise clique would invent three different
actions and destroy the deletion predicate. Four-cell lines require
two deletions, which the capacity-two retained formulation handles but
an ordinary one-cover-per-vertex model does not.

## 4. Actual witnesses and independent finite paths

The n=6 two-regular board (03,12,04,14,35,25) has exactly the five overfull
lines stated in the primary, forming C5 with distinct shared cells.
The exact values are Boolean2, fractional5/2, integral3. Among the
two-regular square boards in the declared all-triple incidence regime,
no smaller board has a five-cycle COMPONENT. The entire n=6 universe
was not enumerated and no minimality for arbitrary finite cell sets is
claimed. Separated translates give independent copies on finite integer
boards, but their bounding squares contain empty rows and columns.

The n=5 board (24,01,23,34,01) has one actual concurrent three-way cell,
and needs one deletion. The pairwise-clique prediction of two is false.
The board (01,02,34,13,24) has an attached triangle with a perfect
matching, so no fractional gap occurs. The attachment must travel with
an odd-cycle charge. The all-triple theorem and the exceptional-cell
completion retain both boundaries explicitly.

The independent executable enumerates the full row-pair product and
filters column degrees, covering every2137 two-regular board n=2,...,5.
Direct increasing-cardinality deletion search is compared with a
separate residual-capacity branch optimizer on all of them. Exactly1320
enter the triple regime; exactly16 n=5 boards have exceptional cells.
For the eligible boards it independently computes ordinary matchings,
Boolean unions and independence. Fractional matching is solved through
half-integer primal and vertex-cover dual enumeration; equality of the
two exact feasible values certifies optimality, without assuming a
numerical LP output or a finite-to-infinite integrality extrapolation.

The same independent primal/dual comparison covers all844 labelled
subcubic graphs on at most five vertices. A literal cell-cover dynamic
program checks the integral formula independently of maximum matching.
These abstract incidence systems are not called geometric realizations.
Five actual boards also pay all24 singleton direction-block orders and
the full fourth-direction exceptional-cell branch. The complete engine
passes12767 always-active gates, normal and optimized raw LF outputs
agree, and the certificate regenerates unchanged. The unbounded
statements use the analytic arguments above and in the primary.
