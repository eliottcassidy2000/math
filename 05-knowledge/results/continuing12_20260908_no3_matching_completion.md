# Matching completes a finite-direction repair regime

**Status: PROVED ANALYTICALLY + FINITE-EXACT controls; INDEPENDENTLY AUDITED.** This extends the continuing11 three-direction repair mechanism to any finite family of directions, gives an exact nonbipartite matching regime, and identifies its first geometric failure boundary. It does not assert a new extremal no-three-in-line bound or a new random-board asymptotic coefficient.

## 1. Inheritance and the representation that changes the problem

The closest proved suppliers are [continuing10 repair](continuing10_20260907_no3_repair.md), which gives integral capacity-two bipartite flow, and [continuing11 three directions](continuing11_20260908_no3_three_directions.md), which gives the safe-original-cell charge and the Boolean dual compiler. Its canonical hostile is the three-direction triangle with `beta=1 < tau_LP=3/2 < tau=2`. The least-used sidecar is how many **original overfull lines** meet at the same original cell.

The live concept board is: original-cell loss; overfull-line incidence; graph edge covers; fractional matching; disjoint charge regions; and concurrency. Targeted repository searches found the inherited two-family capacity-two matching but no earlier no3 edge-cover completion, blossom regime, or safe-cell block theorem. No external priority claim is made for graph matching or edge-cover identities.

The new map reverses the usual choice of vertices. Overfull lines become vertices, and an original cell meeting two such lines becomes an edge. A cell meeting three lines is a genuine hyperedge and must not be replaced by its three pairwise intersections. Keeping that incidence distinction completes a large exact regime while explaining why the unrestricted three-family relaxation remains nonintegral.

## 2. Arbitrarily many directions admit disjoint original-cell block charges

Let `B` be any finite set of distinct cells and `S` a finite family of directions. Write `tau_S(B)` for the minimum number of original cells deleted so that every line in those directions retains at most two cells. Let

`U_S(B)={p in B: every S-line through p has original occupancy at most two}`.

For any further finite direction family `T`,

\[
\boxed{\tau_{S\cup T}(B)\ge\tau_S(B)+\tau_T(B\cap U_S(B)).} \tag{1}
\]

All occupancies defining `U_S` are measured on the original `B`. Any deletion set clearing `S` still clears `S` after its deletions inside `U_S` are undone, because those cells touch no originally overfull `S`-line. Thus a repair of `S union T` deletes at least `tau_S(B)` outside `U_S`. Its deletions inside `U_S` must repair the restricted board for `T`, proving the additional independent charge.

More generally, partition a finite direction set into ordered nonempty blocks `T_1,...,T_k`, and put `S_j=T_1 union ... union T_j`, `U_0=B`, `U_j=U_(S_j)(B)`. Iterating (1) gives

\[
\boxed{\tau_{S_k}(B)\ge\sum_{j=1}^k\tau_{T_j}(B\cap U_{j-1}).} \tag{2}
\]

The charge at step `j` can be realized using only cells in `U_{j-1} setminus U_j`: cells in `U_j` lie on no originally overfull line of `T_j`, hence on no overfull `T_j`-line of the restricted board either. These strata are disjoint. For blocks of one direction the local cost is a sum of line excesses; for blocks of two directions it is the inherited exact unit-capacity bipartite flow. No probabilistic independence is needed. The bound is not claimed exact, and no order-independent optimum or greedy optimality is asserted.

There is also a useful dual form. If `Q` is any selected collection of old `S`-lines, let `N(Q)` be its incident original cells and `c(Q)=|N(Q)|−2|Q|`. Every repair gives

\[
\boxed{\tau_{S\cup T}(B)\ge c(Q)+\tau_T(B\setminus N(Q)).} \tag{3}
\]

At most `2|Q|` retained cells lie in `N(Q)`, and the outside cells must satisfy the new directions. Optimizing (3) with a singleton new direction recovers continuing11's Boolean compiler; a two-direction new block uses a genuine integral flow on the residual cells. This retains more information than a scalar old excess while keeping the two charges disjoint.

Any final no-three-in-line set `A`, even after adaptive insertions or moves, satisfies `|B setminus A|>=tau_S(B)`. Thus (1)–(3) retain the original-cell-loss interpretation and hold for every fixed skeleton and row order. This report makes no further probabilistic assertion.

## 3. General finite-direction graph regime: exact capacity-two b-matching

Let `V` be the collection of all originally overfull lines among the chosen directions. For each original cell `p`, let `I(p) subset V` be the lines containing it. Assume

\[
|I(p)|\le2\quad\text{for every original cell}. \tag{4}
\]

Cells with `I(p)` empty may always be retained. Make a graph with one vertex of capacity two for each member of `V`. A cell with two incidences becomes the edge joining those vertices. A cell with one incidence becomes an edge to its own new private leaf of capacity one. Each edge has capacity one and represents exactly one original cell.

Then retained feasible cell sets **containing all cells of empty incidence** correspond bijectively to graph b-matchings with these vertex capacities, together with those free cells. An arbitrary feasible retained set may omit free cells, but every optimal retained set contains them because adjoining them preserves feasibility. Consequently

\[
\boxed{\tau_S(B)=|\{p:I(p)\ne\varnothing\}|-
\max\{|M|:M\text{ is such an integral b-matching}\}.} \tag{5}
\]

The same construction is exact for the fractional relaxation. The proof is literal: capacity at an original line vertex says precisely that at most two of that line's cells are retained; a private leaf imposes no additional condition beyond the unit cell capacity. Originally safe lines impose no further restriction on any retained subset. This is a graph problem even when there are many directions. It need not be a bipartite flow problem; general matching retains the parity constraints that bipartite incidence would lose.

## 4. The all-triple regime has an exact closed matching formula

Assume in addition that every overfull line has **exactly three** original cells. Remove the private leaves from the graph above and call the resulting line-intersection graph `Gamma`. Its edges retain their original-cell labels. Put `v=|V|`, let `nu(Gamma)` be its maximum matching size, and let `alpha(Gamma)` be its independence number. The graph has maximum degree at most three; an isolated vertex has three private cells.

Let `beta_S=max_Q(|N(Q)|−2|Q|)` be the Boolean line-selection bound, allowing any selected line collection in the chosen directions. Let `tau_LP` be the minimum fractional deletion cost with `0<=x_p<=1` and each line's required deletion constraint. Then

\[
\boxed{\tau_S(B)=v-\nu(\Gamma),\quad
\tau_{\rm LP}(B)=v-\nu_f(\Gamma),\quad
\beta_S(B)=\alpha(\Gamma),} \tag{6}
\]

where `nu_f` is the maximum fractional matching value with edge weights nonnegative and incident sum at most one at every vertex.

**Integral identity.** A deletion must cover every overfull-line vertex at least once. Shared-cell deletions are graph edges; private-cell deletions cover one vertex. A maximum matching covers `2nu` vertices using `nu` cells. Every remaining vertex can be covered by one further incident or private cell, giving `v−nu` deletions. No edge joins two unmatched vertices, so these choices do not accidentally reduce the declared count.

Conversely take an inclusion-minimal cell cover. Its selected shared edges form disjoint stars: an edge with both endpoints already covered by other selected edges or private choices could be removed. A selected private cell cannot coexist with another selected incidence at its vertex. Choosing one edge from each nontrivial star gives a matching. If the cover has `k` cells, it has exactly `v−k` such stars, so `nu>=v−k` and `k>=v−nu`.

**Boolean identity.** A selected originally safe line can be removed from `Q` without decreasing the Boolean value, so it suffices to select vertices of `Gamma`. Since each selected line has three cells and no cell meets three overfull lines,

`|N(Q)|−2|Q|=|Q|−|E(Gamma[Q])|`.

An independent set achieves its cardinality. From any other `Q`, deleting at most one endpoint per induced edge leaves an independent set of size at least `|Q|−|E(Gamma[Q])|`. Maximization therefore gives exactly `alpha`.

**Fractional identity and an exact flow compiler.** First drop the deletion bounds `x_p<=1` without changing the optimum: capping every nonnegative `x_p` at one preserves each triple's constraint `sum x_p>=1` and can only decrease the objective. If a triple contains a capped variable, that variable alone still contributes one; otherwise its sum is unchanged. The dual of this equivalent formulation assigns `y_a>=0` to overfull-line vertices, with `y_a+y_b<=1` on shared edges and `y_a<=1` for private cells. The latter bound also follows at every nonisolated vertex from any incident edge and nonnegativity; isolated vertices have private cells. Thus its feasible set is precisely `0<=y<=1`, `y_a+y_b<=1`. Substituting `z=1−y` gives the fractional vertex-cover dual of maximum fractional matching. Finite LP duality gives `tau_LP=v−nu_f`.

Alternatively make the bipartite double cover of `Gamma`, with `a_L b_R` and `b_L a_R` for every edge. If its maximum matching has size `M`, then

\[
\nu_f=M/2,\qquad \tau_{\rm LP}=v-M/2. \tag{7}
\]

Any fractional matching gives a flow of value twice its weight in this double cover. Conversely symmetrizing an integral double-cover matching gives a feasible original fractional matching of value `M/2`. The inherited integral bipartite flow proves equality. A minimum double-cover vertex cover symmetrizes to a fractional vertex cover of the same value, yielding a matching fractional deletion dual. This compiles the precise fractional gap without an external numerical LP solver.

The scalar Boolean bound and the LP need not agree with integral repair, but the exact nonbipartite matching in (6) completes this entire declared regime. The failure is not an unrestricted impossibility of graph methods after a third direction is added.

## 5. Odd cycles, an actual five-cycle, and unbounded additive gaps

For a connected component `Gamma=C_(2j+1)`, its incident cells are disjoint from those of every other graph component. The formulas give

`beta=j`, `tau_LP=j+1/2`, `tau=j+1`.

Equivalently, summing the `2j+1` triple-deletion constraints and using that a cell covers at most two of them gives the rounded inequality `number of incident deletions >= j+1`. A matching attains it. This is a genuine integer parity charge. It may be summed over **separate graph components**, because their original-cell supports are disjoint. Counting arbitrary odd cycles without their attachment structure is not valid.

The continuing11 triangle is the `j=1` case. A new two-regular `6 by 6` board has zero-based row column-pairs

\[
(03,12,04,14,35,25). \tag{8}
\]

For directions `d=c−r`, `s=c+r`, `q=c−2r`, its only overfull lines are

`d=0, d=1, s=7, q=0, q=−5`.

Their adjacency cycle is `d0—q(−5)—s7—d1—q0—d0`, with distinct original cells at its intersections. Every one of these lines has three cells. Thus

\[
\boxed{\tau_2=2,\quad J_3=0,\quad\beta_3=2,
\quad\tau_{\rm LP}=5/2,\quad\tau_3=3.} \tag{9}
\]

The original safe-set and Boolean bounds miss the additional required deletion, while the matching theorem recovers it exactly. Exhaustive smaller-board controls establish that `n=6` is the smallest square dimension supporting a five-cycle component in this all-triple, at-most-two-incidences regime; the `n=6` search itself was bounded and not exhaustive.

There are unbounded additive gaps already for these **fixed three directions** on finite integer boards. Take `k` copies of (8), translating copy `j` by `(100j,300j)`. Their line labels shift by `200j,400j,100j` in the three respective families, larger than the entire label ranges of one copy, so no line is shared across copies. Hence `Gamma` is `k` disjoint five-cycles and

`tau_3=3k`, `tau_LP=5k/2`, `beta_3=2k`.

Every occupied row and column still contains two cells. There are empty rows and columns in the bounding square, however; these are not saturated `n by n` boards and are not a random-board lower-density construction.

## 6. Sharp hostiles and what must remain in the state

**A third incidence is a real hyperedge.** The two-regular `5 by 5` board

`(24,01,23,34,01)`

has exactly three overfull triples: `d=0`, `s=4`, `q=−2`. All three meet at the **same** original cell `(2,2)`. Deleting that cell repairs them all, so `tau_3=1`. Their naive pairwise line-intersection graph is a triangle, whose ordinary edge-cover formula would predict two. Hypothesis (4) is therefore load-bearing; the first failed implication is replacing one three-way deletion choice by three different two-way choices. This is dimension-minimal for the fixed three directions, since slope-two lines have at most two grid cells for `n<=4`.

**The presence of an odd cycle is not enough.** The two-regular board

`(01,02,34,13,24)`

lies in the theorem's regime but its line graph is a triangle with one pendant vertex. It has a perfect matching and `beta=tau_LP=tau=2`. The attachment absorbs the apparent parity defect. A local odd-cycle count without boundary incidence would overcharge; (6) retains the full graph.

**Ordinary edge covers require triple lines.** If a line has four cells it needs two deletions, so a one-cover-per-vertex model is invalid. Formula (5) still retains its capacity-two b-matching exactly when (4) holds. The all-triple hypothesis is used only for the closed formulas (6), not for the general graph correspondence.

The successful connection is therefore: arbitrary finite direction constraints -> original overfull-line incidences -> capacity-two graph b-matching when cells have at most two incidences -> exact matching/edge-cover formulas when each overfull line is a triple. It preserves the complete deletion predicate. The discarded three-way incidences and larger line demands are exhibited separately, rather than silently imposed as filters on a claimed general theorem.

## 7. Exact branching only on genuine hyperedge cells

There is a paid extension beyond hypothesis (4). Put

`H={p:|I(p)|>=3}`, `h=|H|`, `C={p:1<=|I(p)|<=2}`.

Retain the complete incidence set of every cell in `H`. For each subset `R subset H` proposed to be retained, reject the branch if some overfull line contains more than two cells of `R`. Otherwise give each old line vertex residual capacity

`b_a(R)=2−|R intersection a|`.

Make the graph from cells of `C` exactly as in Section 3, with these residual capacities and private-leaf capacities one. Let `mu(R)` be its maximum integral b-matching size. Then for **every** original finite board and finite direction set,

\[
\boxed{\tau_S(B)=h+|C|-\max_{R\subseteq H\ \mathrm{valid}}
\bigl(|R|+\mu(R)\bigr).} \tag{10}
\]

Every retained cell set fixes exactly one such `R`; the remaining capacity constraints are precisely its residual graph b-matching. Conversely every valid branch and b-matching gives a feasible retained set, together with all cells of empty incidence. This proves (10). Thus there are at most `2^h` exact graph optimization calls, branching only on genuine higher-incidence original cells rather than on every cell. No bound on `h` for random or saturated boards is claimed. The graph optimization is integral; a fractional replacement would reintroduce the already demonstrated parity gap.

For the concurrent-triple hostile, the shared central cell is exactly the exceptional set `H`. Its original three-way action survives the branch, so (10) recovers the one-deletion solution that a pairwise clique replacement loses.

## 8. Declared finite universe and stopping boundary

The exploratory probe first exhausted all 2,137 two-regular boards with `n=2..5`. Exactly 1,320 satisfy the all-triple, at-most-two-incidences regime (counts `1,6,73,1240`). Six `n=5` boards have a triangle component; none has a five-cycle. It then used random seed `12092026`, up to 20,000 pairs of shuffled permutations per dimension starting at `n=6`, stopping on the first five-cycle component. It stopped at draw 11,850: 4,384 pairs were simple two-regular boards, 2,098 entered the declared regime, and the first five-cycle was (8). These draws were a bounded discovery device, not a uniform uncolored-board sample or an exhaustive `n=6` census.

The [standalone producer](../../04-computation/continuing12_20260908_no3_matching_completion.py) passes **58,304 always-active exact gates**. It verifies the full exceptional-cell branching formula on all 2,137 smaller boards, without filtering by the graph regime. Sixteen `n=5` boards have genuine higher-incidence cells. It separately verifies (6) on every one of the 1,320 eligible boards by comparing exact triple-hitting deletion search with matching, original Boolean union enumeration with independence, and literal fractional deletion primal/dual certificates from the double cover.

An independent abstract bank contains all **844 labelled subcubic graphs on zero through five vertices**, with counts `1,1,2,8,64,768`. Each graph produces an exact triple-incidence system using a shared cell per edge and `3−degree` private cells per vertex. The same independent deletion and primal/dual checks apply. These incidence systems are not all claimed realizable by geometric grid lines.

The geometric control bank retains the five explicit positive/hostile boards above, testing every retained subset against the residual-capacity branch predicate: **8,192 subsets** in total. It also verifies their complete four-direction branch calculations. For the 97 complete boards of sizes two through four and the five named controls, it checks all 24 singleton block orders and all six ordered partitions into two blocks of two directions. These are **3,060** block inequalities with disjoint supporting strata. There are another **32** two-new-direction residual-dual checks on the named boards. The separated-copy control checks one through four literal translated five-cycle boards and confirms that occupied rows and columns each have two cells.

Normal and optimized Python runs have byte-identical LF stdout and regenerated certificates. The [certificate](continuing12_20260908_no3_matching_completion_certificate.json) contains the complete universe counts, graph matchings, fractional primal/dual vectors, exceptional retained choices and all gate counts. The [stored transcript](continuing12_20260908_no3_matching_completion.out) records the bounded scope. The implementation uses explicit exceptions, exact integers and rational fractions, no numerical LP solver, and no imported mathematical producer. Its finite graph optimization uses exact bounded dynamic programming; this is not a timing claim about the unbounded reduction.

Reproduce after filing by running the standalone source normally and with `-O`. An outside source writes the certificate beside itself; a source under `04-computation` writes it under `05-knowledge/results`. The unbounded statements above use the proofs, not finite extrapolation.

The matching regime itself is now exact, and (10) supplies an exact extension parameterized by the number of genuine higher-incidence cells. The next structural frontier is to control that branch parameter or exploit the overlap among its exceptional incidence sets. No small parameter bound, improvement in random-board mean repair, or uniform probabilistic estimate is claimed here.

## Independent acceptance

The [independent audit](continuing12_20260908_no3_matching_audit.md) accepts
this scoped result and its analytic proof, with 12767 exact gates per
normal and optimized pass. Its producer-independent controls retain the
declared positive and hostile boundaries. The parent replays both modes
against frozen raw LF outputs and certificates.
