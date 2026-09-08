# Cell separators localize the genuine hyperedge repair problem

**Status: PROVED ANALYTICALLY + FINITE-EXACT controls; INDEPENDENTLY AUDITED.** This is an exact reduction for any finite original board and any finite set of directions. It localizes continuing12's exceptional-cell branching, and gives connected integer-geometric families with arbitrarily many exceptional cells but no local branching. It asserts no random-board frequency, no asymptotic saturated-board construction, and no improvement to the extremal no-three-in-line constant.

## 1. Inheritance and precise new coordinate

The closest supplier is [continuing12 matching completion, Sections 3 and 7](continuing12_20260908_no3_matching_completion.md): cells incident to at most two originally overfull lines compile to integral capacity-two graph b-matching; preserving all other cells gives an exact bound of `2^|H|` graph calls. The triangle requires integral nonbipartite matching; a concurrent triple requires one genuine three-way cell; and a triangle with an attachment prevents charging odd cycles without their boundary state. Those are the inherited positive and hostile controls. The least-used sidecar is the topology of the **bipartite original line--cell incidence graph**, before projecting cells to edges.

The live board is: original-cell loss, literal line constraints, exceptional concurrency, articulation cells, weighted matching, and inseparable joint choices. A targeted search of the previous no3 result reports found no earlier separator-localized repair statement. No external priority claim is made for block-cut trees, tree dynamic programming, or weighted graph matching.

The new source-to-target map is

`original repair -> complete line constraints in incidence bags -> two-state original-cell separators -> weighted integral graph optimization inside each bag`.

The retained/deleted decision of a shared original cell is preserved exactly. A scalar local repair cost would lose how that cell interacts with the next bag; the two conditional costs are the required sidecar. The remaining exponential parameter counts genuine higher incidence **inside one inseparable bag**, rather than across the entire board.

## 2. The correct bags keep every line constraint whole

Let `V` be all originally overfull lines in the selected directions, and `P` all original cells lying on at least one such line. Define the bipartite graph `X` with vertices `V disjoint-union P` and one edge for each literal line--cell incidence. Cells outside `P` can and will all be retained in an optimum; arbitrary feasible retained sets need not contain them.

Take the usual biconnected blocks of `X`, including every bridge as a two-vertex block, and its block-cut forest. Merge all blocks connected through an articulation vertex belonging to `V` (a **line vertex**). Do not merge through articulation cell vertices. The merged blocks are called bags. The remaining articulation cell vertices link bags in a forest.

This construction has four properties:

1. Each incidence edge belongs to exactly one bag.
2. Each original line vertex belongs to exactly one bag, and that bag contains every cell of the line.
3. A nonseparator cell belongs to exactly one bag. A separator cell can belong to several bags, with one shared binary state.
4. The bipartite bag--separator graph is a forest.

Indeed the block-cut graph is a forest. Contracting each connected group of blocks joined through line articulation nodes preserves that property. All blocks containing a given line are in the same group, which proves complete line ownership. Blocks can share only articulation vertices; the only uncontracted ones are cells. These facts also cover disconnected incidence graphs and bridge-only components. If there are no overfull lines, the repair cost is zero and the forest is empty.

For a bag `B`, let `I_B(p)` be the original line vertices owned by `B` containing `p`, and set

\[
h_B=|\{p\in B:|I_B(p)|\ge3\}|. \tag{1}
\]

This is a local incidence count. A cell with many incidences globally may have only one or two in each bag.

It is not valid to obtain bags simply by deleting all articulation cells and taking the remaining connected components. Two such cells can lie on a common cyclic block; the corresponding quotient can have a cycle. Section 5 gives an actual three-direction geometric example of this failure. The block-cut construction retains the cyclic block whole.

## 3. Exact two-state recurrence and weighted matching compiler

Root each component of the bag--separator forest at a bag. Write `s_p=1` when original cell `p` is retained and `s_p=0` when it is deleted. For a separator cell `p` with parent bag, let its message be

\[
C_p(s)=1-s+\sum_{B\text{ child of }p} C_B(s). \tag{2}
\]

Here `C_B(s)` is the minimum deletion cost inside the descendant bag subtree when its parent separator is fixed to state `s`, excluding the deletion cost of that parent cell. Thus (2) pays each separator cell exactly once.

Inside bag `B`, fix its parent separator if present. For every other cell `q` in the bag, give it unary costs

\[
(c_q(0),c_q(1))=
\begin{cases}
(1,0),&q\text{ is not a separator},\\
(C_q(0),C_q(1)),&q\text{ is a child separator}.
\end{cases} \tag{3}
\]

Minimize their sum, subject to the complete owned line constraints

\[
\sum_{p\in a}s_p\le2\qquad(a\text{ owned by }B). \tag{4}
\]

The fixed parent cell has no cost in (3), since its cost is paid upstream. Since each constraint has a unique owner and all overlaps are the declared separator cells, induction from the leaves proves that (2)--(4) equal the original minimum deletion cost. Sum the unconditional root-bag minima over components. The recurrence is finite even for parent state one: retaining one fixed cell and deleting every other descendant cell is feasible.

There is an exact graph compiler for each bag state. Put

`base=sum_q c_q(0)` and `w_q=c_q(0)-c_q(1)`.

Retaining cell `q` earns reward `w_q`. Enumerate only the free cells for which `|I_B(q)|>=3`; a fixed parent cell is already decided. For each assignment, reject it if the fixed retained cells exceed two on any owned line. Otherwise reduce each line's capacity by the number of those fixed retained cells. Every remaining decision cell has one or two incidences in `B` and becomes a unit-capacity graph edge of weight `w_q`, with a private capacity-one leaf when it has only one incidence. Maximize the total weight of an **integral** b-matching, with the residual line capacities in `{0,1,2}`. Subtract this reward and the chosen high-incidence rewards from `base`.

Weights may be negative. Retaining a separator can force additional deletions below it, so an unweighted maximum-cardinality matching is invalid. Optional negative-weight edges may be omitted by the weighted optimizer. No bipartite-flow claim is made: the inherited triangle still requires nonbipartite integral matching.

Every feasible retained set gives exactly one assignment of the local exceptional cells, followed by exactly one feasible residual graph edge set. Conversely every such branch and edge set satisfies all owned original line constraints. This is the same literal cell correspondence as the supplier, now with exact continuation costs rather than unit rewards.

Each nonroot bag has at most two parent states; each state uses at most `2^{h_B}` graph calls. Therefore

\[
\boxed{\text{graph calls}\ \le\ 2\sum_B2^{h_B}.} \tag{5}
\]

There is only one root state, and fixing a high-incidence parent cell saves a factor two, so (5) is a convenient uniform upper bound. The bag construction and message assembly have polynomial overhead in the finite incidence size. If `p=max_B h_B`, this gives `O(number_of_bags * 2^p)` weighted integral graph calls. The previous global method can also be selected; the new call bound need not numerically dominate `2^|H|` on every small input.

For **three directions**, each exceptional cell has exactly three incidences. If such a cell is an articulation vertex, its incidences split between at least two bags, so it contributes to no `h_B`. Otherwise it lies in exactly one bag and contributes there. Thus the three-direction local parameter counts the exceptional cells that are **not articulation cells**, separated by bags. Arbitrarily many exceptional articulation cells cause no exponential branching.

For more than three directions, an articulation cell may still have three or more incidences inside one bag. It must then remain in that bag's branching state. Dropping all articulation cells from the exceptional count is a three-direction corollary, not the general theorem.

## 4. A connected geometric family has unbounded global concurrency and zero local branching

Use fixed slopes `1,-1,2`, with labels `d=c-r`, `s=c+r`, and `q=c-2r`. For every integer `k>=1`, start with centres

\[
P_i=(2i,\,2i+2\lfloor i/2\rfloor),\qquad 0\le i<k. \tag{6}
\]

Select all three direction lines through every centre, identifying repetitions. Consecutive centre differences alternate `(2,2)` and `(2,4)`, hence successive joining slopes alternate one and two. More explicitly their labels are

`d_i=2 floor(i/2), q_i=-2 ceil(i/2), s_i=4i+2 floor(i/2)`.

Each repeated `d` or `q` label occurs at exactly one adjacent pair; the `s` labels are all distinct. There are exactly `2k+1` selected lines, connected through the centres in a chain. Each centre belongs to exactly three selected lines; each selected line currently has either one or two centres.

Fill each selected line to exactly three points with private integer points. Such a filling can ensure that every nonselected line in the three directions is a singleton, and that no new point shares a row or column with a previous point. To see this, for a selected slope `a` and label `ell`, candidates are `(r, ar+ell)` for integer `r`. Excluding existing rows and columns forbids finitely many values of `r` because `a` is nonzero. For each other slope, exclude every already used label and every other selected line label. Distinct slopes make each excluded equality forbid at most one `r`. Choose any remaining integer and iterate. Every centre is installed before filling, so no later private point lands on another selected line or recreates a previous nonselected label. A common translation then puts the board in the nonnegative integer quadrant.

There are `3k+3` private points and `4k+3` original points total. The active incidence graph has `6k+4` vertices and `6k+3` edges and is connected, hence is a tree. Its bags each own one complete selected line; every cell has only one incidence within each bag. Consequently

\[
|H|=k,\qquad h_B=0\text{ for all }B,
\qquad\boxed{\tau=k.} \tag{7}
\]

The upper bound deletes all centres. For the lower bound, the `k` distinct slope-minus-one lines each need a deletion, and their original-cell supports are pairwise disjoint. Thus these are actual connected geometric instances where the old global branch count is `2^k` but the new bound is `4k+2` graph calls. With the declared rooting, the producer uses `4k+1` calls.

Every occupied row and column has one point; the bounding square can contain gaps. This is an arbitrary finite integer board, not a saturated two-regular board. No assertion about positive frequency under a random board law follows from this construction.

Even `k=1` forces a negative separator reward. Root at one line. The central cell's two child lines each cost zero if the cell is deleted and one if it is retained. Including its own deletion cost gives `C_p(0)=1, C_p(1)=2`, so `w_p=-1`. The weighted continuation sidecar is load-bearing in the smallest family member.

## 5. Geometric hostiles retain the failure boundary

**Two inseparable exceptional cells can improve only jointly.** An explicit board is

```
(6,12), (12,12), (9,15), (10,8), (0,0), (8,10),
(7,17), (11,11), (13,19), (15,18), (14,28).
```

Its six overfull lines, with zero-based original point indices, are

```
s=18: {0,3,5}     s=24: {1,2,6}
d=0:  {1,4,7}     d=6:  {0,2,8}
q=-12:{1,3,9}     q=0:  {0,4,10}.
```

Only points zero and one have three incidences. They lie in one bag, whose incidence core is a theta with three internally disjoint paths between those cells; hence `h_B=2`. If their forced **deletion** states are `00,01,10,11`, the exact respective costs are

\[
\boxed{3,4,4,2.} \tag{8}
\]

When neither centre is deleted, points two, three, and four give three pairwise line covers, and no remaining point covers more than two of the six lines. When exactly one centre is deleted, the other three still-uncovered lines have pairwise disjoint available cells, so need three further deletions. When both centres are deleted, every overfull line is repaired. Thus screening out a cell because its unilateral deletion gives no improvement is false. The missing coordinate is the joint choice inside an inseparable incidence block; the localized theorem retains it. The minimum is two, but a greedy one-cell improvement from the `00` state cannot reach it.

**Deleting all articulation cells does not give a forest quotient.** Another actual board is

```
(0,2), (6,8), (2,0), (1,1), (3,11), (9,5),
(5,3), (12,10), (4,6), (8,12), (7,16), (10,22).
```

Its three centres are points zero, one, and two. The base lines `s=2,d=2,q=-4` form a six-cycle alternating these centres and line vertices. Each centre also has a separate overfull arm line, respectively `q=2,s=14,d=-2`. These are exactly the six overfull triples. All three centres are articulation cells. Removing them first leaves the three base-line components; reconnecting those components to the separator cells recreates the six-cycle. The naive quotient therefore is not a tree, despite all separators being genuine articulation cells.

The actual block-cut construction has four bags: the whole cyclic base and the three separate arm lines. All `h_B=0`; the exact cost is three, forced by the three pairwise disjoint arm lines and attained by deleting all centres. Cyclic bags with only local two-way cells remain exact graph problems. Splitting a cyclic bag into independent line costs would lose its attachment information.

The inherited triangle, concurrent triple, triangle with pendant, and larger-line control remain in the finite bank. Together they separate integral from fractional matching, a literal hyperedge from a projected clique, attachment state from cycle counting, and capacity-two retention from the special all-triple edge-cover formula.

**The articulation simplification stops at three directions.** Add slope minus two to the theta construction and include its line through each centre, filling each to a triple with private points. The certificate supplies the literal fifteen-point integer board. The two centres now each have four incidences and are articulation vertices because each has a separate new arm. Both still have three incidences inside the same theta bag, so its `h_B=2`. There are three bags total, with local counts `0,0,2`, and the exact repair cost is two. Testing all three root choices exercises a fixed parent separator of local incidence three. Thus the general theorem must use (1); globally being an articulation cell alone does not permit suppressing its genuine local hyperedge when four directions are present.

**The first small saturated local core occurs in dimension six.** Every one of the 2,137 complete two-regular boards of dimensions two through five has `h_B=0` for the three fixed directions. This includes all sixteen dimension-five boards with genuine global exceptional cells. In dimension six the board

`(01,01,23,45,35,24)`

has one bag with `h_B=1`, at cell `(3,4)`, and exact repair cost four. Its overfull lines are `s=7, d=-1, d=0, d=1, q=-2, q=-1`; the `d=1` line has four cells, and the others are triples. Consequently dimension six is the smallest square dimension admitting a nonzero local parameter in this saturated class. A lexicographic discovery scan stopped at its fourteenth dimension-six board; this is not an exhaustive dimension-six census. This finite frontier is compatible with the connected unsaturated family: the class and incidence geometry remain explicit.

## 6. Scope and decisive next boundary

The theorem works for every finite board and finite direction set, including lines with more than three cells. It preserves the original-cell-loss predicate: if a final no-three-in-line set `A` is obtained using arbitrary adaptive insertions or moves, then `B intersect A` is a feasible retained subset, so `|B setminus A|>=tau(B)` and all computed costs remain valid original losses.

The connected family shows that the raw global exceptional-cell count alone is too coarse to measure repair complexity. The theta witness shows the opposite boundary: one must not erase interacting exceptions merely because every one-cell test fails. The remaining frontier is the genuine high-incidence core inside a bag. No small local parameter under the uniform two-regular board law is claimed. The finite saturated-board bank tests correctness of the reduction on that class; it cannot establish an asymptotic frequency or mean improvement.

## 7. Reproduction and finite controls

The [standalone source](../../04-computation/continuing13_20260908_no3_cell_separators.py) passes **45,383 always-active exact gates** and uses two separate exact routes: whole-instance triple-hitting deletion recursion, and a Tarjan block decomposition with line-block contraction, two-state messages, and weighted residual-capacity graph dynamic programming. It imports no mathematical producer. All gates use explicit exceptions and remain active under optimized Python. The finite graph routine is a bounded exact verifier, not a timing claim about the unbounded oracle reduction.

The complete actual-board universe is all 2,137 simple two-regular boards for `n=2..5`, with counts `1,6,90,2040`. Every board is tested without filtering out higher incidence. Sixteen dimension-five boards have global exceptional cells, all of them articulation cells; twelve boards generate negative continuation rewards. A separate complete abstract universe contains all 109 linear triple-incidence systems on zero through four line vertices, with counts `1,1,2,9,96`: choose any shared cells of size two or three, require each pair of lines to occur together in at most one cell and shared degree at most three, then fill each line to a triple with private cells. These are incidence controls, not all claimed geometric realizations.

Six named saturated boards (the five inherited controls and the new dimension-six core) are checked in three and four directions. The connected geometric chain family is checked for every `k=1..8`, including all claimed counts, connectivity, exact costs, singleton nonselected lines, and one point per occupied row and column. The two-cell theta, cyclic naive-quotient hostile, and actual four-direction boundary are checked independently against their original cells. Every bag-root choice is replayed for the abstract bank, named controls, and constructed geometric examples. Literal line ownership, separator membership, component count, and the forest identity are checked; articulation cells are recomputed independently by deleting each cell vertex and comparing connected-component counts.

The [stored transcript](continuing13_20260908_no3_cell_separators.out) and [certificate](continuing13_20260908_no3_cell_separators_certificate.json) record the exact gate count, declared universes, complete constructed coordinates, negative continuation rewards, local parameters, and all hostile tables. Running the source normally and with `-O` gives identical raw LF stdout and regenerated certificates. A source placed in `04-computation` writes its certificate to `05-knowledge/results`; an outside source writes beside itself. Unbounded claims above follow from the analytic proofs, not from finite extrapolation.

## Independent acceptance

The [independent referee](continuing13_20260908_no3_separators_audit.md) accepts the complete scoped theorem without mathematical repair, with **85,491 always-active exact gates per mode**. Its alternative route finds blocks by recursive vertex deletion and optimizes local states by literal retained-subset enumeration. It independently recovers all 2,137 complete boards, all 109 abstract systems, the constructed coordinates, every tested root choice, the first saturated local core, and the four-direction fixed-parent boundary. Normal and optimized referee runs have identical raw LF stdout and certificates.

The referee also proves a useful sufficient bound: if the ordinary nonbridge biconnected blocks in a bag have cycle ranks `beta_i`, then `h_B<=sum_i(2 beta_i-2)`. Handshaking in each such block proves it, since every vertex has degree at least two and every locally exceptional cell contributes at least one to the total degree excess. In particular incidence cacti have `h_B=0` for every finite direction family. This is a sufficient condition; a subdivided `K4` with branch vertices of line type has `h_B=0` and is not a cactus. The typed line--cell incidence remains finer than untyped cycle complexity.
