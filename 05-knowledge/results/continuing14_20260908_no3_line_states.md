# Line capacity states keep the incidence blocks separate

**Status: PROVED ANALYTICALLY + FINITE-EXACT controls; INDEPENDENTLY AUDITED.** This extends the continuing13 exact original-cell repair reduction. It keeps the ordinary biconnected blocks separate by paying both kinds of separator state: a binary retained/deleted cell state, and a three-valued retained-count state for a line. An explicit weighted ordinary-matching gadget preserves arbitrary continuation costs at a capacity-two line, including nonconvex costs. A connected fixed-three-direction geometric family exhibits an exponential reduction in the inherited graph-call bound. No random-board frequency, extremal no-three-in-line constant, or polynomial algorithm for unrestricted inseparable hyperedge cores is claimed.

## 1. Inheritance: merging a line was an information choice

The closest proved supplier is [continuing13 cell separators](continuing13_20260908_no3_cell_separators.md), including its [independent cycle-excess/cactus audit](continuing13_20260908_no3_separators_audit.md). It merges ordinary incidence blocks through articulation **line** vertices so that every line constraint stays whole, and then transmits two states through articulation **cell** vertices. Its remaining branching parameter sums the genuinely higher-incidence cells across such merged bags.

The inherited triangle requires integral graph matching, the concurrent triple retains one literal three-way cell, and the geometric theta has two exceptional cells that improve only jointly. The [continuing12 matching theorem](continuing12_20260908_no3_matching_completion.md) is the original graph consumer. The least-used sidecar is the number of retained cells contributed to one articulation line by each incident block.

The new map is

`original incidence -> ordinary block-cut forest -> binary cell states plus exact line counts -> weighted degree-cost graph -> ordinary weighted matching`.

It preserves every original cell identity, the capacity-two constraint on every originally overfull line, and the complete continuation cost for each separator state. A single scalar repair value per piece loses the number of cells that the piece leaves on a shared line. Treating line slots as independent linear costs loses nonconvex continuation profiles. Sections 3 and 5 give explicit repairs and hostiles for these two losses.

The live board is: original-cell loss; typed incidence; full separator states; nonconvex degree costs; integral matching; and inseparable hyperedge choices. This is a scoped reduction in the repository's repair problem. No external priority claim is made for block-cut trees, finite-state dynamic programming, or weighted matching.

## 2. Exact recurrence on the unmerged block-cut forest

Let `V` be the originally overfull lines among any finite set of chosen directions, and let `P` be the original cells on at least one line in `V`. Use the literal bipartite incidence graph `X` on `V disjoint-union P`. Cells outside `P` can all be retained in an optimum; an arbitrary feasible retained set may omit them.

Take the ordinary biconnected blocks of `X`, with each bridge included as a two-vertex block. **Do not merge through either type of articulation vertex.** Their block-cut graph is a forest. Each incidence edge has exactly one block owner, while an articulation vertex can belong to several blocks. The state at an articulation cell is its one original retained/deleted decision. The state at an articulation line records the number of its retained incident cells contributed by a specified block, in `{0,1,2}`.

Root each block-cut component at a block. For a separator cell `p` and parent block, let

\[
C_p(s)=1-s+\sum_{B\text{ child of }p} C_B(p,s),
\qquad s\in\{0,1\}. \tag{1}
\]

Here `C_B(p,s)` excludes the deletion cost of its fixed parent cell. A nonseparator cell pays its deletion cost inside its unique block. Thus every original cell is charged exactly once.

For a separator line `a`, let `t` be its parent block's retained contribution. The line message is

\[
C_a(t)=\min_{\substack{r_B\in\{0,1,2\}\\
\sum_B r_B\le2-t}}
\sum_{B\text{ child of }a} C_B(a,r_B),
\qquad t\in\{0,1,2\}. \tag{2}
\]

The block message `C_B(a,r)` requires **exactly** `r` retained incidences at the parent line inside that block. It excludes the other blocks of this line. The minimum in (2) is a capacity-two min-plus convolution, so only the three budget states `0,1,2` are needed regardless of the number of child blocks. A block message for an impossible exact count is infinite; every line message itself is finite, since all child cells can be deleted.

Inside a block, every nonparent cell has unary costs `(1,0)` if internal, or `(C_p(0),C_p(1))` if a child separator. A parent cell is fixed and has no local deletion cost. Each child separator line adds the continuation cost `C_a(k)` when the block retains `k` of its incident cells. A nonseparator line has its ordinary local capacity two. A parent line has its prescribed exact count. All local counts are at most two.

These recurrences are exact. Every original line--cell incidence belongs to one block, so summing the block contributions at a line counts every retained original cell on that line exactly once. The shared cell states guarantee that a cell appearing in several blocks makes one consistent original decision. Conversely states satisfying (1)--(2) and every block predicate glue to a feasible original retained set. Induction over the forest therefore gives the original minimum deletion cost.

The old merged-line construction kept the line constraint whole geometrically. Formula (2) now keeps the same constraint whole **in the state**, permitting a finer decomposition.

## 3. Arbitrary capacity-two degree costs compile to ordinary matching

First fix all local cells incident to at least three line vertices of the block, and fix a parent cell if present. Reject a choice that already exceeds any line's permitted count. Each line then has residual capacity in `{0,1,2}`, and a parent line may additionally require an exact residual count. Every remaining cell has one or two local incidences.

Write a decision cell's costs as `c_p(0),c_p(1)`. Start with the constant `sum c_p(0)` and give retention the reward `w_p=c_p(0)-c_p(1)`. At each line, its child continuation supplies a reward depending on the number of further retained cells. The resulting graph problem maximizes

\[
\sum_{p\text{ retained}}w_p+\sum_a R_a(k_a), \tag{3}
\]

where `0<=k_a<=b_a<=2`, each original cell is used at most once, and there can be one prescribed exact degree at the parent line. The line reward profile can be shifted by already retained fixed cells; its entries are still finite integers. In particular it must not be assumed linear or convex.

Here is an explicit ordinary weighted-matching compiler for any such profile. For a capacity-two line with rewards `R_0,R_1,R_2`, put

\[
L=\min(R_1-R_2,R_0-R_1),\quad
a=R_1-R_2-L,\quad b=R_0-R_2-2L. \tag{4}
\]

Then `b>=a>=0`. Make two slot vertices. Join them by an edge of weight `b`, and join each to **one shared private leaf** by an edge of weight `a`. Charge `-L` on every selected incidence using a slot and add the constant `R_2+2L` to the objective. If `k` slots are occupied externally, the best remaining gadget weight is respectively `b,a,0` for `k=0,1,2`. Indeed two free slots can use their internal edge or only one edge to the shared leaf, and `b>=a`; one free slot can use its leaf edge. Thus the total line reward is exactly

`R_2+2L-L*k+(b,a,0)_k=R_k`.

A capacity-one line needs one slot, the constant `R_0`, and the incidence charge `R_1-R_0`. A capacity-zero line contributes the constant `R_0` and has no slots. These constructions retain all three entries, without a convexification or fractional relaxation.

Each ordinary original cell gets **two private port vertices**, joined by an internal edge. A cell with two line incidences connects one port to each slot of its first line and the other port to each slot of its second line. A cell with one line incidence connects its first port to that line's slots and its second port to a new private free vertex. Put its reward `w_p` on the external edges from the first port, once; add each endpoint line's incidence charge. All internal port-pair edges initially have weight zero.

Require all cell ports to be covered by the matching. This can itself be enforced with a finite integer weight. Add `M` for each covered port to every incident edge: an internal port-pair edge gets `2M`, an external port edge gets `M`, and line-gadget edges get none. Choose `M` greater than the sum of absolute values of all other edge weights. Covering every cell port is feasible by using all internal port-pair edges. Any matching leaving even one port uncovered loses more coverage reward than its entire possible gain elsewhere, so every maximum-weight matching covers all ports.

Consequently each original cell is either unused, with its two ports matched internally, or retained, with both ports matched externally. A half-retained cell is impossible. Shared line slots enforce their capacities. The remaining matching reward, after subtracting coverage bonuses and adding the declared constants, is exactly (3).

To enforce one exact degree `r` at the parent line, first truncate its capacity and profile to `r`. Give each selected incidence there a bonus `N` exceeding the sum of absolute base edge weights. Among matchings covering all cell ports, this maximizes that line's degree first; check whether it reaches `r`, and reject the branch if not. Choose the coverage bonus `M` only **after** this degree bonus, larger than the sum of all resulting absolute noncoverage weights. Thus port coverage remains the first priority, the exact parent degree the second, and the actual objective the third. Subtract both artificial bonuses after decoding. These weights have polynomial bit length in the finite integer input.

The compiler is exact for positive, zero, or negative cell rewards and arbitrary integer line profiles. It uses an ordinary **integral weighted matching** oracle, so the inherited odd-cycle obstruction remains represented.

## 4. The improved parameter and its boundary

For an ordinary unmerged incidence block `B`, set

\[
h_B=|\{p\in B:\ p\text{ has at least three incidences inside }B\}|. \tag{5}
\]

Enumerate only these local high-incidence cells. A parent cell is already fixed and need not be branched. A parent cell has at most two states, a parent line at most three, and the root only one. The recurrence and compiler therefore give

\[
\boxed{\text{ordinary weighted-matching calls}
\ \le\ 3\sum_{B\text{ ordinary block}}2^{h_B}.} \tag{6}
\]

The remaining work per enumerated high-cell assignment and finite-state convolution is polynomial; the total runtime has a polynomial factor times the displayed sum of local branching terms. The new bound need not be numerically smaller on every small input; one can also choose either inherited reduction. Its strict gain is that local high-incidence counts are no longer added together merely because cyclic blocks meet at a line.

For a nonbridge ordinary block with cycle rank `beta_B`, the continuing13 handshaking argument gives `h_B<=2 beta_B-2`: every block vertex has degree at least two, so each counted cell consumes at least one of the `2 beta_B-2` total degree excesses. Bridge blocks have `h_B=0`. In particular a bound `beta_B<=b` on each nonbridge block gives at most `3 * number_of_blocks * 2^(max(0,2b-2))` graph calls, regardless of the total cycle rank of the whole incidence graph. Incidence cacti still require no exponential branching. This is a sufficient condition: cells of local incidence two can compile to a graph even when the untyped block has large cycle rank. The roles of line vertices and cell vertices remain distinct.

The genuine inseparable boundary is retained. The two exceptional cells in the geometric theta occupy one ordinary block and still require a joint choice. The dimension-six saturated control `(01,01,23,45,35,24)` retains an ordinary block with `h_B=1`. The new line state does not remove these genuine local hyperedges or imply a uniform small parameter on arbitrary two-regular boards.

## 5. Exact geometric gains and load-bearing hostiles

### 5.1 Connected theta blocks with one shared original line

Start with the continuing13 eleven-point theta board

```
(6,12),(12,12),(9,15),(10,8),(0,0),(8,10),
(7,17),(11,11),(13,19),(15,18),(14,28).
```

For every integer `k>=1`, take the union of its `k` copies under

\[
(r,c)\longmapsto(1000^j(r+100),\,1000^j(c+100)),
\qquad0\le j<k. \tag{7}
\]

Use the fixed directions of slopes `1,-1,2`. The slope-one line `d=c-r=0` is shared by all copies and contains exactly `3k` original cells. Every other originally overfull line belongs to one copy and is a triple. There are no other cross-copy line coincidences: nonzero slope-one labels are bounded nonzero integers multiplied by different powers of 1000; slope-minus-one labels lie in `[200,256]*1000^j`; slope-two labels lie in a fixed negative interval inside `[-200,-1]*1000^j`. These scaled intervals and the corresponding nonzero signed integer label ranges are pairwise disjoint across distinct `j`. The occupied row and column intervals are disjoint across copies as well.

There are `11k` points and `5k+1` overfull lines. Its incidence graph is connected. Each copy contributes one nonbridge theta block and six bridge blocks for its private cells; the theta blocks meet only at the common line vertex. Thus there are `7k` ordinary blocks, with `k` local values equal to two and `6k` equal to zero. Formula (6) gives at most **`30k`** ordinary matching calls.

By contrast all these blocks merge into one bag under the continuing13 line-contraction rule, with `2k` exceptional cells. Even after its capacity filter the number of valid exceptional retained assignments is

\[
2^k\left(1+k+\binom{k}{2}\right). \tag{8}
\]

There is one exceptional centre of each copy on the common line, so at most two of those `k` centres may be retained; the other `k` exceptional centres may be chosen independently. No other line contains more than one of these exceptional centres. This proves (8), not just the coarser `2^(2k)` upper bound.

The exact original-cell deletion cost is

\[
\boxed{\tau=4k-2.} \tag{9}
\]

To verify it, let `r` be the number of retained cells on the common line contributed by one copy. The global capacity implies `r in {0,1,2}`. Clearing that copy's other five lines while retaining exactly `r` of its three common-line cells costs respectively `4,3,2` deletions. For `r=0`, all three common-line cells are deleted, including its exceptional centre and one shared cell; two remaining lines through the other centre still require a fourth deletion. For `r=1`, the two required common-line deletions cannot clear all three lines through the other centre, so at least one further deletion is needed. For `r=2`, one common-line deletion cannot clear the other five lines, while deleting the two exceptional centres attains two. The upper bounds for the first two states follow by additionally deleting suitable common-line cells from that two-centre repair. Hence the exact per-copy cost is `4-r`; summing with `sum r<=2` gives (9), attained by any allocation of two retained common-line cells.

This also refutes summing the isolated-copy scalar cost two after identifying a shared line: that would give `2k` and omit the original common-line capacity. The exact line state is what permits separation without losing this global obligation.

Every occupied row and column has at most two points, but there are empty rows and columns in the bounding square. These are connected finite integer boards, not saturated square boards and not a random-board frequency construction.

### 5.2 An actual nonconvex continuation profile

The ten-point board

```
(100,100),(101,102),(102,104),(103,103),(104,102),
(100000,100000),(101000,102000),(102000,104000),
(103000,103000),(104000,102000)
```

has two triangle blocks joined at the line `d=0`. That line contains four cells, and the other four overfull lines are triples. For either triangle branch, the minimum deletion costs with exactly `r=0,1,2` of its two common-line cells retained are `(2,2,1)`. Consequently the continuation cost seen at the common line from a parent contribution `t=0,1,2` is

\[
\boxed{C(t)=(1,2,2).} \tag{10}
\]

This profile is not convex: its successive costs are one and then zero. Replacing its entries by independent linear slot costs, or assuming convexity without proof, loses an exact conditional state. The gadget in Section 3 preserves all three entries. The complete board has repair cost three, verified directly on its original cells.

### 5.3 Mandatory cell-port coverage is essential

A cell incident to two lines must use both incidences together. If one line has capacity one and the other capacity zero, the cell cannot be retained. Without the coverage condition, a positive-weight edge from its first port to the available slot can still be selected while its other port is left unmatched. This is a false half-retained original cell. The finite gadget hostile gives that edge weight ten: the unprotected matching obtains ten, while the true retained-cell optimum is zero. The dominant coverage weight forces the internal port pair and repairs the mismatch.

The inherited concurrent triple and theta remain additional controls: they test the literal shared original decision and the inseparable joint choice. No clique replacement of a three-way cell is used anywhere in the new reduction.

## 6. Scope, finite universe, and reproduction

The theorem applies to every finite original board and finite chosen direction set, with arbitrary overfull line lengths. All local capacities are residual capacities derived from the complete original constraints. If a final no-three-in-line set `A` is obtained after adaptive insertions or moves, `B intersect A` is still a feasible original retained set, so `|B setminus A|>=tau(B)`. This remains an original-cell-loss statement.

The [standalone producer](../../04-computation/continuing14_20260908_no3_line_states.py), [stored transcript](continuing14_20260908_no3_line_states.out), and [certificate](continuing14_20260908_no3_line_states_certificate.json) pass **463,582 always-active exact gates per mode**. The graph consumer is NetworkX's integer-weight blossom implementation; the certificate records version 3.5. The verifier compares it with independent literal retained-subset optimization for the degree-cost gadget and with whole-instance triple-hitting deletion for the original board. No mathematical producer is imported. The all-size theorem follows from the analytic recurrence and explicit compiler, not a finite extrapolation.

The complete actual-board universe is all 2,137 simple two-regular boards of dimensions two through five, with counts `1,6,90,2040`, and no incidence filter. The abstract bank is all 109 linear triple-incidence systems on zero through four line vertices, with counts `1,1,2,9,96`: choose shared two- or three-way cells with no repeated line pair and maximum shared line degree three, then fill lines to triples with private cells. They are not all claimed geometric realizations. Every root choice is tested for this abstract bank and for the six inherited geometric controls in three and four directions.

The standalone degree-cost bank exhausts all 125 reward triples in `{-2,-1,0,1,2}^3`, two different original-cell incidence patterns, and four exact-degree modes: 1,000 literal-subset comparisons. The separate half-cell hostile tests mandatory port coverage. Every block kernel used in the board bank is also compared with literal subset optimization, including the prescribed parent-degree cases. The block decomposition is checked by full incidence ownership, the block-cut forest identity, and independent vertex-deletion articulation tests for both node types.

The connected theta family is checked for `k=1..5`, including all original coordinates, line lengths, local and merged exceptional counts, every valid old exceptional assignment, and the exact repair formula. Every root choice is replayed for its first two members. A direct 2,048-retained-subset check on the eleven-point theta recovers the complete local common-line cost table `(4,3,2)`. The ten-point nonconvex example supplies its literal `(1,2,2)` line message. These controls preserve the distinction between complete saturated boards, unsaturated geometric families, abstract incidences, and auxiliary matching gadgets.

All gates raise explicit exceptions and remain active under optimized Python. Normal and `-O` runs have byte-identical raw LF stdout and regenerated certificates. An outside source writes its certificate beside itself; a source placed under `04-computation` writes under `05-knowledge/results`.

The next missing coordinate is now confined to genuine higher-incidence cells inside one ordinary biconnected block. This report does not assert that they can always be removed cheaply, that a feedback-cell set has a particular optimal size, or that geometric two-regular boards avoid such cores. The established improvement is exact line-state separation and its concrete exponential gain over the previous merged-bag parameter.

Independent acceptance: [referee and reproducible controls](continuing14_20260908_no3_line_states_audit.md). Audited status, referee link and the referee-authorized polynomial-overhead-per-branch clarification only; theorem, call bound and computation unchanged. The accepted producer and filed report hashes are recorded in the manifest.
