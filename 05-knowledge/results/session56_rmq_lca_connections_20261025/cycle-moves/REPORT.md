# Path-locked score fibres: sharp triangle/quadrilateral wiring

Status: PROVED elementary graph statements, with established cycle-reversal machinery credited below; FINITE-EXACT certificates at explicit bounds. No novelty claim. No LRC, dyadic-guard, native-complement-metagraph connectivity, Hamilton-count invariance, or shortest-route algorithm claim.

## Inheritance and concept board

The immediate source is the Session 56 tournament-hierarchy packet: its cut-flux identity identifies the full ordered cut-flux array with the labelled tournament score vector, and its fixed-path H=11/H=13 pair differs by a directed residual quadrilateral. This packet explains that collision by identifying the entire fibre operation. The inherited source witness is independently reconstructed here, not assumed from its numerical output.

Board: labelled degree/divergence; directed disagreement circulation; chord refinement; fixed Hamilton path/forest; Hodge cycle coordinates; native-operation and Hamilton-count guards. Anchor: exact legal moves inside a score fibre. Niche: arbitrary forest locks. Wildcard: sharp universal induced-cycle threshold.

## 1. General orientation-fibre theorem

Let G be a finite simple undirected graph on labelled vertices. A score fibre consists of all orientations D of G with one prescribed outdegree at each vertex. Vertices are not relabelled or identified. A legal cycle reversal reverses every edge of a currently directed simple cycle.

**Theorem.** Each nonempty score fibre is connected by reversals of directed cycles whose underlying cycles are induced in G. Consequently, if G has no induced cycle longer than h, cycles of length at most h suffice.

Proof. For orientations D,E in one fibre, take exactly the edges on which they disagree, oriented as in D. At each vertex the number of outgoing disagreement arcs equals the number of incoming ones: their difference is d_D^+(v)-d_E^+(v)=0. This balanced digraph decomposes into edge-disjoint directed simple cycles. Reversing these cycles one at a time transforms D to E.

It remains to implement any directed simple cycle C using induced directed cycles. If C is induced, do so directly. Otherwise choose a chord uv. The two directed paths of C go from u to v and v to u. Exactly one closes to a directed cycle with the current chord orientation. Reverse that smaller cycle first. The chord now points the other way, so the other path closes to a directed cycle; reverse this second cycle. Every edge of C was reversed once, the chord twice, and all other edges not at all. Both smaller cycles have fewer vertices. Induction terminates in induced cycles. Importantly, every intermediate orientation is legal and remains in the same score fibre. This is not merely a formal sum of circulation vectors.

**Constructive bound.** Reversing an r-cycle requires at most r-2 induced-cycle reversals. A chord split has smaller lengths r1,r2 with r1+r2=r+2; the inductive costs add to at most r-2. Thus if the disagreement set has M edges decomposed into q edge-disjoint cycles, the constructed route has at most M-2q moves. For unequal endpoints q≥1, giving the coarser M-2 bound. This is a route bound, not a geodesic formula. Chords outside the disagreement set may be temporarily reversed and restored.

**Universal sharpness.** Let h(G) be the maximum induced-cycle length, with h(G)=0 when G is a forest. For any induced k-cycle C, construct an orientation by directing C cyclically, directing all edges between outside vertices in a total order, and directing every edge from an outside vertex into C. Its only directed simple cycle is C. Any other orientation with the same labelled scores differs by a balanced subdigraph, so its disagreement edges are either empty or exactly C. The entire fibre therefore has exactly two states, related only by the k-cycle reversal. Consequently h(G) is the exact universal maximum-move-length threshold, with 0 interpreted as requiring no moves because every forest fibre is a singleton. In fact every induced-cycle support is indispensable for some fibre. The construction is not constrained to an externally specified score vector.

## 2. Forest locks, and the Hamilton-path specialization

Fix any forest F on the tournament vertex set, with the direction of each of its edges locked. The free edges form G=K_n-E(F). Prescribing full labelled tournament scores is equivalent to prescribing labelled residual scores on G, since the locked outgoing contribution is constant.

**Corollary.** Any two tournaments with the same labelled scores and the same directions on F are connected by directed triangle and quadrilateral reversals avoiding every edge of F.

Proof. If G had an induced cycle C_k, k≥5, then F on those k vertices would contain every non-cycle pair, namely k(k-3)/2≥k edges. A forest on k vertices has at most k-1 edges. Contradiction. Apply the general theorem. No assumption on the maximum degree of F is needed.

In particular, fix the directed Hamilton path P: 0→1→...→n-1. Every intermediate orientation retains this very path, not just the existence of some Hamilton path. The residual graph is K_n minus its n-1 path edges. The permitted moves are its directed induced triangles/quadrilaterals. A residual induced quadrilateral is precisely a four-vertex set on which the locked path has two disjoint edges and no other edges.

For n≤4 this residual graph is a forest, so all its labelled-score fibres are singleton. At n=5, four-cycle moves become genuinely indispensable; the next section gives a sharp example. The conclusion is all-size; the computation is only a bounded audit.

## 3. Exact sharp n=5 fibre and Hamilton-count obstruction

Take the fixed path 0→1→2→3→4. Begin with all arcs forward in the natural order.

- State A reverses (0,3), (1,4), (2,4). Its Hamilton-path count is 11.
- State B reverses (0,4), (1,3), (2,4). Its Hamilton-path count is 13.

Their full labelled scores are (3,2,1,2,2); their residual scores are (2,1,0,1,2). Their entire fixed-path, fixed-score fibre consists of these two states. There are no legal directed residual triangles in either state. The single transition reverses 0→4→1→3→0 in A. This is an induced quadrilateral in the residual graph: its missing diagonals are the locked edges 0-1 and 3-4.

A direct non-enumerative check of fibre size is short. The residual edges are 02,03,04,13,14,24. Residual score of vertex2 is zero, forcing 0→2 and 4→2. The remaining edges form the four-cycle on 0,4,1,3. Each of those vertices must have exactly one outgoing edge within this cycle, forcing one of its two cyclic orientations. Thus there are exactly two states, with no possible intermediate third state.

Therefore degree preservation and fixed-path preservation do not imply preservation of H, or even H modulo8: the values are 11 and13. The example is minimal in n for a nonsingleton fixed-path score fibre, since all n≤4 residual graphs are forests.

## 4. Incidence/Hodge interpretation and the cube boundary

Choose reference orientations for edges and let B be the signed incidence matrix with a +1 at each reference tail and -1 at its head. Represent an orientation by x in {−1,+1}^E. Then

    Bx = 2d^+ − deg_G.

Thus equality of labelled scores is exactly Bx=By, and the difference belongs to ker B. Over the reals,

    R^E = im(B^T) orthogonal-direct-sum ker B.

Scores determine the gradient part; cycle circulation is the forgotten coordinate. In K_n, BB^T=nI−J, so the gradient component is B^T(Bx)/n. Its cycle-space dimension is binom(n,2)−n+1=(n−1)(n−2)/2, equal to the number of nonpath edges after locking a spanning path.

That dimension match is not an orientation-preserving map. Relative to a spanning tree, arbitrary real chord coordinates determine the balancing tree coordinates of a circulation, but these tree coordinates need not stay in the orientation cube, and locking the tree forbids changing them at all. Indeed the correct space for differences with path arcs fixed is ker(B_G) on the residual graph, not the full K_n cycle space. For n≥4 this residual graph is connected and has cycle-space dimension

    binom(n,2)−(n−1)−n+1 = (n−1)(n−4)/2.

At n=3, the full triangle has one circulation dimension, yet the fixed-path residual graph is a single edge and every residual score fibre is singleton. The legal chord-refinement proof supplies the missing cube-feasibility argument; a linear dimension count cannot replace it.

The new moves preserve every property determined solely by the labelled score vector and the locked arcs. They are equivariant under graph relabelling and under reversal of all residual arcs: the latter maps score d to deg_G−d and sends each directed cycle to its reverse. Full tournament converse additionally reverses the locked path. These operation-level facts do not assert that a prescribed mirror-fixed locus is connected by mirror-compatible moves, nor that the native complement metagraph or extra dyadic guards admit the same routes. Mirror-coupled cycle lifting remains a separate obligation.

## 5. Literature scope

The cycle-reversal mechanism is established, not claimed new. Primary source checked 2026-10-25: H.W. Willie Wong and E.G. Tay, *Arc reversals of cycles in orientations of G vertex-multiplications*, arXiv:2108.12605v1, https://arxiv.org/html/2108.12605v1 . Theorem3.1 uses the balanced directed difference graph; the proof of Theorem3.4 uses chord splitting. Its automorphism-based score-list formulation specializes here to the identity map and labelled scores. We give the full elementary argument to make the lock and intermediate-legality constraints explicit.

For standard alpha-orientation terminology and the distinction between existence of a route and optimizing its length: Aichholzer et al., *Flip distances between graph orientations*, https://arxiv.org/abs/1902.06103 . No complexity result is imported for our specific residual-graph family. No identification with a generic toric Markov basis is needed or asserted: the orientation cube and legally directed reversals are the relevant constraints.

## 6. Reproducibility and finite universes

Run `python3 check.py` in this directory. Standard library only; output `certificate.json` and stdout `run.out` are reproducible. Bits enumerate residual edges lexicographically, with bit1 meaning lower label points to higher label.

1. Exhaust all fixed-path orientations for n=1..7. Fibre counts are 1,1,2,8,52,478,5706; orientation counts are 1,1,2,8,64,1024,32768.
2. Independently union all states along every legal induced-cycle reversal and check each labelled-score fibre is connected. Every induced residual cycle is of length3 or4.
3. For every nonrepresentative state, construct and replay a route to its fibre's minimum-mask representative. Replay checks directedness, inducedness, fixed residual support, scores at every step, endpoint, and the M−2q bound. The maximum constructed route lengths for n5,n6,n7 are 2,4,7; these are not diameter claims. Longest examples are retained in JSON.
4. Exhaust every simple labelled graph for n=1..5 and every orientation: 1099 graphs and 59809 orientation states in total. Check all score fibres connected by induced-cycle reversals. For every cyclic graph, construct the two-state longest-induced-cycle witness and verify its exact fibre. The induced C5 graph is a hostile control against an unjustified universal 3-/4-cycle claim.
5. Reconstruct the n5 sharp pair, exhaust its entire fixed-score fibre, check no directed triangle moves, and independently enumerate all vertex permutations for its H counts.

This packet makes no remote edits and does not allocate canon identifiers; the parent session owns auditing and publication.
