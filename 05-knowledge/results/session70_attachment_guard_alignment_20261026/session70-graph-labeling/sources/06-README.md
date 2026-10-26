# Square-sum edge lifetimes, endpoint markings, and parity reconstruction

Session40, 2026-10-23. **PROVED** elementary all-size statements; **FINITE-EXACT** complete square-sum path census for `1 <= n <= 32`. No external priority claim. The ordinary path counts and the endpoint-marked two-path-cover mechanism are inherited; the present additions are the full edge/pair-use ledger, temporal support changes, exact endpoint-conditioned minimum sensor banks through32, and their parity-decoding interpretation.

## 1. What was inherited, and what must not be conflated

The graph `Q_n` has vertices `1,...,n`; distinct vertices `a,b` are adjacent iff `a+b` is a square. Paths below are unoriented Hamiltonian paths, counted once up to reversal.

Closest mechanisms and corrections:

- `05-knowledge/results/session32_square_monster_20261023/square-sum/README.md`: endpoint-aware chain kernels, exact marked-ear/two-path-cover correspondence, all paths through25.
- `05-knowledge/results/collatz_mod6_20260922_w6_square_sum_hamiltonicity.md`: earlier counts through36 and cut obstructions. Its obsolete rule that every graph-degree-two vertex forces both edges was explicitly retracted. A degree-two vertex may be a path endpoint.
- `05-knowledge/results/decoder_prime_square_20260925.md`, Sections4–5: repaired endpoint-aware proof for24 and the compulsory ear `11–25–24`.
- `05-knowledge/results/collatz_mod6_20260922_w6_summand_closure_square_filter.md`: the square-sum founders are `1,2,4`; the missing strict-summand closure elements are `1,4,6`. These are different objects, not three interchangeable components.
- `05-knowledge/results/arithmetic_braids_20260917_summand.md`, Section8: corrected synchronous closure and directional restoration. Adding4 can generate6 but cannot generate1.
- `05-knowledge/results/square_sum_affine_lift_protected_edge_obstruction_20261019.md`: an existing affine-lift architecture cannot preserve an arbitrary specified edge uniformly across its whole output residue range.

The concept board is: vertex birth / path endpoint ownership / paired endpoint-component incidence / edge support over time / parity boundary / reconstruction bank. Counts alone do not retain the middle two coordinates.

## 2. The user's three colors have an exact birth interpretation

Read the supplied `13:,12` as `13–12`. The separate literal `3–12` is not an edge: `3+12=15`, not a square. No correction of that literal to another edge is assumed.

For an edge `a–b`, with `a<b`, its birth is the arrival of vertex `b`. Through25 the supplied first class is exactly the collection of edges born at a traceable order:

- birth15: `1–15,10–15`;
- birth16: `9–16`;
- birth17: `8–17`;
- birth23: `2–23,13–23`;
- birth25: `11–25,24–25`.

The supplied second class is exactly the edges born at nontraceable vertices with two new neighbors: births `13,14,19,20,21,22,24`. The third class is the remainder. It includes the earlier scaffolding and the failed leaf arrival `7–18`.

This is a legitimate three-way edge partition, but **not a proper edge coloring** (two edges incident at23 have the same first color), and **not the forced/optional/unused partition**. At25, first-color `1–15` is unused, first-color `13–23` is optional, and first-color `11–25` is forced. All three usage states already occur in that one color.

The arrival law is inherited:

`d_new(n) = floor(sqrt(2n-1)) - floor(sqrt(n))`.

Components grow to two at2 and three at4, then stay at three through12. The actual founders are1,2,4. The components at12 are the paths

- `2–7–9` (three vertices),
- `12–4–5–11` (four vertices),
- `10–6–3–1–8` (five vertices).

Vertex13 merges two components;14 merges the last two. The graph stays connected for every larger order. There is no later recurrent three-component phase in this graph. The directed strict summand shadow is instead `x -> z` iff `x<z` and `z != 2x`; its faithful two-parent generation carrier is a hypergraph. The square-sum graph retains the two parents of a square-target hyperedge, not its directed target.

## 3. Exact edge usage is nonmonotone

For an existing edge `e`, define `h_n(e)` as the number of Hamiltonian paths of `Q_n` containing it, and `H_n` as the total. At orders with `H_n>0`:

- UNUSED: `h_n(e)=0`;
- FORCED: `h_n(e)=H_n`;
- OPTIONAL: `0<h_n(e)<H_n`.

At orders with no Hamiltonian path we use **NO_PATH**, not “forced”: otherwise the universal quantifier would misleadingly become vacuous. “Used at every observed admissible size” means positive use at every tested order having a path and containing the edge. It is weaker than “in every path,” and neither means forever.

The complete counts at positive orders through32 are:

`n: 1,15,16,17,23,25,26,27,28,29,30,31,32`

`H: 1, 1, 1, 1, 3,10,12,35,52,19,20,349,392`.

The first eight order comparisons requested in the prompt therefore check exactly:16,17,23,25 are possible;18,19,22,24 are impossible.

Important edge histories:

- `1–3`: unused at15,16,17; forced at23; optional at25; unused again at26; optional at27–32.
- `1–15`: forced at15,16,17; unused at23 and25; reactivates at26 and remains optional through32.
- `5–20`: forced at23,25,26,27,28; unused at29,30; reactivates at31.
- `7–9`: forced at15–17; optional at23,25–27; unused at28–31; reactivates at32, in exactly2 of392 paths.
- `9–16`: forced at16,17,23,25,26; optional at27; forced again at28–31; optional at32, in390 of392 paths.

At25, **only `1–15` is unused**;14 edges are forced and17 optional. At32, **all46 edges are optional**. Thus every edge present by25 has been used by32, and no edge present at32 is forced at32. This proves neither permanent future activity nor any all-order prescribed-edge theorem.

A readable full horizon table is in [EDGE_HISTORY.md](EDGE_HISTORY.md); the exact counts, all paths, all joint edge counts, and activation changes are in [receipt.json](receipt.json).

### Structural certificate for the dead edge at25

The receipt contains a path-free contradiction certificate for the supposition `1–15` is used. Since18 is a leaf, fix candidate endpoints `{18,e}` for every `e != 18`, all24 possibilities. Each case applies only these sound rules:

1. the selected degree is at most1 at a marked endpoint and at most2 elsewhere;
2. the available degree must reach that exact target;
3. reaching the target with selected edges deletes all other incident edges;
4. having exactly the target number of available edges forces them;
5. a forced cycle is impossible in a Hamiltonian path.

All24 cases close without branching beyond the explicitly marked endpoint. The small recurring obstruction is the potential forced cycle `1–8–17–19–6–10–15–1`; exceptional endpoints break that cycle but trigger a degree deficit or another forced cycle. Every force retains the endpoint marks. The script records each deletion/force and final contradiction, rather than silently applying the retracted degree-two rule.

## 4. Actual endpoint markings work through25, then need shared incidences

**FINITE-EXACT.** Every unordered endpoint pair supports at most one path through order25. At25 there are ten endpoint pairs, each containing18, with the other endpoint in

`{2,3,4,8,9,10,11,13,22,23}`.

This uniqueness was already visible in the inherited Session32 census. The next order is a sharp failure: at26 there are12 paths but only9 endpoint pairs. The pairs `{18,13}`, `{18,22}`, `{18,26}` each support two paths. The single incidence bit

`“does this path contain edge3–6?”`

separates all three surviving ambiguities. Therefore one extra edge-membership bit is sufficient and necessary at26.

The finite generalization is exact. Supply the actual unordered endpoint pair and a fixed global bank `B` of edges. Observe one membership bit for each edge in `B`. The minimum bank sizes are:

| order | paths | endpoint pairs | largest endpoint fibre | minimum extra bits |
|---|---:|---:|---:|---:|
|25|10|10|1|0|
|26|12|9|2|1|
|27|35|16|4|3|
|28|52|19|5|4|
|29|19|15|2|1|
|30|20|15|2|1|
|31|349|172|8|5|
|32|392|220|9|6|

“Minimum” here means **within this specific bank of individual edge-membership tests, fixed globally across all endpoint pairs**. It is not minimum arbitrary information, not an adaptive query lower bound, and not a theorem about all tournament invariants.

The checker gives exact upper banks and an exhaustive lower-bound search. Independent audit enumerated all smaller bank combinations after merging identical test columns. For32 it rejected all169,911 five-column combinations among31 distinct nonzero tests; for31 it rejected all8,855 four-column combinations among23. Smaller bank sizes fail automatically because any smaller successful bank can be enlarged.

### The exact ambiguity hypergraph

Fix any graph and its Hamiltonian paths. For every distinct pair of paths with the same endpoint pair, take their symmetric difference of edge sets as a hyperedge. A bank reconstructs all paths iff it meets every such hyperedge. Thus the minimum bank is exactly the transversal number of this **realized path-difference hypergraph**.

Proof: two membership vectors coincide exactly when the bank misses their symmetric difference. This is an all-size equivalence, but determining the hypergraph may require knowing the paths. It is a precise reconstruction specification, not an efficient all-size classification algorithm.

Shared incidences matter even without complete recovery. At25, edges `1–24` and `4–5` occur individually in4 and3 paths, yet occur together in **zero** paths. They are disjoint edges, so their incompatibility is not simply a degree-three clash at a shared vertex. Marginal edge frequencies do not retain joint realizability.

## 5. An all-size parity decoder, with its validity gate

Let `G` be any connected graph on at least two vertices, choose a spanning tree `T`, and retain actual distinct endpoints `s,t`. A Hamiltonian path's edge indicator `x` satisfies the binary boundary equation

`B x = 1_s + 1_t` over `F_2`,

where `B` is the vertex-edge incidence matrix. Retain all membership bits on the cotree `E(G) minus E(T)`.

**Theorem (PROVED, standard cycle-space mechanism).** Endpoints and cotree bits uniquely determine the path edge set. Therefore at most `|E|-|V|+1` such bits always suffice.

**Proof.** If two path indicators agree on the cotree, their difference is supported on the tree and has zero boundary. A nonempty forest has a vertex of degree1, so it cannot be an even subgraph. Their difference is empty. Constructively, subtract the cotree boundary from the endpoint boundary and peel tree leaves: the unique tree edge at a leaf is selected iff its remaining boundary charge is1. Update the parent charge and continue. The total boundary is even, so the last root equation holds.

This also gives an iff statement for the larger affine space of **all binary edge sets with the specified boundary**: observing a bank is injective on that space iff the unobserved subgraph is a forest. A cycle in the unobserved subgraph can be toggled invisibly; a forest cannot. Its minimum coordinate bank is exactly the cycle rank. This larger-model hostile need not itself be a pair of Hamiltonian paths, so it is not a lower bound for the restricted path bank. For example, the cycle rank of `Q_32` is15, whereas the exact restricted minimum above is6.

**Validity gate.** Binary boundary controls parity, not degree or connectedness. The decoded candidate must still have degree1 exactly at the two endpoints, degree2 elsewhere, and be connected. These conditions are necessary and sufficient for a Hamiltonian path. A boundary-correct short path omitting other vertices is an immediate hostile. The decoder authenticates a valid candidate; it does not assert that every chosen cotree bit string is valid.

**Tournament transfer.** The same argument applies to the undirected support of the arcs of a fixed, labelled tournament. Actual path endpoints and cotree arc-use bits reconstruct the path's support; then test that all its arcs follow one directed ordering. This does not solve the quotient metagraph by itself: vertex alignment, the original tournament adjacency, converse transport, and the nonpath-arc complement operation still have to be retained. Numerical class signatures do not supply these labels automatically.

The rigorous electrical analogy is Kirchhoff's parity boundary: the two marked endpoints are the odd sources and cotree bits carry the cycle freedom. No physical electromagnetic or topological field-theory identification is needed.

## 6. Exact edge-marked birth polynomial

This extends the inherited marked-ear correspondence to every arrival degree and to edge-use and shared-incidence statistics.

Let a new vertex `v` be added to a graph `G`, with neighbor set `N`. Give every edge a commuting weight `x_e`. Let `H_G` be the sum of edge-monomials of unoriented Hamiltonian paths. Let `E_a(G)` be the analogous polynomial restricted to paths having endpoint `a`. Let `F_ab(G)` count spanning forests consisting of exactly two nonempty paths, with `a` and `b` endpoints in **different components**. Each undirected forest is counted once; singleton components are allowed and their sole vertex counts as an endpoint without a factor2.

Then, for nonempty `G`,

`H_(G+v) = sum_(a in N) x_(va) E_a(G) + sum_({a,b} subset N) x_(va) x_(vb) F_ab(G)`.

For the single-vertex base graph, its one path contributes once to `E_a`; this is the only convention needed to avoid a double endpoint count. The empty graph is handled separately.

**Proof.** A Hamiltonian path gives the new vertex degree1 or2. Deleting a degree-one new vertex leaves an old Hamiltonian path with the required exposed endpoint. Deleting a degree-two new vertex leaves precisely two spanning path components with distinct joining endpoints `a,b`. Conversely each displayed object uniquely reconstructs the path. Joining endpoints of the **same** component would close a cycle and omit the other component, which is why shared component ownership cannot be discarded.

The polynomial is multiaffine. Evaluation at all weights1 gives the count recurrence. Its first derivative in an old edge gives that edge's participation count; a mixed derivative in two old edges gives their joint participation. Differentiating in a new edge separates attachment and merge contributions. Thus the exact same carrier retains endpoint, edge and shared-incidence constraints rather than approximating them by separate marginals.

Selected exact arrival decompositions:

|new order|endpoint-extension terms|two-path merge terms|total|
|15|0|`F_(1,10)=1`|1|
|16|`E_9=1`|0|1|
|17|`E_8=1`|0|1|
|23|0|`F_(2,13)=3`|3|
|24|0|`F_(1,12)=0`|0|
|25|0|`F_(11,24)=10`|10|
|26|`E_10+E_23=1+1`|`F_(10,23)=10`|12|
|27|`E_9+E_22=1+2`|`F_(9,22)=32`|35|
|28|`E_8+E_21=1+4`|`F_(8,21)=47`|52|
|29|0|`F_(7,20)=19`|19|
|30|`E_6+E_19=1+1`|`F_(6,19)=18`|20|
|31|`E_18=20`|`F_(5,18)=329`|349|
|32|`E_4+E_17=35+24`|`F_(4,17)=333`|392|

Here the old graph is always `Q_(n-1)`. The script independently verifies every extension term against the preceding order's endpoint counts. The merge terms are identified exactly by the deletion bijection; they are not represented as an independent enumeration of all two-path forests.

The key loss is visible at25: `Q_24` has no Hamiltonian paths, so all `E_a` vanish, yet `F_(11,24)=10`. Any recurrence retaining only Hamiltonian path totals or exposed endpoints of Hamiltonian paths loses all ten new solutions. A sufficient growing state must retain the two-component joining data (and, for continued growth, the corresponding higher-component forest boundary states).

## 7. All-size boundary and next useful obligation

The strongest proved general reconstruction statement here is the endpoint-plus-cotree decoder, together with degree/connectivity/orientation validation. It is not a small all-size metagraph wiring theorem. The exact minimum bank computations are finite, and their candidate family is explicit.

Unrestricted square-sum interval Hamiltonicity is already established in the inherited literature route. It does not imply every specified edge occurs eventually or at every sufficiently large order. The repository's protected-edge obstruction is concrete: the edge `482–3362`, whose sum is `62^2`, cannot be preserved throughout a full49-residue output range by one common parity-compatible input-pair, intact-block affine-lift architecture. This is an obstruction to that architecture, not a proof the edge is absent from any particular large graph.

The existing absorber theorem in `square_sum_prescribed_edge_absorbers_20261019.md` gives a positive conditional route: if actual host edges in a backbone satisfy Hall's matching condition, disjoint omitted blocks can be inserted while protected edges and endpoints remain intact. Numerical host existence is not host membership or Hall expansion.

A useful next obligation is therefore to carry the **realized difference hypergraph and marked forest attachment states** through an explicit interval construction, rather than extrapolate the finite support table. For tournaments, carry those same kinds of actual owner/endpoint incidences through aligned converse and complement charts, retaining the independent mod8 guard. Shared terminology alone does not provide that transport.

## 8. Reproduce and audit

Run in this directory:

`python3 check_square_incidence.py > receipt.json`

`python3 -O check_square_incidence.py > receipt.optimized.json`

The two receipts are byte-identical. The primary search is direct vertex-order DFS with only sound leaf/isolated-vertex reductions. It stores every path through32, every edge count, every pair count, exact endpoint fibres, minimum-bank search traces, birth terms, and the24 unused-edge contradiction traces. Edge and pair mass identities are checked. The parity decoder is tested on every nontrivial enumerated path.

Independent Session40 audit uses degree-constrained acyclic edge-subset enumeration, not path DFS, through25, and independent brute combinations for sensor minimality through32. Its artifacts are in the session audit packet. The result is finite exactness and all-size elementary proofs, not an unbounded experimental extrapolation.
