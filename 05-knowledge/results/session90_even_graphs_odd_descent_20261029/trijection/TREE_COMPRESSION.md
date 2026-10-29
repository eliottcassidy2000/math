# A natural tree-edge compression behind the marked trijections

**Session90, 2026-10-29. Status: PROVED elementary; finite verification separately stated.**
This is a generalization of the repository's marked path/star constructions, not an assertion of external priority. It is a correspondence of **marked completions**, not a bijection from all unmarked tournament classes of order `n` to those of order `n−1`.

## Inheritance and distinctions

Sources are read at repository commit `3eb49fa63d7ff4e58742b1ae84944f832e87a118`. Relevant paths:

- `05-knowledge/results/trijection_four_nodes_three_edges_20261020.md`, Theorem1: a tournament with a marked directed Hamiltonian path, directed cycle-plus-solo, or out/in-star becomes a smaller tournament modulo respectively `1`, `C_(n−1)`, or `S_(n−1)`. The two star directions give separate sheets. These are three marked species, not three equal unmarked class counts.
- `05-knowledge/results/tournament_graph_trijections_20261020/TRIJECTIONS_AND_TREE_GAUGES.md`, Theorems1–2: fixed oriented spanning trees are cut-switching transversals; arbitrary constraint graphs have cycle consistency conditions and component multiplicities. Its observer, cut-quotient, and Hamiltonian-incidence maps are distinct.
- `05-knowledge/results/session66_root27_aliquot_graph_interfaces_20261026/session66-graph-pairs/REPORT.md`: joint graph layers retain shared matching; independent graph classes can forget it. The minimal converse hostile is `source ⊕ C3` versus its converse.

The live conceptual board is: tree paths; orientation torsors; mark stabilizers; cycle/component budgets; enlarged faithful carriers. Closest mechanism is the fixed-path staircase codec. The underused coordinate is the **edge set of the tree as the new vertex carrier**, rather than selecting an arbitrary vertex to delete.

## 1. Outer-endpoint theorem

Let `F` be any finite undirected tree on `V`, with `|V|=n≥2` and edge set `E`, so `|E|=n−1`. Write

- `D_V = {(u,v): u≠v, {u,v}∉E}` for ordered nonadjacent vertex pairs;
- `D_E = {(e,f): e,f∈E, e≠f}` for ordered distinct tree-edge pairs.

Define `Φ_F(u,v)` to be the first and last edges of the unique `u`–`v` path, in that order.

**Theorem 1.** `Φ_F:D_V→D_E` is a bijection. It commutes with reversal and with every tree isomorphism. Consequently the unordered nonedges of `F` are naturally identified with the unordered pairs of its edges, including the orientation-reversal information.

**Proof.** The path has length at least two, so its first and last edges differ. Conversely, given distinct edges `(e,f)`, delete `e`. Exactly one of its two endpoint components contains `f`; let `u` be the endpoint in the other component. Similarly let `v` be the endpoint of `f` on the side not containing `e`. The unique `u`–`v` path begins with `e` and ends with `f`. These choices are unique and invert the forward map. Swapping the edge order swaps `(u,v)`. A tree isomorphism preserves unique paths, deletions and their components, proving naturality. ∎

An equivalent inverse chooses the unique pair of endpoints, one from each edge, whose distance is maximal. The supplied producer uses this distance definition, while the proof and independent audit use deletion components/endpoint paths.

## 2. Compression of every oriented-tree completion

Fix an orientation `X` of all edges of `F`. Let `Comp(X)` consist of tournaments on `V` that extend these oriented edges. Given `T∈Comp(X)`, define a tournament `C_X(T)` on the **tree edges** by

`e→f in C_X(T)` exactly when `u→v in T`, where `Φ_F(u,v)=(e,f)`.

**Theorem 2.** `C_X` is a bijection

`Comp(X) ↔ Tour(E)`.

It is equivariant for `Aut(X)` acting on vertices on the left and on tree edges on the right. For every completion,

`Aut(T,X) = Stab_(Aut(X))(C_X(T))`

under the induced action. In particular, isomorphism classes of tournaments carrying a mark isomorphic to `X` are exactly `Tour(E)/Aut(X)`.

**Proof.** All arcs on `E` are already prescribed by `X`. Theorem1 assigns every other arc, with its direction, to precisely one arc on the new edge carrier. This has an explicit inverse: install `X`, then use the outer endpoints to restore every remaining arc. Naturality proves equivariance and equality of stabilizers. Normalizing any mark isomorphic to `X` differs only by an element of `Aut(X)`, which gives the final orbit statement. ∎

More globally this is a natural isomorphism of tree-marked species:

`(F, tournament T on V(F)) ↔ (oriented tree X=T|F, tournament on E(F))`.

Every tree isomorphism, including one to a differently labelled tree, transports the edge carrier and commutes with compression and expansion. The fixed-`Aut(X)` quotient is one fibre of this full natural correspondence, not its entire naturality statement.

**No orientation assumption on the tree is needed:** it need not be an out-tree, in-tree, directed path, or have a selected root. The theorem handles every one of its `2^(n−1)` orientations. At `n=2`, the completion and zero-edge tournament are both unique; the oriented edge's automorphism group is trivial. At `n=1`, the trivial analogous statement holds with an empty new carrier, although the endpoint construction has no inputs.

### Exact enumeration formula

For `g∈Aut(X)`, let `o_2(g)` be the number of orbits of its induced action on unordered pairs of tree edges. The action on tree edges is faithful for an oriented tree with at least two vertices. Thus

`|Comp(X)/Aut(X)| = (1/|Aut(X)|) Σ_(g of odd order) 2^o_2(g)`.

An even-order permutation of the edge carrier has a power swapping two of its elements and fixes no tournament. An odd-order permutation reverses no pair and permits one independent choice per unordered-pair orbit. This is the inherited tournament Burnside mechanism on a newly identified carrier, not a new Burnside formula.

`Aut(X)` itself need **not** be odd: an outward star has leaf transpositions. Only the stabilizer of a completed marked tournament is necessarily odd. Confusing these two groups would discard precisely the even elements whose fixed-point sets vanish.

### Recovery of the old special cases

For a path with vertices `v_0,…,v_(n−1)` and edges `e_i={v_i,v_(i+1)}`, the chord `{v_i,v_j}`, `j≥i+2`, maps to `{e_i,e_(j−1)}`. This is the old staircase shift, now derived from paths rather than a coordinate trick.

For an outward star, each tree edge identifies its leaf. The chord between two leaves maps to their two spokes. Hence the smaller tournament is literally the leaf tournament and the frame group is `S_(n−1)`. The inward-star sheet is separate.

A cycle-plus-isolated-vertex frame is not a tree. Its earlier compression uses the directed cycle's successor map to move the solo vertex's bits onto cycle edges. The common abstraction is an equivariant correspondence of free **oriented pair slots**; the concrete map is different. The tree theorem does not erase the triangle-plus-isolate obstruction in the cut-switching chart.

## 3. Why this is not switching normalization

The new map preserves every completion arc once `X` is retained, and preserves the marked automorphism action. It does not quotient by cuts and does not claim that the smaller tournament has the original Hamiltonian-path count, strong components, or score sequence.

A minimal example is the path `0→1→2`. Its two completions are the transitive tournament and the directed triangle, with Hamiltonian-path counts `1` and `3`. Their compressed tournaments each have two vertices and Hamiltonian-path count `1`. The missing statistic is not a lost source bit: both smaller labelled orientations are distinct, and expansion with `X` reconstructs the original. Rather, Hamiltonian counting does not commute with this change of carrier without the decoding rule.

Likewise, dropping `X` is not harmless. The smaller tournament alone does not specify which tree, which orientation of its edges, or which edge action was used. A quotient can be economical for a fixed mark while requiring the whole frame for unrestricted reconstruction.

## 4. A normal form for arbitrary `n−1`-edge constraint graphs

Let `F` now be any spanning graph on `n` vertices, with `m=n−1` edges and `c` components. Its cycle rank is `β=m−n+c=c−1`.

Choose a spanning forest `R⊆F`, and add `c−1` connecting edges `B` between its components, producing a tree `K=R∪B`. Any tournament extending a prescribed orientation of `F` is equivalent to:

1. arbitrary orientations of the `c−1` new bridge edges in `B`;
2. a tournament on `E(K)` under Theorem2;
3. `β=c−1` fixed arc constraints on that edge tournament, one for each removed edge of `F\R`, transported by the outer-endpoint map.

The coordinates are distinct, so the number of states remains `2^binom(n−1,2)`. This is an exact labelled normal form after the forest and bridges are chosen. It is **not** an unmarked equivariant normal form under `Aut(F)` unless these choices and their transporters are retained.

This makes the global balance visible: each disconnected component beyond the first supplies one freely orientable bridge, while each independent cycle supplies one imposed chord condition. Equal counts of the two explain the bit count, not equality of native operations or stabilizers.

For cut-switching normalization the inherited result is sharper: coverage is `2^(−β)`, with `2^(c−1)` representatives in each admissible switching class. A spanning tree is precisely the case where both debts vanish. At order four, `K3 ⊔ K1` has one of each; it is not interchangeable with `P4` or `K1,3` as a complete switching chart.

## 5. Forgetting the tree mark: the exact counting boundary

Summing the marked orbit counts over isomorphism types of oriented trees counts pairs `(T,F)` with an underlying spanning tree, not ordinary tournament classes. The fibre over `[T]` is the number of `Aut(T)`-orbits on spanning trees of `K_n`.

Already at order3, the transitive tournament has three such orbits and the directed triangle has one. The marked total is4, whereas `a(3)=2`. At order4 the marked total is44, whereas `a(4)=4`. The underlying tree has its orientation induced by `T`; it is not an independent extra choice once the pair is fixed.

The weighted groupoid mass has a clean identity, for `n≥2`:

`Σ_[(T,F)] 1/|Aut(T,F)| = n^(n−2) · 2^binom(n,2) / n!`.

This follows by counting all labelled tournaments and all labelled spanning trees, then dividing by `n!`. It does not imply that all summands have equal weight or that forgetting the mark has a constant fibre. The verifier checks the marked totals `1,4,44` and masses `1,4,128/3` at orders2,3,4.

A full A000568 recurrence from these coordinates still owes the mark-forgetting incidence or an exact canonical-mark rule. That obligation is the analogue of retaining alignment in the other graph/tournament correspondences.

## 6. Reproduction and claims boundary

Run:

```
python3 check_tree_compression.py
python3 -O check_tree_compression.py
```

The normal run has **54,552 active checks**. It includes:

- every labelled tree through six vertices, by Prüfer words: `1,3,16,125,1296` at sizes2–6;
- the outer-endpoint inverse and full pair-slot coverage for each;
- completion round trips `1,6,128,8000` at sizes2–5, with a varying prescribed tree orientation;
- all vertex relabelings for every tree through size5, and five deterministic random relabelings per size6 tree;
- every oriented labelled tree through size4, comparing exact marked stabilizers and literal Burnside orbit counts.

The sum `880` printed for size4 is a sum over **labelled oriented tree frames** of their respective orbit counts, not A000568, not a count of unmarked tournament classes, and not a new named sequence claim.

The theorem is an exact frame-based classification and counting transfer. It supplies no asymptotic improvement for the unmarked prime atlas without an additional theorem controlling which marks cover each tournament, how they overlap, and how to choose or certify representatives efficiently. Those are separate algorithmic obligations.
