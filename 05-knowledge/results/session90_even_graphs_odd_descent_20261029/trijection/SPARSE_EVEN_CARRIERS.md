# Faithful tournament carriers with one fewer edge than vertices

**Session90, 2026-10-29. PROVED elementary constructions; FINITE-EXACT tests below.**
No external priority claim is made for these gadget reductions. The contribution here is an explicit, strictly decodable construction connecting the repository's graph-pair, trijection, odd-automorphism, and attachment-port themes with the exact scope stated below.

## 1. The meaning of the edge budget

For a finite graph with `v` vertices, `e` edges and `c` connected components, the first Betti number is

`β = e−v+c`.

Thus `e=v−1` means `β=c−1`, or Euler characteristic1. It does **not** mean a tree unless connectedness is also required. This identity was already emphasized in the October20 trijection reset and is inherited here.

The new question is sharper: can a graph in this apparently sparse family retain a tournament **and its complete automorphism group** without relying on arbitrary canonical labeling? Yes, after increasing the carrier size. In fact the same global edge budget permits an exact encoding of arbitrary simple graphs, so isomorphism on this disconnected family is GI-complete. The details identify where the information lives.

This does not contradict the fixed-size capacity obstruction: a bounded tuple of `n`-vertex graphs with `n−1` edges cannot encode all tournament classes on the same `n` vertices for large `n`. The parent session's capacity note proves the quantitative bounds. Here vertex count grows explicitly.

## 2. A connected, automorphism-exact tournament graph

Let `T` be a tournament on `n≥1` vertices. For each `v`, create three core vertices `c_v,t_v,h_v`, with edges `c_v t_v` and `c_v h_v`. Attach internally disjoint pendant paths with these lengths:

- at `c_v`: `1,2,3`;
- at `t_v`: `1,2,4`;
- at `h_v`: `1,2,5`.

For every tournament arc `u→v`, add the undirected port edge `{t_u,h_v}`. Call the resulting uncoloured simple graph `Q(T)`.

Each original vertex contributes24 graph vertices and23 internal tree edges. Therefore

`v(Q(T)) = 24n`,

`e(Q(T)) = 23n + binom(n,2)`,

`β(Q(T)) = binom(n−1,2)`.

The graph is connected: each gadget is connected, and every pair of gadgets has one port edge.

**Theorem 1.** Every graph isomorphism `Q(T)→Q(T')` is induced uniquely by a tournament isomorphism `T→T'`. In particular

`Aut(Q(T)) ≅ Aut(T)`.

**Proof.** Exactly the core vertices have degree at least4. Their pendant paths are therefore recognizable as the components outside the induced core, each attached at one end. Their three distinct lengths identify the roles `c,t,h`; none can be permuted internally. Each center is adjacent in the core to precisely its own tail and head, recovering the matching. The remaining core edges say exactly which `u→v` arcs occur. Consequently any graph isomorphism uniquely induces and is induced by a tournament isomorphism. ∎

This repairs the even leaf-permutation kernel in the Session66 `9n`-vertex graph codec. That earlier graph was faithful on classes but had automorphism kernel `(S2×S3)^n`. Here unequal pendant paths remove that kernel. The price is a larger but still linear connected codec before padding.

### Homotopy loses exactly the port information

Contract the entire 24-vertex gadget tree for each original vertex. The result is `K_n`, irrespective of the orientation of `T`. Contracting trees in a graph is a homotopy equivalence, so every `Q(T)` of this fixed order is homotopy equivalent to a wedge of `binom(n−1,2)` circles.

Thus the graph is fully faithful before contraction, but the contracted graph and the homotopy type forget the tournament completely. The forgotten information is which tail/head port each connecting edge used. The coincidence of the Betti number with the fixed-path completion bit count does not identify those two coordinate systems or preserve the full attachment operations.

At order3, the core for a transitive tournament has a unique cycle of length5, while the core for the directed triangle has a unique cycle of length9. Both contract to `K3` and both have the homotopy type of a circle. This is a literal small hostile to replacing port attachment data by ordinary homology or cycle rank.

## 3. Rigid tree padding with logarithmic size per type

Let `k≥1`, `ℓ=ceil(log2 k)`, and `L=ℓ+5`. For each `j=0,…,k−1`, define a tree `R_(k,j)` as follows:

1. take the path `p_0,…,p_L`;
2. attach one leaf at `p_2`;
3. for each binary bit `i=0,…,ℓ−1` of `j` that is1, attach one leaf at `p_(3+i)`.

There is no attached leaf at `p_(L−2)`. Each tree has maximum degree3 and

`v(R_(k,j)) = ℓ+7+popcount(j)`.

**Lemma.** These `k` trees are pairwise nonisomorphic and asymmetric.

**Proof.** Deleting all degree-one vertices once leaves precisely the path `p_1,…,p_(L−1)`. Any isomorphism acts on this remaining path by identity or reversal. Reversal would send `p_2`, which has an extra pendant leaf, to `p_(L−2)`, which does not. Therefore the path is fixed pointwise, and the original endpoints and all attached leaves are fixed as well. The bit positions then recover `j`. ∎

For `k=1`, this is the seven-vertex asymmetric tree with arms of lengths1,2,3, familiar from the earlier graph-pair codec. The construction extends that old rigid frame to an explicitly distinguishable family of as many components as needed.

## 4. An exact sparse Royle-even embedding

Set `k=binom(n−1,2)` and define

`S(T) = Q(T) ⊔ R_(k,0) ⊔ … ⊔ R_(k,k−1)`.

For `k=0` no trees are appended. For the empty tournament alone, define `S(∅)=K1`.

**Theorem 2.** For every finite tournament:

1. `e(S(T))=v(S(T))−1`;
2. `T≅T'` iff `S(T)≅S(T')`;
3. every isomorphism extends uniquely, so `Aut(S(T))≅Aut(T)`;
4. `S(T)` is Royle-even;
5. for `n≥1`, it has exactly `k+1` components and size
   `24n + Σ_(j<k)(ceil(log2 k)+7+popcount(j))`, with the sum empty at `k=0`.

In particular its size is `O(n² log(n+1))`.

**Proof.** `Q(T)` has `e−v=k−1`; each appended tree has `e−v=−1`, so the total is−1. The component `Q(T)` is uniquely identified by containing vertices of degree at least4, whereas padding has maximum degree3. The padding components are pairwise distinct and asymmetric, so they contribute no automorphism kernel or component permutations. Apply Theorem1.

For the Royle-even statement, use the correct character. Given any temporary orientation of a graph's present edges, its endpoint-reversal character is

`χ_G(g)=(-1)^(number of present edges whose transported orientation is reversed by g)`.

For `g∈Aut(G)` this is an orientation-independent homomorphism to `{±1}`. A graph is Royle-even when this character is identically+1. It is **not** the sign of the permutation of the edge set: on `K2`, the vertex transposition fixes its only edge but has endpoint character−1.

The automorphism group of a tournament has odd order: an involution would interchange some pair and reverse its arc. Thus every homomorphism `Aut(T)→{±1}` is trivial, proving Royle-evenness. Alternatively, orient each internal rooted pendant path away from its root, orient center–port edges from the center, and orient all port edges tail-to-head. Every tournament-induced graph automorphism preserves these orientations, giving the character directly. Padding is fixed pointwise. ∎

Small sizes for `n=1,…,5` are `24,48,79,125,187`. The first two outputs are trees. Starting at `n=3`, the one cyclic component is balanced by enough disconnected tree components. All fixed-order outputs have the same component count and homotopy type, yet retain every tournament isomorphism distinction and its complete automorphism group.

The equinumerosity between tournament and Royle-even graph classes is a separate, same-order counting statement. `S` is an explicit embedding into a larger sparse subfamily, not that equinumerosity bijection, and it does not claim that every sparse Royle-even graph is an image.

## 5. The global edge condition still permits full graph-isomorphism complexity

For any simple graph `G` on `n≥1` vertices with `m` edges, form a connected tagged cone `A(G)`:

- add a hub `r` adjacent to every original vertex;
- attach pendant paths of lengths `1,2,3,4` to each original vertex;
- attach pendant paths of lengths `1,2,3,5` to the hub;
- retain all original edges of `G`.

All cone vertices have degree at least4, while pendant vertices have degree at most2. The path-length tags recover original vertices and distinguish the hub. Hence every isomorphism is the unique extension of a graph isomorphism of the original `G`, and

`Aut(A(G))≅Aut(G)`.

The counts are

`v(A(G))=11n+12`, `e(A(G))=11n+11+m`, `β(A(G))=m`.

Append the `m` rigid trees `R_(m,j)`, obtaining `B(G)`. Define `B(∅)=K1` separately.

**Theorem 3.** `B` is an explicit polynomial-time isomorphism reduction from arbitrary simple graphs to simple graphs satisfying `e=v−1`, with a strict image decoder and exact automorphism-group preservation. Its output size is `O(n+m log(m+1))`. Consequently isomorphism restricted to these disconnected graphs is GI-complete under polynomial many-one reductions.

**Proof.** The tagged-cone component is the unique component with high-degree core vertices. Recover it, remove its intrinsically recognized tags and hub, and recover `G`. The rigid, distinct padding again contributes no symmetry. The same Euler budget calculation gives `e=v−1`. This is a polynomial reduction from unrestricted GI, while restricted GI is a special case of GI. ∎

This is a complexity classification of the target family, not a new algorithm for GI and not an NP-completeness claim. The tournament construction alone gives a TI reduction, not GI-hardness; the generic graph construction supplies the latter.

### Exact preservation of the Royle character

The generic graph codec does more than preserve the abstract automorphism group:

`χ_(B(G))(lift(g)) = χ_G(g)` for every `g∈Aut(G)`.

Choose all pendant-path orientations away from their roots and all hub spokes outward. Their orientations are preserved by every lift. The padding is pointwise fixed. Only original edges contribute endpoint reversals, exactly as in `G`. Thus `B(G)` is Royle-even if and only if `G` is Royle-even, even when the automorphism group has even order.

The distinction matters: the empty graph on two vertices is Royle-even with automorphism group `S2`; Royle-evenness does not imply odd-order automorphism group. The sparse tournament image is the stricter odd-group subfamily.

## 6. What transfers, what is lost, and where improvement remains

- **Source/target:** tournament classes to one larger sparse Royle-even graph class; arbitrary graph classes to one larger sparse graph class.
- **Preserved:** isomorphism iff, all actual isomorphisms, automorphism groups, and in the graph codec the endpoint-reversal character.
- **Not promised:** Hamiltonian-path count of the output graph, a same-size bijection, minimum carrier size among all possible encodings, or a faster isomorphism algorithm.
- **Loss under contraction:** all orientation choices of `T`, specifically their tail/head attachments, despite unchanged cycle rank and homotopy type.
- **Decisive controls:** `K2` separates edge-set permutation sign from endpoint sign; `TT3`/`C3` separates faithful incidence from identical homotopy; repeated padding components preserve the global edge budget but create unwanted permutation symmetries and fail the strict image test.

The parent session's capacity/padding analysis shows the `O(k log k)` overhead is asymptotically sharp **within core-preserving disjoint padding with odd automorphism group**. It is not a universal lower bound on tournament encodings: another representation could alter the core, use other component structures, or exploit a different global constraint. Retaining this scope prevents a gadget-optimality result from being mistaken for an unrestricted information lower bound.

The useful unifying theory is therefore concrete: global scalar budgets can be made compatible with fully faithful incidence by adding rigid components, while operations that merge attachment ports can still erase the original object. The exact category and allowed operations determine which compression is real.

## 7. Reproduction

Run:

```
python3 check_sparse_carrier.py
python3 -O check_sparse_carrier.py
```

The verifier performs **10,364 active checks**, not optimization-disabled assertions, including **406 character-preserving lifted automorphisms**. It verifies:

- all1,100 labelled tournaments through order5 and all1,100 labelled simple graphs through order5;
- strict codec round trips, Euler budgets, component counts, and deterministically scrambled relabelings;
- an independent generic uncoloured graph-isomorphism backtracker, without calling the decoders, on every source isomorphism class through order4, for distinct-class separation and exact automorphism counts;
- all rigid padding families through32 components, with generic automorphism checks through8 and two extremal codewords per larger family;
- lifted automorphisms and exact endpoint-reversal characters through source order4;
- malformed inputs, including duplicate padding that keeps `e=v−1`, loops, and nonsymmetric adjacency;
- the order3 cycle lengths5/9.

The strict parsers check the full image, including matching/roles, absence of extra core edges, and the precise padding isomorphism multiset. This is stronger than the valid-image-only decoder contract in the earlier Session66 producer. The generic isomorphism backtracker is an audit instrument for the stated finite banks, not a claimed practical general GI solver.
