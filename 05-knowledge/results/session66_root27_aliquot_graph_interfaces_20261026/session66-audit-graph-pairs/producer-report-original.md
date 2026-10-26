# Session 66 — tournaments as two simple graphs: faithful codes and lost alignment

**Status:** Authored elementary proofs below; inherited/cited results explicitly marked; finite exact checks. Ready for independent hostile audit, not a canon promotion or novelty claim. **Source:** read-only repository `eliottcassidy2000/math` at `63f52936f8dd8592d4f2fee869748baad2e1997e`. No source-repository modification or publication by this worker.

## Findings first

- Yes: every tournament isomorphism class has a **faithful, polynomial-size encoding by a pair of simple graphs**. A particularly transparent joint pair uses `7n` vertices. An independent-unlabelled pair can use a single intrinsically decodable `9n`-vertex graph plus an empty graph.
- But **joint isomorphism and two independent graph isomorphisms are different**. Forgetting alignment destroys the `7n` code already on four tournament vertices. Splitting forward/backward arcs relative to a forgotten order fails even on three vertices.
- The useful link to the previous `{7,21}` work is the odd-cycle conflict graph `Ω(T)`, where the inherited identity is `H(T)=I(Ω(T),2)`. It retains the statistic, not the tournament. A genuinely faithful independent pair is `(Γ(T),Ω(T))`, where `Γ` is the reconstruction graph below. This second coordinate is not polynomial-size in general.
- `{7,21}` are forbidden **Hamiltonian-path counts**, not tournament sizes, forbidden primes, or forbidden counts for arbitrary graph pairs. In fact `I(K3,2)=7` and `I(P4,2)=21`.
- The Frobenius/general-permutation defect transports exactly: for the lifted permutation, the number of undirected edge disagreements in either codec is `2Dσ(T)`. The inherited minimum `min_T Dσ(T)=# even vertex cycles of σ` therefore remains a constrained-image graph-edit theorem.
- Canonical transitive clone modules can be compressed exactly **only with the oriented quotient and module sizes retained**. Their undirected path graph alone is not sufficient.

## 1. Inheritance, correction lineage, and live concept board

The bounded inheritance pass used the navigation/guardrails, research protocol, relevant canon, Session62 clone theorem, and Session65 signed-pair result. Source bytes and Git blob IDs are in `SOURCE_PINS.json` and `inputs/`.

| Role | Exact source or object | Surviving mechanism / boundary |
|---|---|---|
| Closest proved scalar mechanism | THM-002, `THM-002-ocf.md` | Directed odd-cycle collections give `H=I(Ω,2)` |
| `{7,21}` theorem | THM-1370, `THM-1370-h-spectrum-omits-7-21-all-n.md`; repaired `rotation_owners_spectrum_20261020e/spectrum.md` | Strong minima plus SCC multiplication exclude both values; complete range remains open |
| Canonical hostile | `source ⊕ C3` versus its converse | Same scalar and natural undirected shadows, different isomorphism classes |
| Corrected near miss | Historical THM-209 calls cycles vertex sets and says proof needed | Current OCF is proved; literal directed cycles modulo cyclic rotation must be counted, with multiplicity on a support |
| Least-used sidecar | Shared vertex identification / port matching / orientation frame | Restores what independent graph quotients discard |
| Clone mechanism | THM-4463, `homogeneous-tournament-contacts-and-pinned-scale-floor.md` | Homogeneous-pair components are oriented paths, exactly maximal transitive modules |
| Signed action | Session65 `GENERAL_PERMUTATION.md` | Negative return on pair orbits, not permutation sign alone, determines minimum defect |

The source research protocol asks for a clean startup/freshness pass; this worker was instead explicitly assigned a **fixed read-only pin**. No startup script that could modify the repository, fetch, rebase, or publication was run. This preserves the assigned reproducibility boundary rather than claiming a fresh remote state.

### What the old `{7,21}` result actually says

Let `H(T)` count Hamiltonian paths as **vertex words**, not modulo automorphisms. Strong components form a dominance order, so paths concatenate and

`H(T)= product_C H(C)`.

The inherited Busch minimum for a strong tournament of order `n>=3` is

```
f(n)=3*5^((n-3)/3), n=0 mod3;
     5^((n-1)/3),   n=1 mod3;
     9*5^((n-5)/3), n=2 mod3.
```

Thus the first minima are `3,5,9,15,25,45,75,125`. The complete strong spectra through six (independently reproduced here) exclude seven and twenty-one; the floor excludes larger strong witnesses. A product equal to seven or twenty-one would require a strong factor seven or twenty-one. This proves the two global omissions. The minimum theorem is cited to [Busch 2006, Corollary 1](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf); the local finite computations are separate evidence.

Correction boundaries: the global range being *all* positive odds except these two is **OPEN**. `H=63` is attained, so a putative `7*3^k` forbidden tower is false. The repository retracted an earlier novelty claim: the 2016 discussion already had the omission argument. The false Leonardo floor and the later false replacement `f(n)=3f(n-2)` are not inherited; the correct recurrence is `f(n+3)=5f(n)`. Paley's seven-vertex tournament has ordinary `H=189` and automorphism order21: neither21 nor7 is forbidden as an automorphism datum.

## 2. Three meanings of “a pair of graphs”

For simple undirected graphs `G1,G2` on a common set `W`:

1. **Labelled pair:** the identities of the vertices are retained.
2. **Joint unlabelled pair:** `(G1,G2)~(G1',G2')` iff one bijection `f:W→W'` is simultaneously an isomorphism of both layers. Layer order is retained.
3. **Independent unlabelled pair:** retain only `([G1],[G2])`; the two isomorphisms may use unrelated bijections. This forgets the cross-layer matching. The graphs need not have equal orders, unless separately required.

Joint pairs are two binary symmetric relations on one carrier (equivalently, each pair has one of four edge-membership states). Independent pairs are strictly coarser. Merely asserting that a class “maps to” a pair asserts no injectivity. Here **faithful** always means an iff on isomorphism classes; where stronger automorphism preservation is available it is stated separately.

### Minimal order-loss hostile

Fix order `0<1<2`. Split each tournament arc into a forward graph `F` and backward graph `B`; these are complementary simple graphs. From the transitive tournament `0→1→2,0→2`:

- reverse adjacent arc `0→1`: the result is still transitive, with `H=1`;
- reverse long arc `0→2`: the result is the directed triangle, with `H=3`.

In both cases `B` is a single edge plus an isolate and `F` its complement. Therefore the pairs are even **jointly** isomorphic. The missing coordinate is the total order, not just the cross-layer identification. Retaining a labelled order makes this split invertible, but a freely chosen order is not an intrinsic class invariant. Canonicalizing the tournament first is a separate computational operation, not a free sidecar.

## 3. Faithful joint pair on `7n` vertices

### Construction

Use the asymmetric seven-vertex tree `R` with one degree-three center and three arms of lengths `1,2,3`. Distinguish its vertices intrinsically by their positions on these unequal arms. Let `t` be the endpoint of the length-one arm and `h` the endpoint of the length-two arm.

For a tournament `T` on vertex set `V`, form `W=V×V(R)` and define:

- `G1(T)` is the disjoint union of the copies `R_v`;
- `G2(T)` has the edge `{t_u,h_v}` exactly when `u→v` in `T`, and no other edges.

Both are simple undirected graphs. `G1` is usually disconnected; `G2` is bipartite and has `5n` unused isolates, with possibly further port isolates. There is no demand that either be connected.

### Theorem J (PROVED, elementary)

`T≅T'` iff the ordered pairs `(G1(T),G2(T))` and `(G1(T'),G2(T'))` are jointly isomorphic. Moreover every joint isomorphism is uniquely induced by a tournament isomorphism, and

`Aut(G1(T),G2(T)) ≅ Aut(T)`.

**Proof.** The center of `R` is its unique degree-three vertex. Removing it gives paths of distinct orders `1,2,3`; every arm is fixed, then every vertex is fixed. Thus `R` is asymmetric. Every `G1` isomorphism permutes its connected components, and has a unique map on each component. In particular it preserves the roles `t,h`. Its component permutation `f` preserves `G2` exactly when `u→v` iff `f(u)→f(v)`. Conversely every tournament isomorphism has this unique lift. QED.

There are `7n` vertices and `6n+C(n,2)` edges counted across layers. Construction and direct decoding take polynomial time (`O(n²)` in an adjacency-table representation), without solving tournament canonicalization. The theorem gives a polynomial reduction to **joint** pair isomorphism, not a new polynomial algorithm for that target.

### Why independent graph classes break this code

Reversing all tournament arcs transposes the incidence matrix. Swapping all tail and head ports is an ordinary graph isomorphism `G2(T)≅G2(T^op)`; all other vertices are isolated in that layer. The first layer is independent of T. Therefore

`([G1(T)],[G2(T)])=([G1(T^op)],[G2(T^op)])`.

At order four take `A=source⊕C3`, whose outdegree multiset is `{3,1,1,1}`, and `A^op=C3⊕sink`, whose multiset is `{0,2,2,2}`. They are nonisomorphic. This is minimal for a converse hostile because the tournament classes of orders at most three are all self-converse. The swap of tails and heads is **not** a frame-graph isomorphism; demanding one common bijection is precisely what blocks it.

## 4. Faithful independent-unlabelled pair via a `9n`-vertex graph

### Construction of `Γ(T)`

For each original vertex `v`, create three core vertices `c_v,t_v,h_v`. Join `c_v` to `t_v` and `h_v`. Attach:

- one new pendant leaf to `c_v`;
- two new pendant leaves to `t_v`;
- three new pendant leaves to `h_v`.

Finally add `{t_u,h_v}` exactly when `u→v`. The result is an ordinary uncoloured simple graph with `9n` vertices and `8n+C(n,2)` edges.

### Theorem I (PROVED, elementary)

`Γ(T)≅Γ(T')` iff `T≅T'`. Thus `T↦([Γ(T)],[empty_(9n)])` is a faithful independent-unlabelled pair, even if both graphs are required to have the same order.

**Proof.** In the construction every leaf has degree one and every core has degree at least three. Hence leaf/core status is intrinsic. Every core's number of leaf neighbours is respectively `1,2,3`, so the center/tail/head types are intrinsic. Each center has precisely one tail and one head core neighbour. Consequently the center-to-port matching is recovered inside this one graph, not via a second graph's labels. Between these matched ports, `{t_u,h_v}` means exactly `u→v`. Every graph isomorphism induces a tournament isomorphism, and every tournament isomorphism extends to a graph isomorphism. QED.

Construction and decoding again take `O(n²)` adjacency-table time. No canonical ordering is used. Given a graph canonical-form oracle, one can canonicalize `Γ(T)`, decode its intrinsically recognized centers in the oracle's order, and obtain a canonical representative of T. This does not assert that such an oracle is free or prove any new complexity separation/equivalence.

**Automorphism boundary.** This second codec is faithful on classes but not automorphism-identical:

`Aut(Γ(T)) ≅ (S2×S3)^n ⋊ Aut(T)`.

The kernel independently permutes the identical tail leaves and head leaves; centers, port types, and their matching are already intrinsic. The quotient action on centers is exactly `Aut(T)`. Hence Γ usually has even-order automorphisms despite T having only odd-order automorphisms. This is not a contradiction: fidelity on isomorphism classes is weaker than equality of automorphism groups.

## 5. The genuine `{7,21}` graph-pair bridge

Let `Ω(T)` have one vertex for each literal directed odd cycle of length at least three, **modulo cyclic rotation only**; two vertices are joined when their cycles share an original vertex. Two distinct directed cycles on the same support remain distinct adjacent vertices.

The inherited OCF says `H(T)=I(Ω(T),2)`. The primary source is [Irving–Omar, Corollary20](https://arxiv.org/html/2412.10572v3), which attributes the result to Grinberg–Stanley. Reversing every cycle gives `Ω(T)≅Ω(T^op)` and preserves cycle intersections.

The pair

`([Γ(T)],[Ω(T)])`

therefore has two different exact jobs:

- Γ reconstructs the tournament isomorphism class;
- Ω exposes Hamiltonian-path count as an independence-polynomial evaluation.

This is a faithful pair of independent simple-graph classes. **Size boundary:** Γ is linear in vertex count; Ω may contain exponentially many cycles (indeed the maximum cycle count can be much larger than exponential in n), so this informative two-coordinate construction is not claimed polynomial-size or polynomial-time. The empty-second-layer version in §4 supplies the polynomial-size existence theorem.

### Image restrictions are the mathematical content

For arbitrary simple graphs, `I(K3,2)=1+3*2=7` and `I(P4,2)=1+4*2+3*4=21`. Consequently neither graph can equal an entire `Ω(T)`. The tournament result is a **realizability restriction on the image of Ω**, not a universal graph theorem. Nor does this forbid K3 or P4 as induced subgraphs or arbitrary components without a separate argument.

Likewise an arbitrary graph pair need not lie in the image of either codec. A faithful recoding transports an already proved exclusion onto a precisely defined image; it does not by itself explain that image or solve the open spectrum-completeness question.

### A tempting natural pair fails minimally

Let `J(T)` be the homogeneous-pair graph from §7. Both Ω and J are unchanged up to isomorphism by converse, so `(Ω,J)` is not a faithful replacement for `(Γ,Ω)`. On `source⊕C3` and its converse, both Ω are `K1`, both J are edgeless four-vertex graphs, and both H equal3. Their different score multisets prove nonisomorphism.

## 6. Symmetry obstruction and exact defect transport

### Same-vertex naturality no-go (PROVED)

There is no relabelling-equivariant encoding of every labelled tournament into **any fixed number of named simple-graph layers on the same original vertex set** with a relabelling-equivariant exact decoder.

**Proof.** On three vertices, the directed cycle is preserved by the cyclic group C3. Equivariance forces every output graph layer to be C3-invariant. C3 is transitive on the three unordered pairs, so each layer is empty or complete. The whole tuple is consequently preserved by S3, including a transposition τ. Equivariance of a decoder would imply `τT=D(τE(T))=D(E(T))=T`, contradiction because τ reverses the cycle. Equivalently `E(T)=E(T^op)` as labelled outputs. QED.

The statement is deliberately scoped: it does **not** forbid an abstract injection of tournament isomorphism classes into same-size graph-class pairs, and does not forbid symmetry-breaking/non-natural choices or enlarged vertex carriers. The codecs above enlarge the carrier and intrinsically distinguish source and target roles. A bare cardinality argument proves none of these stronger properties.

### Exact graph form of the Session65 defect (PROVED from construction)

For any permutation σ of V, lift σ to all gadget copies, retaining every local gadget coordinate. Define

`Dσ(T)=#{ {u,v}: A(u,v) != A(σu,σv) }`.

For the joint code, every frame edge is preserved by this lift. For Γ, every matching and pendant edge is preserved. For each original unordered pair `{u,v}`, exactly one of the two possible port edges `{t_u,h_v}`, `{t_v,h_u}` is present. An orientation disagreement exchanges these two edges, causing two edge-membership disagreements. Different unordered pairs give disjoint port-edge pairs. Therefore

`|E(G2(T)) Δ E(σG2(T))| = |E(Γ(T)) Δ E(σΓ(T))| = 2Dσ(T)`.

The same identity holds using `σ^-1`; the convention for image versus pullback changes neither count. It is an exact labelled action statement, not an optimization over all graph permutations.

Session65 proves `min_T Dσ(T)=e(σ)`, the number of even cycles of σ on **vertices**, by signed return along unordered-pair orbits. Hence the minimum graph symmetric difference **among these encoded tournament graphs and under this prescribed lifted action** is `2e(σ)`. An involution with b swapped vertex pairs has minimum2b in this graph metric. A positive-sign permutation with two transpositions still has minimum4, so sign alone is insufficient. In an actual arithmetic Frobenius application, the original arithmetic set and its actual Frobenius permutation must still be justified and retained.

Finite tournament automorphism groups have odd order, since a nontrivial involution exchanges some vertex pair and reverses its arc. This is reflected exactly in the joint pair's automorphism group. Γ's extra even leaf kernel is harmless to class fidelity, but must not be confused with a new tournament symmetry.

## 7. What clone merger can and cannot simplify

Define a homogeneous pair `{u,v}` by agreement of their arcs to **every exterior vertex**, not equality of their complete adjacency rows. The inherited THM-4463 proves that J(T)'s components are consistently oriented paths, and their vertex sets are precisely maximal transitive modules.

Thus there is a canonical compression:

1. find these modules;
2. retain the tournament quotient Q between them;
3. retain each positive size `m_v`;
4. reconstruct `T=Q[TT_(m_v)]` by the unique transitive order inside each block.

This is lossless up to isomorphism. Every automorphism of T induces a size-preserving automorphism of Q, and each such automorphism lifts uniquely, because a transitive tournament has no nonidentity automorphism. An encoded weighted quotient can therefore replace a full codec if the sizes are encoded as genuine retained data; the undirected path components do not supply the directions between modules.

**Minimal equivalence-relation hostile:** in `TT3`, pairs `{0,1}` and `{1,2}` are homogeneous while `{0,2}` is not. Endpoint-deleted row agreement is not transitive. The lawful unit is the complete canonical path component, not an arbitrary relation obtained by declaring one coordinate equal.

**Higher clone warning:** Session62's exact cut-subset orbit polynomial is controlled by the full odd-order automorphism action on module gaps. Its root multiplicity at−1 is `min_g #cycles(g)`, while the rank of legal fixed symmetric-difference toggles is `# group orbits`. These can differ. Passing to independent graph classes, or merely retaining module sizes, can lose the action and cannot license a labelled gap toggle. Our joint codec preserves the full automorphism action; Γ preserves it as the quotient action on recovered centers. A particular cut/gap labelling is still a separate marked input.

## 8. Verification and boundary of completion

Run from this directory with Python3 (standard library only):

```
python3 check.py > normal.json
python3 -O check.py > optimized.json
cmp normal.json optimized.json
```

All conditions are explicit exceptions, not removable assertions. The checker performs:

- every labelled tournament on orders0–5 (1,100 total): both codec round trips, weighted clone reconstruction, independent comparison of Hamiltonian-path permutation count, subset DP, and cycle-conflict independence recursion;
- every isomorphism class through five, with known counts `1,1,1,2,4,12`: all same-order class pairs for both codecs checked by a generic exact edge-coloured graph-isomorphism backtracker that does **not** call the decoders, plus independently scrambled positive controls;
- every labelled tournament through order6 for the two forbidden-value controls and complete strong spectra at3–6; this is finite evidence, not the proof of the all-order theorem;
- every tournament/permutation through order4: exact two-to-one graph-disagreement identity for both codecs;
- all8 graphs on three vertices for the C3-invariance obstruction;
- minimal order-loss, converse/alignment-loss, false natural-pair, arbitrary-graph7/21, clone-equivalence, and even-permutation defect hostiles.

No all-order tournament enumeration, new graph-isomorphism complexity bound, full spectrum classification, arithmetic canonical tournament, or new proof of a major open conjecture is claimed. The mathematical continuation worth pursuing is **characterizing the image of Ω or of a smaller statistic-bearing codec**, not proving the mere existence of a graph encoding again.
