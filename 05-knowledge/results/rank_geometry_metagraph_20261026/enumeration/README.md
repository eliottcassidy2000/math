# Exact class generation from Hamiltonian charts and canonical deletion

2026-10-26. **PROVED (elementary algorithm correctness), IMPLEMENTED, FINITE-EXACT through order eight.** The primary algorithm constructs every tournament isomorphism class exactly once by adjoining a vertex. It does not use the Burnside sum to generate or count its output. Its finite output agrees with the existing independent reverse counting recurrence. This is an explicit construction and implementation of the canonical construction-path method, not a novelty claim or an asymptotically smaller scalar recurrence.

The counts at orders zero through eight are

```
1, 1, 1, 2, 4, 12, 56, 456, 6880.
```

All-order correctness means that the terminating algorithm has the stated output for every finite input order. The run here stops at eight; it is not an enumeration of infinitely many orders or a practical claim at arbitrarily large orders.

## 1. Inheritance and the information that must survive

- Closest finite class census: [Session66 graph-pairs, §§7–8](../../session66_root27_aliquot_graph_interfaces_20261026/session66-graph-pairs/REPORT.md), with all classes through five, faithful graph codecs, and weighted clone reconstruction.
- Closest aggregate count mechanism: [the production reverse potential](../../lca_ports_reverse_20261025/counting/README.md), then the [incoming forced-skip contraction](../../lca_ports_reverse_20261025/counting/INCOMING_COUNTING.md). Those count without constructing every class; the present catalog is complementary and is not a replacement performance claim.
- Canonical hostile: adjoining a vertex produces **rooted** classes. Their multiplicity over an unrooted tournament is its number of vertex orbits, not its order. The old `04-computation/a000568_fiber_recursion.py` already warns that dividing by n is not exact.
- Corrected near miss: homogeneous pairs do not themselves form an equivalence relation. In TT3 its consecutive pairs are homogeneous, but its endpoints are not. Their connected components, with their transitive order, are the canonical modules.
- Least-used sidecar: a canonical labeling **and all minimizing transports**, so the exact automorphism action and the deletion orbit survive. A group order alone cannot act on extension masks or on weights.

The live board is Hamiltonian chart / rooted deletion / exact automorphism action / weighted transitive core / literal port rank / aggregate count. The Anchor is constructive enumeration; the Niche is exact weighted-core allocation; the Wildcard is whether forgetting the owner action already fails in the smallest symmetric examples. It does.

The general canonical construction-path technique is classical: Brendan D. McKay, [*Isomorph-Free Exhaustive Generation*, J. Algorithms 26 (1998), 306–324](https://doi.org/10.1006/jagm.1997.0898), with the [author's corrected copy](https://users.cecs.anu.edu.au/~bdm/papers/orderly.pdf). Primary publication metadata/abstract and author-copy availability were checked on 2026-10-26. The elementary specialization below is proved here; no external algorithmic complexity result is needed.

## 2. The Hamiltonian chart is an exact finite family

Put a fixed directed spine

```
0 -> 1 -> ... -> n-1.
```

Each of its nonconsecutive pairs receives one free direction. For n>=1 the family therefore has exactly

```
2^((n-1)(n-2)/2)
```

members. Use one empty chart at n=0. A word stores the nonconsecutive pair bits in lexicographic pair order; bit one means i->j for i<j. `unchart(n,word)` decodes the whole tournament, including the spine.

Every tournament has such a presentation. To prove this, insert vertices one at a time into an existing directed path: put v immediately before the first path vertex that v beats, or append v if there is none. Every preceding vertex beats v and v beats its successor, so a Hamiltonian path results. `hamilton_path` implements that argument; `chart` returns the word and the actual transport from chart positions to original vertices.

Two words represent one unmarked tournament class exactly when their decoded tournaments are isomorphic. Thus an all-order, non-Burnside enumeration is already

```
R_n = { canonical(unchart(n,w)) : all chart words w },
A000568(n) = |R_n|.
```

This is a faithful generating family with an explicit duplicate criterion. It is not yet an efficient implementation: order eight has 2,097,152 words. The augmentation below reduces the tested candidates to 54,256 rooted extensions at that order.

The fiber over T has exactly `tc(T)=H(T)/|Aut(T)|` words. Isomorphisms from chart presentations to T amount to Hamiltonian vertex words, and the free Aut(T) action on complete words identifies precisely the same presentations. This is the inherited free-action lemma **LEM-003**, [aut acts freely on ham paths](../../../../01-canon/theorems/LEM-003-aut-acts-freely-on-ham-paths.md), with its [Session66 tc clarification](../../session66_root27_aliquot_graph_interfaces_20261026/session66-tc-scope/REPORT.md). It is not H(T) itself when T has symmetry.

## 3. Exact canonicalization and the transported group

The implementation stores one direction for each unordered pair; it does not forget orientation or identify a tournament with its converse. `canonical(n,code,colors=None)` uses the following finite search.

1. Start with the ordered color classes, or one class when uncolored.
2. Split each cell by the vector of outgoing-neighbor counts into all current cells. Order resulting pieces lexicographically by these vectors and iterate until no split occurs.
3. If a nonsingleton cell remains, take the first such cell. Branch on each vertex v, making `{v}` a first singleton followed by the rest of that cell, and refine again.
4. A discrete partition gives an ordering p, with `p[new]=old`. Retain the least integer adjacency code over all leaves and **every** ordering attaining it.

Every rule is preserved by isomorphisms, including branch sets and their chosen cell. Thus an isomorphism maps the entire search tree bijectively, preserves all leaf codes, and sends the minimum leaves to the minimum leaves. Conversely equal output codes give an explicit isomorphism by composing their leaf orderings. This proves the iff on isomorphism classes. The rule need not produce the lexicographic minimum over all n! vertex orders; it only needs this equivariance and exact decoding.

Fix one minimum ordering p0. For each other minimum p, the map `p0[i] -> p[i]` is an automorphism, because both induce the same labeled tournament. Conversely each automorphism sends p0 to a distinct minimum leaf. These are therefore **all** automorphisms, without multiplication closure or a missing-generator assumption. The code can express them on original labels or canonical labels. Equal initial colors remain separate throughout; the analogous colored statement follows verbatim.

The selected deletion orbit is

```
d(T) = { p[n-1] : p ranges over all minimum orderings }.
```

It is one Aut(T)-orbit, and it transports under every tournament isomorphism. Choosing a labeled last vertex without this orbit would not be invariant.

## 4. Canonical augmentation gives an exact recurrence with no accepted duplicates

Assume `R_(n-1)` contains exactly one labeled representative of every class. For each parent P, let Aut(P) act on subsets S of V(P). For one representative of every subset orbit, form `P+S` by adjoining v with `u->v` iff `u in S`. Accept this extension exactly when `v in d(P+S)`.

**Theorem A (PROVED).** The accepted children contain exactly one representative of every n-vertex tournament class. In particular

```
A_n = sum_[P in R_(n-1)] sum_[S in subsets(V(P))/Aut(P)]
          1_[new vertex belongs to d(P+S)].
```

Here every sum is a finite explicit computation, with `R_0={empty}`. It is an unrestricted structural recurrence. Its state is the actual canonical parent and its group action; it is not a closed recurrence on previous scalar counts.

**Proof.** Root a tournament T at any vertex v in d(T). Its deletion T-v is isomorphic to exactly one parent representative P. Transporting the arcs to v gives a subset S. Different choices of the parent isomorphism differ by an automorphism of P, so they give one subset orbit. Conversely two subset choices over P produce isomorphic rooted children iff the isomorphism restricts to an automorphism of P taking one subset to the other. Thus parent-plus-subset-orbit pairs are exactly rooted tournament classes.

The rooted classes above an unrooted T are exactly its vertex orbits. The acceptance rule selects the single orbit d(T), so it selects exactly one rooted class and hence one extension pair. This proves existence and uniqueness. QED.

`enumerate_classes` treats repeated acceptance as an error. A set insertion is not allowed to conceal an incorrect uniqueness proof. Each output record retains the canonical code, the complete automorphism action, the parent code, and the extension mask; the parent/mask refer to the pre-canonical adjoining operation. Canonical label transports are recomputable by `canonical` when a literal operation needs them.

**Exact multiplicity boundary.** Let r(T) be the number of Aut(T)-orbits on vertices. The unfiltered candidate number is

```
sum_[P] # (subsets(V(P))/Aut(P)) = sum_[T of order n] r(T).
```

Replacing r(T) by n fails already at n=3: the transitive class contributes three vertex orbits and C3 contributes one. Four rooted classes divided by three do not equal the two unrooted classes. Equivalently one may weight each rooted extension by `1/r(P+S)`, but the canonical acceptance recurrence uses only zero/one weights.

## 5. Weighted primitive cores: exact, with the action and every merger retained

Here “primitive core” means **pair-free**, not an indecomposable/prime tournament in every modular-decomposition convention.

The inherited **THM-4463**, [homogeneous tournament contacts and pinned scale floor, §2](../../../../01-canon/theorems/THM-4463-homogeneous-tournament-contacts-and-pinned-scale-floor.md), identifies homogeneous-pair graph components with maximal transitive modules. `clone_core` computes these literal blocks, their unique internal orders, the quotient Q, and the positive size vector w. It explicitly checks homogeneous external arcs and transitivity. Flattening its ordered blocks gives the literal reconstruction

```
T = Q[TT_(w_v)]
```

under the returned label transport. Its quotient is pair-free: a homogeneous pair of quotient vertices would join the corresponding transitive modules into a larger transitive module. Conversely, expanding a pair-free Q creates no homogeneous pair between distinct blocks: such a pair would force agreement against every other quotient block, contradicting pair-freeness. Thus precisely the supplied blocks are the canonical modules.

Isomorphisms of inflated tournaments are exactly the weight-preserving isomorphisms of their quotient cores, with unique internal lifts because a transitive tournament is rigid. Consequently

```
A_N = sum_(q=1..N) sum_[Q pair-free, |Q|=q]
      # { w in positive integers^q : sum(w)=N } / Aut(Q),   N>=1.
```

The slash means **orbits**, not division by the group order. `weighted_core_allocations` enumerates compositions and takes one least tuple per actual group orbit; it does not call Burnside. Inflation yields each class once. At order eight this splits the 6,880 classes by core size as

| Core order | 1 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---:|---:|---:|---:|---:|---:|---:|
| Classes | 1 | 7 | 26 | 107 | 444 | 1693 | 4602 |

**The remaining cost is explicit.** The all-one allocation at q=N includes every N-vertex pair-free core. This formula does not avoid enumerating those cores, and the present implementation obtains them from the complete class catalog. No recurrence depending only on core counts is asserted.

Two smallest useful hostile controls are retained. First, `TT2[TT1,TT2]` and `TT1[TT3]` both inflate to TT3: arbitrary block partitions are not unique cores. Second, the cyclic triangle with weights `(1,2,3)` and `(1,3,2)` gives distinct six-vertex classes. Their size multisets agree, but C3 has rotations, not all six permutations. Discarding the orientation or the weight-to-vertex incidence merges them incorrectly.

## 6. Ports, graph codecs, and metagraph scope

An ordered port list may repeat vertices. `marked_key` uses the full incidence vector of port occurrences as colors; it therefore preserves both distinct ports and shared endpoints. `port_coordinates(blocks,ports)` retains each port's **block and internal rank**. Retaining only the block is insufficient: in TT3, ports at ranks `(0,2,0)` and `(0,1,0)` give different marked objects despite identical block labels. All these transports are literal and testable. The enumerated A000568 classes are unmarked; no assertion identifies different marked port configurations.

The source/target map is explicit: a marked Hamiltonian chart decodes the entire tournament, then canonicalization forgets its chosen spine and labels while retaining its unmarked class. Reconstruction of the marked object needs the chart word and spine transport; reconstruction after clone compression needs Q, vertex-attached weights, and any port ranks. Merely knowing H, tc, a degree multiset, or a clone-size multiset is not the implemented equivalence test. At order four, a source joined to C3 and its converse already have H=3 and different classes.

The Session66 joint `7n` graph-pair codec can substitute for the literal tournament because it preserves the exact automorphism action. Its independent `9n` graph codec Gamma preserves the class but has the leaf-permutation kernel `(S2 x S3)^n`; one must use the action on recovered centers. Using the full graph automorphism order as if it were Aut(T) changes the chart fiber by `12^n` and is wrong. This is precisely the [Session66 leaf-kernel/tc scope](../../session66_root27_aliquot_graph_interfaces_20261026/session66-tc-scope/REPORT.md), not a new codec claim.

The full family of chart words **does suffice**, with the exact equivalence test above. This proof uses neither black-core connectivity nor any asserted classification by an abstract merged-metagraph shape. It therefore cannot promote such a shape-only classifier: that would require a separate injectivity/reachability theorem with the chosen blue/black operations and their retained owners stated. The code retains the class owner and actual group action rather than assuming that missing theorem.

### Full path inventory and exact operations

`hamilton_paths` lists every literal Hamiltonian vertex word. For any retained ordered port list L, it consequently gives the complete inventory

```
{ (chart(T,P), positions_of_L_in_P) : P a Hamiltonian path of T }.
```

The multiplicity of each inventory entry is exactly `|Aut(T,L)|`, where every port occurrence is fixed. Indeed equal entries induce an automorphism between the two vertex orders, and equality of port-position lists says that this automorphism fixes L. Conversely every such automorphism preserves the entry and acts freely on complete vertex words. Thus the distinct inventory has `H(T)/|Aut(T,L)|` entries. This is the exact marked version of the chart fiber, and retains coincident endpoints. Reframing along any listed path transports every port by the inverse path permutation, with no lost vertex identification.

Two explicitly implemented operations are exact on this framed space:

1. K flips every nonconsecutive arc and fixes all spine arcs: its chart word is the bitwise complement. This is an involution at the **retained** spine P; ports keep their original vertices and their frame positions.
2. R takes the converse tournament and reverses the spine. Its free-pair word is reflected by `(i,j)->(n-1-j,n-1-i)` without bit complementation; a port position i becomes `n-1-i`.

They commute as chart operations. The checker compares both to literal original-label arc edits, rather than checking only a word-level tautology, and compares the corresponding marked canonical objects. It also finds an explicit source class whose different spines give different unmarked K-target classes. Therefore K does **not** descend to a function of the unmarked class alone. A class transition relation can retain all these choices, but choosing one without its spine loses information.

The first retained K hostile is the order-four source code34, with out-neighbor bitmasks `(4,1,10,3)`. Paths `(0,2,3,1)`, `(1,0,2,3)`, `(2,3,1,0)`, and `(3,1,0,2)` give canonical K-target codes9,0,34,40 respectively: all four order-four classes occur. These codes use this producer's canonical convention, not necessarily the globally least upper-pair integer.

The complete re-presentation relation is reconstructible from this inventory: two chart words lie over the same class iff a transported Hamiltonian path relates their full adjacency tables. This supplies all presentations and their literal port maps, not a proof that a smaller selected set of blue/black moves connects them. [The information lane](../information/RESULT.md) gives the corresponding class-mass accounting and distinct scalar-information boundaries.

## 7. Executable scope and receipts

Standard-library Python, from this folder:

```powershell
python enumerate.py --max-n 8 --catalog catalog.json
python check.py > normal.json
python -O check.py > optimized.json
```

The catalog option writes actual canonical codes, actions, parent records and chart words. The committed producer and receipts do not require a stored catalog; it is reproducibly generated on demand.

`check.py` performs the following exact universes, with non-removable exception gates:

- All 1,100 labeled tournaments of orders 0–5: compare against an independent full-permutation matrix canonicalization, including the **entire** automorphism group, false-merge/false-split detection, and label transport.
- All 1,552 relabelings of the generated representatives through five: canonical invariance.
- All 33,869 Hamiltonian chart words of orders 0–7: literal round trips, onto all classes, and every fiber exactly H/Aut.
- All 7,413 classes of orders 0–8: literal automorphism preservation; Hamiltonian subset-DP counts; chart mass; labeled orbit-stabilizer mass; the rooted-extension identity; exact clone reconstruction and symmetry; and one-to-one weighted-core allocation back to the complete class set.
- All 1,335 literal Hamiltonian paths of every generated class through six: full word multiplicities, exact K/R arc transport, commuting involutions, and repeated-port transport. All 388 ordered port inventories `(u,v,u)` through five (including the empty control) check the exact pointwise-stabilizer multiplicity.
- Orders 0–8: independent agreement with the already published integer reverse recurrence. That oracle is a verification route, not part of the non-Burnside class generator.
- Minimal hostiles for division by n, an H collision, a nonmaximal clone partition, an orientation-forgetting weight multiset, and repeated/internal-rank ports.

At n=8 the generator considers 54,256 rooted candidates, visits 58,520 discrete canonical-search leaves, and accepts 6,880 classes. A direct chart scan would have 2,097,152 inputs; an all-labeled scan would have 268,435,456 inputs. These are measured work counts for different candidate spaces, **not** a fair timing or asymptotic performance comparison with established tournament generators or the aggregate reverse count engine.

Worst-case canonicalization is finite factorial search; a coarse bound for the next layer is `O(A_(n-1) 2^(n-1) n! poly(n))`, plus exact group/subset processing, and storing all output classes is inherently proportional to that output. The refinement helps greatly in this bounded run but receives no polynomial-time claim. The structural recurrence is complete; a smaller aggregate recurrence derived only from metagraph family types remains a separate target.

The normal and optimized receipts are byte-identical; SHA256 of each is `ec87ed60f39a3146127a1cf2454ae7e5efb53aeaaa985ef05303e8ecd95cfffd`.
