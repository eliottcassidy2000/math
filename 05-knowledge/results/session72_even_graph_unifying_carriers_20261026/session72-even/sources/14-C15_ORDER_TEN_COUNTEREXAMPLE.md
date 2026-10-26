# A ten-vertex counterexample to the mixed-prime symmetry refinement

Status: **PROVED below; independently check before promotion. FINITE-EXACT verification included.**

Let `H=< (0 1 2)(3 4 5 6 7) > <= S_10`, fixing vertices 8 and 9. Thus `H` is cyclic of order 15, with orbit sizes `3,5,1,1`. Write `tau_T(H)` for the number of tournament isomorphism types whose automorphism group contains a conjugate of `H`, and `tau_E(H)` for the analogous count of automorphism/sign-even graph types. Then

**Theorem.**

`tau_T(H)=32`, while `tau_E(H)=31`.                              (1)

Consequently HYP-9340(a), as stated for every odd-order permutation subgroup, is false. This leaves intact the proved Sylow statement in THM-4721, the new p-group refinement, and the relative-even normalizer identity. The discrepancy is precisely a globally odd graph, `K_{5,5}`, which remains even relative to the normalizer of this chosen `H`.

## 1. Inheritance, source freshness, and exact universe

The inherited source is

`01-canon/theorems/THM-4721-twisted-burnside-on-odd-fixed-point-sets-sylow-refinement-of-even-graphs-and-tournaments-components-walsh-basis.md`,

with hypothesis

`05-knowledge/hypotheses/HYP-9340-refined-equinumerosity-for-every-odd-order-symmetry-and-what-a-tournaments-even-graphs-bijection-can-be.md`.

Before computing, we inspected incoming commit `d5c80ed23f810a0f78387a2a5ee5a3ba81291328` and its inherited script/output

- `04-computation/experiments/even_graphs_20261019/subgroups_fix.py`;
- `04-computation/experiments/even_graphs_20261019/subgroups_fix_n9.out`.

Those verify the disjoint `5+3` action at order eight (`2=2`) and `5+3+1` at order nine (`8=8`). They do not check the present ten-vertex action. No earlier computation was rerun unnecessarily.

Set `A={0,1,2}`, `B={3,4,5,6,7}`, and `R={8,9}`. Since 3 and 5 are coprime, `H` contains the cycles on `A` and `B` separately. There are nine orbits of `H` on unordered pairs:

- one within `A`;
- two within `B`;
- one between `A` and `B`;
- four between a nontrivial orbit and a vertex of `R`;
- one within `R`.

Hence there are exactly `2^9=512` invariant labelled graphs and exactly 512 invariant labelled tournaments. The final counts are full `S_10`-isomorphism-type counts, not a silent substitution of normalizer orbits for types.

## 2. Tournament count: 32 full isomorphism types

On `A`, invariance forces a directed triangle; its two orientations are isomorphic. On `B`, invariance gives a circulant tournament on five vertices. Its connection set chooses one element from each of `{1,4}` and `{2,3}`. The four choices are all multiplier-equivalent: multiplication of `{1,2}` by `1,2,3,4` gives all four choices. Let `Q_5` be the fixed representative with arcs `i -> i+1,i+2 mod 5`.

All arcs between distinct `H`-orbits are uniform. Thus every invariant tournament is isomorphic to a substitution into a four-vertex quotient tournament, with modules `C_3,Q_5,{8},{9}`. There are 64 quotient presentations. The only possible extra identifications are settled by the following elementary lemma.

**Prime-module lemma.** A strongly connected prime tournament module of order at least three cannot overlap another module without one containing the other. Consequently two distinct such prime modules are disjoint.

*Proof.* Suppose prime strong module `U` overlaps a module `W`, with all three sets `U∩W`, `U\W`, `W\U` nonempty. The intersection is a proper module in `T[U]`: each vertex of `U\W` lies outside `W` and sees its intersection uniformly. Primeness forces the intersection to be a singleton `{u}`. Fix `w in W\U`. For each `x in U\W`, the module properties give

`x -> u iff x -> w iff u -> w`.

The last statement does not depend on `x`, so `u` is a source or a sink in `T[U]`, contradicting strong connectivity. Thus modules cannot overlap. If two prime modules were properly nested, the smaller would be a nontrivial proper module of the larger, also impossible. QED.

Both `C_3` and `Q_5` are prime and strong. For `Q_5`, a module of order three or four would force an outside vertex to dominate or be dominated by at least three vertices, impossible because every vertex has indegree and outdegree two. No pair is a module: rotation reduces to the pairs `{0,1}` and `{0,2}`, distinguished respectively by vertices 3 and 1. Strong connectivity follows directly from the directed five-cycle.

The existing `C_3` and `Q_5` modules occupy eight vertices. By the lemma, any other prime strong module must be disjoint from both, unless identical to one of them. Only two vertices remain. Therefore the `C_3` module and the `Q_5` module are each intrinsically unique. Every isomorphism preserves them and can only interchange the two remaining singleton vertices. Conversely such a quotient isomorphism lifts.

The interchange of the two singleton vertices fixes no quotient tournament, because it reverses the arc joining them. Therefore the 64 quotient presentations occur in pairs, proving

`tau_T(H)=64/2=32`.

This proves the complete isomorphism classification; no assumption about normalizer orbits is needed.

## 3. Even-graph reduction: both nontrivial orbit blocks must be independent

Any automorphism of an induced `H`-orbit graph extends by the identity outside the orbit, since all cross-orbit adjacency is uniform. An odd internal automorphism is therefore an odd automorphism of the full graph.

Within `A`, the options are the empty graph and `K_3`; a transposition is odd for `K_3`.

Within `B`, the four options are the empty graph, a five-cycle, its complement (another five-cycle), and `K_5`. Each nonempty choice has an odd automorphism:

- in either five-cycle, reflection `i -> -i mod 5` reverses an odd number of edges (one of its two transposed vertex pairs is an edge);
- in `K_5`, any transposition reverses an odd number of edges.

Thus every full-even invariant graph has `A` and `B` independent. This eliminates exactly 448 of the 512 labelled invariant graphs. The remaining 64 graphs are the blowups of an arbitrary four-vertex graph with independent vertex weights `(3,5,1,1)`.

## 4. The relative quotient count is 32

Temporarily preserve the roles of the 3-block and 5-block, while allowing the two singleton vertices to be interchanged. Permutations internal to either independent block are even, so the only possible odd role-preserving automorphism is the interchange of the singleton vertices.

It is an automorphism precisely when the singletons have equal adjacency to each of the two large blocks. Among the 64 quotient masks:

- 16 are fixed by this interchange: one free bit for the large-block edge, two free common-adjacency bits, and the singleton edge bit;
- the remaining 48 form 24 two-element orbits;
- among the 16 fixed masks, precisely eight have the singleton edge present and are odd relative to this interchange;
- the other eight fixed masks are relatively even.

Thus the relatively-even quotient types number

`24+8=32`.                                                    (2)

The ordinary weighted quotient has `(64+16)/2=40` types.

## 5. Full isomorphism and sign: exactly one additional obstruction

Two vertices of a graph are twins if they have the same neighbors outside their pair. Twin classes are intrinsic under isomorphism. Each class induces either a clique or an independent set.

In a relatively-even candidate above, there can be no nontrivial clique twin class. A root cannot be a true twin of a vertex in either large independent block, since a second vertex in that block distinguishes them. The only remaining possibility would be the two singleton roots being adjacent true twins, and that was excluded in §4.

Therefore all nontrivial twin classes are independent. Permutations within such classes are even. A full automorphism can acquire a new sign only by permuting twin classes of the same cardinality.

The original 3- and 5-blocks each lie entirely within one twin class. There are two cases.

1. **They merge into one class.** That class has size at least eight and is uniquely large; all remaining vertices are roots.
2. **They lie in distinct classes.** Their class sizes are
   `3+r_A` and `5+r_B`, where `r_A,r_B>=0` and `r_A+r_B<=2`.
   All other classes contain only roots, so have size at most two. The two large class sizes can coincide only when `r_A=2,r_B=0`, giving sizes five and five.

Except in that final case, class sizes distinguish the original large-block roles. Additional permutations inside classes contribute no sign, and any interchange of singleton classes is already the interchange tested in §4. Also, relocating a root within its independent twin class leaves its quotient adjacency unchanged. Hence two relatively-even quotient types cannot become newly isomorphic after forgetting roles.

In the exceptional case, both roots join the 3-block to make an independent class of size five. If the two size-five classes were nonadjacent, they would merge into a single twin class. Therefore, as distinct twin classes, they are joined completely: the graph is exactly `K_{5,5}`. It has only one weighted quotient presentation type. Swapping its two parts reverses all 25 cross edges in a block-compatible reference orientation, so its sign is `(-1)^25=-1`.

Thus exactly one of the 32 relatively-even types is globally odd. The remaining 31 are full-even and pairwise nonisomorphic. This proves

`tau_E(H)=32-1=31`,

and completes the proof of (1).

## 6. Explicit missed graph and why the old theorem survives

For the fixed labeling, the exceptional graph has bipartition

`{0,1,2,8,9} | {3,4,5,6,7}`.

It is fixed by `H`. A concrete odd automorphism interchanging its parts is

`[3,4,5,0,1,2,8,9,6,7]`,

meaning that entry `i` is the image of vertex `i`. It is a product of five transpositions and reverses an odd number of graph edges.

The normalizer of the chosen `H` preserves its 3-orbit and 5-orbit individually, so it cannot contain this part swap. Relative to that normalizer, the graph is even. This is the missing Frattini step outside the ambient Sylow setting: normalizer-evenness need not imply global evenness. THM-4721 itself already identifies that logical danger; at this mixed-prime ten-vertex action, its compensation no longer balances in the full type count.

This is the first failure found beyond the inherited checked range; the inherited data covers all listed subgroup actions through nine vertices. We do not independently claim that its subgroup list is a formally proved exhaustive classification of every odd-order subgroup of `S_9`; “minimal” should be used only with that provenance qualified.

## 7. Reproduction and exact checks

Files:

- `c15_n10_counterexample.py`;
- `c15_n10_counterexample.json`;
- `c15_n10_counterexample_optimized.json`;
- `inherited_subgroups_fix.py`, recovered for provenance from the incoming commit, not rerun.

Run:

`python3 c15_n10_counterexample.py > c15_n10_counterexample.json`

`python3 -O c15_n10_counterexample.py > c15_n10_counterexample_optimized.json`

`cmp c15_n10_counterexample.json c15_n10_counterexample_optimized.json`

The verifier requires only Python's standard library. It:

1. constructs `H`, its nine pair orbits, and all 512 invariant labelled graphs;
2. checks an explicit odd automorphism for each of the 448 excluded graphs;
3. uses the intrinsic full twin quotient for the 64 remaining graphs, enumerates every quotient automorphism, lifts it to ten vertices, and compares direct edge-inversion sign with weighted quotient sign;
4. independently enumerates all 512 invariant labelled tournaments from oriented `H`-orbits;
5. classifies them using full directed-isomorphism backtracking with equitable color refinement, rather than normalizer equivalence;
6. checks every pair of the 64 normalized quotient presentations, verifying that the full isomorphism relation is exactly identity or singleton interchange;
7. prints the complete offending `K_{5,5}` witness.

Result: **6,326 explicit checks, PASS.** Normal and optimized outputs agree byte for byte. Graph count 31; tournament count 32; exactly one relatively-even type is globally odd. The full theorem proof is independent of the computation; the computation tests both its positive and hostile branches.

## 8. Research mechanism

The earlier sign-aware augmentation packet exposed the missing coordinate: parity on root-moving automorphisms is not carried by a root stabilizer. Here that warning produces a concrete failure. Two fixed vertices migrate into the 3-orbit's false-twin class, balancing it with the 5-orbit and creating a new odd part swap. Prime strong tournament modules cannot perform the analogous migration, because they are intrinsic and unique. This asymmetric behavior is exactly the one-unit discrepancy.

The correct survivors are the p-group refinement, the ambient Sylow refinement, and normalizer-relative signed identities. A single bijection preserving the presence of every odd-order symmetry subgroup is impossible, because the two finite fibres in (1) have different sizes.
