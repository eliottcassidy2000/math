# Independent correspondence audit

Result: PASS, with the stated observer/quotient scopes. No source files changed.

## Proof audit

1. Complement and observer reversal each reverse all tournament arcs. Their correlated product fixes every arc. This is a precise action identity, not a claim that either graph or tournament is self-converse.
2. Every nonconstant edge-bit assignment on a complete graph with n>=3 has a vertex incident to differently colored edges. For the chosen mixed triple, moving c across adjacent a,b reverses its directed two-step path. The two adjacent-swap cyclic-count changes therefore differ by exactly ±2; all remaining vertices' contributions agree. This proves the four-order witness and the full observer-stability classification, independently of finite census. At n=2 the classification is vacuous apart from the empty/complete graphs, as required.
3. The ordinary and affine relabeling actions differ, for the first adjacent transposition, by exactly one coordinate flip. A function invariant under both actions is invariant under that flip; conjugation by graph permutations supplies every edge flip. Their translations act transitively on the cube, proving constancy. Complement/converse can only coarsen these equivalences.
4. A disconnected graph–tournament incidence relation would give a nonconstant component-label function invariant under both actions, contradicting item3. The observer-stability theorem proves the unique degree-one left class is the empty/complete pair. Neither assertion identifies the bipartite incidence graph with the user's native metagraph.
5. Full Zeckendorf digit strings order vertices exactly: at the highest differing digit F_k, every legal lower tail is at most F_k−1. Finite tags do not retain that observer on an unrestricted infinite address set: three identically tagged addresses support the one-edge order-three hostile. This does not forbid finite tags from identifying a particular finite inventory, and does not weaken independent electrical guards.

## Independent computations

All five source-manifest files were verified byte-for-byte by SHA256 and size. The original checker was replayed normally and with `-O`; both outputs equal the frozen source receipts. A separate direct adjacency-matrix canonicalizer, not the source affine-mask implementation, exhausts n=2..5 and recovers counts:

- graph classes: 2,4,11,34;
- graph-complement classes: 1,2,6,18;
- tournament classes: 1,2,4,12;
- tournament-converse classes: 1,2,3,10;
- quotient-incidence edges: 1,3,14,120.

A separate direct three-vertex orientation test verifies the four local ±2 differences. The source's 1,090 nontrivial labeled-graph witness checks and full incidence connectedness controls also pass under both optimization modes.

These are finite controls of independently reviewed all-size proofs. They do not prove a new unlabeled bijection, resolve the Hamiltonian spectrum, or establish all-size native metagraph wiring.
