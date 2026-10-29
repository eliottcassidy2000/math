# Session90 independent mathematical and computational audit

**Verdict: PASS in the scopes below.** The independently implemented checks run identically in normal Python and `python3 -O`. There are **1,993,328 explicit runtime checks per mode** across nine audit programs. This does not mean nine independent proofs of every claim; dependencies and finite universes are stated separately. Primary-literature priority is not claimed.

Source repository snapshot: `3eb49fa63d7ff4e58742b1ae84944f832e87a118`.

## 1. Definitions and corrected scope

Royle sign is the parity of reversed **endpoints** of graph edges under an automorphism. It is not vertex-degree parity or the sign of the induced edge permutation. K2 already separates these definitions. Decoration changes the stabilizer on which this character is tested. A graph which is odd without marks can be relatively even after marking.

The tree-boundary map concerns completions of a retained oriented tree. It is not a bijection between all unmarked tournaments of consecutive sizes. An oriented tree can have an even automorphism group: an outward three-leaf star has group S3. Only the stabilizer after tournament completion is necessarily odd.

The aligned-pair/self-converse class bijection changes stabilizers to anti-involution centralizers. Marked mass, ordinary mass, and uniformly labelled sampling are different targets. The singleton uses a signed anti-involution whose underlying permutation is identity.

## 2. Decorated Walsh correspondence and integral lattice

The arbitrary finite-decoration theorem is correct: partial Walsh transformation gives signed orbit lines, and precisely the joint-stabilizer-even lines survive. No averaging division by the group order is needed. The Sylow normalizer/Frattini proof extends because the decoration contributes the same fixed-point multiplier to both traces. For n>=5, a globally fixed decoration supplies the Sylow3/5 Hall obstruction. That hypothesis is sufficient, not a classification of every possible decoration set.

`independent_decorated_counts.py` exhausts the labelled universes through n=4, with an independent action implementation. It obtains pair and relative-even graph counts

`1,2,8,88,2752`.

Requiring the underlying graph itself to be even incorrectly gives

`1,2,6,60,1092`.

`independent_sylow_decorations.py` uses literal S5 actions and independently checks full-stabilizer versus normalizer evenness. The source Sylow3/5 families have512 and8 classes and are disjoint. Their target families intersect in exactly4 classes, so the target union has516 against520 source classes. This is one exact Hall deficit, not the full bipartite matching deficiency.

`independent_fourier_lattice.py` builds the complete matrix by signed graph-orbit summation and builds its inverse numerator directly from labelled Fourier coefficients. It verifies both products `MN=NM=2^m I`, exact Bareiss determinants, and the weighted group-product formula. For n=0,1,2,3 the determinant valuations are0,0,2,116 and the cokernel exponents are1,1,2,8. The proof that all Smith factors divide2^m is valid. In a general affine binary action, the source 2-adic stabilizer term must be retained; its omission is valid only in the tournament case.

## 3. Tree-boundary compression and sparse graph carriers

The first/last-edge map is a genuine natural bijection. Its inverse chooses the outer endpoint of each edge relative to the other, so it retains ordered-pair reversal information. The induced completion codec and the marked stabilizer equality follow. Faithfulness of an oriented tree's edge action holds, including the two-vertex case where its oriented-edge stabilizer is trivial.

`independent_tree_boundary.py` checks every labelled tree through seven vertices, totaling18,249 trees, with reversal, inverse, and relabeling naturality. It separately reconstructs every labelled tournament for every tree through four vertices,1,051 round trips. This uses Prüfer enumeration and deletion-component inverse, independently of the producer's distance-based implementation.

For an arbitrary n−1-edge constraint graph, forest-plus-bridges supplies c−1 free bridge orientations and c−1 transported chord constraints. Their counts cancel, but their operations do not. The retained forest/bridge choice is essential for equivariance.

The sparse encoder proof is valid. Its high-degree core vertices and distinct path-length tags recover the original roles intrinsically. The added caterpillar trees have maximum degree3, so they cannot be confused with the tagged core. Leaf pruning leaves a path whose reflection is destroyed by the forced near-end leaf and the absent opposite leaf. Thus padding components are rigid and pairwise nonisomorphic.

`independent_sparse_audit.py` imports the encoder as its input, but uses its own structural decoder and its own center/rooted-tree automorphism recursion. It recovers all source graph/tournament classes through order4 and independently confirms8,256 padding trees across every padding family k<=128 are rigid and pairwise distinct.

`independent_sign_codec.py` checks all1,100 labelled simple graphs through order5 and all4,374 lifted automorphisms. The generic cone/tag encoder preserves the **entire Royle sign character**, including1,411 negative-sign maps. The all-order proof orients new hub/tag edges canonically and equivariantly; only the original graph edges contribute a character. Padding is fixed pointwise. Exact automorphism recovery then gives the claimed evenness equivalence.

The generic graph reduction is GI-complete for the disconnected e=v−1 target restriction. The tournament-only reduction is a faithful tournament-isomorphism reduction; it does not itself establish GI-hardness. Neither gives a polynomial-time solver for either isomorphism problem.

## 4. Capacity and Euler-repair lower bounds

The capacity proof is sound. From the labelled tournament lower bound and labelled sparse-graph upper bound, fixed k same-size outputs fail eventually. The finite13/28/43 thresholds are first failures of those **loose inequalities** through500, not first actual class-count failures. Indeed at n=3 there is only one unmarked two-edge simple graph class but two tournament classes, so any fixed-length tuple of such same-size graph classes already fails there.

For output order N, `k(N−1) log2(eN/2) >= binom(n,2)−log2(n!)` gives

`liminf N log2(n)/n² >= 1/(4k)`.

The variable-total-size bound is valid for fixed tuple length. Marks or external bit strings cannot be free extra storage.

The Euler identity yields t=1+sum(beta_i−1) over cyclic components, so retaining a connected rank-beta core requires at least beta additional tree components. Odd automorphism order forces every tree component to be rigid and every component type distinct. Counting tree types by rooted plane trees gives total size at least `sum_i log4(3i+1)`, hence Omega(beta log beta). The construction matches this order only within the stated disjoint core-preserving odd-symmetry contract. No universal lower bound against all encodings is implied.

## 5. Anti-involution transport and its infinite family

The unique odd-order square root is correct, including in the ambient permutation group: an odd-order square root is a power of its square, and cyclewise rooting uses only odd cycles. The conjugator `sqrt(ts)` is exact and equivariant. The polar factorization A=KC is a unique set factorization; K need not be a subgroup. Flat pair-local transport is equivalent to ambient commutativity of K. A chosen global base supplies a flat alternative even when pair-local transports are nonflat.

`independent_holonomy.py` independently constructs C3, Paley7 and C3[C3,C3,C3], and searches the full symmetric groups to confirm their complete automorphism groups. Every triangle of the27 anti-involutions of the nine-vertex example is checked. Exactly11,664 triangles are nonflat. The fixed-base holonomies lie in its C3 centralizer; when all bases vary, seven different permutations appear across the conjugate centralizers. The script independently verifies the order27 twisted operation's Latin laws, automorphic inverse law and **left Bol** convention, and finds nonassociativity.

`independent_small_selfconverse.py` does **not** use the pair lift or its action. It builds ordinary tournament banks using previously audited, SHA256-pinned S87/S88 code, then searches literal permutations for every automorphism and every anti-involution. The ordinary banks have1,2,12,456 classes at orders1,3,5,7; filtering yields1,2,8,88 self-converse classes. All99 are flat. Together with the independently checked order9 witness, this validates minimum odd nonflat order9 without sharing the current pair-action implementation.

The iterated-C3 proof also passes inspection. Same-top-block neighborhood distances are even and at most s−1; cross-block distance is2s. This intrinsically recovers the full wreath group. The centralizer recurrence and formulas A=N C² and i=N C follow. There are exactly C anti-involutions fixing each named vertex. The zero-prefix embedding of the order9 witness stays anti across the module boundary and proves nonflatness at every level k>=2. The growing fixed-root mark information bound is for choosing an individual actual mark, not for describing the family or answering a mark-insensitive decision.

## 6. Coprime descent, path observer, and weighted quotient

The p-group fixed-orbit proof is correct by prime-to-p orbit and transporter cardinalities. The explicit involution transporter repair is valid. Its application must retain the joint relational object before descent; separate marginal classes do not supply their relative alignment.

The folded Hamilton-path recurrence is correct and intentionally computes h_tau rather than H. Internal mirror-pair arcs are never queried. The mean `E[h_tau]=n!` follows because each candidate half-word imposes n independent fair constraints on distinct matrix-coordinate blocks.

`independent_odd_descent.py` supplies a separate nonabelian control at the nine-vertex wreath tournament. All512 subsets give24 ordinary orbits,32 fixed subsets and16 fixed orbit classes, equal to16 centralizer orbits. It checks432 repaired transporters. Literal9! Hamilton words yield H=3159, h_tau=33 and11 blue presentations.

The same audit implements the centralizer state quotient independently and compares **each orbit's intermediate total with literal partial-word enumeration**, not another recurrence. It confirms64 states→24 orbits and96 arcs→32 weighted entries, with4 singleton and20 size-three orbits. Initializing each orbit with mass1 incorrectly returns19. Initializing with orbit size returns33.

The proof of the weighted quotient is sound even for a supplied subgroup H<=C, but the final blue normalization must still use the **full** |C|. Division by |H| instead counts H-orbits. This clarification was requested in the producer note.

The marked sampler is correctly normalized. After acceptance probability1/i(W) and uniform vertex relabeling, each labelled W has the same unnormalized probability `n! / (2^(n²)(2n+1)!)`. Finding i(W) remains a charged input cost.

## 7. Reproduction and dependency inventory

Run each command with both `python3` and `python3 -O`; JSON outputs agree byte-for-byte:

| Script | Explicit checks | Dependency |
|---|---:|---|
| independent_decorated_counts.py |3614|standalone|
| independent_sylow_decorations.py |2085|independent_decorated_counts.py|
| independent_fourier_lattice.py |239079|independent_decorated_counts.py|
| independent_tree_boundary.py |1614504|standalone|
| independent_sparse_audit.py |80397|producer encoder as input; independent decoder/AHU|
| independent_sign_codec.py |8748|producer encoder through independent_sparse_audit.py|
| independent_holonomy.py |42119|standalone|
| independent_small_selfconverse.py |1750|SHA-pinned S87/S88 ordinary generator; independent_holonomy.py permutation primitives|
| independent_odd_descent.py |1032|independent_holonomy.py permutation primitives|

`independent_small_selfconverse.py` finds the repository above its path, or accepts `--baseline-root REPO`. `independent_sparse_audit.py` defaults to the sibling `trijection/check_sparse_carrier.py` and accepts a positional producer path. No network, optional numerical package, or assertion-elision dependency is used.

Generated JSON receipts are optional local outputs, not required checked-in fixtures. The audit does not independently compute Monster/Baby Monster groups, prove external priority, or resolve any claimed major open problem. Those literature bridges remain exactly as scoped in the source-linked notes.
