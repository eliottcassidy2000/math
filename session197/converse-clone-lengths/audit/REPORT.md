# Independent review: finite/all converse clone lengths and stabilized symmetry

**PASS, no correction required.** The canonical transitive-module argument gives an all-order finite/all decision for self-converse clone lengths, the sharp threshold and the stated automorphism stabilization. The unit-clock consequences and positive-binary classification use their named reviewed algebraic dependencies at the correct scope.

## Frozen source

Directory: `session197-orbit-wildcard/converse-clone-lengths/`.

- PROOF.md: `c3516df7522a100e692665f2d33e2427ff76a1dc83d0e873f072bbff6c66d890`.
- verify.cpp: `cb14a6b38cc5d2688a5c11ad0a20f151f7e97c0fa44ed19ad84f011148a7c5f9`.
- Both frozen receipts and copied optimized/unoptimized reruns: `923653d2b805f90c8e2b21cc1dc8c56ebcf3855abf604403a798ed531941af1c`.

Original files remain unchanged. Source and independent optimized/debug runs all pass. The proof dependencies are the general unit-clock reduction, the marked binary classification (SHA `65102468958b198f8538022288238c40d4e3ca33684cf0d419bc86b0993681a9`) and the marked negative-base rigidity theorem (SHA `bea7513641935a9ad0b5af15551dac626ec981667bddc62dbae054bed3339b9d`). The latter retains the separate nonattainment disclaimer at amplitude242.

## Mathematical audit

1. **Maximal transitive blocks are canonical.** Two intersecting modules have a modular union, and a triangle in the union would either lie inside one transitive module or use two vertices uniformly seen by a third. Both are impossible. Thus intersecting transitive modules have a transitive union; maximal blocks partition the vertices. This includes nested blocks and singleton sources. A nontransitive quotient has at least three blocks.
2. **Cloning cannot create a new merge.** The old root block expands transitively at the root's old position. A competing transitive module intersecting it can first be united with the whole block, then contracted to a strictly larger old transitive module. A competing disjoint block contracts unchanged. Both contradict old maximality. This also justifies that after contraction of the whole root block, the marked quotient vertex has a singleton maximal transitive block.
3. **Weighted quotient maps lift bijectively.** Every graph isomorphism or anti-isomorphism preserves canonical blocks. On a transitive chain, the order-preserving map is unique, and the order-reversing map is unique. Consequently weighted quotient automorphisms and anti-isomorphisms lift to exactly all actual vertex maps, without an extra internal permutation factor.
4. **The finite/all criterion and cutoff are exact.** A root-fixed converse map of the contracted graph extends at every positive length. Outside the finite collision set, a converse map must preserve the root block because its weight is unique. At a collision, the exact weighted quotient test decides the case. In a nonpersistent case there are at least three nonempty blocks, giving w_j<=n−l−1 and m<=n−2l. This proves the n−1 threshold for n>=2; n=1 is separately transitive. The cyclic three-block examples attain the finite bound at every n>=4.
5. **A single small self-converse input is insufficient.** The stated order-four source has H_m=3·2^m−1 and is self-converse exactly at m=1,2. Its cycle count derivation includes the even two-cycle boundary correctly. It gives a concrete hostile to dropping persistence or the threshold.
6. **The stabilizer is the right group.** Any automorphism fixing one vertex of a canonical transitive block preserves that block and then fixes its whole chain pointwise. Thus the original root stabilizer extends to every clone graph and equals the stabilizer of any vertex in the expanded block. Outside collision lengths it is the full automorphism group. This is an identification of actions, not only a divisibility statement.
7. **Unit clocks force trivial root stabilizer.** The named factorial theorem makes b and A integral and gives b−1 dividing A. Consecutive positive clone counts are coprime. The fixed stabilizer acts freely on each Hamilton-word set, so its order divides both counts and is1. Full rigidity follows at every noncollision length. The claim does not assert that all isolated early outputs are rigid or that outputs from different sources are distinct.
8. **Amplitude renormalization is correct.** If the old root block has length l, the contracted clone parameter is j=m+l−1 and its amplitude is A/b^(l−1). This is initially only an eventual identity on positive lengths. Finite-exponential uniqueness is applied afterward, so no deleted or negative length is silently introduced. The marked negative-base theorem then gives exactly the stated b=3, odd-square amplitude and lower bound.
9. **The positive-binary j=0 step is decisive and valid.** The marked classification gives a single-backward-edge core whose internal transitive block has size2j+1 and contains the root. The contracted root block is already singleton, so j=0. Equal outside padding is preserved. Restoring the old chain yields exactly Tr_d→C3[Tr_l,K1,K1]→Tr_d, with the marked root at any chain position and A=2^(l−1).
10. **All rooted and labelled counts are correct.** For l>1 the source is rigid and the l chain positions are distinct rooted classes. For l=1 the central triangle has three automorphisms, its three possible roots form one orbit, and each rooted stabilizer is trivial. Thus the rooted labelled count is l·n!, the prescribed-root-label count l·(n−1)!, and the unmarked labelled count is n!/3 for l=1 and n! otherwise. Equal padding does not create automorphisms between distinct SCC positions. The root-fixed anti criterion is exactly an odd chain with its central root. The modulo-seven and necessary power-of-two exponent claims follow as stated, without an infinite-prime claim.

## Independent controls

The standalone independent `check.cpp` imports no source code. It identifies transitive modules by explicit triangle exclusion and uniform outside incidences, then forms maximal blocks by union-find over all such subsets. It separately checks the union lemma. It enumerates all degree-compatible permutations by Cartesian products of degree groups, rather than using the source's incremental map backtracking.

It regenerates all 5,405 labelled rooted tournaments through order five and all 32,110 corresponding clone graphs through m=n+1. For every clone graph it constructs the complete sets of actual automorphisms and anti-isomorphisms and compares them **map for map** with the weighted quotient lifts. There are 64,220 exact map-fibre comparisons. It checks all root stabilizers, every noncollision full group, all finite/all decisions, the strict cutoff, and all 5,405 two-clone block partitions.

The independent inventory exactly matches the source:

- 1,853 persistent rooted sources.
- 3,552 finite cases.
- Finite-length set sizes: zero for1,664 cases, one for1,600, two for288.

This finite histogram is not asserted to bound the possible cardinality at larger orders.

Additional independent positive-binary controls use l=1,...,6 and equal padding d=0,1,2. They check 189 actual clone Hamilton counts with arbitrary chain-root positions, all full automorphism groups, 63 rooted isomorphism classes, and the central-root criterion. The source separately replays its sharp-cutoff and negative-binary hostile controls.

Independent checker SHA: `1b8556ca5d287eb5166a916b0536c3f361d0a2aa496537386a03d2c12124afae`.
Independent optimized/debug receipt SHA: `af2465860174e834df38216607354a77ff0bd1ccc7a0bf88108a1ead0698c1c9`.

## Separate optional consequence

The author's proposed, not-yet-included converse-orbit corollary is also correct: for fixed l,d the l rooted classes in the persistent positive-binary family form ceil(l/2) orbits under converse. Reversal pairs chain positions k and l+1−k. There is one fixed rooted class for odd l and none for even l. The l=1 triangle case is already a single rooted class and satisfies the same statement. The independent checker verifies these orbit counts in the padded controls, but this paragraph does not modify the frozen source theorem.

## Reproduction and scope

Compile either standalone C++17 checker normally and in an assertion-enabled unoptimized build; run with no arguments. Both output deterministic receipts. The independent checker needs no bank or external file.

The all-order claims concern one fixed root repeatedly replaced by a nonempty transitive chain. No conclusion about adaptive graph reconstruction, a single arbitrary self-converse source, or unrestricted Hamilton-path spectrum closure is inferred.
