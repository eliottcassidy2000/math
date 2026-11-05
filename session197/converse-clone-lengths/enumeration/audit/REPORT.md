# Independent review: exact binary-family orbit enumeration and sampling measures

**PASS, no correction required.** The four class generating functions, automorphism mass, exact normalization fibres and conditional sampling discrepancy all follow at every order for the classified persistent self-converse positive-binary family. They are not statements about all tournament classes.

## Frozen source

Directory: `session197-orbit-wildcard/binary-orbit-enumeration/`.

- PROOF.md: `a7bddbae7fa55a44bd0f847b477934a72f82428a313ee24681f9335b1e056d19`.
- verify.cpp: `f34d61e189ebf603bb15029b8aa116d7efbd7efbb6d9b200f4677e490940d6a8`.
- Both source receipts and copied optimized/debug reruns: `2897101be769b4b981c38ee10580592790bd3594e108e2330b6089581666a622`.

Original files remain unchanged. Named dependencies are the reviewed converse-clone-length classification (SHA `c3516df7522a100e692665f2d33e2427ff76a1dc83d0e873f072bbff6c66d890`) and converse mass/fibre theorem (SHA `07f4b20dde9f9d23b7711f65bd484d7468bdbef3cc6da0cc2888aa91d402a620`). Their roles are kept distinct: the first establishes completeness of the physical family; the second identifies the normalization fibre, not the ordinary class weight.

## All-order proof audit

The pair (l,d) determines the unique nontrivial SCC and equal outside padding, so distinct pairs give different unmarked classes. Summing weights 1, l, ceil(l/2), and [l odd] against z^(l+2+2d) gives exactly U,R,C,F. Converse reverses chain positions. The identity R+F=2C is the correct order-two orbit relation. The stated floor and triangular-number coefficient formulas follow separately for even and odd n.

At l=1 the central triangle has three automorphisms, its three possible marked roots form one orbit, and each root stabilizer is trivial. At l>1 the graph is rigid. Hence the unmarked mass correction is exactly2/3 on odd orders n>=3. All sequences are explicitly zero below order three, avoiding a spurious negative mass at n=1. The labelled rooted, prescribed-root-label and unmarked labelled formulas use the appropriate different stabilizers.

For the fixed-word formula, there is one increasing Hamilton word. Every other word uses the unique backward edge and is determined by a subset S of the internal l-chain placed before that edge, with the complement afterward. Reversal followed by reflection sends S to the reflected complement. Fixed subsets therefore choose one member of every mirrored pair; a fixed central internal vertex makes this impossible when l is odd. This proves exactly1 fixed word for odd l, and2^(l/2)+1 for even l. The source's equivalent central-edge argument is valid, and padding is forced on both ends.

The anti-involution centralizer is trivial. For l>1 this follows from rigidity. For l=1 reflection conjugates each triangle rotation to its inverse, so only the identity commutes. The reviewed fibre theorem then gives exactly the fixed-word count per unmarked class, without an additional Aut factor or a multiplicity from the other conjugate anti-involutions.

Summing the fibres gives S_(2s+1)=s and S_(2s)=2^s+s−3 with the stated domains. The displayed rational generating function agrees with the independent decomposition

    S(z)=U(z)+2z^4/[(1−2z^2)(1−z^2)].

The total normalized universe has the named theorem's size, so division by it gives the unconditional probability of representing this family. This is correctly distinguished from uniform unmarked-class sampling.

Conditioned on this family, exactly d=0 is strong. At even order the strong class has fibre2^(s−1)+1, while the total is2^s+s−3. Thus the two conditional probabilities and limits0 versus1/2 are exact. At odd order all fibres are1, so both measures give1/s. The conclusion is a genuine sampling discrepancy for this family and does not claim a universal distributional law.

## Independent controls

The standalone independent `check.cpp` imports no source code. It constructs all 42 family graphs through order14. Unlike the source's full Hamilton-word DFS, it counts all Hamilton paths by an endpoint subset DP and enumerates fixed words directly from one chosen vertex in each reflected pair, followed by the uniquely determined reflected half. It handles the odd central fixed vertex separately.

Across these graphs it verifies:

- The full Hamilton-count sum10,950.
- Exactly282 fixed words, from1,926 admissible first-half path candidates.
- Exact sets of normalized adjacency matrices, with no duplicates inside a class and disjoint sets across different classes of the same order.
- Every normalized graph's consecutive Hamilton edges and reflection anti-isomorphism.
- Full automorphism groups by a separate incremental adjacency backtracking algorithm, their anti centralizers, and all rooted/converse orbit counts.
- Every per-order class count, rooted count, mass and normalized-fibre total.

The independent optimized and assertion-enabled unoptimized outputs are identical.

- check.cpp SHA: `d4682b0cdbaaf7244ee794e2770b7d38e4a9ef1476a5a870823d027b500d62d5`.
- independent-normal.txt and independent-debug.txt SHA: `9b9baf256c64c6513f7492f30258c1baed070cf8c54e0d6def8e48267bfb55ac`.

Compile either standalone checker as C++17 and run without arguments. There are no input-file dependencies. The finite controls supplement the explicit generating-function and subset/reflection proofs; they are not extrapolated to establish the all-order claims.

Conclusion: all claimed enumeration formulas and the measure-sensitive sampling example pass, with the stated nonconstant-family and order boundaries.
