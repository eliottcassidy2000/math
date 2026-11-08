# Independent adversarial audit: Session227 physical-state and graft images

Verdict: **PASS** for the exact final inputs bound by the hashes below. No unresolved mathematical defect was found. This audit does not claim Chan–Pak completeness, scalar coverage, or literature priority. It did not edit the producer's files or publish anything.

## 1. Proof review and exact boundaries

### Physical closed cone

The state convention is exactly Session225's `(D,E,U,d)`, with `H(empty)=1`, ordered two-path bins, one empty unmarked path, and a nonempty rooted bin beginning at the root. The four inequalities are valid for all nonempty rooted tournaments, including singleton and sink-root sources:

- `E>=0` is counting positivity.
- Deleting a terminal root is an injection from the E paths into deletion paths, giving `E<=d`.
- Inserting the root at a deterministically selected legal position injects deletion paths into full paths, giving `d<=D`. The empty deletion at a singleton has one insertion.
- Splitting every full Hamilton path not ending at the root immediately before the root injects it into two-bin covers with nonsingleton rooted bin. Exactly d covers have singleton rooted bin, giving `U-d>=D-E`.

The proposed rays form a unimodular basis. Their unique coefficients are precisely `(E,d-E,D-d,U-D-d+E)`. Singleton and source edge realize the first two rays. The marked cyclic-block family has states `(2^m+1,1,2^m+1,1)`, whose normalization tends to the third ray. Source-rooted transitive tournaments have states `(1,0,2^(n-1),1)` and tend projectively to the fourth ray. Consequently the closed conic hull is exactly the four-facet cone, and its dual consists of the nonnegative combinations of these facets.

An explanatory gap in the initial draft was repaired: being nonphysical by itself does not exclude membership in the finite conic hull. The required argument is that every nonzero finite nonnegative combination of physical states has `d>0`, whereas rays three and four have `d=0`. Thus closure really is necessary. This is a clarification of a true result, not a counterexample to it.

The result concerns homogeneous inequalities over all orders. It does not classify fixed-order, affine, nonlinear, or discrete constraints.

### Affine integer lattice and nonrealizability

Rédei parity supplies odd D and d. The singleton-graft count `D-E+U` is also odd, so `U=E mod 2`. The five stated physical witnesses give a difference determinant of magnitude eight and generate `(2,0,0,0)`, `(0,-1,1,0)`, `(0,0,2,0)`, `(0,0,0,2)`. These are exactly the difference lattice of the stated parity coset. Thus the affine integer hull and its index eight are correct.

The hostile `(1,0,6,1)` satisfies both exact hull descriptions but is not physical. Insertion monotonicity makes every tournament containing a directed triangle have at least three Hamilton paths. A tournament without a directed triangle is transitive. For a transitive rooted tournament the rooted cover count is exactly `2^outdegree(root)`, not six. This establishes all-order nonrealizability, rather than a bounded-search failure.

### Cyclic-block gadgets and determinant signs

For `X -> u -> v -> X`, with X a nonempty transitive m-block, the Hamilton paths divide into the `2^m` choices with consecutive `u,v`, and the single path beginning at v and ending at u. Thus `H=2^m+1`.

For two bins, placing u and v in the same bin contributes `2(3^m+2^m-1)`, and placing them in separate bins contributes `2^(m+1)`. The subtraction at zero X-members in the joint bin is essential: that induced graph is a directed edge with H=1, not the extension of `2^m+1` to m=0. The draft handles this boundary correctly. The total is `Q=2*3^m+2^(m+2)-2`.

For either singleton root, `s=e=q=1`, `b=2^(m+1)`, and

    Delta = 2(3^m+2^m-1-4^m).

It vanishes at m=1 and is negative for all m>=2 by the displayed induction. Marking u gives `U=1+2^m`; marking v is not claimed to have that same U.

For any X root, deletion gives the smaller cyclic block when m>=2. Substitution yields

    Delta = (2^m+4)3^(m-1)+3*2^m > 0.

The endpoint counts at index i are `s=2^(m-1-i)`, `e=2^i`. All three blocks are nonempty, so the graph is strong. Combining with Session226's proved determinant factorization gives full rank of either sign at every order >=4, with the triangle boundary singular. No all-order inference is being made from a finite gadget census.

### Transported hulls, exact inverse, and scalar congruence

For invertible K, linear homeomorphism transports the closed cone exactly. Integer-affine combinations transport the affine lattice exactly. The indices are `8*abs(det K)` and, for a fixed chain, eight times the product of determinant magnitudes. A specified chain has actual graph realizations for every physical source; no formal source vector is silently substituted for a graph.

The four displayed B4 inverse formulas are correct on the whole rational space. For integer output `(A,B,C,F)`, integrality of the recovered d forces `A+3F=0 mod 8`; this condition also clears every other denominator. Hence the integer image lattice has index eight. Applying the recovered source parity conditions gives physical affine-hull index 64. The vector `(145,31,159,37)=K(1,0,6,1)` is in both transported hulls but not in the physical graft-state image by injectivity. Its own realizability as some unrelated tournament state is deliberately not asserted.

The scalar-projection corollary is also correct. The first row maps the affine basepoint to H. Its increments on the four difference generators are `2H,Q,2(Q-H),2(2H-Q)`, with gcd `gcd(2H,Q)`. Thus the scalar affine hull is exactly `H+gcd(2H,Q) Z`. This is the strongest affine congruence, not a coverage statement. Odd m>=3 in the cyclic-block family gives H and Q divisible by three while the matrix is invertible, so these full-rank gadgets still force scalar divisibility by three. The first-row argument itself does not require invertibility, although the surrounding image-cone discussion does.

### Graph-level certificate and explicit collision

The module certificate is genuinely an equivalence. Necessity follows directly from the graft's cross-arcs. Conversely, deleting the candidate module and marking its unique externally dominated vertex reproduces every induced arc and cross-arc. Properness/nonempty complement, root containment, rooted gadget isomorphism, and the unique dominated exterior vertex are all retained. It works for singleton and nonstrong gadgets, does not require an invertible matrix, and does not claim module uniqueness. Repeating the certificate in reverse order recognizes a prescribed chain only when a full successful module sequence is retained.

The two supplied six-vertex graphs are valid strong tournaments with state `(37,7,43,9)`. The good graph has exactly the claimed rooted B4 certificate, module `{0,2,3,4}` and old root 1; its complement is the source-rooted edge `1 -> 5`. The bad graph has no four-vertex module containing 0 at all. Our independent receipt supplies a mixed-orientation exterior witness for each of the ten candidate subsets, before any gadget-isomorphism test. Their common numerical inverse `(1,0,2,1)` is physical. This therefore separates numerical inversion, physical source-state existence, and actual rooted graph-image membership.

The initial awkward sentence about an absent candidate was replaced by the precise stronger claim above. No graph minimality is claimed or proved.

## 2. Inheritance and novelty audit

Read and checked the Session225 graft proof, the Session226 determinant/inverse proof and its independent audit, and Session226 root-descent Section 1. Their conventions, singleton exception, and logical boundaries are preserved.

The earlier `adaptive_metagraph_reformulation_20261103/ADAPTIVE.md` Section 3 and `adaptive_state_induction_20261103/BLOCK_SWITCH.md` B2 were read through the repository's committed objects. They already give the orientation-dual endpoint-defect inequality and positive four-coordinate system. The new draft correctly treats these inequalities and observer dimensions as inherited, rather than claiming priority. `nonzero_preparation_20261103/algebra/ONE_EXIT_AND_PERSISTENCE.md` Section 5 already states the cyclic transitive-block scalar count `1+2^m`; the new draft acknowledges this as well.

The older ADAPTIVE Section 4 H37 hostile uses different rooted states and separates arbitrary root switching. It does not assert the graph-image module obstruction proved here. Rédei parity is catalogued by `01-canon/theorems/THM-001-redei.md` and the current core-paper route. No literature-priority search or assertion is made.

## 3. Independent exact computation

`independent_check.py` imports no producer code. It generates tournaments as directed edge sets, recursively enumerates every simple directed path, assembles ordered path-cover counts over subsets, and computes determinants by the Leibniz formula. This is independent of the producer's forward subset DP and elimination routines. Every check uses an explicit exception, so optimization cannot disable it.

Coverage:

- All 5,405 labeled rooted presentations at orders 1..5, with every root and no graph filters: cone inequalities, affine parity, the H=1 power-of-two boundary, determinant factorization, and scalar affine congruence generators.
- All roots of transitive tournaments at orders 1..10.
- Every root of the cyclic family for m=1..9: H, Q, endpoint/deletion counts, Delta formulas/signs, determinant factorization, strongness, ray-sequence state, and the odd-m divisibility-three barrier.
- Five actual affine-lattice witnesses, determinant eight, and four exact cone rays.
- 1,140 literal graft outputs and graph-certificate positive controls: every rooted gadget through order four on the four physical spanning sources.
- 625 exact rational inverse and integer-image lattice probes for B4.
- The transported nonphysical source, the actual strong state-collision pair, and all ten separate module-obstruction witnesses for the bad graph.

Reproduction, in this audit directory:

    python3 -I -B independent_check.py > independent.json
    python3 -I -B -O independent_check.py > independent.optimized.json
    cmp independent.json independent.optimized.json

Both isolated runs passed all 32,425 explicit checks with byte-identical JSON. An AST inspection also confirmed that the checker contains no assert guards. Finite tests corroborate the proofs; they do not establish the universal claims on their own.

## 4. Producer replay and final binding

The final producer's normal and optimized isolated runs both passed. Their outputs are byte-identical to each other and to both final producer receipts. Its inherited verifier was compared byte-for-byte to the committed Session225 verifier and is unchanged.

Reproduction from this audit directory:

    python3 -I -B ../session227-graft-images/verify.py > producer.normal.json
    python3 -I -B -O ../session227-graft-images/verify.py > producer.optimized.json
    cmp producer.normal.json producer.optimized.json
    cmp producer.normal.json ../session227-graft-images/CHECK.json

The final producer counts are: 5,405 profiles, 75 rooted cyclic-family statistics, 18 direct-permutation family-statistics controls, 11 transitive-ray controls, 6,561 B4 integer-image grid points, 75 scalar-lattice gadget projections, 5,625 literal graft pairs, and 75 literal two-step chains. Its exhaustive rooted-at-zero graph-image census tests 1,024 order-five and 32,768 order-six tournaments, finding exactly 24 and 240 B4 images. Its recognizer and comparison image bank are built by separate module-test and direct-graft/relabeling paths.

All final mathematical edits, including the scalar affine-lattice corollary and precise verification universes, were inspected. No further mathematical correction is required. The producer's final proof header reads 'PROVED, independently audited with FINITE-EXACT controls.' The recorded PROOF.md hash below already binds that promoted version; no further status change is needed.

### Audited input hashes (SHA-256)


- PROOF.md: 6afd19677753eec9041f472e7d9a19791c13d1bbb564e38fb7923d7b3d0ffb10
- INHERITANCE.md: ca25b41930ba44b7355a26d0fe47eb1ce3b030b8432dab42e272c6f892931a56
- verify.py: ba4c47f6eb20b2f8b4eaff43e9d182dce3e1c55f8a78c9b15473fc000e407c51
- inherited_verify.py: 9e834801a3dcaeb74a22a9fd7787d2bdf8c6c1b49e594f7177f0031bdfe088d5
- CHECK.json: f1261e359f966e9f6fa0de406f7edbdc2869aa5acc3bad3e46437203d2a832d3
- CHECK.optimized.json: f1261e359f966e9f6fa0de406f7edbdc2869aa5acc3bad3e46437203d2a832d3

These six inputs are frozen in the separate `audited_inputs.sha256` manifest. Audit artifacts are `independent_check.py`, `independent.json`, `independent.optimized.json`, `producer.normal.json`, `producer.optimized.json`, this report, and the input manifest.
