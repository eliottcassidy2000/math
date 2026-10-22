# Independent Session 26 cut, classification, and spectral review

Status: PASS. Reviewed 2026-10-22. Reviewer: independent audit_cut_partition_spectrum lane. No repository writes were made by this reviewer.

## Exact sources bound by this review

- `session26-tournament/APEX_CLASSIFICATION_AND_SPECTRAL_LAW.md`
  SHA256: `82d01a9acb14fab5d9f53acffc8202ce81beb2f7b55c4dbc8826a291f3a3e9f7`
- `session26-tournament/CUT_ADDENDUM.md`
  SHA256: `05d8be14ff381ff5dc403ca2d2d50ca5199b9d34d2c1dcd832e9cd803bebe244`

Both complete files were read directly, and their hashes were computed from the files. This review binds these source versions, not merely the earlier mathematical messages. Session 25's modular-completeness dependency is explicitly inherited, not independently re-audited here. The separate observer proof and H13 witness are outside this review's file-binding scope.

## Mathematical findings

The cut/Ferrers identity, hook-distance interpretation, minimum and maximum values and equality cases, source/sink stripping, exclusion of H=13, and dyadic orbit bound are correct. The special low-length weight plateau is handled correctly. The Ferrers statistic is not misidentified with the tableau hook-length formula, and the equal-H shapes are not incorrectly asserted to be different unrooted tournament classes.

The deletable-apex classification is correct. Triangle intersection characterizes the deletable vertices. With two apices, their insertion cuts into the common remaining chain give precisely the triangles between those cuts; strong connectivity forces the cyclic blow-up C3(1,1,n-2). The three-apex strong exception is C3. Trivial apex stabilizers and odd automorphism-group order imply rigidity in strong orders at least four. The ordinary and self-converse class formulas then follow with the exceptional double-apex multiplicity correctly subtracted.

The binary product identity is exact, including the adjacent-10 correction and backward-pair remainder. Each backward geometric sum is less than one half. The leading-prefix coupling uses disjoint source coordinates, while correctly allowing dependence between the full binary integers. Bounded product density yields the stated uniform CDF error. Both CDF branches, variance 55/2304, and finite mean 9*2^(m-4)+1/2 for m at least three check out.

The exact W(X) decomposition and cancellation identity are correct, with valid length cutoffs. They yield W(X)=(log 2)X+O(sqrt(X)). The ordinary-core duplicate correction, self-converse bound, and converse-merged count are correct. Independently, reverse-complement word orbits suffice for the distinct-value bound without requiring class enumeration. The lower density of missing odds is a lower bound, not an asserted natural-density limit.

The optional Dirichlet consequence is valid: partial summation gives the main term c*s/(s-1), and the counting remainder integral is holomorphic for Re(s)>1/2. It concerns presentation-weighted or class-weighted strong-core spectra, not distinct H values.

## Essential scope gate

PASS: the final classification/spectral source explicitly limits all-order finite class-counting asymptotics and Dirichlet series to strong cores. Arbitrary transitive padding would produce infinitely many classes at fixed realized H. The distinct-value conclusion nevertheless applies to all feedback-vertex-at-most-one tournaments, because padding preserves H. The final source states this distinction clearly and repeatedly.

## Independent finite controls

Pure JavaScript checks conducted independently during the shell outage exhausted all strong words of lengths 2 through 12. They checked the exact binary-product identity, backward-remainder bound, minimum sequence, absence of H=13, and the dyadic excluded family. All checks passed.

A separate canonical-isomorphism enumeration exhausted pointed words for tournament orders 3 through 7, canonicalizing adjacency matrices over permutations within degree classes. Ordinary class counts were 2,4,9,21,48; strong counts 1,1,3,7,15; all self-converse counts 2,2,5,3,10; strong self-converse counts 1,1,3,1,5. Every strong order at least four had one or two deletable apices, exactly one class had two, and its pointed-word multiplicity equaled the apex count. All matched the symbolic formulas.

These controls support the proofs but are not substitutes for the all-order arguments. No literature-priority claim is assessed or made. No full tournament-spectrum closure or general metagraph statement follows.
