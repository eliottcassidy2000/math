# Independent completeness and finite-census review for C=3

Status: **PASS**, no correction. The exact source snapshot below has been independently regenerated through order twelve, with every retained graph matched up to isomorphism. This is not extrapolation from a finite pattern: the new strong-extension completeness argument and the stated minimal-counterexample reduction are reviewed separately from the numerical inventory.

## 1. Frozen source binding

All source files are under `../session197-lowcycle/`:

- `C3_COMPLETE.md`: SHA256 `a919deb88bc78a84b8d0cf40199fc7966051541f1e7027be3666bb9f0713b444`.
- `verify_c3_complete.cpp`: SHA256 `e49f0b9ff6e29f7b7260d6516ccb0fafe03ea4bbcf06dccc0664b7fddead806f`.
- `verify_c3_complete.py`: SHA256 `88096aede3e208dfa66a45f02ae576d7459c8322841663b8627e630e24547799`.
- `c3_complete_certificate.json`: SHA256 `5461864f55486797ef825342841ef30e5258f90cc74b194087573b5fb2ce341c`.

All four source hashes were checked again after the audit. The Python wrapper was copied into an audit-only directory before replay, since it writes its certificate. Normal and Python `-O` runs reproduce the frozen certificate byte for byte. No source file was modified.

## 2. Independent enumeration and counting

The new checker `check_c3_complete.cpp` does not use the source's insertion lower-bound pruning or its capped Hamilton-cycle DFS.

For each retained parent A on n vertices, it computes the all-start Hamilton-path endpoint matrix M by bottom-up subset dynamic programming. M_(s,t) is the number of spanning directed paths from s to t, capped at four. For a new vertex whose out-neighbor set is O, every new Hamilton cycle, started at that vertex, corresponds uniquely to an old Hamilton path starting in O and ending outside O. Thus its count is

    sum_(s in O, t outside O) M_(s,t),

with the same exact distinction between zero, one, two, three and at least four. Capping individual nonnegative summands at four cannot change a total below four or turn a total at least four into one below four.

Every nonempty proper O is tested. No insertion-count filter is used. The implementation computes these directed cut sums by the recurrence

    cut(U+v)=cut(U)+sum_t M_(v,t)
              -sum_(u in U)(M_(v,u)+M_(u,v)).

The endpoint DP packs independent start counts into four-bit lanes, with an explicitly tested saturating addition. This is an implementation optimization, not a different mathematical assumption. For each parent its known cycle count is also checked from the endpoint matrix by summing closing endpoint pairs and dividing by n.

The independent complete generation used:

- 9,567 parent endpoint matrices.
- All 15,611,794 nonempty proper incidence masks across the generations.
- 212,712 retained labelled extensions.
- 27,945 final isomorphism classes across orders four through twelve.

Its complete C1/C2/C3 inventory agrees with every source entry. In particular the C3 counts at orders five through twelve are

    1,11,38,124,390,1193,3574,10537.

Every C3 representative has an indegree-one or outdegree-one vertex. The independent checker tests the disjunction; it does not incorrectly require both orientations.

## 3. Independent isomorphism handling

The source uses all actual Hamilton cycles and rotations as canonical coordinates. The new checker instead starts with outdegree classes, repeatedly refines them by their out-neighbor counts in every current class, and exhaustively permutes within every final class. The full adjacency encoding is minimized over those permutations.

The refinement is isomorphism-invariant and never merges old classes. Every isomorphism respects the final classes, while equality of full adjacency encodings supplies an actual isomorphism. This is exact canonicalization, not degree-signature deduplication. The maximum permutation workload for any graph in this entire audit was 120; the explicit safety ceiling of one million was never approached.

The source C++ verifier was run with an audit-only CSV output. Every one of its 27,945 records was recanonicalized by this independent method and matched to the new bank, including its cycle count. No source record was missing and no two source records were isomorphic. Equality of the two complete graph sets was checked in both directions, not just equality of their histograms.

The source record stream has SHA256 `1d0159168c03fcd011e90caab5f0fe1dd170362ba414d2f74730fe9d644b480a`, agreeing with the frozen source certificate. The CSV is reproducible and is not required as a separately published input bank.

## 4. Primitive implementation controls

`check_c3_primitives.cpp` supplies additional independent controls:

- 100,000 packed saturating-addition tests against scalar lane-by-lane arithmetic.
- Literal Hamilton-word permutation counts for every one of the 33,867 labelled tournaments through order six, checking the entire endpoint matrix against the capped subset DP.
- 31,668 directed-cut cap controls.
- Full unrestricted permutation canonicalization of every labelled tournament through order five, confirming both invariance and absence of false identifications in the refinement canonicalizer. Its unlabelled counts are 1,1,2,4,12.

Both the complete independent census and these primitive controls were compiled and run at `-O2` and at `-O0 -D_GLIBCXX_ASSERTIONS`. Their respective receipts agree byte for byte.

## 5. Why this finite bank is complete for the theorem

The strong-extension completeness argument is sound. A mixed-incidence vertex extends a strong induced tournament strongly, and inserting that named vertex in a Hamilton-cycle gap injects old cycles into new cycles. If there is no mixed vertex, the outside vertices split into those dominating S and those dominated by S. Strongness forces a dominated-to-dominating arc, giving the exact two-vertex ear S -> q -> p -> S. Its Hamilton cycles are precisely the H(S) old Hamilton paths, so their count is at least |S| C(S). Repeating the extension proves induced cycle monotonicity.

For a current strong S of order at least four inside a tournament of at most three cycles, a two-vertex ear would already create at least four cycles. Therefore the entire remaining extension can be made by mixed singleton steps. Every strong tournament of order at least four contains a strong four-vertex induced graph: from a directed triangle either use a mixed vertex, or use two triangle vertices and the two-vertex ear. The unique strong four-vertex class and its count one are correct. This justifies the seed, every generation step, and the absence of a missed two-vertex-ear branch in the census.

The all-order theorem also requires the minimal-bridge-free-counterexample bound twelve. Its new covering argument is valid, assuming the previously independently reviewed C1/C2 classifications and separating-vertex product stated as dependencies: a bridge-free C3 counterexample is 2-strong; deleting each vertex exposes a degree-one vertex in at least one orientation; vertices of outdegree two or indegree two then cover all vertices by their respective neighborhoods. Each degree-two side together with its neighbors has at most six vertices, including the empty-side case, so the total is at most twelve. The bridge-to-degree-one implication uses the exact endpoint-one lemma and correctly allows weak pieces and either orientation.

The companion structural review owns the detailed audit of those already reviewed C1/C2 and endpoint prerequisites. The present independent finite bank excludes every bounded counterexample supplied by that reduction. No stabilization premise, unrestricted path-spectrum conclusion, or all-order C3 class-count generating function is inferred from the numerical row.

## 6. Audit artifacts and reproduction

- `check_c3_complete.cpp`: SHA256 `f5d974967398108a01b98e298bbc1a266db41e6e1b94b9abd52331ce5ccb2946`.
- `c3-independent.txt`: SHA256 `6175f8f9e3a3daa477c255ba2aba7c1ce91668417fd90f4af9592b4df04355ad`.
- `check_c3_primitives.cpp`: SHA256 `546a0dab77866a7a79ab13c8815e046d5288397d6fb340495cc38aa534cfc839`.
- `c3-primitives.txt`: SHA256 `8a40fc632496088af29495519f5516f1031e903e76eba0c68ab8623e1a7e3592`.

First compile the named source C++ verifier and run it with a fresh audit CSV filename. Compile `check_c3_complete.cpp`, then run it with two arguments: a fresh independent CSV filename and the source CSV filename. Compile and run `check_c3_primitives.cpp` in the same directory as the independent checker, which it includes. No network, external graph library, precomputed graph bank or repository modification is needed.
