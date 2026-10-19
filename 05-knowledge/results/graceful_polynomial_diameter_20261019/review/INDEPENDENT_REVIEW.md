# Independent review: polynomial rank-staircase reconfiguration diameter

2026-10-19. **PASS.** The universal finite-matroid theorem and its graceful upper-tail corollary follow from the stated staircase hypothesis and the correctly cited strict-rank theorem. No mathematical correction is required.

Reviewed proof: [the proof](../POLYNOMIAL_STAIRCASE_DIAMETER.md) (reviewed input hash; publication-only changes are recorded in SOURCE_PROVENANCE.json), SHA256 `57322dd240c0107cfc93db519259a2591296eb873484c8f55de6b9cffe2f5414`.

## Exact conclusions approved

For s disjoint color classes satisfying rank(A_i) >= s-i+1, independent transversals have diameter at most s^2-s+1 under single-coordinate moves and double-coordinate moves that preserve the full selected-set closure. If the maximal tight suffix leaves a colors outside it and consecutive tight blocks have sizes b_j, the sharper bound is g(a)+sum c(b_j), with g(0)=0, g(n)=n^2-n+1 for n>0, c(1)=1, and c(b)=b^2-3b+3 for b>=2.

For graceful upper-tail completion fibers, the total free system is tight. The diameter is at most r^2-3r+3 for r>=3, at most one for r=1,2, and zero for r=0. This is a sufficient diameter bound, not an optimality assertion. The theorem neither proves one-color connectivity nor prescribes the abstract tree type.

## Primary source independently verified

D. Fon-Der-Flaass, *Local equivalence of transversals in matroids*, Electronic Journal of Combinatorics 3(1), R24 (1996), DOI https://doi.org/10.37236/1248 . Primary publication: https://www.combinatorics.org/ojs/index.php/eljc/article/view/v3i1r24 . Primary PDF: https://www.combinatorics.org/ojs/index.php/eljc/article/download/v3i1r24/pdf/ .

I independently retrieved and read the primary PDF's indexed full text on 2026-10-19. Theorem 1, printed page 2, assumes strict union rank for every nonempty subcollection and concludes single-coordinate connectivity. The paragraph following the proof of Theorem 2, printed page 3, explicitly attributes the n^2-n+1 bound to the proof of Theorem 1. It does not require the separate thickness hypothesis of Theorem 2. The paper defines transversals as independent tuples of actual matroid elements, so distinguishable parallel realizations are permitted. Screenshot retrieval failed with a cache miss; no saved PDF image is asserted. The primary text retrieval was successful and unambiguous.

## Proof audit

1. **Existence and tight chain.** Reverse greedy selection works because the already selected s-i elements cannot span A_i. If I is tight and i=min I, the chain |I|=rank(A_I)>=rank(A_i)>=s-i+1>=|I| forces I to be the entire suffix. This proves nestedness for all tight sets, not merely among suffixes already examined.
2. **Fixed contraction flat.** A state's |K| independent selected elements in a tight suffix K span cl(A_K), because that flat has rank |K|. Thus contraction is independent of the chosen state. Rank in the quotient is rank(A_(K union J))-|K|. Equality at |J| would make a strictly larger tight set. Consequently the outside prefix satisfies the exact strict-rank source hypothesis.
3. **Consecutive blocks.** Their quotient rank is the difference of consecutive tight ranks, namely b. Contraction lowers rank(A_i) by at most the inner rank, so the first class has quotient rank at least b and hence exactly b. Equality in the weak Hall inequality for a nonempty subset of the remaining block colors would produce an intermediate tight suffix, which does not exist. Loops created by contraction may be discarded without altering these statements or any state.
4. **Lift of each partial step.** A partial selection has quotient rank b-1, while the first color has rank b. A buffer element outside its closure therefore exists. Replacing the changed partial coordinate and possibly the buffer yields another basis of the same block flat. Together with the unchanged inner basis it spans exactly the fixed outer tight flat. Keeping all outer coordinates unchanged therefore preserves both full independence and the full selected-set closure. This validates actual width-two moves, rather than just abstract quotient moves.
5. **No missing last buffer step.** If the partial path is nonempty, its final vertex is the target partial selection. The target buffer belongs to an independent target block over the already completed inner flat, so it can be installed on that last step. This remains valid with the current outer coordinates because the tight outer flat does not change. If the partial path has length zero, only the buffer may differ, requiring at most one move. Thus c(b)=max(1,g(b-1)) is correct.
6. **Completed inner coordinates.** Processing inside out never changes a completed inner coordinate. The proof does not assume that an initially correct buffer must remain fixed; the known target-preservation hostile is respected.
7. **Sum bound.** For positive x,y, g(x+y)-g(x)-g(y)=2xy-1>=1, while c(b)<=g(b). Zero prefix is harmless because g(0)=0. This proves the general g(s) bound. For the all-tight corollary, c(x+y)-c(x)-c(y)=2xy-3 for x,y>=2, and c(x+1)-c(x)-c(1)=2x-3 for x>=2. Thus all nonsingleton parts and then singleton parts can be merged. The exceptional pair (1,1) is handled by the all-singleton case, where r<=c(r) for r>=3. The r=2 bound is obtained directly by one width-two move, not by the block-sum estimate, which can be two.
8. **Graphic specialization.** An acyclic fixed set F has m-r edges, so its contraction has r+1 components and free rank r. The difference-d graph on labels 0,...,m has d residue-chain components. Adding F can only decrease the count, giving quotient rank at least r+1-d. Difference 1 is connected, making the whole free system tight. Distinct parallel edge realizations remain distinct coordinates.

## Independent finite checks

`independent_controls.py` imports no producer code. It computes binary spans as explicit sets, constructs the complete allowed-move state graph, and checks every state pair in each bounded system. Normal and optimized runs agree byte-for-byte on **202,884 checks**, **1,844 systems**, and **81,490 unordered state pairs**:

- 1,568 exhaustive F_2^3 class systems with a fixed first-color basis; all 2- or 3-element second-color subsets and all 1- or 2-element third-color subsets of the seven nonzero vectors.
- 120 deterministic systems with loops, distinguishable parallel copies, and possible strict outside prefixes in ambient dimensions three through five.
- Seven prescribed multiple-tight-block systems.
- Every nonempty graceful upper-tail fiber for labels 0,...,m through m=6, yielding 149 fibers; every pair of states is checked.
- Dynamic programming over every possible composition type through total size 500 verifies both block-sum inequalities independently of their algebraic proof.

The finite checker verifies tight-set suffix structure, block ranks, strict quotient Hall conditions, connectivity, and both diameter bounds. The universal theorem rests on the proof and primary source, not extrapolation from these controls.

The producer's normal and optimized checker were also rerun independently. Both reproduce `controls.json` exactly: **29,727 checks**, **2,158 constructed moves**, 200 binary systems, 585 graceful upper-tail fibers through m=7, and the target-preservation hostile. The producer's final target-buffer insertion and closure checks are active under optimization.

## Algorithmic and novelty scope

The proof is a short-path existence theorem with a blockwise lifting construction. The supplied implementation uses BFS for finite strict-rank subproblems. No general polynomial-time oracle implementation has been certified. The new step is the tight-suffix staircase decomposition and its lifted bound; the strict-rank connectivity and n^2-n+1 bound are attributed to Fon-Der-Flaass. No global literature-priority claim is certified by this review.

No files inside `math/` were edited, and no Git write was performed.
