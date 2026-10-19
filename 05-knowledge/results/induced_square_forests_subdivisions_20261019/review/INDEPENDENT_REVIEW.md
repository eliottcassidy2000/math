# Independent review: induced square-sum forests and subdivisions

Date: 2026-10-19. Verdict: **PASS** for the mathematical statements and bounds reviewed below. No external priority claim. This review records the pre-publication source audit. Publication updated status labels, dependency links, and package-relative commands only; SOURCE_PROVENANCE.json records every change.

## Audited conclusions

1. **Finite forests at near-quadratic height: PROVED.** For every eta>0, every n-vertex forest admits an induced square-sum realization with largest label at most C_eta n^(2+eta). The star lower bound makes the worst-case exponent exactly 2. The statement does not prove O(n^2) with an absolute constant.
2. **Explicit quartic bound: PROVED.** The stated height 2(64(n+1))^4 follows from tau(d)<=9d^(1/4). The exact fourth-power prime-factor product is 127401984/25025 < 9^4; the listed maximizing exponents are correct.
3. **Countable forests and exact finite seed extension: PROVED.** The repaired parent-before-child enumeration by inserting finite ancestor paths works even with infinitely many vertices in a level. An arbitrary finite induced seed, including a triangle, is retained when each added vertex has at most one required previous neighbor. The uniform finite bound is O_eta(max(H,n^(2+eta))). An auxiliary vertex may increase the bookkeeping count from n to n+1 without changing this asymptotic claim.
4. **Induced two-vertex bridge: PROVED.** The large-parameter construction adds exactly A-w-v-B between distinct old labels, preserving any old A-B edge. With A=B it attaches a triangle. The square-gap checks cover every old vertex and the collision bounds hold. Consequently finite connected unicyclic graphs have polynomial-height realizations, and countable cacti are qualitatively induced-realizable.
5. **Generic simultaneous two-subdivision bound: PROVED.** The common-window argument gives O_eta(K^(4+eta)), K=n+2m. Its codegree differences are within 5R^4; the mapping r to |s| is strictly increasing in the chosen window; each collision removes at most one parameter. The countable extension requires the stated finite-stage fair vertex/edge schedule, rather than an already-infinite seed.
6. **Mod-seven strengthening: PROVED.** `../SUBDIVISION_MOD7_CONSTRUCTION.md` gives an absolute O(K^4) bound for the exact two-subdivision of every finite simple graph. It audits all six pair types among original, V and W vertices and proves injectivity. This bound uses no divisor estimate. The coordinating reviewer independently checked its inequalities.

## Main proof audit

At a tree insertion with k previous vertices, at most k occupied labels and (k-1)lambda common-neighbor labels are forbidden among the parent's neighbors. This overcounts harmlessly: the parent itself is not its own neighbor, and exclusion sets may overlap. Since k<=n-1, the displayed strict minimum-degree inequality suffices. The parent-before-child ordering ensures exactly one earlier required neighbor.

The forest reduction must introduce and later delete an auxiliary root. Connecting components directly would leave unwanted edges; the text explicitly rejects that incorrect shortcut.

For interval labels [M,2M], the arithmetic difference satisfies 1<=|v-b|<=M. Each positive difference-of-squares factor pair determines one candidate, so the tau bound is valid. Counting square roots in the closed interval loses at most one, and removing a loop loses at most one more. The lower root-interval length occurs at b=2M. The weaker sqrt(M)/4-2 bound is conservative.

For the asymptotic choice M=C n^(2+eta), the exponent comparison is 1+(2+eta)epsilon < 1+eta/2, exactly the stated epsilon<eta/[2(2+eta)]. Constants can absorb all small n. To formalize n^(2+o(1)), define the worst-case minimum height over the finite set of n-vertex forests; its logarithm divided by log n tends to 2 by the upper bounds and the star lower bound.

For a seed of height H, all parent labels are <=2M and positive when M>=H, so the same root-interval lower bound applies. Differences are now at most 2M. With M=C_eta max(H,n^(2+eta)), the degree dominates nD(2M) uniformly in H. A finite auxiliary root attached once to the seed may be used to produce new independent components; delete it after construction. Empty seeds are covered by the original forest result.

The fixed pair 1,2 has no common positive square-sum neighbor, since distinct positive squares cannot differ by 1. This is a correct obstruction to arbitrary exact two-neighbor attachments, and is not contradicted by adding a path with two internal vertices.

## Exact independent controls

Two standard-library scripts use explicit exception checks, never optimization-disabled assertions. Ordinary and `python -O` outputs agree byte for byte.

### Forest and arithmetic controls

`independent_exact_review.py`, receipt `exact_results.json`: **6,316,664 always-active checks**.

- tau(d)^4<=9^4 d for all 1<=d<=1,000,000
- Small-prime maxima through exponent 100, with monotonicity checks certifying the remaining tail
- Direct common-neighbor sets versus independent factor-pair reconstruction on [M,2M] for every M=1,...,80 and M in {100,127,128,255,256,511,512,1000}
- Every increasing rooted forest parent array through order eight: 46,233 arrays in total, including every unlabeled forest shape in that range
- Paths, stars, binary trees, edgeless forests and combs through order 256
- Sixty exact one-parent extensions of the fixed triangle seed 482,3362,6242
- A hostile to the older ordinary greedy path: labels at positions 0 and 4 are 1 and 15, adding the unwanted square-sum edge 1+15=16

### Bridge and modular controls

`verify_mod7_subdivision.py`, receipt `mod7_results.json`: **1,961,781 always-active checks**.

- Every labeled simple graph through order five: 1,099 graphs total
- Complete graphs of orders six through twenty-four
- One hundred deterministic randomized edge lists, orders and orientations on at most forty original vertices
- 2,920 arbitrary-seed bridge cases, including equal endpoints for attached triangles
- Every pair of output labels is compared with the exact intended adjacency relation; distinctness, residue classes and numerical height are checked separately

Combined total: **8,278,445** always-active checks. These finite controls validate implementations and boundary examples; they do not replace the proofs.

Reproduce from the staging parent directory:

    python3 -B review/independent_exact_review.py > /tmp/forest-review.json
    python3 -B -O review/independent_exact_review.py > /tmp/forest-review-optimized.json
    cmp /tmp/forest-review.json /tmp/forest-review-optimized.json
    python3 -B review/verify_mod7_subdivision.py > /tmp/mod7-review.json
    python3 -B -O review/verify_mod7_subdivision.py > /tmp/mod7-review-optimized.json
    cmp /tmp/mod7-review.json /tmp/mod7-review-optimized.json

## Inheritance and literature boundary

The nearest repository theorem is `math/05-knowledge/results/square_sum_elementary_sparse_path_cycle_bounds_20261019.md`. It explicitly permits extra square-sum edges in its O(n^2) tree construction and gives the inherited two-vertex bridge identity. The present induced-forest theorem and simultaneous modular induced subdivision are not stated there.

The older `math/05-knowledge/results/zenodo_square_20260925.md` audits a particular six-vertex Pythagorean path and identifies its possible extra edge. It is relevant evidence that the nonedge predicate is substantive, but does not contain general induced-forest universality or its quantitative bound.

A targeted local search across square-sum and induced-forest/tree terminology located no duplicate of the present claims. This is not an exhaustive priority investigation.

Targeted public searches are substantially polluted by unrelated graph-labeling definitions:

- Ajitha, Arumugam and Germina, *On square sum graphs*, DOI [10.1080/09728600.2009.12088866](https://www.tandfonline.com/doi/pdf/10.1080/09728600.2009.12088866), defines injective edge weights f(u)^2+f(v)^2. Its tree results concern that different predicate.
- Sonchhatra and Ghodasara, *Some new sum perfect square graphs*, [original paper](https://ijpam.eu/contents/2017-113-3/11/11.pdf), defines injective edge weights (f(u)+f(v))^2. It likewise does not assert that f(u)+f(v) itself is square exactly on edges.

These primary definitions were checked. No matching external induced-universality theorem was identified in this limited search; external novelty remains unclaimed.

## Quantitative next target

The O(n^2) forest question remains open within this packet. The promising obstruction to address is aggregate blocking, rather than the maximum divisor count. For an interval host, the exact identity

    sum_{y != b} codeg(b,y) = sum_{x in N(b)} (deg(x)-1)

shows that average codegree over all labels is O(1), since both degrees are O(sqrt(M)) and the interval has order M labels. However this average alone does not control the adaptively chosen prior labels or edge-conditioned next candidates. A valid strengthening would maintain a uniform O(n) bound for each active parent's total blocked neighborhood, or prove a distribution-preserving randomized embedding statement. Neither is supplied here. Merely substituting average codegree for maximum codegree in the greedy lemma is invalid.

The mod-seven subdivision result succeeds through an additional exact residue invariant plus simultaneous placement. It is a genuine way around the divisor loss in that different construction, not a proof of the O(n^2) forest target.

## Exposition notes, not mathematical blockers

The countable-subdivision paragraph can begin directly with its fair finite-stage enumeration. The preliminary phrase about choosing all original labels first is immediately rejected by the paragraph and is unnecessarily confusing. Proof-candidate status labels can be promoted by the coordinating owner after this review. No mathematical corrections remain outstanding in the reviewed statements.
