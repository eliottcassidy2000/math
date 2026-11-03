# Independent audit of the Session145 bridges

Date: 2026-11-03 UTC. **Verdict: PASS for the final publication scope.** This is an independent mathematical review and exact computational reproduction, not a formal proof-assistant verification or a claim of historical novelty.

The original discovery snapshot was preserved. Two subsequent clarifications were reviewed: the knot paragraph now explicitly distinguishes arc-only tangles from tangles containing closed components; the fence model now explicitly excludes collinear overlaps. Both changes are correct. The tangle correction has new exact controls. See [CORRECTION_LINEAGE.md](CORRECTION_LINEAGE.md); no identified scope issue remains unresolved in the final publication revision.

## 1. Magic labels and the fixed-placement family

The following proofs are accepted:

- Every real 3-by-3 magic square has center equal to the mean of its nine entries. For the nine stated decimal labels this is 555, absent from the input. Consequently no rearrangement works, and no nonconstant affine relabeling repairs the obstruction.
- In the prescribed rank layout, the exact decomposition is `528+111U+(567/2)V+(81/2)V²`. U and V have all eight line sums zero. V² contributes two on seven lines and zero on the anti-diagonal, giving the exact 81 defect.
- The space of real magic arrays is three-dimensional, with basis 1,U,V. Thus the complete family of value transformations on these nine fixed inputs that works in this fixed arrangement is `c+pU+qV`. This is not a classification for arbitrary arrangements or all real inputs.
- Strict increase is exactly p>0 and q>2p. Positivity is exactly c>|p|+|q|. Injectivity is exactly p,q nonzero with p/q outside {−2,−1,−1/2,1/2,1,2}. Constants and coordinate projections are repeated-value magic arrays, not bijections.
- The orthogonal projection of the raw grid to the magic space is `555+111U+(567/2)V`, uniquely minimizing squared error, with error 6561/2.
- The increasing rank map and trace-address cubic map are distinct bijections. They preserve the stated three-fiber partitions, but transported digit rotation fails to preserve the numerical magic predicate.
- The degree-seven interpolants and the nine-node annihilator give the complete fixed-trace polynomial classification. The polynomial identity cancels exactly to the exceptional cubic. A second independent nonconstant mode appears first at degree seven; through degree six the family is constant plus that cubic.

Independent exact linear algebra confirms the six independent magic constraints, projection formula, and degree filtration. A full **9!=362,880** permutation census, without assuming the center theorem, finds exactly eight target arrangements and identifies them with the D4 orbit. **2,025 rational parameter pairs** independently check injectivity, order, and positivity boundaries. The checker imports no discovery implementation.

All discovery checks also reproduce normally and under Python -O, including nine result-mutation rejections, eleven mathematical hostile examples, and 625 integer parameter pairs. Result mutations are integrity controls, not a claim of testing every possible algorithm mutation.

The three archived inherited sources were compared directly with their blobs at repository revision `7b76044ca13e6684d1746923285b83dbf83828d2`. The elliptic-curve realization and point counts remain inherited from the separate Frobenius packet; they were not independently re-established here. No acoustics, health, curve-level action, or canonical identification with a separate 126-dimensional carrier is inferred.

## 2. Fence and product geometry

The primary [fence page](https://erich-friedman.github.io/packing/fence/) and [product page](https://erich-friedman.github.io/packing/product/) were independently opened. Their objectives and intersection rules differ. The fence model permits supported endpoints in another segment's interior and imposes a unit area cap. The endpoint-only product model does not impose that area cap. Both pages describe best-known configurations, not universal optimality certificates.

### Fence pseudoforests

The n/4 bound and equality case are accepted under the final explicit no-overlap model. Splitting all junctions gives minimum degree at least two. In each pseudoforest component, E≤V and 2E≥2V force every degree to be two. No through-fence junction survives, so every cycle edge is an original unit fence. Triangles satisfy the equilateral area bound; cycles with at least four edges use the field-area cap. Nested cycle regions are correctly subtracted from containing fields.

Equality requires four-edge cycles with area one. A unit quadrilateral reaches that area only as a square, and positive nested components reduce the containing field's area. Precisely disjoint, non-nested unit squares attain n/4. The seven-fence adjacent-square example has area two but cycle rank two, so it is a valid hostile to an unrestricted extension.

### Endpoint-only product cacti

The exact balanced-cycle formula is accepted for all cactus graphs and all pseudoforests, with at least one bounded face. Cycle edges are disjoint and remain unit edges. The classical polygonal isoperimetric inequality bounds each cycle's disk; nested fields have no greater area. Strict increase of `a(x)=x cot(pi/x)/4` exhausts the edge budget, and strict concavity of log a(x) balances cycle lengths at each fixed cycle count. Disjoint regular polygons attain every balanced candidate.

The written derivative identities were independently checked symbolically. The cited [Bogosel author-hosted paper](https://beniamin-bogosel.github.io/pdfs/DiscreteIsopShorter_Rev.pdf) was inspected. The no-bounded-face convention is essential: otherwise an empty product would change small-n cases.

### Exact twelve-edge witness

Independent radical arithmetic confirms twelve original unit segments, eleven vertices, endpoint-only intersections, no crossings or overlaps, and strict convexity of each six-corner field. A fan-triangulation area calculation, separate from the discovery half-edge traversal, gives product `(89+26√13)/16`, strictly greater than the cactus maximum `6+3√3`. The common length-two side consists of two original unit edges. The graph is theta(5,5,2), with cycle rank two, not a cactus.

Each field has area greater than one, so this is not a fence witness. The five-fence/two-triangle T-junction fixture independently has cycle edge lengths 1,1,1/3 and seven split edges from five original fences. It refutes a proposed unit-cycle argument for general fence cacti, not the unproved n/4 inequality itself. Only the endpoint-only cactus fence corollary is established.

The discovery geometry checker reproduces **4,014 gates**, **1,099 graphs**, and **1,885 partitions** under normal and optimized Python. Its partition comparisons use rational intervals, not arbitrary floating-point tolerances. Its face traversal is used for its connected/non-nested fixtures, not asserted to implement arbitrary nested face geometry.

Five inherited repository sources were checked directly at the named revision. The quoted THM-4509 constants and EMPIRICAL/OPEN distinctions match those sources; its older computer-assisted convex-field bound was not re-proved in this audit.

## 3. Bounded theta optimality

**Accepted only for the written class:** legal endpoint-only theta(5,5,2), the two five-edge paths forming a simple outer polygon, the two-edge path separating its two bounded fields, and the open straight terminal chord lying inside that outer polygon.

The chord closes each outer five-edge chain to a simple hexagon with side lengths (1,1,1,1,1,ell). Their areas sum to the outer area, independently of the initial separator's bending. Fixed-side maximality gives total area≤2g(ell); AM-GM gives product≤g(ell)². The two-edge separator implies ell≤2. The half-angle branch, monotonicity, and polynomial derivative identity are correct, giving the exact upper bound `(89+26√13)/16`.

Equality requires ell=2, forcing a straight separator with its intermediate vertex at the midpoint. Fixed-side cyclic uniqueness gives the two reflected maximal caps. Independent exact tests confirm the derivative identity coefficient by coefficient, the valid branch, the ell=2 algebraic solution, the exact area, and cocircularity of all six witness corners, with circle center (0,D). All **15 independent gates** and **177 discovery gates** agree under normal/-O execution.

The [published Kouřimská–Skuppin–Springborn chapter](https://link.springer.com/chapter/10.1007/978-3-662-50447-5_5) directly supplies fixed-side maximality as Theorem 1.5 and cyclic existence/uniqueness as Theorem 1.1. Its [arXiv version](https://arxiv.org/pdf/1506.08069) labels these Theorem B and Theorem 1. These are cited classical inputs.

A chord leaving a nonconvex polygon, a different separating theta path, and unrestricted twelve-edge optimization are outside the theorem.

## 4. Sparse Ramsey graphs, matching loops, knots, and 1729

### Restricted graph bounds

The four sparse non-green graph classes have chromatic bounds 2,3,4,5, respectively. Their induced-subgraph edge inequalities justify the recursive degree deletion; the argument is not only an average-degree assertion about the original graph. The stated K2/K3/K4/K5 blocks attain the corresponding thresholds while respecting the red/blue triangle restrictions. These are exact all-order restricted maxima, not unrestricted Ramsey bounds.

For one surviving monomial in a product of r attachment observations, its support lies in a union of r linear forests. Every induced subgraph has fewer than r|S| edges, hence is (2r−1)-degenerate and 2r-colorable. If that support is the actual non-green graph, the claimed `n≤2r(k−1)` consequence follows. A polynomial's aggregate support, a signed coefficient, or the full source tournament cannot be substituted for that one monomial's support. The necessary r≥6 conclusion is not sufficient construction or a carrier-dimension lower bound.

### Matching quotient and knot boundary

The splice bijection proves `F_r(z)=∏(z+2j)` for all r. The doubled edge in P∪Q is retained correctly. Choosing Q=P separates different named matchings, proving exactly `(2r−1)!!` contextual states when every matching closure is permitted. This is not a linear dimension or a theorem for planar-only or direction-compatible closures.

The final tangle paragraph correctly restricts the matching-only component count to arc-only tangles. With m pre-existing closed components, it retains m and uses m+c(P,Q). The independent checker adds disjoint circles and verifies this correction by binary incidence ranks. No knot type, embedding, over/under information, twist, or framing is recovered from endpoint matching alone.

The Fox-3 braid hostile is accepted. Full crossing/closure systems independently give three colorings for the one-crossing closure and nine for the three-crossing closure. The algebraic Reidemeister identities and braid relation are checked, while their topological interpretation remains classical input. Fox coloring is not asserted to be complete.

### Paley-11

The forward/reverse endpoint split 628/1101 is correct. The square-affine group acts freely on Hamiltonian words and transitively on arcs. Independent enumeration of the fixed-first-arc chart yields **1,729** affine orbits, whose full expansion gives **95,095** words. Every one of the 55 directed arcs occurs in **17,290** words. Direct rotation canonicalization yields **5,505** directed cycles from 60,555 closing words. Thus 1729 counts the orbit/endpoint-channel observer, not a specified-arc incidence or knot.

### Reproduction and external bounds

The discovery checker has **106,874** gates. Its separate independent checker originally had 112,856; the repaired publication version has **112,989**, adding 76 closed-component equality controls and 57 extra-circle hostiles. Both revised normal/-O outputs were rerun and reproduce the final receipts exactly.

The present independent checker additionally covers all matching closures through fourteen ports, **2,208,196** ordered pairs of six-vertex linear forests, full Fox diagram systems, direct Paley orbits, and **496** matching-plus-closed-component cases. These finite checks supplement, rather than replace, the all-order written arguments.

The [Radziszowski revision-18 survey](https://www.cs.rit.edu/~spr/ElJC/ejcram18.pdf) supports the lower bounds 45 and 61 and upper bound 57. The [current Dauth author repository](https://github.com/Lachy-Dauth/ramsey-neighbourhood-counting) reports the computer-assisted/preprint upper bound 86, superseding 91. Hence the source-audited intervals are **45≤R(3,3,5)≤57** and **61≤R(3,3,6)≤86**. The external 86 certificate was not rerun here; direct PDF fetches failed but the author's repository table was accessible. No claim about absence of newer unpublished work is made.

## Reproduction and provenance

See [README.md](README.md) for self-contained package-relative commands. `SOURCE_EVIDENCE.json` retains exact original/final target hashes and direct revision-blob evidence. `AUDIT_RECEIPTS.json` records the discovery replays and narrow amendment checks. The three independent checkers import no discovery implementations; their normal and optimized outputs match. `MANIFEST.json` and `verify_package.py` define and verify the exact artifact whitelist.
