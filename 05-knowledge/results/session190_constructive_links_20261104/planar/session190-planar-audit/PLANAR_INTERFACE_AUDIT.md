# Independent audit of short planar Hamiltonian interfaces

## Verdict

The constructive results pass: exact gluing, unique completion through five ports, realizable unit profiles, the six-port parity kernel, deletion-criticality with the stated interior-size restriction, the square-sum conclusions at their stated evidence level, and the complete-stellation hostile.

All three requested scope clarifications are incorporated in the frozen proof: universal-response/profile equivalence in Theorem 3 explicitly has k<=5; the universal-vertex path/cycle bijection explicitly has |G|>=2; and the linear-dimension statement is a lower bound on a factorization of the formal response operator or a linear/bilinear gluing representation. It is not a lower bound on arbitrary decoding of a real linear encoding restricted to integer profiles. None changes the constructive probes or finite conclusions.

## 1. Exact interface and the five-port theorem

The gluing argument is correct. With a private vertex on each side, neither restriction of a spanning cycle can already contain a cycle. Every positive-edge component is a path with boundary endpoints, and terminal degrees add to two. Degree-two boundary vertices have no exterior edge in the selected cycle, so their detailed positions along paths need not be retained. Contracting the paths leaves the union of the two endpoint matchings. Treating duplicate matching edges as distinct is correct: on two endpoints they form one alternating two-edge cycle, corresponding to an ordinary cycle after restoring the two original paths.

The count S_k=C_(k+1) is the two-colored Motzkin identity. Removing the 2^k patterns with no degree-one terminal gives L_k. Such a no-endpoint pattern cannot contain a spanning acyclic fragment with a nonempty interior.

For k<=5 the endpoint set has size two or four. The two-point matching is compatible with itself; among the two noncrossing four-point matchings, only the other matching completes one cycle. Complementing the degrees and making this matching change is an involution, so K is a symmetric permutation matrix. The dimensions 1,6,26,100 for k=2,3,4,5 are exact.

This does not identify coloring extension data with Hamiltonian data. The four-port hostile is valid: two length-two paths on each side give either two disjoint four-cycles or one eight-cycle, depending only on the pairing.

## 2. Physical unit probes

The collar construction is valid for every live state. Assigning each degree-two port to the preceding degree-one port partitions the degree-two ports into disjoint cyclic intervals. The collar chains occupy disjoint intervals; their tips remain in the endpoints' cyclic order. A noncrossing matching of these tips can therefore be drawn inside the collar without crossing the chains. Degree-zero ports can be bypassed in the collar and remain isolated.

After subdividing every raw edge, each private vertex has graph degree two. A feasible fragment must select both of its incident edges. Consequently every edge of the probe is forced, and the whole graph is a forest with exactly the desired degree vector and endpoint pairing. This proves a unit *count* profile, not merely existence of a fragment. There is at least one subdividing vertex because the state has at least one matched pair.

There are z+m raw edges when z ports have degree two and 2m ports have degree one. Since z<=k-2m and m>=1, there are at most k-1 private vertices and at most 2k-1 vertices total. No boundary-boundary edge remains, avoiding collisions with an arbitrary patch's boundary edges.

Through five ports these probes recover each profile coefficient individually. This proves necessity of equal profiles for equality of all exterior Hamiltonian-count responses. It concerns arbitrary simple disk exteriors, not only triangulated exteriors.

### Precise dimension scope

The rank statement is correct over every field. The safe lower-bound statement is: any factorization of the formal response operator through a d-dimensional vector space has d>=L_k; hence any linear/bilinear gluing representation has at least L_k coordinates. Rank gives this directly.

A bare assertion about any linear map on realizable integer profiles plus unrestricted decoding is stronger and not established by rank. For example, real coefficients linearly independent over the rationals encode a finite integer vector injectively in one real linear coordinate; a non-linear decoder could recover it. Likewise, existence of the unit-profile graphs alone does not establish closure of graph-realizable profiles under all formal linear combinations. The frozen producer proof adopts the operator-factorization wording.

## 3. Six-port block

The five displayed noncrossing perfect matchings are complete. Direct independent connectivity calculation gives the displayed K2 disjoint-union K3 matrix. Its determinant is (-1)*2=-2. The K2 block has rank two in every characteristic; the K3 adjacency block has rational rank three and characteristic-two rank two. Therefore the full six-port ranks are 365 over the rationals and 364 over F2.

The remaining 360 states form a permutation block, consisting of 180 transpositions. In fact the full determinant is -2, strengthening the stated absolute-value assertion. There are 362 rows of sum one and three of sum two.

The sum of the three K3 basis vectors is exactly the nonzero F2-kernel vector. The prose correctly treats this as a formal linear combination, not an unproved single-graph realization, and correctly distinguishes loss of injectivity of the full modular response operator from the simpler fact that an individual nonzero cycle count may be even.

## 4. Deletion and hypotraceability scope

Two private vertices on each side are sufficient for every one-vertex deletion to leave both sides nonempty. Boundary deletion must occur on both sides, with a shorter ordered boundary, as stated. Under this condition the same gluing theorem tests non-Hamiltonicity and every deletion card. The explicit warning about a side becoming empty is necessary and correct: an already completed cycle then requires a new state or a separate direct test.

For |G|>=2, adding a universal vertex z bijects unoriented Hamiltonian paths of G with unoriented Hamiltonian cycles of G+z. For testing hypotraceability, only deletions of original vertices are relevant; the z-deletion is not an extra condition. This construction may destroy planarity, so the augmented boundary requires unrestricted matchings. At six boundary vertices the ambient count is 64+240+180+15=499 and the live count is 435. No planar 100-state conclusion is being claimed for square-sum graphs.

The raw singleton exception is why the |G|>=2 proviso was requested: a singleton has a trivial spanning path, whereas its cone K2 has no simple cycle. Actual separator applications in this packet have more vertices and are unaffected.

## 5. Square-sum conclusions

The minimum-degree arguments for hypohamiltonian graphs and for hypotraceable graphs of order at least three are correct. In the degree-two hypotraceable case, deleting one neighbor leaves the degree-two vertex as a forced path endpoint; restoring the neighbor by appending it creates a spanning path in the original graph, a contradiction.

The four low-degree witness bands through order 31 are valid. For vertex 16, only 25 and 36 can occur as square sums within the stated interval, giving at most the neighbors 9 and 20. The smaller bands are immediate by the same calculation.

I checked the primary author-attributed statement on [OEIS A090461](https://oeis.org/A090461): Gerbicz's January 21, 2018 entry states path existence for every n>=25 and cycle existence for every n>=32. The packet explicitly cites rather than independently verifies that infinite construction. Conditional on this cited result, the all-n>=3 exclusion of both hypo-properties is valid. The Q2 exception is accurately retained under the singleton-path/raw-definition convention.

The six finite nontraceable graphs and all 124 deletion cards were independently recomputed here with an all-start visited-subset endpoint dynamic program, not the producer's path backtracking or its chain-kernel method. The exact restoring sets agree:

- Q18,Q19,Q20,Q21: {18}.
- Q22: {2,18}.
- Q24: {11,18,22,24}.

Thus the finite negative decisions have three algorithmically different checks, while the infinite claim remains citation-dependent.

## 6. Complete face-stellation hostile

For an n-vertex simple triangulation, Euler gives 2n-4 faces. Stellation adds that many pairwise nonadjacent degree-three vertices, yielding N=3n-4. Every new vertex is a legitimate reducible singleton: its three mutually adjacent neighbors use three distinct colors.

A cycle can use at most n independent new vertices, and a path at most n+1. This gives failure of Hamiltonicity for n>=5 and traceability for n>=6. After deleting a new vertex the independent side has size 2n-5; comparison with n and n+1 gives the universal deletion-card thresholds n>=6 and n>=7 respectively. Old-vertex deletion only strengthens the obstruction. The exact abundance (2N-4)/3 is correct.

The explicit 11-vertex triangular-bipyramid stellation is valid. Independent subset-DP cycle counts give zero cycles in the original graph and after each old-vertex deletion, but six cycles after each of the six new-vertex deletions. The displayed deletion-10 cycle uses every required vertex and edge. The original graph is traceable, as expected at n=5. This is a concrete and useful separation between coloring reducibility and Hamiltonian preservation, not a claim of a new graph family.

## 7. Reproduction

`independent_planar_audit.py` imports no producer code and checks:

- all 498 live states through k=6 using a Catalan-interval matching generator;
- all compatibility entries with a disjoint-set connectivity test;
- all 498 unit probes by independently constructed paths and exhaustive edge subsets;
- the determinant by the Leibniz expansion and F2 ranks by a separate elimination routine;
- all 189 stored physical graph cycle counts by a memoized visited-subset recurrence;
- all six square-sum base graphs and 124 deletion cards by all-start subset dynamic programming;
- all eleven stellation deletion-cycle counts and traceability.

Output: `INDEPENDENT_RESULTS.json`. Exact source hashes are recorded there. Original normal and optimized verifier replays also match saved outputs byte-for-byte, as do the original kernel and stellation replays. No installations or forbidden files were used.
