# Independent audit of Session33 geometric and tournament addenda

Verdict: PASS, with no mathematical corrections required. This report and its new checker/receipts are separate from the prior frozen audit package, which was left unchanged.

## Exact artifacts and reproduction

Run `python3 observer-audit/check_addendum_independent.py` and `python3 -O observer-audit/check_addendum_independent.py`. Outputs `addendum_independent.json` and `addendum_independent.optimized.json` are byte-identical. Producer replays `plane_producer_replay.json` and `addendum_producer_replay.json` match their saved originals exactly.

The checker imports neither producer. Its geometry engine uses all 3^7 signed difference coefficient vectors in integer-scaled coordinates, rather than the producer's 8,128 rational point pairs. The common denominator is408850. It proves all2,187 signed sums are zero exactly for the zero coefficient vector, and unit exactly for support size1; there are14 signed unit differences. Thus all128 subset positions are distinct, and every and only cube adjacency is unit distance. This is an exact integer certificate, not a tolerance test or probabilistic search.

The K8 Gram matrix is independently eliminated over rational numbers, with seven nonzero pivots and determinant1/16. Tournament profile checks use subset Hamiltonian permutations plus least-vertex unordered set partitions, not the producer permutation-and-cut formula. Triangle censuses use the score complement identity, rather than the producer's triple test. Context word counts use a finite automaton rather than star-slot enumeration; permanents use Ryser inclusion-exclusion rather than matching permutations.

## Geometric proof audit

The inductive embedding proof is sound for every finite dimension d. Distinct old subset sums imply that an old signed sum w with coefficients in{-1,0,1} can vanish only if its coefficient vector is zero: the positive and negative supports would otherwise be distinct subsets with the same sum. After adding a unit vector v, every genuinely new difference has the form w+v or w-v.

- Collisions with w nonzero restrict v to at most two circle points in total. If w=0, v cannot be zero, so no collision arises.
- For an unwanted unit distance with w nonzero, expanding the norm gives w·v=−|w|²/2 or +|w|²/2. Each line meets the unit circle in at most two points. No such equation is identically true on the circle because w is nonzero.
- If w=0, the two new subset vertices differ only in the new coordinate, so this is exactly a desired cube edge. Old differences remain unchanged.

There are finitely many old signed sums, hence finitely many excluded directions. The circle has infinitely many rational points (indeed they are dense through the rational parameterization), so a rational next direction exists. The base case d=1 is immediate. Finite induction proves the stated theorem; no infinite simultaneous embedding or uniform coordinate-size bound is asserted or needed.

For d=2^r−1, the nonzero syndrome columns give a locally bijective graph homomorphism to K_(2^r): every source vertex has exactly one neighbour mapping to each other quotient vertex. The source remains bipartite and has edges, so its chromatic number is exactly2. Quotient chromatic numbers are unbounded with r. For K8, seven vectors obtained by subtracting one point from the other seven would have Gram matrix (I+J)/2, of rank7; hence eight mutually unit-separated points need dimension at least7 and cannot lie in R².

The coloring direction is correct: a graph homomorphism G->Q pulls Q-colorings back to G, proving chi(G)<=chi(Q), not the opposite. A homomorphic section Q->G would give the opposite inequality, but K8 cannot map into a bipartite graph because its triangles cannot map to closed odd walks in a loopless bipartite target. Thus locally bijective quotients do not provide lower-bound-preserving maps for the Hadwiger–Nelson problem. The assertion concerns geometric unit-distance representations, not topological planarity. Nothing here improves a lower or upper bound for the chromatic number of the full plane.

The novelty framing is appropriately limited: graph-homomorphism monotonicity, graph covers, Hamming syndromes and finite generic-direction avoidance are classical elementary mechanisms. The report establishes the joint exact witness and its negative-transfer interpretation, not priority for a general theorem.

## Tournament addendum audit

A1 is a disjoint exhaustive classification: an (n−2)-path cover has one 3-vertex component or two 2-vertex components. A transitive triple gives1 Hamiltonian path and a cyclic triple3; each 4-set gives3 unordered pairings, each oriented uniquely. Therefore pc(n−2)=C(n,3)+2c3+3C(n,4). The cases pc(n)=1 and pc(n−1)=C(n,2) have the correct unordered normalization.

A2 is correctly extracted from E(u)=sum_d F(d)(u−1/2)^(d−1): the u^(n−3) coefficient receives exactly the top three F entries, giving the displayed expression and lattice increment2(n−2)!. The independent checker verifies A1 and A2 on all1,096 labelled tournaments of orders3..5; its generic binomial-shift extraction agrees with the producer's simplified three-term expression.

At order7, E differences are even of degree at most4 after cancellation of the common leading coefficient. Equality F1,F2 gives a double root at u=1/2 and hence by parity at u=−1/2. The result is a scalar multiple of (u²−1/4)². A2 identifies that scalar as240Δc3. Expanding at u=z+1/2 yields240Δc3·z²(1+z)², so the profile differences are ΔF1=ΔF2=0, ΔF3=240Δc3, ΔF4=480Δc3, ΔF5=240Δc3, with higher differences zero. In particular ΔR2=960Δc3. These statements are necessary algebraic conditions only. They neither prove existence nor rule out a hypothetical deeper twin, and the note correctly says so.

For the Hall hostile, the two matrices have the same labelled row margins(1,1,2) and column margins(2,1,1). Their permanents are0 and1, independently certified by Ryser. B_bad fails Hall on its first two rows; B_good has exactly the diagonal matching. With the specified transitive shores, both quotient score vectors are(3,3,2,2,4,4,3). Independent complete word-coefficient vectors, indexed c=0,1,...,7 and normalized to0 at c=0, are:

- bad: (0,145,150,50,0,0,0,0)
- good: (0,141,152,60,6,0,0,0)

Thus top coefficients0 versus6 agree with r!per(B), and same labelled scores do not even determine its support. Equal triangle counts follow from the score identity; this fact does not recover the missing higher-order interface information. Hall's condition is exactly the top coefficient's positivity criterion, not a claim that Hall subset coverage determines all lower coefficients.

The r<=2 margin-minimality proof is valid: for2×2 matrices, nonzero rows and columns suffice for a perfect matching, and permanent2 occurs only for the all-one matrix. Labelled margins distinguish all those cases. The independent checker exhausts all18 matrices at r=1,2. Therefore r=3 is minimal for THIS row/column-margin compression, not proven minimal for arbitrary quotient-score obstructions.

No source note, original audit, repository file, or Git state was modified. No further blocker remains on these addenda.

### Concrete prior-canon provenance found after proof review

`01-canon/theorems/THM-432-unit-distance-product-crossover-hamming-tie.md`, section(A), explicitly states the generic relative-angle Minkowski-sum realization of G Cartesian-product H as an induced plane unit-distance graph, with distinct points and no accidental diagonal unit distances, and labels it classical. Iterating K2 gives every Qd. Hence the all-cube embedding proposition is an inherited special case of this already-recorded mechanism. The current finite-direction proof is still valid and self-contained; the rational-coordinate strengthening is explicitly supplied by its own argument. Parent was advised to cite this exact filename because THM432 has an unrelated namespace collision. This provenance update does not change the PASS verdict or the scope of the joint Hamming-cover/negative-color-transfer witness. No claim of new general embedding theory is justified.
