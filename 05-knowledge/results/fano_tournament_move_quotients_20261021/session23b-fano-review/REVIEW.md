# Independent audit: tetrahedral Fano orientation and triangle-move quotients

Verdict: PASS. No mathematical correction is required in the frozen packet. The finite claims were independently recomputed over their complete stated universes. All seven source hashes match the producer manifest before and after the audit; both producer scripts reproduce their receipts byte-for-byte in normal and optimized Python. This review changes no source packet or repository file.

## Evidence and independence

`independent_audit.py` does not import either producer script. Its global generator traverses all 2,097,152 labelled seven-vertex orientations in Gray order, updating the seven outdegrees after each single arc change. It selects exactly 2640 regular tournaments, independently of the producer's outdegree-constrained recursive generator. For every selected tournament it tests all 5040 possible vertex orders directly against six arc predicates, rather than using Hamilton-path dynamic programming.

The resulting H classes have sizes 240, 1680, 720 for H=189,171,175. Full relabeling orbits certify that these are precisely three isomorphism classes, not merely three invariant values. The script constructs every cyclic-triangle neighbor of every tournament and checks all 36,960 directed edge incidences. Each individual vertex has the claimed quotient row. The 128 fixed-Fano orientations are independently constructed from XOR-zero triples. A breadth-first search from each of their 128 vertices in the full graph verifies all 16,384 ordered cube-pair distances, not only their adjacent edges.

`receipt.json` and `optimized-receipt.json` are byte-identical independent-audit results. The script also replays the producer only after its own enumeration and checks. `small_circle_witness.py` and its normal/optimized identical receipts give an additional exact rational geometry certificate. All checks use explicit exceptions, so optimization does not remove them. No random sample, floating-point tolerance, numerical eigensolver, or external Python dependency is used.

## Tetrahedral incidence and the orientation seam

The reduced binary incidence columns of K4 are exactly labels 1 through 6. Its four triangular faces and the three opposite-edge pairs with the formal point 7 are exactly the seven XOR-zero triples. Pair coverage is exact. The full Fano automorphism group has 168 elements; the stabilizer of point 7 has 24, agreeing with the faithful S4 edge action. The identification is combinatorial and retains neither a sphere nor a Euclidean metric.

Using the explicit coherent oriented boundary faces (123),(032),(013),(021) gives fixed normal-coordinate bits m1=0,m2=1,m4=0,m7=1. The eight completions are exactly 66,70,82,86,98,102,114,118; all have H=171. For every completion, both members of an opposite-edge pair have identical directions to the other four edge-vertices and exactly two common outneighbors in the full tournament. Thus the stated structural obstruction to the doubly regular Paley tournament is valid, independently of its Hamilton count. Swapping members of the opposite pairs provides the claimed isomorphisms among the eight completions.

The code parity argument is also valid. Dotting the syndrome with 7 selects precisely normals 1,2,4,7. All Paley code words therefore have odd parity on those four face bits, while the coherent boundary has even parity. Flipping one face is necessary and sufficient to remove this parity obstruction. The three remaining normal columns 3,5,6 have rank two and one relation, so each compatible syndrome equation has two solutions. The verified repair pairs are:

- Face 1: 83 and 103.
- Face 2: 84 and 96.
- Face 4: 74 and 126.
- Face 7: 6 and 50.

The 16/112 Hamming-coset split and regular7 classification are inherited, not new claims: the local September 17 `arithmetic_braids2_20260917_fano_code.md` gives the code and H census, and September 22 `collatz_mod6_20260922_paley_fano_octonion_design.md`, S17–S18, gives the full regular census and contained Fano designs. The new checks agree with both.

## Cube, syndrome quotient, and exact-state distinction

The seven Fano lines partition the 21 pairs, so their independent cyclic orientations form a literal Q7. The independent full graph confirms that each cube vertex has exactly seven internal and seven external triangle-reversal neighbors, and no additional cube adjacencies. A non-line cyclic triangle changes one arc in each of three distinct line blocks, leaving that presentation family. If two presentations differ in k line coordinates, their tournaments differ in 3k arcs. A triangle reversal changes only three arcs, proving the lower distance bound k; reversing those k whole lines realizes it. This establishes the isometric claim structurally as well as computationally.

The syndrome map has sixteen states per fiber and a local neighbor bijection onto K8. The parity-enriched map has eight states per fiber and the exact crown graph K8,8 minus its diagonal matching. Neither quotient should be mistaken for the original state graph. Every base triangle has nonzero weight-three codeword holonomy; exactly seven support types occur. In particular, at length three there are 42 words returning to the Paley locus and no words returning to their original cube state.

The formulas for the locus-return generating function, H expectation, and geometric first-hitting time with mean seven are correct. They describe the specified observable or locus, not full-state mixing; the warning about cube bipartiteness is essential and correct.

## Full regular-tournament graph

With source-to-target rows and class order (189,171,175), the full quotient is exactly

    ((0,14,0),(2,9,3),(0,7,7)).

All 2640 vertices satisfy this, so the quotient is equitable rather than an average transition table. The degree is fourteen because every regular7 tournament has 35−7*C(3,2)=14 cyclic triples. Reversing one preserves all scores and distinct triples give distinct neighbors.

Independent trace/trace-square/determinant calculations give characteristic polynomial X³−16X²+14X+196, hence (X−14)(X²−2X−14). The stationary class masses are (1,7,3)/11 and stationary mean is 1911/11. The first-step equations reduce to 5a−3b=14 and −7a+7b=14, giving Paley first-hit means a=10,b=12. Actual full-state distributions from a fixed Paley vertex, not merely quotient powers, verify the H-weighted series through length twelve. The exact equitable quotient proves the all-length rational identity with numerator 189−630z−1470z² and denominator (1−14z)(1−2z−14z²). The numerator shares no denominator factor, so all three modes survive.

Operation scope is crucial: this graph reverses one currently directed triangle. It is not the original nonpath-complement operation, a single-arc operation, or a tournament-order/marked-path/renewal generating function. No operation-intertwining map to those carriers is supplied or implied. No larger-order equitability claim follows from n=7.

## Anti-symmetry and the two Fano sheets

The statement that every anti-action exchanges the two Paley Fano planes needs exactly two contained designs, not merely one decomposition. The audit explicitly enumerates all C(14,7)=3432 seven-block subsets of the fourteen cyclic triples and tests exact coverage of all 21 pairs. Exactly two pass; they are block-disjoint. Each of the 21 independently enumerated Paley anti-automorphisms exchanges these two. This also agrees with S17–S18 of the inherited September 22 note.

The general sign proof is sound. An anti-automorphism has at most one fixed point. An odd nontrivial cycle would yield an odd power fixing at least two vertices, impossible for an anti-action. A cycle of length divisible by four would yield an even power interchanging an opposite pair, impossible for an ordinary tournament automorphism. Thus all nontrivial cycles have length 2 modulo 4. At odd order 3 modulo 4 there is one fixed point and an odd number of nontrivial cycles, so the permutation is odd. Every Fano collineation is even, as seen either from transvection generators or the complete 168-element check. Therefore no anti-automorphism of any seven-vertex tournament preserves a single labelled Fano plane.

All stated ordinary/extended subset stabilizer counts and all seven marked-path rows also pass independently. The distinction between an ordered path and its underlying subset is valid. The Fano Property B observation in the separate palette packet is compatible, but its acyclic-color predicate is different from distinguishing color and proper graph color.

## Real-linear obstruction versus selected small circles

The determinant proof is correct under the standard faithful rank-three/projective representation meaning. Normalize binary-labeled points 1,2,4 to e1,e2,e3 and joining points 3,5 to (1,1,0),(1,0,1). The incidences put point 7 at (1,1,1), then point 6 at (0,1,1). The required last collinearity of 3,5,6 has determinant −2. This is impossible over R, but vanishes in characteristic two. Great-circle incidence has the same obstruction. Antipodal point identifications cannot evade it faithfully, since they would give identical incident great circles to two points with distinct Fano incidences.

The selected-small-circle caveat is genuinely different and can be made entirely explicit. For t=1,...,7, set

    P_t=(2t,2t²,t⁴+t²−1)/(t⁴+t²+1).

These are seven distinct rational points on the unit sphere. Exact Fraction calculations verify that all 35 triple planes are noncentral and contain no fourth labelled point. Select only the seven planes corresponding to XOR-zero triples. Their intersections with the sphere are seven small circles realizing precisely the designated point-circle incidence. Additional unlabelled circle intersections are not Fano points. This is an imposed incidence structure, not a canonical metric realization or linear representation.

The cited Baez section was opened and supports the signed Fano multiplication/background, not the new move graph. The cited Oxley PDF URL was reachable but its extracted text was unusably encoded; this audit therefore does not independently certify its exact Figure 8/Proposition 5.3 attribution. The elementary determinant and deletion proofs do not depend on that citation extraction.

## Publication recommendation

Accept the frozen packet with its inherited/finite-exact/proved labels intact. For maximum clarity, add a sentence explicitly citing September 22 S17–S18 (or this bounded exact-cover certificate) for “exactly two contained Fano planes,” and explicitly state that triangle reversal is not the original nonpath-complement operation. These are explanatory strengthening suggestions, not detected mathematical errors. Keep the rational small-circle witness separate from any claim about the user's given geometric feature points.
