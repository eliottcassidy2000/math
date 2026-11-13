# Independent review of the pairing/wiring boundary

**Verdict: PASS for the corrected note and the stated orders 3,4,5.** The first-failure assertion is restricted to the nondegenerate n≥3 complement-line tower.

Initially reviewed snapshots (the first note snapshot contained the code-convention error below):

- `PAIRING_WIRING.md`: SHA-256 `9f66cdf708cc3daab0bb0bf4f1674c52cc2a26459ed3028dca8d7b60f039a9f2`.
- `audit/check_small.py`: SHA-256 `ea71adf577950bbda301271355f1893418d3017c77f4297bb330459d601d3406`.
- Different, previously frozen generator `metagraph/replay.py`: SHA-256 `aa656c169826ff90e26de26b6bc315a03994d70ade72810629ac6cc912a8c821`.

## Correction lineage: two canonical label gauges were confused

The root and figure agent independently identified a genuine prose error after the first review. The initial note claimed **unrestricted minimum over relabelings**, but its table and proof used the producer's **score-sorted canonical codes**. The initial review, SHA-256 `9be24baed5ddee4afa7300dffe7b9d976cfa78850f83a27d139032a83719a9e0`, failed to reconcile these gauges and accepted that inaccurate description. This was a review failure, not a new mathematical counterexample.

The producer canon sorts vertices into ascending outdegree buckets and minimizes only inside those buckets. The independent root checker instead enumerates full relabeling orbits and chooses the unrestricted minimum. Both are faithful class identifiers, but their numeric labels cannot be mixed. Recomputing the minimum over all120 permutations for every order-five representative gives:

| producer score-sorted code | root unrestricted code |
|---:|---:|
|0|0|
|2|2|
|4|4|
|8|8|
|32|5|
|40|40|
|68|41|
|76|76|
|97|12|
|258|10|

The additional ordinary chiral partners map256→11 and64→9. Hence the root's ordinary converse pairs are {2,11} with H3 and {4,9} with H5. Under the corrected root labels TT0 has sole neighbor10, while the H3 self-converse node5 has sole neighbor41. The correction relabels the table and adjacency statements; it changes no graph, invariant, pairing count, or group conclusion. This review's mathematical proof below now uses the root's unrestricted labels.

Final repaired note SHA-256: `aa26589421cadaab94f587fb6ba3de541bd815868249d5fa12e3601ae580e7a7`. A fresh post-repair pass compared all10 table rows with the unrestricted-orbit checker, recomputed the full12-class ordinary code transport, and matched all100 merged adjacency entries after transport: **115 additional always-active checks passed**. The ordinary-pair prose and explicit nondegenerate n≥3 qualifier were checked as well. No generator or checker changed.

## Proof and typing review

The weighted-vertex-rigidity proof is correct. At order5 the only repeated dart degrees are1 at unrestricted codes0,5 and9 at codes8,10. The degree1 nodes have unique neighbors of degrees9 and13, respectively, so each is fixed. Adjacency to the fixed transitive node then fixes both degree9 nodes. Other degrees are unique. This proves triviality of the **vertex** automorphism group of A, without any H assumption. It leaves parallel-edge or dart permutations fixing vertices intact.

Consequently every diagram-compatible ordinary-class permutation acts inside each converse orbit. The eight singleton orbits are fixed; each of the two size-two orbits can be independently flipped. This gives exactly C2×C2, even without imposing H preservation or involutivity in advance. All four permutations are involutions and preserve H because both classes in each nontrivial orbit are converse. The fixed-count histogram is {8:1,10:2,12:1}, and the unique eight-fixed element is global converse. This classification uses the note's explicit compatibility definition (R-equivariance and the induced weighted node action). It does not claim arbitrary independent pair flips have a preferred lift to every marked-path structure.

Global converse does have the claimed physical lift. For fixed path order, reverse the vertex order and all arcs. The resulting tiling is reflection of the free-pair bit array. Reflection preserves the free-pair set, commutes with its all-bit complement, and preserves whether a tiling is reflection-fixed. Thus it induces ordinary converse, respects the native dart-pairing involution, and preserves blue/black. The graph-guard minimum therefore has a literal presentation lift.

The unguarded H-preserving/converse-equivariant optimum4 at n5 pairs both H9 self-converse classes, both H15 self-converse classes, and both chiral pairs. The first exchange breaks adjacency to TT5, and the second exchanges degrees5 and3. At n4 degrees1,2,5 already force vertex rigidity and leave two fixed self-converse classes, so no free compatible pairing exists. At n3 the single edge can swap its two self-converse endpoints, changing H1↔H3; hence free diagram compatibility and optimal H-preserving compatibility have different first obstructions.

## Independent checks of the new checker block

The new enumeration is complete: H-fibre permutations exhaust every H-preserving permutation; the explicit square and converse-commutator tests retain exactly the requested involutions. Converse commutation makes the induced merged permutation well-defined. Testing A retains literal multiplicities and doubled loops. The additional B test suffices for both colors because the black matrix is A−B. Separately, degree-group permutations exhaust all possible weighted graph vertex automorphisms, without an H restriction.

A focused execution imported the new checker and ran `census(n)` only for n=3,4,5: **3,597 internal always-active gates** passed. An additional **308 review checks** compared the full histograms and automorphism counts, and tested every fixed-path tiling at these orders for the physical reflection formula, ordinary-converse projection, complement commutation and color transport. The frozen generator had previously obtained the same results independently:

| n | unguarded fixed histogram | weighted and colored guard histogram | weighted vertex automorphisms |
|---:|---|---|---:|
|3|{2:1}|{2:1}|2|
|4|{2:1,4:1}|{2:1,4:1}|1|
|5|{4:1,6:4,8:6,10:4,12:1}|{8:1,10:2,12:1}|1|

No n>5 run was made for this review. No frozen proof or generator was edited. Root owns the final full normal/optimized receipts and final repository checkpoint.
