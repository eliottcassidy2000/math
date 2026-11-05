# Independent audit: complete C3 two-source reduction

**Result: PASS. No correction to the frozen source packet is required.**

Audit completed 2026-11-05. This review covers the all-order physical union, its exact scalar forms, the finite source inverse, and the supplied implementation. It does not infer universal odd coverage, uniqueness of a presentation, or a returning interface.

## 1. Snapshot and dependency scope

The reviewed source is `session197-hybrid-constructions/c3-two-source/`. The complete source manifest was checked before copying and again at closeout. Source files were not edited. Both Python modes were run on a separate copied snapshot because the verifier writes adjacent output files.

The six named structural dependencies are included in the source manifest. In particular, the review uses the corrected ROOT_AMALGAMATION proof, which allows its bundled suffix to be the singleton closer. The complete C3 degree-one theorem is an inherited, independently established all-order theorem with a justified finite reduction. The new packet does not try to prove that theorem from its own finite controls.

The dependency chain for the C3 bridge and start-three classification includes the exact published sources at commit `6ed2d7a0f1e1bfb26b9cc865ddc9f35dd0ad4567`, paths `session197/lowcycle/C3_COMPLETE.md` and `session197/lowcycle/root-amalgamation/START_THREE.md`. Their bytes agree with the named copies.

## 2. Mathematical review

### Arbitrary selected-entry graft

The scalar identity is valid with an arbitrary nonempty right graph, including a singleton. All backward cross arcs leave its marked root, so a Hamilton word has at most one such seam. In the one-seam term the only bad forward join is the stated singleton-root reserved path followed from a different selected old endpoint. This is exactly the ordered distinct-entry correction K times the deletion count. The cycle product needs a nontrivial right graph, as stated; it is not extended to the singleton case.

The injection cutting a Hamilton word before its first selected old vertex proves R_A-K_A >= h, not merely a numerical inequality on tested graphs. Together with E_A <= h and the two ending-root cover injections, this verifies both displayed graft lower bounds.

### Universal left response rows

All four rows were checked directly from their path-cover decompositions. In particular, the two directions of the K correction must remain separate until summed; no diagonal pair is included. The A/B clone-pair endpoint and two-cover cases, the six reserved-path cases for C, and the sole-outneighbor pair identity for T account for every term of the four rows. These identities hold for arbitrary rooted input graphs; source assumptions are only needed for the C1 interpretation and the selected start sum.

The mixed-pair classification produces A or B when the selected start-one vertex is one of the two old consecutive vertices. The third-vertex case forces the additional sole-outneighbor step and gives C. Deleting the first two vertices of the complete start-one triple gives T. Weak remainders are legitimate sources. A/B require nontrivial sources, whereas C/T include their singleton inputs.

### Triangle amalgamation and boundaries

The ordered root-amalgamation formula, applied twice to the four triangle incidences, gives exactly M0 through M3. Its endpoint formula has the explicit singleton exception; the scalar and cover identities do not need that exception. The transpose/converse relation is exact.

For two nontrivial sources every permitted Hamilton cycle traverses the two ordered source segments and the triangle, yielding exactly three choices. With one source singleton the cycle counts are 3-j or j. Thus precisely the two stated one-singleton boundaries are admitted. Both sources singleton are excluded. The graph construction and the scalar evaluation enforce the same boundaries.

### All-order completeness

After choosing an outdegree-one orientation using the inherited theorem, the SCC product produces the two stated cases. Recursive ordinary Type-I peeling terminates because it strictly decreases graph order. At a mixed/triple terminal branch, bundling all peeled material and the final closer leaves exactly one selected-entry boundary. The bundle contains both an entry root and the closer, so it is nontrivial. The cycle product then forces its ending count to be one; its converse is an actual source, not a formal vector.

At a Type-II terminal branch the later bundle is only required to be nonempty. If there is no ordinary suffix it is the singleton closer. This is explicitly handled by the corrected dependency and by the present proof. The two shared-root sides are still nontrivial. Their unique endpoint paths and their ordered cross arcs give an actual J_j graph. The direct initial start-one plus triangle-tail case supplies the j=0 singleton boundary, with its converse supplying j=3.

This is a terminating, surjective graph decomposition. It makes no injectivity claim and does not rely on agreement of finite graph inventories for the all-order step.

### Source cone, cutoff and inverse

For a nontrivial source the root has an outneighbor. Inserting it just before the first outneighbor in each deletion Hamilton word injects into non-ending-root words, proving D-E >= d. Separately, cutting each full Hamilton word before the root injects into the U covers. Its intersection with the d covers whose reserved path is just the root has size exactly E. This gives U >= D+d-E, and hence U >= 2d. These are different injections; the latter is the cut-cover image used in the overlap argument.

The nonnegative coordinates (D-E-d,U-D,E,d) sum to U. Direct substitution gives the four stated nonnegative P matrices. Their entry bounds, together with d <= U/2 for P0/P3, establish the Type-II product bounds. The two boundary cases give H >= 6U. The graft inequalities and the left scalar lower bounds give the Type-I cutoffs, including the separately permitted singleton left sources. Therefore every source in every allowed representation obeys U <= floor(H/3).

The source oracle includes weak sources and retains physical recipes. Its U coordinate strictly increases. Merging equal profiles preserves every future profile and every evaluated scalar, without identifying the underlying graphs. Restricting to the cutoff therefore loses no representation. Exhaustive admissible ordered source-pair evaluation is an exact membership algorithm for the strong C3 path spectrum, with actual witness reconstruction.

The stated complexity is deliberately numeric-polynomial: O(N^4) profile states and O(N^8) pair arithmetic, with O(log N)-bit fixed-degree intermediate values. It is not a polynomial-in-bit-length claim. All singleton and small-target cases remain governed by admissibility, rather than by extrapolating inequalities to an empty or invalid graph.

## 3. Independent numerical controls

The standalone `check.cpp` imports none of the author's Python implementation. It uses forward subset Hamilton-path inventories on both orientations to compute endpoint counts and ordered path-cover responses, and literal adjacency constructions for all templates. It independently canonicalizes strong graphs by their Hamilton cycles and source graphs by their unique marked starting Hamilton word.

The optimized and assertion-enabled unoptimized runs produce identical receipts:

- 21,620 universal response checks: all 5,405 rooted labelled tournaments through order five, with all four left templates.
- 29,667 arbitrary selected-entry graft checks, including singleton right graphs for the scalar formula and only nontrivial right graphs for the cycle formula.
- 1,694 physical source-cone checks.
- Complete rooted start-one graph-class counts through order seven: 1, 1, 2, 5, 13, 34, 89.
- 2,022 literal admissible template presentations of output order at most ten.
- Exact agreement of the resulting graph union, including converses, with all C3 classes at orders five through ten: 1, 11, 38, 124, 390, 1193, totaling 1,757 classes.
- Every one of the 6,040 labelled C3 graphs through order six belongs to that union.
- A separately generated source bank at U <= 167 has 811 states. Its gadget constants use a direct binary formula and reverse-complement response, with separately reconstructed graph recipes.
- 5,252,040 admissible ordered pair tests were performed without the source implementation's target pruning. They reproduce exactly all 217 supported counts through 501.
- All 217 retained witnesses were physically reconstructed and independently checked for their exact H and C=3. Their largest order is eleven.

The unlabelled-through-seven input is the complete uncapped tournament bank, not a cap-selected source inventory. The C3 graph fixture is the restriction of the previously independently matched complete C<=3 census to C=3 and order at most ten. The expected scalar support is used only for final set comparison, not to generate any source, candidate pair, or witness.

An initial local census-fixture preparation accidentally retained C1/C2 rows as well as C3. The independent checker rejected the fixture immediately. The final fixture explicitly filters both C=3 and order <=10. This was a reviewer input-preparation correction, not a defect in the frozen source or its mathematics.

## 4. Source replays and reproduction

Both copied source runs, `python3 verify.py` and `python3 -O verify.py`, reproduce the frozen validation bytes. The regenerated SUPPORT.json also exactly matches the frozen source support. The independent optimized and debug receipts agree byte for byte.

Run the independent checker from the directory containing its three fixtures:

    g++ -std=c++20 -O2 check.cpp -o check
    ./check > independent-normal.json
    g++ -std=c++20 -O0 -D_GLIBCXX_ASSERTIONS check.cpp -o check-debug
    ./check-debug > independent-debug.json

Required independent files are `check.cpp`, `unlabelled-through7.csv`, `c3-through10.csv`, and `expected-support.txt`. No Python source module, external package, network, or repository checkout is needed by this checker. The two independent receipts and the two source replay receipts should accompany this report. Local compiled binaries are not part of the audit deliverable.

Fixture provenance:

- `session197-prime-atoms/source-bank-order7.json`, SHA-256 `a441288c50604918fd20e1a3a6e87a8eb85ddb63f27c81b618a64cbf2646a4b3`. Complete tournament counts at orders one through seven: 1,1,2,4,12,56,456. The CSV records the order followed by its adjacency rows.
- `session197-bridges/c3-source-records.csv`, SHA-256 `1d0159168c03fcd011e90caab5f0fe1dd170362ba414d2f74730fe9d644b480a`. Restrict rows to order <=10 and C=3; retain the order,C,adjacency-row format. This census was independently matched before the present review.
- `expected-support.txt` records the sorted finite support in the frozen validation. It is not an input to the source-bank generation or physical graph construction.

## 5. Frozen file binding

The following hashes bind this review. A PASS applies to these bytes and the stated dependency scope; later edits require a fresh comparison.

### Source packet

- `PROOF.md`: `103a0c347881c7641f17d53bb8796187ac80f057a90af0a3acf8e6f89032c64d`
- `core.py`: `dc32dc1ef3647b694015cb4ff1ef99415c162bfbaf58d3b545175b09f0378042`
- `source_oracle.py`: `33cfe6cea487b17da44e083d7e81bd860dfeddcab3096deb6d135191e5537ac0`
- `c3.py`: `564034e130625d8bfd06ac04f7f365c14c93cd22e9ef8f33da3f5330a0825bb9`
- `verify.py`: `5fcd2e5ec2ae66c3e5e0d8c2d7312fd112608f7cf13fb4a59ee12d8d3f97c986`
- `SUPPORT.json`: `fc8f2eb1bd64634aa34927302105583914bc13812b1433fb954bd3917e2dd943`
- `VALIDATION.json`: `2f05815a209685996bcdb9256cf523d74843cb82e246ec870eafb13335f1c84f`
- `VALIDATION-normal.json`: `2f05815a209685996bcdb9256cf523d74843cb82e246ec870eafb13335f1c84f`
- `VALIDATION-optimized.json`: `2f05815a209685996bcdb9256cf523d74843cb82e246ec870eafb13335f1c84f`
- `dependencies/BILINEAR_C2.md`: `5eab9deb2871819d8a1b5684c4987f2f4da807adc55063a30675d531c5ad3813`
- `dependencies/C3_COMPLETE.md`: `a919deb88bc78a84b8d0cf40199fc7966051541f1e7027be3666bb9f0713b444`
- `dependencies/ROOTED.md`: `ba7e632518219eeb8d86b6db7da386a88e98317ea417c9b22e649278c132b32e`
- `dependencies/ROOT_AMALGAMATION.md`: `b773da29cfb8298f3df4b99a741dccad7618d1bc84165f69e9066d86b819f072`
- `dependencies/START_THREE.md`: `5b0010444d0339741c1ad3ece83e3a30511457d9c3e9f74f037af18510887753`
- `dependencies/START_TWO.md`: `a746d2e33b0fb78a98252afa602dc8985c25f718225c2f9ce37923fd20b0ada3`

### Independent checker, fixtures, and replay receipts

- `check.cpp`: `a6073ec6885cede0904b894b073bed397ff3b01ad1a5e3b90f568a7cc8fa7f84`
- `unlabelled-through7.csv`: `11d184c47fbccca0c6d7bddd45fc52a60ccc557d3989c323f36114e0d6af365c`
- `c3-through10.csv`: `5b571b5d2ebf7edaf6af94363caf88328088f984d840e7bf5878b813b63cfe80`
- `expected-support.txt`: `c2115240d601f7f9618d2dc1d3d99318dd02b7e87bc1389fd8e9a4d1e2bd9d1f`
- `independent-normal.json`: `160398039e6c8b9ee691759f64930604dabc2baa2edff2a898ec2e6a4368bdec`
- `independent-debug.json`: `160398039e6c8b9ee691759f64930604dabc2baa2edff2a898ec2e6a4368bdec`
- `source-normal.json`: `2f05815a209685996bcdb9256cf523d74843cb82e246ec870eafb13335f1c84f`
- `source-optimized.json`: `2f05815a209685996bcdb9256cf523d74843cb82e246ec870eafb13335f1c84f`

## 6. Conclusion

The frozen complete C3 physical reduction and exact finite inverse pass independent review. No source correction or unresolved blocker remains. The finite support through 501 remains a finite result; neither it nor the all-order reduction proves the unrestricted {7,21} spectrum conjecture.
