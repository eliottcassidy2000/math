# Independent review of C3 automorphisms and labelled enumeration

**PASS, no correction.** The group action, canonical symmetry recursion and labelled orbit mass are sound. This is a separate review from `session197-start-three-audit/REPORT.md`, which covers the underlying start-three grammar, exceptional subclass and unlabelled C3 enumeration.

## Frozen source

Under `session197-lowcycle/`:

- `C3_AUTOMORPHISMS.md`: 0f38f6d2989275bb6f7271cea0c8aafb94c05d99c7d4ca73760025859eaf944f.
- `verify_c3_automorphisms.py`: e7538f50f723f03c84674c19263a734acd836fc4533185404ade9c8cffb774a2.
- `c3_automorphisms_certificate.json`: 26adef3f570fd06b8b5868d3e26d2c7c80f95a1d4f1f7f1b1b073ec66254def1.

The verifier's named C++ census dependency is `verify_c3_complete.cpp`, SHA e49f0b9ff6e29f7b7260d6516ccb0fafe03ea4bbcf06dccc0664b7fddead806f. The rooted word helper is the pinned `verify_start_two.py` from the preceding audit. The source was copied before replay, since it writes a certificate beside itself. Both normal and Python -O runs match the frozen certificate exactly. No original source was edited.

## 1. Full automorphism group and rooted stabilizers

The complete C3 theorem supplies a degree-one vertex in at least one orientation. There are one or two such vertices in that orientation, and every automorphism fixes them individually; exchanging a pair would reverse its mutual arc. Fix one such vertex v. If an automorphism fixes one Hamilton cycle as a cyclic object, its fixed v determines a position in that cycle and forces every vertex to be fixed. Thus every cycle stabilizer is trivial. The group acts freely on three cycles, so its order divides three. This proves precisely the two cases 1 and 3, not merely a bound on an unspecified subgroup.

For a rooted graph with three start paths, its root stabilizer likewise acts freely on those three words, so it has order 1 or 3. The source verifier counts nontrivial rooted stabilizers only among C3 graphs; this is valid. A nontrivial order-three root stabilizer would also act freely on the graph's Hamilton cycles, while 1<=C<=s_r=3. Hence C must equal three. No C1/C2 rooted symmetry case is silently omitted.

Rooted start-one inputs have trivial stabilizer even at the triangle boundary. The unrooted triangle's rotations are therefore relevant only for a separate unmarked tail triangle or a whole-triangle triple entry, not for the initial marked start-one component.

## 2. Symmetry recursion

In the single-entry Type-I branch, every output automorphism fixes the new root and the ordered tail, and restricts to exactly an automorphism of the marked input. Conversely every marked-input automorphism extends. The mixed-pair input is rigid. A triple-entry input has symmetry only at the triangle, and then the mandatory nonempty tail leaves all three rotations intact. These give K R and z^3(K-z).

In Type II the initial rooted input is rigid. A tail triangle rotates if and only if its incidence to the new root is uniform. For a terminal triangle only the all-zero subset is admissible; for a nonterminal triangle both all-zero and all-one are admissible. The ordered singleton positions and SCC order are fixed. At total tail size k this gives (k-2)2^(k-3), so the tail series is z^3/(1-2z)^2. This proves

    R=z^5(2-4z+z^2)/Delta^2.

For exceptional roots, only the inherited one-entry symmetry and triangle triple-entry branch remain. The unique order-four mixed-pair input is rigid, and all triangle-tail inputs fail the exceptional condition. Thus

    R_E=K R_E+z^3(K-z)=z^5/Delta.

The distinctions between full groups and marked stabilizers are maintained at every step.

## 3. Unmarking and labelled mass

Every automorphism fixes each degree-one vertex individually, so marking such a vertex does not reduce the full group. The two-degree-one insertion/deletion bijection restricts and extends automorphisms exactly: the new vertex is fixed, and its incidences are determined by the old marked vertex. This holds both for the entire degree-one class and for the subclass with no opposite degree-one orientation.

It follows that the same marked/unmarked correction as in START_THREE is valid on the order-three automorphism class. Adding the converse-only orientation gives

    S_C3=z^5(1-2z)/Delta^2.

A rigid unlabelled class contributes n! labelled graphs; a class with three automorphisms contributes n!/3. Therefore the labelled count is exactly n!(c_n-2s_n/3). The coefficients 40 and 6000 at orders five and six are correct, and the source's normalized labelled generating function follows by subtracting two thirds of S_C3 from C3.

## 4. Independent controls

`check.cpp` uses the already independently certified complete low-cycle bank through order twelve, `session197-bridges/c3-source-records.csv`, SHA1d0159168c03fcd011e90caab5f0fe1dd170362ba414d2f74730fe9d644b480a. Its completeness and graph-by-graph comparison are documented in `session197-bridges/C3_COMPLETE_CENSUS_REVIEW.md`.

Unlike the source's candidate generator extracted from three path words, the new checker enumerates every automorphism by isomorphism-invariant degree refinement and full compatible-bijection backtracking. It counts start paths by reverse subset DP, explicitly forms each eligible root orbit, and counts its full stabilizer. It verifies that every nontrivial start-three stabilizer lies in a C3 graph and that every C3 full group has order one or three. It tests all degree-one vertices are fixed and all elements of each nontrivial C3 group have cube identity.

Across all 27,945 graphs it examines 37,019 automorphisms. The sum of stabilizer sizes over all start-three root orbits is 63,865. All three generating-function sequences match:

- Nonrigid C3 classes at orders five through twelve: 1,4,13,40,120,354,1031,2972.
- Start-three root classes with stabilizer three: 2,8,27,86,265,798,2362,6898.
- Exceptional rooted classes with stabilizer three: 1,3,8,21,55,144,377,987.

Direct orbit mass gives labelled counts 40,6000,147840,3924480,112492800,3472761600,115226496000,4098178022400. Optimized and unoptimized C++ with library assertions produce identical receipts. The prior start-three audit additionally checked literal Hamilton words and all full vertex permutations for every tournament through order seven; the copied source independently checks the complete labelled universe through six.

## Artifacts

- `check.cpp`: d0409c03c017c0c05b585c8fc8bcdf5629ae34ca3fd7012f3b3339f1ae80c3e6.
- `independent-normal.txt`: 955e30544aa12ffbf4cb1dcf1da0cefcd1d0b73281befb18f04f4d058721bfdd.
- `source-normal.json`: 26adef3f570fd06b8b5868d3e26d2c7c80f95a1d4f1f7f1b1b073ec66254def1.

No assertion about automorphism groups at higher cycle counts is inferred, and no blocker remains.
