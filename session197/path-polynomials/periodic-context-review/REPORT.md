# Independent review: periodic whole-module context dichotomy

**PASS, no correction required.** The exact scope is a fixed finite word of rooted tournament contexts, each replacing its root by the entire current nonempty tournament as a uniform module. Repeating such a word is either exactly geometric, or has the stated superexponential lower bound. This does not constrain exposed-port births, adaptive words, root changes, graph rebuilding, or the unrestricted spectrum.

## Frozen source and dependencies

Original directory: `session197-odd-cycle-path-polynomials/`. All source files remain unchanged.

- PERIODIC_CONTEXT_DICHOTOMY.md: `421d2b57adbd891e79fa773848eb42efcafbb738a7276fa91078e7c120f2c4da`.
- verify_periodic_context_dichotomy.py: `b9ba608cd8ec3ca8666c1bf2c0def97f5496ead6bcd6ed2c9e17170dc9755cf0`.
- Both frozen receipts and both copied-source reruns: `152e647712779b3a77cfbf8a97a2582f7bf41d6a2700a4f4a087eb0908c56eb1`.
- UNIVERSAL_CONTEXTS.md: `8615845e7bf06a6cdb1925fdbd8198a77b1a8864aa9c22704b0459678ff2ed51`.
- Its companion PROOF.md: `b88c86ac0afa64d3a1bf399c56fa7f8aad046fb93424c70a9ace8d164708183b`.
- verify_extension.py: `3b537a9cbfea66aabe8b1d4cf59c84de45cc41d2ea6f80d053df8c083b1cad76`.
- Companion verify.py: `10c2a35ce3f506aefb031645e81cffc1a2fb4d73bf9bbb38a163ccbf8e35dfcb`.

Normal and Python `-O` source replays agree exactly. The later adaptive supplement is not part of this review and is not a dependency.

## Proof audit

1. **Separate skeleton positivity is available.** The named universal theorem proves more than positivity of the sum or physical outputs. Each packing contribution is a nonnegative composition of X, Delta and delta on the extended central-factorial basis. The recurrences for the two parity kernels have nonnegative integer coefficients. Thus it is legitimate to retain two individual skeletons and discard all others entrywise on virtual columns.
2. **The selected triangle has the claimed coefficient and ordering.** The empty exterior packing contributes X^N. A directed root triangle has two exterior vertices, one marked gap, odd assigned-run parity, and N−2 unused exterior vertices. Its contribution is exactly X^(N−1) delta, with coefficient1. The shifts act first and X afterward; the proof does not commute them. Their sum is X^(N−2)(X^2+X delta), with N>=2. These are distinct packings, including for the smallest triangle context.
3. **The physical profile starts in the required cone.** Odd-cycle packing polynomials have nonnegative monomial coefficients and the graph's parity. Repeated multiplication by X is positive in the extended basis, so physical profiles are nonnegative central-factorial vectors. Nonempty graphs have no constant term. Evaluation at1 selects precisely their index0 coefficient, equal to H. Hence every input dominates its own H times the proper parity vacuum, even when higher coordinates or Hankel behavior are arbitrary.
4. **Cross-parity domination is valid.** Each context shifts parity by its number of exterior vertices. Selecting X^N or X^(N−2)D has the same shift. Matrix multiplication is therefore correctly typed even for an odd-size period. The proof uses the two appropriate alternating matrices, not an unjustified fixed parity-return matrix. Every individual column is finite, so no infinite-series convergence is involved.
5. **Every active period has the three required transitions.** Both D parity blocks have upward weight1, positive downward weight for k>=1, and diagonal at least(k+1)^2. X has a positive diagonal in both cross-parity directions. Choose one D for the required move and all other factors for their diagonal moves; after a move, those remaining diagonals are taken at the new index. This yields the stated period entries in either starting parity.
6. **The path bound handles all boundaries.** For 0<=k<=floor(m/2), climb for k periods, stay for m−2k, and descend for k. This gives h_m>=h(k+1)^(2m−4k). Taking k=floor(m/4) gives the claimed h(1+floor(m/4))^m, including m=0 and m=1,2,3. Its m-th root tends to infinity.
7. **The passive alternative is exact.** No root triangle means every cross edge is I→O. The substituted graph is literally the ordered join Q[I]→T→Q[O]. Empty sides have H=1. All counts multiply, allowing cyclic exterior components and the empty period. Thus h_m=h c^m with a positive integer c. If this sequence is nonconstant and eventually satisfies an affine recurrence, subtracting consecutive equations forces multiplier c and additive term0. Active growth cannot satisfy any fixed affine recurrence because such a recurrence has at most exponential absolute growth.

No step assumes that a virtual basis state is a realizable tournament, that every positive operator is physical, or that every construction is a whole-module substitution.

## Independent numerical controls

`check_physical.cpp` imports none of the source implementation. It independently builds actual transitive-clone graphs, enumerates directed odd cycles and disjoint cycle packings for their complete bin polynomials, and checks evaluation at1 against a separate endpoint Hamilton-path DP. It then obtains virtual-column images by ordinary linear combinations of those actual monomial images, rather than using the author's cyclic-skeleton compiler or central-factorial transition rules.

- All285 rooted labelled contexts through order4.
- 1,710 virtual odd/even columns, including the exact triangular remainder inequality.
- 810 passive virtual columns with exact full-profile multiplication.
- 396 distinct actual bin-profile graphs, each with an independent endpoint-H comparison.
- All63 literal repeated-context counts in the source's21 histories, including parity-switching active periods, a passive exterior triangle, identity and empty periods, and three different nonempty seeds.

The independent C++ optimized and assertion-enabled unoptimized outputs are identical and reproduce every source physical history.

`check_transitions.py` is separately written using integer polynomial arithmetic and no source imports. It checks124 exact D/X band identities through index30, all2,160 parity/index controls for D/X words of length at most6 containing D, and560 longer iteration states. The long controls deliberately begin with nonphysical positive higher-coordinate vectors as well as positive vacuum coefficients. Every climb/stay/descend bound is checked, not only the chosen floor(m/4) bound. Normal and `-O` outputs agree.

Independent hashes:

- check_physical.cpp: `dff397c065647d9aec6c6c6aedd3b4dc5705fb3a5a38c141930549a1c625642c`.
- independent-normal.json and independent-debug.json: `e8a6a505aea6bf103bed7d878442511c850bcd6b2bba39a03b4a2ccb0896d209`.
- check_transitions.py: `4285122c0253a58d6493a5809535a120e086c1e8d531aaa513392dd63731f04e`.
- transitions-normal.json and transitions-optimized.json: `0b570f8fa4a3a325f55b4ccf92904715dd083deed35b28aa194305a8756e69f9`.

## Reproduction

Run the unchanged source checker with its two included executable dependencies, normally and with Python `-O`. Compile `check_physical.cpp` as C++17 and save its output as `independent-normal.json`. Save a source replay as `source-normal.json`. Then run `python check_transitions.py` and its `-O` form; the script also compares all physical history arrays. All independent mathematical computations use only standard libraries.

Conclusion: the fixed whole-module periodic dichotomy and the exclusion of nonconstant affine clocks with nonzero additive term are sound at precisely the stated scope.
