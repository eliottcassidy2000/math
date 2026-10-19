# Neglected mechanisms and proof connections — 19 October 2026

Status: completed scoped proofs and exact controls, with cited inputs identified in each packet. Major conjectures remain open. This session follows [descent and global obstructions](descent_and_global_obstructions_session_20261019.md), with the user's emphasis on recently underexplored repository topics.

## The concept board

1. Induced arithmetic graphs: preserve nonedges as well as the prescribed square sums.
2. Graceful exchange fibers: retain the assignment of colors to fixed components, then change the allowed operation.
3. Oriented-tree parity: distinguish universal evenness from a useful odd-count existence certificate.
4. Ternary Pythagorean counts: retain the wall's spatial digits rather than fitting a recurrence to a long plateau.
5. Mahler carry words: strict real-tail safety and a positive ordinary integer source are different requirements.
6. Matroid tight sets: use their rank structure to understand whether recursive reconfiguration really needs exponential length.

## Proof packets

- [Induced square-sum forests](induced_square_forests_subdivisions_20261019/INDUCED_SQUARE_FORESTS.md): every n-vertex forest has an induced square-sum realization at height n^(2+o(1)); the exponent 2 is forced by stars. Exact finite seeds can be retained under leaf attachment. Induced cycles, countable forests and countable two-subdivisions have explicit constructions.
- [Simultaneous modular subdivisions](induced_square_forests_subdivisions_20261019/SUBDIVISION_MOD7_CONSTRUCTION.md): every finite graph's exact two-subdivision embeds induced at height O(K^4), K its subdivided vertex count. Residues3,5,6 modulo 7 remove within-family edges; adjacent-square gaps remove unintended cross-family edges. No unit-distance or chromatic conclusion is transported.
- [Factorial graceful fibers](graceful_factorial_fibers_20261019/GRACEFUL_FACTORIAL_FIBERS.md): explicit frozen upper tails have exactly k! one-color components, each a k-cube. Two-color moves give exact distance k−cycles(relative permutation)+ceil(f/2), with f the fixed-position sign errors, and diameter k−1 for k>=2. This changes tree shape and does not prove graceful existence for every abstract tree.
- [Pythagorean digit counts](pythagorean_counts_natural_boundary_20261019/TERNARY_COUNTS_NATURAL_BOUNDARY.md): c_n=2b_n+a_(n+1), strict discrepancy below 2/3, and the count tower recovers the wall digits. The regularized generating function has the unit circle as a natural boundary by the cited Pólya–Carlson theorem. No polynomial-coefficient linear recurrence holds, even after any fixed arithmetic-progression decimation. Transcendence of the density itself is not claimed.

- [Square-scaled pair extension](induced_square_pair_extension_20261019/SCALED_PAIR_EXTENSION_AND_TWO_DEGENERATE.md): every finite two-degenerate graph occurs induced. Every exact one-subdivision has labels below 4096n^4. The arbitrary-seed extension preserves adjacency after common square scaling and is sharp at two prescribed neighbors; seed 1,2,3 obstructs three via Fermat's area-one theorem.
- [Polynomial graceful reconfiguration](graceful_polynomial_diameter_20261019/POLYNOMIAL_STAIRCASE_DIAMETER.md): the rank-staircase diameter is at most s²−s+1. Arbitrary graceful upper tails with r >=3 free colors have diameter at most r²−3r+3 under two-color moves, replacing the former exponential sufficient bound. The proof uses tight suffixes and the cited strict-Hall transversal theorem. It is not a proof of graceful labeling for every abstract tree.
- [Twin-signature parity](twin_signature_parity_20261019/TWIN_SIGNATURE_PARITY.md): a minimal nine-vertex tree refutes HYP-9173(b)'s hereditary-symmetry converse; the rooted orientation-capacity recurrence detects the automorphism-forced-even branch exactly. A restricted tournament host universe has additional cancellations.
- [Two leaf-pair cancellations](twin_signature_parity_20261019/TWO_CHERRY_PARITY.md): two disjoint pendant leaf-pairs force even counts in every tournament, including asymmetric orientations. The retained isolated vertices supply the cancellation missed by the original automorphism filter.
- [Complete nine-tree classification](twin_signature_parity_20261019/NINE_VERTEX_GOOD_TREES.md): exactly 36 of 47 types are good; 12 are not hereditarily symmetric. Every nine-tree embedding parity is affine in orientation bits. In the good cases, the odd/even value is selected by specific leaf edges whose deletion is P8 or D8. The all-order classification remains open.
- [Padding complex](tournament_padding_complex_20261019/PADDING_COMPLEX.md): vertex deletion and universal-source adjunction satisfy d²=0 and dh+hd=1. This gives exact padding kernels, rank formulas and graph-embedding identities, while retaining the isolated-vertex information.
- [Mahler template closure](mahler_thue_morse_source_exclusion_20261019/THUE_MORSE_FINITE_SOURCE_EXCLUSION.md): every finite-prefix non-erasing binary block image of Thue–Morse misses each positive ordinary source, with an exact finite mismatch certificate. The [analytic extension](mahler_thue_morse_source_exclusion_20261019/THUE_MORSE_SCALAR_MAHLER_EXTENSION.md) gives rational sources for commuting blocks and transcendental sources for noncommuting blocks, with the stated uniform irrationality-exponent bound. Genuine strict-safe examples remain; arbitrary safe-state entry is open.

- [Sharp induced-book height](square_sum_book_height_20261019/BOOK_HEIGHT_ASYMPTOTIC.md): for a book with m pages, log H(m) ~ log(m) log log(m)/log(2). An explicit all-m construction matches the divisor-function upper obstruction. Equivalently, the largest induced book in [1,N] has logarithmic size asymptotic to log(2) log(N)/log log(N). The proof is elementary, with classical ingredients attributed.
- [Linear partition routing](graceful_linear_partition_routes_20261019/LINEAR_PARTITION_DIAMETER.md): the fully tight partition-matroid subclass has diameter at most 2r−3 for r>=3, with an explicit O(N+r²) actual-element router. Its graceful application requires the entire available contact graph to simplify to a tree. No optimal-diameter claim or new nonstar arithmetic realization is asserted.
- [Induced bipartite cleanup](square_sum_bipartite_bridge_20261019/SQUARE_SUM_BIPARTITE_CLEANUP.md): opposite translations after square scaling preserve all cross edges and remove all within-part edges. Combined with the cited Dujella–Elsholtz theorem, this gives induced K3,m for every m with varying anchors; no common total-height estimate is imported. A separate complete finite scan proves that the first K3,3 in the interval square-sum graph occurs at N=673, with an induced witness.

## Why these connections are useful

The square-sum construction separates two jobs cleanly: a small finite field suppresses entire classes of unintended edges, while an Archimedean size gap handles the remaining pairs. The graceful construction exposes an invisible permutation, then measures exactly how a two-coordinate operation repairs it. The ternary result turns an exact finite count into a recoverable boundary code and an analytic obstruction. None of these connections relies on a numerical coincidence between isolated constants.

## Current proof targets

- Remove the subpower loss in induced forest height: is a universal O(n²) bound possible?
- Narrow the quantitative gap for general two-degenerate graphs: books require the proved superpolynomial height, while the universal iterative construction currently has a much larger bound.
- The exponential graceful upper-tail bound is now replaced by a quadratic bound. Determine its sharp order and identify subclasses admitting linear-length paths.
- Extend the complete nine-tree parity classification to structural all-order families, retaining the padding and leaf-cut coordinates.
- Extend the finite-source carry obstruction beyond the specified substitution templates without discarding ordinary-integer entry.

## Precise differences between neighboring problems

- Forests have induced square-sum embeddings at near-quadratic height. Two-degenerate graphs all embed, but a universal polynomial height bound for that larger class is obstructed by large books (many common neighbors of one edge); the book-height theorem gives the sharp superpolynomial rate for that family.
- One-subdivision loses odd-cycle parity; two-subdivision preserves it. Neither preserves chromatic number or Euclidean unit-distance geometry.
- A fixed graceful upper tail can have factorially many one-color components and still have short two-color routes. The operation is part of the theorem.
- A parity-rigid graph can have constant zero. Only constant one gives a Rédei existence certificate. The nine-tree classification includes both kinds.
- A safe symbolic carry word can have a rational negative or transcendental 2-adic source. Safety alone cannot supply a positive ordinary source.

## Focused inherited repair

Incoming THM-4716 omitted the modular inverse in its changed mechanical-word index. The corrected index is p^(-1)(q−1−t) modulo q; the adjacent nonwrapping swap conclusion is retained. Only that local statement was audited. The session also repaired the HYP-9173(b) status and a downstream bullet that had treated its formerly open converse as proved; unrelated hypothesis parts remain unchanged.

## A useful distinction about repeated pairs

The quartic one-subdivision formula treats a simple graph: at most one new vertex is assigned to each original pair. Repeatedly requesting new common neighbors of the same pair creates a book. That repetition is exactly where divisor multiplicity imposes superpolynomial height. Thus the two construction theorems are compatible, and the numerical size of the labels carries genuine arithmetic information.

## Incoming concurrent work

The session preserved the external Farey/Mandelbrot and subsequent Sturmian-circle commits. THM-4716's local index repair is noted above. Incoming THM-4717 and the Gilbreath/duality hypotheses were recorded as external work with their stated scopes; this session did not independently audit those larger claims or use them as proved dependencies. The common language of circles, carries and symbolic words is not asserted to transport ordinary integer entry.

## Further leads retained, without new proof claims

- Glaudo–Kravitz, *Reconstructing a set from its subset sums: 2-torsion-free groups* (Discrete Analysis, 2024), supplies an established reconstruction theorem relevant to the repository's signed divisor kernels. [Primary paper](https://discreteanalysisjournal.com/article/125856-reconstructing-a-set-from-its-subset-sums-2-torsion-free-groups). Translating its full arithmetic consequences remains a separate task; no prime-generator completeness result follows here.
- A blanket exclusion of all automatic Mahler carry codes was not justified by the available scalar theorem. The completed result retains the exact binary block/Thue–Morse hypothesis.
- The original 482,3362,6242 elliptic common-neighbor construction was located in existing AP/congruent-number work and was not republished as a new result.

## Verification and publication

Each packet includes its proof, declared finite universe, runnable controls, and independent review where applicable. Normal and optimized Python runs retain active checks and matching receipts. Computations support the exact implementations; they do not substitute for the all-parameter proofs. Classical inputs are cited, and no exhaustive historical-priority claim is made.

The session was published as ordinary commits directly to main, preserving concurrent external commits. The final publication ledger verifies the final file bytes and lineage. The current-assignment links route here without replacing unrelated frontier statements.
