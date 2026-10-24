# Independent audit: integer move lattice versus binary packet connectivity

24 October 2026. **PASS**, with the native-source and dyadic-bank antecedents
retained explicitly. No producer file, checkout or Git state was changed.

## Findings

The proof in `move-lattices/THEOREM.md` is correct on its declared complete
28-source, common-target, capacity-one carrier. Its central distinction is
essential: all four- and six-slot trades span the entire integer charge
kernel, but those moves do not connect every binary fiber.

Independent ordinary and optimized runs each pass **894,666** explicit checks.
Their JSON receipts are byte-identical. Each separately invokes the producer
checker in the corresponding mode, which passes **10,184** checks; those
producer replay receipts are also byte-identical.

## Independent algorithms

1. **Actual labels.** All 28 five-coordinate source labels are re-read from
   the original electrical receipt, rather than accepted from search output.
   They are distinct, excluding any nontrivial two-slot trade.

2. **Complete move census.** Instead of the producer's equal-sum subset-pair
   enumeration, choose every changed support and every balanced sign partition
   with its first index positive. For support four this examines 20,475
   supports; for support six, 376,740. A support-parity sieve leaves 1,262
   and 23,636 supports respectively. Exact cyclic tests give **27 four-slot**
   and **1,917 six-slot** moves. Their full sets agree with the ordered list
   used solely to identify the saved integer witness's columns.

3. **Rank and invariant.** An independent rational reduced row-echelon
   calculation gives rank **23** for the four-slot moves. Every such move
   preserves the exact integer occupation sum b1638+b1974.

4. **Lattice certificate.** All 27 displayed basis columns are checked as
   exact integer combinations of the first 48 actual four/six moves, in all
   **28** coordinates, including the coordinate reconstructed from the zero
   margin. The supplied column matrix is upper triangular with diagonal
   128,2,2,2,2 followed by 22 ones. Its determinant is therefore 2048 directly,
   independently of the producer's Bareiss algorithm. Finite subgroup closure
   based at the first source, not its last, gives all 2048 values of the
   authenticated group Q.

   Let B be the certified basis lattice, M the move lattice and K the exact
   zero-margin charge kernel, all in Z^27. The certificates prove
   B contained in M contained in K and [Z^27:B]=[Z^27:K]=2048. Hence B=M=K.
   No unverified HNF/SNF computation is needed. Since the marked charge image
   is the whole known Q, the quotient type is (C2)^4+C128 as well.

5. **Binary fibers.** Rather than join every pair at Hamming distance at most
   six, enumerate actual local removal/replacement operations of one, two or
   three selected slots, keyed by exact charge, then use union-find. This
   checks 1,188 local replacements at cardinality two, 29,358 at cardinality
   three and 415,710 at cardinality four. The full component distributions
   agree with the producer:

   - k=2: 378 states, 351 nonempty fibers, all connected.
   - k=3: 3,276 states, 1,646 nonempty fibers, all connected.
   - k=4: 20,475 states, 2,048 nonempty fibers; 1,851 connected,
     175 with two components, 21 with three, and one with four.

   The replacement construction includes every admissible event of support
   at most six: remove its changed selected indices and insert its changed
   unselected indices. Conversely each enumerated replacement stays in the
   binary fiber and has support at most six. Thus it checks the intended
   legal-event graph, not merely a convenient subset of its edges.

6. **Exact hostile.** The full k=4, Q=(parity 7, cyclic 15) fiber consists
   exactly of the two V-selected source sets

       {176,1425,1549,2087},  {1394,1923,2837,3369}.

   They are disjoint and have Hamming distance eight. All 28 source choices
   were allowed during enumeration, so the other 20 supplied slots cannot
   buffer a six-slot route at these fixed margins. The direct eight-slot
   event is legal. This is a sharp bottleneck in this exact fiber, not a
   lower bound after adding outside sources or changing pins/margins.

## Proof-level review and native antecedents

The lattice index argument has both necessary inclusions and the correct
ambient zero-margin lattice. Full lattice generation implies no further
additive integer/modular invariant, but does not supply feasible binary
contexts for its coefficients. The theorem correctly does not infer a
six-slot Markov basis from its lattice generators.

The general legal-state vector-space argument is also correct: event
differences span exactly the vectors whose coefficient sum is zero on each
connected component, using paths to component roots. The underlying graph
here is the legal assignment graph, not the original electrical metagraph.
The square-zero warning is valid: multiplying a formal trade by a resource
already present on only one side can kill that monomial and leave a relation
to zero, which is not a two-state native move.

The native completeness/injectivity input is separately established by
`audit-compensation/AUDIT.md` and its receipts, inspected for this review.
It checks all 1,128,960 owner permutations across the 28 sources, finds 57
full pinned paths with distinct normalized resources, and retains actual
owner maps, source matrices, kernel lifts and complement/converse closure.
The sole extra path is source 1923 to target 1559; its zero required margin
excludes it from every accepted assignment, leaving exactly 56 allowed U/V
records. The original marked all-dyadic theorem and its separate fresh
filtration replay supply guard completeness. This audit does not duplicate
that million-permutation/native reconstruction or silently replace it by
charge arithmetic.

No odd-primary guard, arbitrary per-port character, common-host amalgam,
canonical full-atlas rewiring, global six-slot connectivity or all-order native
supply theorem follows. The original per-slot X mod 8 phase is separately
unchanged because the two actual targets have identical X.

## Reproduction

    python3 -B audit.py --repo REPO --producer PRODUCER --supply REPO/05-knowledge/results/session45_injective_multisource_moves_20261023/multisource-trades/target_pair_inventory.json --out normal.json
    python3 -O -B audit.py --repo REPO --producer PRODUCER --supply REPO/05-knowledge/results/session45_injective_multisource_moves_20261023/multisource-trades/target_pair_inventory.json --out optimized.json
    cmp normal.json optimized.json
    cmp producer.normal.json producer.optimized.json

PRODUCER is the directory containing the reviewed `check.py` and
`integer_basis_witness.json`. All replay outputs stay beside the requested
audit output. SOURCE_PINS.json records the reviewed theorem and independent
native antecedents as well as direct arithmetic inputs; MANIFEST.json hashes
the final audit artifacts, excluding itself.
