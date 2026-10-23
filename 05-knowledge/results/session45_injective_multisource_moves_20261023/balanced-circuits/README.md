# Large primitive native cancellation circuits

Read THEOREM.md for proof, native scope, the sharp sixteen-slot boundary in the fixed twenty-option supply, and the Boolean-assembly obstruction.

From the repository root, run:

    python3 -B 05-knowledge/results/session45_injective_multisource_moves_20261023/balanced-circuits/check.py --repo .
    python3 -O -B 05-knowledge/results/session45_injective_multisource_moves_20261023/balanced-circuits/check.py --repo .

The maximum-support and Boolean projection checkers locate their inputs beside themselves and may be run from any working directory:

    python3 -B 05-knowledge/results/session45_injective_multisource_moves_20261023/balanced-circuits/max_support.py
    python3 -B 05-knowledge/results/session45_injective_multisource_moves_20261023/balanced-circuits/boolean_boundary.py

Repeat each with -O for the saved optimized controls. Their saved normal and optimized receipts are identical. The independent audit uses different graph, subset, maximal-support and Boolean algorithms.

image_witnesses.json retains seven full integer image certificates and their corrections. The optional solve_images.py accepts --repo and --out; use a separate output path and compare it with the frozen witness file. Both ordinary and optimized regeneration were completed and matched byte-for-byte; regeneration.json records the digest. The verifier does not trust the elimination algorithm.

search.py is optional deterministic discovery provenance. It writes candidates.json; it is not part of routine read-only validation. The supplied cases are independently authenticated by the checkers.

The current guard is nodewise and dyadic. Chosen ordinary per-slot X phases are also retained, but the packet does not preserve every possible ordinary per-port phase, odd-primary character, or an unqualified common-host embedding. Original tiling capacity belongs to normalized resources; owner lifts remain explicit witnesses.
