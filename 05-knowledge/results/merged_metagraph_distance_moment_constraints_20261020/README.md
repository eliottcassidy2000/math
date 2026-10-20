# Merged-metagraph distance moments

This packet proves necessary matrix-valued constraints for the Hamming-distance line layers of the converse-merged Hamiltonian-path tiling metagraph.

Main results:

- Walsh transforms of the distance matrices are integer Gram matrices, with reflection-refined rank bounds.
- Alternating even/odd-distance counts factor through a nonzero signed tiling-parity charge at every merged node.
- The parity charge equals a sign-corrected signed Rédei count divided by the automorphism count, with the appropriate converse-sheet multiplicity.
- A block covariance identity quantifies exactly what is lost by compressing single-tiling moves to nodes.
- An explicit order-five rewiring preserves all layer degrees, diagonal distance distributions, entry parities, pair totals, and the entire complement graph, but fails a two-by-two positive-semidefinite test.

These are necessary conditions, not a complete realization criterion or a unique-arrangement theorem.

Read `PROOF.md` for the definitions, proofs, exact scope, and failure boundaries. Run:

    python3 verify.py

Requirements: Python 3, NumPy, and SymPy. The program is independent of repository computation modules and runs from any working directory. It exhaustively checks every fixed-path tiling at n=3,4,5,6, compares direct distance counts to direct Walsh sums, checks the signed-Rédei identity on all 49 merged nodes, and reproduces the impossible rewiring.

Outputs are `audit.json`, `audit_n6.json`, `signed_redei_check.json`, `hostile.json`, and `hostile.out`. `verify.out` records the successful run. Mathematical dependency fingerprints are in `DEPENDENCIES.json`; these document inherited statements and distinguish them from the new derived constraints.
