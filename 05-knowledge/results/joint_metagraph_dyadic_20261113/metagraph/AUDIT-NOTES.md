# Local audit correction lineage

The root audit identified a wording error in the unpublished proof candidate.

- Before, PROOF.md SHA-256: `0473918265fc5e5a8f72ea97b41d70ec92faaf959b5f94193cd44fa9ca4028b5`.
- Incorrect sentence: “it is not uniform over physical paths in one representative unless the orbit and converse factors are included.”
- Correct statement: uniform physical Hamiltonian-path sampling in one ordinary representative induces exactly the same **normalized merged-node transition law**. The automorphism action is free on spanning paths, hence all path orbits have the same size; converse-sheet duplication also cancels in the ratio. Raw adjacency counts and actual dart identities still require those factors.
- Mechanism: confusing the normalization-sensitive raw count with its factor-canceling probability ratio. No matrix computation or numerical claim changed.
- After correction, stronger equal-weight hostile and candidate-status clarification, PROOF.md SHA-256: `1eb44268e5bbbfc1f374f97281cf4def000d6dab6c3bcd72352412d17874d9db`.

The root independently checked the all-order triangle twist, small native matrices, weighted Gram identity and two-step return. The final Krylov/source audit and status promotion belong to the root; this file does not claim them completed. The candidate remains outside canon.
