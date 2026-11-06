# Two explicit colored repair contexts

Status: **PROVED, independently audited in the stated scope.**

For each of the two physical parent/child core pairs in PROOF.md, an arbitrary two-colored tournament body T satisfies H(child) >= 2 H(parent) + H(T). Equality holds exactly when every color-1 vertex beats every color-0 vertex. On a transitive spine this permits only words 1*0*, so these two promoted repairs do not give a free binary encoder. Nontransitive ordered color blocks do attain equality.

This is a classification of these two displayed operations, not all two-root contexts. It does not establish an eventual odd-spectrum threshold, a globally missing count, or strong outputs. The B construction always uses its full eight-vertex child; a sink deletion from the monochromatic construction is not valid for mixed colors.

## Contents and provenance

- PROOF.md is the audited proof with a status-only promotion. REVIEWED_PROOF.md preserves the exact reviewed input, and STATUS_PROMOTION.diff records the complete change.
- verify_colored_identity.py, NORMAL.json, OPTIMIZED.json, colored_tensor.cpp, COLORED_TENSOR.txt, and PROMOTED_TWO_ROOTS.json are unchanged producer artifacts. REVIEWED_INPUTS.json preserves their original input manifest; its PROOF.md digest refers to REVIEWED_PROOF.md after the documented promotion.
- audit/ contains all 15 entries of the independent audit manifest and SHA256SUMS itself, unchanged. REVIEW.md states the independent checks and scope. The audit directory is standalone.
- The proof includes the cut identity, strict positivity, and the fixed-core gap bijection. No external tournament bank or unproved completeness statement is a proof dependency. Earlier repair packets explain the origin of the rows, but are not needed to reconstruct or check this result.

The producer checked 198 literal parent/child pairs. The independent audit recomputed all 51,120 core orders, both convolution factorizations, 1,099 colored cut inputs, and 204 literal graph pairs. The finite core exhaustion certifies an all-order identity; the finite body controls are additional checks, not an extrapolation.

## Reproduction

From this directory:

    python3 verify_colored_identity.py > replay.json
    python3 -O verify_colored_identity.py > replay_optimized.json
    cmp replay.json NORMAL.json
    cmp replay.json replay_optimized.json
    g++ -O3 -std=c++17 colored_tensor.cpp -o /tmp/colored_tensor
    /tmp/colored_tensor > ordered_replay.txt
    cmp ordered_replay.txt COLORED_TENSOR.txt

From audit/:

    sha256sum -c SHA256SUMS
    python3 check_tensors.py > tensor_replay.json
    python3 check_bodies.py > body_replay.json
    cmp tensor_replay.json TENSORS_NORMAL.json
    cmp body_replay.json BODIES_NORMAL.json

The Python verifiers use explicit failure checks, including under optimization. Running the producer writes the redundant VERIFICATION.json next to its script.
