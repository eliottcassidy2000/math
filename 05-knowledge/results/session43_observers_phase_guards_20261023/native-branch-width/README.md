# Native branch width versus shared-state transfer

This extends Section5 of the Session43 native observer packet, following the root synthesis's regular-host and affine-cover sharpening. Read THEOREM.md for all proofs and scope restrictions.

On explicit regular tournaments of order3t, one genuine contiguous-path family has3*2^(t-1) paths and requires exactly2^(t-1) affine subgroup cosets in ANY exact union representation. Yet a3-state shared-phase transfer counts it exactly, with a period6 endpoint correction.

The prime-p version has p*|E|^(t-1) paths, exact coset width|E|^(t-1), and a p-state transfer. This is a representation-specific lower bound paired with a compact native compiler, not a computational hardness theorem and not a count of all host Hamilton paths.

Run from this directory:

    python -B check.py > receipt.json
    python -B -O check.py > receipt.optimized.json
    cmp receipt.json receipt.optimized.json

No dependencies, floating-point decisions, assertions, network calls, source mutations, or Git writes are used. The receipts retain actual-host checks, endpoint normalization, literal-path hashes, exhaustive small affine-coset checks, and integer period6 clocks.
