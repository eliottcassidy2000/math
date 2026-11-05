# Reviewed unbounded weak deletion ambiguity

For every m, the [proof](UNBOUNDED_WEAK_AMBIGUITY.md) constructs at least m rooted graphs of the same order and weak response triple, but with distinct root-deletion counts. The [independent audit](audit/REPORT.md) passes without correction.

This is a separate corollary of the balanced power-two gadget law and the weak-source decomposition. It does not claim an exhaustive list of deletion values, a higher interface-rank requirement, or scalar coverage.

Run the source wrapper normally or with -O:

    python3 verify_weak_ambiguity.py
    python3 -O verify_weak_ambiguity.py

Each run compiles the standalone C++17 counter in a temporary directory in optimized and debug modes, compares their results, and writes weak_ambiguity_certificate.json. It directly counts10 graphs through order16 and separately checks1,275 symbolic recipes through m50. Use a copy if preserving a pristine checkout matters.

From audit/, run:

    python3 independent.py
    python3 -O independent.py

The exact copied source receipt and replay layout are included. This independent construction and recursive counter verifies the10 graphs, adds literal permutation controls, and checks all5,913 gadget permutations through m7. The proof dependencies are included; earlier full inverse audits are published in the sibling backward-profile-inverse/ and weak-start-one-oracle/ checkpoints. All executable dependencies use standard libraries. FILES.sha256 binds every released byte.
