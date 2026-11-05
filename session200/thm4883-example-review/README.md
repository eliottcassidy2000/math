# Scoped review of one THM-4883 example

**VERIFIED COUNTEREXAMPLE:** the unqualified example “every J-born source over a strong source” does not always give all odd C2 offsets. The exact T4 offset identity and its gcd-qualified residue statement survive.

See [REVIEW.md](REVIEW.md) for the proof, repair, source citations and limits. This standalone packet does not modify the original theorem, its proof, or any index. It is not a full audit of THM-4883.

Reproduce from this directory:

    python3 verify.py > replay-normal.json
    python3 -O verify.py > replay-optimized.json
    cmp replay-normal.json replay-optimized.json
    cmp replay-normal.json VALIDATION_NORMAL.json
    python3 audit/independent.py > audit-replay-normal.json
    python3 -O audit/independent.py > audit-replay-optimized.json
    cmp audit-replay-normal.json audit-replay-optimized.json
    cmp audit-replay-normal.json audit/INDEPENDENT_NORMAL.json
    sha256sum -c MANIFEST.sha256

The portable standard-library verifier builds and counts the actual sources and both C2 grafts independently by subset dynamic programming. It checks all gates with explicit exceptions, including under Python optimization.

The second audit uses direct permutation and path enumeration on the published adjacency rows, with no subset-DP code or verifier imports.
