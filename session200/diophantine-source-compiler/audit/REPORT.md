# Independent compiler audit

Status: PASS for the frozen proof SHA256 `b22e1701e4656ee690426174efb9aa59f6827631d65abc1b2f2a83cbd4fd97f2` and verifier SHA256 `6b389a90c57f9bdcd21c6cac11d42ff7e7fea01da7b67a3c7b869b98e3ab5282`.

The audit read the physical-source and bilinear inheritance and checked the source-growth cutoff, quadratic one-hot encoding, induction forcing integer state coordinates over every characteristic-zero field, actual C1/N/M graph reconstruction, variable/equation counts, binary-input NP verification, and target-dependent large-prime exactness. The finite split source solution set introduces no additional field-extension witnesses. This establishes a family of systems generated from each target, not one fixed existential polynomial defining the entire spectrum. Universal odd coverage remains open.

Independent `check_graph_compiler.py` uses separately written graph births, Hamilton-path/cover subset dynamic programming, and graft construction. It checked all 13 source words of length at most two and all 169 ordered word pairs for each of K, N, M, with graphs through order ten. All source response coordinates and root-start counts matched. K decoded to C1 and N/M to C2. The Boolean-plus-orthogonality selector systems have exactly four solutions over each of F2,F3,F5,F7.

The author's verifier was rerun normally and with Python optimization. Receipts `compiler-replay-normal.json` and `compiler-replay-optimized.json` are byte-identical. Independent normal/optimized graph receipts are also byte-identical. Every audit gate raises explicitly rather than using removable Python assertions.

Before freeze, the author incorporated the requested fixed-formula caveat, provisional status, and separate balanced-source finite-field upper bound 6*4^L. No unresolved correctness issue was found. Audit checks inherit the complete source and C2 classification theorems; finite enumerations do not independently reprove those classifications.

The final proof promotes only its first-line status from PROVISIONAL to
PROVED after this PASS. All mathematical text is unchanged.

Reproduce from this packet directory:

    python3 audit/independent.py
    python3 -O audit/independent.py
    python3 verify.py
    python3 -O verify.py
