# Independent audit

Reviewed `NOTE.md`, `check.py`, and `independent_matrix_check.py` on 2026-11-03. No substantive defect was found.

The least-significant-digit sign recurrence matches all four cases of the specified Cayley–Dickson doubling convention. Its six-state reachability/minimality argument is correct. Signed-XOR permutation composition gives exponent at most four and the stated exact cycle multiplicities at every fixed level without an associativity assumption. The local-field valuation proof gives precisely the stated 26·3^j clock periods and state counts. Marked-site occupation gives exactly k triggers per clock lap, yielding the sharp 26k, 52k, and 104k families. Three-row observability has determinant two modulo three.

The tournament encoding is injective for the prescribed family with fixed marker and decoding convention; the full profile is constant. The note correctly distinguishes synchronous changes, externally scheduled class-flip paths, literal physical transitions, and complement-line metagraph operations. It does not claim an autonomous unmarked tournament realization. The inherited strong-wrapper result was not independently re-proved in this audit.

Normal and optimized runs of both provided programs passed byte-identically:

- Main output SHA-256: `7498ba766d5ab18ff398ef2237b5ee44716dd00a31c87d6818f8aefbfd1057d2`.
- Matrix output SHA-256: `d7e20dc4ab8cb8f1be59a983de164165faeaffe28bd9337a6271255655efcd05`.
- Audited note SHA-256: `ac2fcc4b7a0c12994f3de6efa5edf1d44a89c297df434bffef321b4d696b482a`.
- Audited main checker SHA-256: `5ed5257f3076936440d2e49b3158690958b0977fde7753f82a6b27e5774e77d3`.
- Audited matrix checker SHA-256: `43c439c0f1f7c1f577bd04eb6337a18498680196e312ac4b16696f044a4ae602`.

`independent_audit_check.py` imports neither existing checker. It exhausts all signed-XOR maps through level three, checks the exponent-four and cycle-count formulas, directly evolves literal guard-coordinate clocks through width 81, verifies uniform local occupation and exact four-lap first returns, and checks nontrivial valuation examples through width 27. These finite controls do not replace the proofs.

Optional scope polish: state n≥1 explicitly when choosing the e1-triggered four-lap family. The general exponent bound itself remains valid at level zero.

Reproduce the independent controls with:

    python3 independent_audit_check.py > independent_audit_results.json
    python3 -O independent_audit_check.py > independent_audit_results.optimized.json
    cmp independent_audit_results.json independent_audit_results.optimized.json
