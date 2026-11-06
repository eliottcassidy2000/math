# An arbitrary-module equal-row switch

Status: **PROVED, independently audited in the stated scope.**

Two explicit physical contexts of an arbitrary tournament module X have the same padding row (B,e,V,d) while their Hamilton and marked-start counts both differ by 2(p1(X)+p2(X)), where pj counts ordered nonempty j-path covers. For nonempty X this gives two distinct endpoint choices at one actual row; empty X has zero gap. The full row and all seven profile coordinates have explicit exact cover-polynomial formulas.

The contexts are an incidence-bundle reversal after the relabeling independently recorded in audit/REVIEW.md. Selectively reversing only some incidences is a separate construction and is not covered by this theorem. Both outputs retain a global sink; no strong-output claim is made.

This is a positive arbitrary-body identity, not a fiber-widening or coverage theorem. Nesting need not preserve the row. The moment-kernel direction in PROOF.md is a necessary arithmetic condition for varying the module, not a realization theorem or a bounded-endpoint menu.

## Evidence and dependencies

- PROOF.md is the corrected, reviewed proof with status-only promotion. ORIGINAL_PROOF.md preserves the earlier prose; EMPTY_SCOPE_CLARIFICATION.diff qualifies the distinct-endpoint statement to nonempty X. REVIEWED_PROOF.md retains the exact corrected audit input; STATUS_PROMOTION.diff records the later status change only.
- module_tensor.py and MODULE_TENSOR.json give the complete fixed-core certificate. The optimized tensor replay is identical. verify_module.cpp and the normal/debug receipts check 152 literal graphs independently of those gap tensors.
- audit/ contains every one of the 19 frozen audit entries and the manifest itself, unchanged. The independent checker directly partitions core vertices into two labeled paths instead of using the producer's auxiliary merged separators. It checks all arrays, the empty boundary, the integer moment kernel, 176 literal context graphs, 56 padding cases and 24 four-core children.
- dependencies/ freezes published proofs for the 3B+2s physical readout and O/Z padding corollaries, together with the rooted amalgamation identity used there. PROVENANCE.json gives exact repository commits, paths and blob identities. These corollary dependencies are not needed for the self-contained module tensor theorem itself. Where the frozen row-normalization proof refers to its dependency folder, the corresponding three proofs are located together here under their displayed names; its FREE_ADDITIVE_DIFFERENCE and ROOTED_ORDERED_AMALGAMATION names are preserved.

## Reproduction

From this directory:

    python3 module_tensor.py > replay.tensors.json
    cmp replay.tensors.json MODULE_TENSOR.json
    g++ -O3 -std=c++17 verify_module.cpp -o /tmp/module_pair
    /tmp/module_pair > replay.json
    cmp replay.json MODULE_NORMAL.json

From audit/, follow REVIEW.md for independent normal/debug replays and run:

    sha256sum -c SHA256SUMS

The finite core enumeration certifies the all-order identity. The literal graph controls provide additional checks, not an extrapolation or a complete spectrum census.
