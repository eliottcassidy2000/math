# An infinite four-choice endpoint fiber

Status: **PROVED, independently audited in the stated scope.**

For every k>=0, four explicit tournaments on 12+2k vertices share one padding row (B,e,V,d). Their marked-start values differ from the first by {0,8,14,22} times 2^k. Their Hamilton counts have exactly the same offsets; L,R,e,d and beta remain equal.

The construction selectively reverses the incidences between one retained core vertex and a common transitive spine. The four middle words are 01010,01110,10001,10101, each wrapped by k leading ones and k trailing zeros. The all-order proof diagonalizes the exact source/sink actions on colored path-cover polynomials and checks every fixed-core eigencomponent. It is not an extrapolation from sample scales.

There are four choices at each scale, not a growing menu. The starts are large and the offsets do not supply adjacent odd targets through the ordinary readout. This establishes no interval coverage, F(8,137), or eventual odd-spectrum threshold. Every output retains a global sink; strong outputs are not claimed.

## Proof and evidence

- PROOF.md is the audited proof with status-only promotion. REVIEWED_PROOF.md preserves the frozen original; STATUS_PROMOTION.diff is the complete change.
- derive.py and TENSORS.json contain the exact colored and diagonal-coordinate tensors and all seed eigencomponents. The optimized replay is byte-identical.
- check.cpp, NORMAL.txt and DEBUG.txt independently control the formulas on all four outputs at k=0,1,2,3, through order 18. The all-order proof is the finite operator certificate, not those 16 graph counts.
- audit/ retains all 17 entries of the independent manifest and SHA256SUMS unchanged. The audit independently verifies all 70 monomials through degree 4, every profile/eigencomponent, 12,600 source/sink controls including nontransitive bodies, and all 16 literal graphs through order 18. The zero-eigenvalue term is present only at k=0.
- dependencies/ freezes the published physical 3B+2s readout used only in the stated corollary, with exact immutable provenance. The main four-choice theorem is self-contained and does not depend on that readout.

The separate exploratory selective-word scans and failed attempts to increase the choice count are not part of this reviewed release. No iteration beyond the displayed wrapping law is asserted.

## Reproduction

From this directory:

    python3 derive.py > replay.tensors.json
    cmp replay.tensors.json TENSORS.json
    g++ -O3 -std=c++17 check.cpp -o /tmp/four_choice
    /tmp/four_choice > replay.txt
    cmp replay.txt NORMAL.txt

From audit/:

    sha256sum -c SHA256SUMS
    python3 check_algebra.py > algebra.replay.json
    cmp algebra.replay.json ALGEBRA_NORMAL.json
    g++ -O3 -std=c++17 check_graphs.cpp -o /tmp/four_choice_independent
    /tmp/four_choice_independent > graphs.replay.json
    cmp graphs.replay.json GRAPHS_NORMAL.json

See audit/REVIEW.md for the complete optimized/debug commands. Checks use explicit failures and remain active under optimization.
