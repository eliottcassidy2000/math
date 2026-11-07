# Portable audited three-core release

Every published file is at most131072 bytes. The4,000,008-byte recipe certificate is stored only as19 ordered base64 parts of a gzip stream. No large original is included. certificate/MANIFEST.json pins each part, the compressed stream, and the reconstructed original.

From any directory:

    python3 /path/to/this/package/replay.py

This reconstructs the certificate in a temporary directory, verifies its original hash, runs producer verification, compiles and runs the independent audit, and verifies the standalone exception graphs. Python3 and a C++17 compiler are required; no third-party packages are needed. The independent audit includes full finite coverage enumeration and literal small-graph controls. A successful replay matches all published receipts and removes its temporary files.

For a standalone reconstructed certificate:

    python3 /path/to/this/package/reconstruct.py /external/path/RECIPES.bin

The destination must not already exist. Keep reconstructed large files outside this publication package. Existing producer commands in PROOF.md apply to a working directory containing the reconstructed certificate, not directly to the chunk-only package.

The only mathematical-document change from frozen input is the status promotion in audit/STATUS.diff. audit/FROZEN_PROOF.md preserves the reviewed text; audit/FROZEN_SHA256SUMS pins the complete original packet, including the binary certificate represented losslessly here. REVIEW.md and INPUT_FREEZE.md retain original audit provenance and paths as historical context. No tail-completeness claim is introduced.
