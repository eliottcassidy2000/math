# Independent ambient lossless-package review

PASS. The separately named ambient package restores only ambient.normal.json and ambient.optimized.json. Each complete 280,317-byte output was independently decoded and compared byte for byte with its untouched original. The paired normal and optimized audit receipts agree exactly.

The AMBIENT_SHA256SUMS allowlist and every mathematical source pin in both originals were verified against their actual files. The source pins resolve to the sibling even_cyclic_mod16_ambient packet. The local package, its outputs, local review and local manifests are outside this review and remain unchanged.

- ambient.normal.json: 280,317 bytes; SHA-256 0a21a17128b096012d1eb1269ff82f7d23913d09045868e3c69ed05ba6d3a305.
- ambient.optimized.json: 280,317 bytes; SHA-256 0a21a17128b096012d1eb1269ff82f7d23913d09045868e3c69ed05ba6d3a305.

The reassembler validates all compressed chunks, restored lengths and hashes, unique local output paths and pre-existing outputs before its only write loop. Independent isolated tests in normal and optimized Python cover fresh restoration, idempotence, first/late chunk corruption, malformed final raw hash/size, duplicate/nonlocal output names, a nonlocal chunk name, and a conflicting last output. A separate late-conflict test retains an already-correct first output byte for byte without rewriting it. Every failure leaves absent outputs absent and conflicts untouched.

The checker is self-contained and uses paths relative to its own repository directory. Isolated controls copy only AMBIENT_REASSEMBLY.md, AMBIENT_REASSEMBLY.json, ambient_reassemble.py and the declared chunks, then run from an unrelated working directory. To reproduce: run python ambient_reassemble.py, then python ambient_lossless_check.py and python -O ambient_lossless_check.py. The mathematical source/dependency files pinned by the originals must also be available.

This is only a lossless packaging review. It does not revise mathematical claims or replace complete physical-name, coefficient or query data with summaries.

Exact incremental allowlist: AMBIENT_LOSSLESS_REVIEW.md, ambient_lossless_check.py, ambient_lossless.normal.json, ambient_lossless.optimized.json, ambient_lossless.sha256.
