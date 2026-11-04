# Independent lossless-package review: five_profile_types

PASS. Every restored output was independently decoded and compared byte for byte with its untouched original. Both normal and optimized checker modes pass with identical receipts.

The package contains 1 compressed groups in 4 chunks, restoring 2 complete files and 3,016,756 total raw bytes.
- normal.json: 1,508,378 bytes; SHA-256 7b9cba1e667b794a050adb5ab412f52e5ad12880de280e27779f2faa660c6e90.
- optimized.json: 1,508,378 bytes; SHA-256 7b9cba1e667b794a050adb5ab412f52e5ad12880de280e27779f2faa660c6e90.

Existing mathematical manifests and the final source/dependency pins in the complete original receipts were checked against the actual files. Mathematical sources, data, prior audits and existing manifests were not changed. The exact verified paths and hashes are retained in the paired lossless receipts.

The reassembler validates every chunk, decoded length/hash, unique local output path and existing output across all groups before its sole output-write loop. This was reviewed and independently exercised in both Python modes. Isolated tests cover fresh restoration, idempotence, first-chunk and late-group corruption, a conflicting final existing output, a late conflict after an already-valid first output, malformed final raw hash and size, duplicate output names, nonlocal outputs and nonlocal chunk names. Every failure left absent outputs absent; the valid early file was not rewritten and the conflicting final file was untouched.

The isolated restoration tests copy only this package’s metadata, reassembler and declared chunks and run from an unrelated working directory. They require no mathematical code imports. To rerun the complete checker in a repository checkout, restore the complete records first with `python reassemble.py`, then run `python lossless_check.py` and `python -O lossless_check.py`; original source/dependency files named by the mathematical receipts must also be available for their integrity checks.

This is a packaging review, not a replacement mathematical proof or audit. It preserves all original bytes, including exact coefficients, physical names, zero entries and complete query records. It does not replace those records by hashes or summaries.

Exact incremental five-file allowlist: LOSSLESS_REVIEW.md, lossless_check.py, lossless.normal.json, lossless.optimized.json, lossless.sha256.
