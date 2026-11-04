# Independent lossless signature-bank packaging review

**PASS.** The five base64 chunks decode to the exact original131,072-byte `h_vectors.bin`, not a digest-only proxy, selected-coordinate extract or sorted/unlabelled summary. Mathematical sources and literal physical response data were not edited during this review.

The independent checker verifies every chunk's byte count and SHA256, decodes the concatenated base64 with strict validation, decompresses gzip and compares every restored byte with the original bank. The raw SHA256 is182784ab88f16b360ad12371eb03991bbd8d8a8dcfde7c1ec92c9d0e76bb31f9. The gzip stream has115,054 bytes and SHA25658ceaabc8c1349d27a625025bad9ef552d6057cbb099da6cbb3ad2ab6a926ca6.

The schema is checked on all32,768 unsigned little-endian words. Both unused high bits are zero; each of the thirty retained bits agrees with the exact named row-major h_ij convention, using the independently decoded quotient edge code. Thus the restored bank preserves physical coordinate identity and the complete signature for every quotient, including zero coordinates. This bank is not a full response census of the2^21 physical states.

The reassembler validates all chunks, all decoded output sizes and hashes, duplicate/nonlocal names, and every pre-existing output before writing anything. It refuses to replace a different existing file. The review tested isolated reconstruction and idempotent replay, corruption of the first and last chunk, a conflicting existing output, a second-group late conflict, a bad decoded hash, a nonlocal output, a duplicate output and a nonlocal chunk. Every failure left original conflicting data unchanged and created no partial output.

All eighteen normal/optimized scenarios pass, and paired packaging receipts match exactly. Run `python packaging_check.py > packaging.json` and `python -O packaging_check.py > packaging.optimized.json` from this audit directory. The separate `PACKAGING_SHA256SUMS` binds this report, checker and both receipts without changing the mathematical auditor's manifest.
