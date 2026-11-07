# Portable release review

PASS. The only producer proof edit is its independent-audit status promotion, preserved as STATUS.diff against FROZEN_PROOF.md. README.md additionally links directly to REPLAY.md and explains that the original binary is absent until reconstruction; this packaging-only change is included in STATUS.diff and FROZEN_README.md preserves its original. The original binary certificate is represented losslessly by nineteen ordered base64 chunks of one gzip stream, with per-part, compressed-stream and original-file SHA-256 verification. The reconstructed bytes match the exact audited original certificate.

Full replay was run from /tmp using replay.py. Producer verification, the independent full construction scan, all independent literal controls, and standalone exception counts match the retained receipts exactly. REPLAY.json records this successful run.

Every file listed in PUBLICATION_FILES.txt is at most131072 bytes. No original large certificate, compiled executable or cache file is included. The largest chunk is120000 bytes. Reconstruction occurs in a temporary or explicitly chosen external path. All code resolves package inputs relative to __file__ or takes explicit command-line paths.

The historical input manifest remains FROZEN_SHA256SUMS; the current publication manifest is SHA256SUMS. No Git or stage edits were performed.
