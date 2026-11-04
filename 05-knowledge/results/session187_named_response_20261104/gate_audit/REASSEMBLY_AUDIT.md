# Independent audit of lossless gate-data packaging

**PASS.** This audits the local publication package, not remote GitHub delivery.

All23 named parts part000 through part022 were copied, together with REASSEMBLY.json, REASSEMBLY.md, reassemble.py and the unchanged original verifier dependencies, into a fresh workspace audit directory. No original data file was copied into that directory before reassembly.

The published reassemble.py reconstructed exactly2149227 bytes with SHA256

    09420f040e1294d429ab29d1cbd47e34b12717b695e8a93dafcbfe1d439bbfc7

This equals the original full_response.json byte for byte and the independently pinned hash from the earlier gate audit. The complete39-packet JSON is preserved; no named query or coefficient was dropped.

An independent decode used Python's zlib gzip wrapper rather than the publisher's gzip.decompress. It checked ordered unique chunk names, all23 chunk SHA256 values, complete base64 decoding, gzip EOF, absence of trailing unused data, exact output length, final SHA256 and byte equality with the script's result. Two hostile controls deliberately corrupted only the audit copy: a changed chunk byte and a wrong final manifest hash. Both were rejected before modifying the already reconstructed output; the audit inputs were restored afterward. Normal and optimized independent runs produce identical receipts.

The unchanged original verify.py was then run against the reconstructed JSON normally and under python -O. Both outputs match each other and the original gate receipt byte for byte. They certify39 packets,221208 retained paired-battery count slots,39763 retained nonzero coefficients and714 author literal count controls, with the original proof scope unchanged.

The reassembly script, manifest and instructions are consistent. The script writes only full_response.json in its own directory after checking every hash and the final length. The actual manifest contains only the23 expected local basename paths. No packages were installed and no files outside the permitted workspace were changed.

Receipts: reassembly.receipt.json, reassembly.receipt.optimized.json, reassembled_author.normal.json and reassembled_author.optimized.json. The full independent gate-mathematics review remains REPORT.md; packaging does not add any mathematical reconstruction or collision claim.

To reproduce from a checkout, run python reassemble.py in the gate directory, then python verify.py and python -O verify.py. A separate audit can compare the reconstructed SHA256 to the value above before running either verifier.
