# Lossless ordered-block receipt package

Run `python reassemble_ordered.py` in this directory to restore `ordered.json`, `ordered.optimized.json`. Each restored file has 148161 bytes. The records are preserved exactly, not replaced by summary hashes.

The reassembler verifies every chunk and restored output against the manifest before writing anything. It refuses to replace a differing existing output. Normal and optimized modes must produce identical files.

Before running `verify_ordered.py` from a fresh checkout, run `python ../simultaneous_response/reassemble.py` to restore its published full named banks. Then compare the verifier output to `ordered.json`. The separate six-name modulo-four bank is not part of this ordered-block package.
