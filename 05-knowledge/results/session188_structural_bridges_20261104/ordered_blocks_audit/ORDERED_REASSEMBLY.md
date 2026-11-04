# Lossless ordered-block receipt package

Run `python reassemble_ordered.py` in this directory to restore `producer-replay.json`. Each restored file has 148161 bytes. The records are preserved exactly, not replaced by summary hashes.

The reassembler verifies every chunk and restored output against the manifest before writing anything. It refuses to replace a differing existing output. Normal and optimized modes must produce identical files.
