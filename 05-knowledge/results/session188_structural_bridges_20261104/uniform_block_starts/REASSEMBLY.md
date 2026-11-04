# Restore complete uniform-start records

Run `python reassemble.py` from any directory. Both complete 362,393-byte original receipts are restored exactly, including all physical start/query counts, arbitrary-third controls, all sixteen root families and the invariant-blindness example. No exact count is replaced by a digest.

The reassembler validates every chunk, both output lengths and SHA256 hashes, and every existing output before writing. Differing existing files are never overwritten. Mathematical sources and original records are unchanged.
