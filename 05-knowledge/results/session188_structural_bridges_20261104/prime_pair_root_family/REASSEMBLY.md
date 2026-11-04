# Restore complete marked-pair root-family records

Run `python reassemble.py` from any directory, then `sha256sum -c SHA256SUMS` in this directory. Both complete 272,915-byte original receipts are restored exactly, including all 260 root-containing families, physical XOR targets, complete selected query values, branch data, and fixed batteries. No exact value is replaced by a digest.

The reassembler validates every chunk, both output lengths and hashes, and every existing output before writing. Differing existing outputs are never overwritten. The mathematical sources and original records are unchanged. Replaying verify.py requires the sibling arbitrary_third_root_family helper and the mathematical dependencies pinned in its receipt.
