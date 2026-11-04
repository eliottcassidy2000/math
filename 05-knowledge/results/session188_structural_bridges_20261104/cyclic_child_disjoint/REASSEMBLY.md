# Restore complete cyclic-child quartic records

Run `python reassemble.py` from any directory, then `sha256sum -c SHA256SUMS` in this directory. Both complete 201,379-byte original receipts are restored exactly, including every retained physical query, finite-difference value, kernel control, and root-containing family. Summary counts and digests do not replace exact data.

The reassembler validates every chunk, output length and hash, and every existing output before writing. Differing existing files are never overwritten. Mathematical sources and original receipts are unchanged. Replaying verify.py requires the adjacent arbitrary_third_root_family helper and the other proof dependencies pinned in the receipt.
