# Restore full staircase-quotient records

Run `python reassemble.py` from any directory, then `sha256sum -c SHA256SUMS` in this directory. Both complete 103,600-byte original receipts are restored exactly, including all retained positive and zero-boundary quartic batteries and every actual query value. No cubic/quartic effect is omitted, and no full record is replaced by a digest.

The reassembler validates all chunks, output hashes and sizes, and existing output contents before writing, and refuses conflicting files. The mathematical source and original receipts are unchanged. The verifier requires the sibling marked_pair_quotient helper and the proof dependencies named in its receipt.
