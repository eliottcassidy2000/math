# Restore complete ordered cyclic module audit records

Run `python reassemble.py` from any directory, then verify this directory's existing source/data hash manifest where provided. The following original files are restored byte for byte:

- normal.json: 907,980 complete bytes.
- optimized.json: 907,980 complete bytes.
- independent.normal.json: 794,943 complete bytes.
- independent.optimized.json: 794,943 complete bytes.

All exact coefficients, physical labels, zeros and complete retained query records are preserved. Digests and summaries do not replace data. Every chunk, raw output length/hash and existing output is validated before any writes; differing existing files are never overwritten. Original mathematical sources and receipts are unchanged. The verifier's documented sibling dependencies still apply.
