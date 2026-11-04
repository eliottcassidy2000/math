# Restore complete six-block endpoint and physical response records

Run `python reassemble.py` from any directory, then verify this directory's existing source/data hash manifest where provided. The following original files are restored byte for byte:

- endpoint_table.tsv: 1,163,178 complete bytes.
- normal.json: 472,749 complete bytes.
- optimized.json: 472,749 complete bytes.

All exact coefficients, physical labels, zeros and complete retained query records are preserved. Digests and summaries do not replace data. Every chunk, raw output length/hash and existing output is validated before any writes; differing existing files are never overwritten. Original mathematical sources and receipts are unchanged. The verifier's documented sibling dependencies still apply.
