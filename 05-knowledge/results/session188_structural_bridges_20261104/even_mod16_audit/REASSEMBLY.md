# Restore complete even mod16 audit records

Run `python reassemble.py` from any directory, then verify this directory's existing source/data hash manifest where provided. The following original files are restored byte for byte:

- local.normal.json: 2,480,334 complete bytes.
- local.optimized.json: 2,480,334 complete bytes.

All exact coefficients, physical labels, zeros and complete retained query records are preserved. Digests and summaries do not replace data. Every chunk, raw output length/hash and existing output is validated before any writes; differing existing files are never overwritten. Original mathematical sources and receipts are unchanged. The verifier's documented sibling dependencies still apply.
