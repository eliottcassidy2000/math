# Restore the complete continuous-cone certificate

Run `python reassemble.py` from any directory, then `sha256sum -c manifest.sha256` in this directory. Both complete 140,629-byte receipts restore exactly, retaining every cubic coefficient, including zeros, at all stated orders and cover coordinates. No coefficient is replaced by a digest or a summary.

The reassembler validates all chunks, output lengths and hashes, and existing outputs before writing; differing existing files are never overwritten. The mathematical source and original receipts are unchanged. The verifier uses only the Python standard library.
