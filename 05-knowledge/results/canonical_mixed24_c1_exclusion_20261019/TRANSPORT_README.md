# Lossless certificate transport

The large gzip certificate `certificates.jsonl.gz` is stored in this repository as ordered raw-byte parts, not as one gzip blob. No data has been omitted, recompressed, or changed. The original mathematical proof, verifier sources, receipts, and frozen manifest remain unchanged.

Before running the original audit or replay commands, run from this directory:

```sh
python3 reconstruct_transport.py
```

The script verifies every part and the complete SHA256 before creating the original path. It accepts an already matching original and refuses to overwrite a mismatching existing file. `python3 reconstruct_transport.py --check-only` checks the parts without writing the original. The original verifier commands then work unchanged, including their SHA256 and file-size checks.

- Original path: `certificates.jsonl.gz`
- Original size: 18660613 bytes
- Original SHA256: `84b5d2172b7e2f3e1282db98d41452046b34e9be17578041e2b452a20c6265ba`
- Parts: 3, each at most 8 MiB

[Transport manifest](TRANSPORT_MANIFEST.json) records order, size, and SHA256 of every part. [Reconstruction script](reconstruct_transport.py) uses only Python's standard library. This is a transport adaptation, not a new mathematical certificate.
