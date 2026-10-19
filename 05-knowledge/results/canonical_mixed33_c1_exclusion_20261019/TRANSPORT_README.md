# Lossless certificate transport

The large gzip certificate `canonical331-stage/certificates.jsonl.gz` is stored in this repository as ordered raw-byte parts, not as one gzip blob. No data has been omitted, recompressed, or changed. The original mathematical proof, verifier sources, receipts, and frozen manifest remain unchanged.

Before running the original audit or replay commands, run from this directory:

```sh
python3 reconstruct_transport.py
```

The script verifies every part and the complete SHA256 before creating the original path. It accepts an already matching original and refuses to overwrite a mismatching existing file. `python3 reconstruct_transport.py --check-only` checks the parts without writing the original. The original verifier commands then work unchanged, including their SHA256 and file-size checks.

- Original path: `canonical331-stage/certificates.jsonl.gz`
- Original size: 14372221 bytes
- Original SHA256: `3c65b259d662e4d3490cd745c11239bdaad9513d50872be1199dd97ad32fdafb`
- Parts: 2, each at most 8 MiB

[Transport manifest](TRANSPORT_MANIFEST.json) records order, size, and SHA256 of every part. [Reconstruction script](reconstruct_transport.py) uses only Python's standard library. This is a transport adaptation, not a new mathematical certificate.
