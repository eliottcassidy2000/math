# Complete equal-profile response data

Run `python reassemble.py` before either the producer verifier or the adjacent audit. This restores `receipt.json`, `receipt.optimized.json` and `full_response.json` byte for byte. The full response table retains all 2,116 named coordinates, including zeros, for each of its three nine-name sources, together with exact integer interpolation values. All cubics and quartics remain present.

All chunk and output lengths and hashes are verified before any write. The reassembler refuses to replace differing existing files. `manifest.sha256` checks the restored original packet; summary receipts do not replace the full records.
