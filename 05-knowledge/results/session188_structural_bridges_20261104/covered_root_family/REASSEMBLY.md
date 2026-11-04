# Complete root-family records

Run `python reassemble.py` to restore both complete verifier receipts. Every original family, physical target row, selected query, value, and retained endpoint battery is preserved byte for byte. The original `manifest.sha256` checks the restored files.

All chunks and outputs are length- and hash-checked before any write. Differing existing outputs are not replaced. The producer imports the already published sibling `covered_response/verify.py`; the new independent finite audit imports no producer code.
