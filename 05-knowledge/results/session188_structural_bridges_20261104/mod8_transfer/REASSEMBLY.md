# Restore complete modulo-eight records

Run `python reassemble.py` from any directory, then `sha256sum -c manifest.sha256` in this directory. This reconstructs both complete 470,986-byte producer receipts, including all 1,683 profile records, all 570 literal ambient substitutions, and every retained physical port-core profile and count. Nothing is replaced by a summary hash.

The reassembler validates all chunks, all output hashes and sizes, and all existing output contents before writing. It refuses to replace differing existing files. Source proofs and verifiers are unchanged. Replaying the mathematical verifier requires the sibling `covered_response` and `response` packets.
