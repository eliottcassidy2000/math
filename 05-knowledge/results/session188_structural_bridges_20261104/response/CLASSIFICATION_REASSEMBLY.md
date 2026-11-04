# Restore the exact order-five classification

The two author outputs and two independent replay copies are byte-identical, each159,614 bytes. Their complete100-coordinate vectors and all1,024 named tournament records are stored losslessly in `classification.json.gz.b64`.

Run `python3 reassemble_classification.py` from this directory. It restores `classification.json`, `classification.optimized.json`, and the review's two producer-replay copies in `../response_review/`, verifying every size and SHA-256 before writing. It refuses to replace different existing bytes. Run `sha256sum -c SHA256SUMS` from `../response_review/` after reassembly.

Then reproduce the result with `python3 classify_five.py` and `python3 -O classify_five.py`, or run `python3 independent.py` and `python3 -O independent.py` in `../response_review/`. The original mathematical verifiers and review are unchanged. This is storage packaging, not a change in the finite universe or scope.
