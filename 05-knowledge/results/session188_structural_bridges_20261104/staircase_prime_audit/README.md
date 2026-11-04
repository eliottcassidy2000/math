# Staircase-prime independent audit

Run `python audit.py > normal.json` and `python -O audit.py > optimized.json`, then compare the outputs. The producer and its runtime dependency must be sibling directories named `staircase_prime_quotient` and `marked_pair_quotient`; the audit itself imports neither.

Publication allowlist: REVIEW.md, README.md, audit.py, normal.json, optimized.json, PACKAGING_REVIEW.md, packaging_audit.py, packaging.normal.json, packaging.optimized.json, SHA256SUMS. Producer source and dependency pins are in both receipts. No downloaded source images or private scratch files are needed.

The packaging audit is independent and checks lossless restoration and validation-before-write negative controls. Run `python packaging_audit.py` and `python -O packaging_audit.py`; complete original producer receipts must be restored first.
