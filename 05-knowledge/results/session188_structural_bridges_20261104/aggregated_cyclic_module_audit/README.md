# Independent aggregate audit

Run `python audit.py` or `python -O audit.py`. The checker reads the frozen aggregated_cyclic_module packet, imports no producer code, and compiles its independent predecessor-DP counter using the existing compiler. No installations are needed.

Exact seven-file allowlist: REVIEW.md, README.md, audit.py, literal.cpp, normal.json, optimized.json, SHA256SUMS.
