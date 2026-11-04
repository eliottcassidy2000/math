# Marked-pair root-family independent audit

Run `python audit.py > normal.json` and `python -O audit.py > optimized.json`, then compare. A standard C++17 compiler is used in a temporary directory; no software is installed. The audit imports no producer code.

Publication allowlist: REVIEW.md, README.md, audit.py, literal.cpp, normal.json, optimized.json, SHA256SUMS. Producer replay requires sibling prime_pair_root_family and arbitrary_third_root_family/verify.py. The mathematical proof dependencies are pinned in both receipts.
