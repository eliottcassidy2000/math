# Uniform-block starts independent audit

Run `python audit.py > normal.json` and `python -O audit.py > optimized.json`, then compare. A standard C++17 compiler is used for the independent literal controls in a temporary directory. No software is installed; the audit imports no producer code.

Publication allowlist: REVIEW.md, README.md, audit.py, literal.cpp, normal.json, optimized.json, SHA256SUMS. Producer replay requires sibling uniform_block_starts; the mathematical proof dependencies are pinned in both receipts.
