# Alternating-hub positivity independent audit

Run `python audit.py > normal.json` and `python -O audit.py > optimized.json`, then compare. A standard C++17 compiler runs in a temporary directory; no software is installed. The audit imports no producer code.

Mathematical audit allowlist: REVIEW.md, README.md, audit.py, literal.cpp, normal.json, optimized.json, SHA256SUMS. Producer source and final receipts must exist in sibling cyclic_child_disjoint. Exact proof dependencies are pinned in both audit receipts. Separate lossless-packaging review files may be added under their own packaging manifest.
