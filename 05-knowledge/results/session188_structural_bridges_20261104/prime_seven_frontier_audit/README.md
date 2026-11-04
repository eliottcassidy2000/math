# Independent order-seven frontier integration audit

Read REVIEW.md for the exact proof scope and pins. Reproduce with `python audit.py > normal.json` and `python -O audit.py > optimized.json`, then compare. The auditor reads the frozen sibling integration receipt and reviewed dependencies, independently reconstructs original modules/strong trees, recomputes profiles by unordered partitions, and recounts physical graphs by forward path extension. No producer algorithms are imported and nothing is installed.

Exact allowlist: REVIEW.md, README.md, audit.py, forward.cpp, normal.json, optimized.json, SHA256SUMS.
