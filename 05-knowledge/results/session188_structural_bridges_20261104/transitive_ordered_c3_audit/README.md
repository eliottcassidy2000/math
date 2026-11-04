# Transitive ordered-C3 independent audit

Frozen publication allowlist: `REVIEW.md`, `README.md`, `audit.py`, `weighted.cpp`, `literal.cpp`, `normal.json`, `optimized.json`, `SHA256SUMS`.

Run `python audit.py > normal.json`, `python -O audit.py > optimized.json`, then compare. The script reads the frozen producer receipts, reconstructs the 23-term polynomial independently by exterior run words, and compiles its own two predecessor count engines into a temporary directory. It imports no producer implementation.

PASS: 598 weighted actual-query counts, 144 full-name counts, nine exact quartic comparisons, and all 32 whole-child family controls. The review separately checks the all-order Stirling-ratio proof and the broader original-module restoration argument. No arbitrary strong ambient lift is inferred.
