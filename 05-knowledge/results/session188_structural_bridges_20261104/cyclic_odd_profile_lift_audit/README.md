# Odd-core cyclic lift independent audit

Frozen publication allowlist: `REVIEW.md`, `README.md`, `audit.py`, `count.cpp`, `normal.json`, `optimized.json`, `SHA256SUMS`.

Run `python audit.py > normal.json`, `python -O audit.py > optimized.json`, and compare. Only the standard Python library and the available C++17 compiler are used. The audit imports no producer implementation and independently recomputes every retained physical Hamilton and two-cover count.

PASS: full32×64 canonical table, complete profile/mask transfer, eleven larger-child/ambient cases, and64 root-family/context records including16 strict increases. The proof review separately checks all-cover homogeneity, odd-profile cancellation, the physical outer lift and common negative dual query. No exact collision or unrestricted prime-quotient theorem is inferred.
