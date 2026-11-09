# Session231 packing/carry lane

Independent audit: PASS with minor corrections applied ([audit](../ribbon_audit/AUDIT.md)). Chan–Pak completeness remains OPEN. The attempted positive short-block strong-state reset search found no eligible blocks; two precise constructive obstructions and an explicit order kernel survive.

- REPORT.md: arbitrary-cutwidth sparse-arc vertex kernel and exact counting bound. Qualitative sparsity is inherited from THM-4878; explicit kernel and uniform all-cutwidth bound are the scoped extension.
- RIBBON_RESET.md: all-order finite-menu bounded-digit obstruction on arbitrary ribbon prefixes. Strong invariant subsets and unbounded-reach apex operations remain outside its scope.
- SPECTRAL_RESETS.md: exact short-block spectral tests and their all-time implication for repeated growing affine orbits, plus a counterexample to a tempting parity proof.
- verify.py / VERIFICATION.json:2,366 small ordered tournament controls.
- search_resets.py / RESET_SCAN.jsonl: span2 all-letter blocks through length6.
- search_span3.py / SPAN3_SCAN.jsonl: span3 all-letter blocks through length3.
- search_long_resets.py / LONG_SCAN.jsonl: span2 nonempty-letter blocks through length9.
- parity_probe.py / PARITY_PROBE.json: discovery of the odd-constant-term counterexample.

Run scripts with PYTHONDONTWRITEBYTECODE=1. Spectral scripts use MATH_REPO (defaulting to this session's complete checkout) for inherited exact matrices. INPUT_SHA256SUMS pins the inherited source files. The replay uses the inherited, hash-pinned matrix sources in the repository and does not modify them.
