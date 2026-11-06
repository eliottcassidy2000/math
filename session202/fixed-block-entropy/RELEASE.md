# Reviewed release

The independently reviewed input proof is preserved in `originals/PROOF.md`, SHA-256 `96aac520797bfc6f8aef07ca9959de4cc8d51c87378b711cf92530ca5f542620`. The released proof changes only its status sentence after the separate review returned PASS. Audit files are in `audit/`.

This release contains the fixed-strong-block shared-root entropy theorem, its physical prime-atom benchmark, literal controls and independently compiled small seeds. The separate common-content and target-67 candidates are excluded. No eventual-coverage threshold is claimed.

Reproduce author controls with `python3 verify.py` and `python3 verify_seeds.py`, and their `python3 -O` forms. The two exploratory matrix scripts are not theorem dependencies. See the audit report for the independent reproduction command.
