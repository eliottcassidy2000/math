# Session208: Balanced order14 witness

An explicit strong tournament on14 vertices has H=11390625=15^6. Dominance joins give the improved global construction rate upper bound7/3. This is neither a global optimality theorem nor a minimum-order classification.

Run `python3 verify.py` and `python3 independent_check.py` with assertions enabled. Output should match `verify.out` and `independent.out`; RESULT.json is deterministically regenerated. PROOF.md gives the complete code convention, inference, limits and optional discovery replay. Verification is fast; the capped discovery replay takes longer.

`DEPENDENCIES.json` pins the preceding oracle publication. `REVIEW.md` records independent audit. Verify frozen payloads using `sha256sum -c SHA256SUMS`.
