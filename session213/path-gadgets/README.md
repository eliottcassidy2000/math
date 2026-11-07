# Binary-comb path barrier

The universal Hamiltonian-cycle binary comb has a density-zero set of **total** Hamiltonian-path counts: at most O(X^{log_3 2}) distinct values up to X. This does not rule out the wider two-apex language and does not settle Chan-Pak.

`PROOF.md` gives an independent odd-cycle-collection proof and an exact same-count/different-child hostile. The concurrent context-transform packet proves a stronger bound using explicit Hamiltonian paths.

Reproduce with `python3 verify.py` or `python3 -O verify.py`; output must match `RESULT.json`. All 2,040 combinations through 8 coding bits pass, including 120 independent odd-cycle-formula checks, 56 complete-permutation checks, and 22,138 explicit two-cycle certificates.

Publication inventory: `README.md`, `PROOF.md`, `verify.py`, `RESULT.json`, `SOURCES.json`, `SHA256SUMS`. The exploratory `probe.py`, byte-identical optimized receipt and bytecode are not release files.
