# Prime-factor flow deliverable

- `PROOF.md`: self-contained LTE exponent-box theorem, monotone valuation/LCM pruning, strict multiplicative cut inequalities, two verified support certificates, and exact finite decision reduction.
- `check_primeflow.py`: standard-library reproducibility checks; run `python check_primeflow.py`.
- `RESULTS.json`: generated exact-check counts and machine-readable rejection traces.

All 34,200 LTE comparisons, 58,482 capacity checks, 600 LCM pattern checks, and eight cut-bookkeeping checks passed. The latter use even perfect numbers to check the algebra and do not supply odd perfect number examples.

Known/open boundary: this gives sound fixed-support certificates, not a bound on the support itself. The example is already excluded by existing prime-count bounds. No historical novelty or global odd-perfect-number result is claimed.
