# Positive approximate construction of actual source slopes

Unreviewed proof candidate. Exact target occupancy remains OPEN.

- `PROOF.md`: actual-index bit reversal, positive Z-run errors, explicit slope
  approximation and discrepancy, elementary joint-grid bounds with rounding,
  and the exact error transferred to the J form.
- `verify.py`: portable exact rational controls; no external dependencies.
- `SOURCES.md`: standard Hammersley naming only; no sharp literature estimate
  is imported.
- `VALIDATION-normal.json`, `VALIDATION-optimized.json`: matching PASS reports.

The construction chooses A then u=2^(k+1)-rev_k(A), always u>=2. It guarantees
an approximate slope, not an exact prescribed bilinear count.
