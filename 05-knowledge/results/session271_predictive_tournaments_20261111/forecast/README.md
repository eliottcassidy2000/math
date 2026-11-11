# Giant-index tournament forecasts

**Scope:** predict ordinary unlabelled tournament counts accurately, without
exact-count computation. Model bounds are proved in [PROOF.md](PROOF.md);
decimal evaluations are high-precision, not interval-certified.

Representative predictions (rounded for readability):

| n | A000568(n) prediction | log10 of model relative-error upper bound |
|---:|---:|---:|
| 300 | 5.122742023739228 × 10^12886 | −172.46398… |
| 1,000 | 7.554225356196997 × 10^147796 | −592.33429… |
| 10,000 | 1.509893164724477 × 10^15014335 | −6007.87304… |
| 1,000,000 | 1.190087511215613 × 10^150509281608 | −602041.26433… |
| 10^9 | 6.012301390415699 × 10^150514989115770076 | −602059963.60096… |
| 10^12 | 2.039399725886988 × 10^150514997820274377090934 | −602059991291.23539… |

The error column concerns the unrounded identity model 2^binom(n,2)/n!, not
rounding of these displayed mantissas. It does not certify exact integer digits.
The full JSONL also covers10^30,10^100 and10^1000, retaining giant exponents.

Methods implemented: logGamma identity model; four-term Stirling with remainder;
finite-support Burnside corrections; defect-layer polynomial coefficients;
adaptive truncation; local-growth-ratio forecast; two precision levels.

Quantitative structural comparison: the k-disjoint-3-cycle sector is enhanced over a naive
Poisson approximation by64^binom(k,2), apart from finite-population depletion.
The inherited fixed defect layers have explicit finite polynomials, but this is a
classification of symmetry corrections, not all recursive tournament types.

Artifacts:

- `predict.py`, `forecasts.jsonl`: giant-index forecasts and separate errors.
- `compare_methods.py`, `method_comparisons.json`: method errors and adaptive
  cutoffs at6 scales and6 requested precision levels.
- `defect_coefficients.json`: complete cycle inventories for k=1,…,6 (defects2k=2,…,12).
- `validate.py`, `validation.json`:159 pair-orbit/deficit controls and independent
  Stirling checks; no exact sequence cache is consulted.
- `PROOF.md`: definitions, derivations, rigorous envelopes and boundaries.

Reproduce from this directory:

    python3 predict.py --output forecasts.jsonl
    python3 validate.py
    python3 compare_methods.py

Requires Python3 and mpmath (tested1.3.0). No network or large-count files needed.

For other code, import `tail_log(n,3)` to obtain the natural logarithm of a
rigorous symbolic upper bound on δ_n=a(n)n!/2^binom(n,2)−1. Set mpmath precision
before calling; the bound's decimal evaluation is not outward rounded.
`forecast(n)` manages its own precision and gives a much sharper evaluated
upper bound after retaining support through12.

Provenance: the defect-polynomial expansion and all-3 leading terms are
inherited from
`01-canon/theorems/THM-4845-a000568-as-a-gaussian-average-gcd-covariance-cycle-index-exponential-and-the-meromorphic-normalized-series-with-poles-at-powers-of-four.md`,
which explicitly credits
`05-knowledge/results/session98_tournament_growth_20261029/limits/ALL_POLES_CORE_EXPANSION.md`.
See PROOF.md for scope and source availability. This deliverable contributes
the practical support-envelope/adaptive prediction implementation and gas
comparison, rather than claiming the inherited transseries as new.
