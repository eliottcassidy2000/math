# Independent review: higher-order bottleneck clocks

19 October 2026

## Verdict

**The two stated asymptotic formulas and their explicit error bounds are correct under the stated hypotheses.** The corrected clock controls the discrete dynamics uniformly in the small parameter; the proof does not identify a discrete orbit with a differential equation. The quartic remainder estimate is uniform for the full interval `0 <= lambda <= 1`.

No mathematical blocker remains. The initial finiteness wording and numerical-certification implementation were hardened during review. Historical priority has not been established and is appropriately not claimed.

## Analytic checks

1. **Clock defect.** With `s(t)=q(t)-1`, the bounds `|s(t)| <= 5/16`, `|R(t)| <= F t^2`, `A^2 <= 4F`, and `F <= 1/32` give exactly the two displayed integral/logarithm estimates. Their sum is `5057/660 < 8`. Summation gives `8 sum F = 8(y_N+a) <= 18a`. Positivity and compactness justify finiteness in the general lemma.

2. **All even degrees.** Put `c=5/32`. The bounds for `r c^(r-1)`, `r(r-1)c^(r-2)`, and `r^2 a^(r-2)` hold at `r=2`. Their ratios upon replacing `r` by `r+2` are at most `2c^2`, `6c^2`, and `4a^2`, respectively, all less than one. Thus the finite-degree software checks are backed by a genuine all-degree argument. The overshoot is at most `2a^r <= a/4`.

3. **Pure-power endpoints.** The endpoint logarithm lies in `[0,(r/2)log(5/4)]`. The two omitted full-line tails total at most `2a^(1-r)/(r-1)`, even though the actual upper endpoint is the first overshoot rather than exactly `a`. These estimates yield the claimed absolute error bound.

4. **Quartic endpoints and expansion.** Here `f(a) <= 2a^4+a^6 <= min(a/4,1/32)`. Since `f` increases on the positive half-line, the overshoot integral is at most `(y_N-a)/f(a) <= 1`. The logarithmic endpoint term is at most `3log(5/4)`. The exact resolvent remainder is nonnegative and at most `lambda^2 <= 1`, since its denominator is at least `y^12`. Its integral is at most `2a`. The two tail bounds are `2/(3a^3)` and `2/a`; their sum with the discrete error gives the stated uniform constant.

5. **Coefficients.** Substitution gives `(2/r)B(1/r,1-1/r) = (2pi/r)csc(pi/r)` for the pure-power integral. For the correction, scaling contributes `epsilon^(-1/4)` and the remaining integral is `(1/2)B(7/4,1/4) = 3pi sqrt(2)/8`. The correction has a minus sign. These identities agree with the [Euler beta integral](https://dlmf.nist.gov/5.12) and [gamma recurrence/reflection formulas](https://dlmf.nist.gov/5.5).

## Independent numerical certificates

`independent_checks.py` uses an exact grid with denominator `2^192`. It evaluates each endpoint polynomial using integer arithmetic and rounds only with exact floor/ceiling division. Monotonicity of the entire step map preserves the enclosure. Every step checks its domain, enclosure order, and the threshold decision; undecidable threshold comparisons fail explicitly.

All ten reported counts agree exactly with the author's verifier:

`85, 387, 888, 8767, 72452, 8761, 72433, 8756, 72415, 59230`.

Four additional checks pass:

- `r=2, epsilon=2^-6, lambda=0`: `N=13`.
- `r=4, epsilon=2^-12, lambda=1`: `N=887`.
- `r=4, epsilon=2^-16, lambda=1/3`: `N=8763`.
- `r=4, epsilon=2^-18, lambda=1`: `N=25376`.

Normal and optimized Python runs produce identical independent receipts. The receipts include exact threshold margins and hashes identifying the reviewed source files. These finite controls supplement the proof; they do not establish its parameter-uniform statements.

The revised author verifier uses directed nonnegative repeated products for even powers. This avoids relying on the qualified C-implementation guarantee for two-argument `Decimal.power` in the [Python documentation](https://docs.python.org/3/library/decimal.html#decimal.Context.power). Its transcendental constant comparisons remain appropriately labeled numerical controls rather than rigorous interval certificates.

## Literature and scope

The limited primary-source check confirms substantial prior context:

- [Kuehn, *Scaling of Saddle-Node Bifurcations: Degeneracies and Rapid Quantitative Changes*, Proposition 2.1](https://arxiv.org/abs/0807.1546) gives the higher-order continuous-time exponent `-(m-1)/m`.
- [Fontich and Sardanyes, *General scaling law in the saddle-node bifurcation: a complex phase space study* (2008)](https://doi.org/10.1088/1751-8113/41/1/015102) studies degenerate analytic one-dimensional flows and their scaling laws. Kuehn explicitly credits their integral calculation.
- [Duarte, Januario, Martins and Sardanyes (2012)](https://repositorio.ipl.pt/entities/publication/2be5889e-c076-405a-afa6-e4185a7b5d70) studies inverse-square-root delay for one-dimensional discrete maps.
- [Brockmoeller, Scherz and Srkalovic, *Pi in the Mandelbrot set everywhere*, Theorem 5](https://arxiv.org/html/2505.07138v1) explicitly includes an `O(1)` escape-time estimate in its satellite setting, alongside its discussion of earlier quadratic results.

These sources warrant the note's conservative attribution. This check does not establish novelty of the corrected-clock argument, exact polynomial families, or uniform constants. The arithmetic conclusion should remain the stated one: retaining only a leading relative asymptotic can lose an unbounded next contribution; it proves no prime or Collatz claim.

## Reproduction

Run `python independent_checks.py` from this directory, with the reviewed note and its `verification.json` in the parent directory. Run again with `python -O independent_checks.py` to confirm that validation does not depend on removable Python assertions. The checked output is in `independent-verification.json`. The parent directory contains the matching normal and optimized author-verifier receipts.
