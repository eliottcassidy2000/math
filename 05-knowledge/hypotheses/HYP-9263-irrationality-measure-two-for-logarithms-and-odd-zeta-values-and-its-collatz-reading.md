---
id: HYP-9263
title: "Owner's hypothesis (2026-10-14): the irrationality measure of every odd zeta value zeta(2k+1), k >= 1, equals 2, and so does the irrationality measure of every logarithm -- log alpha for algebraic alpha > 0, alpha != 1, and ratios log alpha / log beta of multiplicatively independent algebraic numbers, in particular theta = log2 3; Collatz reading through THM-4633: mu(theta) = 2 is the strongest possible Diophantine input for cycles and forces a nontrivial cycle to have of order sqrt(m) odd elements for its minimum m, a scale the computed continued fraction already attains"
status: "OPEN (owner's hypothesis, filed 2026-10-14 by mac-mini-2026-10-14-irrmeasure). Literature: Liu-Jiang-Zhang, arXiv:2610.10192 (7-8 Oct 2026, unrefereed) claim mu(log r) = 2 for every positive rational r != 1 and for nonzero rational linear combinations of such logarithms, by a separated-weight determinant argument; this does not cover ratios such as theta = log 3/log 2 and is not verified here. Known unconditional bounds: mu(log 2) <= 3.57455 (Marcovecchio 2009), mu(log 3) <= 5.125 (Salikhov 2007; 5.1163051 Wu-Wang 2014; 5.116201 Bondareva-Luchin-Salikhov), mu(pi) <= 7.1032 (Zeilberger-Zudilin 2020), mu(zeta(3)) <= 5.513891 (Rhin-Viola 2001); zeta(5) is not known to be irrational (Ball-Rivoal: infinitely many odd zeta values are irrational; Zudilin: one of zeta(5), zeta(7), zeta(9), zeta(11) is), so for zeta(5) the hypothesis presupposes irrationality. Heuristic support: Khinchin's theorem (almost every real has mu = 2) and the standard expectation that 'natural' constants are not exceptionally well approximable; the continued fraction of theta through 20,000 partial quotients is Khinchin-typical: geometric mean 2.6724 (Khinchin 2.6855), Levy constant 1.1819 (pi^2/(12 ln 2) = 1.1866), largest partial quotient 24,477 at index 18,287, and log a_(n+1)/log q_n below 0.001 for n > 1000 (log23_khinchin.py). Concurrent: opus-2026-10-14-S31's HYP-9262 (transfer-function side of the same prompt) also types the Liu-Jiang-Zhang preprint as not reaching log 3/log 2."
session: mac-mini-2026-10-14-irrmeasure
related:
  - THM-4633 (irrationality-measure parametrisation of cycle bounds; the convergent ledger)
  - THM-4634 (parity-number dictionary)
  - J. Liu, K. Jiang, P. Zhang, Irrationality exponents of logarithms of positive rational numbers, arXiv:2610.10192 (2026)
  - R. Marcovecchio, The Rhin-Viola method for log 2, Acta Arith. 139 (2009); V. Kh. Salikhov, On the irrationality measure of log 3, Dokl. Math. 76 (2007); G. Rhin, C. Viola, The group structure for zeta(3), Acta Arith. 97 (2001); D. Zeilberger, W. Zudilin, The irrationality measure of pi is at most 7.103205334137..., Moscow J. Comb. Number Theory 9 (2020)
---

# HYP-9263 — irrationality measure 2 for logarithms and odd zeta values

## Statement

(a) `μ(ζ(2k+1)) = 2` for every `k ≥ 1`. (b) `μ(log α) = 2` for every algebraic `α > 0`, `α ≠ 1`, and
`μ(log α / log β) = 2` for multiplicatively independent algebraic `α, β > 0`; in particular `μ(log_2 3) = 2`.

## Where it stands

* `μ = 2` is the generic value (Khinchin) and the value for all algebraic irrationals (Roth); no
  "natural" transcendental constant is known to have `μ > 2`, and none of the listed ones is known to
  have `μ = 2` unconditionally. Liu–Jiang–Zhang (October 2026) claim (b) for logarithms of positive
  rationals and their rational linear combinations; if correct this covers `log 2` and `log 3` but not
  the ratio `θ = log 3/log 2`, whose Diophantine exponent is governed by linear forms `|q ln 3 - p ln 2|`.
* For `θ` the known effective information is a linear-form exponent (Baker/Laurent–Mignotte–Nesterenko
  type, or Padé-type three-term bounds); the exponent `μ(θ) = 2` would say `|q ln 3 - p ln 2| > q^(-1-ε)`.

## Collatz reading (THM-4633)

* A Syracuse cycle with minimum `m` gives `0 < K/k - θ < 1/(3 m ln 2)`. Under `μ(θ) = 2` with an effective
  constant `c`, every cycle with minimum `≥ M_0` has `k > (3 c M_0 ln 2)^(1/2)`: of order `sqrt(M_0)` odd
  elements. The computed ledger realises exactly this scale at `M_0 = 2^71` (first surviving
  denominator `7.2·10^10 ≈ 1.5·sqrt(2^71)`, the semiconvergent of THM-4633(5); Hercher's bound
  `1.375·10^11 ≈ 2.8·sqrt(2^71)`), because the partial quotients of `θ` are small so far. So the
  hypothesis is *consistent with* and *numerically visible in* the ledger, but adds no effective
  exclusion beyond it, and it says nothing about divergent trajectories (irrational parity numbers).
* (a) has no bearing on Collatz; it is recorded because the owner states both parts together.

## Evidence (FINITE-EXACT, 2026-10-14)

`log23_khinchin.py`: 20,000 partial quotients of `θ` (120+24,000 digits): geometric mean `2.6724` against Khinchin's
`2.6855`, `log q_n / n = 1.1819` against Lévy's `1.1866`, largest partial quotient `24,477` (index 18,287), ten largest
all below `25,000`, and `log a_(n+1)/log q_n < 5·10^(-3)` for `n > 1000` (`311` indices in `(1000, 20000)` exceed
`10^(-3)`, the largest `0.00453` at `n = 1179`, the last at `n = 5113`, so `< 10^(-3)` for `5113 < n < 20000`;
audit E, `audit_E_khinchin.py`; an earlier "`< 10^(-3)` beyond `n = 1000`" here was false, MISTAKE-605). So `θ` is
Khinchin-typical in this range;
`μ(θ) = 2` is the generic expectation it conforms to.

## Decisive tests

* For (b) with `θ`: extend the continued fraction of `θ` to millions of partial quotients and compare the
  growth of `log a_(n+1)/log q_n` with Khinchin's law; a single large partial quotient proves nothing,
  but `a_(n+1) ≥ q_n^ε` along a subsequence would refute `μ = 2`.
* For the Liu–Jiang–Zhang claim: independent verification of the interpolation theorem for
  logarithmic jets (their geometric core); the repo should not cite `μ(log 2) = 2` as known.
