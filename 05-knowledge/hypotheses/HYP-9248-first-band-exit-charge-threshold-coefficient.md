---
id: HYP-9248
title: "The first half-band exit charge of the Collatz odd map has a finite threshold coefficient: the limsup over odd n>1 of charge/(E(n)^2-E(h)^2), attained along drop-one glide-record exits, is conjectured to be the random-model constant c* = (ln 2/2) max f/KL = 4.6407 and the supremum to lie below 5; coefficient 3 is REFUTED (least witness 56,924,955), coefficient 4 is REFUTED (1,008,932,249,296,231), coefficient 5 is the least integer candidate (max observed ratio 4.278); heuristically every first-exit potential must grow at least like 4.64 E^2"
status: >
  OPEN (the value of the supremum and the sufficiency of coefficient 5). REFUTED, FINITE-EXACT: the
  coefficient-3 hypothesis of three_problems_marked_defects_20261012a.md section 3 (least counterexample
  56,924,955 by exhaustive scan of all odd n < 2^33 with exact big-integer recheck of the 24 orbits that
  left int64; 37 counterexamples below 2^33; nine among the 34 glide-record holders to 2.6e18) and the
  coefficient-4 hypothesis (witness 1,008,932,249,296,231, ratio 415/97). HEURISTIC (random model, exact
  optimization of a two-parameter exponential tilt): the threshold c* = 4.6407 and the per-exit price 9.28 E.
  Independently audited 2026-10-12 (witnesses, scans and the tilt optimum reproduced; the mechanism retyped
  HEURISTIC and the supremum weakened to a limsup; MISTAKE-593).
session: opus-2026-10-12-S24
source: 05-knowledge/results/collatz_band_budget_hostile_20261012.md
scripts:
  - 04-computation/experiments/collatz_band_budget_hostile_20261012.py (+ 05-knowledge/results/collatz_band_budget_hostile_20261012.json)
related:
  - 05-knowledge/results/three_problems_marked_defects_20261012a.md (section 3, the hypothesis (3))
  - 05-knowledge/results/collatz_paid_exit_127_20261012.md (the exit compiler whose quadratic invoice uses coefficient 3)
  - 05-knowledge/results/collatz_pascal_tower_lyapunov_20261005.md (the no-descent rate 3^(h*-1) = 0.9465 per odd step, which is exp(-KL) at the pure tilt)
  - THM-4514 (Beatty structure of no-descent words)
---

# HYP-9248 — the threshold coefficient of the first half-band exit charge

## Setting

For odd `n > 1` iterate `U(x) = (3x+1)/2^a`, `a = v_2(3x+1)`, until the first state `h` with `2(h+1) <= n+1`. The **charge** of this prefix is its number of odd letters `a`; `E(n) = floor(log2(n+1))`. The paid-band hypothesis asked for

    charge <= c [E(n)^2 - E(h)^2]                           (3_c)

for every odd `n > 1`, with `c = 3`; induction then gives `o(n) <= c(E(n)^2-1)` for the total odd charge.

## Statements

1. **(REFUTED, c = 3 and c = 4; FINITE-EXACT.)** The least counterexample to `(3_3)` is `n = 56,924,955` (a glide-record holder): first exit after 198 odd steps at `h = 25,168,801`; `E(n) = 25`, `E(h) = 24`, charge 151 against budget 147. Below `2^33` there are exactly 37 counterexamples; the largest ratio below `2^33` is `176/51 = 3.451` at `n = 96,883,183`. Nine of the 34 glide-record holders to `2.6e18` fail `(3_3)`, always with drop one; `n = 1,008,932,249,296,231` fails `(3_4)` with ratio `415/97 = 4.278`. No tested `n` fails `(3_5)`.
2. **(Conjecture.)** The limsup of `charge/(E(n)^2-E(h)^2)` over odd `n`, attained along drop-one exits, equals `c* = (ln 2/2) max_nu f(nu)/KL(nu||P) = 4.6407`, the maximum over letter laws `nu` on `{1,2,...}` with mean at most `log2 3`, where `P(a) = 2^-a`, `f` is the odd-letter mass and `KL` the relative entropy; the maximizer has `f = 0.7420`, `KL = 0.05541`, letter law `(0.638, 0.222, 0.090, 0.031, ...)`, against the pure no-descent tilt `f = 0.7304`, `KL = 0.05498`, `exp(-KL) = 0.9465 = 3^(h*-1)`; and the supremum over all odd `n` lies below 5. In particular `(3_5)` is conjectured to hold for every odd `n`; the model leaves a margin of only `0.36` in the limit, and fluctuations of a few percent at the scales where records are known.
3. **(Mechanism; HEURISTIC: large-deviation reasoning in the random letter model, not proved even there.)** Among the `N/2` odd starts below `N`, about `N exp(-k KL(nu))` have a prefix of length `k` and empirical law `nu` that stays above `n/2`; the longest such prefix with odd fraction `f` has `k = ln N/KL`, charge `f ln N/KL = (f/KL)(ln 2) E(n)`, that is `9.28 E(n)` at the optimum; the budget of a drop-one exit is `3(2E-1) ~ 6E`. So the first-exit charge and the per-band budget are both linear in `E`, and the coefficient must exceed `9.28/2 = 4.64`. Any potential `Phi(E)` that pays first exits must satisfy `Phi(E) - Phi(E-1) >= 9.28 E`, hence `Phi(E) >= 4.64 E^2`: the quadratic shape of the invoice is forced and so is its coefficient; the total bound it then yields, `o(n) <= 4.64 E(n)^2`, is far above the conjectured linear truth, which is the price of paying at the first exit.

## Evidence for the value of the supremum

Observed `charge/(2E(n)-1)` along the drop-one glide-record failures is `3.08, 3.18, 3.16, 3.44, 3.69, 4.28, 3.81, 3.88, 3.95` at `E = 25, 25, 39, 40, 44, 49, 56, 57, 61` (not monotone; maximum 4.28 at `E = 49`), all below the model's limit, which the longest prefix approaches slowly because it is `ln N/KL` only up to a first-passage correction of order `ln ln N`; the odd fractions of the record prefixes (0.65-0.76) scatter around the tilted law's 0.74. Fluctuations of the charge at fixed length are of relative order `k^(-1/2) ~ 4%`, so ratios up to about `4.8` are not excluded at the scales where glide records are known.

## Cheapest decisive tests

* Extend the exhaustive scan beyond `2^33` with 128-bit arithmetic, or scan the known glide-record list beyond `2^61` (the Oliveira e Silva records): a ratio above 5 refutes statement 2's integer clause.
* Fit the drop-one record ratios against `c* - a/E` to test the approach from below; a plateau below 4.64 at `E` near 100 would argue for a smaller supremum and a sharper letter-law model.
