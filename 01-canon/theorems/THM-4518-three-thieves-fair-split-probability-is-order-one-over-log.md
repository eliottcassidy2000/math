---
id: THM-4518
title: "Three thieves and the planar walk: for a uniformly random binary word of shape (3m, 3x) with x/m -> rho in (0,1), the probability that some cut position gives three consecutive fair arcs (each with x ones) is Theta(1/log m): E[N_3] = m C(m,x)^3/C(3m,3x) -> sqrt(3)/(2 pi rho(1-rho)), the exact second moment E[N_3^2] = (m/C(3m,3x)) sum_d sum_t C(d,t)^3 C(m-d,x-t)^3 ~ (log m)/(pi sqrt3 rho(1-rho)), the lower bound is E[N]^2/E[N^2], and the upper bound follows from the return count ~ c log m of the planar imbalance walk after its first visit; with j = 2 always (IVT) and j >= 4 transient (E[N_j] -> 0 like m^(-(j-3)/2)), the one-arc-per-thief necklace problem is Polya's trichotomy"
status: >
  PROVED (Theorem 1 of the note: two-sided order 1/log m; the local-limit lemma
  for hypergeometric prefix counts is derived from the binomial local limit
  theorem with a uniform remainder, CITED as standard; Hoeffding for sampling
  without replacement CITED) + FINITE-EXACT / EMPIRICAL for the constants
  (P(N_3 >= 1) log m = 2.13..2.25 at m = 10^3..10^5; E[N_3^2]/log m -> 0.735 at
  rho = 1/2, computed exactly to m = 10^4 as 1.05, 1.01, 0.98, 0.95, 0.93 from m
  = 10^2). NOT independently audited. Closes the open item of THM-4515 section 5.
  Setting: w uniform of shape (K,X) = (3m,3x); a cut at r in Z/m is fair if the
  three windows [r+im, r+(i+1)m) each hold x ones; N = # fair cuts in [0,m).
  Mechanism of the upper bound: with W(r) = (a_2 - a_1, a_3 - a_2)(r) the
  difference walk of the three windows' prefix counts (independent
  hypergeometric processes given the window counts), a fair cut is a visit of
  W to the point z = (x - x_1, x - x_2); after the first visit the suffix
  processes are fresh (Lemma 1.1), and three independent hypergeometric prefix
  counts coincide at time u with probability >= c_3/u on u <= m/(C_1 log m)
  (Lemma 2.3, from the local limit theorem and Hoeffding's good event), so
  E[N | first visit at u < m/2] >= (c_3/2) log m; since E[N] <= 1.2 this forces
  P(first visit < m/2) <= C/log m, and rotation by m/2 doubles it to P(N >= 1).
  Lower bound: Cauchy-Schwarz with the exact second moment (two fair cuts at 0
  and d force the six window counts t, x-t, t, x-t, t, x-t).
  NOT claimed: the constant c in P(N_3 >= 1) ~ c/log m (about 2.2 at rho =
  1/2); anything about orbit-generated (non-uniform) words.
source: collatz-necklace-20260929 session (mac-mini), part 2, 2026-09-29; owner's request "prove the O(1/log m) upper bound for three thieves"
depends_on:
  - 01-canon/theorems/THM-4515-fair-consecutive-splits-of-cycle-necklaces-are-circulant.md (the exact expectation and the trichotomy statement)
related:
  - 05-knowledge/results/collatz_thieves_20260929_three_thieves_log_law.md (the full proof)
  - 05-knowledge/results/collatz_necklace_20260929_fair_splits_power_clocks_basins.md (section 1.8)
note: 05-knowledge/results/collatz_thieves_20260929_three_thieves_log_law.md
scripts:
  - 04-computation/experiments/collatz_necklace_20260929_fairsplit_asymptotics.py
  - 04-computation/experiments/collatz_necklace_20260929_fairsplit_loglaw.py
outputs:
  - 05-knowledge/results/collatz_necklace_20260929_fairsplit_asymptotics.out
  - 05-knowledge/results/collatz_necklace_20260929_fairsplit_loglaw.out
hash_basis: raw LF bytes (hashes in THM-4515)
audit: NOT independently audited; the proof is elementary given the two cited classical inputs, and its quantitative content (exact expectation, exact second moment) is scripted.
---

# THM-4518 -- three thieves: the fair three-way split of a random necklace has probability Θ(1/log m)

**PROVED (order), EMPIRICAL (constant); not independently audited.** Full proof: [collatz_thieves_20260929_three_thieves_log_law](../../05-knowledge/results/collatz_thieves_20260929_three_thieves_log_law.md).

## Statement

Fix `rho in (0,1)`, let `x/m -> rho`, and let `w` be a uniformly random word of shape `(3m, 3x)`. Let `N` be the number of cut positions `r in [0,m)` at which the three consecutive windows of length `m` each contain `x` ones. Then for `m >= m_0(rho)`

`c_1/log m <= P(N >= 1) <= c_2/log m`,

while `E[N] = m binom(m,x)^3/binom(3m,3x) -> sqrt 3/(2 pi rho(1-rho))` and `E[N^2] ~ (log m)/(pi sqrt 3 rho(1-rho))`.

## Proof in one paragraph

Given the window counts `(x_1,x_2,x_3)` the three prefix-count processes `a_i` are independent hypergeometric processes, and the cut at `r` is fair iff `(a_2 - a_1, a_3 - a_2)(r) = (x - x_1, x - x_2)`: a visit of a planar lattice walk to a fixed point. *Upper bound.* Let `tau` be the first fair cut. Conditionally on `tau = u` (an event of the prefix sigma-field) the suffixes of the three windows are fresh uniform arrangements, so the cuts `u + u'` are fair exactly when the three suffix prefix-counts coincide at time `u'`. On Hoeffding's good event (all window counts and all prefix counts within `sqrt(m log m)` of their means) the suffix parameters are balanced enough that the three hypergeometric laws at time `u'` overlap on a window of width `~ sqrt(u')` where each has mass `>= c_0/sqrt(u')` (binomial local limit theorem applied to the conditional representation of the hypergeometric law), giving coincidence probability `>= c_3/u'` for `u_0 <= u' <= m/(C_1 log m)`; summing, `E[N | tau = u] >= (c_3/2) log m` for `u < m/2`. As `E[N] <= 1.2`, `P(tau < m/2) <= 2.4/(c_3 log m) + 18/m`, and rotating the word by `m/2` shows `P(N >= 1) <= 2 P(tau < m/2)`. *Lower bound.* `P(N >= 1) >= E[N]^2/E[N^2]` with the exact second moment `E[N^2] = (m/binom(3m,3x)) sum_{d<m} sum_t binom(d,t)^3 binom(m-d,x-t)^3` (two fair cuts at `0` and `d` force the six window counts `t, x-t, t, x-t, t, x-t`), which the local limit theorem evaluates as `(1+o(1)) (m/(2 pi sqrt3 rho(1-rho))) sum_d 1/(d(m-d)) ~ (log m)/(pi sqrt 3 rho(1-rho))`. ∎

## Numbers and the trichotomy

Monte Carlo (`rho = 1/2`): `P(N >= 1) log m = 2.13, 2.17, 2.22, 2.19, 2.25` at `m = 10^3, 3 10^3, 10^4, 3 10^4, 10^5`; `E[N | N >= 1] = 3.5 .. 5.4 ~ 0.47 log m`; exact `E[N^2]/log m = 1.05, 1.01, 0.98, 0.95, 0.93` at `m = 10^2 .. 10^4`. The size-biased mean `E[N^2]/E[N] = E[N | fair at 0] ~ 0.735 log m` exceeds the ordinary conditional mean, which is why the second moment gives only the lower bound and the return-count argument is needed for the upper bound. With THM-4515: two thieves always (a one-dimensional bridge has a forced zero), three thieves with probability `~ 2.2/log m` (planar, recurrent but marginal), four or more essentially never (transient, `E[N_j] ~ sqrt j m^{(3-j)/2}(2 pi rho(1-rho))^{-(j-1)/2}`). Collatz reading: a hypothetical `3x+1` cycle with `gcd(K,X) = 3` admits the cyclotomic (Eisenstein-type) factorisation of its clock through a fair 3-split with probability about `2.2/log(K/3)` under the uniform-necklace model.
