---
id: THM-4495
title: "The no-descent residue count W_k = |Bad_k| of Collatz has exact order 2^(hk) k^(-3/2): a self-contained Spitzer identity k W_k = sum_n B_n W_(k-n) (binomial tails B_n) plus an elementary convolution bound; hence the strategy-cube distance delta_k (THM-4479), the periodic-deletion price and its chain (THM-4485), the expanding-necklace count, and the gamma = 1 dip count of THM-4487 are all Theta(2^(hk) k^(-3/2))"
status: >
  PROVED (self-contained, elementary) + FINITE-EXACT. Setting: T the
  Collatz shortcut, h = h(log_3 2) = 0.9499555, rho = log_3 2; for a residue
  r mod 2^k with parity word w and odd counts o_j, M_j = 3^(o_j)/2^j.
  W_k = #{w in {0,1}^k : 3^(o_j) > 2^j for all 1 <= j <= k} (= |Bad_k| of
  THM-4479, the residues with no k-step descent), B_n = sum_(j : 3^j > 2^n) C(n, j),
  N_k = # necklaces of length k with more than k rho ones.
  (A) Identity: k W_k = sum_(n=1)^k B_n W_(k-n), W_0 = 1, i.e.
  sum_k W_k t^k = exp(sum_n B_n t^n/n). Proved in four elementary steps
  (minimum decomposition, reversal, ladder blocks, rotation averaging)
  using only that log_2 3 is irrational (no ties); it is Spitzer's 1956
  combinatorial lemma for two letters. Checked against a ballot DP for
  k <= 300; the recurrence is integral to k = 3000; W_1..W_11 are
  THM-4479's |Bad_k| table.
  (B) Order: 0.26 * 2^(hk) k^(-3/2) <= B_k/k <= W_k <= 545 * 2^(hk) k^(-3/2)
  for every k >= 1 (lower bound from the n = k term of (A) and the Stirling
  bound C(k, j_0) >= 2^(k h(j_0/k))/sqrt(8 k p(1-p)); upper bound from the
  exponential of a series with coefficients g_n = B_n 2^(-hn)/n <= 2.05 n^(-3/2)
  through the convolution lemma (a*b)_k <= 2^(3/2)(A beta + B alpha) k^(-3/2),
  giving w_k <= C e^(2 sqrt2 sigma) k^(-3/2), sigma = sum g_n = 1.9738).
  Also N_k >= C(k, j_0)/k >= 0.26 * 2^(hk) k^(-3/2). The normalised B_k and
  N_k oscillate with the fractional part of k rho in a window whose ends
  differ by rho/(1-rho) = 1.709 and have no limit; W_k k^(3/2) 2^(-hk) is
  observed in [9.6648, 11.0517] for 500 <= k <= 3000.
  (C) Consequences. (C1) THM-4479's sandwich N_k <= delta_k <= |Bad_k| gives
  delta_k = Theta(2^(hk) k^(-3/2)) with explicit constants: the sharp rate
  -(1-h) of THM-4479 becomes a sharp order, and THM-4479's 2^(hk)/(3k^2)
  lost k^(1/2) only through the crude binomial bound. (C2) THM-4485's chain
  N <= nu <= FVS_(log_3 2)(k) <= FVS^odd <= delta_k puts the cycle-packing
  number and the periodic-deletion prices in the same window. (C3) On both
  sheets D_b(X, 1) = #{n <= X : T_b^j(n) >= n, j <= floor(log_2 n)} =
  Theta(X^h (log X)^(-3/2)): THM-4487's bracket log^(-3/2) .. log^(+1)
  closes at its lower end (plus sheet: positive words are non-dippers and
  the rest inject into positive words by prepending two odd letters; minus
  sheet: non-dippers are positive words, and 11w is a non-dipper).
  FINITE-EXACT: sum_(t<T) W_t = 281, 2903, 31730, 367698 for T = 12, 16, 20, 24
  are exactly THM-4487's brute-force counts D_+-(2^T, 1) on both sheets.
  NOT claimed: anything about Collatz orbits; THM-4476's X^(h+eps) keeps its
  eps; the polynomial factor at gamma < 1 (expected (log X)^(-1/2)); sharp
  constants (observed W_k k^(3/2) 2^(-hk) is about 10-11, N_k's about 1.1-2).
source: collatz-exponent-atlas-20260926 session (opus), 2026-09-26, integrating the procgen lane's THM-4479/THM-4485 (one exponent 1 - h) with THM-4487 (dip spectrum) after reading the day's incoming work; the identity was found by asking for the exact order of |Bad_k| and proved combinatorially.
depends_on:
  - 01-canon/theorems/THM-4479-strategy-cube-distance-to-provability-sharp.md (the sandwich N_k <= delta_k <= |Bad_k|)
  - 01-canon/theorems/THM-4485-periodic-edit-price-feedback-sets.md (the chain delta_k >= FVS^odd >= FVS >= nu >= N)
  - 01-canon/theorems/THM-4487-dip-spectrum-entropy-curve-and-sharpness-of-thin-divergence.md (the gamma = 1 dip count)
related:
  - 05-knowledge/results/collatz_nodescent_order_20260926_spitzer_ballot.md (full note with proofs and controls)
  - 01-canon/theorems/THM-4476-thin-divergent-orbits-reciprocal-sums-finite.md (1 - h as the Chernoff rate)
  - 05-knowledge/results/collatz_dipspectrum_20260926_entropy_curve.md (THM-4487's note; its gamma = 1 counts are reproduced exactly)
script: 04-computation/experiments/collatz_nodescent_order_20260926.py
output: 04-computation/experiments/collatz_nodescent_order_20260926.out
script_sha256: 9e3d674295ae19206c7ff5a3c9e053b372e9dc25273c4e8015f787299a1b148a
output_sha256: 9f70677d1971708c6952ea21ada16abfb67a896cc243ba4e83104e06fa45a7e1
hash_basis: raw LF bytes
audit: >
  Self-audited: the four steps of the identity re-derived with attention to
  ties (none, by irrationality); the identity verified exactly against an
  independent ballot DP (k <= 300) and by integrality (k <= 3000); the
  binomial bounds, the geometric tail and h(j_0/n) <= h checked at every
  n <= 3000; the convolution induction re-checked; the dip-count corollary
  checked against THM-4487's brute-force counts to 2^24 (exact equality).
  Independent audit not yet performed.
---

# THM-4495 -- the no-descent count has exact order 2^(hk) k^(-3/2)

**PROVED (self-contained) + FINITE-EXACT.** Full note:
[collatz_nodescent_order_20260926_spitzer_ballot](../../05-knowledge/results/collatz_nodescent_order_20260926_spitzer_ballot.md).

## 1. Statement

With `W_k`, `B_n`, `N_k` as in the status block:

```text
k W_k = sum_(n=1)^(k) B_n W_(k-n)   (W_0 = 1),            [identity]
0.26 * 2^(hk) k^(-3/2) <= B_k/k <= W_k <= 545 * 2^(hk) k^(-3/2),   N_k >= 0.26 * 2^(hk) k^(-3/2),   [order]
delta_k, FVS_(log_3 2)(k), FVS^odd, nu, N_k, |Bad_k| = Theta(2^(hk) k^(-3/2)),   D_(+-)(X, 1) = Theta(X^h (log X)^(-3/2)).
```

## 2. Why the identity holds (one paragraph)

Positive words (all prefix sums `S_j = o_j log_2 3 - j > 0`) split uniquely at
their minimum into a positive min-ending word and a positive word, so
`W = 1/(1 - P)`; reversal identifies min-ending positive words with
first-passage words (`S_j < 0` before `S_m > 0`), so `W = 1/(1 - L)`; words
ending at their maximum are concatenations of first-passage words cut at
the ladder epochs, so `[t^n] log W = sum_r E_(n,r)/r`; and for a word `x`
with positive total, extended periodically, the rotations that end at their
maximum are exactly those starting at a record residue, each with `r(x)`
ladder epochs, so `sum_i [rot_i x ends at max]/r(rot_i x) = 1` and
`B_n = n sum_r E_(n,r)/r`. No ties anywhere because `log_2 3` is irrational.

## 3. Why the order holds (one paragraph)

`B_n 2^(-hn) = Theta(n^(-1/2))` by Stirling and the geometric tail beyond
`j_0 = floor(n rho) + 1` (ratio `(1-rho)/rho`). The recurrence's `n = k` term
gives `W_k >= B_k/k`. For the upper bound, `sum_k W_k 2^(-hk) u^k = exp(G(u))`
with `G` having nonnegative coefficients `<= 2.05 n^(-3/2)`; the split-at-`k/2`
convolution lemma bounds the `m`-fold convolution by `C m (2^(3/2) sigma)^(m-1) k^(-3/2)`,
and the `1/m!` of the exponential sums the series to `C e^(2 sqrt2 sigma) k^(-3/2)`.

## 4. What it changes upstream

* THM-4479 (procgen): `delta_k = Theta(2^(hk) k^(-3/2))`, not only
  `(1/k) log_2(delta_k/2^(k-1)) -> -(1-h)`; the necklace lower bound improves
  from `2^(hk)/(3k^2)` to `0.26 * 2^(hk) k^(-3/2)`.
* THM-4485 (procgen): the whole chain `N <= nu <= FVS <= FVS^odd <= delta_k <= |Bad_k|`
  lies in a constant-factor window; the `1.35` of the refuted Golomb analogue
  is a point inside it.
* THM-4487 (this lane): the `gamma = 1` count is `Theta(X^h (log X)^(-3/2))`
  on both sheets, and `|Bad_k|` is computable for every `k` from binomial
  tails alone (no enumeration).
* Reading: `1 - h = I(g)` (Chernoff rate, THM-4476 section 1.8) and
  `k^(-3/2)` (ballot correction of a zero-drift walk pinned at bounded
  height) are the common exponent and polynomial of four prices. Collatz
  itself is untouched.
