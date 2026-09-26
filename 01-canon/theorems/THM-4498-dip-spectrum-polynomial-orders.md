---
id: THM-4498
title: "The dip spectrum with its polynomial orders: D_b(X, gamma) = Theta(X^(h(gamma/log_2 3)) (log X)^(-1/2)) for log_4 3 < gamma < 1, Theta(X^(h*) (log X)^(-3/2)) at gamma = 1, and Theta(X) for 0 < gamma <= log_4 3"
status: >
  PROVED (elementary; CITED Hoeffding 1963 Theorem 4 for sampling without
  replacement) + FINITE-EXACT controls + INDEPENDENTLY AUDITED (2026-09-26:
  HAS GAPS in the proof text, repaired; theorem true as stated; see the
  audit field).
  Setting: b = +-1, T_b(x) = x/2 or (3x+b)/2, alpha = log_2 3, h the binary
  entropy, D_b(X, gamma) = #{n <= X : T_b^i(n) >= n^gamma for 0 <= i <= floor(log_2 n)}.
  (a) For log_4 3 < gamma < 1 and rho = gamma/alpha:
  c_1(gamma) X^(h(rho)) (log_2 X)^(-1/2) <= D_b(X, gamma) <= c_2(gamma) X^(h(rho)) (log_2 X)^(-1/2)
  for X >= X_0(gamma). Upper bound: a geometric binomial tail replaces
  THM-4487's (t+1) max_o C(t,o) (note, section 1.2b). Lower bound: prepend K
  odd letters to a word with ceil(rho t') odd letters whose final rise
  (largest terminal-block sum) is at most a constant; Hoeffding without
  replacement shows at least half of the arrangements qualify, so a block
  contributes >= C(t', j)/2 = Theta(2^(t h(rho)) t^(-1/2)) non-dippers.
  (b) gamma = 1: Theta(X^(h*) (log_2 X)^(-3/2)), THM-4495 (Corollary C3).
  (c) 0 < gamma <= log_4 3: Theta(X), positive lower density: a quarter of
  the tails v, i.e. a fraction 2^(-K)/4 (K >= 40) of all words, work
  (Hoeffding for independent letters plus the symmetry of the binomial for
  the endpoint; multiplicative carry of THM-4487's Theorem 4 for
  gamma <= log_2(3/2)).
  So the polynomial factor along the entropy curve E(gamma) = h(max(1/2, gamma/alpha))
  is 1 up to Korec's endpoint, (log X)^(-1/2) between Korec and Terras, and
  (log X)^(-3/2) at Terras's endpoint; the constant of (a) blows up like
  1/(1 - (1-rho)/rho) as gamma decreases to log_4 3. FINITE-EXACT: exact DP
  word counts W_t(gamma) t^(1/2) 2^(-tE) to t = 3000 (settled for gamma >= 0.94,
  converging slowly near log_4 3, as the geometric ratio predicts), and the
  brute-force orbit counts of THM-4487 to 2^24 renormalised; the audit's
  exact DP to t = 10^4 gives decade exponents 0.467/0.467/0.485 at
  gamma = 0.91/0.94/0.97 against the claimed 1/2. NOT claimed: sharp
  constants; anything about Collatz orbits.
source: collatz-exponent-atlas-20260926 session (opus), 2026-09-26; completing THM-4487's polynomial bracket after THM-4495 settled gamma = 1.
depends_on:
  - 01-canon/theorems/THM-4487-dip-spectrum-entropy-curve-and-sharpness-of-thin-divergence.md (Theorem 1's exponent, the Terras bijection, the carry bounds, Theorem 4's multiplicative carry)
  - 01-canon/theorems/THM-4495-no-descent-count-exact-order-spitzer-ballot.md (the gamma = 1 order)
related:
  - 05-knowledge/results/collatz_dipspectrum_20260926_entropy_curve.md (proof: section 1.6 and 1.2b)
  - 05-knowledge/results/collatz_thin_20260925_thin_divergent_orbits.md (Lemma 1.4b: the same geometric tail turns THM-4476's eps into (log X)^(0.014+eps))
script: 04-computation/experiments/collatz_dipspectrum_20260926_orders.py
output: 04-computation/experiments/collatz_dipspectrum_20260926_orders.out
script_sha256: ef5662a616822d622625d34f7065ce88462ce195d8354b9e73b24c337af61601
output_sha256: 1633e5e4c6f114a492f7e8617274457ce8ebae60fb175269feab85946dd5db27
hash_basis: raw LF bytes
audit: >
  Self-audited: the two Hoeffding applications (without replacement for
  (a), independent for (c)) re-derived with the block sums and the margins
  written out; the prefix-sum inequality for w = 1^K v checked.
  Independent audit (separate agent, 2026-09-26;
  04-computation/experiments/collatz_dipspectrum_20260926_orders_audit.py ->
  05-knowledge/results/collatz_dipspectrum_20260926_orders_audit.out, sha256
  1baf8ba25c86487b0beb1edda424711642caec5d30b1273aef1e10d0168dc280 /
  ab999da6380dd0aa6514f88ac4de1743897fa83fac572e4fd2bc64d64f706978): 1.2b
  CONFIRMED, constant uniform on compact subsets of (log_4 3, 1], blowing up
  like rho/(2rho-1) with t_0 > 4/(2rho-1); (c) and Remark (i) CONFIRMED
  (delta' >= 29, K >= 40 recomputed; the CLT is unnecessary, binomial
  symmetry gives >= 1/4 for every t'; multiplicative carry checked on
  integers at gamma = 1/2, both sheets). (a) lower bound: ERROR in step (vii)
  as written: with gamma - 1 < 0 the worst n of the block is 2^t, and
  "4 2^((gamma-1)(t+1)) >= 4 n^(gamma-1)" is false for every n in
  [2^t, 2^(t+1)). Repair: the stated condition (c+0.585)K >= delta_0 + c + 2
  gives S_i(w) >= -ct + c + 2 > -ct + 2, hence M_i n >= 4 n^gamma; theorem
  unaffected, text corrected. Hoeffding 1963 Thm 4 correctly cited (exact
  hypergeometric tails <= 0.62 of the bound); union bound, prefix-sum
  inequality, margin, carry (needs only gamma >= log_2(3/2)), Stirling and
  the |h'| <= 1 step (t' >= 28) CONFIRMED. Controls: the float DP's "ties
  impossible" was false: exact ties at o = 0, j = (1-gamma)t were admitted
  by the float comparison at gamma = 0.82, 0.85, 0.94, 0.97; five cells of
  table (1) changed by at most 1.2% (repaired: exact tie test at o = 0,
  output regenerated). Exact DP to t = 10^4: decade exponents
  0.467/0.467/0.485 at gamma = 0.91/0.94/0.97 support t^(-1/2); at 0.82 the
  normalised count (10.1 at 10^4) is 0.91-0.96 of the exact binomial tail,
  which converges slowly to 11.9. Wording corrected: "a fifth of all words"
  is 2^(-K)/4 of all words; "positive density" is positive lower density.
  Verdict: HAS GAPS (theorem true as stated; one false inequality in the
  proof of (a) with a one-line repair; control table slightly wrong in five
  cells); all corrections applied.
---

# THM-4498 -- the dip spectrum with its polynomial orders

**PROVED + FINITE-EXACT + INDEPENDENTLY AUDITED (gaps in the proof text repaired).** Proof: section 1.6 (with 1.2b) of
[collatz_dipspectrum_20260926_entropy_curve](../../05-knowledge/results/collatz_dipspectrum_20260926_entropy_curve.md).

## 1. Statement

```text
D_b(X, gamma) = Theta( X^(h(gamma/alpha)) (log X)^(-1/2) )   for log_4 3 < gamma < 1,
D_b(X, 1)     = Theta( X^(h*) (log X)^(-3/2) )               (THM-4495),
D_b(X, gamma) = Theta( X )                                   for 0 < gamma <= log_4 3.
```

## 2. Mechanism

Upper bounds are binomial tails with a geometric ratio `(1-rho)/rho < 1`
(this is exactly what fails at Korec's endpoint, where the ratio is `1`).
Lower bounds prepend a block of odd letters to a word whose terminal
blocks never sum to more than a constant; Hoeffding's inequality (without
replacement for a fixed number of odd letters, independent for the
uniform case) shows that a constant fraction of arrangements have this
property, and the prepended block supplies the margin that the carries and
the range of `n` in a dyadic block require. At `gamma = 1` the barrier is
not receding and the count is a ballot count (THM-4495).

## 3. What it changes

THM-4487's bracket `[log^(-3/2), log^(+1)]` is now an exact order at every
`gamma`. The same geometric tail is Lemma 1.4b of the thin-divergence note,
where it turns THM-4476's `X^(h*+eps)` into `X^(h*) (log X)^(0.014+eps)`.
