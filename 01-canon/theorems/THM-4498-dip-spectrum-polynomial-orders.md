---
id: THM-4498
title: "The dip spectrum with its polynomial orders: D_b(X, gamma) = Theta(X^(h(gamma/log_2 3)) (log X)^(-1/2)) for log_4 3 < gamma < 1, Theta(X^(h*) (log X)^(-3/2)) at gamma = 1, and Theta(X) for 0 < gamma <= log_4 3"
status: >
  PROVED (elementary; CITED Hoeffding 1963 Theorem 4 for sampling without
  replacement and the central limit theorem) + FINITE-EXACT controls.
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
  (c) 0 < gamma <= log_4 3: Theta(X): a fifth of all words work (Hoeffding
  for independent letters plus the CLT for the endpoint; multiplicative
  carry of THM-4487's Theorem 4 for gamma <= log_2(3/2)).
  So the polynomial factor along the entropy curve E(gamma) = h(max(1/2, gamma/alpha))
  is 1 up to Korec's endpoint, (log X)^(-1/2) between Korec and Terras, and
  (log X)^(-3/2) at Terras's endpoint; the constant of (a) blows up like
  1/(1 - (1-rho)/rho) as gamma decreases to log_4 3. FINITE-EXACT: exact DP
  word counts W_t(gamma) t^(1/2) 2^(-tE) to t = 3000 (settled for gamma >= 0.94,
  converging slowly near log_4 3, as the geometric ratio predicts), and the
  brute-force orbit counts of THM-4487 to 2^24 renormalised. NOT claimed:
  sharp constants; anything about Collatz orbits.
source: collatz-exponent-atlas-20260926 session (opus), 2026-09-26; completing THM-4487's polynomial bracket after THM-4495 settled gamma = 1.
depends_on:
  - 01-canon/theorems/THM-4487-dip-spectrum-entropy-curve-and-sharpness-of-thin-divergence.md (Theorem 1's exponent, the Terras bijection, the carry bounds, Theorem 4's multiplicative carry)
  - 01-canon/theorems/THM-4495-no-descent-count-exact-order-spitzer-ballot.md (the gamma = 1 order)
related:
  - 05-knowledge/results/collatz_dipspectrum_20260926_entropy_curve.md (proof: section 1.6 and 1.2b)
  - 05-knowledge/results/collatz_thin_20260925_thin_divergent_orbits.md (Lemma 1.4b: the same geometric tail turns THM-4476's eps into (log X)^(0.014+eps))
script: 04-computation/experiments/collatz_dipspectrum_20260926_orders.py
output: 04-computation/experiments/collatz_dipspectrum_20260926_orders.out
script_sha256: 87f31934e76b6c958038ce46b4d26cd1de111e602c9ec08418d97275d8109225
output_sha256: 270c93827ae4e7e2da91bcc49e43505093664d6d850540b6fd9bee3727b2873f
hash_basis: raw LF bytes
audit: >
  Self-audited: the two Hoeffding applications (without replacement for
  (a), independent for (c)) re-derived with the block sums and the margins
  written out; the prefix-sum inequality for w = 1^K v and the carry step
  checked; the exact DP counts confirm the t^(-1/2) window for gamma >= 0.94
  and the slow convergence predicted near log_4 3. Independent audit not
  yet performed.
---

# THM-4498 -- the dip spectrum with its polynomial orders

**PROVED + FINITE-EXACT.** Proof: section 1.6 (with 1.2b) of
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
