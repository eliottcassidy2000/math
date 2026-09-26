---
id: THM-4499
title: "Thin divergence is o(X^(h*)): every non-eventually-periodic orbit of x -> x/2, (3x+b)/2 (b odd), and every injective invariant set, has N(X) <= K X^(h*) (log_2 X)^a for every a > lambda*/h* - 3/2 = -0.9862; the Terras exponent h* = h(log_3 2) is not attained"
status: >
  PROVED (elementary; rests on THM-4476's recursion (R) and THM-4495's
  Spitzer/convolution apparatus) + FINITE-EXACT controls. Setting as in
  THM-4476: T_b(x) = x/2 (x even), (3x+b)/2 (x odd), (x_i) an orbit in
  Z\{0} that is not eventually periodic, N(X) = #{i : |x_i| <= X},
  h* = h(log_3 2) = 0.9499555, lambda* = log_2(rho*/(1-rho*))/log_2 3 =
  0.488077 (the tilt -E'(1) of THM-4487).
  Theorem: for every a > a* := lambda*/h* - 3/2 = -0.986211 there is
  K(a, |b|) with N(X) <= K X^(h*) (log_2 X)^a for all X >= 2; in particular
  N(X) = o(X^(h*)). Same for every T_b-invariant set on which T_b is
  injective (THM-4476, Cor. 4).
  Mechanism. (1) Lemma M (moving-barrier ballot bound): the number M_k(y)
  of words of length k whose walk S_i = o_i log_2 3 - i stays above -y is
  at most D_s 2^(hk) k^(-3/2) 2^(lambda* y) e^(sy) for every s > 0. Proof:
  decompose a word at its minimum (M_k(y) = W_k + sum_m N_m(y) W_(k-m),
  N_m(y) = negative words of length m with total > -y; exact identity,
  checked to k = 60); tilt to the zero-drift Bernoulli(log_3 2) law
  (Q(w) = 2^(-hm) 2^(lambda* S_m)) and pay e^(sy) for the endpoint; bound
  E_Q[e^(sS_m); all partial sums < 0] by THM-4495's convolution argument
  applied to the WEIGHTED Spitzer identity n P_n = sum B_m P_(n-m) (same
  bijective proof with multiplicative letter weights; checked exactly for
  rational weights, both signs, n <= 40), with the local bound
  E_Q[e^(sS_n); S_n < 0] <= 3.2 n^(-1/2)/(1 - e^(-s log_2 3)); assemble with
  THM-4495's W_(k-m) <= 545 2^(h(k-m)) (k-m)^(-3/2) and the convolution
  lemma. (2) Lemma 1.4c: #F_b(X, theta) <= 2|b| X^(0.585+theta) +
  4 D_s X^(h*) (log_2 X)^(-3/2) 2^(lambda* theta log_2 X) e^(s theta log_2 X)
  (THM-4476's counting lemma with the ballot factor). (3) THM-4476's
  bootstrap with theta_X = (1+eta) log_2 log_2 X/(h* log_2 X): the L
  dippers per landing point cost c_1 h* = 1 + eta and reappear as
  (log X)^(lambda* c_1); the ballot factor contributes -3/2; hence a*.
  History: THM-4476 gave X^(h*+eps); addendum 1.6b of its note (same day)
  gave (log X)^(0.014+eps) via a geometric binomial tail; this theorem
  adds the ballot factor of THM-4495 at a moving barrier. Observed:
  M_k(y) ~ (y+1) c 2^(hk) k^(-3/2) 2^(lambda* y) for y << k^(1/2), so the
  e^(sy) slack is not the truth and a* is the exponent the linear law
  would also give. NOT claimed: anything excluding divergent orbits;
  N(X) = O(log X) (expected for an actual divergent orbit) is far away,
  since the method uses one free window per element and a landing
  multiplicity L; a multiplicity L^(1/2) would give a* = lambda*/(2h*) - 3/2.
source: collatz-exponent-atlas-20260926 session (opus), 2026-09-26; found by asking whether THM-4495's ballot factor survives the moving barrier of THM-4476's recursion. No priority claimed.
depends_on:
  - 01-canon/theorems/THM-4476-thin-divergent-orbits-reciprocal-sums-finite.md (reductions, Terras bijection, carry bound, recursion (R), Cor. 4)
  - 01-canon/theorems/THM-4495-no-descent-count-exact-order-spitzer-ballot.md (Spitzer identity, convolution lemma, W_k <= 545 2^(hk) k^(-3/2))
related:
  - 05-knowledge/results/collatz_thin_20260926_little_o_thin_divergence.md (full proof and controls)
  - 05-knowledge/results/collatz_thin_20260925_thin_divergent_orbits.md (THM-4476's note; addendum 1.6b superseded)
  - 01-canon/theorems/THM-4487-dip-spectrum-entropy-curve-and-sharpness-of-thin-divergence.md (the tilt lambda*; the sharpness of the counting lemma, which this theorem does not contradict: the ballot factor is a polynomial, not an exponent)
script: 04-computation/experiments/collatz_thin_20260926_movingbarrier.py
output: 04-computation/experiments/collatz_thin_20260926_movingbarrier.out
script_sha256: da2d12874ab368dbff10136fea23e780687633e27b92c38eaf3d3963df05722e
output_sha256: 0b2c11a427dd218cf4edef6a21d0823aa9ad9c41303dbe91383380c8584878b6
hash_basis: raw LF bytes
audit: >
  Self-audited: the decomposition (D) re-derived and checked exactly; the
  weighted Spitzer identity checked exactly with rational weights; the
  tilt identity, the endpoint weight, the local binomial bound and both
  applications of the convolution lemma re-derived; the bootstrap
  bookkeeping written out for negative a (log_2 Y in [L/2, L]).
  Independent audit not yet performed.
---

# THM-4499 -- thin divergence is o(X^(h*))

**PROVED + FINITE-EXACT.** Full note:
[collatz_thin_20260926_little_o_thin_divergence](../../05-knowledge/results/collatz_thin_20260926_little_o_thin_divergence.md).

## 1. Statement

```text
N(X) <= K(a, |b|) X^(h*) (log_2 X)^a      for every a > lambda*/h* - 3/2 = -0.9862,   X >= 2,
```

for every non-eventually-periodic `T_b`-orbit and every injective
`T_b`-invariant set. In particular `N(X) = o(X^(h*))`.

## 2. The three steps

1. **Moving-barrier ballot bound.** Words staying above `-y` number at
   most `D_s 2^(hk) k^(-3/2) 2^(lambda* y) e^(sy)`: decompose at the
   minimum into a negative word with total `> -y` and a positive word;
   count the negative words by the weighted Spitzer identity and the
   convolution lemma of THM-4495; the tilt pays `2^(lambda* y)`.
2. **Counting lemma with the ballot factor.** THM-4476's no-dip set has
   `#F_b(X, theta) <= 4 D_s X^(h*) (log_2 X)^(-3/2) 2^(lambda* theta log_2 X) e^(s theta log_2 X)` plus the small-element term.
3. **Bootstrap with a moving `theta`.** `theta_X = (1 + eta) log_2 log_2 X/(h* log_2 X)`
   makes the dipper term geometric; the no-dip term is
   `X^(h*) (log_2 X)^(lambda*(1+eta)/h* - 3/2 + o(1))`.

## 3. What it says and does not say

The exponent `h*` of THM-4476 is an upper bound that no divergent orbit
attains; the count is smaller by `(log X)^(0.986 - eps)`. Divergent orbits
are not excluded, and the expected `O(log X)` is out of reach of the
one-window method. The sharpness of THM-4487 (the no-dip set at a fixed
`theta` has exponent `h(rho) > h*`) is untouched: this theorem lives in the
polynomial factor at `theta -> 0`.
