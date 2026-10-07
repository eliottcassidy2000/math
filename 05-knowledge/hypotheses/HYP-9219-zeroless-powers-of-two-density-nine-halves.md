---
id: HYP-9219
title: "Zeroless powers of two: Z_k = c (9/2)^k (1 + O(0.4^k)) with c = 0.8876940431151482645, i.e. the 2-adic zeroless-digit measure has a continuous density and the 5-adic exponent set of zeroless trailing digits has Hausdorff dimension exactly log_5(9/2); sufficient: max_t |P_m(zeta_(2^m)^t)| = O(rho^m) for some rho < 9/2"
status: >
  OPEN. NUMERICAL: (2/9)^k Z_k converges to c to 19 digits by k = 40; the 2-adic density's sup and inf stabilize at 1.141346047
  and 0.8870365582 (no change in 10 digits from resolution 2^18 to 2^25); the largest conjugate grows like rho^m with
  rho ~ 2.30, 2.64, 2.87 at m = 10, 20, 30. PROVED brackets: growth in [4.47848, 4.52386], dimension in [0.93156, 0.93782] (THM-4580).
source: mac-mini-2026-10-07-oaimath3 (zeroless reader), 05-knowledge/results/oai3_two_orbits_twos_and_threes_20261007.md
related:
  - 01-canon/theorems/THM-4580-zeroless-powers-of-two-the-5-adic-tree-and-verification-to-1-1e11.md
---

# HYP-9219 — the zeroless density is continuous

**Statement.** Let `μ` be the law of `Σ a_i 10^i` on `Z_2`, with `a_i` i.i.d. uniform on `{1, …, 9}`. Then `μ` has a continuous density with respect to Haar measure. Equivalently, `Z_k = c (9/2)^k (1 + O(0.4^k))`, and the exponent set has dimension `log_5(9/2) = 0.934536`.

**Sufficient condition.** For the cyclotomic units of THM-4580 (3), a bound `max_t |P_m(ζ_(2^m)^t)| = O(ρ^m)` with `ρ < 9/2`.

**Why it may be hard.** Large conjugates need `5^i t mod 2^(m−i)` to be small for many `i` in a row. Bounding those returns looks like a two-logarithm (2 against 5) problem of Baker type. It is the decimal twin of the `2` against `3` problems behind Collatz and Mahler's 3/2 problem.
