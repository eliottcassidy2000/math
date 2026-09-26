---
id: HYP-9139
title: "The second moment of the Collatz hub weight grows polynomially: M2(L) = E_Haar[W_L^2] = poly(L), equivalently the ladder cross-terms gamma_k = <Z_k o tau, Z_k> stay bounded (tau(u) = 3u+1)"
status: >
  OPEN HYPOTHESIS; FINITE-EXACT for L <= 16 (M2(16) = 762.55, M2/L^3
  decreasing to 0.186) and for ||Z_k||^2, k <= 15 (increments
  gamma_(k-1)/2 = 0.5, 0.407, 0.382, 0.362, 0.357 at k = 3, 6, 9, 12, 15,
  slowly decreasing). With THM-4477 (Theorem CS) it would give the price
  lower bound delta_L >= 2^(-(2(1-h)+o(1))L) = 2^(-0.1001L), the
  moment-method limit (THM-4477, Theorem H). The certified class
  majorant (theta_15 = 1.0312629) cannot reach theta = 1, because
  class-level Cauchy-Schwarz does not see that the true cross term stays
  bounded.
source: collatz-procgen-20260922 session, Cauchy-Schwarz lane (2026-09-25), open item (O1)
depends_on:
  - 01-canon/theorems/THM-4477-cauchy-schwarz-price-bound-and-am-fair-criticality.md
---

# HYP-9139 -- the hub weight has polynomially bounded second moment

**Statement.** `M2(L) = O(L^C)` for some `C`. The data suggest `M2 ≈ 0.06 L^3`,
i.e. `||Z_k||^2` grows linearly with slope near `0.35`. Equivalently,
`gamma_k = <Z_k∘tau, Z_k>` is bounded (numerically `gamma_15 = 0.713`).

**Equivalent forms (lane note §7).**
* The level energies `E_n = ||pi||_n^2 - ||pi||_(n-1)^2` of the Syracuse
  law `pi` are bounded (`E_1 = 2/3`, `E_2 = 10/21`, `E_15 = 0.4708`, slowly
  rising toward about 0.476).
* The ladder autocorrelations of the backward tree decay geometrically,
  via the exact Poisson-kernel identity
  `||pi||_n^2 = sum P_(1/4)(xi) |nu_hat(xi)|^2`.

**Why it matters.** It is the exact point where `q = 3`'s AM-fairness
(`g(2) = 1`, neutral diagonal) meets the sibling dependence of the Collatz
tree: `3x+1` and `x` are linked. Proving it closes the gap between the
proved exponent 0.1445 and the moment-method limit 0.1001. It is a
statement about 3-adic equidistribution of the Syracuse random variable in
`L^2`, the second-moment cousin of Tao's Proposition 1.17.
