---
id: HYP-9137
title: "The sharp price of provable descent: delta_L = 2^(-(1-h)L + o(L)), h = h(log_3 2), matching the undecided density up to subexponential factors"
status: >
  PROVED by THM-4478-critical-tube-affine-capacity-sharp-provability-price,
  2026-09-26, independently audited. The lower bound is
  2^(-(1-h)L-O(sqrt(L log L))), even for lower natural flip density.
  Together with THM-4475 this proves the stated sharp exponential rate,
  not constant-factor or polynomial-factor comparability. Collatz is OPEN.
  Historical route: THM-4475 gives delta_L <= 2 rho_L <= 2^(1-(1-h)L) (upper
  bound, explicit construction) and delta_L >= rho_L/W_L =
  2^(-(0.7737+o(1))L) (lower bound). The gap is the factor W_L ~ lambda^L,
  lambda = (3+sqrt 13)/4, which bounds how many bad orbits a single flip at
  a "hub" v can serve. Evidence: the greedy G_L pays 0.86-0.92 rho_L, and
  the brute-force hub weight over v <= 20000 is 80.6 against the tree
  bound 99.3 at L = 8.
  UPDATE 2026-09-25 (THM-4477): the proved lower bound is now
  2^(-(0.1445+o(1))L), via Cauchy-Schwarz with a certified second-moment
  majorant. Any distribution-only argument is stuck at 2(1-h) = 0.1001
  (moment-method limit, forced by g(2) = 1). The harmonic bound
  H_L = E[1_Bad/max W_L] ~ 3.2 rho_L/M2(L) (L <= 16) is the proposed route
  to the sharp exponent. It needs quenched independence of 2-adic badness
  and 3-adic hub weight (condition (SD) of the Cauchy-Schwarz note).
source: collatz-procgen-20260922 session, price lane (2026-09-25)
depends_on:
  - 01-canon/theorems/THM-4475-price-of-provable-descent-tends-to-zero.md
  - 01-canon/theorems/THM-4478-critical-tube-affine-capacity-sharp-provability-price.md
---

# HYP-9137 -- the sharp price of provable descent

**PROVED, 2026-09-26:** [THM-4478](../../01-canon/theorems/THM-4478-critical-tube-affine-capacity-sharp-provability-price.md).
Critical growth bands retain `2^(hL-o(L))` words. At fixed time and slope,
their affine offsets confine actual integer ancestors to an interval of
length at most `L/3`. A first-hit capacity cut then gives
`delta_L >= 2^(-(1-h)L-O(sqrt(L log L)))`. This bypasses distribution-only
moment estimates. The original proposed necessities below were sufficient
routes, not necessary conditions: no quenched-independence theorem is used.

**Statement.** `delta_L >= 2^(-(1-h)L - o(L))`, so that with THM-4475,
`delta_L = 2^(-(1-h)L + o(L))`: the price of bounded-lookahead provability
equals Collatz's undecided density up to subexponential factors.

**Refutation form.** A family of members of `P_L` with flip density
`2^(-(1-h+epsilon)L)` for a fixed `epsilon > 0`. Such a family would flip at
rare "hubs" whose backward trees carry exponentially many bad orbits.

**Historical proposed proof route (superseded as a requirement).**
* Hubs are rare: the number of `Bad_L` orbits served by one flip is
  `O(poly L)` on average. This is a second-moment or 3-adic equidistribution
  statement about backward trees.
* Flips at great height do not finish the job: after a flip at height `R`,
  the orbit is still at height about `R/2`.

**Why it matters.** This makes "the price of provability is the
undecided density" an exact law, `delta_L = rho_L^(1+o(1))`. That law ties
the pairing ladder (THM-4470/4475), the strategy cube (THM-4474, Theorem A)
and the exceptional-set dimension `h(log_3 2)` into one statement.
