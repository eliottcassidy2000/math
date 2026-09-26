---
id: HYP-9137
title: "The sharp price of provable descent: delta_L = 2^(-(1-h)L + o(L)), h = h(log_3 2), i.e. every L-step-provable pairing must flip a density comparable to Collatz's undecided density"
status: >
  OPEN HYPOTHESIS. THM-4475 gives delta_L <= 2 rho_L <= 2^(1-(1-h)L) (upper
  bound, explicit construction) and delta_L >= rho_L/W_L =
  2^(-(0.7737+o(1))L) (lower bound). The gap is the factor W_L ~ lambda^L,
  lambda = (3+sqrt 13)/4, which bounds how many bad orbits a single flip at
  a "hub" v can serve. Evidence: the greedy G_L pays 0.86-0.92 rho_L, and
  the brute-force hub weight over v <= 20000 is 80.6 against the tree
  bound 99.3 at L = 8.
source: collatz-procgen-20260922 session, price lane (2026-09-25)
depends_on:
  - 01-canon/theorems/THM-4475-price-of-provable-descent-tends-to-zero.md
---

# HYP-9137 -- the sharp price of provable descent

**Statement.** `delta_L >= 2^(-(1-h)L - o(L))`, so that with THM-4475,
`delta_L = 2^(-(1-h)L + o(L))`: the price of bounded-lookahead provability
equals Collatz's undecided density up to subexponential factors.

**Refutation form.** A family of members of `P_L` with flip density
`2^(-(1-h+epsilon)L)` for a fixed `epsilon > 0`. Such a family would flip at
rare "hubs" whose backward trees carry exponentially many bad orbits.

**What a proof needs.**
* Hubs are rare: the number of `Bad_L` orbits served by one flip is
  `O(poly L)` on average. This is a second-moment or 3-adic equidistribution
  statement about backward trees.
* Flips at great height do not finish the job: after a flip at height `R`,
  the orbit is still at height about `R/2`.

**Why it matters.** It would make "the price of provability is the
undecided density" an exact law, `delta_L = rho_L^(1+o(1))`. That law ties
the pairing ladder (THM-4470/4475), the strategy cube (THM-4474, Theorem A)
and the exceptional-set dimension `h(log_3 2)` into one statement.
