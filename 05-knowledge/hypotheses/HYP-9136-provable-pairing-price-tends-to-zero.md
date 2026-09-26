---
id: HYP-9136
title: "The price of provable descent tends to zero: the minimal flip density delta_L of an L-step-descent pairing near Collatz tends to 0 as L grows"
status: >
  SETTLED -- PROVED and promoted as THM-4475 (2026-09-25). For every L >= 8
  an explicit greedy member G_L of P_L (flipping one free pair per
  Collatz-undecided n: its own pair, the first 3 (mod 4) up-image below 2n,
  or the pair of T(n)) is a tree with flip density <= 2 rho_L <=
  2^(1-(1-h)L), where 1 - h = 0.050044. So delta_L -> 0 at the conjectured
  rate, as an upper bound. The lower bound delta_L >= 2^(-0.774L) is proved;
  the sharp exponent is HYP-9137. The secondary claim (delta_2 = 0.2907...
  as a limit) remains OPEN. Earlier evidence: the landing => down class
  (L = 2) needs at least 0.2908 (CP-SAT window optimum); window designs
  flip 14.5% (L = 4) and 5.9% (L = 5).
source: collatz-procgen-20260922 session, brackets/pairings lane candidates (P1, P2), 2026-09-24
depends_on:
  - 01-canon/theorems/THM-4470-collatz-pairing-ladder-am-fair-and-defect-blind.md
---

# HYP-9136 -- the price of provable descent tends to zero

**Statement.** `delta_L -> 0` as `L -> infinity`. Secondary claim:
`delta_2 = lim_X` (window optimum on `[1, X]`) exists, and equals `0.2907...`.

**Refutation form.** A uniform `delta > 0` such that every infinite pairing
with bounded-lookahead descent differs from Collatz on a set of density at
least `delta`.

**Why it matters, and what it cannot do.**
* Collatz sits at flip density **zero** from a divergent member (THM-4470 §4).
  HYP-9136 would put it at density zero from *provable* trees as well: the
  conjecture would then sit on the boundary between two dense regions of the
  pairing cube.
* This does not help prove Collatz. Density-zero modifications change the
  truth value (DEFECT). HYP-9136 is a statement about the *geometry* of the
  pairing family: it measures how far bounded-lookahead certificates are
  from Collatz's own rising runs, the 2-adic neighbourhood of `-1`.
* Collatz's own failures of `P_4` are `n = 7, 11, 15 (mod 16)` (density
  `3/16`). The price must be paid mainly on the rising runs.
