---
id: HYP-9374
title: "Mixed-order clone depth equals the candidate minimum: for a regular action of a finite group A of odd order p^e m (p not dividing m) on itself, the order of vanishing at t=-1 of the subset-orbit polynomial P_A(t) reduced modulo p equals min({|A|/o' : p does not divide E_{o'}/p^e} union {|A|/o' - p^j + p^(j-1) : p does not divide e_{j,o'}/p^(j-1)}), where o' runs over divisors of m, E_{o'} counts elements whose p'-part has order o', and e_{j,o'} counts elements of order p^j o'. The inequality >= is THM-4764 Theorem B''; equality would fail only by a collision of candidate orders from different p'-parts with cancelling leading coefficients."
status: >
  OPEN (klein-2026-10-25). PROVED as a lower bound (THM-4764, Theorem B''), PROVED with equality for p-groups
  (m = 1; THM-4764, Theorem B) and for primes not dividing |A| (Theorem B', where the candidate set is
  {|A|/o : l does not divide e_o}). FINITE-EXACT equality on all 22 tested regular actions, including the mixed-order
  groups C15, C21, C45, C7 x| C3 and C3 x C5^2 at every prime dividing their order (clone-normal.json, field
  candidate_lower_bound). Direct products A = S x H reduce to d_p = p^e m/o* - p^J + p^(J-1) with
  o* = max{o' : p does not divide e^H_{o'}} and J = J(S), absent collisions; C7 x| C3 needs the general form
  (E_3 = 14 has 7-adic valuation 1 = e, giving candidate 7, the actual depth).
session: klein-2026-10-25 (codex/session-tournament-twist-unify)
source: 05-knowledge/results/clone_depth_modular_twist_average_20261025/README.md
scripts:
  - 05-knowledge/results/clone_depth_modular_twist_average_20261025/check_clone_depth.py (function mixed_candidates)
related:
  - 01-canon/theorems/THM-4764-modular-cancellation-depth-of-odd-regular-gap-actions-and-the-fifty-vertex-realization.md
  - 05-knowledge/results/session62_agent_theme_bridges_20261025/clone/README.md
---

# HYP-9374 — Mixed-order clone depth equals the candidate minimum

## Why it should be true

Modulo `p`, `P_A` is a sum over `p'`-parts `o'` of blocks `(E_{o'}/|A|) s'^{|A|/o'}·(unit)` plus `μ=1` terms
`(e_{j,o'}/φ(p^j)) (p-1) s'^{|A|/o'-p^j+1} g_j(s')`, `s' = 1 + t^{o'} = o' s (1 + …)`; all `μ ≥ 2` terms vanish
(THM-4764, proof of B″). The surviving terms have the listed `s`-orders with unit leading coefficients. Equality
fails only if two surviving terms from different `o'` share an order and their leading coefficients cancel
modulo `p`. Orders of distinct blocks differ by at least `|A|(1/o'_1 - 1/o'_2)`, which exceeds the correction
`p^j - p^(j-1) < p^e` whenever `o'_1 | o'_2` (ratio at least 3), so collisions need incomparable `p'`-orders
with close reciprocals.

## Cheapest decisive test

Search odd-order groups with two incomparable `p'`-element orders `o'_1 < o'_2` such that
`|A|/o'_1 - p^j + p^(j-1) = |A|/o'_2` for some `j ≤ e`, compute the two leading coefficients modulo `p`, and
compare the actual depth with the candidate minimum. The smallest such `|A|` is unknown to this session.

## What a proof must retain

The `p`-adic valuations of `E_{o'}` and `e_{j,o'}` (Frobenius' theorem supplies `p^e | E_{o'}`), the exact
unit `u_j = C(p^j,p^(j-1))/p`, and the leading coefficient `o'^{|A|/o'}` of `s'`. A proof by "distinct orders"
alone is insufficient: in `C21` modulo 3 the order-21 and order-7 contributions to `[s^3]` are `760` and `98`
and cancel (`858 ≡ 0`), the mechanism that makes mixed order different from the `p`-group case.
