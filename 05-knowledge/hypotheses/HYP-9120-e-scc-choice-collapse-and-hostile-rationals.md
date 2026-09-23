---
id: HYP-9120
title: "E-SCC via choice: the E-graph exceptional sets are thin sets of rationals over the other prime, and an Applegate--Lagarias induction closes"
status: >
  OPEN HYPOTHESIS with FINITE-EXACT support and PROVED partial lemmas.
  (a) Q1: no positive integer n>1 lies in the forward exceptional set
  Bad_inf(E) in Z_2; (b) Q2: no integer m>1 lies in the backward exceptional
  set in Z_3^x; (c) structure: every rational point of Bad_inf(E) is negative
  with a power-of-3 denominator near -1 (observed range [-1.6,-1]), every rational point of the
  backward set is 1, 1/2 or a positive dyadic rational, and both sets are far
  thinner than the no-choice set (dimension 0.95); dimension OPEN; exact threads to 2^64 favour 0 (power law about m^1.6); Bad_inf(E) is infinite (Theorem F: -1-2^i/3^alpha(i) hostile for all i>=3). PROVED: -1, -13/9, 1, 1/2 are hostile; the 1-escape
  lemma; descent by at most two reverse moves (factor <= 8/9, depending only on m mod 27)
  for every m not 1 or 14 mod 27, Lean-checked (q2_descent_off_1_and_14), so Q2 reduces
  to the two hostile neighbourhoods; Q2 verified to 7.87e17; escaping 1/2 costs at least 8/3. FINITE-EXACT: exceptional counts 908 mod 2^36 and 52 mod 3^17; no
  positive integer below 1.5e8 is exceptional at precision 2^36.
source: collatz-procgen-20260922 (mac-mini); E is Le--Smith's Loosened Collatz Graph (arXiv 2109.01180), where Q1/Q2/SCC are not stated
depends_on:
  - 05-knowledge/results/collatz_procgen_20260922_choice_ladder.md
related:
  - 05-knowledge/results/collatz_mod6_20260917_extended_collatz_scc.md
  - 05-knowledge/results/collatz_mod6_20260917_three_adic_g_map.md
---

# HYP-9120 -- choice collapse and hostile rationals for E-SCC

See the lane note for definitions and data. The decisive next tests are:

1. Prove membership of the whole family `-p/3^j` near `-1` by the
   carry-lower-bound argument used for `-13/9`, or find a member that
   escapes.
2. Push the backward thread computation to `3^25` with a tree-restricted
   DFS and check that no new rational family appears.
3. Assemble one see-saw induction step, in the manner of
   Applegate--Lagarias, for Q2 near `1` and `1/2`, using the 1-escape lemma
   and a bounded-lookahead descent factor.

Refutation form: an integer `n>1` whose class is exceptional at every
precision, or an exceptional thread converging to a positive integer.

## Q2 status after the loops lane (same session)

* The thread of `1` is settled for every depth `<=6001`, since the 1-escape
  loops exist for `s<=6000`.
* The `1/2` escape price is exactly `2c(k-1)/3`, in `(1,2)`. A one-round
  see-saw cannot close. Chains of hostile transfers raise the value by at
  most about `m^0.104` (known thread prices), and the endgame is governed
  by the **base-3 digits of `2^K`**: an Erdős-ternary-type transversality,
  the same shape as the divergence-half target of the synthesis.

## Dimension-lane update (same session)

* **PROVED:** Theorem F makes `Bad_inf(E)` infinite, with `-1` an
  accumulation point. Theorem P, the perturbation lemma, generates the
  hostile rationals.
* **FINITE-EXACT:** 382 hostile points `-p/3^j` with `j<=27` and
  `1<=|x|<3/2` are certified. The counts per `j` stay bounded (11--29),
  and the points form a perturbation tree of depth at most 7.
* Q2 is verified below `7.87*10^17`.
* The growth evidence favours dimension `0`. Conjecture G (a finite seed
  generates every hostile rational) is verified for `15<=j<=27`.

See [exceptional_dimension](../results/collatz_procgen_20260922_exceptional_dimension.md).

