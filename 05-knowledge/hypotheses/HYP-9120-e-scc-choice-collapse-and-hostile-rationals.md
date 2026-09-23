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
  thinner than the no-choice set (dimension 0.95); dimension 0 versus about 0.1 is OPEN. PROVED: -1, -13/9, 1, 1/2 are hostile; the 1-escape
  lemma; descent by at most two reverse moves (factor <= 8/9, depending only on m mod 27)
  for every m not 1 or 14 mod 27, Lean-checked (q2_descent_off_1_and_14), so Q2 reduces
  to the two hostile neighbourhoods; Q2 verified to 2.02e13; escaping 1/2 costs at least 8/3. FINITE-EXACT: exceptional counts 908 mod 2^36 and 52 mod 3^17; no
  positive integer below 1.5e8 is exceptional at precision 2^36.
source: collatz-procgen-20260922 (mac-mini)
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
