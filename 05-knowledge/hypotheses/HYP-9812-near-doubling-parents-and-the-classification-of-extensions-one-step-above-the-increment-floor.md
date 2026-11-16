---
id: HYP-9812
title: "After THM-5262: (a) near-doubling parents - there is an absolute constant C such that every count h >= 5 has a realiser T and a vertex x with |H(T) - 2H(T-x)| <= C (complete to h = 223 with C = 5: defect 1 for 79 counts, 3 for 27, 5 for 41 and 125; order-10 realisers give at most 7 up to 3357); (b) a nontrivial one-vertex extension has 2H(T+x) - 3H(T) = 3 only if H(T) = 2 phi_a phi_b - 1 with phi_0 = 1, phi_a = 2^a + 1 (true for every tournament of order at most 9); (c) the basin densities of the signed descent exist and t rho(t) is a function F of the 3-adic half-count with F(n) = F(2n)/2 + (3/4) F(r(n))"
status: >
  OPEN (CONJECTURE), mac-mini-2026-11-15-descent. (a) FINITE-EXACT to h = 223 (complete parent sets), NUMERICAL to
  3357 (upper bounds from order-10 realisers: 1,351 of 1,567 counts with defect 1, 194 with 3, 18 with 5, and 399,
  1285, 2005, 2569 with 7). If true it is a redundant binary expansion of the half-count in which every digit is a
  vertex; it does NOT by itself give Chan-Pak, because the parent has to be found before the child. (b)
  FINITE-EXACT (extensions of every tournament of order at most 9); the constructions for every phi_a phi_b are PROVED (THM-5262 (C)); the
  first test that separates products of two from products of three such numbers is 149 -> 225 at order 11.
  (c) NUMERICAL (X = 2*10^7); the class means mod 9 forced by the functional equation agree to three digits.
  Chan-Pak remains OPEN.
session: mac-mini-2026-11-15-descent
date: 2026-11-15
source: 05-knowledge/results/chanpak_signed_descent_lifts_20261115/README.md
scripts:
  - 05-knowledge/results/chanpak_signed_descent_lifts_20261115/scripts/lift_par.c, lift_near.c, lift_gap.c,
    parents_closed.py, dbasin.py
related:
  - THM-5262, HYP-9469, THM-4876, six_tutte_collatz_20261116/COUNTS.md
---

# HYP-9812 - near-doubling parents, one step above the floor, and the basin function

Source: [the note](../results/chanpak_signed_descent_lifts_20261115/README.md), sections 3-5.

**(a)** Every count `h >= 5` has a realiser `T` and a vertex `x` with `|H(T) - 2 H(T - x)| <= C`.
Falsifier: a count all of whose parents are far from `h/2`. Data: `out/parents_closed_223.txt`,
`out/near_doubling_summary.txt`.

**(b)** If `x` is a new vertex with `H(T + x) > H(T) >= 5` then `2H(T + x) - 3H(T) >= 3`, with equality only for
`H(T) = 2 phi_a phi_b - 1`. (The inequality `>= 1` is the forced-shuffle floor; the case of equality 1 is
`3 -> 5` alone through order 10.) Falsifier: a pair `(2m - 1, 3m)` with `m` not of that form, first candidate
`m = 75`.

**(c)** For every `t` with half-count prime to 3 the limit `rho(t)` of the fraction of odd numbers below `X` whose
signed descent passes through `t` exists, and `t rho(t) -> F(n)` along 3-adic classes, where `F` solves
`F(n) = F(2n)/2 + (3/4) F(r(n))` on the 3-adic units with mean `3 / (ln 2 (1/2 + (1/2) log_2 (4/3)))`.
