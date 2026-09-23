---
id: HYP-9130
title: "2-adic zero estimate for cube Pade lattices (implies HYP-9127)"
status: >
  OPEN; FINITE-EXACT support. For some w in (1.4265, 2) and infinitely many
  E, the lattice {A in Z^(E+1) : A Theta_3 = B mod x^(ceil(wE)), deg B <= E}
  contains a vector of height 2^(o(E)) whose first remainder coefficient is
  not divisible by 2^10. Here Theta_3(x) = sum x^(k^3). By Corollary F3 of
  the cube-theta note, this implies HYP-9127 (irrationality of
  sum (2^10/3^9)^(k^3) in Q_2). Evidence: true for every tested E <= 1000
  (diagonal forms, v_2 of the leading remainder <= 3) and E <= 300 (Siegel
  forms, odd leading remainder). A purely combinatorial statement about
  integer Pade-type forms of the 0/1 cube series.
source: collatz-procgen-20260923 cube-theta lane (candidate HC-CT1)
depends_on:
  - 05-knowledge/results/collatz_procgen_20260923_cube_theta.md
  - 05-knowledge/hypotheses/HYP-9127-cube-swap-cubic-theta.md
---

# HYP-9130 -- the combinatorial route to the cube-swap number

**Refutation form:** for every `w > 1.4265` the lattices eventually contain
no low-height vector with a small-2-power leading remainder.

**Companion conditional (PROVED reduction, Theorem U):** a uniform S-unit gap
`U(1/3)` also implies HYP-9127. The n-term abc conjecture does not.
