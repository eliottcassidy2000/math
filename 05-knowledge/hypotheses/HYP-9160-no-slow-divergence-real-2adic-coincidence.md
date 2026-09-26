---
id: HYP-9160
title: "No slow divergence: every divergent Syracuse orbit has summable reciprocals, equivalently the real and 2-adic values of the Bernstein series of a non-periodic word never coincide at a positive integer"
status: >
  OPEN HYPOTHESIS. Let m_0 = n, m_(l+1) = (3 m_l + b)/2^(v_(l+1)), b = +-1,
  be a Syracuse orbit on either sheet and of either sign, with halving word
  d_L = v_1 + ... + v_L. Proposition 6 of the companion note gives the exact
  identity m_L 2^(d_L)/3^L = n prod_(l<L) (1 + b/(3 m_l)) and the real
  value R(d) = sum_l 2^(d_l)/3^(l+1) = b n (prod_l (1 + b/(3 m_l)) - 1).
  Hypothesis: if the orbit diverges then sum_l 1/m_l < infinity. On the
  minus sheet (positive n, b = -1; equivalently negative n on the plus
  sheet) this says c = n prod (1 - 1/(3 m_l)) > 0, i.e. R(d) < n, so that a
  word d that is not eventually periodic never has R(d) = R_2(d) = n with
  R_2(d) the 2-adic value (Bernstein: n = -b R_2(d)). On the plus sheet it
  says every divergent positive orbit has R(d) < infinity. Implied by "no
  divergent orbit" on the relevant sheet and by nothing weaker that is
  known. Holds on words of bounded critical discrepancy (they carry no
  integer orbit at all: capacity plus ordered carry, in-house). A slow
  divergence m_l ~ l^alpha, alpha <= 1, has discrepancy of order log l,
  the first open regime.
source: collatz-squares-doubles-20260925 session (opus), squares/doubles transport foundry, REALVAL x VAL card
depends_on:
  - 05-knowledge/results/collatz_sqdbl_20260925_squares_doubles_foundry.md (Proposition 6)
related:
  - 05-knowledge/results/collatz_procgen_20260922_hard_class.md (Proposition T, two places)
  - 05-knowledge/results/collatz_mod6_20260922_counterexample_portrait.md (S4, Eliahou identity on cycles)
  - 05-knowledge/hypotheses/HYP-9127-cube-swap-word-cubic-theta.md (the smallest open Periodicity instance)
---

# HYP-9160 -- no slow divergence (real and 2-adic values never coincide off periodic words)

**Refutation form.** A divergent Syracuse orbit (on `3n+1` or `3n-1`, of
either sign) whose odd iterates satisfy `sum_l 1/m_l = infinity`; for example
polynomial growth `m_l ~ l^alpha` with `alpha <= 1`. On the minus sheet such
an orbit would be a word `d` with `R(d) = R_2(d) = n`: the real and the 2-adic
value of the same digit series, both equal to the starting integer.

**Why it matters.** The synthesis of the procedural session asks for a
mechanism that sees the sign and the drift together. Proposition 6 shows the
sign enters the divergence half exactly once, as the inequality `R(d) <= n`
on the minus sheet, with equality iff `sum 1/m_l = infinity`. Full-rate
divergence (`c > 0`) is invisible to the real place; slow divergence is the
only corner where the real place gives an equality. This hypothesis isolates
that corner. It is strictly weaker than the Periodicity Conjecture on the
sheet.

**Evidence.** FINITE-EXACT (probe P6): on `5n +- 1`, where divergence is
expected, the orbits of `n = 37, 47` (plus) and `n = 9, 11` (minus) have
`c_L` stabilised to twelve digits by `L = 100` (`37.374..., 47.839...,
8.528..., 10.660...`) and `sum 1/m_l` between `0.05` and `0.27`: full rate.
Every `3n-1` orbit of odd `n <= 2000` enters a cycle, with `R(d) = n` exactly
(closed form), as the hypothesis requires for eventually periodic words.

**Cheapest partial.** Extend the bounded-critical-discrepancy mechanism to
discrepancy `<= C log l`; that regime contains every orbit with
`m_l = O(l^alpha)`. No decisive finite test exists.

**Non-consequences.** Proving it would not prove Collatz or `3n-1`; it would
remove one of the two divergence types on each sheet.
