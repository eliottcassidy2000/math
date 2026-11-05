---
id: HYP-9475
title: "Signed one-apex responses. (A) Pre-registered P1-P6 (commits 0e5ba011be, c6d7e35dc4): P1 small-witness coverage of Session 198's class; P2 both-odd slopes give every even integer; P3 even denominators below 3/4 confine positive values to one class modulo 4; P4 above 3/4 the parity hypothesis is unnecessary; P5 one sign missing outside (1/2, 1); P6 the cyclic-corrected clock K_cyc = K_lin - b. Outcomes: P1 supported, P2 CONFIRMED then PROVED (THM-4882 (3)(ii)), P3 CONFIRMED then PROVED, P4 REFUTED (the negative side is confined), P5 CONFIRMED, P6 CONFIRMED. (B) Conjecture: two-apex blocks carry two independent drift functionals (degree-1 and degree-2 counters of Session 194), and for a mixed block with a nonzero degree-1 drift the degree-1 functional attains every integer of its parity on strong two-apex words. (C) Conjecture: mixing two clocks with coprime drops bounds the witness size of Session 198's recipe by a polynomial in the target."
status: (A) P1 supported (small-witness), P2 PROVED, P3 PROVED, P4 REFUTED, P5 CONFIRMED, P6 CONFIRMED (all outcomes recorded in the source note, section 5); (B), (C) OPEN, untested
session: opus-2026-11-05-S80
date: 2026-11-05
related: [THM-4882, THM-4878, THM-4876, HYP-9468, HYP-9470, session198/rational-response-coverage/PROOF.md, session194/two-apex/PROOF.md]
script: 04-computation/experiments/signed_slope_probe_20261105.cpp; 04-computation/experiments/signed_clock_recipe_20261105.py; 04-computation/experiments/signed_slope_prereg_20261105.md
---

# HYP-9475 - signed responses: pre-registered predictions and the two-apex clock conjectures

## (A) Pre-registered predictions and outcomes

See `signed_slope_prereg_20261105.md` (P1-P5 written before the probe, P6 before the recipe run) and the source note
`05-knowledge/results/signed_responses_slope_classification_clock_barrier_20261105.md`, section 5. P4 ("above 3/4
every odd target is attained for a odd, b even") was the one guess made against the structure, with low confidence;
the data refuted it (the negative side is confined to one class modulo 4) and the proof of THM-4882 (3)(iii) explains
why: the even-m branch is eventually positive above 3/4, so negative values come only from the odd-m branch, which
has a fixed residue modulo 4.

## (B) Two-apex drift functionals (FINITE-EXACT data; the naive form withdrawn, the lattice form OPEN)

For a block d of the two-apex grammar G(2,0), Session 194 proved H(u d^k) = A + Jk + Kk^2 + (D + Ek) 2^(Lk) + G 3^(Lk)
with A, J, K depending on the prefix u. The first draft conjectured that the degree-1 drift J(u) is onto the integers
of a parity. The test (`two_apex_drift_functionals_20261105.py` on the 200-state automaton dump of THM-4878: the six
coefficients solved exactly from k = 3..8 for all 21,845 prefixes of length <= 7, the form verified at k = 9, 10 with
no failure) shows the coefficients are rational, not integral: for d = 12, K = 0 for every prefix and J(u) = (6t + 4)/9
for an integer t(u), with 192 distinct values in [-758/9, 556/9] and every value of that progression attained in the
window [-56/9, 58/9]; for d = 120, K = 13/49 and J has denominator 4459 (2,439 distinct values, none in [-60, 60]/4459);
for d = 0312, K = 8/15 and J has denominator 375 with numerators 2 mod 4 (675 values). Appending the block to the
prefix leaves J unchanged and shifts the constant term A by J, so the two-apex analogue of Session 198's signed
response is the pair (A(u), J(u)): the constant term is the programmable response and the drift is its clock
constant, which here depends on the seed. Revised conjecture: for d = 12, {J(u)} = {(6t + 4)/9 : t in an interval
growing linearly with |u|} and the constant terms A(u) over prefixes with a fixed J form a full arithmetic
progression of difference J; a modular seed theorem for two-apex profiles would then give exact coverage of a signed
response of G(2,0) with the ternary radix present. Untested beyond the data above.

## (C) Witness size (OPEN)

Session 198's recipe reaches a target r with h = O(|r|) repetitions of a period that multiplies U by 2^k, so the
witness one-apex tournament has about k h chain vertices and U with that many binary digits. Conjecture: for two
slopes with coprime drops K, K' (e.g. 2 and 4 from 2/3 and 6/7 do not qualify; 2 and 18 from 4/5 and 8/11 do not
either; coprime pairs exist for other slopes, e.g. drops 2 and 11 would), alternating the two clocks reaches every
target with O(log |r|) periods, i.e. a witness of order O(log |r|), matching the lower bound of S79's counting lemma
in spirit (a signed, not a positive, count). Untested.
