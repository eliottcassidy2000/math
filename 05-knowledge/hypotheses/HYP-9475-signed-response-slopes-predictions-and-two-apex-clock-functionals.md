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

## (B) Two-apex drift functionals (OPEN)

For a block d of the two-apex grammar G(2,0), Session 194 proved H(u d^k) = A + Jk + Kk^2 + (D + Ek) 2^(Lk) + G 3^(Lk)
with A, J, K depending on the prefix u. The eigenvalue-1 generalized left eigenspace of M_d therefore carries drift
functionals of degree 1 and 2. Conjecture: for every mixed block d of length <= 6 the degree-1 functional (the
coefficient J as a linear functional of the 17-dimensional prefix profile) is a signed response in the sense of
THM-4882 (3): onto the integers of its parity on strong two-apex words, by a modular seed theorem for two-apex
profiles. Test: compute the left Jordan chains of M_d from the G(2,0) automaton dump of THM-4878 and scan prefixes.

## (C) Witness size (OPEN)

Session 198's recipe reaches a target r with h = O(|r|) repetitions of a period that multiplies U by 2^k, so the
witness one-apex tournament has about k h chain vertices and U with that many binary digits. Conjecture: for two
slopes with coprime drops K, K' (e.g. 2 and 4 from 2/3 and 6/7 do not qualify; 2 and 18 from 4/5 and 8/11 do not
either; coprime pairs exist for other slopes, e.g. drops 2 and 11 would), alternating the two clocks reaches every
target with O(log |r|) periods, i.e. a witness of order O(log |r|), matching the lower bound of S79's counting lemma
in spirit (a signed, not a positive, count). Untested.
