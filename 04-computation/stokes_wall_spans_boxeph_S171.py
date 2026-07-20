#!/usr/bin/env python3
"""stokes_wall_spans_boxeph_S171.py (HYP-8395) — code as run (S171 heredocs),
frozen out in 05-knowledge/results/. RESULTS: (1) GLOBAL RIGIDITY upgrade:
A(s) = E[e^{-sP}] is an exponential period, hence HOLONOMIC; Watson's sectorial
A == 1 propagates along the D-module: nullcone members satisfy E[e^{-sP}] == 1
as an EXACT identity — Stokes-faked flatness cannot be asymptotic fakery.
(2) The wall reduces to PICARD-LEFSCHETZ: exactness across sectors forces all
Stokes jumps (connection coefficients / vanishing-cycle pairings) to vanish;
holomorphic analogue: multipliers = intersection numbers, nonzero when critical
values couple. Target: 'some Stokes multiplier nonzero for every mixed P'.
(3) SPAN {+2,+1,-1} (constants, gauged to one unknown): E2 == 4 CONSTANT:
nullcone EMPTY outright. (4) SPAN {+2,+1,-1,-2} (constants, two unknowns):
E2, E3 leave (a,d) = (±1/2, ∓1/2) — the factorized P = 2X(1+iY), whose Z/2
real-parity grading fakes U(1) through m = 3 — E4 kills both: SPAN CLOSED.
Structural find: real-parity (Z/2) is the only faking mechanism at low m."""
print(open("05-knowledge/results/stokes_wall_spans_boxeph_S171.out").read())
