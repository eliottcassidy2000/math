#!/usr/bin/env python3
"""
klein-2026-07-05-S134 (HYP-4099) - THE PEEL TOWER: rho-parametric window thresholds.

Generalized window lemma (output radius rho): base margin beta at one point + top above
  threshold(k) = rho * B / (beta - rho),   beta = 1/(k+1)  (citation floor at k runners)
=> margin rho for all k+1 runners at some t in the window.

Tower for the DICHOTOMY loose side (rho = 2/25, 12-runner target): peel the 12-base
top-down; each level either (a) top big -> window step, (b) repeat -> citation bump,
(c) small spread -> kps band_margin_reciprocal (margin a/(a+b) >= rho iff b <= a*(1-rho)/rho),
(d) the wedge -> finite/CRT checks (mac-mini). Termination: threshold(k) <= B ends the tower.

Print the threshold table for rho = 2/25 and rho = 1/14, the termination level, and the
a-priori LEAF SPREAD BOUND = product of thresholds down to termination -- the explicit
height bound for the finite-check domain.
"""
from fractions import Fraction as F

for rho in (F(2,25), F(1,14)):
    print(f"\n=== rho = {rho} ===")
    print(f"{'k (base size)':>14} {'beta=1/(k+1)':>13} {'threshold/B':>12} {'terminates?':>12}")
    prod = F(1)
    for k in range(12, 1, -1):
        beta = F(1, k+1)
        if beta <= rho:
            print(f"{k:>14} {str(beta):>13} {'(beta<=rho)':>12} {'CITE-DIRECT':>12}")
            break
        thr = rho / (beta - rho)   # threshold = thr * B
        term = thr <= 1
        print(f"{k:>14} {str(beta):>13} {float(thr):>12.3f} {('YES' if term else ''):>12}")
        if term:
            break
        prod *= thr
    print(f"leaf spread bound (product of thresholds above): {float(prod):.1f}")
    # band-loose criterion at this rho: spread b/a <= (1-rho)/rho
    print(f"kps band-loose spread bound (1-rho)/rho = {float((1-rho)/rho):.2f}")
