#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
LRC(14) ANGLE B -- explicit wide/dissociated closure for the L_y crux (HYP-2607).

Goal: show primitive E (gcd=1) with wide / dissociated structure has L_y(E) < cap_k,
with EXPLICIT constants, to combine with the bounded-spread finite check.

Three pieces:
  Step 1  L_y^inf = sum_r y_r C(6,r) (1-r/7)^{k-1}   (FULLY-INDEPENDENT limit; exponent k-1
          because e=0 is pinned in sector 0 and the OTHER k-1 generators are the free uniforms).
          Result (exact Fraction):
             k=8 :  40573/823543   ~ 0.049266
             k=9 :  2616535/17294403 ~ 0.151294
             k=10:  8830611/40353607 ~ 0.218831
          All << cap_k (margins ~0.33).
  Step 2  Fourier / BV stranger-peel error bound:
             J(A,E) = (1-r/7)^{k-1} + sum over NONTRIVIAL integer relations  sum_{e!=0} m_e e = 0
                                        of  prod_e c_{m_e}(A),
          where c_m(A) is the Fourier coeff of the avoid-indicator of |A|=r sectors,
             |c_m(A)| <= r|sin(pi m/7)|/(pi|m|) <= r/(pi|m|),  c_m=0 if 7|m,  c_0 = 1 - r/7.
          Single-peel form (E = E' + {w}, V = total variation of F=prod_{E'} avoid-indicator):
             |J(A,E) - (1-r/7) J(A,E')| <= r*V/(6 w),   V <= 14*sum(E').
          Propagation to L_y: |L_y(E) - L_inf| <= sum_{r>=1} |y_r| C(6,r) eps_r.
  Step 3  Reconciliation with scale-invariance (THM-531) and the partition:
          dilated APs are NOT primitive (gcd=D); dividing by gcd returns consec (bounded spread),
          so they live in the FINITE CHECK. Genuinely-wide primitive sets cap out empirically at
          L_y ~ 0.18 (k=8), 0.28 (k=9), 0.38 (k=10), all with margin >= 0.20 to cap.

NOTE (honest gap): the dissociation peel rigorously covers sets with a far dissociated stranger,
driving L_y -> L_inf << cap. It does NOT by itself cover wide-but-STRUCTURED primitive sets
(partial APs whose every element sits in a small relation). Those are handled empirically here
(<cap with large margin) but a self-contained proof for that pocket is still open.

kind-pasteur-2026-06-19-S12.
"""
import sys
from fractions import Fraction as F
from math import comb, sin, pi
sys.path.insert(0, '.')

CAPS = {8: F(3815, 10000), 9: F(49426, 100000), 10: F(6044, 10000)}


def g_vals(k):
    if k == 8:
        return [F((t-1)*(t-2)*(t-4)*(t-5), 40) for t in range(7)]
    return [F(-(t-2)*(t-3)*(t-6), 36) for t in range(7)]


def y_coeffs(k):
    """y_r from g(t)=sum_r y_r C(t,r) via inverse binomial (forward difference at 0)."""
    g = g_vals(k)
    return [sum((-1)**(r-t) * comb(r, t) * g[t] for t in range(r+1)) for r in range(7)]


def L_inf(k):
    y = y_coeffs(k)
    return sum(y[r] * comb(6, r) * (1 - F(r, 7))**(k-1) for r in range(7))


def error_budget(k):
    """Per-A error eps allowed so that |L_y - L_inf| < cap - L_inf."""
    y = y_coeffs(k)
    Li = L_inf(k)
    slack = CAPS[k] - Li
    coef = sum(abs(y[r]) * comb(6, r) for r in range(1, 7))
    return Li, slack, coef, slack / coef


def peel_bound(r, sumEprime, w):
    """Rigorous single-peel error: |J(A,E)-(1-r/7)J(A,E')| <= r*V/(6 w), V<=14*sum(E')."""
    V = 14 * sumEprime
    return F(r * V, 6 * w)


if __name__ == "__main__":
    print("=== LRC(14) Angle B: dissociation closure schema ===\n")
    print("Step 1 -- fully-independent limit  L_inf = sum_r y_r C(6,r)(1-r/7)^(k-1):")
    for k in (8, 9, 10):
        Li = L_inf(k)
        print(f"  k={k}: L_inf = {Li} = {float(Li):.6f}   cap = {float(CAPS[k]):.5f}   "
              f"margin = {float(CAPS[k]-Li):.6f}")
    print("\nStep 2 -- L_y error budget (need per-A error eps below this for dissociated branch):")
    for k in (8, 9, 10):
        Li, slack, coef, eps = error_budget(k)
        print(f"  k={k}: slack={float(slack):.5f}  sum|y_r|C(6,r)(r>=1)={float(coef):.3f}  "
              f"=> eps < {float(eps):.5f}")
        print(f"        peel gives eps <= r*14*sum(E')/(6 w); for r<=6, sum(E')~S, "
              f"need w > {float(6*6*14/(6*eps)):.0f}*S/?  (separation grows ~S/eps)")
    print("\nFourier coeff bound used: |c_m(A)| <= r|sin(pi m/7)|/(pi|m|), c_m=0 if 7|m. "
          "Per-m magnitudes (single sector):")
    print("  ", [round(sin(pi*m/7)/(pi*m), 5) if m % 7 else 0.0 for m in range(1, 8)])
