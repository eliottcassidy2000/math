#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ADVERSARIAL VERIFICATION of the three "outside-the-box" LRC(14) closers:
  C1 Selberg-Beurling bandlimited per-factor majorant
  C2 Banaszczyk transference / 6-fold L^2 (Cauchy-Schwarz) convergence
  C3 Lam-Leung visibility / dilation reduction of unbounded->bounded

All numeric work in exact Fractions where the quantity is rational; the Fourier
kernel K(n) is a finite signed sum of products of EXACT rational ctilde values
(ctilde_T(0)=1-|T|/7) times sine-kernel pieces, so per-relation K(n) is computed
in exact closed form using the cyclotomic field Q(zeta_7) -> real via cos(2pi/7).

We test the SPECIFIC quantitative claims of each candidate, not vague intuition.
"""
import sys, itertools, math
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True) if hasattr(sys.stdout,'reconfigure') else None

# ---------- exact geometric meas(S7) (ground truth, all rational) ----------
def measS7(E):
    E=sorted(set(E)); bps=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(0,7*abs(e)+1): bps.add(F(m,7*e))
    bps=sorted(b for b in bps if 0<=b<=1); total=F(0)
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2; secs=set()
        for e in E:
            y=(e*xm)%1; secs.add((y.numerator*7)//y.denominator)
        if len(secs)==7: total+=x1-x0
    return total

def M7(k):
    """K(0): iid limit, exact rational."""
    s=F(0)
    for t in range(0,7):
        s += F((-1)**t * math.comb(6,t)) * F(7-t,7)**(k-1)
    return s

# ---------- exact per-coordinate sine kernel |s_T(n)| = |sin(pi n/7)|/(pi|n|) ----------
# We do NOT need the exact complex K(n) for the convergence tests below; we need the
# ENVELOPE behaviour and the actual SIGNED correction (= measS7 - M7), both available.

if __name__=="__main__":
    print("="*80)
    print("STEP 0: reconcile conventions — exact meas(S7), M7(k), correction delta_k at AP")
    print("="*80)
    caps={8:F(2243,5880), 9:None, 10:None}
    cap_float={8:0.38153,9:0.49426,10:0.6044}
    for k in [8,9,10]:
        E=list(range(k))
        g=measS7(E); m7=M7(k); delta=g-m7
        print(f"  k={k} AP: meas(S7)={g}={float(g):.5f}  M7={m7}={float(m7):.5f}  "
              f"delta(correction)={float(delta):.5f}")
        capk = caps[k] if caps[k] is not None else cap_float[k]
        print(f"        cap_{k}={float(capk):.5f}  margin meas->cap = {float(capk)-float(g):.5f}  "
              f"margin delta->(cap-M7) = {float(capk)-float(m7)-float(delta):.5f}")
    print()
