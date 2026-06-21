#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_routeB_Ly_vs_U4_order_opus_0620s2.py   (opus-2026-06-20-S2)

Sharpen Route B.  TWO functionals both bound p_0 (the lonely measure):
  L_y (THM-534 dual): weights phi^Ly; the OPERATIVE certificate that CLOSES k=8,9,10
                      (L_y(consec)<=cap).  k=8: phi=(1,0,0,1/10,0,0,1).
  U4  (THM-556 Bonf): weights phi^U4=(1,0,0,0,0,1,5); a LOOSER bound, U4 does NOT
                      close (U4(consec_8)=0.480 > cap_8=0.381).
So the REAL extremality target = L_y-extremality (the one that closes), and Route B
must be run on L_y, not U4.

This script, for BOTH functionals and k=8,9,10:
  (1) confirm convexity of the weight vector phi (2nd differences).
  (2) confirm consec is the bank-max.
  (3) THE CENTRAL TEST: is phi a NONNEGATIVE combination of the SHIFTED-CONVEX-CORNER
      basis  c_a(t) = (t-a)_+  PLUS the reflected corner  d_b(t)=(b-t)_+  ?
      i.e. phi(t) = alpha + beta*t + sum_a w_a (t-a)_+ + sum_b v_b (a... )  -- find the
      decomposition.  If phi = const + linear + (nonneg)*upper-corners + (nonneg)*lower-corners,
      then U4(E) = const + linear*E[N] + sum_a w_a E[(N-a)_+] + sum_b v_b E[(b-N)_+].
      We then test, CUT BY CUT, whether consec maximizes each *active* corner functional.
  (4) Since E[N] is NOT constant, find the EXACT linear-coefficient and see if the
      linear part fights or helps (consec has near-min E[N]; if linear coeff>0 it hurts).
"""
import sys, itertools
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)

def dist_p(E):
    E = sorted(set(E))
    bps = set([F(0), F(1)])
    for e in E:
        if e == 0: continue
        for a in range(0, 7*e+1):
            bps.add(F(a, 7*e))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    p = [F(0)]*7
    for i in range(len(bps)-1):
        lo, hi = bps[i], bps[i+1]
        if hi == lo: continue
        mid = (lo+hi)/2
        hit = set()
        for e in E:
            v = (e*mid) % 1
            hit.add((v.numerator*7)//v.denominator)
        t = sum(1 for j in range(1,7) if j not in hit)
        p[t] += (hi-lo)
    return p

def phi_Ly(k):
    if k==8: return [F((t-1)*(t-2)*(t-4)*(t-5),40) for t in range(7)]
    if k in (9,10): return [F(-(t-2)*(t-3)*(t-6),36) for t in range(7)]
    return [F((t-3)*(t-4),12) for t in range(7)]
PHI_U4 = [F(1),F(0),F(0),F(0),F(0),F(1),F(5)]

def apply_phi(p, phi): return sum(p[t]*phi[t] for t in range(7))
def EN(p): return sum(t*p[t] for t in range(7))

def upper_corner(p,a): return sum((t-a)*p[t] for t in range(a+1,7))   # E[(N-a)_+]
def lower_corner(p,b): return sum((b-t)*p[t] for t in range(0,b))     # E[(b-N)_+]

def decompose_convex(phi):
    """Write phi(t)=alpha+beta*t+sum_{a=1..5} w_a (t-a)_+  on t=0..6.
       For a convex sequence this gives w_a = 2nd-difference at a >=0 (and is unique
       given alpha=phi(0), beta=phi(1)-phi(0))."""
    alpha=phi[0]; beta=phi[1]-phi[0]
    # w_a = Delta^2 phi(a) for a=1..5
    w={a: phi[a-1]-2*phi[a]+phi[a+1] for a in range(1,6)}
    # verify reconstruction
    def rec(t):
        s=alpha+beta*t
        for a in range(1,6):
            if t>a: s+=w[a]*(t-a)
        return s
    ok=all(rec(t)==phi[t] for t in range(7))
    return alpha,beta,w,ok

def consec(k): return list(range(k))

if __name__=="__main__":
    for name,getphi in [("L_y",phi_Ly),("U4",lambda k:PHI_U4)]:
        print("="*78); print(f"FUNCTIONAL {name}"); print("="*78)
        for k in (8,9,10):
            phi=getphi(k)
            d2=[phi[t-1]-2*phi[t]+phi[t+1] for t in range(1,6)]
            alpha,beta,w,ok=decompose_convex(phi)
            print(f"\n--- k={k}  phi={[str(x) for x in phi]} ---")
            print(f"  2nd diffs (t=1..5) = {[str(x) for x in d2]}  "
                  f"{'CONVEX' if all(x>=0 for x in d2) else 'NOT convex'}")
            print(f"  decomposition: phi(t)= {alpha} + ({beta})*t + sum_a w_a (t-a)_+   exact={ok}")
            print(f"     linear coeff beta = {beta} ({'>0 hurts consec (low mean)' if beta>0 else '<=0 helps' })")
            print(f"     corner weights w_a (a=1..5) = {{ {', '.join(f'{a}:{w[a]}' for a in range(1,6))} }}")
            print(f"     => active upper-corners: {[a for a in range(1,6) if w[a]!=0]}")

            C=consec(k); pc=dist_p(C); Uc=apply_phi(pc,phi)
            cornersC={a:upper_corner(pc,a) for a in range(1,6)}
            # bank
            span=k+4 if k==8 else k+3
            bank=[[0]+list(r) for r in itertools.combinations(range(1,span+1),k-1)]
            beats=0; worst=None
            # cut-by-cut: for each active corner a (w_a>0), does consec maximize E[(N-a)_+]?
            cornermax={a:0 for a in range(1,6) if w[a]!=0}
            beta_hurts_but_total_wins=0
            for E in bank:
                p=dist_p(E); U=apply_phi(p,phi)
                if U>Uc+F(1,10**15): beats+=1
                for a in cornermax:
                    if upper_corner(p,a)>cornersC[a]+F(1,10**15): cornermax[a]+=1
            print(f"  consec U={float(Uc):.6f}, beats={beats}/{len(bank)} "
                  f"(0 = consec is bank-max)")
            print(f"  consec corners E[(N-a)_+]: "
                  f"{ {a:float(cornersC[a]) for a in range(1,6) if w[a]!=0} }")
            print(f"  # shapes EXCEEDING consec on each active corner a: {cornermax}")
            print(f"     => if all zero, consec maximizes every active corner -> "
                  f"with beta and the corners, decompose the win.")
