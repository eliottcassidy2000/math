#!/usr/bin/env python3
"""exact_width_kps_S128c75.py -- kind-pasteur S128 cont.75 (part 2).
CONFIRMING THE EXACT WIDTH 1/21.
The continuum profile is linear on [0,1/4]: the d=3 tooth has centre 13/12 - (7/2)u, so its
LEFT edge is 1 - (7/2)u, and the free piece [0, 1-(7/2)u] gives
        F_inf(u) = 1 - (7/2) u        on [0,1/4].
Threshold 1/6 is crossed at 1 - (7/2)u = 1/6, i.e. u = 5/21.
Grid measurement put the run end at ~0.28576, and 6/21 = 2/7 = 0.285714.  So the per-run
width should be EXACTLY 1/21, total bad 2/21.  Verify against a fine exact grid."""
import sys
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
DS=(1,2,3)
def Finf(u):
    cuts=[]
    for d in DS:
        g=(-d*u)%1
        c=F(7,6)*g-F(1,12)
        a=max(c-F(1,12),F(0)); b=min(c+F(1,12),F(1))
        if b>a: cuts.append((a,b))
    cuts.sort(); cur=F(0); L=F(0)
    for a,b in cuts:
        if a>cur and a-cur>L: L=a-cur
        cur=max(cur,b)
    if 1-cur>L: L=1-cur
    return L
print("### the linear branch: is F_inf(u) = 1 - (7/2)u on [0,1/4]? ###")
ok=True
for n in range(0,26):
    u=F(n,100)
    pred=1-F(7,2)*u
    got=Finf(u)
    if got!=pred: ok=False; print("   mismatch at u=%s: got %s want %s"%(u,got,pred))
print("  exact match on u = 0, 1/100, ..., 25/100 : %s"%ok)
print()
print("### the endpoints of the bad run ###")
print("  entry: 1 - (7/2)u = 1/6  =>  u = 5/21 = %s = %.7f"%(F(5,21),5/21))
print("  F_inf(5/21)      = %s"%Finf(F(5,21)))
print("  F_inf just below = %s (u=5/21 - 1/10000)"%Finf(F(5,21)-F(1,10000)))
print("  F_inf just above = %s (u=5/21 + 1/10000)"%Finf(F(5,21)+F(1,10000)))
print()
print("  exit candidate 2/7 = 6/21 = %.7f"%(2/7))
print("  F_inf(2/7)       = %s = %.7f"%(Finf(F(2,7)),float(Finf(F(2,7)))))
print("  F_inf just below = %.7f"%float(Finf(F(2,7)-F(1,10000))))
print("  F_inf just above = %.7f"%float(Finf(F(2,7)+F(1,10000))))
print()
print("### exact bad-set measure on a fine grid ###")
N=44100
bad=sum(1 for t in range(N) if Finf(F(t,N))<=F(1,6))
print("  grid N=%d : bad measure = %.7f"%(N,bad/N))
print("  2/21 = %.7f  ; difference = %.2e"%(2/21,abs(bad/N-2/21)))
print()
print("### VERDICT ###")
print("  per-run width  = 1/21 = %.7f  (CONVERGES; finite k1 values rise toward it)"%(1/21))
print("  total bad      = 2/21 = %.7f"%(2/21))
print("  S(P) minimum   = 0.164")
print("  counting needs total < S(P): %s   margin = %.5f"%(2/21<0.164,0.164-2/21))
print("DONE")
