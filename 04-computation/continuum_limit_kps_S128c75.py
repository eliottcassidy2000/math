#!/usr/bin/env python3
"""continuum_limit_kps_S128c75.py -- kind-pasteur S128 cont.75.
THE CONTINUUM LIMIT F_inf(u), settling the width growth of THM-1146(IV)/(V).

DERIVATION.  The k1-gap at index j is [j/k1 + 1/(14k1), (j+1)/k1 - 1/(14k1)], of length
G = 6/(7k1).  Normalise it to [0,1] by dividing by G.  For k_i = k1 + d_i:
  * a tooth of k_i has half-width 1/(14 k_i), normalised k1/(12 k_i) -> 1/12
    (so full width -> 1/6);
  * its centre is at ((k1 s - j d_i)/(k1 k_i) - 1/(14 k1)) * 7k1/6  ->  (7/6) g_i - 1/12,
    where g_i = frac(-d_i u) and u = j/k1;
  * the threshold 1/(7 k4) normalises to k1/(6 k4) -> 1/6.
So in the limit: THREE teeth of width 1/6 at centres (7/6)*frac(-d_i u) - 1/12 inside [0,1],
and BAD means the longest surviving piece is <= 1/6.  One parameter, u in [0,1].
Compute the bad set exactly enough to settle whether the per-run width converges.
PRINT DATA ONLY."""
import sys
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
DS=(1,2,3)
def Finf(u):
    """longest surviving piece in the limit, u a Fraction"""
    cuts=[]
    for d in DS:
        g=(-d*u)%1
        c=F(7,6)*g-F(1,12)
        a=c-F(1,12); b=c+F(1,12)
        a=max(a,F(0)); b=min(b,F(1))
        if b>a: cuts.append((a,b))
    cuts.sort(); cur=F(0); L=F(0)
    for a,b in cuts:
        if a>cur and a-cur>L: L=a-cur
        cur=max(cur,b)
    if 1-cur>L: L=1-cur
    return L
print("### F_inf sampled: where is it below the 1/6 threshold? ###")
N=20160   # divisible by lots; exact rational grid
thr=F(1,6)
bad=[]
for t in range(N):
    u=F(t,N)
    if Finf(u)<=thr: bad.append(t)
print("  grid points: %d ; bad points: %d ; bad measure ~ %.6f"%(N,len(bad),len(bad)/N))
runs=[]
if bad:
    cur=[bad[0]]
    for x in bad[1:]:
        if x==cur[-1]+1: cur.append(x)
        else: runs.append(cur); cur=[x]
    runs.append(cur)
    if len(runs)>1 and bad[0]==0 and bad[-1]==N-1:
        runs[0]=runs[-1]+runs[0]; runs.pop()
print("  runs: %d ; run widths (as fractions of the period): %s"%(
    len(runs),[round(len(r)/N,6) for r in runs]))
print("  run centres: %s"%[round((r[0]+r[-1])/2/N,5) for r in runs])
print()
print("### convergence check: finite k1 per-run fraction vs the limit ###")
print("  the finite values from THM-1146(IV) were:")
print("    k1 = 157  0.0382 | 207  0.0386 | 257  0.0428 | 307  0.0423 | 357  0.0448 | 407  0.0467")
if runs:
    per=max(len(r)/N for r in runs)
    print("  CONTINUUM per-run width  = %.6f"%per)
    print("  CONTINUUM total bad      = %.6f"%(len(bad)/N))
    print("  finite values are %s the limit"%("BELOW and rising toward" if per>0.0467 else "ABOVE -- limit already exceeded"))
print()
print("### the load-bearing comparison ###")
tot=len(bad)/N
print("  total bad measure (continuum)      = %.6f"%tot)
print("  S(P) minimum measure (atlas)       = 0.164")
print("  counting argument needs total < S  : %s"%("HOLDS" if tot<0.164 else "FAILS"))
print("  margin                             = %.6f"%(0.164-tot))
print()
print("### where exactly does F_inf dip?  profile around the minimum ###")
print("   u        F_inf(u)    vs 1/6")
for num in [0,1,2,3,4,5,6,7,8,9,10,11,12]:
    u=F(num,48)
    v=Finf(u)
    print("  %-8s %-11.6f %s"%(str(u),float(v),"BAD" if v<=thr else "ok"))
print("DONE")
