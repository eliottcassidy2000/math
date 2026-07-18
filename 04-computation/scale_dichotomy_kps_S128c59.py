#!/usr/bin/env python3
"""scale_dichotomy_kps_S128c59.py -- kind-pasteur S128 cont.59.
THE DICHOTOMY that replaces the refuted non-covering lemma.
Defeating every modulus in [15,Q] forces the killers to be HUGE.  Huge killers are
SCALE-SEPARATED from a core inside {1..12}, and then the two blocks decouple:
   t chosen in the core-safe set S ; u = k1*t mod 1 equidistributes inside S ;
   so M -> min( M(core) , max_u min(||u||,||2u||) ).
Test the convergence on {2..12} u {V, 2V} as V grows.  PRINT DATA ONLY."""
import sys
from fractions import Fraction as F
from math import gcd
sys.stdout.reconfigure(line_buffering=True)
def nd(x):
    fx=x-int(x); fx=fx+1 if fx<0 else fx; return min(fx,1-fx)
def Mexact(V,qhi):
    best=(F(0),None)
    for q in range(2,qhi+1):
        for a in range(1,q//2+1):
            if gcd(a,q)!=1: continue
            t=F(a,q); mn=min(nd(v*t) for v in V)
            if mn>best[0]: best=(mn,t)
    return best
CORE=list(range(2,13))
print("### the two block measures the dichotomy predicts ###")
mc,tc=Mexact(CORE,4000)
print("  M(core {2..12})            = %-10s (%.5f)  at t=%s"%(mc,float(mc),tc))
mk,tk=Mexact([1,2],4000)
print("  max_u min(||u||,||2u||)    = %-10s (%.5f)  at u=%s"%(mk,float(mk),tk))
pred=min(mc,mk)
print("  PREDICTED limit min(.,.)   = %-10s (%.5f)   vs threshold 1/14 = %.5f"%(pred,float(pred),1/14))
print()
print("### convergence: M({2..12} u {V,2V}) as V grows ###")
print("   V      M exact          value     vs 1/14   vs predicted limit")
for V in [157,200,260,400,650,1000,1600]:
    Vs=CORE+[V,2*V]
    if len(set(Vs))!=13: continue
    m,t=Mexact(Vs,max(1200,3*V))
    print("  %-6d %-16s %.5f   %-8s %s"%(V,m,float(m),"OK" if m>=F(1,14) else "BELOW",
        "= limit" if m==pred else ("%.3f x limit"%(float(m)/float(pred)))))
print()
print("### the adversarial families themselves: are they scale-separated enough? ###")
print("  the criterion-defeating pairs found at Q=18,20,22 were (182,271),(161,1274),(881,1274)")
for k1,k2 in [(182,271),(161,1274),(881,1274)]:
    Vs=CORE+[k1,k2]
    if len(set(Vs))!=13: continue
    m,t=Mexact(Vs,4000)
    print("  K=(%-4d,%-4d) ratio=%.2f  M=%-12s (%.5f)  >=1/14: %s"%(
        k1,k2,max(k1,k2)/min(k1,k2),m,float(m),m>=F(1,14)))
print()
print("### how fast does the min defeating-killer size grow?  extend the search ###")
def la(r,q):
    r%=q; return min(r,q-r)
def killing(P,k1,k2,q):
    thr=-(-q//14)
    for a in range(1,q):
        if all(la(p*a,q)>=thr for p in P) and la(k1*a,q)>=thr and la(k2*a,q)>=thr: return False
    return True
def legal(P,k1,k2):
    Vv=sorted(P+[k1,k2])
    if len(set(Vv))!=13: return False
    if not all(any(v%q==0 for v in Vv) for q in range(2,15)): return False
    for i,v in enumerate(Vv):
        if not any(j!=i and v<=13*Vv[j] for j in range(len(Vv))): return False
    return any(a>13*b for a in Vv for b in Vv)
print("  Q    min max(k1,k2) defeating [15,Q]")
for Q in [16,17,18,19,20,21,22,23,24]:
    best=None
    for kmax in range(157,12000):
        found=False
        for other in range(157,kmax+1):
            k1,k2=other,kmax
            if k1==k2: continue
            if not legal(CORE,k1,k2): continue
            if all(killing(CORE,k1,k2,q) for q in range(15,Q+1)): found=True; break
        if found: best=(kmax,other); break
    print("  %-4d %s"%(Q,("%d  (pair %d,%d)"%(best[0],best[1],best[0])) if best else "> 12000"))
print("DONE")
