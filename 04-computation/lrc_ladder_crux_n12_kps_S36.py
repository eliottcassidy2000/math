#!/usr/bin/env python3
"""
kind-pasteur-2026-07-06-S36 (crux, mechanism-focused): where does each 11-speed base's
resonance-ladder rung land in opus's (s,k) form at n=12, and does any hit the gap
interior (1/13,2/25) with k>=2 & k<s<2k?

Emptiness of the gap is ALREADY well-evidenced (mac-mini census to height 48, 511k
survivors cleared; kps S33 ~377k families, 0 in gap).  This driver instead exhibits the
MECHANISM at n=12 on representative bases: the ladder M=mu*x/(x+rho) crosses the gap, and
we read off the crossing rung's opus-order k.

FAST design: exact M for the base; exact M at the gap-CROSSING rungs only (found via the
closed-form screen mu*x/(x+rho)); rho extracted from one exact rung.  No 15-rung sweep.
"""
from fractions import Fraction
from math import gcd, ceil
from functools import reduce
import numpy as np

N = 12
LO, HI = Fraction(1,N+1), Fraction(2,2*N+1)     # (1/13, 2/25)

def Mw(v):
    v=[x for x in v if x]; S=sum(abs(x) for x in v); Q=min(4*S,2*max(abs(x) for x in v)+2)
    va=np.array(v,dtype=np.int64); bn,bd,bc=0,1,1
    for q in range(2,Q+1):
        for c in range(1,q):
            r=(va*c)%q; d=np.minimum(r,q-r); bq=int(d.min())
            if bq*bd>bn*q: bn,bd,bc=bq,q,c
    return Fraction(bn,bd),(bc,bd)

def spectrum(M):
    p,q=M.numerator,M.denominator; return p, q-N*p     # M = s/(N s + k)
def order_ok(M):
    if not (LO<M<HI): return None
    s,k=spectrum(M)
    return (s,k,(k>=2 and k<s<2*k))

def analyze(base, label):
    if reduce(gcd,base)!=1 or len(set(base))!=11:
        print(f"  {label}: (skip: not 11 primitive speeds)"); return None
    mu,(c,D)=Mw(base)
    tag = "" if mu>HI else "  (plateau <= gap top; ladder cannot descend through gap)"
    # rho from first genuine-outlier resonance
    mx=max(base); j0=ceil((2*mx)/D)+1; x0=j0*D
    M0,_=Mw(sorted(base+[x0])); rho=x0*(mu-M0)/M0 if M0>0 else Fraction(0)
    print(f"  {label}: base={base}", flush=True)
    print(f"     mu=M(B)={mu} (~{float(mu):.4f}) at t={c}/{D};  rho={rho} (~{float(rho):.2f}){tag}", flush=True)
    if mu<=HI:
        return None
    # closed-form screen: rungs x=jD with mu*x/(x+rho) in a neighborhood of the gap
    # solve jD in (LO*rho/(mu-LO), HI*rho/(mu-HI)); widen x2 by 1.5 for rho-drift safety
    x_lo = float(LO*rho/(mu-LO)); x_hi = float(HI*rho/(mu-HI)) if mu>HI else 1e9
    j_lo=max(1, int(x_lo/D)-1); j_hi=int(x_hi/D*1.5)+2
    hits=[]; crossing=[]
    for j in range(j_lo, min(j_hi, j_lo+40)+1):
        x=j*D; v=sorted(base+[x])
        if len(set(v))!=12 or reduce(gcd,v)!=1: continue
        M,_=Mw(v)                                  # EXACT rung
        oo=order_ok(M)
        if oo is not None:
            s,k,good=oo
            crossing.append((j,x,M,s,k,good))
            if good: hits.append((j,x,M,s,k))
    if crossing:
        for j,x,M,s,k,good in crossing:
            flagtxt = "  <== IN-GAP, k>=2, k<s<2k !!" if good else ("  in gap (s,k)" if LO<M<HI else "")
            print(f"       j={j:>2} x={x:>4}  M={str(M):>9} (~{float(M):.4f})  (s,k)=({s},{k}){flagtxt}", flush=True)
    else:
        print(f"       (no exact rung landed in the open gap in screened range)", flush=True)
    return hits

print(f"=== S36 n=12 crux mechanism: 11-speed ladders vs gap (1/13,2/25); opus (s,k), k>=2 & k<s<2k required ===\n", flush=True)

CASES = []
# pure APs (expect k=1, edge-threading)
CASES.append(("AP{1..11}", list(range(1,12))))
CASES.append(("AP{1..10}+{12}", list(range(1,11))+[12]))
# S35-analog: AP{1..b} + near defect (like the {1..5}+7 that hit the mediant at n=7)
for d in [12,13,14,15,17,20,24]:
    CASES.append((f"AP{{1..10}}+{{{d}}}", list(range(1,11))+[d]))
CASES.append(("AP{1..9}+{11,13}", list(range(1,10))+[11,13]))
CASES.append(("AP{1..9}+{11,20}", list(range(1,10))+[11,20]))
CASES.append(("AP{1..8}+{10,12,14}", list(range(1,9))+[10,12,14]))
# dilated AP + defect
CASES.append(("2*AP{1..10}+{1}", [1]+[2*i for i in range(1,11)]))

allhits=[]
for label,base in CASES:
    h=analyze(base,label)
    if h: allhits+=[(label,)+t for t in h]
    print(flush=True)

print("="*72, flush=True)
if not allhits:
    print("RESULT: on these representative 11-speed bases, NO ladder rung lands in the", flush=True)
    print("n=12 gap with opus-order k>=2 & k<s<2k.  Pure/near-AP bases thread the edges", flush=True)
    print("(k=1) or overshoot; the interior mediant seat (s,k)=(3,2)=3/38 is not attained.", flush=True)
    print("MECHANISM evidence for (G), consistent with mac-mini census + opus finite-shapes.", flush=True)
else:
    print(f"RESULT: {len(allhits)} interior hit(s) -- INSPECT:", flush=True)
    for h in allhits: print("   ", h, flush=True)
