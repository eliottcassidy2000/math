#!/usr/bin/env python3
"""
lrc14_sector_perT_macmini_0618s7b.py  (mac-mini-2026-06-18-S7, ANGLE B re-examination)

CAREFUL re-test of per-T extremality.  s7b claimed "AP is IE-extremiser for every T-class"
but it only tested ONE representative T per rotation-class.  The blockword run showed
meas(B_T) is NOT rotation-invariant and the t=1 term has AP NOT the minimizer.  Resolve:

  meas(S7(E)) = sum_{M subset Z/7} (-1)^|M| meas(B_M(E)),  M = missed residue set.
  AP maximizes meas(S7) AGGREGATE (verified). Question: is the extremality TRUE per-M,
  or only after the IE sum?  Determine EXACTLY which M-blocks AP wins/loses, and whether
  the LOSSES are dominated by the WINS (a cancellation/majorization structure).

Also: corr(E) can be NEGATIVE (dissociated). So the right inequality is meas(S7(E)) <=
meas(S7(AP_k)), i.e. AP is the GLOBAL MAX. Reframe the proof target as:
   for every primitive E with min=0, |E|=k:   sum_M (-1)^|M| meas(B_M(E)) <= sum_M (-1)^|M| meas(B_M(AP)).

STRATEGY: test the SIGNED per-M difference D_M(E) = (-1)^|M|(meas(B_M(AP)) - meas(B_M(E))).
If sum_M D_M >= 0 always (that's the claim), is each D_M >= 0 (strong) or is there
cancellation? Report the histogram of signs and the dominant contributors.
"""
import sys, itertools
from fractions import Fraction as F
from math import comb, gcd
sys.stdout.reconfigure(line_buffering=True)

def measB_M(E, M):
    E=sorted(set(e for e in E if e!=0)); Ms=set(M)
    if not E: return F(1) if 0 not in Ms else F(0)
    if 0 in Ms: return F(0)  # e=0 pins residue 0; if 0 missed -> impossible (0 always hit)
    bps=set([F(0),F(1)])
    for e in E:
        for m in range(0,7*e+1): bps.add(F(m,7*e))
    bps=sorted(bps); total=F(0)
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2; ok=True
        for e in E:
            if (int(7*e*xm)%7) in Ms: ok=False; break
        if ok: total+=x1-x0
    return total

def measS7_via_M(E):
    s=F(0)
    for r in range(0,8):
        for M in itertools.combinations(range(7),r):
            s+=F((-1)**r)*measB_M(E,M)
    return s

def gen_shapes(k,maxE):
    out=[]
    for rest in itertools.combinations(range(1,maxE+1),k-1):
        E=(0,)+rest; g=0
        for e in E: g=gcd(g,e)
        if g!=1: continue
        out.append(E)
    return out

print("="*92)
print("Per-M signed-difference analysis vs AP. D_M(E)=(-1)^|M|(B_M(AP)-B_M(E)).")
print("Claim: sum_M D_M(E) = meas(S7(AP))-meas(S7(E)) >= 0 for all primitive E.")
print("Question: is each D_M>=0, or is there sign cancellation?")
print("="*92)

for k in [8]:
    maxE=k+4
    AP=tuple(range(k))
    BM_AP={M:measB_M(AP,M) for r in range(8) for M in itertools.combinations(range(7),r)}
    shapes=gen_shapes(k,maxE)
    # Look at the shapes closest to AP (largest meas(S7)) and a few far ones
    data=[]
    for E in shapes:
        s7=measS7_via_M(E); data.append((s7,E))
    data.sort(reverse=True)
    apS7=measS7_via_M(AP)
    print(f"\nk={k}: meas(S7(AP))={float(apS7):.6f}. Top-5 competitor shapes (closest to AP):")
    # the AP itself should be top; show next ones
    cnt=0
    for s7,E in data:
        if E==AP: continue
        cnt+=1
        # per-M sign analysis
        allpos=True; npos=0; nneg=0; maxneg=F(0); negM=None
        for r in range(8):
            for M in itertools.combinations(range(7),r):
                d=F((-1)**r)*(BM_AP[M]-measB_M(E,M))
                if d>0: npos+=1
                elif d<0:
                    nneg+=1; allpos=False
                    if d<maxneg: maxneg=d; negM=M
        diff=apS7-s7
        print(f"  E={E}: meas(S7)={float(s7):.5f} (AP-E={float(diff):+.5f}) "
              f"perM: {npos} pos, {nneg} neg; each_D>=0? {allpos}"
              + ("" if allpos else f"  worst neg D_M={float(maxneg):+.5f} at M={negM}"))
        if cnt>=6: break

print()
print("="*92)
print("Per-M extremality of the SINGLE-term: which M does AP win/lose? (t=|M|=1 detail)")
print("AP loses on t=1 (min not AP). But t=1 enters with sign -1; the LOSS is when AP has")
print("LARGER B_M (so -(larger) is worse). Check: does AP MAXIMIZE meas(B_M) for ALL t>=2?")
print("="*92)
for k in [8]:
    maxE=k+4; AP=tuple(range(k)); shapes=gen_shapes(k,maxE)
    for t in range(1,7):
        # take consecutive block M={0..t-1}? but 0 in M -> 0. Use M={1..t}
        win_all=True; details=[]
        for M in itertools.combinations(range(7),t):
            if 0 in M: continue
            apv=measB_M(AP,M)
            # is AP the max over shapes?
            mx=max(measB_M(E,M) for E in shapes)
            mn=min(measB_M(E,M) for E in shapes)
            sign=(-1)**t
            # IE-correct extremal: even t want MAX, odd t want MIN
            if sign>0:
                ok=(apv==mx); details.append(('MAX',M,ok))
            else:
                ok=(apv==mn); details.append(('MIN',M,ok))
            if not ok: win_all=False
        nok=sum(1 for _,_,o in details if o); ntot=len(details)
        print(f"  t={t} (sign {(-1)**t:+d}, want {'MAX' if (-1)**t>0 else 'MIN'}): "
              f"AP is IE-extremal on {nok}/{ntot} blocks M (0 notin M)")
print("\nDONE.")
