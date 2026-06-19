#!/usr/bin/env python3
"""
lrc14_wsb_pN6_proof_stage5_kps-S9-wf.py  (kind-pasteur-2026-06-19-S9)

Stage 5.  Focus: the cleanly-provable pieces of the empty-sector law that consec extremizes.

LEMMA-CANDIDATE A (P(N=6) bound).  For ANY cluster E (0 in E, |E|=k, positive distinct speeds):
    P(N_E = 6) = meas{ x : frac(e x) in [0,1/7)  for all e in E }  <=  1/(7 * max(E)),
with EQUALITY iff E = {0,1,...,k-1} up to the binding speed, where max(E)=k-1 is the MINIMUM
possible (k distinct nonneg integers incl 0). So consec UNIQUELY maximizes P(N=6) = 1/(7(k-1)).

  PROOF SKETCH to verify:  Let M=max(E).  The event requires frac(M x) in [0,1/7), i.e.
  x in union_{j=0}^{M-1} [j/M, j/M + 1/(7M)).  That's measure 1/7 already (M arcs).  But ALSO
  e=1 in E? not necessarily.  The cleanest: the event S = {x: all frac(ex)<1/7} is a union of
  intervals; on the component containing 0 it is [0, 1/(7M)) ... need all e, binding is largest.
  Actually for a GENERAL E the bound is meas(S) <= 1/(7M)*(#components).  Test the exact claim
  and find the right RHS (is it 1/(7M), or something involving the whole set?).

LEMMA-CANDIDATE B (p_0 = meas(S7) extremality).  Test whether consec maximizes meas(S7) itself
  among same-k over WIDE spans -- if it holds universally we don't even need L_y.  If it FAILS
  at some span, record the first failure (that pins WHY the L_y combination is needed).

This stage uses a FASTER law/p6 engine and pushes the box wide for k=8 with early reporting.
"""
import sys, itertools
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)

def measAllInArc(E,a,b):
    """meas{x in [0,1): frac(e x) in [a,b) for all e in E\\{0}} as exact union-of-intervals measure."""
    E=sorted(set(e for e in E if e!=0)); bps=set([F(0),F(1)])
    for e in E:
        for t in (a,b):
            for m in range(e): bps.add((t+m)/e)
    bps=sorted(z for z in bps if F(0)<=z<F(1)); tot=F(0)
    for i in range(len(bps)):
        x0=bps[i]; x1=bps[i+1] if i+1<len(bps) else F(1)
        if x1<=x0: continue
        xm=(x0+x1)/2
        if all(a<=(e*xm)%1<b for e in E): tot+=x1-x0
    return tot

print("="*78)
print("LEMMA-CANDIDATE A:  P(N=6)(E) = meas{all frac(ex)<1/7} <= 1/(7 max(E)), equality iff consec")
print("="*78)
print("exhaustive small-k, all primitive E with maxE<=B:")
for k in (4,5,6,7,8):
    for B in (k+2, k+5):
        worst_ratio=(F(0),None)  # P(N=6) / (1/(7 max E))
        n_over=0; cnt=0; maxval=(F(0),None)
        for r in itertools.combinations(range(1,B+1),k-1):
            E=[0]+list(r); cnt+=1
            m=measAllInArc(E,F(0),F(1,7))
            M=max(E); bound=F(1,7*M)
            if m>bound+F(1,10**12): n_over+=1
            ratio=m/bound
            if ratio>worst_ratio[0]: worst_ratio=(ratio,E)
            if m>maxval[0]: maxval=(m,E)
        print(f"  k={k} maxE<={B}: {cnt} sets; P(N=6)>1/(7maxE) violations={n_over}; "
              f"max ratio={float(worst_ratio[0]):.4f} (=1 at {worst_ratio[1]}); "
              f"global max P(N=6)={float(maxval[0]):.6f} at {maxval[1]}")

print()
print("="*78)
print("structural: the all-in-[0,1/7) set near 0 is EXACTLY [0, 1/(7 max E)) -- verify, & list comps")
print("="*78)
def arc_components(E,a,b):
    E=sorted(set(e for e in E if e!=0)); bps=set([F(0),F(1)])
    for e in E:
        for t in (a,b):
            for m in range(e): bps.add((t+m)/e)
    bps=sorted(z for z in bps if F(0)<=z<F(1)); out=[]
    for i in range(len(bps)):
        x0=bps[i]; x1=bps[i+1] if i+1<len(bps) else F(1)
        if x1<=x0: continue
        xm=(x0+x1)/2
        if all(a<=(e*xm)%1<b for e in E): out.append((x0,x1))
    out.sort(); mg=[]
    for a2,b2 in out:
        if mg and a2<=mg[-1][1]: mg[-1]=(mg[-1][0],max(mg[-1][1],b2))
        else: mg.append((a2,b2))
    return mg
for E in [[0,1,2,3,4,5,6,7],[0,2,3,5,7],[0,1,3,7,11],[0,3,6,9,12]]:
    comps=arc_components(E,F(0),F(1,7))
    M=max(E)
    near0 = comps[0] if comps else None
    print(f"  E={E} maxE={M}: #components={len(comps)}, total={sum(b-a for a,b in comps)}"
          f"  component@0={near0}  (expect [0,1/(7M))=[0,{F(1,7*M)}))")

print()
print("="*78)
print("LEMMA-CANDIDATE B: consec maximizes meas(S7) over same-k?  -- where does it FIRST fail?")
print("="*78)
def p0_measS7(E):
    """meas(S7) = meas{ all 7 sectors {0..6} hit by orbit {frac(e x)} }, e in E (0 in E => sector0)."""
    E=sorted(set(E)); bps=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for j in range(7):
            c=F(j,7)
            for m in range(e): bps.add((c+m)/e)
    bps=sorted(z for z in bps if F(0)<=z<F(1)); tot=F(0)
    for i in range(len(bps)):
        x0=bps[i]; x1=bps[i+1] if i+1<len(bps) else F(1)
        if x1<=x0: continue
        xm=(x0+x1)/2
        hit=set()
        for e in E:
            fr=(e*xm)%1; hit.add((fr.numerator*7)//fr.denominator)
        if len(hit)==7: tot+=x1-x0
    return tot
for k in (8,):
    pc=p0_measS7(list(range(k)))
    print(f"  k={k}  meas(S7)(consec)={pc}={float(pc):.6f}")
    for B in (12,14,16):
        over=0; cnt=0; worst=(F(0),None)
        for r in itertools.combinations(range(1,B+1),k-1):
            E=[0]+list(r); cnt+=1
            m=p0_measS7(E)
            if m>pc+F(1,10**12):
                over+=1
                if m>worst[0]: worst=(m,E)
        msg=f"  maxE<={B}: {cnt} sets, meas(S7)>consec: {over}"
        if worst[1]: msg+=f"  WORST={float(worst[0]):.6f} at {worst[1]}"
        print(msg)
print("\nDONE stage 5.")
