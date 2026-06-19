#!/usr/bin/env python3
"""
lrc14_wsb_lemmaA_components_stage8_kps-S9-wf.py  (kind-pasteur-2026-06-19-S9)

Stage 8.  RIGOROUS structure of LEMMA A: P(N=6)(E) = meas(G_E) <= 1/(7(k-1)) for primitive E,
where G_E = {x in [0,1): frac(ex) in [0,1/7) for all e in E}.

PER-COMPONENT BOUND (claimed, to verify):  every connected component of G_E has length
<= 1/(7 max(E)).  PROOF: on a component, the winding integers n_e = floor(ex) are constant; the
constraint frac(ex)<1/7 reads x in [n_e/e, n_e/e + 1/(7e)).  Intersecting over e in E gives an
interval of length <= min_e 1/(7e) = 1/(7 max(E)).  => meas(G_E) <= #components / (7 max(E)).

So LEMMA A reduces to a COMPONENT-COUNT bound:  #components(G_E) <= max(E)/(k-1).
For consec: max(E)=k-1, #components=1, ratio = 1/(k-1) ... wait that gives bound 1/(7(k-1)) with
1 component.  For a general primitive E: we need  #comp / max(E) <= 1/(k-1), i.e.
#comp <= max(E)/(k-1).  Verify this combinatorial bound (it's a statement about how many x-cells
admit all-k frac in [0,1/7)).

GOALS
 (G1) verify per-component length <= 1/(7 maxE) exactly (all components, many E).
 (G2) verify the component-count bound #comp(G_E) <= max(E)/(k-1) for primitive E; find tightness.
 (G3) the AP-confinement generalization: meas(G_E) <= meas(G_{consec_k}) directly (Lemma A's claim)
      -- exhaustive primitive, plus the structure of WHICH E approach equality.
 (G4) connect to S7: is meas(S7) controlled by a similar 'all components short' argument? (the
      S7 set is the all-7-sectors-hit set; its complement = union of A_j; harder.  Just document
      whether an analogous per-component bound exists for the binding sector A_1 (sector 1 empty).)
"""
import sys, itertools
from math import gcd
from functools import reduce
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
def gcd_all(E): return reduce(gcd,[e for e in E if e!=0],0)

def GE_components(E):
    """components of {x: frac(ex) in [0,1/7) for all e in E\\{0}} as (lo,hi) merged intervals."""
    E=sorted(set(e for e in E if e!=0)); a,b=F(0),F(1,7); bps=set([F(0),F(1)])
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

print("="*78)
print("(G1) per-component length <= 1/(7 maxE)?  (the rigorous building block)")
print("="*78)
viol=0; tested=0; maxlen_ratio=(F(0),None)
for k in (4,5,6,7,8):
    for r in itertools.combinations(range(1,13),k-1):
        E=[0]+list(r)
        if gcd_all(E)!=1: continue
        M=max(E); comps=GE_components(E); tested+=1
        for lo,hi in comps:
            L=hi-lo; ratio=L*7*M
            if L>F(1,7*M)+F(1,10**14): viol+=1
            if ratio>maxlen_ratio[0]: maxlen_ratio=(ratio,E,(lo,hi))
print(f"  tested {tested} primitive E (k=4..8, maxE<=12); per-comp length>1/(7maxE) violations: {viol}")
print(f"  max (component length)*(7 maxE) = {float(maxlen_ratio[0]):.4f} (should be <=1)  @ {maxlen_ratio[1]}")

print()
print("="*78)
print("(G2) component-count bound:  #comp(G_E) <= max(E)/(k-1) for primitive E?")
print("="*78)
viol=0; tested=0; worst=(F(0),None)
for k in (4,5,6,7,8):
    for r in itertools.combinations(range(1,13),k-1):
        E=[0]+list(r)
        if gcd_all(E)!=1: continue
        M=max(E); nc=len(GE_components(E)); tested+=1
        ratio=F(nc*(k-1), M)
        if F(nc, 1) > F(M, k-1)+F(1,10**9): viol+=1
        if ratio>worst[0]: worst=(ratio, E, nc, M)
print(f"  tested {tested} primitive E; #comp > maxE/(k-1) violations: {viol}")
if worst[1]:
    r,E,nc,M=worst
    print(f"  tightest: E={E} #comp={nc} maxE={M} k-1={len(E)-1} -> #comp*(k-1)/maxE={float(r):.4f} (<=1)")

print()
print("="*78)
print("(G3) LEMMA A direct: meas(G_E) <= meas(G_consec)=1/(7(k-1)); equality structure")
print("="*78)
for k in (5,6,7,8):
    cc=F(1,7*(k-1)); over=0; cnt=0; eqsets=[]; near=[]
    for r in itertools.combinations(range(1,14),k-1):
        E=[0]+list(r)
        if gcd_all(E)!=1: continue
        cnt+=1; m=sum(hi-lo for lo,hi in GE_components(E))
        if m>cc+F(1,10**12): over+=1
        if m==cc: eqsets.append(E)
        near.append((m,E))
    near.sort(reverse=True)
    print(f"  k={k}: {cnt} primitive sets; over 1/(7(k-1)): {over}; #equality={len(eqsets)}")
    print(f"        top-3 meas(G_E): "+", ".join(f"{float(m):.5f}@{E}" for m,E in near[:3]))

print()
print("="*78)
print("(G4) sector-1-empty set A_1={x: no frac(ex) in [1/7,2/7)} -- analogous per-comp structure?")
print("="*78)
# A_j = sector j empty.  meas(S7)=1-meas(union A_j).  The binding single-sector measure.
def Aj_meas(E,j):
    E=sorted(set(e for e in E if e!=0)); a,b=F(j,7),F(j+1,7); bps=set([F(0),F(1)])
    for e in E:
        for t in (a,b):
            for m in range(e): bps.add((t+m)/e)
    bps=sorted(z for z in bps if F(0)<=z<F(1)); tot=F(0)
    for i in range(len(bps)):
        x0=bps[i]; x1=bps[i+1] if i+1<len(bps) else F(1)
        if x1<=x0: continue
        xm=(x0+x1)/2
        if all(not (a<=(e*xm)%1<b) for e in E): tot+=x1-x0
    return tot
print("  meas(A_j) (sector j empty) for consec_8, j=1..6 (sym?):")
E=list(range(8))
print("   ", [f"{j}:{float(Aj_meas(E,j)):.4f}" for j in range(1,7)])
print("  meas(S7)=1-meas(union)=", float(F(481,1470)))
print("\nDONE stage 8.")
