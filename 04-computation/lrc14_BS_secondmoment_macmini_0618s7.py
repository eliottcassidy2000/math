#!/usr/bin/env python3
"""
lrc14_BS_secondmoment_macmini_0618s7.py  (mac-mini-2026-06-18-S7) ANGLE A pivot

RIGOROUS, EXACT-RATIONAL upper bound on meas(S7(E)) using ONLY pairwise empty-sector data
(no minorants -> sidesteps the real-analyticity impossibility wall that blocks the IE-minorant
and Selberg-minorant routes).

S7 fails  iff  some sector j in {1..6} is empty.  Let A_j = {x: sector j empty},
P_j = meas(A_j),  P_{jl} = meas(A_j cap A_l).  Then meas(S7) = 1 - meas(union_j A_j), and we
LOWER-bound the union measure (=> UPPER-bound meas(S7)) by classical second-moment inequalities
that need only P_j and P_{jl}, all EXACT rationals (breakpoint method):

 (CS)  Cauchy-Schwarz / Chung-Erdos:  meas(union) >= (sum_j P_j)^2 / sum_{j,l} P_{jl}.
 (KOU) Kounias/Hunter (best degree-2 Bonferroni over a spanning tree):
        meas(union) >= max_{spanning tree tau} ( sum_j P_j - sum_{(j,l) in tau} P_{jl} ).
 (DAW) Dawson-Sankoff (optimal 2-moment bound):  with S1=sum P_j, S2=sum_{j<l} P_{jl},
        meas(union) >= 2 S1/(t+1) - 2 S2/(t(t+1)),  t = 1 + floor(2 S2/S1)  (t in 1..m).
We take the BEST (max) of these as the union lower bound L_union, then meas(S7) <= 1 - L_union.

All exact.  Compare to cap_k for consec (extremiser) and several shapes, k=8..13.
If 1 - L_union <= cap_k for the dangerous rows, this single closed-form certificate closes them.
"""
import sys, itertools, math
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)

def meas_empty(E, J):
    """meas{x: for every e in E (e!=0), frac(e x) avoids EVERY sector in J}."""
    E=sorted(set(E)); Jset=set(J)
    bps=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(0,7*e+1): bps.add(F(m,7*e))
    bps=sorted(bps); total=F(0)
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2; empty=True
        for e in E:
            if e==0: continue
            y=(e*xm)%1
            if int(y*7) in Jset: empty=False; break
        if empty: total+=x1-x0
    return total

def measS7_geom(E):
    E=sorted(set(E)); bps=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(0,7*e+1): bps.add(F(m,7*e))
    bps=sorted(bps); total=F(0)
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2; secs=set()
        for e in E:
            y=(e*xm)%1; secs.add(int(y*7))
        if len(secs)==7: total+=x1-x0
    return total

def M7(k):
    s=F(0)
    for t in range(0,7):
        s += F((-1)**t * math.comb(6,t)) * F(7-t,7)**(k-1)
    return s

cap_exact={8:F(2243,5880)}
cap_float={8:0.38146,9:0.4943,10:0.6044,11:0.7253,12:0.8571,13:1.0}

def pair_data(E):
    secs=list(range(1,7))
    P={j: meas_empty(E,[j]) for j in secs}
    Pjl={}
    for j,l in itertools.combinations(secs,2):
        Pjl[(j,l)]=meas_empty(E,[j,l])
    return P,Pjl

def union_CS(P,Pjl,secs):
    S1=sum(P.values(),F(0))
    if S1==0: return F(0)
    den=F(0)
    for j in secs: den+=P[j]
    for j,l in itertools.combinations(secs,2): den+=2*Pjl[(j,l)]
    if den==0: return F(0)
    return S1*S1/den

def union_kounias(P,Pjl,secs):
    # Hunter's bound: union >= sum P_j - sum over a spanning tree of P_{jl}, maximize over trees
    # => use Kruskal: max-weight spanning tree on edges weighted P_{jl}, subtract its weight...
    # Hunter: union >= sum P_j - max-weight-spanning-tree(P_{jl}). To maximize the lower bound we
    # SUBTRACT the MIN total, i.e. min-weight spanning tree? Hunter uses ANY tree; best is min tree.
    S1=sum(P[j] for j in secs)
    edges=sorted(((Pjl[(min(a,b),max(a,b))],a,b) for a,b in itertools.combinations(secs,2)))
    # minimum spanning tree to subtract least:
    parent={j:j for j in secs}
    def find(x):
        while parent[x]!=x: parent[x]=parent[parent[x]]; x=parent[x]
        return x
    w=F(0); cnt=0
    for pw,a,b in edges:
        ra,rb=find(a),find(b)
        if ra!=rb:
            parent[ra]=rb; w+=pw; cnt+=1
            if cnt==len(secs)-1: break
    return S1-w

def union_dawson_sankoff(P,Pjl,secs):
    S1=sum(P[j] for j in secs)
    S2=sum(Pjl[(min(a,b),max(a,b))] for a,b in itertools.combinations(secs,2))
    if S1==0: return F(0)
    # optimal integer t in 1..m
    best=F(0)
    m=len(secs)
    for t in range(1,m+1):
        val=F(2,1)*S1/F(t+1) - F(2,1)*S2/F(t*(t+1))
        if val>best: best=val
    return best

shapes_by_k={
 8:[("consec{0..7}",list(range(8))),
    ("perf{0,2,..,9}",[0,2,3,4,5,6,7,9]),
    ("dissoc 2^i",[0,1,3,7,15,31,63,127]),
    ("Sidon",[0,1,3,7,12,20,30,44]),
    ("spread 0..3,40..43",[0,1,2,3,40,41,42,43]),
    ("generic",[0,5,13,27,41,58,79,97])],
 9:[("consec{0..8}",list(range(9)))],
 10:[("consec{0..9}",list(range(10)))],
 11:[("consec{0..10}",list(range(11)))],
 12:[("consec{0..11}",list(range(12)))],
 13:[("consec{0..12}",list(range(13)))],
}

print("="*104)
print("ANGLE A pivot: rigorous exact-rational meas(S7) <= 1 - max(CS,Kounias,DawsonSankoff) union bound")
print("="*104)
# NOTE on direction: CS (Chung-Erdos) and Dawson-Sankoff are LOWER bounds on the union measure,
# hence give UPPER bounds 1-union on meas(S7). Hunter/Kounias is an UPPER bound on the union ->
# WRONG direction for us (would give a lower bound on meas(S7)); we keep it only as a sanity ceiling
# on the union, and DO NOT use it for the meas(S7) upper bound.
print(f"{'shape':<22}{'k':>3}{'meas(S7)':>10}{'unionLB':>9}{'U=1-LB':>9}{'cap_k':>9}{'slack':>9}{'<=cap?':>8}")
print("-"*104)
secs=list(range(1,7))
for k in sorted(shapes_by_k):
    capk = float(cap_exact[8]) if k==8 else cap_float[k]
    for name,E in shapes_by_k[k]:
        s7=float(measS7_geom(E))
        P,Pjl=pair_data(E)
        uCS=union_CS(P,Pjl,secs); uDS=union_dawson_sankoff(P,Pjl,secs)
        uLB=max(uCS,uDS)               # best (largest) VALID lower bound on the union
        U=F(1)-uLB                     # rigorous upper bound on meas(S7)
        flag="OK" if float(U)<=capk else "OVER"
        print(f"{name:<22}{k:>3}{s7:>10.4f}{float(uLB):>9.4f}{float(U):>9.4f}{capk:>9.4f}{capk-float(U):>9.4f}{flag:>8}")
print("-"*104)
print("U = 1 - max(Chung-Erdos, Dawson-Sankoff) = rigorous exact-rational upper bound on meas(S7).")
print("(Hunter/Kounias is an upper bd on the union -> wrong direction, excluded.)")
print("If U <= cap_k for consec (extremiser) at the dangerous rows, the pairwise certificate closes them.")
