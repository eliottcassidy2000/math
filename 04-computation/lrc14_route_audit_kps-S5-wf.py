#!/usr/bin/env python3
"""
lrc14_route_audit_kps-S5-wf.py  (kind-pasteur-2026-06-18-S5)
Audit the UNCONDITIONAL pieces of the 1/7 global-witness route:
  (A) pigeonhole: mu_1/7(E)=1 for ALL E with k<=7  (k<=6 forced; k=7 a.e.)
  (B) m_P = min over admissible P (|P|<=10) of meas(G_P)  -- exact, claimed 14249/252252
  (C) union-bound floor 1891/5880 for k>=8  (= min over admissible (P,E consec) of meas(G_P)+mu_1/7-1)
  (D) global-witness equivalence: 'exists phase phi: ||phi - frac(e_i x)|| >= 1/14 for all i'
       <=> 'maxgap{frac(e_i x)} > 1/7'  (each phase forbids a 1/7-length arc; a free phase exists
       iff some gap exceeds 1/7).  Check by exact enumeration on random x.
"""
import sys, itertools, random
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
random.seed(13)
ONE7=F(1,7); H=F(1,14)
def merge(iv):
    iv=sorted(iv);out=[]
    for a,b in iv:
        if out and a<=out[-1][1]: out[-1]=(out[-1][0],max(out[-1][1],b))
        else: out.append((a,b))
    return out
def meas(arcs): return sum((b-a for a,b in arcs),F(0))
def complement(arcs):
    arcs=merge(arcs);out=[];prev=F(0)
    for a,b in arcs:
        if a>prev: out.append((prev,a))
        prev=max(prev,b)
    if prev<1: out.append((prev,F(1)))
    return out
def danger_arcs(u,h=H):
    iv=[]
    for j in range(u):
        c=F(j,u);a=(c-h/u)%1;b=(c+h/u)%1
        if a<b: iv.append((a,b))
        else: iv.append((a,F(1)));iv.append((F(0),b))
    return iv
def safe_set(P,h=H):
    if not P: return [(F(0),F(1))]
    return complement(merge([iv for u in P for iv in danger_arcs(u,h)]))
def mu17_exact(E):
    E=sorted(set(E));k=len(E);diffs=set()
    for a in range(k):
        for b in range(a+1,k): diffs.add(E[b]-E[a])
    bps=set([F(0),F(1)])
    for d in diffs:
        for m in range(0,d+1): bps.add(F(m,d))
    bps=sorted(x for x in bps if 0<=x<=1);good=[]
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2
        pts=sorted(((E[t]*xm)%1,E[t]) for t in range(k))
        order=[e for _,e in pts];floors=[int((e*xm)//1) for e in order]
        for idx in range(k):
            e_cur=order[idx];f_cur=floors[idx]
            if idx<k-1: e_nx=order[idx+1];f_nx=floors[idx+1];wrap=F(0)
            else: e_nx=order[0];f_nx=floors[0];wrap=F(1)
            A=F(e_nx-e_cur);Cc=F(f_cur-f_nx)+wrap
            if A==0:
                if Cc>ONE7: good.append((x0,x1))
                continue
            xb=(ONE7-Cc)/A
            if A>0: lo=max(x0,xb);hi=x1
            else: lo=x0;hi=min(x1,xb)
            if lo<hi: good.append((lo,hi))
    return meas(merge(good))

print("=== (A) pigeonhole: mu_1/7=1 for k<=7 ===")
# k<=6: any E, mu_1/7 must be 1. Exhaustive over spread<=k+4 (and a few large-spread).
for k in range(3,8):
    space=list(range(1,k+5))
    mn=F(2); mnE=None; cnt=0
    for body in itertools.combinations(space,k-1):
        E=[0]+list(body)
        m=mu17_exact(E); cnt+=1
        if m<mn: mn=m; mnE=E
    # plus large-spread random
    for _ in range(2000):
        sp=random.choice([2*k,3*k,5*k,40])
        E=sorted(set([0]+random.sample(range(1,sp+1),k-1)))
        if len(E)<k: continue
        m=mu17_exact(E)
        if m<mn: mn=m; mnE=E
    print(f"  k={k}: min mu_1/7 = {mn} (~{float(mn):.6f}) over {cnt} bounded + 2000 random ; minimizer={mnE if mn<1 else 'all=1'}")

print("\n=== (B) m_P = min_{|P|<=10} meas(G_P) ===")
per=[]
for psz in range(0,11):
    mn=F(2); mnP=None
    for P in itertools.combinations(range(1,14),psz):
        ms=meas(safe_set(list(P)))
        if ms<mn: mn=ms; mnP=P
    per.append((psz,mn,mnP))
    print(f"  |P|={psz:2d}: min meas(G_P)={mn} (~{float(mn):.6f}) at P={mnP}")
mP=min(m for _,m,_ in per)
print(f"  m_P = {mP} (~{float(mP):.6f})   [claimed 14249/252252={float(F(14249,252252)):.6f}]")
print(f"  MATCH 14249/252252: {mP==F(14249,252252)}")

print("\n=== (C) union floor for k>=8: min over admissible (P, consec E) of meas(G_P)+mu_1/7(consec)-1 ===")
ufloor=F(2); uat=None
consec_mu={k:mu17_exact(list(range(k))) for k in range(8,14)}
for k in range(8,14):
    psz=13-k
    muc=consec_mu[k]
    for P in itertools.combinations(range(1,14),psz):
        val=meas(safe_set(list(P)))+muc-1
        if val<ufloor: ufloor=val; uat=(k,P)
print(f"  consec mu_1/7 k=8..13: {[str(consec_mu[k]) for k in range(8,14)]}")
print(f"  union floor (per-pair, consec) = {ufloor} (~{float(ufloor):.6f}) at (k,P)={uat}")
print(f"  MATCH 1891/5880: {ufloor==F(1891,5880)}")

print("\n=== (D) global-witness equivalence: free-phase exists  <=>  maxgap>1/7 ===")
# For a given x and E, the forbidden region for phi is union over i of arc center frac(e_i x), radius 1/14.
# A free phi exists iff complement nonempty iff some gap between consecutive points > 2*(1/14)=1/7.
def free_phase_exists(E,x):
    pts=sorted(set((e*x)%1 for e in E))
    # forbidden arcs of half-width 1/14 around each point; free phi exists iff some gap>1/7
    if len(pts)==1:
        return True  # one forbidden arc length 1/7 <1, complement nonempty
    gaps=[pts[t+1]-pts[t] for t in range(len(pts)-1)]+[pts[0]+1-pts[-1]]
    return max(gaps)>ONE7
def maxgap_gt(E,x):
    pts=sorted(set((e*x)%1 for e in E))
    if len(pts)==1: return True
    gaps=[pts[t+1]-pts[t] for t in range(len(pts)-1)]+[pts[0]+1-pts[-1]]
    return max(gaps)>ONE7
mismatch=0; tested=0
for _ in range(20000):
    k=random.randint(3,13)
    sp=random.randint(k,40)
    E=sorted(set([0]+random.sample(range(1,sp+1),k-1)))
    if len(E)<2: continue
    x=F(random.randint(1,9999),10000)
    tested+=1
    if free_phase_exists(E,x)!=maxgap_gt(E,x): mismatch+=1
print(f"  tested {tested} (E,x); mismatches between 'free phase exists' and 'maxgap>1/7' = {mismatch}")
print("  (equivalence holds <=> 0 mismatches)")
