#!/usr/bin/env python3
"""
lrc14_rho17_floor_kps-S5-wf.py  (kind-pasteur-2026-06-18-S5)
Confirm the COMBINED rho*_{1/7} floor c0' over admissible (P,E), |P|+|E|=13, k=|E|>=3.
  - k<=7: rho*_1/7 = meas(G_P) (since mu_1/7=1 => Good=whole circle) >= m_P.
  - k>=8: rho*_1/7 >= meas(G_P)+mu_1/7(E)-1 >= 1891/5880  (granting mu_1/7(E)>=mu_1/7(consec_k)>=thr_k).
  => combined c0' = min(m_P, 1891/5880) = m_P = 14249/252252.
We directly compute rho*_1/7 exactly on a broad admissible sample and confirm:
  (i) min observed rho*_1/7 = m_P, attained at k=3 (or wherever Good=circle);
  (ii) no rho*_1/7 below m_P anywhere;
  (iii) the union LB is a valid lower bound (rho* >= union LB) on the sample.
"""
import sys, itertools, random
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
random.seed(99)
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
def good17_arcs(E):
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
    return merge(good)
def intersect(A,B):
    A=merge(A);B=merge(B);out=[];i=j=0
    while i<len(A) and j<len(B):
        lo=max(A[i][0],B[j][0]); hi=min(A[i][1],B[j][1])
        if lo<hi: out.append((lo,hi))
        if A[i][1]<B[j][1]: i+=1
        else: j+=1
    return out

mP=F(14249,252252)
print(f"m_P = {mP} (~{float(mP):.6f})")
minrho=F(2); minat=None; n_below_mP=0; union_violations=0; tested=0
for trial in range(6000):
    k=random.randint(3,13); psz=13-k
    P=tuple(sorted(random.sample(range(1,14),psz))) if psz>0 else ()
    sp=random.choice([k, k+2, 2*k, 3*k])
    E=sorted(set([0]+random.sample(range(1,sp+1),k-1)))
    if len(E)<k: continue
    GP=safe_set(list(P)); GOOD=good17_arcs(E)
    rho=meas(intersect(GP,GOOD)); mgp=meas(GP)
    tested+=1
    if rho<minrho: minrho=rho; minat=(k,P,tuple(E))
    if rho<mP: n_below_mP+=1
    # union LB check
    ulb=mgp+meas(GOOD)-1
    if rho < ulb - F(1,10**9): union_violations+=1
print(f"tested {tested} admissible (P,E)")
print(f"min rho*_1/7 observed = {minrho} (~{float(minrho):.6f}) at (k,P,E)={minat}")
print(f"# with rho*_1/7 < m_P : {n_below_mP}")
print(f"# union-LB violations (rho < meas(GP)+meas(GOOD)-1): {union_violations}")
print()
# Confirm k<=7 branch: rho = meas(GP) exactly (Good=whole circle)
print("k<=7 spot check (Good should be whole circle => rho=meas(GP)):")
for k in [3,5,7]:
    psz=13-k
    P=tuple(sorted(random.sample(range(1,14),psz)))
    E=sorted(set([0]+random.sample(range(1,3*k+1),k-1)))
    while len(E)<k: E=sorted(set([0]+random.sample(range(1,3*k+1),k-1)))
    GP=safe_set(list(P)); GOOD=good17_arcs(E)
    print(f"  k={k} P={P}: meas(Good)={meas(GOOD)} (==1?{meas(GOOD)==1}); rho={meas(intersect(GP,GOOD))} meas(GP)={meas(GP)} eq:{meas(intersect(GP,GOOD))==meas(GP)}")
print()
print(f"CONCLUSION: combined c0' = min(m_P, 1891/5880) = m_P = {mP} (~{float(mP):.6f}), if min==m_P and 0 below.")
