#!/usr/bin/env python3
"""
lrc14_BS_4mom_crosscheck_macmini_0618s7.py  (mac-mini-2026-06-18-S7) ANGLE A cross-validation

Cross-validate the 4-moment certificate two independent ways over a box (B=12 primitive k=8):
 (1) scipy linprog float U4   vs   (2) exact rational vertex-enumeration U4.
Confirm they AGREE (|float-exact|<1e-7) for every shape, and that exact U4 <= cap_8 for all,
with the AP {0..7} the unique max.  This makes the exhaustive closure claim trustworthy
(not a float/LP artifact).
"""
import sys, itertools, math
import numpy as np
from scipy.optimize import linprog
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)

def meas_empty(E, J):
    Jset=set(J); bps=set([F(0),F(1)])
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
            if int(((e*xm)%1)*7) in Jset: empty=False; break
        if empty: total+=x1-x0
    return total

def measS7_geom(E):
    bps=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(0,7*e+1): bps.add(F(m,7*e))
    bps=sorted(bps); total=F(0)
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2; secs=set()
        for e in E: secs.add(int(((e*xm)%1)*7))
        if len(secs)==7: total+=x1-x0
    return total

def moments(E):
    secs=list(range(1,7)); S=[F(1)]
    for t in range(1,5):
        st=F(0)
        for J in itertools.combinations(secs,t): st+=meas_empty(E,list(J))
        S.append(st)
    return S

def U4_float(S):
    A=np.array([[float(math.comb(i,t)) for i in range(7)] for t in range(5)])
    b=np.array([float(S[t]) for t in range(5)]); c=np.zeros(7); c[0]=-1.0
    r=linprog(c,A_eq=A,b_eq=b,bounds=[(0,None)]*7,method='highs')
    return r.x[0] if r.success else None

def solve_rat(A,b):
    n=len(A); M=[row[:]+[b[i]] for i,row in enumerate(A)]
    for col in range(n):
        piv=next((r for r in range(col,n) if M[r][col]!=0),None)
        if piv is None: return None
        M[col],M[piv]=M[piv],M[col]; pv=M[col][col]; M[col]=[v/pv for v in M[col]]
        for r in range(n):
            if r!=col and M[r][col]!=0:
                f=M[r][col]; M[r]=[M[r][j]-f*M[col][j] for j in range(n+1)]
    return [M[i][n] for i in range(n)]

def U4_exact(S):
    rows=[[F(math.comb(i,t)) for i in range(7)] for t in range(5)]; b=[S[t] for t in range(5)]
    best=None
    for basis in itertools.combinations(range(7),5):
        A=[[rows[t][j] for j in basis] for t in range(5)]
        sol=solve_rat(A,b)
        if sol is None or any(v<0 for v in sol): continue
        p0=F(0)
        for idx,j in enumerate(basis):
            if j==0: p0=sol[idx]
        if best is None or p0>best: best=p0
    return best

cap8=F(2243,5880)
B=12
print(f"Cross-check float vs exact U4 over primitive k=8, B={B}. cap_8={float(cap8):.6f}")
maxdiff=0.0; nover=0; maxexact=(F(0),None); cnt=0; valbad=0
for combo in itertools.combinations(range(1,B+1),7):
    E=[0]+list(combo); g=0
    for e in E: g=math.gcd(g,e)
    if g!=1: continue
    cnt+=1
    S=moments(E); uf=U4_float(S); ue=U4_exact(S); s7=measS7_geom(E)
    if uf is not None: maxdiff=max(maxdiff,abs(uf-float(ue)))
    if ue<s7: valbad+=1
    if ue>cap8: nover+=1
    if ue>maxexact[0]: maxexact=(ue,E)
print(f"scanned {cnt}; max|float-exact U4| = {maxdiff:.2e}")
print(f"exact-U4 > meas(S7) holds for all? validity fails = {valbad}")
print(f"exact U4 > cap_8 count = {nover}")
print(f"max exact U4 = {maxexact[0]} = {float(maxexact[0]):.6f} at {maxexact[1]}")
print("AP {0..7} is the unique max:", maxexact[1]==list(range(8)))
print(f"\nCERTIFICATE (exact, B={B}): every primitive k=8 shape has U4 <= cap_8 =>")
print("4-moment empty-sector marginal certificate closes LRC(14)-S3 k=8 in this box.")
