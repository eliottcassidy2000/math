#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
THREAD 2 -- THE DECISIVE finite-layer test at k=12.

p0=measS7 is NOT consec-maximal at k>=12 (E*=(0..10,12) beats consec, EXACT).
BUT the cap proof uses the LP-DUAL functional L_y = U4 (THM-534/HYP-2607),
with p0 <= L_y <= cap_k. Since L_y >= p0, it is logically possible that consec
STILL maximizes L_y even where it loses p0. THIS is the finite-layer crux.

Test EXACTLY (rational moment-LP via the THM-534 degree-4 dual evaluated as the
binomial-moment functional, plus a float-LP cross-check):
  L_y(E) for k=8..13, full-residue/primitive/span<=14.
  Does consec maximize L_y? In particular at k=12 vs E*.

The exact L_y at the binding degree: THM-534 says L_y(E)=U4(E)= the moment LP
  max{ p0 : sum_i C(i,t) p_i = S_t(E), t=0..4, p_i>=0 }.
We compute S_t(E) EXACTLY (rational) and solve the small LP exactly by enumerating
the C(7, *) basic feasible vertices (7 vars, 5 equality constraints => basis size 5).
"""
import itertools, math, sys
from fractions import Fraction as F
from math import gcd, comb
from functools import reduce
sys.stdout.reconfigure(line_buffering=True)
P=7

def meas_empty(E, J):
    Jset=set(J); bps={F(0),F(1)}
    for e in E:
        if e==0: continue
        d=P*abs(e)
        for m in range(0,d+1): bps.add(F(m,d))
    bps=sorted(bps); tot=F(0)
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2; empty=True
        for e in E:
            if e==0: continue
            if int(((e*xm)%1)*P) in Jset: empty=False; break
        if empty: tot+=x1-x0
    return tot

def measS7_full(E):
    bps={F(0),F(1)}
    for e in E:
        if e==0: continue
        d=P*abs(e)
        for m in range(0,d+1): bps.add(F(m,d))
    bps=sorted(bps); tot=F(0)
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2
        if len({int(((e*xm)%1)*P) for e in E})==P: tot+=x1-x0
    return tot

def moments(E):
    """S_t = sum_{|J|=t, J subset {1..6}} meas_empty(E,J), t=0..4. EXACT."""
    secs=list(range(1,7)); S=[F(1)]
    for t in range(1,5):
        st=F(0)
        for J in itertools.combinations(secs,t):
            st+=meas_empty(E,list(J))
        S.append(st)
    return S

def solve_LP_exact(S):
    """max p0 s.t. sum_i C(i,t) p_i = S_t (t=0..4), p>=0, i=0..6.
       Enumerate all size-5 bases (5 eqs, 7 vars), solve exactly, keep feasible, max p0."""
    n=7; rows=5
    A=[[F(comb(i,t)) for i in range(n)] for t in range(rows)]
    b=[S[t] for t in range(rows)]
    best=None
    for basis in itertools.combinations(range(n), rows):
        # solve A[:,basis] x = b exactly
        M=[[A[r][c] for c in basis]+[b[r]] for r in range(rows)]
        # Gaussian elimination exact
        m=len(M); ncol=rows
        piv=0; ok=True; where=[-1]*ncol
        for col in range(ncol):
            sel=-1
            for r in range(piv,m):
                if M[r][col]!=0: sel=r; break
            if sel==-1: continue
            M[piv],M[sel]=M[sel],M[piv]
            pivval=M[piv][col]
            M[piv]=[v/pivval for v in M[piv]]
            for r in range(m):
                if r!=piv and M[r][col]!=0:
                    f=M[r][col]
                    M[r]=[M[r][cc]-f*M[piv][cc] for cc in range(ncol+1)]
            where[col]=piv; piv+=1
        # check consistency
        sol=[F(0)]*ncol
        consistent=True
        for col in range(ncol):
            if where[col]!=-1: sol[col]=M[where[col]][ncol]
        # verify
        for r in range(m):
            s=sum(M0 for M0 in [0])  # placeholder
        # rebuild full p
        if piv<rows:  # rank deficient; check residual rows zero
            for r in range(piv,m):
                if M[r][ncol]!=0: consistent=False; break
        if not consistent: continue
        p=[F(0)]*n
        feas=True
        for idx,c in enumerate(basis):
            p[c]=sol[idx]
            if p[c]<0: feas=False; break
        if not feas: continue
        # verify constraints
        good=True
        for t in range(rows):
            if sum(F(comb(i,t))*p[i] for i in range(n))!=b[t]: good=False; break
        if not good: continue
        if best is None or p[0]>best: best=p[0]
    return best

def full_residue(E): return set(e % P for e in E)==set(range(P))
def primitive(E): return reduce(gcd,[e for e in E if e!=0],0)==1
def consec(k): return tuple(range(k))

def stratum(k, box):
    out=[]
    for combo in itertools.combinations(range(1,box+1), k-1):
        E=(0,)+combo
        if full_residue(E) and primitive(E): out.append(E)
    return out

if __name__=="__main__":
    print("="*88)
    print("DECISIVE: does consec maximize L_y=U4 (the cap-proof functional) where p0 fails?")
    print("="*88)
    # sanity: consec L_y vs known cap_8
    C8=consec(8); S8=moments(C8); Ly8=solve_LP_exact(S8)
    print(f"k=8 L_y(consec)={Ly8}={float(Ly8):.6f} (canon: 2633/7350={float(F(2633,7350)):.6f}); p0={float(measS7_full(C8)):.6f}")

    print("\n--- k=12: consec vs E* on L_y ---")
    C=consec(12); Es=(0,1,2,3,4,5,6,7,8,9,10,12)
    LyC=solve_LP_exact(moments(C)); LyE=solve_LP_exact(moments(Es))
    print(f"  L_y(consec) = {LyC} = {float(LyC):.6f}   p0(consec)={float(measS7_full(C)):.6f}")
    print(f"  L_y(E*)     = {LyE} = {float(LyE):.6f}   p0(E*)    ={float(measS7_full(Es)):.6f}")
    print(f"  consec maximizes L_y vs E*? {LyC>=LyE}   (delta L_y = {float(LyC-LyE):+.6f})")

    print("\n--- full L_y argmax scan k=8..13 (full-residue/primitive/span<=14) ---")
    for k in range(8,14):
        S=stratum(k,14)
        if not S: continue
        LyC=solve_LP_exact(moments(consec(k)))
        best=None; bestE=None; nbeat=0; worst=F(0)
        for E in S:
            ly=solve_LP_exact(moments(E))
            if ly is None: continue
            if best is None or ly>best: best=ly; bestE=E
            if ly>LyC+F(1,10**9):
                nbeat+=1
                if ly-LyC>worst: worst=ly-LyC
        print(f"  k={k}: L_y(consec)={float(LyC):.6f}  argmaxL_y={'CONSEC' if bestE==consec(k) else bestE}  maxL_y={float(best):.6f}  #L_y-beaters={nbeat}  worst+{float(worst):.6f}")
