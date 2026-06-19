#!/usr/bin/env python3
"""
lrc14_BS_4mom_exact_AP_macmini_0618s7.py  (mac-mini-2026-06-18-S7) ANGLE A exact certificate

EXACT-RATIONAL value of the 4-moment certificate U4 for the AP (consec k=8), the unique tight
extremiser.  The moment LP  max p_0  s.t.  sum_{i=0}^6 C(i,t) p_i = S_t (t=0..4), p_i>=0  has its
optimum at a basic feasible solution: 5 equality constraints => 5 basic variables, 2 at zero.
We enumerate all C(7,5)=21 choices of basis, solve the 5x5 rational system, keep feasible (p>=0)
vertices, and take max p_0.  All EXACT (fractions).  Confirms U4(consec8) and U4 < cap_8 exactly.

Then prints the exact binomial moments S_1..S_4 for consec8 and the certifying distribution.
"""
import sys, itertools, math
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
        for e in E:
            secs.add(int(((e*xm)%1)*7))
        if len(secs)==7: total+=x1-x0
    return total

def binom_moments_exact(E, upto=4):
    secs=list(range(1,7)); S=[F(1)]
    for t in range(1,upto+1):
        st=F(0)
        for J in itertools.combinations(secs,t):
            st+=meas_empty(E,list(J))
        S.append(st)
    return S

def solve_rational(A, b):
    # Gaussian elimination on exact Fractions; A is nxn list of lists, b length n. Return x or None.
    n=len(A); M=[row[:]+[b[i]] for i,row in enumerate(A)]
    for col in range(n):
        piv=None
        for r in range(col,n):
            if M[r][col]!=0: piv=r; break
        if piv is None: return None
        M[col],M[piv]=M[piv],M[col]
        pv=M[col][col]
        M[col]=[v/pv for v in M[col]]
        for r in range(n):
            if r!=col and M[r][col]!=0:
                f=M[r][col]; M[r]=[M[r][j]-f*M[col][j] for j in range(n+1)]
    return [M[i][n] for i in range(n)]

def U4_exact(S):
    # variables p_0..p_6 (7). constraints t=0..4 (5 eqns). basis = 5 of 7 vars.
    rows=[[F(math.comb(i,t)) for i in range(7)] for t in range(5)]
    bvec=[S[t] for t in range(5)]
    best=None; bestsol=None
    for basis in itertools.combinations(range(7),5):
        A=[[rows[t][j] for j in basis] for t in range(5)]
        sol=solve_rational(A,bvec)
        if sol is None: continue
        if any(v<0 for v in sol): continue
        p=[F(0)]*7
        for idx,j in enumerate(basis): p[j]=sol[idx]
        if best is None or p[0]>best:
            best=p[0]; bestsol=p
    return best, bestsol

if __name__=="__main__":
    cap8=F(2243,5880)
    print("EXACT 4-moment certificate U4 for the AP consec k=8 (the tight extremiser):")
    for name,E in [("consec8",list(range(8))),("dilAP",[0,2,4,6,8,10,12,14]),("consec9",list(range(9)))]:
        S=binom_moments_exact(E,4)
        s7=measS7_geom(E)
        u4,sol=U4_exact(S)
        cap = cap8 if len(E)==8 else None
        print(f"\n{name}  E={E}")
        print(f"  binomial moments S1..S4 = {[str(S[t]) for t in range(1,5)]}")
        print(f"  meas(S7) = {s7} = {float(s7):.6f}")
        print(f"  U4 (exact) = {u4} = {float(u4):.6f}")
        print(f"  certifying count-distribution p_0..p_6 = {[str(x) for x in sol]}")
        if cap is not None:
            print(f"  cap_8 = {cap} = {float(cap):.6f}")
            print(f"  U4 <= cap_8 ?  {u4 <= cap}   (slack = {cap-u4} = {float(cap-u4):+.6f})")
            print(f"  *** EXACT RATIONAL: U4(consec8) {'<' if u4<cap else '>='} cap_8 ***" if name=="consec8" else "")
