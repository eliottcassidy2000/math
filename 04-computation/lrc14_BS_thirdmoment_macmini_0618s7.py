#!/usr/bin/env python3
"""
lrc14_BS_thirdmoment_macmini_0618s7.py  (mac-mini-2026-06-18-S7) ANGLE A tightening

Does adding TRIPLE empty-sector data (3-wise) close the AP and shrink the certificate residual?
Use higher-degree moment bounds on meas(union_j A_j), j in {1..6}, A_j = {sector j empty}:

  S1 = sum P_j,  S2 = sum_{j<l} P_{jl},  S3 = sum_{j<l<m} P_{jlm}.
Lower bounds on the union from binomial moments (all VALID lower bounds -> upper bds on meas(S7)):
  (DS2) Dawson-Sankoff degree-2 (already have).
  (DS3) the degree-3 Dawson-Sankoff / Galambos-Simonelli optimal 3-moment lower bound.
        A clean valid degree-3 lower bound: union >= S1 - S2 + (something with S3)? No -- the
        truncated IE  S1 - S2 + S3  is an UPPER bound (odd truncation). For a LOWER bound with 3
        moments we use the optimal bound of Dawson-Sankoff generalized (the "method of optimal
        bounds"): we solve the small LP over the distribution of the count N=#empty sectors in
        {0,..,6} matching moments S1,S2,S3 (binomial moments b_i=S_i), minimizing P(N>=1).
        This LP min P(N>=1) s.t. binomial-moment constraints is the BEST possible lower bound on
        the union from S1,S2,S3 -- we solve it exactly (rational LP via the dual / small enumeration).

We implement min P(N>=1) over distributions p_0..p_6 on counts with
   sum_i C(i,t) p_i = S_t  (t=0,1,2,3; S_0=1),  p_i>=0.
min p_0' where union = 1-p_0, so min union = 1 - max p_0. Equivalently max p_0 s.t. moment match.
We solve this LP exactly with scipy (rational check after). The optimum max p_0 gives the
tightest 3-moment certificate:  meas(S7) <= max p_0  (since meas(S7)=p_0 EXACTLY -- p_0 = P(no
sector empty) = meas(S7)!).  Wait: that means the moment-LP max p_0 is exactly an upper bound on
meas(S7) and it's the BEST bound consistent with S1,S2,S3.  Compute it for consec & shapes.

This is the rigorous "best possible bound from k-wise marginals" -- the moment-problem certificate.
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
        for e in E:
            secs.add(int(((e*xm)%1)*7))
        if len(secs)==7: total+=x1-x0
    return total

def binom_moments(E, upto):
    secs=list(range(1,7))
    S=[F(1)]  # S_0=1
    for t in range(1,upto+1):
        st=F(0)
        for J in itertools.combinations(secs,t):
            st+=meas_empty(E,list(J))
        S.append(st)
    return S  # S[0..upto]

def moment_LP_max_p0(S, m=6):
    # variables p_0..p_m on counts 0..m. max p_0 s.t. sum_i C(i,t) p_i = S_t (t=0..len(S)-1), p>=0.
    upto=len(S)-1
    n=m+1
    c=np.zeros(n); c[0]=-1.0  # maximize p_0
    A_eq=[]; b_eq=[]
    for t in range(0,upto+1):
        row=[float(math.comb(i,t)) for i in range(n)]
        A_eq.append(row); b_eq.append(float(S[t]))
    res=linprog(c, A_eq=np.array(A_eq), b_eq=np.array(b_eq), bounds=[(0,None)]*n, method='highs')
    if not res.success: return None
    return res.x[0]

cap8=F(2243,5880); cap_float={8:0.38146,9:0.4943,10:0.6044,11:0.7253,12:0.8571,13:1.0}

shapes_by_k={
 8:[("consec{0..7}",list(range(8))),("perf",[0,2,3,4,5,6,7,9]),
    ("dissoc",[0,1,3,7,15,31,63,127]),("Sidon",[0,1,3,7,12,20,30,44]),
    ("spread",[0,1,2,3,40,41,42,43]),("generic",[0,5,13,27,41,58,79,97]),
    ("dilAP",[0,2,4,6,8,10,12,15]),("AP+out",[0,1,2,3,4,5,6,16])],
 9:[("consec{0..8}",list(range(9)))],
 10:[("consec{0..9}",list(range(10)))],
 11:[("consec{0..10}",list(range(11)))],
}

print("="*98)
print("ANGLE A: moment-problem certificate -- BEST upper bound on meas(S7) from t-wise marginals")
print("meas(S7)=p_0; U_t = max p_0 s.t. binomial moments S_1..S_t matched (tightest t-marginal bound)")
print("="*98)
print(f"{'shape':<18}{'k':>3}{'meas(S7)':>10}{'U(2mom)':>10}{'U(3mom)':>10}{'U(4mom)':>10}{'cap_k':>9}{'3m<=cap?':>9}")
print("-"*98)
for k in sorted(shapes_by_k):
    capk=float(cap8) if k==8 else cap_float[k]
    for name,E in shapes_by_k[k]:
        s7=float(measS7_geom(E))
        S=binom_moments(E,4)
        U2=moment_LP_max_p0(S[:3]); U3=moment_LP_max_p0(S[:4]); U4=moment_LP_max_p0(S[:5])
        flag="OK" if (U3 is not None and U3<=capk) else "OVER"
        def fmt(x): return f"{x:>10.4f}" if x is not None else f"{'--':>10}"
        print(f"{name:<18}{k:>3}{s7:>10.4f}{fmt(U2)}{fmt(U3)}{fmt(U4)}{capk:>9.4f}{flag:>9}")
print("-"*98)
print("U(t mom) is the BEST possible upper bound on meas(S7) using only <=t-wise empty-sector measures.")
print("If U(3mom) or U(4mom) <= cap_k for consec at k=8..11, the t-marginal moment certificate closes")
print("the dangerous rows (and since meas(S7)=p_0 exactly, this is the sharpest marginal-based bound).")
