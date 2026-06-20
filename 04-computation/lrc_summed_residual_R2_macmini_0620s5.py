#!/usr/bin/env python3
"""
lrc_summed_residual_R2_macmini_0620s5.py  (mac-mini-2026-06-20-S5)

Test the HYP-2692 redirect: is the SUMMED two-far residual R_2 the controlling quantity?
For E=B u F (B bounded core, F far runners), the Newton/boundary decomposition is
   p0(E) = P_r(B) + sum_{s>=1} R_s,   R_s = sum_{|S|=s, S subset F}[Delta_S(B) - Phi_s(B)].
P_r(B) = decorrelated boundary value (safe, margin grows). The corrections R_s.
Compute R_1, R_2, R_3 and the total R = p0 - P_r for true-wide rows with growing F; check whether
R_2 dominates and whether the total stays within the cap margin.
"""
import itertools, math, sys
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
def sector_of(p): return int((p%1)*7)
def measS7(E):
    E=sorted(set(E)); b=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(0,7*e+1): b.add(F(m,7*e))
    b=sorted(b); tot=F(0)
    for i in range(len(b)-1):
        x0,x1=b[i],b[i+1]
        if x1<=x0: continue
        if len(set(sector_of(e*((x0+x1)/2)) for e in E))==7: tot+=x1-x0
    return tot
def miss_profile(B):
    B=sorted(set(B)); b=set([F(0),F(1)])
    for e in B:
        if e==0: continue
        for m in range(0,7*e+1): b.add(F(m,7*e))
    b=sorted(b); prof=[F(0)]*7
    for i in range(len(b)-1):
        x0,x1=b[i],b[i+1]
        if x1<=x0: continue
        t=7-len(set(sector_of(e*((x0+x1)/2)) for e in B))
        if t<=6: prof[t]+=x1-x0
    return prof
def DeltaS(B,S):
    S=list(S); t=F(0)
    for r in range(len(S)+1):
        for T in itertools.combinations(S,r):
            t+=(-1)**(len(S)-len(T))*measS7(list(B)+list(T))
    return t
def Phi_s(B,s):
    prof=miss_profile(B)
    # Phi_s = 7^-s sum_t (-1)^(s-t) t! Stirling2(s,t) p_t
    from math import comb
    def stir2(n,k): return sum((-1)**(k-j)*comb(k,j)*j**n for j in range(k+1))//math.factorial(k)
    return sum(F((-1)**(s-t)*math.factorial(t)*stir2(s,t), 7**s)*prof[t] for t in range(1,s+1))
def c_t(t,r): return sum((-1)**i*math.comb(t,i)*(1-F(i,7))**r for i in range(t+1))
def P_r(B,r):
    prof=miss_profile(B); return sum(prof[t]*c_t(t,r) for t in range(7))

caps={8:F(2243,5880),9:F(1979,4004),10:F(4,7),11:F(5,7),12:F(6,7)}
print(f"{'E (B u F)':<42}{'p0':>8}{'P_r(B)':>8}{'R_1':>8}{'R_2':>8}{'R_3':>8}{'R_tot':>8}{'cap':>7}{'marg':>7}")
cases=[
 ((0,4,6,8,10,12,14),(15,16)),
 ((0,4,6,8,10,12,14),(15,16,17)),
 ((0,2,4,6,8,10,12),(15,16,17,18)),
 ((0,1,2,3,4),(15,16,17,18,19)),
 ((0,1,2,3),(15,16,17,18,19,20)),
]
for B,Fr in cases:
    r=len(Fr); E=tuple(sorted(B+Fr)); k=len(E)
    p0=measS7(E); Pr=P_r(B,r)
    Rs={}
    for s in (1,2,3):
        if s>r: Rs[s]=F(0); continue
        phis=Phi_s(B,s)
        Rs[s]=sum(DeltaS(B,S)-phis for S in itertools.combinations(Fr,s))
    Rtot=p0-Pr; cap=caps.get(k,F(1))
    print(f"{str(E):<42}{float(p0):>8.4f}{float(Pr):>8.4f}{float(Rs[1]):>8.4f}{float(Rs[2]):>8.4f}{float(Rs[3]):>8.4f}{float(Rtot):>8.4f}{float(cap):>7.4f}{float(cap-p0):>7.4f}")
print("\nR_s = summed order-s residual. If R_2 dominates and R_tot stays << margin, the HYP-2692")
print("redirect (height-weighted summed-R_2 bound) is the right lever; P_r(B) is the safe main term.")
