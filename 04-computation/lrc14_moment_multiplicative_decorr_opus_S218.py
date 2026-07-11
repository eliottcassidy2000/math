"""
opus-2026-07-11-S218 (part 2): the moment-LP multiplicative-decorrelation route for the ungapped wide half.

THM-534 (PROVED per-E): p0(E) = meas(S7(E)) <= L_y(E) = sum_r y_r S_r(E),  S_r(E)=sum_{|A|=r} meas{avoid A},
meas{avoid A}=meas{x: no e in E lands in the r/7-region union_{j in A}[j/7,(j+1)/7)}, A subset {1..6}.
k=9 dual (deg 3): L_y = 1 - (13/18)S1 + (4/9)S2 - (1/6)S3.

KEY STRUCTURE (this script tests): peeling w=max E,
  meas_{E'u w}{avoid A} = (1-r/7)*meas_{E'}{avoid A} + delta_A,   |delta_A| <= V_A(E')/(c w)   (THM-700-type)
=> MULTIPLICATIVE recursion (each far element damps by (1-r/7)<1), so S_r(E) -> S_r^iid = C(6,r)(1-r/7)^{k-1}
geometrically, and L_y(E) -> L_y^inf < cap for WIDE E. Unlike the ADDITIVE p0 tautology, this DAMPS the
per-step discrepancy geometrically. Tests: (1) L_y(E)>=p0(E) [dual holds]; (2) L_y(E)<cap for wide;
(3) the (1-r/7) multiplicative ratio; (4) L_y(E)-L_y^inf small for wide (the decorrelation).
"""
from fractions import Fraction as F
from itertools import combinations
from math import comb

def breakpoints(E):
    Es=[abs(e) for e in E if e!=0]; bps={F(0),F(1)}
    for e in Es:
        for m in range(0,7*e+1): bps.add(F(m,7*e))
    return sorted(b for b in bps if 0<=b<1)

def avoid_measure(E, A):
    """meas{ x : no e in E has floor(7 frac(e x)) in A }, A subset of {0..6}."""
    Aset=set(A); pts=breakpoints(E)+[F(1)]; tot=F(0)
    for i in range(len(pts)-1):
        a,b=pts[i],pts[i+1]
        if b<=a: continue
        mid=(a+b)/2; hit=False
        for e in E:
            if int(((e*mid)%1)*7) in Aset: hit=True; break
        if not hit: tot+=(b-a)
    return tot

def p0_measure(E):
    pts=breakpoints(E)+[F(1)]; tot=F(0)
    for i in range(len(pts)-1):
        a,b=pts[i],pts[i+1]
        if b<=a: continue
        mid=(a+b)/2; occ=set(int(((e*mid)%1)*7) for e in E)
        if len(occ)==7: tot+=(b-a)
    return tot

def Sr(E, r):
    # sum over r-subsets A of {1..6} (sector 0 always hit) of meas{avoid A}
    return sum(avoid_measure(E, A) for A in combinations(range(1,7), r))

def Ly_k9(E):
    S1=Sr(E,1); S2=Sr(E,2); S3=Sr(E,3)
    return 1 - F(13,18)*S1 + F(4,9)*S2 - F(1,6)*S3, (S1,S2,S3)

def Ly_inf(k):
    # iid: meas{avoid A}=(1-r/7)^{k-1}; S_r^iid = C(6,r)(1-r/7)^{k-1}
    S1=6*F(6,7)**(k-1); S2=15*F(5,7)**(k-1); S3=20*F(4,7)**(k-1)
    return 1 - F(13,18)*S1 + F(4,9)*S2 - F(1,6)*S3

cap9=F(1979,4004)
print(f"cap_9={float(cap9):.5f}   L_y^inf(9)={float(Ly_inf(9)):.5f}\n")
fams={
 "consec_9           ":[0,1,2,3,4,5,6,7,8],
 "wide 3-cluster     ":[0,1,2,30,31,32,60,61,62],
 "wide 3-cluster x10 ":[0,1,2,100,101,102,200,201,202],
 "gapped {0..7,500}  ":[0,1,2,3,4,5,6,7,500],
 "2-scale {0..3,300} ":[0,1,2,3,300,301,302,303,600],
 "wide dissociated   ":[0,1,4,9,16,25,37,50,64],
}
print(f"{'family':>20} {'p0(E)':>8} {'L_y(E)':>8} {'L_y>=p0':>8} {'L_y<cap':>8} {'L_y-Lyinf':>10}")
lyinf=Ly_inf(9)
for name,E in fams.items():
    p0=p0_measure(E); ly,_=Ly_k9(E)
    print(f"{name:>20} {float(p0):>8.5f} {float(ly):>8.5f} {str(ly>=p0):>8} {str(ly<cap9):>8} {float(ly-lyinf):>+10.5f}")

print("\n=== multiplicative recursion: meas_{E' u w}{avoid A} / meas_{E'}{avoid A} ~ (1-r/7)? ===")
E=[0,1,2,30,31,32,60,61,62]; w=max(E); Ep=[e for e in E if e!=w]
for r in [1,2,3]:
    ratios=[]
    for A in list(combinations(range(1,7),r))[:6]:
        mE=avoid_measure(E,A); mEp=avoid_measure(Ep,A)
        if mEp>0: ratios.append(float(mE/mEp))
    avg=sum(ratios)/len(ratios) if ratios else float('nan')
    print(f"  r={r}: avg ratio = {avg:.4f}   target (1-r/7) = {float(1-F(r,7)):.4f}   (n={len(ratios)} subsets)")

print("\n=== L_y decorrelation under peeling (does L_y stay < cap and approach L_y^inf)? ===")
cur=E[:]
while len(cur)>=6:
    ly,(S1,S2,S3)=Ly_k9(cur) if len(cur)==9 else (None,(0,0,0))
    p0=p0_measure(cur)
    print(f"  {str(cur):<44} p0={float(p0):.5f}" + (f"  L_y={float(ly):.5f}" if ly is not None else ""))
    if len(cur)<=5: break
    cur=[e for e in cur if e!=max(cur)]
