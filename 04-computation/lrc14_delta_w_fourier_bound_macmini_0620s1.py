#!/usr/bin/env python3
"""
lrc14_delta_w_fourier_bound_macmini_0620s1.py  (mac-mini-2026-06-20-S1)

THE ONE OPEN CONSTANT, made rigorous: a convergent explicit bound on the far-element
deviation Delta_w = p0(E'∪{w}) - Phi(E'),  Phi(E')=p0(E')+(1/7)p1(E')  (HYP-2642 recursion).

DERIVATION (peel the single far element w => 1-D discrepancy, NOT the divergent multi-D lattice):
 p0(E'∪{w}) = p0(E') + Σ_{j=1..6} meas{B_j ∩ {frac(wx)∈sector_j}},  B_j={x: E' misses EXACTLY sector j}.
 (w is one point per x, fills at most one missed sector.)  Hence
   Delta_w = Σ_j [ meas{B_j ∩ {frac(wx)∈sector_j}} - (1/7)meas(B_j) ]
           = Σ_j Σ_{n≠0} shat_j(n) * 1hat_{B_j}(-n w),
 shat_j(n)=|sin(pi n/7)|/(pi|n|) magnitude (VANISHES at 7|n), 1hat_{B_j}(m): |·|<=#arcs(B_j)/(pi|m|).
 => |Delta_w| <= kappa * (Σ_j #arcs(B_j)) / (pi^2 * w) = C(E')/w,
    kappa = 2 Σ_{n>=1, 7∤n} |sin(pi n/7)|/n^2.

VERIFY (exact rational Delta_w vs the bound C(E')/w), and report kappa, C(E'), the implied
cutoff ratio R = C(E')/(margin*max(E')) [w > R*max(E') => |Delta_w|<margin].
"""
import sys, itertools, math
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
SEV=F(1,7)

def sector_of(p): return int((p%1)*7)
def measS7(E):
    E=sorted(set(E)); bps=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(0,7*e+1): bps.add(F(m,7*e))
    bps=sorted(bps); tot=F(0)
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2
        if len(set(sector_of(e*xm) for e in E))==7: tot+=x1-x0
    return tot
def p1_and_Bj_arcs(E):
    """p1(E)=meas{exactly 1 sector missed}; and #arcs of each B_j (E misses exactly sector j)."""
    E=sorted(set(E)); bps=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(0,7*e+1): bps.add(F(m,7*e))
    bps=sorted(bps)
    p1=F(0); inBj=[False]*7; arccount=[0]*7
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2; secs=set(sector_of(e*xm) for e in E)
        missed=[j for j in range(7) if j not in secs]
        if len(missed)==1:
            j=missed[0]; p1+=x1-x0
            if not inBj[j]: arccount[j]+=1; inBj[j]=True
            for jj in range(7):
                if jj!=j: inBj[jj]=False
        else:
            for jj in range(7): inBj[jj]=False
    return p1, arccount

def delta_w(Ep, w):
    E=sorted(set(list(Ep)+[w]))
    p0E=measS7(Ep); p0Ew=measS7(E); p1E,_=p1_and_Bj_arcs(Ep)
    return p0Ew - (p0E + SEV*p1E)

# kappa = 2 sum_{n>=1, 7 not | n} sin(pi n/7)/n^2  (numeric, high precision)
kappa=2*sum(abs(math.sin(math.pi*n/7))/n**2 for n in range(1,200000) if n%7!=0)
print(f"kappa = 2*sum_{{n>=1,7∤n}} |sin(pi n/7)|/n^2 = {kappa:.6f}")
print("="*92)
print(f"{'E_prime':<30}{'w':>5}{'Delta_w':>12}{'|Dw|*w':>10}{'sum#arcs':>9}{'C(Eprime)':>11}{'bound C/w':>11}{'OK?':>5}")
print("="*92)
cases=[
  ([0,1,2,3,4,5,6,7],[10,20,50,140]),         # consec core + far
  ([0,1,2,4,6,7,8,10],[12,24,60,120]),         # codex B13 leader-ish core
  ([0,3,5,16,28,30,33],[35,70,200]),           # third-pocket-ish core
  ([0,1,2,4,8,12,16,20],[24,48,120]),          # dyadic-block m=4 core
]
maxC=0.0
for Ep,ws in cases:
    p1,arcs=p1_and_Bj_arcs(Ep); V=sum(arcs); C=kappa*V/math.pi**2
    for w in ws:
        dw=delta_w(Ep,w); val=abs(float(dw))*w; ok = val<=C+1e-9
        print(f"{str(Ep):<30}{w:>5}{float(dw):>12.6f}{val:>10.4f}{V:>9}{C:>11.4f}{C/w:>11.5f}{'yes' if ok else 'NO!':>5}")
    maxC=max(maxC,C)
print("="*92)
print("If |Dw|*w <= C(Eprime) holds in every row, the bound |Delta_w|<=C(E')/w is verified.")
print("Cutoff: |Delta_w|<margin once w > C(E')/margin. k=9 margin=cap_9-Q(8)=129643/980980~0.132.")
print(f"  e.g. consec_8 core: C={kappa*sum(p1_and_Bj_arcs([0,1,2,3,4,5,6,7])[1])/math.pi**2:.3f}, "
      f"cutoff w>C/0.132 ~ {kappa*sum(p1_and_Bj_arcs([0,1,2,3,4,5,6,7])[1])/math.pi**2/0.132:.0f}")
print("\nDONE.")
