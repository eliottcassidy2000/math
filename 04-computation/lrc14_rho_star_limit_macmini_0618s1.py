#!/usr/bin/env python3
"""
lrc14_rho_star_limit_macmini_0618s1.py  (mac-mini-2026-06-18-S1)

Validate the SINGLE-VARIABLE reformulation of kind-pasteur's good-period density
rho*(Delta,P) for the LRC(14) S3 k>=3 residual (THM-526 / OPEN-Q-108 / HYP-2581d).

CLUSTER L = {w0 + d_i : d_i in Delta}, Vmax = w0 + s (s=max Delta), co-offsets e_i = s - d_i.
P = small part (subset of 1..13). thresh = 1/(7 Vmax).

FINITE-w0 (kind-pasteur): rho_K = (#good Vmax-ruler periods)/Vmax, where period j is GOOD if
   the safe set of A=S\{Vmax} restricted to the ruler arc I_j=((14j+1)/(14Vmax),(14j+13)/(14Vmax))
   contains a sub-arc of width > thresh.

LIMIT (my claim): good_period j  <=>  x:=center(I_j) in G_P  AND  the k points
   {frac(e_i * x)} have max circular gap > 2/7  (e_k=0 from Vmax; teeth width 1/7 each).
   => rho*_limit = meas{ x in [0,1): x in G_P and maxgap{frac(e_i x)} > 2/7 }.

GOAL: show rho_K(w0) -> rho*_limit as w0->inf (validates the reformulation), for several (Delta,P).
"""
import sys
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
H = F(1,14)

def danger_arcs(u, h=H):
    iv=[]
    for j in range(u):
        c=F(j,u); a=(c-h/u)%1; b=(c+h/u)%1
        if a<b: iv.append((a,b))
        else: iv.append((a,F(1))); iv.append((F(0),b))
    return iv

def merge(iv):
    iv=sorted(iv); out=[]
    for a,b in iv:
        if out and a<=out[-1][1]: out[-1]=(out[-1][0],max(out[-1][1],b))
        else: out.append((a,b))
    return out

def safe_set(A, h=H):
    dz=merge([iv for u in A for iv in danger_arcs(u,h)])
    safe=[]; prev=F(0)
    for a,b in dz:
        if a>prev: safe.append((prev,a))
        prev=max(prev,b)
    if prev<1: safe.append((prev,F(1)))
    return safe

def meas(arcs): return float(sum(b-a for a,b in arcs))

def rho_K_finite(P, Delta, w0):
    """Efficient: precompute safe(A) once, then count ruler periods of Vmax that
    contain a safe(A) sub-arc of length > thresh."""
    s=max(Delta); Vmax=w0+s
    L=[w0+d for d in Delta]
    S=sorted(set(list(P)+L))
    A=[u for u in S if u!=Vmax]
    thresh=F(1,7*Vmax)
    safeA=safe_set(A)                     # sorted disjoint arcs, computed ONCE
    good=set()
    for (lo,hi) in safeA:
        # which ruler periods does this safe(A) arc overlap, and does it leave >thresh inside one?
        # ruler period j spans tau in [j/Vmax, (j+1)/Vmax); safe sub-window inside it is
        # I_j=((14j+1)/(14Vmax),(14j+13)/(14Vmax)). Intersect [lo,hi] with each I_j it touches.
        jlo=int(lo*Vmax); jhi=int(hi*Vmax)+1
        for j in range(max(0,jlo), min(Vmax, jhi+1)):
            ilo=F(14*j+1,14*Vmax); ihi=F(14*j+13,14*Vmax)
            x0=max(lo,ilo); x1=min(hi,ihi)
            if x1-x0 > thresh:
                good.add(j)
    return F(len(good),Vmax)

def maxgap(points):
    pts=sorted(p%1.0 for p in points)
    if not pts: return 1.0
    g=0.0
    for i in range(len(pts)-1):
        g=max(g, pts[i+1]-pts[i])
    g=max(g, (pts[0]+1.0)-pts[-1])
    return g

def rho_star_limit(P, Delta, N=200000):
    s=max(Delta); e=[s-d for d in Delta]
    GPf=[(float(a),float(b)) for a,b in safe_set(P)]
    cnt=0
    for t in range(N):
        x=(t+0.5)/N
        ok=False
        for a,b in GPf:
            if a<=x<b: ok=True; break
        if not ok: continue
        if maxgap([(ei*x)%1.0 for ei in e]) > 2.0/7.0:
            cnt+=1
    return cnt/N, meas(GPf)

CASES = [
  ("k=3 P={1..10} Delta={0,1,2}", list(range(1,11)), [0,1,2]),
  ("k=4 P={1..9}  Delta={0,1,2,3}", list(range(1,10)), [0,1,2,3]),
  ("k=4 P={1..9}  Delta={0,2,5,9}", list(range(1,10)), [0,2,5,9]),
  ("k=5 P={1..8}  Delta={0,1,2,3,4}", list(range(1,9)), [0,1,2,3,4]),
  ("k=6 P={1..7}  Delta={0,1,2,3,4,5}", list(range(1,8)), [0,1,2,3,4,5]),
]
print("="*86)
print("VALIDATE: rho_K(finite w0)  vs  rho*_limit (single-variable reformulation)")
print("="*86)
for name,P,Delta in CASES:
    rl, gpm = rho_star_limit(P,Delta, N=120000)
    print(f"\n{name}:  meas(G_P)={gpm:.5f}")
    print(f"    rho*_limit (my reformulation) = {rl:.5f}")
    for w0 in (101, 401, 1009, 3001):
        rk = float(rho_K_finite(P,Delta,w0))
        print(f"    rho_K(w0={w0:5d}) = {rk:.5f}   (Vmax={w0+max(Delta)})")
print("\nDONE.")
