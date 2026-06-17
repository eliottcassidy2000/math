#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Verified number-facts for the triangular/square/Fibonacci/prime ↔ Alcuin-boat synthesis.
kind-pasteur-2026-06-16-S2.  All facts below VERIFIED here; the tournament synthesis
(Φ₃ gates ideal-gas big-boat AND forbidden gaps; the boat dichotomy is a doubling threshold)
is in HYP-2553 + the reflection.  Builds on THM-519 (Alcuin big-boat ⟺ Ω edgeless ⟺ H=3^{α₁}),
S599v ({7,21} phantom volumes = Φ₃(2),3Φ₃(2)), THM-486 (Pisano), QR-resonance.
"""
import sys; sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None
from math import isqrt

def heron_area_sq16(a,b,c): return (a+b+c)*(-a+b+c)*(a-b+c)*(a+b-c)  # = 16·Area²

# (1) exactly 5 integer triangles with Area = Perimeter (m=1); all have inradius r=2
def area_eq_m_perimeter(m, amax=130, bmax=170):
    out=[]
    for a in range(1,amax):
        for b in range(a,bmax):
            for c in range(b,a+b):
                A2=heron_area_sq16(a,b,c)
                if A2<=0 or A2%16: continue
                A=isqrt(A2//16)
                if A*A!=A2//16: continue
                if A==m*(a+b+c): out.append((a,b,c,A,a+b+c, A/((a+b+c)/2)))  # last = inradius
    return out
print("(1) Area = 1·Perimeter Heronian triangles (a,b,c,Area,Perim,inradius):")
for t in area_eq_m_perimeter(1): print("   ",t)
print("    => exactly 5; ALL inradius r=2 (Area=r·s, Area=P=2s ⟹ r=2).")

# (2) Fib ∩ Tri = {1,3,21,55} = Fibonacci prefix of (2n-1)(4n-1)=T_{4n-2}
fib=set(); a,b=1,1
while a<10**13: fib.add(a); a,b=b,a+b
tri=set(k*(k+1)//2 for k in range(2*10**6))
print("(2) Fib ∩ Tri (>0):", sorted(t for t in tri if t in fib and t>0)[:10])
print("    (2n-1)(4n-1)=T_{4n-2}; factor as (2n-1)·(4n-1):")
for n in range(0,6):
    v=(2*n-1)*(4*n-1); m=4*n-2; T=m*(m+1)//2 if m>=0 else 1
    print(f"    n={n}: {2*n-1}·{4*n-1}={v}=T_{m}  Fib={v in fib}   (n=2: 3·7=21=Φ3(1)·Φ3(2))")

# (3) odd² = 8·T_k + 1
print("(3) (2k+1)² = 8·T_k + 1:", all((2*k+1)**2==8*(k*(k+1)//2)+1 for k in range(3000)),
      "(octal odd-square = octal T_k followed by digit 1)")

# (4) marble 2×6 grid up to its Z2×Z2 symmetry = Pólya/Burnside (same machinery as A000568)
cells=[(r,c) for r in range(2) for c in range(6)]
syms=[lambda p:p, lambda p:(p[0],5-p[1]), lambda p:(1-p[0],p[1]), lambda p:(1-p[0],5-p[1])]
def burnside(k):
    tot=0
    for g in syms:
        seen=set(); cyc=[]
        for p in cells:
            if p in seen: continue
            l=0;q=p
            while q not in seen: seen.add(q); q=g(q); l+=1
            cyc.append(l)
        dp=[0]*13; dp[0]=1
        for L in cyc:
            for s in range(12,L-1,-1): dp[s]+=dp[s-L]
        tot+=dp[k]
    return tot//4
print("(4) 2×6 grid, k marbles up to Z2×Z2, k=0..6:", [burnside(k) for k in range(7)],
      "→ {1,3,21,55} is the prefix; 135 breaks both Fib & Tri.")

# (5) cyclotomic Φ3 thread
print("(5) Φ3(x)=x²+x+1: Φ3(1)=3 (big-boat/ideal-gas unit H=3^α₁), Φ3(2)=7 (forbidden phantom),",
      "Φ3(1)·Φ3(2)=21 (forbidden, =T_6, =Fib∩Tri[4]).")

# (6) Fibonacci entry point α(p) | p−(5/p)  (Euler criterion = QR of 5 = Paley/THM-486 structure)
def leg(a,p):
    a%=p
    return 0 if a==0 else (-1 if pow(a,(p-1)//2,p)==p-1 else 1)
def rank_app(p):
    a,b=1,1; m=2
    while True:
        m+=1; a,b=b,(a+b)%p
        if b%p==0: return m
print("(6) Fibonacci rank-of-apparition α(p) | p−(5/p):")
for p in (2,3,5,7,11,13,17,19,23,29,31):
    r=rank_app(p); L=leg(5,p)
    ok = "p=5" if p==5 else str((p-L)%r==0)
    print(f"    p={p:2d}: α={r:2d}  (5/p)={L:+d}  p-(5/p)={p-L:2d}  α|p-(5/p)={ok}")
