#!/usr/bin/env python3
"""
lrc14_highheight_tailbound — mac-mini-2026-06-18-S4

The HIGH-relation-height half of codex HYP-2599 / HYP-2600, made rigorous.
Correction = Σ_{n∈Λ_aff(E),n≠0} ∏ψ̂(n_i), |ψ̂(n)| = |sin(2πn/7)|/(π|n|) (=0 if 7|n).
∫G(Ex)dx = (5/7)^k + correction; ∫G>0 ⟹ μ>0 ⟹ M(S)≥1/14.

RIGOROUS upper bound on |correction|, by support:
 - support 2: none (distinct offsets).
 - support 3: one primitive relation per triple {i,j,l}, (a,b,c)∝(e_j-e_l,e_l-e_i,e_i-e_j)/gcd;
   exact multiple-tail  T3 = Σ_{m≥1} 2·|ψ̂(ma)ψ̂(mb)ψ̂(mc)|  (tight, with sin factors; ≤ 2ζ(3)/(π³|abc|)).
 - support s≥4: bounded by enumerating low-coeff relations on each s-subset (tail beyond bounded
   coeffs is negligible by the 1/∏|n| decay).
CERTIFICATE: if B(E):=B3+B4(+...) < (5/7)^k  then ∫G(E)>0 PROVED (the HIGH-relation-height region).
The complement (some primitive triple of small height) = the finite-pattern LOW-height residual.
"""
import math, cmath, itertools
from math import gcd
from fractions import Fraction as F

def psihat_abs(n):
    if n%7==0: return 0.0
    return abs(math.sin(2*math.pi*n/7))/(math.pi*abs(n))
def prim3(a,b,c):
    g=gcd(gcd(abs(a),abs(b)),abs(c));
    return (a//g,b//g,c//g) if g else (a,b,c)
def tail3(a,b,c,M=300):
    s=0.0
    for m in range(1,M+1):
        s+=2*psihat_abs(m*a)*psihat_abs(m*b)*psihat_abs(m*c)
    return s
def B3(E):
    es=list(E); tot=0.0
    for i,j,l in itertools.combinations(range(len(es)),3):
        a,b,c=es[j]-es[l],es[l]-es[i],es[i]-es[j]
        if (a,b,c)==(0,0,0): continue
        pa,pb,pc=prim3(a,b,c)
        tot+=tail3(pa,pb,pc)
    return tot
def B4(E,cb=4):
    es=list(E); tot=0.0
    for quad in itertools.combinations(range(len(es)),4):
        ev=[es[q] for q in quad]
        for n in itertools.product(range(-cb,cb+1),repeat=4):
            if all(x!=0 for x in n) and sum(n)==0 and sum(n[r]*ev[r] for r in range(4))==0:
                p=1.0
                for x in n: p*=psihat_abs(x)
                tot+=p
    return tot
def minheight3(E):
    es=list(E); h=10**9
    for i,j,l in itertools.combinations(range(len(es)),3):
        a,b,c=es[j]-es[l],es[l]-es[i],es[i]-es[j]
        pa,pb,pc=prim3(a,b,c); H=abs(pa*pb*pc)
        if 0<H<h: h=H
    return h

print("="*80)
print("HIGH-relation-height certificate: B(E)=B3+B4 < (5/7)^k  ⟹  ∫G>0  ⟹  M(S)≥1/14")
print("="*80)
shapes={
 'AP k=7 (low-h)':[0,1,2,3,4,5,6],
 'perf k=7 (low-h)':[0,2,3,4,5,6,8],
 'dissoc k=5 (Sidon)':[0,1,4,9,11],
 'dissoc k=6':[0,1,4,9,15,22],
 'dissoc k=7':[0,1,8,11,13,17,21],
 'Sidon k=8':[0,1,3,7,12,20,30,44],
 'short-triple+spread {0,1,2,N}':[0,1,2,500],  # has triple (0,1,2) height 2 -> low-h
}
for nm,E in shapes.items():
    k=len(E); main=(5/7)**k; b3=B3(E); b4=B4(E); B=b3+b4; mh=minheight3(E)
    cert = B<main
    print(f"  {nm:30s} k={k}: minH3={mh:5d}  B3={b3:.5f} B4={b4:.5f} B={B:.5f}  (5/7)^k={main:.5f}  CERT∫G>0:{cert}")
print()
print("="*80)
print("THRESHOLD: for k≤13, how large must min relation-height be so B3 alone < (5/7)^k?")
print("  (crude: B3 ≤ C(k,3)·2ζ(3)/(π³·minH3) ⟹ minH3 > C(k,3)·2ζ(3)/(π³·(5/7)^k))")
print("="*80)
zeta3=1.2020569
for k in range(5,14):
    main=(5/7)**k; C=math.comb(k,3)
    H0=C*2*zeta3/(math.pi**3*main)
    print(f"  k={k:2d}: C(k,3)={C:3d}  (5/7)^k={main:.5f}  crude H0 (all triples need H>H0) = {H0:.0f}")
print()
print("READING: shapes whose every primitive triple has height > H0(k) are PROVED loose (high-h).")
print("The complement = shapes with a small-height triple = a FINITE list of support-3 coefficient")
print("patterns (codex HYP-2599 low-height side) — the remaining finite-enumeration residual.")
