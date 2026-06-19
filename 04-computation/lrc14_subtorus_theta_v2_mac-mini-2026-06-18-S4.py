#!/usr/bin/env python3
"""
lrc14_subtorus_theta_v2 — mac-mini-2026-06-18-S4  (efficient: support-based lattice sum)

Verify codex's subtorus Fourier identity and the KEY connection (bounded-k escapes THM-504):
  int_0^1 G(Ex) dx = (5/7)^k + Σ_{n∈Λ_aff(E),n≠0} ∏ ψ̂(n_i),
  G(Ex)=Σ_gaps(gap-2/7)_+ (empty 2/7-window measure), ψ̂(0)=5/7, ψ̂(n)=-(1-e^{-4πin/7})/(2πin).
Correction dominated by SUPPORT-3 relations (every triple gives one primitive (a,b,c)∝
(e_j-e_l,e_l-e_i,e_i-e_j)), each summed over multiples m (tail ~ 1/m^3, CONVERGENT), plus support-4.
Bounded k=|E|≤13 ⟹ finitely many primitive triples (≤C(k,3)) ⟹ correction converges ABSOLUTELY,
unlike THM-504's speed-side |T|≥3 DIVERGENCE (unbounded speeds).
"""
import math, cmath, itertools
from math import gcd
from fractions import Fraction as F

def psihat(n):
    if n==0: return 5/7
    return -(1-cmath.exp(-2j*math.pi*2*n/7))/(2j*math.pi*n)

def intG_grid(E, Q=200000, c=2/7):
    es=sorted(set(E)); s=0.0
    for a in range(Q):
        x=a/Q; pts=sorted((e*x)%1 for e in es); aug=pts+[pts[0]+1]
        s+=sum((aug[i+1]-aug[i]-c) for i in range(len(pts)) if aug[i+1]-aug[i]>c)
    return s/Q

def prim(v):
    g=0
    for x in v: g=gcd(g,x)
    return tuple(x//g for x in v) if g else v

def correction(E, Mmult=200, c4bound=6):
    """support-3 (all triples x all multiples) + support-4 (bounded coeffs)."""
    es=list(E); k=len(es); total=0.0
    # support-3
    seen3=set()
    for i,j,l in itertools.combinations(range(k),3):
        a,b,cc = es[j]-es[l], es[l]-es[i], es[i]-es[j]
        if (a,b,cc)==(0,0,0): continue
        pa,pb,pc=prim((a,b,cc))
        key=(i,j,l,pa,pb,pc)
        if key in seen3: continue
        seen3.add(key)
        for m in range(1,Mmult+1):
            for sgn in (1,-1):
                t=psihat(sgn*m*pa)*psihat(sgn*m*pb)*psihat(sgn*m*pc)
                total+=t.real
    # support-4: for each 4-subset, enumerate small relations (Σn=0, Σn e=0), exact support 4
    for quad in itertools.combinations(range(k),4):
        ev=[es[q] for q in quad]
        for n in itertools.product(range(-c4bound,c4bound+1),repeat=4):
            if all(x!=0 for x in n) and sum(n)==0 and sum(n[r]*ev[r] for r in range(4))==0:
                t=1.0
                for x in n: t*=psihat(x)
                total+=t.real if isinstance(t,complex) else t
    return total

print("="*76)
print("VERIFY: (5/7)^k + correction(support 3+4)  vs  int G(Ex)dx")
print("="*76)
tests=[("AP k=5",[0,1,2,3,4]),("AP k=6",[0,1,2,3,4,5]),("AP k=7",[0,1,2,3,4,5,6]),
       ("dissoc k=5",[0,1,4,9,11]),("perforated k=7",[0,2,3,4,5,6,8]),
       ("dissoc k=6",[0,1,5,11,13,18])]
for name,E in tests:
    k=len(E); main=(5/7)**k; corr=correction(E); ig=intG_grid(E)
    print(f"  {name:16s}: (5/7)^{k}={main:.5f} + corr={corr:+.5f} = {main+corr:.5f}   intG(grid)={ig:.5f}   |corr|/main={abs(corr)/main:.2f}")
print()
print("="*76)
print("support-3 multiple-tail convergence (primitive (1,-2,1), the AP triple):")
print("="*76)
acc=0.0
for m in range(1,60):
    t=2*(psihat(m*1)*psihat(m*-2)*psihat(m*1)).real  # +-m
    acc+=t
    if m in (1,3,10,30,59): print(f"   m≤{m:2d}: partial={acc:+.6f}  |term_m|={abs(t):.2e}")
print("  => converges (~1/m^3). Finitely many primitive triples (k≤13) ⟹ absolute convergence,")
print("     ESCAPING THM-504's |T|≥3 divergence (which summed over UNBOUNDED speeds).")
print()
print("="*76)
print("FLOOR readout: is |correction| < (5/7)^k? (then int G>0 ⟹ μ>0 ⟹ M≥1/14)")
print("="*76)
for name,E in tests:
    k=len(E); main=(5/7)**k; corr=correction(E)
    ok = (main+corr)>0
    print(f"  {name:16s}: intG≈{main+corr:+.5f} {'>0 ✓' if ok else '≤0 ✗ (corr overwhelms main; need the hard-indicator μ, not this minorant)'}")
