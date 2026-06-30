#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""What the genus represents: Eisenstein/cusp dims + the Z_p core-gap landscape (klein-S11).

Weight-2 on Gamma_0(N): dim M_2 = dim Eisenstein + dim S_2 (cusp forms).
  dim Eisenstein_2(Gamma_0(N)) = (#cusps) - 1   [the BULK / boundary space]
  dim S_2(Gamma_0(N))          = genus(X_0(N))  [the OBSTRUCTION / global space]
Claim: across the LRC(2p) family the BULK is constant (=3, since 4 cusps) while the OBSTRUCTION = genus
grows 0,0,1,2,2 -- the genus IS the dimension of the missing (global, cusp-form) modes.

Also: the Z_p core-gap landscape gap(O)=min_{k!=0}|sum_{x in O} w^{kx}|^2 (the apex cyclotomic gap),
its distinct values / min / complement-symmetry, across p, to see the qualitative shift.
"""
from math import gcd, cos, pi, prod
import itertools, cmath

def pf(n):
    f=[];d=2;m=n
    while d*d<=m:
        while m%d==0:f.append(d);m//=d
        d+=1
    if m>1:f.append(m)
    return f
def psi(n):
    r=n
    for p in set(pf(n)): r=r*(p+1)//p
    return r
def phi(m):
    r=m
    for p in set(pf(m)): r=r//p*(p-1)
    return r
def divisors(n): return [d for d in range(1,n+1) if n%d==0]
def ncusps(n): return sum(phi(gcd(d,n//d)) for d in divisors(n))
def nu2(n): return 0 if n%4==0 else sum(1 for x in range(n) if (x*x+1)%n==0)
def nu3(n): return 0 if n%9==0 else sum(1 for x in range(n) if (x*x+x+1)%n==0)
def genus(n): return round(1+psi(n)/12-nu2(n)/4-nu3(n)/3-ncusps(n)/2)

print("="*80)
print(" Weight-2 dims on Gamma_0(2p): BULK (Eisenstein) vs OBSTRUCTION (cusp forms = genus)")
print("="*80)
print(f" {'N=2p':>5} {'apex':>4} {'#cusps':>7} {'dimEis=cusps-1':>14} {'dimCusp=genus':>14} {'note':>22}")
notes={6:'LRC(6) SOLVED',10:'apex5 not Mers/Heeg',14:'LRC(14) FIRST HARD',22:'apex11',26:'apex13'}
for N in [6,10,14,22,26]:
    p=max(pf(N)); c=ncusps(N); g=genus(N)
    print(f" {N:>5} {p:>4} {c:>7} {c-1:>14} {g:>14} {notes.get(N,''):>22}")
print(" => BULK/Eisenstein is CONSTANT (=3, four cusps) across the family; OBSTRUCTION/cusp = genus")
print("    grows 0,0,1,2,2. The genus is the DIMENSION of the missing global (cusp-form) modes.")

print("\n"+"="*80)
print(" Z_p core-gap landscape  gap(O)=min_{k!=0}|sum_{x in O} w_p^{kx}|^2  (the apex cyclotomic gap)")
print("="*80)
for p in [3,5,7,11,13]:
    w=cmath.exp(2j*cmath.pi/p)
    vals=set()
    minnz=1e9; argmin=None
    for sz in range(1,p+1):
        for O in itertools.combinations(range(p),sz):
            g=min(abs(sum(w**((k*x)%p) for x in O))**2 for k in range(1,p))
            gr=round(g,5); vals.add(gr)
            if gr>1e-9 and gr<minnz: minnz=gr; argmin=O
    doublet=2+2*cos(2*pi*((p-1)//2)/p)  # min doublet gap ~ angle nearest pi
    # exact min doublet: min over d of 2+2cos(2 pi * round(p/2)*?), just min over k of 2+2cos(2pi k/p)
    dmin=min(2+2*cos(2*pi*k/p) for k in range(1,p))
    print(f"  p={p:>2} (X_0(2p) genus {genus(2*p)}): {len(vals)} distinct gap values; "
          f"min nonzero={minnz:.5f} at O={argmin}; doublet gap={dmin:.5f}")
    print(f"        distinct gaps = {sorted(vals)}")
print("\n NOTE: the doublet gap = min_k 2+2cos(2pi k/p) DECREASES with p (0.198 at 7) independent of genus;")
print(" the genus distinction is GLOBAL (the cusp-form/obstruction dim), not the local gap-landscape min.")
print(" Cross-check (persistence test for hardness=genus): genus 0 vs 1 is about the GLOBAL obstruction.")
