#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
CENSUS + CARATHEODORY-TOEPLITZ FLAT-EXTENSION through the RIEMANN-ZETA / MODULAR-MULTIPLICATION lens.

kind-pasteur-2026-07-01-S7. Thesis to test computationally:
 (1) For the AP core {1..N-1}, the covering-min MAXIMIZERS (= the flat-extension ATOMS, HYP-3789)
     are exactly {a/N : a in (Z/N)*} -- PURE MODULAR MULTIPLICATION, count = phi(N).
 (2) The flat-extension moment matrix of that atomic measure is the RAMANUJAN-SUM Toeplitz matrix
     [c_N(j-l)], c_N(k)=sum_{a in (Z/N)*} e(ak/N).  So the Caratheodory-Toeplitz certificate for the
     AP lonely set IS a Ramanujan-sum matrix: rank = phi(N) (flat), spectrum = the divisor structure.
 (3) RAMANUJAN <-> ZETA: sum_{N>=1} c_N(k)/N^s = sigma_{1-s}(k)/zeta(s).  This is the exact place the
     Riemann zeta function enters -- the unit-restricted (modular) density THM-503(4) left UNTESTED.
     THM-503 proved L (archimedean) has NO Euler product; but the UNIT-RESTRICTED atomic moments DO
     (c_N multiplicative) and carry 1/zeta(s). Refines THM-503 without contradicting it.
 (4) The census target meas(L_C) at band 1/14 = collar volume around the phi(N) atoms; each collar
     width is set by the binding slopes = a cotangent/sawtooth (Dedekind) term.
 (5) The near-tight 11-core extremizer {1..13}\{6,10}: its atoms/denominators = a PERTURBED (Z/N)*.
"""
import sys, itertools
from fractions import Fraction as Fr
from math import gcd, cos, pi, prod
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass

def nrm(x):  # ||x|| distance to nearest integer, exact for Fraction
    f = x - (x.numerator//x.denominator)
    return f if f <= Fr(1,2) else 1-f
def g(S, t):  # min_v ||v t||
    return min(nrm(v*t) for v in S)
def totient_units(N): return [a for a in range(1,N) if gcd(a,N)==1]

def find_M_and_atoms(S, Dmax=None):
    """M(S)=max_t g and its argmax rationals. Maximizers sit at t=a/b with b | (v+-v') or 2v."""
    if Dmax is None: Dmax = 2*max(S)+2
    dens=set()
    for v in S: dens.add(2*v)
    for v in S:
        for w in S:
            dens.add(v+w);
            if v!=w: dens.add(abs(v-w))
    dens={d for d in dens if d>0 and d<=Dmax}
    best=Fr(0); atoms=[]
    for b in sorted(dens):
        for a in range(1,b):
            if gcd(a,b)!=1: continue
            t=Fr(a,b); val=g(S,t)
            if val>best: best=val; atoms=[t]
            elif val==best: atoms.append(t)
    atoms=sorted(set(atoms))
    return best, atoms

def ramanujan(N,k):  # c_N(k) exact integer via multiplicativity-free direct sum of cos
    s=0
    for a in totient_units(N):
        # e(ak/N)+conj is real; sum over units is real => sum cos(2 pi a k /N)
        s+=cos(2*pi*a*k/N)
    return s
def ramanujan_exact(N,k):  # exact via c_N(k)=mu(N/d) phi(N)/phi(N/d), d=gcd(N,k)
    from sympy import mobius, totient
    d=gcd(N,k); return int(mobius(N//d)*totient(N)//totient(N//d))

print("="*94); print(" PART 1  AP CORE MAXIMIZERS = (Z/N)*  (modular multiplication; #atoms=phi(N))"); print("="*94)
for N in [12,13,14]:
    S=list(range(1,N)); M,atoms=find_M_and_atoms(S)
    denoms=sorted(set(t.denominator for t in atoms))
    units=totient_units(N)
    match = all(t.denominator==N for t in atoms) and sorted(t.numerator for t in atoms)==units
    print(f"  S={{1..{N-1}}}: M={M}={float(M):.5f}; #atoms={len(atoms)}; denoms={denoms}")
    print(f"     atoms={[str(t) for t in atoms]}")
    print(f"     (Z/{N})* = {units} (phi={len(units)});  atoms == units/{N}?  {match}")

print("\n"+"="*94); print(" PART 2  FLAT-EXTENSION MOMENT MATRIX = RAMANUJAN-SUM TOEPLITZ  (rank=phi(N))"); print("="*94)
import numpy as np
for N in [12,14]:
    units=totient_units(N); phi=len(units)
    # moments of atomic measure mu=sum_{a in units} delta_{a/N}: muhat(k)=sum e(-a k/N)=c_N(k)
    def muhat(k): return sum(np.exp(-2j*pi*a*k/N) for a in units)
    for W in [phi, phi+2, phi+4]:
        T=np.array([[muhat(j-l) for l in range(W)] for j in range(W)])
        ev=np.linalg.eigvalsh((T+T.conj().T)/2)
        rk=int(np.sum(np.abs(ev)>1e-8))
        print(f"  N={N} phi={phi}: (W={W}) Toeplitz rank={rk}  (flat if =phi)  min|nonzero eig|={min(abs(e) for e in ev if abs(e)>1e-8):.3f}")
    # confirm muhat(k)=c_N(k) Ramanujan sum
    cvals=[round(muhat(k).real,6) for k in range(N)]
    rvals=[ramanujan_exact(N,k) for k in range(N)]
    print(f"     muhat(k) k=0..{N-1} = {[int(round(x)) for x in cvals]}")
    print(f"     c_{N}(k)  Ramanujan = {rvals}   match={all(abs(cvals[k]-rvals[k])<1e-6 for k in range(N))}")

print("\n"+"="*94); print(" PART 3  RAMANUJAN <-> ZETA:  sum_N c_N(k) N^-s = sigma_{1-s}(k)/zeta(s)"); print("="*94)
def zeta(s,T=200000): return sum(1.0/n**s for n in range(1,T))
def sigma(a,k):  # sigma_a(k)=sum_{d|k} d^a
    return sum(d**a for d in range(1,k+1) if k%d==0)
for s in [2.0, 3.0]:
    z=zeta(s)
    for k in [1,6,12]:
        lhs=sum(ramanujan_exact(N,k)/N**s for N in range(1,4000))
        rhs=sigma(1-s,k)/z  # sigma_{1-s}(k)=sum_{d|k} d^{1-s}
        print(f"  s={s} k={k:2d}:  sum_N c_N(k)/N^s = {lhs:.6f}   sigma_{{1-s}}(k)/zeta(s) = {rhs:.6f}   match={abs(lhs-rhs)<2e-3}")
print("  => the Riemann zeta function enters as the DENOMINATOR of the Ramanujan (unit-restricted) Dirichlet")
print("     series. THM-503: L(archimedean) has NO Euler product; the UNIT-RESTRICTED atomic moments DO.")

print("\n"+"="*94); print(" PART 4  NEAR-TIGHT 11-CORE EXTREMIZER {1..13}\\{6,10}: perturbed (Z/N)* atoms"); print("="*94)
BAND=Fr(1,14)
def safe_arcs(v): return [((Fr(k)+BAND)/v,(Fr(k+1)-BAND)/v) for k in range(v)]
def inter(A,B):
    r=[];i=j=0
    while i<len(A) and j<len(B):
        lo=max(A[i][0],B[j][0]); hi=min(A[i][1],B[j][1])
        if lo<hi: r.append((lo,hi))
        if A[i][1]<B[j][1]: i+=1
        else: j+=1
    return r
def Lset(S):
    a=safe_arcs(S[0])
    for v in S[1:]:
        a=inter(a,safe_arcs(v))
        if not a: return a
    return a
def meas(A): return sum(h-l for l,h in A)
for S in [[1,2,3,4,5,7,8,9,11,12,13], list(range(1,12)), [1,2,3,4,5,6,7,8,9,10,11,12]]:
    M,atoms=find_M_and_atoms(S, Dmax=40)
    Lc=Lset(S); m=meas(Lc)
    denoms=sorted(set(t.denominator for t in atoms))
    print(f"  S={S}")
    print(f"     M={M}={float(M):.5f} (>1/14={float(M)>1/14}); #atoms={len(atoms)} denoms={denoms}")
    print(f"     atoms={[str(t) for t in atoms]}")
    print(f"     meas(L_C @1/14)={m}={float(m):.5f}; #collars={len(Lc)}")

print("\n"+"="*94); print(" PART 5  CENSUS: collar volume = sum over atoms of width set by binding slopes"); print("="*94)
S=[1,2,3,4,5,7,8,9,11,12,13]; Lc=Lset(S)
print(f"  extremizer collars (lo,hi,width) around each maximizer:")
for (lo,hi) in Lc:
    mid=(lo+hi)/2
    binding=[v for v in S if nrm(v*mid)<Fr(1,7)+Fr(1,100)]  # near-binding runners
    print(f"    [{str(lo):>8},{str(hi):>8}] w={float(hi-lo):.5f}  center~{float(mid):.4f}")
print(f"  TOTAL meas={float(meas(Lc)):.5f}=313/9702 ; 1/36={float(Fr(1,36)):.5f} ; margin={float(meas(Lc)*36):.3f}x")
print("DONE.")
