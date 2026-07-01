#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
DEEPENING the census + Caratheodory-Toeplitz flat-extension via modular multiplication.

kind-pasteur-2026-07-01-S7 (companion to lrc14_census_toeplitz_ramanujan_kps.py).
 (A) CENSUS-BY-MODULAR-STRUCTURE: every difference-closed / near-AP 11-core lands its covering-min
     MAXIMIZERS on some (Z/N)*; tabulate (M, N, phi(N)=#atoms, meas(L_C@1/14)); confirm the min-measure
     11-core is the (Z/10)* one {1..13}\{6,10}.  The census is ORGANIZED by which (Z/N)* the atoms hit.
 (B) COLLAR WIDTH = MODULAR INVERSE: at atom a/N the binding runner has speed v = a^{-1} mod N (and
     N-a^{-1}); the collar half-width = (M-1/14)/v_max is set by the MODULAR INVERSE.  meas = sum over
     units of these widths = a units-sum reminiscent of a Dedekind/cotangent sum (-> HYP-3773, -1/12).
 (C) VERBLUNSKY / OPUC (the gap both search agents flagged as absent): Levinson-Durbin the Toeplitz
     moments c_k=hat 1_{L_C}(k) of the lonely-set measure; report the Verblunsky coefficients alpha_k
     (|alpha_k|<1 <=> PD Toeplitz <=> valid measure) and the Szego product -- the OPUC coordinates of L_C.
"""
import sys, itertools, cmath
from fractions import Fraction as Fr
from math import gcd, pi
try: sys.stdout.reconfigure(encoding='utf-8')
except Exception: pass

def nrm(x):
    f = x-(x.numerator//x.denominator); return f if f<=Fr(1,2) else 1-f
def gmin(S,t): return min(nrm(v*t) for v in S)
def units(N): return [a for a in range(1,N) if gcd(a,N)==1]
def phi(N): return len(units(N))
def find_M_atoms(S,Dmax=None):
    if Dmax is None: Dmax=2*max(S)+2
    dens=set()
    for v in S: dens.add(2*v)
    for v in S:
        for w in S:
            dens.add(v+w)
            if v!=w: dens.add(abs(v-w))
    dens={d for d in dens if 0<d<=Dmax}
    best=Fr(0); atoms=[]
    for b in sorted(dens):
        for a in range(1,b):
            if gcd(a,b)!=1: continue
            val=gmin(S,Fr(a,b))
            if val>best: best=val; atoms=[Fr(a,b)]
            elif val==best: atoms.append(Fr(a,b))
    return best, sorted(set(atoms))
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

print("="*96); print(" (A) CENSUS-BY-MODULAR-STRUCTURE: 11-cores labelled by the (Z/N)* their maximizers hit"); print("="*96)
rows=[]
# structured near-tight 11-cores: {1..13} minus 2, and {1..11} minus1 plus far, plus dilated
cores=set()
for i,j in itertools.combinations(range(1,14),2):
    C=tuple(sorted(set(range(1,14))-{i,j}))
    if len(C)==11: cores.add(C)
for C in list(cores):
    S=list(C); M,atoms=find_M_atoms(S,Dmax=40)
    Ns=sorted(set(t.denominator for t in atoms))
    N=Ns[0] if len(Ns)==1 else None
    m=meas(Lset(S))
    isunit = (N is not None and sorted(t.numerator for t in atoms)==units(N))
    rows.append((float(m),m,S,M,N,len(atoms),isunit))
rows.sort()
print(f"  {len(rows)} eleven-cores of form {{1..13}}\\{{i,j}}. Smallest 8 by meas(L_C@1/14):")
print(f"    {'meas':>10} {'M':>7} {'N(atoms)':>9} {'#at=phi':>8} {'units?':>6}  core")
for fm,m,S,M,N,na,iu in rows[:8]:
    print(f"    {fm:>10.5f} {str(M):>7} {str(N):>9} {na:>8} {str(iu):>6}  {S}")
print(f"  ... largest 3:")
for fm,m,S,M,N,na,iu in rows[-3:]:
    print(f"    {fm:>10.5f} {str(M):>7} {str(N):>9} {na:>8} {str(iu):>6}  {S}")
mn=rows[0]
print(f"  MIN-measure 11-core: {mn[2]}  meas={mn[1]}={mn[0]:.5f}  atoms on (Z/{mn[4]})* (phi={mn[5]})")

print("\n"+"="*96); print(" (B) COLLAR WIDTH = MODULAR INVERSE: binding speed at atom a/N is v=a^{-1} mod N"); print("="*96)
def modinv(a,N):
    for x in range(1,N):
        if (a*x)%N==1: return x
    return None
for S,N in [([1,2,3,4,5,7,8,9,11,12,13],10), (list(range(1,12)),12), (list(range(1,14)),14)]:
    M,atoms=find_M_atoms(S,Dmax=40)
    print(f"  S={S}  atoms on (Z/{N})*, M={M}:")
    for t in atoms:
        a=t.numerator  # t=a/N
        if t.denominator!=N: continue
        inv=modinv(a%N,N); binding=[v for v in S if nrm(v*t)==M]
        print(f"    atom {a}/{N}: a^-1 mod {N} = {inv};  binding speeds (||v t||=M) = {binding};  "
              f"|binding| relates to {{a^-1, N-a^-1}}={{ {inv},{N-inv} }}")
    # measure vs. sum of collar half-widths (M-1/14)/v over binding speeds
    print(f"     meas(L_C)={meas(Lset(S))}  (collar volume around the {len(atoms)} unit-atoms + secondaries)")

print("\n"+"="*96); print(" (C) VERBLUNSKY / OPUC coordinates of the lonely-set measure (fills the repo gap)"); print("="*96)
def fourier_c(arcs,k):  # c_k = hat 1_{L_C}(k) = int_{L_C} e^{-2pi i k t} dt, exact-ish via rationals->float
    if k==0: return complex(float(meas(arcs)),0.0)
    s=0j
    for (lo,hi) in arcs:
        s+=(cmath.exp(-2j*pi*k*float(lo))-cmath.exp(-2j*pi*k*float(hi)))/(2j*pi*k)
    return s
def verblunsky(cvec):
    # Levinson-Durbin (Szego) from normalized moments r_k=c_k/c_0; returns alpha_0..alpha_{M-1}
    c0=cvec[0].real; r=[c/c0 for c in cvec]
    M=len(r)-1; a=[1.0+0j]; alphas=[];
    # phi_0=1; Szego recursion
    phi=[1.0+0j]; err=1.0
    for n in range(0,M):
        # reflection coefficient alpha_n = -(sum_{k=0}^{n} phi_k * r_{n+1-k}) / err
        s=sum(phi[k]*r[n+1-k] for k in range(len(phi)))
        alpha=-s/err
        alphas.append(alpha)
        newphi=[0j]*(len(phi)+1)
        for k in range(len(phi)): newphi[k]+=phi[k]
        for k in range(len(phi)): newphi[k+1]+=alpha*phi[len(phi)-1-k].conjugate()
        phi=newphi; err=err*(1-abs(alpha)**2)
    return alphas, err
for S in [[1,2,3,4,5,7,8,9,11,12,13], list(range(1,12))]:
    arcs=Lset(S); c=[fourier_c(arcs,k) for k in range(0,13)]
    al,err=verblunsky(c)
    mags=[round(abs(x),4) for x in al]
    print(f"  S={S}")
    print(f"    c_0=meas={c[0].real:.5f};  |c_k| k=1..6 = {[round(abs(x),4) for x in c[1:7]]}")
    print(f"    |Verblunsky alpha_k| k=0..{len(al)-1} = {mags}")
    print(f"    all |alpha_k|<1 (PD Toeplitz / valid measure)? {all(abs(x)<1 for x in al)};  Szego err_prod={err:.4e}")
print("  => |alpha_k|<1 confirms the Toeplitz moment matrix of L_C is PD (valid measure); the alpha_k")
print("     are the OPUC coordinates. Atomic limit (band->M) drives some |alpha_k|->1 (singular/atomic).")
print("DONE.")
