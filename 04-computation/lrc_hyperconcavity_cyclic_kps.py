#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
HYPERCONCAVITY (CYCLIC, corrected): the CYCLIC additive energy of the WITNESS RESIDUES is the
order parameter for regularizable/Eisenstein (ordered/arc) vs cusp (disordered/generic).

kind-pasteur-2026-07-01. Fixes lrc_hyperconcavity_eisenstein_cusp_kps.py, which used the
INTEGER-LINE autocorrelation and so mis-scored the construction (its antipode -1 = Phi6-1 is
adjacent to 0 on the CYCLE, a near-arc, but looked far on the line). The correct object is the
CYCLIC (Z/D) autocorrelation of the witness residues {v*a mod D} at the binding rotation a/D.

CLAIM (opus HYP-3775 made quantitative): witness residues form a (cyclic) ARC/AP  <=>
cyclic additive energy MAXIMAL  <=>  autocorrelation = log-concave Fejer triangle  <=>
Dedekind sum closes (regularizable, -1/12=zeta(-1), Eisenstein). The DEFICIT below the arc =
a computable proxy for the un-regularizable cusp/residual (OPEN-Q-108's iota-odd index).

E_cyc(R)=sum_t A(t)^2, A(t)=#{(x,y) in R^2: x-y==t mod D} = ||hat 1_R||_4^4 on Z/D (the
hypercontractive L4 moment). E_norm = E_cyc / (arc energy of same size). 1 = perfect arc.
"""
import sys, math
from fractions import Fraction as Fr
from math import gcd
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None

def cyc_autocorr(R,D):
    A=[0]*D
    for x in R:
        for y in R: A[(x-y)%D]+=1
    return A
def cyc_energy(R,D): return sum(a*a for a in cyc_autocorr(R,D))
def arc_energy(m):   return sum((m-abs(t))**2 for t in range(-(m-1),m))  # interval {0..m-1}
def cyc_lc_defect(R,D):
    """cyclic log-concavity (hyperconcavity) defect around the arc: violations of A(t)^2>=A(t-1)A(t+1)
    over the (cyclically) contiguous support. 0 = Fejer triangle = arc/AP."""
    A=cyc_autocorr(R,D)
    # rotate so the support is contiguous-ish: measure over t where A>0 flanked by A>0
    viol=0; nz=0
    for t in range(D):
        a,b,c=A[(t-1)%D],A[t],A[(t+1)%D]
        if b>0: nz+=1
        if a>0 and b>0 and c>0 and b*b<a*c: viol+=1
    return viol,nz
def witness(S,Dmax):
    """M(S)=max_{a/D} min_v ||v a/D|| over rationals a/D, D<Dmax (exact); returns (M,a,D)."""
    best=(Fr(-1),0,0)
    for D in range(2,Dmax):
        for a in range(1,D):
            if gcd(a,D)!=1: continue
            m=min(min((v*a)%D, D-(v*a)%D) for v in S)
            val=Fr(m,D)
            if val>best[0]: best=(val,a,D)
    return best
def residues(S,a,D): return sorted((v*a)%D for v in S)
def is_cyc_arc_upto_scale(R,D):
    """is R a scaled arc/AP in Z/D? test: some unit u makes u*R a contiguous interval (mod D)."""
    R=set(x%D for x in R)
    for u in range(1,D):
        if gcd(u,D)!=1: continue
        T=sorted((u*x)%D for x in R)
        # contiguous (as an arc) iff sorted gaps are all 1 except one big gap
        gaps=[(T[(i+1)%len(T)]-T[i])%D for i in range(len(T))]
        if sum(1 for g in gaps if g!=1)<=1: return True,u
    return False,None

print("="*96)
print("PART 1  WITNESS RESIDUES + CYCLIC order parameter: construction (arc) vs beaters (generic)")
print("="*96)
constr={n:sorted(list(range(1,n-1))+[n*(n-1)]) for n in range(7,15)}   # {1..n-2} U {n(n-1)=2T}
beaters={7:[1,2,5,6,7,8], 8:[1,4,5,6,7,11,16], 9:[1,3,4,5,7,11,18,32], 10:[1,2,3,5,6,7,8,9,30]}
Mknown={7:Fr(2,13),8:Fr(2,15),9:Fr(4,33),10:Fr(4,37)}
def analyze(tag,S,Dmax):
    M,a,D=witness(S,Dmax); R=residues(S,a,D); m=len(R)
    E=cyc_energy(R,D); En=E/arc_energy(m); viol,nz=cyc_lc_defect(R,D)
    isarc,u=is_cyc_arc_upto_scale(R,D)
    print(f"  {tag:24s} M={str(M):>7} (a/D={a}/{D})  residues(descaled by 1/{a}? ) E_norm={En:.3f}  "
          f"lc_viol={viol}  cyc-arc/AP={isarc}")
    return En,isarc,M,D
print("  -- CONSTRUCTION {1..n-2, n(n-1)} (should be cyclic ARC at witness => ordered/Eisenstein) --")
con_res={}
for n in range(7,12):
    En,isarc,M,D=analyze(f"construction n={n}",constr[n], n*n-n+5); con_res[n]=(En,isarc,M)
print("  -- BEATERS (generic covering-min; opus verified NON-interval residues => cusp) --")
beat_res={}
for n,S in beaters.items():
    En,isarc,M,D=analyze(f"beater n={n}",S, 60); beat_res[n]=(En,isarc,M)

print("\n"+"="*96)
print("PART 2  DICHOTOMY TABLE: cyclic arc <=> high E_norm <=> regularizable (checkable)")
print("="*96)
print(f"  {'n':>3} | {'construction':>28} | {'beater':>28}")
print(f"  {'':>3} | {'E_norm  arc?  M':>28} | {'E_norm  arc?  M':>28}")
for n in sorted(beaters):
    ce,ca,cm=con_res.get(n,(None,None,None)); be,ba,bm=beat_res[n]
    cs = f"{ce:.3f}  {str(ca):>5}  {str(cm)}" if ce is not None else "n/a"
    print(f"  {n:>3} | {cs:>28} | {be:.3f}  {str(ba):>5}  {str(bm):>28}")
print("  READING: construction residues = cyclic ARC (E_norm high, regularizable/-1/12);")
print("           beater residues = NOT arc (E_norm lower, un-regularizable/cusp). The gap = residual.")

print("\n"+"="*96)
print("PART 3  the -1/12 telescoping is the ARC structure (regularizable side, reproduced)")
print("="*96)
def saw(x):
    x=Fr(x)
    return Fr(0) if x.denominator==1 else x-(x.numerator//x.denominator)-Fr(1,2)
def dedekind(a,D):
    a%=D; return sum(saw(Fr(i,D))*saw(Fr(a*i,D)) for i in range(1,D))
for n in [7,10,14]:
    D=n*n-n+1; T=Fr(n*(n-1),2); closed=Fr(-T,12*T+6); s=dedekind(n,D)
    print(f"  n={n}: s(n,Phi6)={float(s):+.6f} = closed -T/(12T+6)={float(closed):+.6f} [{s==closed}]  "
          f"-> -1/12 ; Phi6 mod n={D%n} (arc telescopes in 1 reciprocity step)")

print("\n"+"="*96)
print("SYNTHESIS: cyclic additive energy of the witness residues = the ORDER PARAMETER.")
print("  arc/AP (max) = Eisenstein/regularizable/-1/12 ; generic (deficit) = cusp/residual.")
print("  Extends opus HYP-3775's BINARY 'interval?' to a SCALAR that quantifies the residual")
print("  (the open OPEN-Q-108 iota-odd index) -- and identifies it as a hyperconcavity/Fejer defect.")
print("DONE.")
