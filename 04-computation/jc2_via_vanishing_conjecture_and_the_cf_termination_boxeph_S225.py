#!/usr/bin/env python3
"""jc2_via_vanishing_conjecture_and_the_cf_termination_boxeph_S225.py -- boxeph-2026-07-21-S225

HISTORICAL / CORRECTED: read MISTAKE-237 and HYP-8905.  Only two independent
checks survive here: exact Laplacian/Hessian calculations for the monomials
P=(x+iy)^d, and a finite illustration of Lame's Euclidean worst case.  Neither
check proves a descent for JC(2), and neither transfers NC2, GMC, or LRC to the
Jacobian conjecture.  THM-1300 refutes JC(n) for n>=3; JC(2) remains open.

The repo's general symmetric/Vanishing-Conjecture route for a planar Keller
map changes the ambient problem to four variables; the binary calculation
below is only a small subcase.  The four-variable VC route, planar
leading-form/Jelonek geometry, and Newton/continued-fraction descent are three
separate programs until an explicit map and decreasing invariant are proved.
"""
from math import gcd
from fractions import Fraction as F

def sep(t): print("\n"+"="*72+"\n"+t+"\n"+"="*72)
# polynomials in x,y as dict {(i,j): coeff(complex)}
def pt_add(a,b):
    r=dict(a)
    for k,v in b.items(): r[k]=r.get(k,0)+v
    return {k:v for k,v in r.items() if v!=0}
def pt_mul(a,b):
    r={}
    for (i1,j1),c1 in a.items():
        for (i2,j2),c2 in b.items():
            k=(i1+i2,j1+j2); r[k]=r.get(k,0)+c1*c2
    return {k:v for k,v in r.items() if v!=0}
def pt_pow(a,m):
    r={(0,0):1}
    for _ in range(m): r=pt_mul(r,a)
    return r
def dxx(p): return {(i-2,j):c*i*(i-1) for (i,j),c in p.items() if i>=2}
def dyy(p): return {(i,j-2):c*j*(j-1) for (i,j),c in p.items() if j>=2}
def lap(p): return pt_add(dxx(p),dyy(p))
def is_zero(p): return all(abs(v)<1e-9 for v in p.values()) if p else True

print("HISTORICAL / CORRECTED: read MISTAKE-237 and HYP-8905; the three programs are separate.")

# ==========================================================================
sep("P1  exact monomial check: P=(x+iy)^d is harmonic with nilpotent Hessian")
for d in (2,3,4,5):
    P={(d-m,m):(1j)**m*__import__('math').comb(d,m) for m in range(d+1)}  # (x+iy)^d
    lP=lap(P)
    # Hessian: [[Pxx,Pxy],[Pxy,Pyy]]; trace=Delta P; det=Pxx*Pyy - Pxy^2
    Pxx=dxx(P); Pyy=dyy(P)
    Pxy={(i-1,j-1):c*i*j for (i,j),c in P.items() if i>=1 and j>=1}
    det=pt_add(pt_mul(Pxx,Pyy), {k:-v for k,v in pt_mul(Pxy,Pxy).items()})
    print(f"  P=(x+iy)^{d}: Delta P = 0 (harmonic)? {is_zero(lP)} ; trace(Hess)=Delta P =0 and det(Hess)=0? {is_zero(det)} -> NILPOTENT Hessian")
# VC: Delta^m(P^m) = 0 for P=(x+iy)^d
for d,m in [(3,2),(3,3),(2,4),(4,2)]:
    P={(d-k,k):(1j)**k*__import__('math').comb(d,k) for k in range(d+1)}
    q=pt_pow(P,m)
    for _ in range(m): q=lap(q)
    print(f"  P=(x+iy)^{d}: Delta^{m}(P^{m}) = 0 (VC holds, invertible map)? {is_zero(q)}")
print("  => these one-sided monomials pass the displayed identities; this is not the full symmetric case.")

# ==========================================================================
sep("P2  dimension warning: the general symmetric route is not binary")
print("""  The relevant reduction for a general planar Keller map changes the ambient target to four variables.
  Therefore the binary monomial calculation above does not settle that target or JC(2).
  Read the scoped primary-source cards in CORE-PAPERS.md and the correction HYP-8905.""")

# ==========================================================================
sep("P3  analogy boundary: similar moment syntax is not a transfer")
print("""  Laplacian iterates, Gaussian moments, constant terms, and LRC Fourier fibers are differently typed.
  No DvdK bypass or coprime-interval theorem has been mapped to the four-variable VC target.
  A future bridge must state the map, preserved predicate, lost data, and a hostile example.""")

# ==========================================================================
sep("P4  independent toy: Euclidean descent has Fibonacci worst cases")
def euclid_steps(a,b):
    s=0
    while b: a,b=b,a%b; s+=1
    return s
# Lame: the worst-case (most steps) Euclidean pair below N is consecutive Fibonacci
best=(0,None)
N=200
for a in range(1,N):
    for b in range(1,a+1):
        st=euclid_steps(a,b)
        if st>best[0]: best=(st,(a,b))
fib=[1,1]
while fib[-1]<N: fib.append(fib[-1]+fib[-2])
print(f"  max Euclidean steps for pairs < {N}: {best[0]} at {best[1]}  ; consecutive Fibonacci below {N}: {fib[-3:]}")
print(f"  the worst-case pair {best[1]} is consecutive Fibonacci? {best[1][0] in fib and best[1][1] in fib}")
print("  => Lame's theorem: the Euclidean/continued-fraction descent is longest for FIBONACCI (golden CF).")
print("  This finite check illustrates ordinary Euclid only; no JC(2) decreasing invariant is proved here.")

sep("SUMMARY -- exact survivors and honest frontier")
print("""  SURVIVES: exact identities for P=(x+iy)^d and a finite Euclidean/Fibonacci check.
  DOES NOT SURVIVE: a binary-to-general symmetric reduction, a GMC/LRC/NC2 transfer, or a proved
  Newton/continued-fraction descent for Keller maps.
  LIVE SEPARATE PROGRAMS: four-variable symmetric/VC; planar leading-form/Jelonek geometry;
  Newton-polygon descent with an as-yet-unproved decreasing invariant.
  JC(2) remains open.  See MISTAKE-237 and HYP-8905.""")
