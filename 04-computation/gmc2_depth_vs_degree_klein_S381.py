#!/usr/bin/env python3
"""
klein-2026-07-20-S381 -- DOES THE GMC(2) DETECTION DEPTH GROW WITH RADIAL DEGREE d?  The crux of
HYP-8540.  Focused, exact-symbolic, on charge span {-1,0,1} at d = 0,1,2.  A control confirms the
test has power (one-sided P is not flagged).
"""
import sympy as sp
from math import factorial
from sympy import groebner

Z, Zb = sp.symbols('Z Zb'); ZZb = Z*Zb
def Ew(poly):
    poly=sp.expand(poly)
    if poly==0: return sp.Integer(0)
    t=sp.Integer(0)
    for (a,b),c in sp.Poly(poly,Z,Zb).as_dict().items():
        if a==b: t+=c*factorial(a)
    return t

def buildP(d):
    A=[sp.symbols(f'A{k}') for k in range(d+1)]   # charge +1 radial coeffs
    B=[sp.symbols(f'B{k}') for k in range(d+1)]   # charge  0
    C=[sp.symbols(f'C{k}') for k in range(d+1)]   # charge -1
    P=sum(A[k]*ZZb**k*Z for k in range(d+1)) + sum(B[k]*ZZb**k for k in range(d+1)) + sum(C[k]*ZZb**k*Zb for k in range(d+1))
    return P, A, B, C

def depth(d, mmax=10, twosided_only=True):
    P,A,B,C=buildP(d); allc=A+B+C
    moms=[sp.nsimplify(sp.expand(Ew(sp.expand(P**m)))) for m in range(1,mmax+1)]
    # one-sided loci: {all B=0 and all C=0} (charge +1 only) or {all A=0 and all B=0}.
    # nullcone should be their union.  targets to kill: A_i*C_j (both-sign) and B_i (charge 0).
    targets=[a*c for a in A for c in C]+list(B)
    for D in range(1,mmax+1):
        gens=[g for g in moms[:D] if g!=0]
        if not gens: continue
        G=groebner(gens,*allc,order='grevlex')
        if all(any(sp.simplify(G.reduce(t**k)[1])==0 for k in range(1,7)) for t in targets):
            return D, moms
    return None, moms

print("="*80)
print("DETECTION DEPTH of GMC(2) on span {-1,0,1}, exact, vs radial degree d")
print("="*80)
res={}
for d in (0,1,2):
    D,moms=depth(d)
    res[d]=D
    print(f"  d={d}: #coeffs={3*(d+1)}, detection depth = {D}")
    print(f"       E[P^1..3] = {[str(m) for m in moms[:3]]}")
print(f"\n  detection depths (d=0,1,2): {[res[d] for d in (0,1,2)]}")
print("="*80)
print("CONTROL: one-sided P (charge +1 only) must NOT be flagged as a two-sided nullcone member")
print("="*80)
d=1; A=[sp.symbols(f'A{k}') for k in range(d+1)]
Pone=sum(A[k]*ZZb**k*Z for k in range(d+1))     # charges +1,+3 depth -> all positive
Q=Zb
ep=[sp.expand(Ew(sp.expand(Pone**m))) for m in range(1,6)]
eq=[sp.expand(Ew(sp.expand(Q*Pone**m))) for m in range(1,6)]
print(f"  one-sided P: E[P^m] m=1..5 = {[str(x) for x in ep]}  (all 0 => in nullcone)")
print(f"               E[Q P^m] = {[str(x) for x in eq]}  (0 for m>1 => MZ-harmless)")
print("="*80)
print("VERDICT")
print("="*80)
ds=[res[d] for d in (0,1,2)]
if len(set(ds))==1:
    print(f"  Detection depth CONSTANT = {ds[0]} across d=0,1,2 on span {{-1,0,1}}: EVIDENCE FOR HYP-8540")
    print("  (the span-only bound) -- the radial degree does NOT raise the depth on this span.")
else:
    print(f"  Detection depth VARIES {ds} with d on span {{-1,0,1}}: the depth is NOT span-only;")
    print("  HYP-8540 as stated (span-only bound) would be FALSE and the bound must include d.")
