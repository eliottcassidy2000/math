#!/usr/bin/env python3
"""
death-star-2026-07-20-S59z (HYP-8245) -- the Zhao Vanishing Conjecture VERIFICATION
machinery, to demonstrate the witness certificate is FINITE and COMPUTABLE (Lean-able).
VC(4): for homogeneous quartic P, if Delta^m(P^m)=0 for all m>=1 then Delta^m(P^{m+1})=0
for all m>>0. A VC COUNTEREXAMPLE = Hessian-nilpotent homogeneous quartic P with
Delta^m(P^m)=0 for all m (automatic for Hessian-nilpotent, Zhao) but Delta^m(P^{m+1})!=0.
Here: implement Delta, Delta^m(P^m); verify on the isotropic harmonic example (VC holds);
show the machinery + note the true witness needs F's reduction (dimension-dependent).
"""
from fractions import Fraction as Fr
from itertools import product

# multivariate polys over Q(i): coeff = (re, im) as Fr pairs; monomial = tuple of exps
def cadd(a,b): return (a[0]+b[0], a[1]+b[1])
def cmul(a,b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def pmul(A,B):
    R={}
    for ka,ca in A.items():
        for kb,cb in B.items():
            k=tuple(x+y for x,y in zip(ka,kb)); v=cadd(R.get(k,(Fr(0),Fr(0))),cmul(ca,cb))
            if v!=(Fr(0),Fr(0)): R[k]=v
            elif k in R: del R[k]
    return R
def padd(*Ps):
    R={}
    for A in Ps:
        for k,c in A.items():
            v=cadd(R.get(k,(Fr(0),Fr(0))),c)
            if v!=(Fr(0),Fr(0)): R[k]=v
            elif k in R: del R[k]
    return R
def d2(A,i):  # second partial d^2/dx_i^2
    R={}
    for k,c in A.items():
        if k[i]>=2:
            k2=list(k); e=k[i]; k2[i]-=2
            R[tuple(k2)]=cadd(R.get(tuple(k2),(Fr(0),Fr(0))), cmul(c,(Fr(e*(e-1)),Fr(0))))
    return {k:v for k,v in R.items() if v!=(Fr(0),Fr(0))}
def laplacian(A,n): return padd(*[d2(A,i) for i in range(n)])
def ppow(A,m):
    R={tuple([0]*NV):(Fr(1),Fr(0))}
    for _ in range(m): R=pmul(R,A)
    return R
def is_zero(A): return len(A)==0

NV=2
def var(i):
    k=[0]*NV; k[i]=1; return {tuple(k):(Fr(1),Fr(0))}
x0,x1=var(0),var(1)
i_unit=(Fr(0),Fr(1))
L=padd(x0, {tuple([0,1]):i_unit})   # L = x0 + i x1 (isotropic: sum of squares of coeffs = 1+i^2=0)
P=ppow(L,4)                          # P = (x0+i x1)^4, homogeneous quartic
print("=== isotropic harmonic example P=(x0+i x1)^4 (Hessian nilpotent) ===")
# Hessian nilpotent check: H = 12 L^2 a a^T with a=(1,i), a.a=0 => H^2=0 (structural)
print("  Hessian nilpotent: H = 12 L^2 (1,i)(1,i)^T, (1,i).(1,i)=1+i^2=0 => H^2=0. YES (structural)")
print("  VC data Delta^m(P^m) and Delta^m(P^{m+1}):")
for m in range(1,5):
    Pm=ppow(P,m); Pm1=ppow(P,m+1)
    dm=Pm
    for _ in range(m): dm=laplacian(dm,NV)
    dm1=Pm1
    for _ in range(m): dm1=laplacian(dm1,NV)
    print(f"    m={m}: Delta^{m}(P^{m})=0? {is_zero(dm)};  Delta^{m}(P^{{m+1}})=0? {is_zero(dm1)}")
print("  => both always 0 (powers of an isotropic linear form are HARMONIC) => VC HOLDS here.")
print("     x+grad P is a SHEAR (L preserved) = INVERTIBLE: consistent, NOT a JC counterexample.")

print("\n=== THE POINT (verification is finite & computable = Lean-able) ===")
print("  A VC COUNTEREXAMPLE needs a Hessian-nilpotent homogeneous quartic P with")
print("  Delta^m(P^m)=0 for ALL m (automatic, Zhao) BUT Delta^m(P^{m+1}) != 0 for some m.")
print("  Such P comes from a NON-invertible Keller map = the JC counterexample F via")
print("  Yagzhev (cubic-homog) + de Bondt-van den Essen (symmetric quartic) reduction.")
print("  Once P is produced, the certificate is a FINITE set of Delta^m(P^m) computations")
print("  (this machinery) -> directly Lean-formalizable as polynomial-identity checks,")
print("  INDEPENDENT of the VC<=>JC equivalence (the robustness point of the assessment).")
print("  Feasibility gate = the dimension the F-reduction lands in (pending recipe agent).")
