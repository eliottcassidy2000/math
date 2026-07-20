import sympy as sp
from math import factorial
z,zb=sp.symbols('z zb')
def E1(e):
    e=sp.expand(e)
    if e==0: return sp.Integer(0)
    p=sp.Poly(e,z,zb); t=0
    for (a,b),c in zip(p.monoms(),p.coeffs()):
        if a==b: t+=c*factorial(a)
    return sp.expand(t)
def seq(P,M=7): return [E1(sp.expand(P**m)) for m in range(1,M+1)]
print("A MIXED-CHARGE P must fall in one of two cases. Test both structurally.\n")
print("CASE B -- charges of BOTH signs, none zero.  Then for c>0 and d<0 present,")
print("  taking |d| factors of charge c and c factors of charge d gives total charge 0,")
print("  so P^(c+|d|) HAS balanced monomials.  Do they cancel?")
for P in [z+zb, z**2+zb, z+zb**2, z**2+zb**2, z-zb, z**3+zb, 2*z+zb, z**2-zb**2,
          z+zb+z**2, z**3+zb**2, z**2*zb**0+zb**3]:
    s=seq(P,6)
    first=next((i+1 for i,v in enumerate(s) if v!=0), None)
    print(f"   P={str(P):18s} E[P^m] m=1..6 = {s}   first nonzero at m={first}")
print()
print("CASE A -- P contains a CHARGE-0 part P0 (necessarily with E[P0]=0).")
print("  The natural candidates use P0 = z*zb - 1, or (z*zb)^2 - 2*z*zb, etc.")
for P0 in [z*zb-1, (z*zb)**2-2, (z*zb)**2-4*z*zb+2, z*zb-1+ (z*zb)**2-2]:
    for extra in [0, z, zb, z**2, z-zb, z*zb*0+z**3]:
        P=sp.expand(P0+extra)
        s=seq(P,6)
        first=next((i+1 for i,v in enumerate(s) if v!=0), None)
        tag = "  <-- NULLCONE!" if first is None else ""
        print(f"   P0={str(P0):22s} +{str(extra):6s}: E[P^m]={s}  first nz m={first}{tag}")
print()
print("CONCLUSION CHECK: any P above with ALL E[P^m]=0 and MIXED charge would refute")
print("the characterisation.  Scan the 'first nz' column for None on a mixed-charge row.")
