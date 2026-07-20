import sympy as sp
from math import factorial
z1,z1b,z2,z2b=sp.symbols('z1 z1b z2 z2b')
def E(e):
    e=sp.expand(e)
    if e==0: return sp.Integer(0)
    p=sp.Poly(e,z1,z1b,z2,z2b); t=0
    for (a1,b1,a2,b2),c in zip(p.monoms(),p.coeffs()):
        if a1==b1 and a2==b2: t+=c*factorial(a1)*factorial(a2)
    return sp.expand(t)
def parts(e):
    """split by VECTOR charge (a1-b1, a2-b2)"""
    e=sp.expand(e); p=sp.Poly(e,z1,z1b,z2,z2b); d={}
    for (a1,b1,a2,b2),c in zip(p.monoms(),p.coeffs()):
        q=(a1-b1,a2-b2)
        d[q]=d.get(q,0)+c*z1**a1*z1b**b1*z2**a2*z2b**b2
    return d
P=sp.expand((1+z2)*(z2b - z1*z1b))
d=parts(P)
print("WHY n=4 ESCAPES THE n=2 ARGUMENT.")
print("The charge group is Z^(n/2): RANK 1 at n=2, RANK 2 at n=4.\n")
print("Vector-charge decomposition of the 4-term witness P' = (1+Z2)(W2 - Z1W1):")
for q in sorted(d): print(f"   charge {q}:  {sp.expand(d[q])}")
P0=d.get((0,0),0)
print(f"\n   charge-0 part P0 = {sp.expand(P0)}   E[P0] = {E(P0)}  (vanishes: |Z2|^2 and |Z1|^2 cancel)")
print(f"   E[P0^2] = {E(sp.expand(P0**2))}   <- STRICTLY POSITIVE, as at n=2")
cross = 2*E(sp.expand(d.get((0,1),0)*d.get((0,-1),0)))
print(f"   2*E[P_(0,1) * P_(0,-1)] = {cross}   <- the OPPOSITE-CHARGE cross term")
print(f"   sum = {E(sp.expand(P0**2)) + cross}   and directly E[P^2] = {E(sp.expand(P**2))}")
print()
print("SO: at n=2 with SIGN-COHERENT charges there is NO opposite charge available, so")
print("    E[P^2] = E[P0^2] = c^T H c > 0 by Hankel positive-definiteness -- P0 must vanish.")
print("    At n=4 the charges (0,+1) and (0,-1) BOTH occur, and their cross term is")
print("    exactly -E[P0^2].  The extra rank of the charge lattice buys the cancellation.")
print("    THAT is the dimensional mechanism separating GMC(2) from GMC(4).")
