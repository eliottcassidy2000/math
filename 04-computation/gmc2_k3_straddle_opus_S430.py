import sympy as sp
from math import factorial
Z,Zb=sp.symbols('Z Zb')
a1,a3,a5,b1,b3,b5=sp.symbols('a1 a3 a5 b1 b3 b5')
def E(expr):
    e=sp.expand(expr)
    if e==0: return sp.Integer(0)
    p=sp.Poly(e,Z,Zb); t=0
    for (aa,bb),c in zip(p.monoms(),p.coeffs()):
        if aa==bb: t+=c*factorial(aa)
    return sp.expand(t)
print("k=3 STRADDLE (induction step for HYP-8550): charges {+-1,+-3,+-5}")
print("P = a5 Z^5 + a3 Z^3 + a1 Z + b1 Zb + b3 Zb^3 + b5 Zb^5")
print("CLAIM: E[P^{2k}] mod <lower even moments> is a MONOMIAL in the bottom pair (a1,b1).")
P=a5*Z**5+a3*Z**3+a1*Z+b1*Zb+b3*Zb**3+b5*Zb**5
mom=[]
for m in [2,4,6]:
    e=sp.expand(E(sp.expand(P**m)))
    if e!=0: mom.append((m,e))
    print(f"  E[P^{m}] computed ({len(sp.Poly(e,a1,a3,a5,b1,b3,b5).monoms()) if e!=0 else 0} terms)")
print(f"  E[P^2] = {sp.factor(mom[0][1])}")
print()
print("  Bottom-up reduction (this is the tower):")
gens=[a1,a3,a5,b1,b3,b5]
G1=sp.groebner([mom[0][1]],*gens,order='lex')
r4=sp.reduced(mom[1][1],G1)[1]
print(f"  E[P^4] mod <E[P^2]> = {sp.factor(r4)}")
G2=sp.groebner([mom[0][1],mom[1][1]],*gens,order='lex')
r6=sp.reduced(mom[2][1],G2)[1]
print(f"  E[P^6] mod <E[P^2],E[P^4]> = {sp.factor(r6)}")
print()
print("  (k=3 needs moments up to m=? -- checking whether E[P^6] already exposes a bottom monomial")
print("   or whether higher m is needed; the pattern from k=2 was terminate at 2k=6 -> here maybe 2k=? )")
