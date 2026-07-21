import sympy as sp
from math import factorial
Z,Zb=sp.symbols('Z Zb')
al,be,ga,de,w=sp.symbols('alpha beta gamma delta w')
def E(expr):
    e=sp.expand(expr)
    if e==0: return sp.Integer(0)
    p=sp.Poly(e,Z,Zb); t=0
    for (aa,bb),c in zip(p.monoms(),p.coeffs()):
        if aa==bb: t+=c*factorial(aa)
    return sp.expand(t)
print("TWO-STRADDLE witness (the HYP-8470 residual klein left open): charges {-3,-1,+1,+3}")
print("P = alpha Z^3 + gamma Z + beta Zb + delta Zb^3   (h=1 pair (be,ga), h=3 pair (al,de))")
print("Both charge-0 pairings be*ga (=-1+1) and al*de (=+3-3) present -> can they cancel?")
print()
P=al*Z**3+ga*Z+be*Zb+de*Zb**3
mom=[]
for m in range(1,9):
    e=sp.expand(E(sp.expand(P**m)))
    if e!=0: mom.append((m,e))
print("nonzero moments E[P^m]:")
for m,e in mom[:6]:
    print(f"   m={m}: {sp.factor(e)}")
print()
print("BOTTOM-UP RESULTANT TOWER: the low moment fires the bottom shell first.")
print("   E[P^2] = 2(be*ga + ...)?  Check which pairs appear at m=2:")
m2=dict(mom)[2] if 2 in dict(mom) else 0
print(f"   E[P^2] = {sp.expand(m2)}   (the h=1 AND h=3 straddles both at m=2!)")
print()
print("ONE-SIDED locus (to exclude): P single-signed <=> {be=de=0} (all charges>0) or")
print("{al=ga=0} (all <0). One-sided ideal J = <(be or de)*(al or ga)> ~ products across signs.")
# nullcone: all moments vanish; test V(moments) minus one-sided is empty
I=[e for _,e in mom]
# saturate by a genuine-two-sided witness: need a positive charge (al or ga) AND negative (be or de)
print("TEST: is V(moment ideal) cap (two-sided) EMPTY?  Saturate by (ga*be), (ga*de), (al*be), (al*de):")
for sat,lab in [(ga*be,'ga*be'),(al*de,'al*de'),(ga*de,'ga*de'),(al*be,'al*be')]:
    G=sp.groebner(I+[1-w*sat],al,be,ga,de,w,order='grevlex')
    print(f"   saturate by {lab}: V empty (1 in ideal)? {list(G)==[sp.Integer(1)]}")
