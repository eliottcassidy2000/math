import sympy as sp
from math import factorial
z1,z1b,z2,z2b = sp.symbols('z1 z1b z2 z2b')
def E(expr):
    e=sp.expand(expr)
    if e==0: return sp.Integer(0)
    p=sp.Poly(e,z1,z1b,z2,z2b); tot=0
    for (a1,b1,a2,b2),c in zip(p.monoms(),p.coeffs()):
        if a1==b1 and a2==b2: tot+=c*factorial(a1)*factorial(a2)
    return sp.expand(tot)
print("(1) THE 4-TERM COUNTEREXAMPLE  P' = (1+Z2)(W2 - Z1 W1)  = (1+z2)(z2b - z1*z1b)")
Pp = sp.expand((1+z2)*(z2b - z1*z1b))
print("    P' =", Pp, "   terms:", len(sp.Poly(Pp,z1,z1b,z2,z2b).monoms()),
      "  degree:", sp.Poly(Pp,z1,z1b,z2,z2b).total_degree())
ok=True; seq=[]
for m in range(1,13):
    Pm=sp.expand(Pp**m); e1=E(Pm); e2=E(sp.expand(z2*Pm)); seq.append(e2)
    if e1!=0: ok=False
    print(f"    m={m:2d}  E[P'^m]={e1}   E[Z2 P'^m]={e2}    m!={factorial(m)}   "
          f"{'= m!' if e2==factorial(m) else ''}")
print(f"    all E[P'^m]=0 through m=12: {ok}")
print(f"    => {'4-TERM COUNTEREXAMPLE CONFIRMED (6 terms -> 4)' if ok and all(s!=0 for s in seq) else 'FAILS'}")
print()
print("(2) THE n=2 'CANDIDATES' -- GMC(2) IS A THEOREM, so they MUST die at higher m.")
print("    My search only tested m<=4 and flagged ANY nonzero -- too weak, since GMC")
print("    only requires E[QP^m]=0 for m >> 0.  Re-testing properly:")
def E1(expr):
    e=sp.expand(expr)
    if e==0: return sp.Integer(0)
    p=sp.Poly(e,z1,z1b); tot=0
    for (a,b),c in zip(p.monoms(),p.coeffs()):
        if a==b: tot+=c*factorial(a)
    return sp.expand(tot)
cands=[(-z1**3 - z1**2*z1b - z1**2 - z1, z1b),
       (-z1**3 - z1**2 - z1, z1b),
       (-z1**3 + z1**2*z1b - z1**2 - z1, z1b)]
for P,Q in cands:
    ms=[E1(sp.expand(P**m)) for m in range(1,9)]
    qs=[E1(sp.expand(Q*P**m)) for m in range(1,9)]
    print(f"    P={P}")
    print(f"       E[P^m] m=1..8 : {ms}")
    print(f"       E[QP^m] m=1..8: {qs}")
    print(f"       -> {'NOT a counterexample (E[P^m] nonzero somewhere)' if any(x!=0 for x in ms) else ('E[QP^m] dies at large m -> GMC(2) intact' if qs[-1]==0 else 'STILL NONZERO -- inspect!')}")
