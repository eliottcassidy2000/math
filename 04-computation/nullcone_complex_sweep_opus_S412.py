import sympy as sp
from math import factorial
from itertools import product
z,zb=sp.symbols('z zb')
I=sp.I
def E1(e):
    e=sp.expand(e)
    if e==0: return sp.Integer(0)
    p=sp.Poly(e,z,zb); t=0
    for (a,b),c in zip(p.monoms(),p.coeffs()):
        if a==b: t+=c*factorial(a)
    return sp.expand(t)
def chg(e):
    e=sp.expand(e)
    if e==0: return set()
    p=sp.Poly(e,z,zb)
    return {a-b for (a,b),c in zip(p.monoms(),p.coeffs()) if c!=0}
def definite(ch): return bool(ch) and (all(c>=1 for c in ch) or all(c<=-1 for c in ch))
print("CRITICAL RECHECK: my S411 sweep used REAL coefficients {-1,0,1}.")
print("GMC is over C.  Redo with GAUSSIAN-INTEGER coefficients {0,+-1,+-i}.\n")
MON=[(a,b) for a in range(3) for b in range(3) if 1<=a+b<=2]   # deg 1..2, 5 monomials
COEF=[0,1,-1,I,-I]
null=[]; mixed=[]; checked=0
for co in product(COEF,repeat=len(MON)):
    if not any(c!=0 for c in co): continue
    checked+=1
    P=sum(c*z**a*zb**b for c,(a,b) in zip(co,MON) if c!=0)
    if all(E1(sp.expand(P**m))==0 for m in range(1,8)):
        null.append(P)
        if not definite(chg(P)): mixed.append((P,chg(P)))
print(f"  deg<=2, coeffs in {{0,+-1,+-i}}: scanned {checked}")
print(f"  nullcone members            : {len(null)}")
print(f"  NOT charge-definite         : {len(mixed)}")
for P,ch in mixed[:10]: print(f"     *** P={P}   charges {sorted(ch)}")
if not mixed: print("  -> every complex nullcone member is still charge-definite at deg<=2")
print()
# degree 3, restricted support, complex
MON3=[(1,0),(0,1),(2,0),(1,1),(0,2),(2,1),(1,2),(3,0),(0,3)]
null3=[]; mixed3=[]; checked3=0
for co in product([0,1,-1,I],repeat=len(MON3)):
    if not any(c!=0 for c in co): continue
    checked3+=1
    P=sum(c*z**a*zb**b for c,(a,b) in zip(co,MON3) if c!=0)
    ch=chg(P)
    if definite(ch): continue        # only need to test the MIXED ones
    if E1(sp.expand(P))!=0: continue
    if E1(sp.expand(P**2))!=0: continue
    if all(E1(sp.expand(P**m))==0 for m in range(3,7)):
        mixed3.append((P,ch))
print(f"  deg<=3, coeffs in {{0,1,-1,i}}: MIXED-charge P scanned {checked3}")
print(f"  mixed-charge NULLCONE members found: {len(mixed3)}")
for P,ch in mixed3[:10]: print(f"     *** P={P}   charges {sorted(ch)}")
if not mixed3:
    print("  -> NO mixed-charge complex nullcone member.  Conjecture survives over C.")
