import sympy as sp
from math import factorial
z,zb,t=sp.symbols('z zb t')
def E1(e):
    e=sp.expand(e)
    if e==0: return sp.Integer(0)
    p=sp.Poly(e,z,zb); tot=0
    for (a,b),c in zip(p.monoms(),p.coeffs()):
        if a==b: tot+=c*factorial(a)
    return sp.expand(tot)
def run(D,M,span1=True):
    A=sp.symbols(f'a0:{D+1}'); B=sp.symbols(f'b0:{D+1}'); C=sp.symbols(f'c0:{D+1}')
    s=z*zb
    P=sp.expand(z*sum(A[k]*s**k for k in range(D+1))
                + sum(B[k]*s**k for k in range(D+1))
                + zb*sum(C[k]*s**k for k in range(D+1)))
    eqs=[]; Pm=sp.Integer(1)
    for m in range(1,M+1):
        Pm=sp.expand(Pm*P); e=E1(Pm)
        if e!=0: eqs.append(sp.expand(e))
    gens=list(A)+list(B)+list(C)
    print(f"  D={D}: {len(gens)} unknowns, {len(eqs)} equations")
    ok=True
    for i in range(D+1):
        for j in range(D+1):
            G=sp.groebner(eqs+[1-t*A[i]*C[j]],*(gens+[t]),order='grevlex')
            triv=(list(G.exprs)==[sp.Integer(1)])
            print(f"     a{i}*c{j} saturation trivial: {triv}")
            ok = ok and triv
    print(f"  => {'ONE-SIDED PROVED for deg <= '+str(D) if ok else 'NOT all trivial -- inspect'}")
print("PUSH: charge span {-1,0,+1}, degree <= 2 (9 unknowns)")
run(2,8)
