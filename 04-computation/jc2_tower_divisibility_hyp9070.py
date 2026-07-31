"""DECISIVE TEST 2 of HYP-9070: does (L1') really give H^{a-b} | Jac(P_{n-1},H)?
Verify on genuine Jacobian pairs with a-b >= 1, and probe the next order."""
import sympy as sp
x, y = sp.symbols('x y')

def jac(P, Q): return sp.expand(sp.diff(P,x)*sp.diff(Q,y) - sp.diff(P,y)*sp.diff(Q,x))
def hom(F, d):
    F = sp.Poly(sp.expand(F), x, y)
    return sp.expand(sum(c*x**a*y**b for (a,b),c in F.terms() if a+b==d))
def deg(F): return sp.Poly(sp.expand(F), x, y).total_degree()

def analyse(P, Q, name):
    n, m = deg(P), deg(Q)
    g = sp.gcd(n, m); a, b = n//g, m//g
    Pn, Qm = hom(P,n), hom(Q,m)
    # extract H: P_n = c H^a  => H = (P_n)^(1/a) up to scalar; get it by factoring
    fl = sp.factor_list(Pn)[1]
    H = sp.prod([f for f,e in fl if sp.total_degree(f)>0])   # radical (K=1 => single factor)
    print(f"  {name}: deg=({n},{m}) g={g} (a,b)=({a},{b}) H={H} Jac={sp.simplify(jac(P,Q))}")
    # verify L1
    L1 = sp.expand(jac(Pn, hom(Q,m-1)) + jac(hom(P,n-1), Qm))
    print(f"     L1 residual: {sp.simplify(L1)}")
    # the divisibility claim, with a>=b or b>=a handled
    if a >= b:
        target = sp.expand(jac(hom(P,n-1), H)); e = a-b
    else:
        target = sp.expand(jac(H, hom(Q,m-1))); e = b-a
    if e == 0:
        print("     a=b: no divisibility claim"); return
    Hp = sp.expand(H**e)
    q_, r_ = sp.div(sp.Poly(target, x, y), sp.Poly(Hp, x, y)) if target != 0 else (0,0)
    ok = (target == 0) or (sp.simplify(sp.expand(r_.as_expr() if hasattr(r_,'as_expr') else r_)) == 0)
    print(f"     H^{e} | Jac(...) ?  target={'0 (vacuous)' if target==0 else 'nonzero'}  -> {'HOLDS' if ok else '*** FAILS ***'}")

# automorphisms with a-b >= 1
u = sp.expand(x + y**2)
analyse(u, sp.expand(y + u**3), "F=(x+y^2, y+(x+y^2)^3)")           # (a,b)=(1,3)
P2 = sp.expand(x + y**2); Q2 = sp.expand(y + P2**3)
analyse(sp.expand(P2 + Q2**2), Q2, "deeper composite")               # (a,b)=(2,1)
v = sp.expand(y + x**2)
analyse(sp.expand(x + v**4), v, "F=(x+(y+x^2)^4, y+x^2)")            # (a,b)=(4,1)
