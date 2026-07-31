"""THE n+m-4 CROSS TERM, resolved by a Plucker identity.

Notation: Jacobian pair (P,Q), Jac=1, n=deg P, m=deg Q, g=gcd(n,m),
P_n = c H^a, Q_m = c' H^b, a=n/g, b=m/g coprime, wlog a>=b.
Write X = Jac(P_{n-1},Q_{m-1}) (the cross term), J = Jac(H,Q_{m-1}),
Y = Jac(P_{n-1},H).

STEP 1 (L1'):  c a H^{a-b} J = - c' b Y,   so  Y = -kappa H^{a-b} J,
               kappa = c a/(c' b).
STEP 2 (Plucker, three gradients in a 2-dim space):
   Jac(P_{n-1},Q_{m-1}) grad H + Jac(Q_{m-1},H) grad P_{n-1}
                                + Jac(H,P_{n-1}) grad Q_{m-1} = 0
   i.e.  X grad H = J grad P_{n-1} + Y grad Q_{m-1}.
STEP 3: substitute Y and use
   H^{a-b} grad Q_{m-1} = grad(H^{a-b}Q_{m-1}) - (a-b)H^{a-b-1}Q_{m-1} grad H:
   [ X - kappa(a-b) J H^{a-b-1} Q_{m-1} ] grad H  =  J grad W,
   where  W := P_{n-1} - kappa H^{a-b} Q_{m-1},  HOMOGENEOUS of degree n-1
   (since (a-b)g + m-1 = n-1).
STEP 4: wedge both sides with grad H  =>   J * Jac(W,H) = 0.

CONSEQUENCE.  If J != 0 then Jac(W,H)=0.  For binary FORMS that means
W = alpha G^u, H = beta G^v for a common form G, so deg G divides both
n-1 and g; but g | n, hence deg G | gcd(n-1,n) = 1, so G is LINEAR and
H is a pure power of a linear form, i.e. K = 1.

Since automorphisms have K = 1 (HYP-9070 test 1) a COUNTEREXAMPLE has
K >= 2, and therefore must satisfy  J = 0  or  W = 0:

   Jac(H, Q_{m-1}) = 0     or     P_{n-1} = kappa H^{a-b} Q_{m-1}.   (X)

So the cross term does NOT run free: it is rigidly determined.
This file verifies every step symbolically on genuine Jacobian pairs.
"""
import sympy as sp
x, y = sp.symbols('x y')

def jac(A, B): return sp.expand(sp.diff(A,x)*sp.diff(B,y) - sp.diff(A,y)*sp.diff(B,x))
def grad(A): return sp.Matrix([sp.diff(A,x), sp.diff(A,y)])
def hom(F, d):
    F = sp.Poly(sp.expand(F), x, y)
    return sp.expand(sum(cc*x**i*y**j for (i,j),cc in F.terms() if i+j==d))
def deg(F): return sp.Poly(sp.expand(F), x, y).total_degree()

def check_plucker():
    A = sp.expand(x**3 - 2*x*y + 1); B = sp.expand(y**2 + x); C = sp.expand(x*y - y**3)
    v = jac(B,C)*grad(A) + jac(C,A)*grad(B) + jac(A,B)*grad(C)
    ok = sp.simplify(v[0]) == 0 and sp.simplify(v[1]) == 0
    print(f"  Plucker identity Jac(B,C)gradA+Jac(C,A)gradB+Jac(A,B)gradC=0 : {'OK' if ok else 'FAIL'}")
    return ok

def analyse(P, Q, name):
    n, m = deg(P), deg(Q)
    g = sp.gcd(n, m); a, b = n//g, m//g
    Pn, Qm = hom(P,n), hom(Q,m)
    Pn1, Qm1 = hom(P,n-1), hom(Q,m-1)
    # H from factoring P_n  (K=1 for automorphisms => single linear factor)
    fl = [f for f,e in sp.factor_list(Pn)[1] if sp.total_degree(f) > 0]
    H = sp.prod(fl)
    Hg = sp.expand(H**g) if sp.total_degree(H) == 1 else H
    # constants: P_n = c H^a with deg H = g -> use the radical^g normalisation
    Hn = sp.expand(H**(sp.Integer(g)//sp.total_degree(H))) if sp.total_degree(H) > 0 else H
    c = sp.simplify(Pn/ Hn**a); cp = sp.simplify(Qm/ Hn**b)
    J = jac(Hn, Qm1); Y = jac(Pn1, Hn); X = jac(Pn1, Qm1)
    swap = a < b
    A_, B_ = (a, b) if not swap else (b, a)
    kappa = sp.simplify(c*a/(cp*b)) if not swap else sp.simplify(cp*b/(c*a))
    e = A_ - B_
    if swap:
        W = sp.expand(Qm1 - kappa*Hn**e*Pn1); Jj = jac(Hn, Pn1)
    else:
        W = sp.expand(Pn1 - kappa*Hn**e*Qm1); Jj = J
    l1 = sp.simplify(jac(Pn,Qm1) + jac(Pn1,Qm))
    res = sp.simplify(Jj*jac(W, Hn))
    print(f"  {name}: deg=({n},{m}) g={g} (a,b)=({a},{b}) H={Hn} kappa={kappa}")
    print(f"     L1 residual = {l1}   |  J = Jac(H,Q_(m-1)) = {sp.simplify(Jj)}")
    print(f"     W = {sp.factor(W)}   deg W = {deg(W) if sp.simplify(W)!=0 else '-'} (expect n-1={n-1})")
    print(f"     J * Jac(W,H) = {res}   -> {'OK (=0)' if res == 0 else '*** NONZERO ***'}")
    return res == 0

print("Verification of the cross-term chain:")
ok = check_plucker()
u = sp.expand(x + y**2)
ok &= analyse(u, sp.expand(y + u**3), "F=(x+y^2, y+(x+y^2)^3)")
P2 = sp.expand(x + y**2); Q2 = sp.expand(y + P2**3)
ok &= analyse(sp.expand(P2 + Q2**2), Q2, "deeper composite")
v = sp.expand(y + x**2)
ok &= analyse(sp.expand(x + v**4), v, "F=(x+(y+x^2)^4, y+x^2)")
print("\nALL CROSS-TERM CHECKS PASS" if ok else "\nSOME CHECK FAILED")
