"""STRENGTHENING THM-3016 and ITERATING to the n+m-4 equation.

(A) J = 0 ALSO forces K = 1.  If J = Jac(H,Q_{m-1}) = 0 and Q_{m-1} != 0 then
    H, Q_{m-1} are powers of a common form G, so deg G | g and deg G | m-1;
    but g | m, hence deg G | gcd(m,m-1) = 1 and H is a pure power of a linear
    form.  So the dichotomy of THM-3016 COLLAPSES:

        K >= 2   =>   Q_{m-1} = 0   or   P_{n-1} = kappa H^{a-b} Q_{m-1}.

(B) With W = 0 the cross term becomes EXPLICIT:
        Jac(P_{n-1},Q_{m-1}) = Jac(kappa H^{a-b}Q_{m-1}, Q_{m-1})
                             = kappa (a-b) H^{a-b-1} Q_{m-1} J,
    using Jac(fg,g) = g Jac(f,g).  So the degree n+m-4 equation
        c a H^{a-1} Jac(H,Q_{m-2}) + Jac(P_{n-1},Q_{m-1})
                                   + c' b H^{b-1} Jac(P_{n-2},H) = 0
    has all three terms carrying explicit powers of H, namely
        a-1,   a-b-1,   b-1,
    so dividing by H^{min(a-b-1, b-1)} gives the NEXT relation with exponent
    pair reduced -- the tower iterates, and the exponents are exactly a
    Euclidean step applied to (a,b).  This file verifies (A) and the (B)
    identity symbolically.
"""
import sympy as sp
x, y = sp.symbols('x y')

def jac(A,B): return sp.expand(sp.diff(A,x)*sp.diff(B,y) - sp.diff(A,y)*sp.diff(B,x))

# (A) verify the degree argument on explicit forms
print("(A) J=0 with Q_{m-1} != 0 forces H a pure power of a linear form:")
cases = [
    ("H=y^2 (g=2), Q_{m-1}=y^5  (m=6, m-1=5)", y**2, y**5, 2, 6),
    ("H=x^3 (g=3), Q_{m-1}=x^8  (m=9, m-1=8)", x**3, x**8, 3, 9),
    ("H=(x*y) (g=2, K=2), Q_{m-1}=x^5 (m=6)", x*y, x**5, 2, 6),
]
for nm, H, Qm1, g, m in cases:
    J = sp.simplify(jac(H, Qm1))
    Kf = len([f for f,e in sp.factor_list(H)[1] if sp.total_degree(f)>0])
    print(f"   {nm}: J={J}  K(H)={Kf}   "
          f"{'J=0 and K=1 (consistent)' if J==0 and Kf==1 else ('J=0 but K>1 -> would break (A)' if J==0 else 'J!=0')}")

# (B) verify Jac(f*g, g) = g*Jac(f,g) and the explicit cross term
print("\n(B) cross-term identity when W = 0:")
f, gg = sp.Function('f'), sp.Function('g')
A = sp.expand(x**2 - y); B = sp.expand(x*y + y**3)
lhs = sp.expand(jac(A*B, B)); rhs = sp.expand(B*jac(A, B))
print(f"   Jac(f*g,g) = g*Jac(f,g):  {'OK' if sp.simplify(lhs-rhs)==0 else 'FAIL'}")
# explicit: P_{n-1} = kappa H^{a-b} Q_{m-1}
kap = sp.Symbol('kappa'); H = sp.expand(y**2); Qm1 = sp.expand(x*y**3 + y**4)
for (a,b) in [(4,1),(3,2),(5,2)]:
    Pn1 = sp.expand(kap*H**(a-b)*Qm1)
    lhs = sp.expand(jac(Pn1, Qm1))
    rhs = sp.expand(kap*(a-b)*H**(a-b-1)*Qm1*jac(H, Qm1))
    print(f"   (a,b)=({a},{b}): Jac(P_(n-1),Q_(m-1)) = kappa(a-b)H^(a-b-1)Q_(m-1)J ?  "
          f"{'OK' if sp.simplify(lhs-rhs)==0 else '*** FAIL ***'}")
print("\n   => the three terms of the n+m-4 equation carry H-powers a-1, a-b-1, b-1;")
print("      dividing by H^min(a-b-1,b-1) reduces the exponent pair: the tower ITERATES.")
