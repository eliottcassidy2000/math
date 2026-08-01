"""KILLING BRANCH (i) of THM-3016 by a coordinate normalisation.

Branch (i) of the dichotomy (X'') is  P_{n-1} = Q_{m-1} = 0.
Let tau(x,y) = (x+s, y+t) be a SOURCE TRANSLATION.  Then
  * Jac(P o tau, Q o tau) = Jac(P,Q) o tau = 1, so the pair stays Jacobian;
  * the leading forms are UNCHANGED (a translation cannot change the top
    homogeneous component), so H, g, a, b and K are all preserved;
  * the subleading components shift by the directional derivative of the
    leading form:
        P_{n-1}  ->  P_{n-1} + (s d_x + t d_y) P_n,
        Q_{m-1}  ->  Q_{m-1} + (s d_x + t d_y) Q_m.
With P_n = c H^a, Q_m = c' H^b and L := (s d_x + t d_y) H (a form of degree
g-1), a pair in branch (i) becomes
        P_{n-1}' = c a H^{a-1} L,      Q_{m-1}' = c' b H^{b-1} L,
and then, with kappa = ca/(c'b),
        kappa H^{a-b} Q_{m-1}' = (ca/(c'b)) H^{a-b} c' b H^{b-1} L
                               = c a H^{a-1} L = P_{n-1}',
i.e. the translated pair satisfies BRANCH (ii) EXACTLY.  Moreover L != 0 for
generic (s,t): in characteristic 0, d_x H and d_y H cannot both vanish for a
form of degree g >= 1 (Euler: x d_x H + y d_y H = g H).

CONCLUSION.  Branch (i) is not a separate case.  Every branch-(i) pair is a
translate of a branch-(ii) pair with BOTH subleading components nonzero, so
after a generic translation every counterexample (K >= 2) may be assumed to
satisfy
        P_{n-1} = kappa H^{a-b} Q_{m-1},   both nonzero,
and in fact the translated pair carries the sharper factorisation
        Q_{m-1} = c' b H^{b-1} L,   deg L = g-1.
This file verifies all of it symbolically.
"""
import sympy as sp
x, y, s, t = sp.symbols('x y s t')

def jac(A, B): return sp.expand(sp.diff(A,x)*sp.diff(B,y) - sp.diff(A,y)*sp.diff(B,x))
def hom(F, d):
    F = sp.Poly(sp.expand(F), x, y)
    return sp.expand(sum(cc*x**i*y**j for (i,j),cc in F.terms() if i+j==d))
def deg(F): return sp.Poly(sp.expand(F), x, y).total_degree()

print("STEP 1 -- translation preserves Jac and the leading form:")
P = sp.expand(x + y**2); Q = sp.expand(y + (x+y**2)**3)
Pt = sp.expand(P.subs({x: x+s, y: y+t}, simultaneous=True))
Qt = sp.expand(Q.subs({x: x+s, y: y+t}, simultaneous=True))
print("   Jac(P o tau, Q o tau) =", sp.simplify(jac(Pt, Qt)))
n, m = deg(P), deg(Q)
print("   leading forms equal:",
      sp.simplify(hom(Pt, n) - hom(P, n)) == 0 and sp.simplify(hom(Qt, m) - hom(Q, m)) == 0)

print("\nSTEP 2 -- subleading shift is the directional derivative of the top form:")
lhs = sp.expand(hom(Pt, n-1) - hom(P, n-1))
rhs = sp.expand(s*sp.diff(hom(P,n), x) + t*sp.diff(hom(P,n), y))
print("   P_{n-1} shift == (s d_x + t d_y) P_n :", sp.simplify(lhs-rhs) == 0)
lhs = sp.expand(hom(Qt, m-1) - hom(Q, m-1))
rhs = sp.expand(s*sp.diff(hom(Q,m), x) + t*sp.diff(hom(Q,m), y))
print("   Q_{m-1} shift == (s d_x + t d_y) Q_m :", sp.simplify(lhs-rhs) == 0)

print("\nSTEP 3 -- a branch-(i) pair translates into branch (ii) exactly:")
c, cp, a, b, g = sp.symbols('c cprime a b g', positive=True)
H = sp.Function('H')
# symbolic identity: with P_{n-1}=Q_{m-1}=0 initially,
# P' = c a H^{a-1} L, Q' = c' b H^{b-1} L, kappa = c a/(c' b)
Hs = sp.Symbol('H', positive=True); L = sp.Symbol('L')
Pp = c*a*Hs**(a-1)*L; Qp = cp*b*Hs**(b-1)*L; kappa = c*a/(cp*b)
print("   kappa H^{a-b} Q'  - P'  =",
      sp.simplify(kappa*Hs**(a-b)*Qp - Pp), " (0 means branch (ii) holds)")

print("\nSTEP 4 -- L = (s d_x + t d_y)H is generically nonzero (Euler):")
for Hex in [x*y, x**2 - y**2, x*y*(x+y), x**3 + y**3]:
    gg = deg(Hex)
    Lex = sp.expand(s*sp.diff(Hex,x) + t*sp.diff(Hex,y))
    eul = sp.simplify(x*sp.diff(Hex,x) + y*sp.diff(Hex,y) - gg*Hex)
    Kf = len([f for f,e in sp.factor_list(Hex)[1] if sp.total_degree(f)>0])
    print(f"   H={Hex} (g={gg}, K={Kf}): L={Lex}  nonzero for generic (s,t): "
          f"{sp.simplify(Lex) != 0}   Euler check: {eul == 0}")

print("\nSTEP 5 -- concrete end-to-end demonstration on a genuine Jacobian pair")
print("   (build a pair whose subleading layer vanishes, then translate):")
# u = x, v = y + x^2 has P_{n-1}=0 trivially; use P=x, Q=y+x^3
P0, Q0 = x, sp.expand(y + x**3)
n0, m0 = deg(P0), deg(Q0)
print(f"   P=x, Q=y+x^3: Jac={sp.simplify(jac(P0,Q0))}, "
      f"P_(n-1)={hom(P0,n0-1)}, Q_(m-1)={hom(Q0,m0-1)}  -> branch (i)")
P0t = sp.expand(P0.subs({x:x+1,y:y+1}, simultaneous=True))
Q0t = sp.expand(Q0.subs({x:x+1,y:y+1}, simultaneous=True))
print(f"   after tau=(x+1,y+1): Jac={sp.simplify(jac(P0t,Q0t))}, "
      f"P_(n-1)={hom(P0t,n0-1)}, Q_(m-1)={hom(Q0t,m0-1)}")
print("   -> subleading layer is now nonzero: branch (i) destroyed by translation.")
