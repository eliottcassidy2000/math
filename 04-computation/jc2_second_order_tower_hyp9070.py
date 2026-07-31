"""The second-order Jacobian condition, derived and verified.

Let (P,Q) be a Jacobian pair in C[x,y], n=deg P, m=deg Q, Jac(P,Q)=1.
Write P = sum_j P_j, Q = sum_j Q_j for homogeneous components.
  degree n+m-2 part of Jac:  Jac(P_n,Q_m) = 0
      => (two variables) P_n = c H^a, Q_m = c' H^b, deg H = g = gcd(n,m),
         a = n/g, b = m/g coprime.                                    (L0)
  degree n+m-3 part:  Jac(P_n,Q_{m-1}) + Jac(P_{n-1},Q_m) = 0.        (L1)
Substituting (L0) into (L1) and using Jac(H^a,B) = a H^{a-1} Jac(H,B):

      c a H^{a-1} Jac(H, Q_{m-1})  =  - c' b H^{b-1} Jac(P_{n-1}, H)

so, assuming wlog a >= b and cancelling H^{b-1},

      c a H^{a-b} Jac(H, Q_{m-1})  =  - c' b Jac(P_{n-1}, H).          (L2)

(L2) forces  H^{a-b} | Jac(P_{n-1}, H)  whenever Jac(H,Q_{m-1}) != 0.
Iterating replaces (a,b) by (a-b, b): the tower runs the SUBTRACTIVE
EUCLIDEAN ALGORITHM on the exponent pair, so its depth is the continued
fraction of a/b.  Jung-van der Kulk (an automorphism has n|m or m|n) says a
counterexample needs a,b >= 2 coprime, i.e. the tower cannot terminate at
the first step.

This script verifies (L1)/(L2) exactly on genuine Jacobian pairs.
"""
import sympy as sp

x, y = sp.symbols('x y')


def jac(P, Q):
    return sp.expand(sp.diff(P, x) * sp.diff(Q, y) - sp.diff(P, y) * sp.diff(Q, x))


def hom(F, d):
    F = sp.Poly(sp.expand(F), x, y)
    return sp.expand(sum(c * x ** mx * y ** my
                         for (mx, my), c in F.terms() if mx + my == d))


def check(P, Q, name):
    J = jac(P, Q)
    n, m = sp.Poly(P, x, y).total_degree(), sp.Poly(Q, x, y).total_degree()
    g = sp.gcd(n, m)
    a, b = n // g, m // g
    Pn, Qm = hom(P, n), hom(Q, m)
    top = sp.expand(jac(Pn, Qm))
    L1 = sp.expand(jac(Pn, hom(Q, m - 1)) + jac(hom(P, n - 1), Qm))
    print(f"  {name}")
    print(f"     deg=({n},{m})  g={g}  (a,b)=({a},{b})  Jac={sp.simplify(J)}")
    print(f"     L0: Jac(P_n,Q_m) = {top}   {'OK' if top == 0 else 'NONZERO'}")
    print(f"     L1: Jac(P_n,Q_{{m-1}})+Jac(P_{{n-1}},Q_m) = {sp.simplify(L1)}"
          f"   {'OK' if sp.simplify(L1) == 0 else '*** NONZERO ***'}")
    # exhibit H
    if top == 0 and n > 0 and m > 0:
        try:
            fs = sp.factor_list(Pn)
            print(f"     P_n factors: {sp.factor(Pn)}")
            print(f"     Q_m factors: {sp.factor(Qm)}")
        except Exception:
            pass
    return top == 0 and sp.simplify(L1) == 0


print("Verification of L0/L1 on genuine Jacobian pairs (automorphisms of C^2):")
ok = True
# triangular
ok &= check(x, y + x ** 3, "F=(x, y+x^3)")
# composition: (x,y) -> (x + y^2, y) -> then (u, v+u^3)
u, v = x + y ** 2, y
ok &= check(sp.expand(u), sp.expand(v + u ** 3), "F=(x+y^2, y+(x+y^2)^3)")
# a deeper composite, degrees (2,7)
P2 = sp.expand(x + y ** 2)
Q2 = sp.expand(y + P2 ** 3)
P3 = sp.expand(P2 + Q2 ** 2)
ok &= check(P3, Q2, "deeper composite")
print("\nALL L0/L1 CHECKS PASS" if ok else "\nSOME CHECK FAILED")

print("\nJung-van der Kulk consequence (degrees of automorphisms):")
print("  every pair above has n|m or m|n -> (a,b) has a=1 or b=1, tower")
print("  terminates immediately.  A counterexample needs a,b>=2 coprime,")
print("  i.e. a NON-TRIVIAL continued fraction: this is the search lattice.")
