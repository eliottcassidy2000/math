#!/usr/bin/env python3
"""Small exact controls for packet 134; not a finite-range proof.

Universe: symbolic rational functions over QQ in x,t and formal time s;
four specified g in t*QQ[t]. No numerical approximation or random sampling.
Run from the packet directory or from any working directory.
"""
from math import factorial
import sympy as sp

x, t, s, r = sp.symbols("x t s r")
checks = 0


def equal(lhs, rhs, label):
    global checks
    assert sp.cancel(lhs - rhs) == 0, label
    checks += 1


def equal_pair(lhs, rhs, label):
    equal(lhs[0], rhs[0], label + ": x")
    equal(lhs[1], rhs[1], label + ": t")


def j_map(pair, time):
    X, T = pair
    d = sp.cancel(1 - time * X**2 * T)
    return sp.cancel(X / d**2), sp.cancel(T * d**3)


def middle_map(pair, g, time=sp.Integer(1)):
    X, T = pair
    d = sp.cancel(1 - X**2 * T)
    v = sp.cancel(T * d**3)
    p = sp.cancel(X + time * d**2 * g.subs(t, v))
    q = sp.cancel(d + T * p**2)
    return sp.cancel(p / q**2), sp.cancel(T * q**3)


def trunc(poly, order):
    p = sp.Poly(sp.expand(poly), t)
    return sp.Add(*(c * t**m[0] for m, c in p.terms() if m[0] < order))


def source_exp(coordinate, H, order):
    value, term = coordinate, coordinate
    for n in range(1, order + 1):
        term = trunc(sp.diff(term, x) * sp.diff(H, t)
                     - sp.diff(term, t) * sp.diff(H, x), order)
        value += term / factorial(n)
    return trunc(value, order)


def inverse_series(poly, order):
    coeff = [sp.expand(poly).coeff(t, n) for n in range(order)]
    assert coeff[0] == 1
    inv = [sp.Integer(1)]
    for n in range(1, order):
        inv.append(sp.expand(-sum(coeff[k] * inv[n-k]
                                  for k in range(1, n+1))))
    return sp.Add(*(c * t**n for n, c in enumerate(inv)))


print("packet134 exact control")
print("SymPy", sp.__version__)
print("Universe: QQ(x,t), symbolic times; g=0,t,-2t+3t^2,t^3")

J = x**3 * t**2
D = 1 - x**2 * t
X, T = j_map((x, t), s)
equal(sp.diff(X, s), 2 * X**3 * T, "J x ODE")
equal(sp.diff(T, s), -3 * X**2 * T**2, "J t ODE")
equal(sp.diff(X, x) * sp.diff(T, t) - sp.diff(X, t) * sp.diff(T, x),
      1, "J Jacobian")
equal_pair(j_map(j_map((x, t), s), r), j_map((x, t), r+s),
           "J scalar group law")
equal_pair(j_map(j_map((x, t), 1), -1), (x, t), "J inverse")
print("PASS: J differential equations, Jacobian, scalar group law, inverse")

choices = [sp.Integer(0), t, -2*t+3*t**2, t**3]
for g in choices:
    f = sp.integrate(g, t)
    H = f.subs(t, t*D**3)
    P = x + s*D**2*g.subs(t, t*D**3)
    Q = D + t*P**2
    bx, bt = P/Q**2, t*Q**3
    equal(sp.diff(bx, s).subs(s, 0), sp.diff(H, t),
          f"middle x vector field g={g}")
    equal(sp.diff(bt, s).subs(s, 0), -sp.diff(H, x),
          f"middle t vector field g={g}")
    equal(sp.diff(H, x), -6*x*t**2*D**2*g.subs(t, t*D**3),
          f"conjugate Hamiltonian x derivative g={g}")
    first = j_map((x, t), -1)
    second = middle_map(first, g)
    final = j_map(second, 1)
    equal_pair(final, (x+g, t), f"complete chronological cancellation g={g}")
    if g != 0:
        p, q = sp.expand(P.subs(s, 1)), sp.expand(Q.subs(s, 1))
        assert sp.gcd(p, q) == 1
        assert sp.Poly(q, x, t).total_degree() > 0
        assert sp.diff(H, x) != 0
        checks += 3
    else:
        equal_pair(middle_map((x, t), g), (x, t), "zero middle identity")
    print(f"PASS: g={g}; exact final map=(x+({g}),t)")
print("PASS: all nonzero middle maps have nonconstant reduced denominator")

# Independent source-series control: repeated Hamiltonian differentiation
# versus exact rational expansion, with no conjugation used by source_exp.
order = 7
g = t
H = (t**2 / 2).subs(t, t*D**3)
p = sp.expand(x + D**2 * g.subs(t, t*D**3))
q = sp.expand(D + t*p**2)
series_x = trunc(p * inverse_series(sp.expand(q**2), order), order)
series_t = trunc(t*q**3, order)
equal(source_exp(x, H, order), series_x, "middle source exponential x")
equal(source_exp(t, H, order), series_t, "middle source exponential t")
print("PASS: g=t middle rational map matches original source exponential mod t^7")

# Hostile order check at an exact rational point; denominators are nonzero.
# Swapping the outer J signs must not accidentally produce the requested shear.
point = (sp.Rational(1), sp.Rational(1, 10))
wrong = j_map(middle_map(j_map(point, 1), t), -1)
target = (point[0] + point[1], point[1])
assert wrong != target
checks += 1
print("PASS: reversed outer order fails at (x,t)=(1,1/10)")
print(f"ALL {checks} EXACT CHECKS PASSED")
print("Scope: supplementary identities only; the uniform theorem is proved in the note.")
