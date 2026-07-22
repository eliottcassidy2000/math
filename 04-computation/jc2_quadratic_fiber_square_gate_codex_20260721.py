#!/usr/bin/env python3
"""Exact symbolic referee for THM-2071.

This is a check of the displayed coefficient identities and normal form, not
an exhaustive search and not evidence for unrestricted JC(2).

Reproduce from the repository root with:
    python3 04-computation/jc2_quadratic_fiber_square_gate_codex_20260721.py
"""

from sympy import (
    Function,
    Poly,
    Rational,
    S,
    diff,
    expand,
    integrate,
    linsolve,
    simplify,
    symbols,
)


x, y, u, v, t = symbols("x y u v t")


def jacobian(p, q):
    return expand(diff(p, x) * diff(q, y) - diff(p, y) * diff(q, x))


print("THM-2071 exact symbolic referee")
print("scope: coefficient identities + normal form + bounded hostile controls")

# 1. The top fiber coefficient for arbitrary complementary fiber degree.
A = Function("A")(x)
B = Function("B")(x)
C = Function("C")(x)
P = A * y**2 + B * y + C
top_checks = []
for n in range(1, 8):
    qs = [Function(f"q{n}_{j}")(x) for j in range(n + 1)]
    Q = sum(qs[j] * y**j for j in range(n + 1))
    observed = expand(jacobian(P, Q)).coeff(y, n + 1)
    expected = n * diff(A, x) * qs[n] - 2 * A * diff(qs[n], x)
    top_checks.append(simplify(observed - expected) == 0)
assert all(top_checks)
print(f"top-coefficient law n=1..7: PASS {top_checks}")
print("  [y^(n+1)] J(P,Q) = n*A'*q_n - 2*A*q_n'")

# 2. The even-degree target shear cancels exactly and preserves the bracket.
d = symbols("d", nonzero=True)
even_checks = []
for k in range(1, 5):
    lead = expand(P**k).coeff(y, 2 * k)
    bracket_preserved = simplify(jacobian(P, d * P**k)) == 0
    even_checks.append(simplify(lead - A**k) == 0 and bracket_preserved)
assert all(even_checks)
print(f"even top cancellation k=1..4: PASS {even_checks}")

# 3. Exact coefficient system when the reduced complementary degree is one.
L = Function("L")(x)
M = Function("M")(x)
Q1 = L * y + M
J1 = expand(jacobian(P, Q1))
expected_coeffs = {
    2: diff(A, x) * L - 2 * A * diff(L, x),
    1: diff(B, x) * L - 2 * A * diff(M, x) - B * diff(L, x),
    0: diff(C, x) * L - B * diff(M, x),
}
linear_checks = [
    simplify(J1.coeff(y, degree) - expression) == 0
    for degree, expression in expected_coeffs.items()
]
assert all(linear_checks)
print(f"reduced-degree-one coefficient system: PASS {linear_checks}")

# 4. If A is constant, completion of the square gives P=A*Y^2+D(x).
# The odd coefficient chain has
#   q_(j-2)' = (j/(2A))*D'*q_j.
# Thus q_1 is a polynomial in D of degree (n-1)/2, and its antiderivative
# has degree (n+1)/2.  For n>=3 it cannot compose with nonconstant D to an
# affine polynomial in x, as the constant bracket would require.
z = symbols("z")
antiderivative_degrees = []
for n in (3, 5, 7):
    q_odd = 1
    for j in range(n, 1, -2):
        integration_constant = symbols(f"k{n}_{j}")
        q_odd = integrate(Rational(j, 4) * q_odd, z) + integration_constant
    primitive = integrate(q_odd, z)
    assert Poly(q_odd, z).degree() == (n - 1) // 2
    assert Poly(primitive, z).degree() == (n + 1) // 2
    antiderivative_degrees.append(Poly(primitive, z).degree())
print(
    "constant-leading odd-tail obstruction n=3,5,7: PASS "
    f"primitive degrees {antiderivative_degrees}"
)

# 5. A nontrivial exact member of the resulting tame normal form, together
# with the displayed inverse in both directions.
A0 = Rational(3, 2)
ell = Rational(5, 3)
kappa = Rational(7, 4)
d0 = Rational(-2, 5)
m0 = Rational(11, 6)
B0 = 2 * x + 1
H = lambda z: z**2 - 2 * z + Rational(4, 7)
Y = y + B0 / (2 * A0)
P0 = expand(A0 * Y**2 + (kappa / ell) * x + d0)
Q0 = expand(H(P0) + ell * Y + m0)
assert simplify(jacobian(P0, Q0) - kappa) == 0

Yinv = (v - H(u) - m0) / ell
xinv = expand((ell / kappa) * (u - A0 * Yinv**2 - d0))
yinv = expand(Yinv - B0.subs(x, xinv) / (2 * A0))
source_roundtrip = (
    simplify(xinv.subs({u: P0, v: Q0}, simultaneous=True) - x),
    simplify(yinv.subs({u: P0, v: Q0}, simultaneous=True) - y),
)
target_roundtrip = (
    simplify(P0.subs({x: xinv, y: yinv}, simultaneous=True) - u),
    simplify(Q0.subs({x: xinv, y: yinv}, simultaneous=True) - v),
)
assert source_roundtrip == (0, 0)
assert target_roundtrip == (0, 0)
print("explicit tame normal form and two-sided inverse: PASS")

# 6. Hostile bounded controls.  The theorem proves the nonsquare case
# impossible at every degree; these finite systems merely catch sign or
# coefficient mistakes.  The square-but-nonconstant row is a control that the
# square condition is necessary, not advertised as sufficient.
P_nonsquare = x * y**2 + (x + 1) * y + x**2
P_square_nonconstant = x**2 * y**2 + (x + 1) * y + x**2
hostile = []
for p in (P_nonsquare, P_square_nonconstant):
    # linsolve returns EmptySet for inconsistency.  Spell the comparison in a
    # version-stable way through bool(solution).
    coeffs = symbols("h0:24")
    q = sum(coeffs[4 * j + i] * x**i * y**j for j in range(6) for i in range(4))
    residual = Poly(expand(jacobian(p, q) - 1), x, y)
    solution = linsolve(list(residual.coeffs()), coeffs)
    hostile.append(solution is S.EmptySet)
assert all(hostile)
print(f"bounded hostile mate searches (deg_x<=3, deg_y<=5): PASS {hostile}")

print("RESULT: PASS")
