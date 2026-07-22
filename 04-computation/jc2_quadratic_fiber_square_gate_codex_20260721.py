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
    binomial,
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

# 4. Centered coordinates P=z^2+D turn the bracket into
#   U*(D'*partial_z-2*z*partial_x).  Check the complete cubic system and its
# translation back to the original fiber coefficient q_1.
z = symbols("z")
Dcenter = Function("Dcenter")(x)
Ucenter = Function("Ucenter")(x)
E = Function("E")(x)
F = Function("F")(x)
c, e, f = symbols("c e f", nonzero=True)


def centered_operator(q):
    return expand(diff(Dcenter, x) * diff(q, z) - 2 * z * diff(q, x))


Qcubic = c * z**3 + E * z + F
Lcubic = centered_operator(Qcubic)
cubic_checks = (
    simplify(Lcubic.coeff(z, 2) - (3 * c * diff(Dcenter, x) - 2 * diff(E, x)))
    == 0,
    simplify(Lcubic.coeff(z, 1) + 2 * diff(F, x)) == 0,
    simplify(Lcubic.coeff(z, 0) - diff(Dcenter, x) * E) == 0,
)
assert all(cubic_checks)
Enormal = Rational(3, 2) * c * Dcenter + e
Qnormal = c * z**3 + Enormal * z + f
assert simplify(centered_operator(Qnormal) - diff(Dcenter, x) * Enormal) == 0

hcenter = B / (2 * Ucenter)
Doriginal = C - hcenter**2
zoriginal = Ucenter * y + hcenter
Qoriginal = expand(c * zoriginal**3 + (Rational(3, 2) * c * Doriginal + e) * zoriginal + f)
expected_q1 = Ucenter * (Rational(3, 2) * c * C + e) + Rational(3, 8) * c * B**2 / Ucenter
expected_q0 = -Rational(1, 2) * c * hcenter**3 + (Rational(3, 2) * c * C + e) * hcenter + f
assert simplify(Qoriginal.coeff(y, 1) - expected_q1) == 0
assert simplify(Qoriginal.coeff(y, 0) - expected_q0) == 0
print(f"centered cubic system and q_1/q_0 translation: PASS {cubic_checks}")

# 5. The full odd recurrence has
#   2*d(a_(j-1))/dD=(2j+1)*a_j.
# Check its leading coefficients and the nonzero polar coefficient in q_0
# exactly through r=6.  The proof uses the same closed formula for every r.
T, hvar, cvar = symbols("T hvar cvar")
pole_coefficients = []
for r in range(7):
    a = [None] * (r + 1)
    a[r] = c
    for j in range(r, 0, -1):
        integration_constant = symbols(f"lambda_{r}_{j}")
        a[j - 1] = (
            integrate(Rational(2 * j + 1, 2) * a[j], T)
            + integration_constant
        )

    for k in range(r + 1):
        j = r - k
        expected_lead = c * binomial(Rational(2 * r + 1, 2), k)
        assert simplify(Poly(a[j], T).LC() - expected_lead) == 0

    q_at_y_zero = expand(
        sum(a[j].subs(T, cvar - hvar**2) * hvar ** (2 * j + 1) for j in range(r + 1))
    )
    observed_pole = Poly(q_at_y_zero, hvar).coeff_monomial(hvar ** (2 * r + 1))
    alternating_sum = sum(
        (-1) ** k * binomial(Rational(2 * r + 1, 2), k)
        for k in range(r + 1)
    )
    expected_pole = c * (-1) ** r * binomial(2 * r, r) / 4**r
    assert simplify(observed_pole - c * alternating_sum) == 0
    assert simplify(observed_pole - expected_pole) == 0
    pole_coefficients.append(simplify(observed_pole / c))
print(
    "centered odd-tail noncancellation r=0..6: PASS "
    f"coefficients {pole_coefficients}"
)

# 6. A nontrivial exact member of the resulting tame normal form, together
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

# 7. Hostile bounded controls.  The theorem proves both nonconstant-leading
# cases impossible at every degree; these finite systems merely catch sign or
# coefficient mistakes and do not supply the uniform proof.
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
