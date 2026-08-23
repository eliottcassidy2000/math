#!/usr/bin/env python3
"""Independent exact hostile checker for THM-3837.

This checker does not use the theorem companion's universal Taylor formula.
It computes both normal jets directly in the truncated ring K[t]/(t^3),
then reduces coefficient identities in Q(alpha)=Q[a]/(3a^3+7a^2+1).
"""

from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path

import sympy as sp


CHECKS = 0


def require(condition: object, label: str) -> None:
    global CHECKS
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(f"FAILED: {label}")
    CHECKS += 1


a = sp.symbols("a")
r_a = 3 * a**3 + 7 * a**2 + 1


def nf_numerator(expression: object) -> sp.Expr:
    """Return the numerator remainder modulo r_a.

    Vanishing of this remainder is equivalent to vanishing in every
    characteristic-zero field containing the chosen root a, once the
    separately audited denominators are units.
    """

    numerator, _ = sp.together(expression).as_numer_denom()
    parameters = sorted(numerator.free_symbols - {a}, key=str)
    domain = sp.QQ.frac_field(*parameters) if parameters else sp.QQ
    return sp.rem(
        sp.Poly(numerator, a, domain=domain),
        sp.Poly(r_a, a, domain=domain),
    ).as_expr()


def same_nf(lhs: object, rhs: object, label: str) -> None:
    require(nf_numerator(sp.together(lhs - rhs)) == 0, label)  # type: ignore[operator]


def same(lhs: object, rhs: object, label: str) -> None:
    require(sp.cancel(sp.expand(lhs - rhs)) == 0, label)  # type: ignore[operator]


# Length-three coefficient lists model the normal parameter modulo t^3.
Jet = list[sp.Expr]


def jet(value0: object = 0, value1: object = 0, value2: object = 0) -> Jet:
    return [sp.sympify(value0), sp.sympify(value1), sp.sympify(value2)]


def jadd(*terms: Jet) -> Jet:
    return [sp.expand(sum(term[i] for term in terms)) for i in range(3)]


def jscale(c: object, term: Jet) -> Jet:
    return [sp.expand(c * term[i]) for i in range(3)]


def jmul(left: Jet, right: Jet) -> Jet:
    return [
        sp.expand(sum(left[i] * right[n - i] for i in range(n + 1)))
        for n in range(3)
    ]


def jpow(term: Jet, exponent: int) -> Jet:
    result = jet(1)
    for _ in range(exponent):
        result = jmul(result, term)
    return result


def selector_residual(h: Jet, k: Jet, d: Jet) -> Jet:
    quadratic = jadd(jscale(7, jpow(h, 2)), jscale(3, jpow(k, 2)))
    cubic = jadd(
        jscale(3, jpow(h, 3)),
        jscale(7, jmul(jpow(h, 2), k)),
        jpow(k, 3),
    )
    a5 = jmul(quadratic, cubic)
    b3 = jmul(
        jmul(jadd(h, k), jadd(jscale(2, h), k)),
        jadd(jscale(3, h), jscale(-1, k)),
    )
    bracket = jadd(jmul(k, b3), jmul(jpow(h, 2), d))
    return jadd(a5, jscale(-1, jmul(d, bracket)))


z = sp.symbols("z")
r = 3 * z**3 + 7 * z**2 + 1
Qz = 7 * z**2 + 3
bz = sp.expand((z + 1) * (2 * z + 1) * (3 * z - 1))
Uz = 49 * z**2 - 9 * z + 7
collision = 140 * z**2 + 307 * z - 93

# Independent arithmetic and all field-unit seams.
same(bz + Qz, 2 * r, "b(z)+Q(z)=2r(z)")
require(sp.gcd(r, Qz) == 1, "Q is coprime to the cubic")
require(sp.gcd(r, Uz) == 1, "U is coprime to the cubic")
require(sp.gcd(r, collision) == 1, "terminal collision is coprime to the cubic")
require(sp.resultant(r, Qz, z) == 1615, "Q resultant")
require(sp.resultant(r, Uz, z) == 14535, "U resultant")
require(sp.resultant(r, collision, z) == -5817545, "collision resultant")
require(r.subs(z, 0) == 1, "alpha cannot vanish")

Q = 7 * a**2 + 3
B = sp.expand(bz.subs(z, a))
rho = 9 * a**2 + 14 * a
Bprime = sp.diff(bz, z).subs(z, a)
mu = -9 * a - 14
U = Uz.subs(z, a)

same_nf(B, -Q, "B=-Q")
same_nf(B, a * rho, "B=a*rho")
same_nf(mu, -B / a**2, "mu=-B/a^2")

# CRT, row completion, and boundary selector data are checked in the full ring.
x, y = sp.symbols("x y")
q_generic = sp.symbols("q_generic")
u = x
v = x * y - 1
f = sp.expand(u * v)
E0 = 1 + x * y * (y - 1)
same(v - y * u, -1, "u and v are comaximal")
same(E0 * (x - v) + f * (y - 1) ** 2, 1, "CRT Bezout identity")

m_generic = E0 * q_generic - (y - 1) ** 2
C_generic = E0 + a * m_generic
k_generic = x - v + f * q_generic
h_generic = a * k_generic + f
same(C_generic * k_generic - m_generic * h_generic, 1,
     "determinant-one row completion")
same(k_generic.subs(x, 0), 1, "line boundary unit")
same(k_generic.subs(y, 1 / x), x, "hyperbola boundary unit")

def b3_scalar(hh: object, kk: object) -> sp.Expr:
    return sp.expand((hh + kk) * (2 * hh + kk) * (3 * hh - kk))


d_boundary = mu * x**2
def a5_scalar(hh: object, kk: object) -> sp.Expr:
    return sp.expand(
        (7 * hh**2 + 3 * kk**2)
        * (3 * hh**3 + 7 * hh**2 * kk + kk**3)
    )


boundary_bracket = k_generic * b3_scalar(h_generic, k_generic) + h_generic**2 * d_boundary
same_nf(a5_scalar(h_generic, k_generic).subs(x, 0), 0,
        "A5 vanishes on the line")
same_nf(a5_scalar(h_generic, k_generic).subs(y, 1 / x), 0,
        "A5 vanishes on the hyperbola")
same_nf(boundary_bracket.subs(x, 0), B,
        "unselected bracket is a unit on the line")
same_nf(d_boundary.subs(y, 1 / x), mu * x**2,
        "selected d is a unit on the hyperbola")
same_nf(
    boundary_bracket.subs(y, 1 / x),
    0,
    "opposite selector bracket vanishes on hyperbola",
)

# Coefficients for general total-degree-at-most-two profiles.
q0, qx, qy, qxx, qxy, qyy = sp.symbols("q0 qx qy qxx qxy qyy")
d0, dx, dy, dxx, dxy, dyy = sp.symbols("d0 dx dy dxx dxy dyy")
Y, X = sp.symbols("Y X", nonzero=True)

# Direct selected-line chart x=t, y=Y.
q_u = jet(q0 + qy * Y + qyy * Y**2, qx + qxy * Y, qxx)
D_u = jet(d0 + dy * Y + dyy * Y**2, dx + dxy * Y, dxx)
f_u = jet(0, -1, Y)
k_u = jadd(jet(1, 1 - Y, 0), jmul(f_u, q_u))
h_u = jadd(jscale(a, k_u), f_u)
d_u = jadd(jet(0, 0, mu), jmul(f_u, D_u))
E_u = selector_residual(h_u, k_u, d_u)

line_boundary_D = d0 + dy * Y + dyy * Y**2
same_nf(E_u[1], -Q * (line_boundary_D + rho),
        "direct line first contact")

# Direct hyperbola chart v=t, x=X, y=(1+t)/X.
q_v0 = q0 + qx * X + qy / X + qxx * X**2 + qxy + qyy / X**2
q_v1 = qy / X + qxy + 2 * qyy / X**2
D_v0 = d0 + dx * X + dy / X + dxx * X**2 + dxy + dyy / X**2
D_v1 = dy / X + dxy + 2 * dyy / X**2
f_v = jet(0, X, 0)
k_v = jadd(jet(X, -1, 0), jmul(f_v, jet(q_v0, q_v1, qyy / X**2)))
h_v = jadd(jscale(a, k_v), f_v)
d_v = jadd(jet(mu * X**2, 0, 0), jmul(f_v, jet(D_v0, D_v1, dyy / X**2)))
E_v = selector_residual(h_v, k_v, d_v)

same(k_v[1] / X, q_v0 - 1 / X, "hyperbola k_f=q-1/X")
hyperbola_contact_rhs = (
    rho + 2 * B / a**2
    - X / a**2 * (2 * B * q_v0 + Bprime - 2 * B / a)
)
expected_v1 = X**5 * (-mu * a**2) * (D_v0 - hyperbola_contact_rhs)
same_nf(E_v[1], expected_v1, "direct hyperbola first contact")

# Verify that coefficient comparison gives exactly the eight advertised
# independent first-contact equations.  These are the Laurent coefficients
# after the three line equations have been imposed.
line_first = {d0: -rho, dy: 0, dyy: 0}
ledger_expression = sp.expand((D_v0 - hyperbola_contact_rhs).subs(line_first))
ledger = {
    -1: 2 * B * qyy / a**2,
    0: dxy - 2 * rho + 2 * B * (qy - 1) / a**2,
    1: dx + (2 * B * (q0 + qxy) + Bprime - 2 * B / a) / a**2,
    2: dxx + 2 * B * qx / a**2,
    3: 2 * B * qxx / a**2,
}
for power, claimed in ledger.items():
    same_nf(sp.expand(ledger_expression).coeff(X, power), claimed,
            f"complete first-contact Laurent coefficient X^{power}")

first_contact = {
    qxx: 0,
    qyy: 0,
    d0: -rho,
    dy: 0,
    dyy: 0,
    qx: -a**2 * dxx / (2 * B),
    qy: 1 + a - a**2 * dxy / (2 * B),
    q0: (-a**2 * dx - Bprime + 2 * B / a) / (2 * B) - qxy,
}
same_nf(E_u[1].subs(first_contact), 0, "ledger kills line first contact")
same_nf(E_v[1].subs(first_contact), 0, "ledger kills hyperbola first contact")

# The two direct second contacts.  No universal Taylor coefficient is reused.
Lu = dx + 18 * a + 38 + 7 * a**2
Lv = dx + 325 * a - 55 + 147 * a**2
line_second = sp.expand(E_u[2].subs(first_contact))
same_nf(line_second, -Q * (dxy * Y + Lu),
        "direct line second contact peels dxy and Lu")

hyperbola_second = sp.expand(E_v[2].subs(first_contact))
coefficient_X8 = hyperbola_second.coeff(X, 8)
same_nf(coefficient_X8, a**2 * dxx**2 / 4,
        "direct hyperbola X^8 second contact peels dxx")

peeled_second = sp.expand(hyperbola_second.subs({dxx: 0, dxy: 0}))
coefficient_X5 = peeled_second.coeff(X, 5)
same_nf(coefficient_X5, -U * Lv / 3,
        "direct hyperbola X^5 second contact forces Lv")
same_nf(sp.diff(coefficient_X5, qxy), 0,
        "surviving qxy cannot change terminal coefficient")
same(Lv - Lu, collision.subs(z, a), "terminal linear-form difference")

# Hostile controls: perturbing either terminal relation is detected, and a
# consistent first-contact specialization really does survive order one.
hostile_dx = -18 * a - 38 - 7 * a**2
same_nf(Lu.subs(dx, hostile_dx), 0, "positive control: line Lu can vanish")
require(nf_numerator(Lv.subs(dx, hostile_dx)) != 0,
        "hostile control: same dx cannot kill Lv")
require(nf_numerator(line_second + Q * (dxy * Y + Lu) + 1) != 0,
        "hostile control: incorrect line identity rejected")

source = Path(__file__).read_text(encoding="utf-8")
require(not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
        "checker contains no Python assert")

semantic = {
    "theorem": "THM-3837",
    "method": "direct K[t]/(t^3) line and hyperbola jets; no universal Taylor reuse",
    "first_contact_rank": 8,
    "survivors": ["dx", "dxx", "dxy", "qxy"],
    "second_contact": "dxy=0, dxx=0, Lu=Lv=0 impossible",
    "collision_resultant": -5817545,
    "scope": "one orientation; n=1; units 1,x; q,D total degree <=2; selector only",
}
semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()

print("audit=THM-3837-independent-direct-contact")
print("method=truncated_normal_ring_direct;universal_Taylor_not_reused")
print("comaximality=PASS;completion=PASS;boundary_signs=PASS")
print("first_contact=8_independent_equations;4_parameters_survive")
print("second_contact=dxy_then_dxx_then_incompatible_Lu_Lv")
print("resultants=1615,14535,-5817545")
print("scope=one_orientation;n1;boundary_units_1_x;degree_le_2;selector_only")
print(f"CHECKS={CHECKS}")
print(f"semantic_sha256={hashlib.sha256(semantic_blob).hexdigest()}")
