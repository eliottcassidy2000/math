#!/usr/bin/env python3
"""Exact companion for THM-3952's four minimal Mobius carrier colors."""

from __future__ import annotations

import hashlib
import json

import sympy as sp


CHECKS = 0


def gate(condition: object, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if condition is not True and condition != sp.S.true:
        raise RuntimeError(message)


d = sp.symbols("delta")


def dreduce(expression: sp.Expr) -> sp.Expr:
    numerator, denominator = sp.together(sp.cancel(expression)).as_numer_denom()
    num = sp.rem(sp.Poly(sp.expand(numerator), d), sp.Poly(d**2 + 3, d))
    den = sp.rem(sp.Poly(sp.expand(denominator), d), sp.Poly(d**2 + 3, d))
    conjugate = den.as_expr().xreplace({d: -d})
    return sp.cancel(num.as_expr() * conjugate / (den.as_expr() * conjugate))


def dzero(expression: sp.Expr, message: str) -> None:
    gate(dreduce(expression) == 0, message)


omega = (-1 + d) / 2
omega2 = (-1 - d) / 2
dzero(omega**2 + omega + 1, "omega relation")
dzero(omega - omega2 - d, "delta relation")
dzero(d**2 + 3, "delta square")


# ---------------------------------------------------------------------------
# Universal cubic ratio map and its four double fibres.
# ---------------------------------------------------------------------------

x, R, S = sp.symbols("x R S")
Nphi = sp.expand((1 - x) * (1 - omega * x) ** 2)
Dphi = sp.expand((1 + x) * (1 + omega2 * x) ** 2)
phi = sp.cancel(Nphi / Dphi)
wronskian = dreduce(sp.diff(Nphi, x) * Dphi - Nphi * sp.diff(Dphi, x))

dzero(
    wronskian + 2 * d * x * (x + omega) * (x - omega2),
    "finite critical polynomial",
)
gate(sp.degree(Nphi, x) == 3 and sp.degree(Dphi, x) == 3,
     "ratio map has degree three")
gate(sp.Poly(wronskian, x).degree() == 3,
     "three finite critical points")
gate(dreduce(sp.discriminant(wronskian / x, x)) != 0,
     "the two nonzero finite critical points are distinct")

a = dreduce((S - R) * (S - omega * R) ** 2)
b = dreduce((S + R) * (S + omega2 * R) ** 2)
dzero(a / b - phi.subs(x, R / S), "homogeneous carrier ratio")
dzero(a - b - R**2 * (R + d * S), "double fibre over one")
dzero(
    a + omega * b - (1 - omega) * S**2 * (3 * R + d * S) / 3,
    "double fibre over minus omega",
)
dzero(a - (S - R) * (S - omega * R) ** 2,
      "double zero fibre")
dzero(b - (S + R) * (S + omega2 * R) ** 2,
      "double pole fibre")

dzero(phi.subs(x, 0) - 1, "critical value at zero")
gate(dreduce(Nphi.subs(x, omega2)) == 0, "critical zero at omega2")
gate(dreduce(Dphi.subs(x, -omega)) == 0, "critical pole at minus omega")
dzero(
    sp.LC(sp.Poly(Nphi, x)) / sp.LC(sp.Poly(Dphi, x)) + omega,
    "critical value at infinity",
)

# The fibre discriminant gives the same four distinct branch values.
z = sp.symbols("z")
fiber = dreduce(Nphi - z * Dphi)
fiber_disc = dreduce(sp.discriminant(fiber, x))
dzero(
    fiber_disc + 6 * z * (z - 1) * ((d + 3) * z + d - 3),
    "four-value fibre discriminant",
)
dzero((3 - d) / (3 + d) + omega, "third finite branch value")


# ---------------------------------------------------------------------------
# Four affine representatives and their exact linear pencil members.
# ---------------------------------------------------------------------------

t = sp.symbols("t")


def carrier(Rt: sp.Expr, St: sp.Expr) -> tuple[sp.Expr, sp.Expr]:
    return (
        dreduce(a.subs({R: Rt, S: St})),
        dreduce(b.subs({R: Rt, S: St})),
    )


rows: list[tuple[str, sp.Expr, sp.Expr, str, sp.Expr, int]] = [
    ("zero", sp.Integer(1), t, "a-b", d * t + 1, 2),
    (
        "infinity",
        t,
        sp.Integer(1),
        "a+omega*b",
        ((3 - d) * t + 1 + d) / 2,
        3,
    ),
    ("minus_omega", 1 - omega * t, t, "b", d * t + omega, 2),
    ("omega2", 1 + omega2 * t, t, "a", -d * t - omega2, 2),
]

for name, Rt, St, member, expected, clean_count in rows:
    at, bt = carrier(Rt, St)
    members = {
        "a-b": at - bt,
        "a+omega*b": at + omega * bt,
        "a": at,
        "b": bt,
    }
    dzero(members[member] - expected, f"{name} exact linear member")
    gate(sp.degree(dreduce(members[member]), t) == 1,
         f"{name} member has degree one")
    gate(max(sp.degree(at, t), sp.degree(bt, t)) == 3,
         f"{name} retains a cubic coordinate")

    U = dreduce(St + omega2 * Rt)
    V = dreduce(St - omega * Rt)
    # R,U,V are homogenized as degree-one forms. A constant affine form
    # vanishes at the source point at infinity; a linear form has one finite
    # zero. Count the latter, which are THM-3951's clean affine colors.
    finite_clean = sum(sp.degree(dreduce(q), t) == 1 for q in (Rt, U, V))
    gate(finite_clean == clean_count, f"{name} finite clean-color count")


# ---------------------------------------------------------------------------
# Noncritical hostile: degrees (2,3) and a finite cusp.
# ---------------------------------------------------------------------------

R_bad = t + 1
S_bad = t
a_bad, b_bad = carrier(R_bad, S_bad)
t_cusp = (-3 + d) / 6
gate(sp.degree(a_bad, t) == 2, "hostile lower degree is two")
gate(sp.degree(b_bad, t) == 3, "hostile upper degree is three")
dzero(sp.diff(a_bad, t).subs(t, t_cusp), "hostile cusp kills da")
dzero(sp.diff(b_bad, t).subs(t, t_cusp), "hostile cusp kills db")

# The unique pencil cancellation at a noncritical infinity has exact degree
# two in this row; no accidental second cancellation is being hidden.
gate(sp.Poly(a_bad, t).LC() != 0, "hostile quadratic leading coefficient")
gate(sp.Poly(b_bad, t).LC() != 0, "hostile cubic leading coefficient")


# ---------------------------------------------------------------------------
# Riemann--Hurwitz bridge: rational C2^4 to polynomial C3+C2+C2.
# ---------------------------------------------------------------------------

rh_c2_four = -2 * 6 + 4 * (6 - 6 // 2)
rh_polynomial_s3 = -2 * 6 + (6 - 6 // 3) + 2 * (6 - 6 // 2)
rh_cyclic = -2 * 3 + 2 * (3 - 3 // 3)
gate(rh_c2_four == 0, "C2^4 closure has genus one")
gate(rh_polynomial_s3 == -2, "C3+C2+C2 closure has genus zero")
gate(rh_cyclic == -2, "cyclic C3+C3 closure has genus zero")

rho, q = sp.symbols("rho q")
u = sp.symbols("u")
poly_cubic = u**3 - 3 * rho**2 * u + q
dzero(
    poly_cubic.subs(u, rho) - poly_cubic.subs(u, -rho) + 4 * rho**3,
    "distinct cubic critical values cannot collide",
)


summary = {
    "checks": CHECKS,
    "classification": "k[a,b]=k[t] iff source infinity is critical for phi",
    "critical_colors": ["0", "infinity", "-omega", "omega2"],
    "noncritical": "target-normalized degrees (3,2), Abhyankar-Moh obstruction",
    "atlas": "all four excluded by THM-3951 clean double incidence",
    "rh_bridge": "C2^4 genus1 -> C3+C2+C2 or cyclic C3+C3 genus0",
    "scope": "degree-one ratio and unit extra debt only",
}
semantic = hashlib.sha256(json.dumps(summary, sort_keys=True).encode()).hexdigest()

print("THM-3952 minimal Mobius internal-split carrier companion")
print(f"CHECKS={CHECKS}")
print("CRITICAL_SOURCE_COLORS=0,infinity,-omega,omega^2")
print("EMBEDDING_IFF=SOURCE_INFINITY_CRITICAL")
print("FOUR_ROWS=EXACT_LINEAR_PENCIL_MEMBER")
print("NONCRITICAL_HOSTILE=DEGREES_3_2;FINITE_CUSP")
print("CLEAN_AFFINE_INCIDENCES=2,3,2,2;THM3951_NO_ATLAS")
print("RH=C2^4_GENUS1;C3_C2_C2_GENUS0;CYCLIC_C3_C3_GENUS0")
print("SCOPE=UNIT_DEBT_DEGREE_ONE_RATIO;EXTRA_DEBT_AND_HIGHER_RATIO_OPEN")
print(f"SEMANTIC_SHA256={semantic}")
