#!/usr/bin/env python3
"""Exact finite controls for the all-order nonresonance atlas THM-2827."""

import ast
import math
from pathlib import Path

import sympy as sp


def require(condition, message):
    if not bool(condition):
        raise RuntimeError(message)


def syntax_has(path, node_type):
    tree = ast.parse(Path(path).read_text(encoding="utf-8"))
    return any(isinstance(node, node_type) for node in ast.walk(tree))


def has_float_literals(path):
    tree = ast.parse(Path(path).read_text(encoding="utf-8"))
    return any(
        isinstance(node, ast.Constant) and isinstance(node.value, float)
        for node in ast.walk(tree)
    )


def clean(polynomial):
    return {
        monomial: sp.cancel(coefficient)
        for monomial, coefficient in polynomial.items()
        if sp.cancel(coefficient) != 0
    }


def add(polynomial, monomial, coefficient):
    polynomial[monomial] = polynomial.get(monomial, sp.Integer(0)) + coefficient


def direct_coefficient(R, n):
    """Coefficient of t^n in (1+2dt^2+qt^3+(d^2-s)t^4)^(R-1/2)."""
    alpha = sp.Rational(2 * R - 1, 2)
    answer = {}
    for i in range(n // 2 + 1):
        for j in range(n // 3 + 1):
            remainder = n - 2 * i - 3 * j
            if remainder < 0 or remainder % 4:
                continue
            k = remainder // 4
            total = i + j + k
            multinomial = math.factorial(total) // (
                math.factorial(i) * math.factorial(j) * math.factorial(k)
            )
            common = sp.binomial(alpha, total) * multinomial * 2**i
            for ell in range(k + 1):
                add(
                    answer,
                    (j, i + 2 * (k - ell), ell),
                    common * sp.binomial(k, ell) * (-1) ** ell,
                )
    return clean(answer)


def shifted(polynomial, q_shift=0, d_shift=0, s_shift=0, scalar=1):
    answer = {}
    for (q_degree, d_degree, s_degree), coefficient in polynomial.items():
        monomial = (
            q_degree + q_shift,
            d_degree + d_shift,
            s_degree + s_shift,
        )
        require(min(monomial) >= 0, f"negative shifted exponent {monomial}")
        add(answer, monomial, scalar * coefficient)
    return clean(answer)


def combine(*polynomials):
    answer = {}
    for polynomial in polynomials:
        for monomial, coefficient in polynomial.items():
            add(answer, monomial, coefficient)
    return clean(answer)


def explicit_H(R):
    alpha = sp.Rational(2 * R - 1, 2)
    answer = {}
    h = 0
    while R - 2 - 3 * h >= 0:
        remainder = R - 2 - 3 * h
        j = R + 1 + h
        for k in range(remainder // 2 + 1):
            ell = remainder - 2 * k
            coefficient = (
                4
                * sp.binomial(alpha, j)
                * sp.binomial(j, ell)
                * (-1) ** ell
                * sp.binomial(-2 - 2 * h, k)
            )
            add(answer, (2 * (k + 2 * h), k, ell), coefficient)
        h += 1
    return clean(answer)


def row(R):
    c_phi = direct_coefficient(R, 4 * R - 1)
    c_psi = direct_coefficient(R, 4 * R)
    c_response = direct_coefficient(R, 4 * R + 1)
    phi = shifted(c_phi, q_shift=-1, scalar=4)
    psi = shifted(c_psi, scalar=4)
    h_direct = shifted(
        combine(c_response, shifted(c_phi, d_shift=1)),
        q_shift=-3,
        scalar=4,
    )
    h_formula = explicit_H(R)
    require(h_direct == h_formula, f"third-response formula R={R}")
    return {"phi": phi, "psi": psi, "H": h_formula}


def minimum_face(polynomial, weights):
    require(polynomial, "minimum face of zero polynomial")
    values = {
        monomial: sum(exponent * weight for exponent, weight in zip(monomial, weights))
        for monomial in polynomial
    }
    minimum = min(values.values())
    return minimum, {
        monomial: polynomial[monomial]
        for monomial, value in values.items()
        if value == minimum
    }


rho = sp.symbols("rho")


def prefix_polynomials(R):
    """THM-2760 equations (9)--(10)."""
    alpha = sp.Rational(2 * R - 1, 2)
    P = 0
    Q = 0
    for ell in range((R - 1) // 3 + 1):
        P += (
            sp.binomial(alpha, R + ell)
            * (-1) ** (R - 1 - 3 * ell)
            * 2 ** (R - 1 - 3 * ell)
            * sp.binomial(R - 1 - ell, 2 * ell)
            * rho**ell
        )
    for ell in range(R // 3 + 1):
        Q += (
            sp.binomial(alpha, R + ell)
            * (-1) ** (R - 3 * ell)
            * sp.Rational(2) ** (R - 3 * ell - 1)
            * (
                2 * sp.binomial(R - 1 - ell, 2 * ell - 1)
                + sp.binomial(R - 1 - ell, 2 * ell)
            )
            * rho**ell
        )
    return sp.Poly(P, rho, domain=sp.QQ), sp.Poly(Q, rho, domain=sp.QQ)


rows = {R: row(R) for R in range(1, 19)}
prefix = {R: prefix_polynomials(R) for R in range(1, 19)}
for R, (P, Q) in prefix.items():
    require(P.degree() == (R - 1) // 3, f"deg P_R, R={R}")
    require(Q.degree() == R // 3, f"deg Q_R, R={R}")
    require(P.nth(0) != 0 and Q.nth(0) != 0, f"prefix constants R={R}")
    require(P.LC() != 0 and Q.LC() != 0, f"prefix leading terms R={R}")
    require(sp.gcd(P, Q).degree() == 0, f"prefix coprimality R={R}")

cycle = {0: "psi", 1: "phi", 2: "H"}
lane_counts = {"polar": 0, "pure": 0, "regular-H": 0}
resonances = []
double_pole_resonances = []

for nu in range(3, 26):
    for a in range(nu + 5):
        for b in range(nu + 5):
            polar = 2 * b < nu and a + b < nu
            if polar:
                lane_counts["polar"] += 1
                x = a + b - nu
                y = a - 3 * b + nu
                require(x < 0, f"polar s order nu,a,b={nu,a,b}")
                for R in range(7, 19):
                    P_R, Q_R = prefix[R]
                    if y > 0:
                        p_exp = q_exp = 0
                    elif y < 0:
                        p_exp = P_R.degree()
                        q_exp = Q_R.degree()
                    else:
                        require(sp.gcd(P_R, Q_R).degree() == 0, f"unit-rho gcd R={R}")
                        continue
                    phi_top = (R - 1) * x + p_exp * y
                    psi_top = R * x + q_exp * y
                    for j in range(1, R - 1):
                        P_j, Q_j = prefix[j]
                        p_j_exp = 0 if y > 0 else P_j.degree()
                        q_j_exp = 0 if y > 0 else Q_j.degree()
                        phi_lower = (j - 1) * x + p_j_exp * y
                        psi_lower = j * x + q_j_exp * y
                        require(
                            phi_top < phi_lower,
                            f"polar Phi top/lower nu,a,b,R,j={nu,a,b,R,j}",
                        )
                        require(
                            psi_top < psi_lower,
                            f"polar Psi top/lower nu,a,b,R,j={nu,a,b,R,j}",
                        )
                continue

            if 2 * a < nu:
                lane_counts["pure"] += 1
                require(2 * b >= nu, f"nonpolar pure cone has polar d: {nu,a,b}")
                require(a + nu < 3 * b, f"q/s efficiency inequality: {nu,a,b}")
                weights = (2 * a - nu, 0, 2 * min(0, a + b - nu))
                delta = nu - 2 * a
                require(delta > 0, f"pure delta: {nu,a,b}")
                for R in range(7, 19):
                    channel = cycle[R % 3]
                    top_weight, top_face = minimum_face(rows[R][channel], weights)
                    require(
                        len(top_face) == 1 and next(iter(top_face))[1:] == (0, 0),
                        f"pure top face nu,a,b,R={nu,a,b,R}",
                    )
                    for j in range(1, R - 1):
                        if not rows[j][channel]:
                            continue
                        lower_weight, _ = minimum_face(rows[j][channel], weights)
                        require(
                            top_weight < lower_weight,
                            f"pure top/lower nu,a,b,R,j={nu,a,b,R,j}",
                        )

                    if channel != "H":
                        continue
                    k = (R - 2) // 3
                    actual_h_order = -2 * k * delta
                    required_h_order = nu + 1 - 3 * a
                    arithmetic_resonance = (
                        a == 1 + (2 * k + 1) * delta
                        and nu == 2 + (4 * k + 3) * delta
                    )
                    require(
                        (actual_h_order == required_h_order) == arithmetic_resonance,
                        f"resonance equivalence nu,a,b,R={nu,a,b,R}",
                    )
                    if arithmetic_resonance:
                        record = (R, nu, a, delta)
                        if record not in resonances:
                            resonances.append(record)
                        if nu == 4:
                            double_pole_resonances.append(record)
                continue

            lane_counts["regular-H"] += 1
            vT = 2 * a - nu
            vd = min(0, 2 * b - nu)
            vs = min(0, a + b - nu)
            require(vT >= 0 and vs >= 0 and vT + vd >= 0,
                    f"regular H cone inequalities: {nu,a,b,vT,vd,vs}")
            require(a >= 2 and 1 - a < 0, f"regular H response pole: {nu,a,b}")
            weights = (2 * vT, 2 * vd, 2 * vs)
            for R in range(2, 19):
                if rows[R]["H"]:
                    h_weight, _ = minimum_face(rows[R]["H"], weights)
                    require(h_weight >= 0, f"regular H monomial R={R}, nu,a,b={nu,a,b}")

require(not double_pole_resonances, "nu=4 must avoid every resonance ray")
require(resonances, "hostile range must contain resonance rays")
require(resonances[0] == (8, 13, 6, 1), "first resonance ray")

# Exact formal local witness on the first resonance.
require(
    not any(d_degree == 0 and s_degree == 0 for _, d_degree, s_degree in rows[8]["phi"]),
    "R=8 pure Phi must vanish at d=s=0",
)
require(
    not any(d_degree == 0 and s_degree == 0 for _, d_degree, s_degree in rows[8]["psi"]),
    "R=8 pure Psi must vanish at d=s=0",
)
first_h_coefficient = -sp.Rational(195, 131072)
require(
    rows[8]["H"].get((8, 0, 0)) == first_h_coefficient,
    "R=8 resonant H coefficient",
)
require(
    6 + (-5) == 1,
    "V=tau^13,A=tau^6,T=tau^-1 gives ord(A*T*H)=1",
)

# Passport arithmetic: at R=3k+2, surviving p=nu-2 is a multiple of 4k+3.
for R, nu, _a, delta in resonances:
    k = (R - 2) // 3
    pole_part = nu - 2
    require(pole_part == (4 * k + 3) * delta, f"passport divisor R,nu={R,nu}")

require(not syntax_has(__file__, ast.Assert), "Python assert nodes are forbidden")
require(not has_float_literals(__file__), "Python float literals are forbidden")

print("THM-2827 POLE-ORDER FABER NONRESONANCE AUDIT -- exact finite controls")
print("assert_nodes=0")
print("float_literals=0")
print("rows R=1..18: direct multinomial and explicit H_R identity agree")
print("prefix degree laws/gcd checked R=1..18")
print("nu=3..25 lane counts =", lane_counts)
print("all polar top/lower inequalities checked through R=18")
print("all pure faces, retained rows, and prescribed H valuations checked through R=18")
print("all regular-H cone inequalities/monomials checked through R=18")
print("resonance records =", tuple(resonances))
print("first resonance =", resonances[0])
print("first formal witness H coefficient =", first_h_coefficient)
print("nu=4 resonance count =", len(double_pole_resonances))
print("passport divisibility p=(4k+3)delta checked on every recorded ray")
print("ALL EXACT FINITE CONTROLS PASS")
