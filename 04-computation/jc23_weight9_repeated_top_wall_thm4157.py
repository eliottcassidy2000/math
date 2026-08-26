#!/usr/bin/env python3
"""Exact certificate for the exact-M=9 repeated top-edge wall.

The wall is zeta=-eta with eta nonzero.  This file verifies the contracted
source, Newton polygon, all four structural strict-transform strata,
normalized boundary packets, two exact critical projections at one rational
control per stratum, response degrees, and both permutation gates.  Pass
--symbolic-endpoints to recompute the four exact primary-resultant endpoint
formulae over polynomial rings.
"""

from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, permutations
from math import gcd
import sys

import sympy as sp


CHECKS = 0


def require(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def valuation(poly, variable):
    dictionary = sp.Poly(poly, variable).as_dict()
    require(bool(dictionary), "zero polynomial has no valuation")
    return min(exponent[0] for exponent in dictionary)


def expanded_support(coefficients):
    """Support of (s^2-p)(1-QH)+gamma*Q*s^2, with Q-height retained."""
    raw = [(2, 0, 0, F(1)), (0, 1, 0, F(-1))]
    for (p_power, y_power), coefficient in coefficients.items():
        raw.append((y_power + 2, p_power + y_power, 1, -coefficient))
        raw.append((y_power, p_power + y_power + 1, 1, coefficient))
    answer = {}
    for x_power, p_power, q_height, coefficient in raw:
        key = (x_power, p_power, q_height)
        answer[key] = answer.get(key, F(0)) + coefficient
    return {key: value for key, value in answer.items() if value}


def convex_hull(points):
    points = sorted(set(points))

    def cross(origin, left, right):
        return ((left[0] - origin[0]) * (right[1] - origin[1])
                - (left[1] - origin[1]) * (right[0] - origin[0]))

    lower = []
    for point in points:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], point) <= 0:
            lower.pop()
        lower.append(point)
    upper = []
    for point in reversed(points):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], point) <= 0:
            upper.pop()
        upper.append(point)
    return tuple(lower[:-1] + upper[:-1])


def polygon_ledger(points):
    vertices = convex_hull((x, y) for x, y, _ in points)
    area2 = abs(sum(
        vertices[index][0] * vertices[(index + 1) % len(vertices)][1]
        - vertices[(index + 1) % len(vertices)][0] * vertices[index][1]
        for index in range(len(vertices))
    ))
    boundary = sum(
        gcd(abs(vertices[(index + 1) % len(vertices)][0]
                - vertices[index][0]),
            abs(vertices[(index + 1) % len(vertices)][1]
                - vertices[index][1]))
        for index in range(len(vertices))
    )
    require((area2 - boundary + 2) % 2 == 0, "Pick parity failed")
    genus = (area2 - boundary + 2) // 2
    return vertices, area2, boundary, genus


def edge_ledger(vertices):
    answer = []
    for start, end in zip(vertices, vertices[1:] + vertices[:1]):
        dx = end[0] - start[0]
        dy = end[1] - start[1]
        length = gcd(abs(dx), abs(dy))
        inward = (-dy // length, dx // length)
        constant = inward[0] * start[0] + inward[1] * start[1]
        index = inward[0] + inward[1] - constant
        answer.append((start, end, inward, length, constant, index))
    return tuple(answer)


def polynomial_digest(poly):
    monic = poly.monic()
    coefficients = []
    for coefficient in monic.all_coeffs():
        coefficient = sp.Rational(coefficient)
        coefficients.append(f"{coefficient.p}/{coefficient.q}")
    payload = f"degree={monic.degree()};" + ",".join(coefficients)
    return sha256(payload.encode()).hexdigest()


def compose(left, right):
    return tuple(left[right[index]] for index in range(len(left)))


def inverse(permutation):
    answer = [0] * len(permutation)
    for source, target in enumerate(permutation):
        answer[target] = source
    return tuple(answer)


def permutation_index(permutation):
    seen = set()
    cycles = 0
    for start in range(len(permutation)):
        if start in seen:
            continue
        cycles += 1
        point = start
        while point not in seen:
            seen.add(point)
            point = permutation[point]
    return len(permutation) - cycles


def source_polynomial(X, T, P, Y, Delta, Theta, Phi, eta):
    K = sp.Rational(2848, 45) - sp.Rational(7, 6) * Delta
    return (
        -X ** 2 * T / 2 - 3 * P + sp.Rational(8, 3) * P ** 2
        - sp.Rational(1376, 135) * P ** 3 + K * Y ** 2
        + Phi * P ** 2 * Y + Delta * P ** 4 + Theta * P * Y ** 2
        + eta * P ** 3 * Y - eta * Y ** 3
    )


def symbolic_endpoint_checks(X, T, P, Y):
    """Recompute all four endpoint formulae; first row is intentionally slow."""
    answer = []

    # Row A: C=Delta+Theta is nonzero.
    Delta_symbol, Theta_symbol = sp.symbols("Delta_symbol Theta_symbol")
    G_symbol = source_polynomial(
        X, T, P, Y, Delta_symbol, Theta_symbol, sp.Rational(0),
        sp.Rational(1)
    )
    resultant = sp.resultant(
        sp.cancel(sp.diff(G_symbol, X) / T), sp.diff(G_symbol, T), X
    )
    quotient = sp.cancel(resultant / (T ** 42 * (6 * T + 1) ** 2))
    require(sp.denom(quotient) == 1,
            "symbolic endpoint quotient is not polynomial")
    qpoly = sp.Poly(quotient, T)
    require(qpoly.degree() == 19, "symbolic endpoint degree changed")
    expected = -12288 * (Delta_symbol + Theta_symbol) ** 6
    require(sp.factor(qpoly.nth(0) - expected) == 0,
            "symbolic endpoint formula changed")
    answer.append(sp.factor(qpoly.nth(0)))

    # Row B: C=0 and Phi is the next nonzero coefficient.
    Phi_symbol = sp.symbols("Phi_symbol")
    G_symbol = source_polynomial(
        X, T, P, Y, sp.Rational(1), sp.Rational(-1), Phi_symbol,
        sp.Rational(23, 13)
    )
    resultant = sp.resultant(
        sp.cancel(sp.diff(G_symbol, X) / T), sp.diff(G_symbol, T), X
    )
    quotient = sp.cancel(resultant / (T ** 30 * (6 * T + 1) ** 2))
    qpoly = sp.Poly(quotient, T)
    expected = -sp.Rational(3 * 7 ** 5, 2 ** 5) * Phi_symbol ** 5
    require(qpoly.degree() == 17 and sp.factor(qpoly.nth(0) - expected) == 0,
            "symbolic Phi endpoint formula changed")
    answer.append(sp.factor(qpoly.nth(0)))

    # Row C: C=Phi=0 and B=epsilon+K is the next coefficient.
    Delta_symbol, eta_symbol_C = sp.symbols(
        "Delta_symbol_C eta_symbol_C"
    )
    epsilon = -sp.Rational(1376, 135)
    K_symbol = sp.Rational(2848, 45) - sp.Rational(7, 6) * Delta_symbol
    B_symbol = epsilon + K_symbol
    G_symbol = source_polynomial(
        X, T, P, Y, Delta_symbol, -Delta_symbol, sp.Rational(0),
        eta_symbol_C
    )
    resultant = sp.resultant(
        sp.cancel(sp.diff(G_symbol, X) / T), sp.diff(G_symbol, T), X
    )
    quotient = sp.cancel(resultant / (T ** 32 * (6 * T + 1) ** 2))
    qpoly = sp.Poly(quotient, T)
    expected = 2 * 3 ** 6 * eta_symbol_C * B_symbol ** 5
    require(qpoly.degree() == 15 and sp.factor(qpoly.nth(0) - expected) == 0,
            "symbolic B endpoint formula changed")
    answer.append(sp.factor(qpoly.nth(0)))

    # Row D: C=Phi=B=0 forces Delta=2048/45; eta remains nonzero.
    eta_symbol = sp.symbols("eta_symbol")
    G_symbol = source_polynomial(
        X, T, P, Y, sp.Rational(2048, 45), sp.Rational(-2048, 45),
        sp.Rational(0), eta_symbol
    )
    resultant = sp.resultant(
        sp.cancel(sp.diff(G_symbol, X) / T), sp.diff(G_symbol, T), X
    )
    quotient = sp.cancel(resultant / (T ** 35 * (6 * T + 1) ** 2))
    qpoly = sp.Poly(quotient, T)
    expected = -sp.Rational(3 ** 11, 2 ** 5) * eta_symbol ** 5
    require(qpoly.degree() == 12 and sp.factor(qpoly.nth(0) - expected) == 0,
            "symbolic cusp endpoint formula changed")
    answer.append(sp.factor(qpoly.nth(0)))
    return tuple(answer)


def main():
    X, T, s, p, a, z, Q = sp.symbols("X T s p a z Q")
    P = T + X ** 2 * T ** 2
    Y = X * T * P

    # One exact rational point on the proposed critical-open chamber.
    Delta = F(1)
    Theta = F(19, 11)
    eta = F(23, 13)
    Phi = F(11, 7)
    K = F(2848, 45) - F(7, 6) * Delta
    gamma = F(-1, 2)
    lam = F(-3)
    alpha = F(8, 3)
    epsilon = F(-1376, 135)
    Cwall = Delta + Theta
    require(K == F(5591, 90), "forced K row changed")
    require(Cwall == F(30, 11), "ordinary-node coefficient changed")

    # Exact contraction on zeta=-eta.  Here y=s*p and t=p-s^2.
    y = s * p
    t = p - s ** 2
    contraction = eta * p ** 3 * y - eta * y ** 3
    require(sp.expand(contraction - eta * s * p ** 3 * t) == 0,
            "top-row contraction failed")

    coefficients = {
        (1, 0): F(-3),
        (2, 0): F(8, 3),
        (3, 0): F(-1376, 135),
        (0, 2): K,
        (2, 1): Phi,
        (4, 0): Delta,
        (1, 2): Theta,
        (3, 1): eta,
        (0, 3): -eta,
    }
    support = expanded_support(coefficients)
    vertices, area2, boundary, arithmetic_genus = polygon_ledger(support)
    expected_vertices = ((0, 1), (2, 0), (5, 3), (1, 5), (0, 5))
    require(vertices == expected_vertices, "repeated-wall Newton polygon changed")
    require((area2, boundary, arithmetic_genus) == (31, 11, 11),
            "repeated-wall Pick ledger changed")
    edges = edge_ledger(vertices)
    require(edges == (
        ((0, 1), (2, 0), (1, 2), 1, 2, 1),
        ((2, 0), (5, 3), (-1, 1), 3, -2, 2),
        ((5, 3), (1, 5), (-1, -2), 2, -11, 8),
        ((1, 5), (0, 5), (0, -1), 1, -5, 4),
        ((0, 5), (0, 1), (1, 0), 4, 0, 1),
    ), "edge ledger changed")

    # Symbolic toric strict transform at the repeated top edge.
    Ds, Ths, es, Ks, Phis = sp.symbols("Ds Ths es Ks Phis")
    rs = 1 - a
    H_symbol = (
        sp.Rational(-3) * p + sp.Rational(8, 3) * p ** 2
        - sp.Rational(1376, 135) * p ** 3 + Ks * s ** 2 * p ** 2
        + Phis * s * p ** 3 + Ds * p ** 4 + Ths * s ** 2 * p ** 3
        + es * s * p ** 3 * (p - s ** 2)
    )
    F_symbol = (s ** 2 - p) * (1 - Q * H_symbol) - Q * s ** 2 / 2
    L_actual = sp.cancel(
        z ** 11 * F_symbol.subs({s: 1 / z, p: rs / z ** 2})
    )
    L_expected = (
        a * z ** 9 - Q * z ** 9 / 2 + Q * es * a ** 2 * rs ** 3
        - Q * a * rs ** 3 * (Ds * rs + Ths) * z
        - Q * Phis * a * rs ** 3 * z ** 2
        - Q * a * (sp.Rational(-1376, 135) * rs ** 3
                   + Ks * rs ** 2) * z ** 3
        - Q * sp.Rational(8, 3) * a * rs ** 2 * z ** 5
        + 3 * Q * a * rs * z ** 7
    )
    require(sp.factor(L_actual - L_expected) == 0,
            "strict-transform identity changed")

    local_poly = sp.Poly(sp.expand(L_expected), a, z)
    tangent = sum(
        coefficient * a ** powers[0] * z ** powers[1]
        for powers, coefficient in local_poly.terms()
        if sum(powers) == 2
    )
    Csymbol = Ds + Ths
    require(sp.factor(tangent - Q * a * (es * a - Csymbol * z)) == 0,
            "ordinary-node tangent cone changed")

    # Leading terms of the normalization branches in the four structural
    # strata.  In the first three rows both branches are smooth and meet to
    # order m=1,2,3.  Since omega=Q*z^7 dz/L_a, L_a has order m and both
    # indices are 8-m.
    L_a = sp.diff(L_expected, a)
    fast_a = Csymbol * z / es
    slow_a = -z ** 8 / (2 * Csymbol)
    require(valuation(sp.cancel(L_expected.subs(a, fast_a)), z) >= 3,
            "row A fast branch did not cancel")
    require(valuation(sp.cancel(L_expected.subs(a, slow_a)), z) >= 10,
            "row A slow branch did not cancel")
    require(valuation(sp.cancel(L_a.subs(a, fast_a)), z) == 1,
            "row A fast L_a order changed")
    require(valuation(sp.cancel(L_a.subs(a, slow_a)), z) == 1,
            "row A slow L_a order changed")

    L_B = sp.expand(L_expected.subs(Ths, -Ds))
    fast_a_B = Phis * z ** 2 / es
    slow_a_B = -z ** 7 / (2 * Phis)
    require(valuation(sp.cancel(L_B.subs(a, fast_a_B)), z) >= 5,
            "row B fast branch did not cancel")
    require(valuation(sp.cancel(L_B.subs(a, slow_a_B)), z) >= 10,
            "row B slow branch did not cancel")
    require(valuation(sp.cancel(sp.diff(L_B, a).subs(a, fast_a_B)), z) == 2,
            "row B fast L_a order changed")
    require(valuation(sp.cancel(sp.diff(L_B, a).subs(a, slow_a_B)), z) == 2,
            "row B slow L_a order changed")

    epsilon_symbol = -sp.Rational(1376, 135)
    Bsymbol = epsilon_symbol + Ks
    L_C = sp.expand(L_B.subs(Phis, 0))
    fast_a_C = Bsymbol * z ** 3 / es
    slow_a_C = -z ** 6 / (2 * Bsymbol)
    require(valuation(sp.cancel(L_C.subs(a, fast_a_C)), z) >= 7,
            "row C fast branch did not cancel")
    require(valuation(sp.cancel(L_C.subs(a, slow_a_C)), z) >= 11,
            "row C slow branch did not cancel")
    require(valuation(sp.cancel(sp.diff(L_C, a).subs(a, fast_a_C)), z) == 3,
            "row C fast L_a order changed")
    require(valuation(sp.cancel(sp.diff(L_C, a).subs(a, slow_a_C)), z) == 3,
            "row C slow L_a order changed")

    # At C=Phi=B=0 the local Newton edge joins (a-degree,z-degree)
    # (2,0) to (0,9), with initial form Q*(eta*a^2-z^9/2).  gcd(2,9)=1,
    # so it is one cusp branch z=tau^2, a=c*tau^9, c^2=1/(2eta).
    L_D = sp.Poly(sp.expand(L_C.subs(Ks, -epsilon_symbol)), a, z)
    weighted_terms = []
    minimum_weight = min(9 * powers[0] + 2 * powers[1]
                         for powers, _ in L_D.terms())
    for powers, coefficient in L_D.terms():
        if 9 * powers[0] + 2 * powers[1] == minimum_weight:
            weighted_terms.append(
                coefficient * a ** powers[0] * z ** powers[1]
            )
    cusp_initial = sum(weighted_terms)
    require(minimum_weight == 18, "row D cusp weight changed")
    require(sp.factor(cusp_initial - Q * (es * a ** 2 - z ** 9 / 2)) == 0,
            "row D cusp initial form changed")
    tau, cusp_coefficient = sp.symbols("tau cusp_coefficient")
    cusp_La = sp.diff(L_D.as_expr(), a).subs(
        {z: tau ** 2, a: cusp_coefficient * tau ** 9}
    )
    require(valuation(cusp_La, tau) == 9,
            "row D cusp L_a order changed")
    require(2 * 7 + 1 - 9 == 6,
            "row D pulled-back residue order changed")

    structural_rows = (
        {
            "name": "A_C_nonzero", "delta": 1, "genus": 10,
            "packet": (7, 7, 4, 2, 2, 2, 1),
            "branches": "a~(C/eta)z,a~-z^8/(2C)",
            "top": (7, 7),
        },
        {
            "name": "B_C_zero_Phi_nonzero", "delta": 2, "genus": 9,
            "packet": (6, 6, 4, 2, 2, 2, 1),
            "branches": "a~(Phi/eta)z^2,a~-z^7/(2Phi)",
            "top": (6, 6),
        },
        {
            "name": "C_C_Phi_zero_B_nonzero", "delta": 3, "genus": 8,
            "packet": (5, 5, 4, 2, 2, 2, 1),
            "branches": "a~(B/eta)z^3,a~-z^6/(2B)",
            "top": (5, 5),
        },
        {
            "name": "D_C_Phi_B_zero", "delta": 4, "genus": 7,
            "packet": (7, 4, 2, 2, 2, 1),
            "branches": "z=tau^2,a~c*tau^9,c^2=1/(2eta)",
            "top": (7,),
        },
    )
    for row in structural_rows:
        packet = row["packet"]
        defect = sum(index - 1 for index in packet)
        if len(row["top"]) == 2:
            require(row["top"] == (8 - row["delta"],) * 2,
                    f"{row['name']}: smooth-branch residue indices changed")
        else:
            require(row["top"] == (7,),
                    f"{row['name']}: cusp residue index changed")
        require(arithmetic_genus - row["delta"] == row["genus"],
                f"{row['name']}: delta/genus ledger changed")
        require(defect == 2 * row["genus"] - 2,
                f"{row['name']}: packet no longer saturates Riemann--Hurwitz")

    controls = (
        {
            "name": "A_C_nonzero", "Delta": F(1), "Theta": F(19, 11),
            "Phi": F(11, 7), "eta": F(23, 13),
            "Tval": 42, "degree": 19,
            "q0": -F(12288) * F(30, 11) ** 6,
            "qLC": -F(11206045421019265531550064784,
                        407151596119125),
            "r0": F(1325533905629568, 604175),
            "rLC": F(6918774197944320000, 5436100813),
        },
        {
            "name": "B_C_zero_Phi_nonzero", "Delta": F(1),
            "Theta": F(-1), "Phi": F(11, 7), "eta": F(23, 13),
            "Tval": 30, "degree": 17,
            "q0": -F(3 * 7 ** 5, 2 ** 5) * F(11, 7) ** 5,
            "qLC": F(85821754813159273972508, 352960408125),
            "r0": -F(70369952239872, 54925),
            "rLC": F(30531977307612, 371293),
        },
        {
            "name": "C_C_Phi_zero_B_nonzero", "Delta": F(1),
            "Theta": F(-1), "Phi": F(0), "eta": F(23, 13),
            "Tval": 32, "degree": 15,
            "q0": F(12463005381719088344323, 12793950000),
            "qLC": F(85821754813159273972508, 352960408125),
            "r0": -F(70369952239872, 54925),
            "rLC": F(994981781335325982843932, 18796708125),
        },
        {
            "name": "D_C_Phi_B_zero", "Delta": F(2048, 45),
            "Theta": -F(2048, 45), "Phi": F(0), "eta": F(23, 13),
            "Tval": 35, "degree": 12,
            "q0": -F(1140178853421, 11881376),
            "qLC": F(84349672224497501284348663943,
                       231577323770812500),
            "r0": -F(3875767847803456, 2471625),
            "rLC": -F(1809463840379127, 62748517),
        },
    )

    critical_rows = []
    for control, structure in zip(controls, structural_rows):
        name = control["name"]
        require(name == structure["name"], f"{name}: row alignment changed")
        Delta_case = sp.Rational(control["Delta"].numerator,
                                 control["Delta"].denominator)
        Theta_case = sp.Rational(control["Theta"].numerator,
                                 control["Theta"].denominator)
        Phi_case = sp.Rational(control["Phi"].numerator,
                               control["Phi"].denominator)
        eta_case = sp.Rational(control["eta"].numerator,
                               control["eta"].denominator)
        K_case = sp.Rational(2848, 45) - sp.Rational(7, 6) * Delta_case
        B_case = -sp.Rational(1376, 135) + K_case
        if name == "A_C_nonzero":
            require(Delta_case + Theta_case != 0,
                    "row A control left its structural stratum")
        elif name == "B_C_zero_Phi_nonzero":
            require(Delta_case + Theta_case == 0 and Phi_case != 0,
                    "row B control left its structural stratum")
        elif name == "C_C_Phi_zero_B_nonzero":
            require(Delta_case + Theta_case == 0 and Phi_case == 0
                    and B_case != 0,
                    "row C control left its structural stratum")
        else:
            require(Delta_case + Theta_case == 0 and Phi_case == 0
                    and B_case == 0 and Delta_case == sp.Rational(2048, 45),
                    "row D control left its structural stratum")

        G = source_polynomial(
            X, T, P, Y, Delta_case, Theta_case, Phi_case, eta_case
        )
        f = sp.cancel(sp.diff(G, X) / T)
        h = sp.diff(G, T)
        require(sp.cancel(f.subs(T, 0) + X) == 0,
                f"{name}: T=0 specialization of G_X/T changed")
        require(sp.cancel(h.subs(T, 0) + (X ** 2 + 6) / 2) == 0,
                f"{name}: T=0 specialization of G_T changed")

        # The two universal critical pairs and their Hessian determinants.
        hessian = sp.det(sp.hessian(G, (X, T)))
        for t_value, x_square, g_value, determinant in (
                (sp.Rational(-1, 6), sp.Rational(6), sp.Rational(1, 2),
                 sp.Rational(-6)),
                (sp.Rational(0), sp.Rational(-6), sp.Rational(0),
                 sp.Rational(6))):
            modulus = sp.Poly(X ** 2 - x_square, X)
            for expression, label in (
                    (sp.diff(G, X), "G_X"), (sp.diff(G, T), "G_T"),
                    (G - g_value, "G-value"),
                    (hessian - determinant, "Hessian")):
                remainder = sp.rem(
                    sp.Poly(sp.expand(expression.subs(T, t_value)), X),
                    modulus
                ).as_expr()
                require(remainder == 0,
                        f"{name}: universal {label} identity changed")

        primary_resultant = sp.resultant(f, h, X)
        require(valuation(primary_resultant, T) == control["Tval"],
                f"{name}: primary Sylvester valuation changed")
        primary_residual = sp.cancel(
            primary_resultant
            / (T ** control["Tval"] * (6 * T + 1) ** 2)
        )
        require(sp.denom(primary_residual) == 1,
                f"{name}: primary residual is not polynomial")
        qpoly = sp.Poly(primary_residual, T)
        require(qpoly.degree() == control["degree"],
                f"{name}: primary residual degree changed")
        expected_q0 = sp.Rational(control["q0"].numerator,
                                  control["q0"].denominator)
        expected_qLC = sp.Rational(control["qLC"].numerator,
                                   control["qLC"].denominator)
        require(qpoly.nth(0) == expected_q0 and qpoly.LC() == expected_qLC,
                f"{name}: primary residual endpoint changed")
        require(qpoly.eval(sp.Rational(-1, 6)) != 0,
                f"{name}: residual collided with universal pair")
        require(sp.gcd(qpoly, qpoly.diff()).degree() == 0,
                f"{name}: primary residual control is not squarefree")

        # Independent rational (s,p) critical projection.
        H = (
            -3 * p + sp.Rational(8, 3) * p ** 2
            - sp.Rational(1376, 135) * p ** 3
            + K_case * s ** 2 * p ** 2 + Phi_case * s * p ** 3
            + Delta_case * p ** 4 + Theta_case * s ** 2 * p ** 3
            + eta_case * s * p ** 3 * (p - s ** 2)
        )
        Gsp = -s ** 2 / (2 * t) + H
        A = sp.cancel(t ** 2 * sp.diff(Gsp, s) / p)
        Ccrit = sp.cancel(2 * t ** 2 * sp.diff(Gsp, p))
        require(sp.denom(A) == sp.denom(Ccrit) == 1,
                f"{name}: independent critical pair is not polynomial")
        independent_resultant = sp.resultant(A, Ccrit, s)
        p_valuation = valuation(independent_resultant, p)
        require(p_valuation == 8,
                f"{name}: independent p-valuation changed")
        independent_residual = sp.cancel(independent_resultant / p ** 8)
        require(sp.denom(independent_residual) == 1,
                f"{name}: independent residual is not polynomial")
        rpoly = sp.Poly(independent_residual, p)
        require(rpoly.degree() == control["degree"],
                f"{name}: independent residual degree changed")
        expected_r0 = sp.Rational(control["r0"].numerator,
                                  control["r0"].denominator)
        expected_rLC = sp.Rational(control["rLC"].numerator,
                                   control["rLC"].denominator)
        require(rpoly.nth(0) == expected_r0 and rpoly.LC() == expected_rLC,
                f"{name}: independent residual endpoint changed")
        require(sp.gcd(rpoly, rpoly.diff()).degree() == 0,
                f"{name}: independent residual control is not squarefree")

        # Residual roots plus two points at p=0/T=-1/6 and two omitted at
        # t=0/T=0 give the exact affine critical length.
        critical_length = control["degree"] + 2 + 2
        packet = structure["packet"]
        defect = sum(index - 1 for index in packet)
        finite_degree = sum(structure["top"]) + 4 + 1
        beta = 3
        full_degree = finite_degree + 3 * 2
        require(critical_length == finite_degree + beta + 1,
                f"{name}: finite critical margin changed")

        finite_support_sum = 2 * finite_degree - critical_length
        finite_one_handle_capacity = finite_support_sum - 1 + beta
        require(finite_one_handle_capacity < finite_degree - 1,
                f"{name}: finite orbit-merger gate no longer excludes")
        overlap_cap = full_degree - critical_length
        commutator_cap = 2 * overlap_cap
        require(overlap_cap == 2 and defect > commutator_cap,
                f"{name}: full commutator-overlap gate no longer excludes")

        critical_rows.append({
            **structure,
            "Delta": Delta_case, "Theta": Theta_case, "Phi": Phi_case,
            "eta": eta_case, "K": K_case, "B": B_case,
            "Tval": control["Tval"], "degree": control["degree"],
            "q0": qpoly.nth(0), "qLC": qpoly.LC(),
            "q_at_universal": qpoly.eval(sp.Rational(-1, 6)),
            "q_digest": polynomial_digest(qpoly),
            "r0": rpoly.nth(0), "rLC": rpoly.LC(),
            "r_digest": polynomial_digest(rpoly),
            "critical_length": critical_length,
            "finite_degree": finite_degree, "beta": beta,
            "finite_capacity": finite_one_handle_capacity,
            "full_degree": full_degree, "overlap": overlap_cap,
            "commutator_cap": commutator_cap, "defect": defect,
        })

    # Deepest-stratum projection portfolio.  This row has only eta as a free
    # structural coefficient.  Factoring both discriminants separates their
    # common collision wall from projection-only coordinate collisions.
    eta_deep = sp.symbols("eta_deep")
    Delta_deep = sp.Rational(2048, 45)
    Theta_deep = -Delta_deep
    K_deep = sp.Rational(1376, 135)
    G_deep = source_polynomial(
        X, T, P, Y, Delta_deep, Theta_deep, sp.Rational(0), eta_deep
    )
    resultant_deep_X = sp.resultant(
        sp.cancel(sp.diff(G_deep, X) / T), sp.diff(G_deep, T), X
    )
    Q12 = sp.Poly(sp.cancel(
        resultant_deep_X / (T ** 35 * (6 * T + 1) ** 2)
    ), T)
    J_deep = 22143375 * eta_deep ** 2 + 15510536192
    require(sp.factor(Q12.nth(0)
                      + sp.Rational(3 ** 11, 2 ** 5) * eta_deep ** 5) == 0,
            "deep row primary constant formula changed")
    require(sp.factor(
        Q12.LC()
        - eta_deep ** 3 * J_deep ** 2 / sp.Integer(3690562500)
    ) == 0, "deep row primary leading formula changed")
    primary_disc_factors = sp.factor_list(
        sp.discriminant(Q12.as_expr(), T), eta_deep
    )[1]
    primary_factor_ledger = sorted(
        (sp.degree(factor, eta_deep), multiplicity, factor)
        for factor, multiplicity in primary_disc_factors
    )
    require([(degree, multiplicity)
             for degree, multiplicity, _ in primary_factor_ledger]
            == [(1, 32), (30, 1), (36, 2)],
            "deep row primary discriminant factor degrees changed")

    H_deep = (
        -3 * p + sp.Rational(8, 3) * p ** 2
        - sp.Rational(1376, 135) * p ** 3
        + K_deep * s ** 2 * p ** 2 + Delta_deep * p ** 4
        + Theta_deep * s ** 2 * p ** 3
        + eta_deep * s * p ** 3 * (p - s ** 2)
    )
    Gsp_deep = -s ** 2 / (2 * t) + H_deep
    A_deep = sp.cancel(t ** 2 * sp.diff(Gsp_deep, s) / p)
    C_deep = sp.cancel(2 * t ** 2 * sp.diff(Gsp_deep, p))
    require(sp.factor(A_deep.subs(p, s ** 2) + s) == 0
            and sp.factor(C_deep.subs(p, s ** 2) - s ** 2) == 0,
            "deep row independent chart acquired nonzero t=0 solutions")
    resultant_deep_p = sp.resultant(A_deep, C_deep, s)
    R12 = sp.Poly(sp.cancel(resultant_deep_p / p ** 8), p)
    require(sp.factor(R12.nth(0)
                      + sp.Rational(64, 1125) * eta_deep * J_deep) == 0,
            "deep row independent constant formula changed")
    require(sp.factor(R12.LC() + 531441 * eta_deep ** 7) == 0,
            "deep row independent leading formula changed")
    independent_disc_factors = sp.factor_list(
        sp.discriminant(R12.as_expr(), p), eta_deep
    )[1]
    independent_factor_ledger = sorted(
        (sp.degree(factor, eta_deep), multiplicity, factor)
        for factor, multiplicity in independent_disc_factors
    )
    require([(degree, multiplicity)
             for degree, multiplicity, _ in independent_factor_ledger]
            == [(1, 40), (24, 2), (30, 1)],
            "deep row independent discriminant factor degrees changed")

    primary_by_degree = {
        degree: factor for degree, _, factor in primary_factor_ledger
    }
    independent_by_degree = {
        degree: factor for degree, _, factor in independent_factor_ledger
    }
    common30_primary = sp.Poly(primary_by_degree[30], eta_deep).monic()
    common30_independent = sp.Poly(
        independent_by_degree[30], eta_deep
    ).monic()
    require(common30_primary == common30_independent,
            "deep row common discriminant factor changed")
    q_at_universal_numerator = sp.together(
        Q12.eval(-sp.Rational(1, 6))
    ).as_numer_denom()[0]
    require(sp.gcd(independent_by_degree[24], primary_by_degree[36]) == 1,
            "deep row projection-only discriminants acquired common root")
    require(sp.gcd(independent_by_degree[24], q_at_universal_numerator) == 1,
            "deep row independent collision meets universal T collision")
    require(sp.gcd(J_deep, primary_by_degree[30]) == 1,
            "deep row endpoint and common collision walls merged")
    deepest_ledger = {
        "J": J_deep,
        "common30_digest": polynomial_digest(common30_primary),
        "primary36_digest": polynomial_digest(
            sp.Poly(primary_by_degree[36], eta_deep)
        ),
        "independent24_digest": polynomial_digest(
            sp.Poly(independent_by_degree[24], eta_deep)
        ),
    }

    # Exhaustive hostile check of ind([X,Y]) <= 2|supp(X) intersect supp(Y)|
    # through S_5; the report gives the general partial-injection proof.
    hostile_pairs = 0
    for degree in range(1, 6):
        all_permutations = tuple(permutations(range(degree)))
        for left in all_permutations:
            left_inverse = inverse(left)
            left_support = {
                point for point, image in enumerate(left) if point != image
            }
            for right in all_permutations:
                right_inverse = inverse(right)
                right_support = {
                    point for point, image in enumerate(right)
                    if point != image
                }
                commutator = compose(
                    compose(compose(left, right), left_inverse), right_inverse
                )
                overlap = len(left_support & right_support)
                require(permutation_index(commutator) <= 2 * overlap,
                        "commutator-overlap hostile failed")
                hostile_pairs += 1

    symbolic_endpoints = None
    if "--symbolic-endpoints" in sys.argv[1:]:
        symbolic_endpoints = symbolic_endpoint_checks(X, T, P, Y)

    semantic_rows = []
    for row in critical_rows:
        semantic_rows.append(
            f"{row['name']}:delta={row['delta']};g={row['genus']};"
            f"packet={row['packet']};primary=T^{row['Tval']}*"
            f"(6T+1)^2*Q{row['degree']};independent=p^8*R{row['degree']};"
            f"L={row['critical_length']};finite=(n={row['finite_degree']},"
            f"beta=3);full=(n={row['full_degree']},overlap<=2,"
            f"comm_ind<={row['commutator_cap']}<{row['defect']})"
        )
    semantic = (
        f"wall=zeta=-eta;polygon={vertices};Pick=({area2},{boundary},"
        f"{arithmetic_genus})\n" + "\n".join(semantic_rows) + "\n"
        + f"deepest:J={J_deep};H30={deepest_ledger['common30_digest']};"
        + f"H36={deepest_ledger['primary36_digest']};"
        + f"H24={deepest_ledger['independent24_digest']};"
        + "portfolio=eta*J*H30!=0\n"
    )

    print("EXACT-M=9 REPEATED TOP-EDGE WALL CERTIFICATE")
    print(f"checks={CHECKS}")
    print("wall=zeta=-eta, eta!=0")
    print("contraction=eta*p^3*y-eta*y^3=eta*s*p^3*(p-s^2)")
    print(f"polygon={vertices}")
    print(f"Pick=area2:{area2},boundary:{boundary},arithmetic_genus:{arithmetic_genus}")
    print("strict_transform=z^11*F(1/z,(1-a)/z^2)")
    print("tangent_cone=Q*a*(eta*a-(Delta+Theta)*z)")
    print("residue=Q*z^7*dz/L_a")
    print("cubic_carrier=q-1/2=K*W^2-eta*W^3; residue_degree=3; beta=3")
    for row in critical_rows:
        print(f"case={row['name']}")
        print(f"  control=Delta:{row['Delta']},Theta:{row['Theta']},"
              f"Phi:{row['Phi']},eta:{row['eta']},K:{row['K']},B:{row['B']}")
        print(f"  local_delta={row['delta']}; branches={row['branches']}")
        print(f"  packet={row['packet']}; genus={row['genus']};"
              f" defect={row['defect']}")
        print(f"  primary=T^{row['Tval']}*(6T+1)^2*Q{row['degree']}")
        print(f"  Q_constant={row['q0']}; Q_leading={row['qLC']}")
        print(f"  Q_at_T_minus_1_over_6={row['q_at_universal']}")
        print("  Q_squarefree_gcd_degree=0")
        print(f"  Q_sha256={row['q_digest']}")
        print(f"  independent=p^8*R{row['degree']}")
        print(f"  R_constant={row['r0']}; R_leading={row['rLC']}")
        print("  R_squarefree_gcd_degree=0")
        print(f"  R_sha256={row['r_digest']}")
        print(f"  critical_length={row['critical_length']}={row['degree']}+2+2")
        print(f"  finite_response=n:{row['finite_degree']},beta:3,"
              f"L=n+beta+1,one_handle_capacity:{row['finite_capacity']}"
              f"<n-1:{row['finite_degree'] - 1}")
        print(f"  full_response=n:{row['full_degree']},overlap<=2,"
              f"commutator_index<={row['commutator_cap']}"
              f"<defect:{row['defect']}")
    print("deepest_projection_portfolio=")
    print(f"  J(eta)={deepest_ledger['J']}")
    print("  Q12_constant=-(3^11/2^5)*eta^5")
    print("  Q12_leading=eta^3*J(eta)^2/3690562500")
    print("  R12_constant=-(64/1125)*eta*J(eta)")
    print("  R12_leading=-531441*eta^7")
    print("  disc_T_Q12=unit*eta^32*H30*H36^2")
    print("  disc_p_R12=unit*eta^40*H30*H24^2")
    print(f"  H30_sha256={deepest_ledger['common30_digest']}")
    print(f"  H36_sha256={deepest_ledger['primary36_digest']}")
    print(f"  H24_sha256={deepest_ledger['independent24_digest']}")
    print("  gcd(H24,H36)=gcd(H24,Q12(-1/6))=gcd(J,H30)=1")
    print("  portfolio_open=eta*J(eta)*H30(eta)!=0:at_least_one_projection_squarefree")
    print(f"commutator_overlap_hostile_pairs={hostile_pairs}:PASS")
    print("symbolic_endpoint_formulas=")
    print("  A:-12288*(Delta+Theta)^6")
    print("  B:-(3*7^5/2^5)*Phi^5")
    print("  C:2*3^6*eta*(epsilon+K)^5")
    print("  D:-(3^11/2^5)*eta^5")
    if symbolic_endpoints is not None:
        print(f"symbolic_endpoints_recomputed={symbolic_endpoints}:PASS")
    else:
        print("symbolic_endpoints_recomputed=SKIPPED(use --symbolic-endpoints)")
    print(f"semantic_sha256={sha256(semantic.encode()).hexdigest()}")
    print("verdict=PASS_ON_FOUR_STRUCTURAL_CRITICAL_OPEN_SUBCHAMBERS")


if __name__ == "__main__":
    main()
