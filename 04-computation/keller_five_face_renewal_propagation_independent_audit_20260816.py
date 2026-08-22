#!/usr/bin/env python3
"""Independent hostile audit of five-face renewal propagation.

This implementation does not import the candidate companion.  It starts from
the inverse-chart numerators frozen in THM-3495, extracts both THM-3513 hybrid
initial systems, checks the nonmonic cubic norm with a multiplication matrix
and a resultant, and then audits the abstract A(e,m) face algebra.

Scope: a full A(e,m) packet propagates to Q=L^e Norm(P) only after Q is known
to be polynomial.  This script does not infer that polynomiality, an image
equation, irreducibility, a fibre degree, or an all-level Keller law.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass

import sympy as sp


EXPECTED_SEMANTIC_SHA256 = "9de2b0a149105263ee1b3a1fba01424f9c7ff274368c689cdc5737fb340bf804"


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def require_zero(expression: sp.Expr, label: str) -> None:
    require(sp.cancel(expression) == 0, (label, sp.factor(expression)))


def weighted_initial(
    expression: sp.Expr,
    variables: tuple[sp.Symbol, ...],
    weights: tuple[int, ...],
    *,
    take_minimum: bool,
) -> tuple[int, sp.Expr]:
    polynomial = sp.Poly(sp.expand(expression), *variables, domain=sp.QQ)
    rows = [
        (sum(int(power) * weight for power, weight in zip(monomial, weights)), monomial, coefficient)
        for monomial, coefficient in polynomial.terms()
    ]
    extreme = min(row[0] for row in rows) if take_minimum else max(row[0] for row in rows)
    face = sum(
        coefficient
        * sp.prod(variable**int(power) for variable, power in zip(variables, monomial))
        for row_weight, monomial, coefficient in rows
        if row_weight == extreme
    )
    return extreme, sp.expand(face)


def zero_mod_cubic(expression: sp.Expr, cubic: sp.Expr, q: sp.Symbol) -> bool:
    numerator = sp.together(expression).as_numer_denom()[0]
    other_symbols = tuple(sorted(numerator.free_symbols.union(cubic.free_symbols) - {q}, key=str))
    field = sp.QQ.frac_field(*other_symbols)
    remainder = sp.rem(sp.Poly(numerator, q, domain=field), sp.Poly(cubic, q, domain=field))
    return remainder.is_zero


@dataclass(frozen=True)
class Packet:
    e: int
    m: int
    a: int
    b: int
    h: int
    rho: int
    gamma_min: int
    beta_min: int


def packet(e: int, m: int) -> Packet:
    require(e >= 0 and m >= 0 and m % 3 == 0, ("packet-domain", e, m))
    a = 2 * e - 4 * m // 3
    b = 2 * e - 2 * m // 3
    h = e - 2 * m // 3
    require(min(a, b, h, 3 * e - 2 * m, e - m) >= 0, ("packet-exponents", e, m))
    require(a == 2 * h and b == e + h and b % 2 == 0, ("common-monomial", e, m))
    return Packet(
        e=e,
        m=m,
        a=a,
        b=b,
        h=h,
        rho=a - 3 * b,
        gamma_min=-8 * e + 2 * m,
        beta_min=-5 * e + 2 * m,
    )


def transform(source: Packet) -> Packet:
    return packet(7 * source.e - 2 * source.m, 3 * source.e - 2 * source.m)


def scalar_step(source: Packet, top_coefficient: tuple[int, int]) -> dict[str, tuple[int, int]]:
    """Prime exponents for coefficients, represented as 2^u 3^v."""

    u, v = top_coefficient
    gamma = (3 * u + source.rho, 3 * v)
    top = (3 * u + source.rho, 3 * v + 3 * source.e - 3 * source.rho)
    return {"top": top, "gamma": gamma}


def face_ledger(value: Packet) -> dict[str, object]:
    return {
        "state": (value.e, value.m),
        "max_lambda": (value.e, value.m),
        "min_lambda": (3 * value.e - 2 * value.m, value.e),
        "min_beta": (
            value.beta_min,
            3 * value.e - 2 * value.m,
            value.e - value.m,
            2 * value.m // 3,
            value.m // 3,
        ),
        "max_z": (value.a, value.b),
        "min_gamma": (value.gamma_min, value.e, value.h),
        "hybrid_minima": (value.gamma_min - value.b, value.gamma_min - 3 * value.b),
    }


def main() -> None:
    # ------------------------------------------------------------------
    # Route I: reconstruct both inverse hybrid limits from THM-3495's
    # literal inverse-chart numerators, independently of the candidate.
    # ------------------------------------------------------------------
    A, B, C, q = sp.symbols("A B C q")
    variables = (A, B, C)

    L = 27 * A**2 * C**2 - 18 * A * B * C + 16 * A + B**3 * C - B**2
    T = 4 - 3 * B * C
    S = 27 * A * C**2 - 9 * B * C + 8
    K = 9 * A * C - B
    Y0 = 81 * A * B * C**2 - 72 * A * C - 15 * B**2 * C + 16 * B
    A1 = 27 * A * B * C**2 + 54 * A * C - 9 * B**2 * C + 2 * B
    A2 = 27 * A * B**2 * C**2 + 18 * A * B * C - 48 * A - 9 * B**3 * C + 10 * B**2
    Z0 = (
        -2916 * A**3 * C**4
        + 2916 * A**2 * B * C**3
        - 4536 * A**2 * C**2
        + 621 * A * B**3 * C**3
        - 1026 * A * B**2 * C**2
        + 504 * A * B * C
        + 64 * A
        - 207 * B**4 * C**2
        + 454 * B**3 * C
        - 256 * B**2
    )
    invariants = {"L": L, "T": T, "S": S, "K": K, "Y0": Y0, "A1": A1, "A2": A2, "Z0": Z0}

    top = {name: weighted_initial(expr, variables, (0, 0, 3), take_minimum=False) for name, expr in invariants.items()}
    top_expected = {
        "L": (6, 27 * A**2 * C**2),
        "T": (3, -3 * B * C),
        "S": (6, 27 * A * C**2),
        "K": (3, 9 * A * C),
        "Y0": (6, 81 * A * B * C**2),
        "A1": (6, 27 * A * B * C**2),
        "A2": (6, 27 * A * B**2 * C**2),
        "Z0": (12, -2916 * A**3 * C**4),
    }
    for name, (expected_weight, expected_face) in top_expected.items():
        actual_weight, actual_face = top[name]
        require(actual_weight == expected_weight, ("top-initial-weight", name, actual_weight))
        require_zero(actual_face - expected_face, f"top-initial-face-{name}")
    top_cubic = 27 * A**2 * C * q**3 - 2
    top_y = -3 * top["K"][1] * top["L"][1] * q**2 / (2 * top["S"][1])
    top_z = top["Z0"][1] / (8 * top["S"][1])
    require(zero_mod_cubic(q * top_y + 1, top_cubic, q), "top-y=-1/q")
    require(zero_mod_cubic(q**3 * top_z + C, top_cubic, q), "top-z=-C/q^3")
    require_zero(27 * A**2 * C**2 * q**3 - 2 * C - C * top_cubic, "top-residual-cubic")

    D = 27 * A**2 * C + B**3
    gamma = {name: weighted_initial(expr, variables, (1, -1, -5), take_minimum=True) for name, expr in invariants.items()}
    gamma_expected = {
        "L": (-8, C * D),
        "T": (-6, -3 * B * C),
        "S": (-9, 27 * A * C**2),
        "K": (-4, 9 * A * C),
        "Y0": (-10, 81 * A * B * C**2),
        "A1": (-10, 27 * A * B * C**2),
        "A2": (-11, 27 * A * B**2 * C**2),
        "Z0": (-17, -2916 * A**3 * C**4 + 621 * A * B**3 * C**3),
    }
    for name, (expected_weight, expected_face) in gamma_expected.items():
        actual_weight, actual_face = gamma[name]
        require(actual_weight == expected_weight, ("gamma-initial-weight", name, actual_weight))
        require_zero(actual_face - expected_face, f"gamma-initial-face-{name}")
    gamma_cubic = D * q**3 - 3 * B * q - 2
    gamma_y = (gamma["Y0"][1] - 3 * gamma["K"][1] * gamma["L"][1] * q**2) / (2 * gamma["S"][1])
    gamma_z = (
        gamma["Z0"][1]
        + 6 * gamma["L"][1] * gamma["A1"][1] * q
        - 9 * gamma["L"][1] * gamma["A2"][1] * q**2
    ) / (8 * gamma["S"][1])
    require(zero_mod_cubic(q * gamma_y + 1, gamma_cubic, q), "gamma-y=-1/q")
    require(zero_mod_cubic(q**3 * gamma_z + C, gamma_cubic, q), "gamma-z=-C/q^3")
    require_zero(C * D * q**3 - 3 * B * C * q - 2 * C - C * gamma_cubic, "gamma-residual-cubic")

    # Two independent nonmonic checks.  The determinant is the algebra norm
    # of q in the monic quotient; the resultant still contains the leading D.
    multiplication_by_q = sp.Matrix([[0, 0, sp.Rational(2, 1) / D], [1, 0, 3 * B / D], [0, 1, 0]])
    require_zero(multiplication_by_q.det() - 2 / D, "nonmonic-companion-determinant")
    require_zero(sp.resultant(gamma_cubic, q, q) - 2, "nonmonic-resultant")
    top_multiplication_by_q = sp.Matrix(
        [[0, 0, sp.Rational(2, 1) / (27 * A**2 * C)], [1, 0, 0], [0, 1, 0]]
    )
    require_zero(top_multiplication_by_q.det() - 2 / (27 * A**2 * C), "top-companion-determinant")

    # ------------------------------------------------------------------
    # Route II: symbolic exponent identities and hostile equality tests.
    # ------------------------------------------------------------------
    e, m = sp.symbols("e m")
    a = 2 * e - sp.Rational(4, 3) * m
    b = 2 * e - sp.Rational(2, 3) * m
    h = e - sp.Rational(2, 3) * m
    rho = a - 3 * b
    e_next = 7 * e - 2 * m
    m_next = 3 * e - 2 * m
    a_next = 2 * e_next - sp.Rational(4, 3) * m_next
    b_next = 2 * e_next - sp.Rational(2, 3) * m_next
    h_next = e_next - sp.Rational(2, 3) * m_next
    gamma_next = -8 * e_next + 2 * m_next

    symbolic_identities = {
        "common-x": 2 * h - a,
        "common-z": e + h - b,
        "top-A": 2 * e - 2 * rho - a_next,
        "top-C": 2 * e + 3 * b - rho - b_next,
        "top-scale": 6 * e - 3 * a + 18 * b - 3 * b_next,
        "gamma-C": e + 3 * b - e_next,
        "gamma-D": e - rho - h_next,
        "gamma-scale": -8 * e + 3 * a - 24 * b - gamma_next,
    }
    for label, expression in symbolic_identities.items():
        require_zero(expression, label)

    # Equality in either hybrid bound forces both endpoint equalities.  This
    # finite hostile grid checks the signs independently of the symbolic work.
    for gamma_gap in range(17):
        for z_gap in range(17):
            delta6_gap = gamma_gap + z_gap
            delta8_gap = gamma_gap + 3 * z_gap
            require((delta6_gap == 0) == (gamma_gap == z_gap == 0), ("delta6-equality", gamma_gap, z_gap))
            require((delta8_gap == 0) == (gamma_gap == z_gap == 0), ("delta8-equality", gamma_gap, z_gap))

    # Exhaust admissible small packets and the literal binomial min-gamma
    # face.  Its only term with k=b is ell=h, equal to x^a z^b.
    admissible_count = 0
    for e0 in range(1, 121):
        for m0 in range(0, e0 + 1, 3):
            source = packet(e0, m0)
            admissible_count += 1
            gamma_terms = []
            for ell in range(source.h + 1):
                monomial = (2 * ell, 3 * (source.h - ell), source.e + ell)
                gamma_value = monomial[0] - monomial[1] - 5 * monomial[2]
                require(gamma_value == source.gamma_min, ("gamma-face-weight", source, ell))
                gamma_terms.append(monomial)
            common = [monomial for monomial in gamma_terms if monomial[2] == source.b]
            require(common == [(source.a, 0, source.b)], ("face-intersection", source, common))

    # ------------------------------------------------------------------
    # Route III: exact rung states and coefficient prime exponents.
    # ------------------------------------------------------------------
    l_packet = packet(1, 0)
    h_packet = transform(l_packet)
    j = transform(h_packet)
    g = transform(j)
    r5 = transform(g)
    r6 = transform(r5)
    require((g.e, g.m) == (271, 99), g)
    require((r5.e, r5.m) == (1699, 615), r5)
    require((r6.e, r6.m) == (10663, 3867), r6)

    # Two earlier exact normalizations provide independent scalar controls:
    # L*Norm(L)=H/2^6 and L^7*Norm(H)=J/2^35.  Canonical extraction gives
    # H top/gamma coefficients 2^2*3^24 and 2^2*3^9, and J gives
    # 2^15*3^171 and 2^15*3^72.
    l_to_h_over_64 = scalar_step(l_packet, (0, 3))
    h_exact = {name: (uv[0] + 6, uv[1]) for name, uv in l_to_h_over_64.items()}
    require(h_exact == {"top": (2, 24), "gamma": (2, 9)}, h_exact)
    h_to_j_over_2_35 = scalar_step(h_packet, h_exact["top"])
    j_exact = {name: (uv[0] + 35, uv[1]) for name, uv in h_to_j_over_2_35.items()}
    require(j_exact == {"top": (15, 171), "gamma": (15, 72)}, j_exact)

    j_to_g = scalar_step(j, j_exact["top"])
    g_to_r5 = scalar_step(g, j_to_g["top"])
    r5_to_r6 = scalar_step(r5, g_to_r5["top"])
    require(j_to_g == {"top": (-117, 1128), "gamma": (-117, 513)}, j_to_g)
    require(g_to_r5 == {"top": (-1369, 7251), "gamma": (-1369, 3384)}, g_to_r5)
    require(r5_to_r6 == {"top": (-10493, 46008), "gamma": (-10493, 21753)}, r5_to_r6)
    require(g_to_r5["gamma"][1] + 3 * r5.h == g_to_r5["top"][1], "R5-common-coefficient")
    require(r5_to_r6["gamma"][1] + 3 * r6.h == r5_to_r6["top"][1], "R6-common-coefficient")

    # The wrong monic treatment would use Norm(q)=2 and leave D exponent e,
    # already contradicting the calibrated J->G value h_G=205.
    require(j.e != g.h and j.e - j.rho == g.h, ("nonmonic-hostile", j.e, g.h))

    ledgers = {"G": face_ledger(g), "R5": face_ledger(r5), "R6": face_ledger(r6)}
    require(ledgers["R5"]["max_z"] == (2578, 2988), ledgers["R5"])
    require(ledgers["R5"]["min_gamma"] == (-12362, 1699, 1289), ledgers["R5"])
    require(ledgers["R6"]["max_z"] == (16170, 18748), ledgers["R6"])
    require(ledgers["R6"]["min_gamma"] == (-77570, 10663, 8085), ledgers["R6"])

    determinants = (
        g.e * r5.m - g.m * r5.e,
        r5.e * r6.m - r5.m * r6.e,
    )
    require(determinants == (-1536, 12288), determinants)
    require(determinants[1] == -8 * determinants[0], determinants)

    semantic = {
        "inverse_top": {name: (weight, str(face)) for name, (weight, face) in top.items()},
        "inverse_gamma": {name: (weight, str(face)) for name, (weight, face) in gamma.items()},
        "norm_q": ("2/(27*A^2*C)", "2/D"),
        "symbolic_identities": tuple(sorted(symbolic_identities)),
        "admissible_packets": admissible_count,
        "ledgers": ledgers,
        "scalars": {
            "L_to_H_over_64": l_to_h_over_64,
            "H_exact": h_exact,
            "H_to_J_over_2_35": h_to_j_over_2_35,
            "J_exact": j_exact,
            "J_to_G": j_to_g,
            "G_to_R5": g_to_r5,
            "R5_to_R6": r5_to_r6,
        },
        "determinants": determinants,
        "boundary": "renewal propagation is conditional on polynomiality; the next finite-sheet/polynomiality gate is open",
    }
    semantic_hash = hashlib.sha256(
        json.dumps(semantic, sort_keys=True, separators=(",", ":"), default=str).encode("ascii")
    ).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_PINNED":
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256, ("semantic-hash", semantic_hash))

    print("KELLER FIVE-FACE RENEWAL PROPAGATION -- INDEPENDENT HOSTILE AUDIT")
    print("implementation=candidate-independent; starts from THM-3495 inverse numerators")
    print("inverse_top=(A,B,C*s^3), w=q/s: q^3=2/(27*A^2*C), inverse~(q/s,-s/q,-C*s^6/q^3)")
    print("inverse_gamma=(A*t,B/t,C/t^5), w=q*t: D*q^3-3*B*q-2=0, inverse~(q*t,-1/(q*t),-C/(q^3*t^8))")
    print("nonmonic_vieta=det(mult_q)=2/D; resultant(D*q^3-3*B*q-2,q)=2; leading D must be divided")
    print("hybrid_gaps=delta6_gap=gamma_gap+z_gap; delta8_gap=gamma_gap+3*z_gap; equality iff both gaps vanish")
    print(f"admissible_packet_intersection_hostiles={admissible_count}; e<=120; all PASS")
    print(f"normalization_calibration=H:{h_exact}; J:{j_exact}; matches canonical extracted coefficients")
    print(f"J_to_G_scalars_2_3={j_to_g}")
    print(f"R5_state={(r5.e, r5.m)}; max_z={(r5.a, r5.b)}; min_gamma={(r5.gamma_min, r5.e, r5.h)}")
    print("R5_renewal_scalars=top 3^7251/2^1369; gamma 3^3384/2^1369")
    print(f"R6_state={(r6.e, r6.m)}; max_z={(r6.a, r6.b)}; min_gamma={(r6.gamma_min, r6.e, r6.h)}")
    print("R6_renewal_scalars=top 3^46008/2^10493; gamma 3^21753/2^10493")
    print(f"cassini_determinants={determinants}")
    print("polynomiality_gate=sufficient for renewal once the source packet is complete; not implied by face asymptotics")
    print("boundary=THM-3506+3513 give polynomial R5; THM-3521 gives polynomial R6; next finite-sheet/polynomiality gate remains open")
    print(f"semantic_sha256={semantic_hash}")
    print("VERDICT=SOUND; promotion recommended at fixed conditional scope")
    print("STATUS=PASS")


if __name__ == "__main__":
    main()
