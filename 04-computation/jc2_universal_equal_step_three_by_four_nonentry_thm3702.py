#!/usr/bin/env python3
"""Exact symbolic companion for the THM-3702 all-degree AP peel.

The proof is structural.  This companion checks the six weight fibres, every
differential compression identity, the integrated middle row, and an
independent source-chart realization of the homogeneous bracket formula.
"""

from __future__ import annotations

import ast
import hashlib
from itertools import combinations
from pathlib import Path

import sympy as sp


CHECKS = 0


def gate(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def bracket_coefficient(
    left_weight: int,
    left: sp.Expr,
    right_weight: int,
    right: sp.Expr,
    coordinate: sp.Symbol,
) -> sp.Expr:
    return sp.expand(
        right_weight * sp.diff(left, coordinate) * right
        - left_weight * left * sp.diff(right, coordinate)
    )


def source_bracket(
    left: sp.Expr, right: sp.Expr, x: sp.Symbol, z: sp.Symbol
) -> sp.Expr:
    return sp.expand(
        sp.diff(left, x) * sp.diff(right, z)
        - sp.diff(left, z) * sp.diff(right, x)
    )


def buckets(left: tuple[int, ...], right: tuple[int, ...]):
    result: dict[int, list[tuple[int, int]]] = {}
    for left_weight in left:
        for right_weight in right:
            result.setdefault(left_weight + right_weight + 1, []).append(
                (left_weight, right_weight)
            )
    return {weight: tuple(addresses) for weight, addresses in sorted(result.items())}


def main() -> None:
    window = (-5, -2, 1, 4)
    hard_support = (-5, -2, 1)
    expected_hard_buckets = {
        -9: ((-5, -5),),
        -6: ((-5, -2), (-2, -5)),
        -3: ((-5, 1), (-2, -2), (1, -5)),
        0: ((-5, 4), (-2, 1), (1, -2)),
        3: ((-2, 4), (1, 1)),
        6: ((1, 4),),
    }
    gate(buckets(hard_support, window) == expected_hard_buckets, "hard bucket table")

    three_supports = tuple(combinations(window, 3))
    gate(len(three_supports) == 4, "three-subset count")
    for support in three_supports:
        table = buckets(support, window)
        if support != hard_support:
            gate(table[9] == ((4, 4),), "easy case lost top singleton")

    b = sp.symbols("b")
    alpha, kappa, lam, mu, nu = sp.symbols(
        "alpha kappa lambda mu nu", nonzero=True
    )
    f = sp.Function("f")(b)
    p = sp.Function("p")(b)
    r = sp.Function("r")(b)
    g = sp.Function("g")(b)
    q = sp.Function("q")(b)
    s = sp.Function("s")(b)
    t = sp.Function("t")(b)
    k = sp.Function("k")(b)

    W = lambda u, left, v, right: bracket_coefficient(u, left, v, right, b)

    # Low endpoint and adjacent-row compression after g=lambda*f.
    gate(W(-5, f, -5, lam * f) == 0, "low endpoint proportional control")
    low_row = W(-5, f, -2, q) + W(-2, p, -5, lam * f)
    gate(
        sp.simplify(low_row - W(-5, f, -2, q - lam * p)) == 0,
        "low adjacent compression",
    )
    gate(
        W(-5, alpha * k**5, -2, kappa * k**2) == 0,
        "coprime common-power control",
    )

    # High endpoint and adjacent-row integration.
    gate(W(1, r, 4, mu * r**4) == 0, "high endpoint fourth-power control")
    high_row = W(-2, p, 4, mu * r**4) + W(1, r, 1, s)
    gate(
        sp.simplify(high_row - W(1, r, 1, s - 4 * mu * p * r**3)) == 0,
        "high adjacent compression",
    )
    s_specialized = 4 * mu * p * r**3 + nu * r
    gate(
        sp.simplify(
            W(-2, p, 4, mu * r**4) + W(1, r, 1, s_specialized)
        )
        == 0,
        "high adjacent integrated solution",
    )

    # The remaining nonscalar middle row is one exact derivative after the
    # two endpoint peels and the two adjacent-row integrations.
    f_specialized = alpha * k**5
    g_specialized = lam * f_specialized
    q_specialized = lam * p + kappa * k**2
    middle_row = (
        W(-5, f_specialized, 1, s_specialized)
        + W(-2, p, -2, q_specialized)
        + W(1, r, -5, g_specialized)
    )
    primitive = (
        5 * alpha * (nu - lam) * k * r
        + 20 * alpha * mu * k * p * r**3
        - 2 * kappa * p / k**2
    )
    gate(
        sp.simplify(middle_row - k**4 * sp.diff(primitive, b)) == 0,
        "middle-row derivative identity",
    )

    # Local arm control.  The integrated identity gives p,q order at least
    # two once k has positive order.  Generic polynomial jets verify that
    # every scalar address then has a visible local-uniformizer factor.
    u = sp.symbols("u")
    jet_symbols = sp.symbols("k0:3 p0:3 r0:3")
    k_jet = u * sum(jet_symbols[index] * u**index for index in range(3))
    p_jet = u**2 * sum(jet_symbols[3 + index] * u**index for index in range(3))
    r_jet = sum(jet_symbols[6 + index] * u**index for index in range(3))
    f_jet = alpha * k_jet**5
    q_jet = lam * p_jet + kappa * k_jet**2
    t_jet = mu * r_jet**4
    s_jet = 4 * mu * p_jet * r_jet**3 + nu * r_jet
    W_u = lambda a, left, c, right: bracket_coefficient(a, left, c, right, u)
    scalar_terms = (
        W_u(-5, f_jet, 4, t_jet),
        W_u(-2, p_jet, 1, s_jet),
        W_u(1, r_jet, -2, q_jet),
    )
    for term in scalar_terms:
        gate(sp.expand(term).subs(u, 0) == 0, "scalar address survived arm")

    # Independent hostile verification of the bracket formula on the literal
    # collision source chart for one negative/positive weight pair.
    x, z = sp.symbols("x z")
    b_source = 1 - x**2 * z
    c_source = x
    h_source = 1 - b_source**2
    e_source = sp.expand(h_source / c_source**2)
    left_source = sp.expand(c_source * e_source**3)  # c^-5 h^3
    right_source = sp.expand(c_source**4 * b_source)
    formula_source = sp.expand(
        c_source ** (-5 + 4 + 1)
        * bracket_coefficient(-5, (1 - b**2) ** 3, 4, b, b).subs(b, b_source)
    )
    gate(e_source == sp.expand(z * (2 - x**2 * z)), "source e chart")
    gate(
        source_bracket(left_source, right_source, x, z) == formula_source,
        "hostile source-bracket verification",
    )
    gate(source_bracket(x, z, x, z) == 1, "ambient orientation")

    source = Path(__file__).read_text(encoding="utf-8")
    gate(
        not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
        "inactive Python assert found",
    )

    semantic_rows = [
        f"hard_buckets={expected_hard_buckets}",
        f"supports={three_supports}",
        f"low={sp.factor(low_row)}",
        f"high={sp.factor(high_row)}",
        f"middle={sp.factor(middle_row)}",
        f"primitive={primitive}",
        "arm_scalar=" + ";".join(str(sp.factor(term)) for term in scalar_terms),
        f"source_formula={formula_source}",
    ]
    semantic_hash = hashlib.sha256("\n".join(semantic_rows).encode()).hexdigest()

    print("theorem=THM-3702-universal-equal-step-three-by-four-Danielewski-nonentry")
    print("supports=P:-5,-2,1;Q:-5,-2,1,4;bucket_sizes=1,2,3,3,2,1")
    print("endpoint_peel=g=lambda*f;t=mu*r^4")
    print("adjacent_rows=q-lambda*p=kappa*k^2;f=alpha*k^5;s=4mu*p*r^3+nu*r")
    print("middle_row=k^4*d_db(5alpha(nu-lambda)kr+20alpha*mu*kpr^3-2kappa*p/k^2)")
    print("arm_gate=ord(k)>=1=>ord(p),ord(q)>=2=>scalar_row_vanishes")
    print("controls=source_bracket:PASS;common_power:PASS;arm_jets:PASS")
    print(f"semantic_sha256={semantic_hash}")
    print(f"CHECKS={CHECKS}")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
