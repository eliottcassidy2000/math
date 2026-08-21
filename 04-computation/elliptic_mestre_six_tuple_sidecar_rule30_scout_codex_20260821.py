#!/usr/bin/env python3
"""Exact square-completion sidecar for six-point elliptic searches.

For a monic sextic p with vanishing x^5 coefficient, this derives over Q the
condition making p(x-t)p(x+t) differ from a square by a polynomial of degree
at most four.  For nonzero t it is an equivalence.  An ICARM
Mestre--Fermigier tuple is a positive control; a zero-sum tuple is a hostile;
and the first three Rule 30 right-edge integers supply a deliberately modest
symmetric sampler.  No rank or search enrichment is inferred.
"""

from __future__ import annotations

from fractions import Fraction
import hashlib
import importlib.util
import json
from pathlib import Path

import sympy as sp


DEPENDENCY = Path(__file__).with_name("rule30_right_edge_odometer_thm3458.py")
DEPENDENCY_SHA256 = "8b9a6d029419079f5507d3b153fc43af760e846ddbc93d199392ff4da81640ec"
KNOWN_TUPLE = (348, -600, -216, 492, 876, -900)
KNOWN_SHIFT = Fraction(2326, 23)
HOSTILE_TUPLE = (25, 71, 135, 149, 221, -601)


def require(condition: bool, label: object) -> None:
    if not condition:
        raise RuntimeError(f"requirement failed: {label}")


def load_dependency():
    payload = DEPENDENCY.read_bytes()
    require(hashlib.sha256(payload).hexdigest() == DEPENDENCY_SHA256,
            "THM-3458 dependency hash")
    spec = importlib.util.spec_from_file_location("rule30_thm3458", DEPENDENCY)
    require(spec is not None and spec.loader is not None, "dependency loader")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def square_completion(poly: sp.Expr, variable: sp.Symbol) -> tuple[sp.Expr, sp.Expr]:
    """Return the monic degree-six polynomial part g and r=g^2-poly."""
    coefficients = sp.symbols("u0:6")
    root = variable**6 + sum(
        coefficients[index] * variable ** (5 - index) for index in range(6)
    )
    remainder = sp.Poly(sp.expand(poly - root**2), variable)
    solved: dict[sp.Symbol, sp.Expr] = {}
    for degree in range(11, 5, -1):
        equation = sp.expand(remainder.coeff_monomial(variable**degree).subs(solved))
        unknown = coefficients[11 - degree]
        solutions = sp.solve(equation, unknown)
        require(len(solutions) == 1, (degree, "unique square-completion coefficient"))
        solved[unknown] = solutions[0]
    root = sp.expand(root.subs(solved))
    return root, sp.expand(root**2 - poly)


def tuple_record(values: tuple[int, ...], shift: Fraction) -> tuple[object, ...]:
    x = sp.symbols("x")
    t_value = sp.Rational(shift.numerator, shift.denominator)
    sextic = sp.Poly(sp.prod(x - value for value in values), x)
    coefficients = tuple(int(item) for item in sextic.all_coeffs())
    require(len(coefficients) == 7 and coefficients[1] == -sum(values),
            (values, "sextic coefficients"))
    c4, c3, c1 = coefficients[2], coefficients[3], coefficients[5]
    sidecar = 2 * c1 - c3 * c4
    product = sp.expand(sextic.as_expr().subs(x, x - t_value)
                        * sextic.as_expr().subs(x, x + t_value))
    root, quartic = square_completion(product, x)
    degree = sp.Poly(quartic, x).degree()
    base_x = tuple(sp.Rational(value) + sign * t_value
                   for value in values for sign in (-1, 1))
    require(len(set(base_x)) == 12, (values, "distinct base abscissae"))
    for x_value in base_x:
        require(sp.expand(quartic.subs(x, x_value) - root.subs(x, x_value)**2) == 0,
                (values, x_value, "base point"))
    squarefree = sp.gcd(sp.Poly(quartic, x), sp.Poly(sp.diff(quartic, x), x)).degree() == 0
    base_y_nonzero = all(root.subs(x, x_value) != 0 for x_value in base_x)
    coefficient_payload = tuple(str(item) for item in sp.Poly(quartic, x).all_coeffs())
    quartic_sha256 = hashlib.sha256(
        json.dumps(coefficient_payload, separators=(",", ":")).encode("ascii")
    ).hexdigest()
    return (
        values,
        str(shift),
        sum(values),
        sidecar,
        degree,
        len(base_x),
        squarefree,
        base_y_nonzero,
        quartic_sha256,
    )


def main() -> None:
    x, t = sp.symbols("x t")
    c4, c3, c2, c1, c0 = sp.symbols("c4 c3 c2 c1 c0")
    sextic = x**6 + c4 * x**4 + c3 * x**3 + c2 * x**2 + c1 * x + c0
    product = sp.expand(sextic.subs(x, x - t) * sextic.subs(x, x + t))
    _, remainder = square_completion(product, x)
    remainder_poly = sp.Poly(remainder, x)
    require(remainder_poly.degree() == 5, "generic remainder degree")
    x5 = sp.factor(remainder_poly.coeff_monomial(x**5))
    require(x5 == -12 * t**2 * (2 * c1 - c3 * c4), "codimension-one sidecar")

    known = tuple_record(KNOWN_TUPLE, KNOWN_SHIFT)
    require(known[2:8] == (0, 0, 4, 12, True, True),
            "ICARM curve 211 positive control")

    module = load_dependency()
    edge = tuple(module.direct_rows(2))
    require(edge == (1, 7, 25), "Rule 30 right-edge source")
    symmetric_tuple = edge + tuple(-value for value in edge)
    symmetric = tuple_record(symmetric_tuple, Fraction(2))
    require(symmetric[2:7] == (0, 0, 4, 12, True),
            "symmetric Rule 30 sampler")

    hostile = tuple_record(HOSTILE_TUPLE, Fraction(11))
    require(hostile[2] == 0 and hostile[3] != 0 and hostile[4:7] == (5, 12, True),
            "zero-sum hostile")

    semantic = (
        DEPENDENCY_SHA256,
        str(x5),
        known,
        symmetric,
        hostile,
    )
    semantic_sha256 = hashlib.sha256(
        json.dumps(semantic, separators=(",", ":")).encode("ascii")
    ).hexdigest()

    print("ELLIPTIC_MESTRE_SIX_TUPLE_SIDECAR_RULE30_SCOUT_20260821")
    print("status=PROVED_SYMBOLIC_IDENTITY+FINITE_EXACT_CONTROLS;no_rank_claim")
    print(f"dependency_sha256={DEPENDENCY_SHA256}")
    print("identity=[x^5](g^2-p(x-t)p(x+t))=-12*t^2*(2*c1-c3*c4)")
    print("universe=Q;for_t_nonzero:degree_remainder_le_4_iff_e1=0_and_2*e5=e2*e3")
    print("icarm_curve_211_mestre_fermigier_control=" + repr(known).replace(" ", ""))
    print("rule30_symmetric_sampler=" + repr(symmetric).replace(" ", ""))
    print("zero_sum_hostile=" + repr(hostile).replace(" ", ""))
    print(f"semantic_sha256={semantic_sha256}")
    print("scope=quartic_and_twelve_base_abscissae_only;no_independence;no_rank_enrichment")
    print("reproduce=python3_or_python3_-O_04-computation/elliptic_mestre_six_tuple_sidecar_rule30_scout_codex_20260821.py")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
