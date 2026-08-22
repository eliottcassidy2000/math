#!/usr/bin/env python3
"""Exact two-clutch tangent audit of the THM-3306/3309 critical deck.

The fixed slice ``C=c+x, d=k=1`` has a transverse linear-subresultant base
ideal and a nonsplit quadratic common-root deck.  This scout releases the two
physical affine clutch slopes by

    C=c+d*x,       E'=k,

and differentiates the *generic* subresultant row in the ``d`` and ``k``
directions.  In each THM-3212 accessory field it solves the exact implicit
tangent equations for the moving base point ``(x,c)``.  It then checks whether
the two clutch tangents are independent and whether the fixed deck units needed
for algebraic etale and formal persistence remain present.

This is local algebraic geometry at the declared degree-36 points.  It does
not construct a Keller mate, an inverse, a global component, or a
JC(2)/DC(2) consequence.
"""

from __future__ import annotations

import ast
from hashlib import sha256
from pathlib import Path

import sympy as sp
from sympy.polys.domains import QQ


ROOT = Path(__file__).resolve().parents[1]
DEPENDENCIES = {
    "01-canon/theorems/THM-3306-affine-c-critical-section-square-discriminant-and-transverse-base-locus.md":
        "514b5c07a70ea0dd0020857af3278008c3df2a03af6826784d96d059a3b26111",
    "01-canon/theorems/THM-3309-exceptional-quadratic-deck-passport-and-gradient-unimodularity-obstruction.md":
        "9d16b3880af8b5c3fadab68b8654bbb3159d9e1e641eaf67a960aec9c57a3f95",
}
EXPECTED_DELTA_DIGESTS = {
    "4111": "6b66fda17d31412a565a23649197ca94e9372c7368f01a31f60da61714b328fc",
    "3211": "fa7ce86608c083d3e8e722cf5de425a9a3d3ae3ddbdbc5e829beb6e10fb7187a",
}
EXPECTED_TANGENT_PACKETS = {
    "4111": (
        "0a36343b73059ad745d5efa1e9d665c7b268ce678c983a9b584ab7730ea5e158",
        "33078177dcd13534621da151a1bdea693945851da8c1efe73bfc294fbda467e7",
        "3f0a0d1d6385099c25f6fa42aaeb11be990d59dba2260c56f9c0c1395f3b63ec",
    ),
    "3211": (
        "1ddc493057a13c1401da3ec78a761ae37fa330be18d4f40c54e201aaf9b0d31b",
        "3923cbc795af7af9e0aa9145e55e55c54f9113c03911ab05cbf0a615cfab8c5c",
        "1aba7e3b235a75045ab90bbeaeb5e58c00ef8f0af2cc23854ead0ca4fb1b65d1",
    ),
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def lf_hash(path: Path) -> str:
    payload = path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return sha256(payload).hexdigest()


def element_text(element) -> str:
    return "|".join(
        f"{monomial}:{coefficient}"
        for monomial, coefficient in sorted(element.to_dict().items())
    )


def element_digest(element) -> str:
    return sha256(element_text(element).encode("utf-8")).hexdigest()


def packet_digest(elements) -> str:
    return sha256("\n--\n".join(element_text(item) for item in elements).encode("utf-8")).hexdigest()


def specialize_sparse(poly, base_values, polynomial_ring):
    result = polynomial_ring.zero
    for monomial, coefficient in poly.terms():
        term = polynomial_ring.convert(coefficient)
        for base, exponent in zip(base_values, monomial):
            if exponent:
                term *= base**exponent
        result += term
    return result


def inverse_mod(element, modulus):
    inverse, _, common = element.gcdex(modulus)
    require(common.degree() == 0, "requested quotient element is not a unit")
    return (inverse / common.LC) % modulus


def evaluate_c_mod(poly, value, modulus):
    x_generator = modulus.ring.gens[0]
    result = modulus.ring.zero
    for (x_degree, c_degree), coefficient in poly.to_dict().items():
        result += (
            modulus.ring.ground_new(coefficient)
            * x_generator**x_degree
            * value**c_degree
        )
    return result % modulus


def derive_generic_rows():
    y, response_A, C, V, derivative_V, d, k = sp.symbols("y A C V Vp d k")
    L = y**2 + y + C * V
    R1 = sp.expand(2 * L * (2 * y + 1) + V * response_A)
    R2 = sp.expand(
        V**3 * k + V**2 * y + L * (-derivative_V * y + 2 * V**2 * d)
    )
    rows = sp.subresultants(R1, R2, y)
    require(
        tuple(sp.Poly(row, y).degree() for row in rows) == (3, 3, 2, 1, 0),
        "generic cubic subresultant profile changed",
    )
    linear = sp.Poly(rows[-2], y)
    raw_a = sp.Poly(linear.nth(1), response_A, C, V, derivative_V, d, k)
    raw_b = sp.Poly(linear.nth(0), response_A, C, V, derivative_V, d, k)
    raw_resultant = sp.Poly(
        sp.Poly(rows[-1], y).nth(0), response_A, C, V, derivative_V, d, k
    )
    quadratic = sp.Poly(rows[-3], y)
    raw_quadratic = tuple(
        sp.Poly(quadratic.nth(degree), response_A, C, V, derivative_V, d, k)
        for degree in (2, 1, 0)
    )
    require(raw_a.degree(C) == 0 and raw_b.degree(C) == 1, "generic C degrees")
    require(sp.Poly(R1, y).LC() == 4, "projective y-infinity hostile")
    return {
        "a": raw_a,
        "b": raw_b,
        "resultant": raw_resultant,
        "quadratic": raw_quadratic,
        "A_symbol": response_A,
        "C_symbol": C,
        "d_symbol": d,
        "k_symbol": k,
    }


def response_case(name: str):
    u, x, c = sp.symbols("u x c")
    if name == "4111":
        accessory = sp.Poly(100 * u**3 + 244 * u**2 + 237 * u + 44, u, domain=QQ)
        exponent_a, exponent_b = 4, 1
    elif name == "3211":
        accessory = sp.Poly(75 * u**3 - 89 * u**2 - 31 * u + 61, u, domain=QQ)
        exponent_a, exponent_b = 3, 2
    else:
        raise RuntimeError(f"unknown response case {name}")

    field = QQ.alg_field_from_poly(accessory, alias="u")
    alpha = field.ext
    ring = field.poly_ring(x, c)
    X, CP = ring.gens
    if name == "4111":
        accessory_v = (8 * alpha**2 + 9 * alpha + 8) / 7
        shift = 5 * (alpha + 1) / 7
        A0 = 80 * accessory_v**2 * (alpha + 1) / 343
        extras = (9, 0)
    else:
        accessory_v = (24 * alpha**2 - 16 * alpha - 16) / 21
        shift = (5 * alpha - 4) / 7
        A0 = 9 * accessory_v**2 * (5 * alpha - 4) / 343
        extras = (6, 3)

    gamma = -7 * A0
    quadratic = X**2 - alpha * X + accessory_v
    Dsource = X**exponent_a * (X - 1) ** exponent_b * quadratic
    T = X * (X - 1) * quadratic
    Esource = (
        exponent_a * (X - 1) * quadratic
        + exponent_b * X * quadratic
        + X * (X - 1) * (2 * X - alpha)
    ) / 7
    S = X + shift
    V = 4 * S * Dsource * T**2 / gamma**2
    Aresponse = 2 * S * Esource * T / gamma
    boundary = S**3 * T**8 * X**extras[0] * (X - 1) ** extras[1]
    require((V.degree(X), Aresponse.degree(X), boundary.degree(X)) == (16, 8, 44), name)
    require(2 * V * Aresponse.diff(X) - Aresponse * V.diff(X) == 2 * V, name)
    return {
        "name": name,
        "field": field,
        "ring": ring,
        "X": X,
        "CP": CP,
        "V": V,
        "Aresponse": Aresponse,
        "g": S * T,
        "boundary": boundary,
    }


def total_parameter_derivative(raw, parameter: str, generic, case, values):
    ring = case["ring"]
    X = case["X"]
    if parameter == "d":
        # C=c+d*x, so the physical d derivative has the chain-rule C term.
        direct = specialize_sparse(raw.diff(generic["d_symbol"]), values, ring)
        through_C = specialize_sparse(raw.diff(generic["C_symbol"]), values, ring)
        return direct + X * through_C
    if parameter == "k":
        return specialize_sparse(raw.diff(generic["k_symbol"]), values, ring)
    raise RuntimeError(f"unknown clutch parameter {parameter}")


def central_parameter_difference(raw, parameter: str, case):
    """Independent exact derivative at (d,k)=(1,1).

    After C=c+d*x, the physical d degree is at most two; the k degree is at
    most one.  Thus the central difference between parameter values zero and
    two equals the derivative at one in either direction.
    """

    ring = case["ring"]
    X = case["X"]
    CP = case["CP"]
    V = case["V"]
    derivative_V = V.diff(X)
    Aresponse = case["Aresponse"]
    if parameter == "d":
        low = (Aresponse, CP, V, derivative_V, ring.zero, ring.one)
        high = (Aresponse, CP + 2 * X, V, derivative_V, 2 * ring.one, ring.one)
    elif parameter == "k":
        low = (Aresponse, CP + X, V, derivative_V, ring.one, ring.zero)
        high = (Aresponse, CP + X, V, derivative_V, ring.one, 2 * ring.one)
    else:
        raise RuntimeError(f"unknown clutch parameter {parameter}")
    return (specialize_sparse(raw, high, ring) - specialize_sparse(raw, low, ring)) / 2


def audit_case(name: str, generic) -> tuple[str, ...]:
    case = response_case(name)
    field = case["field"]
    ring = case["ring"]
    X = case["X"]
    CP = case["CP"]
    V = case["V"]
    Aresponse = case["Aresponse"]
    derivative_V = V.diff(X)
    boundary = case["boundary"]

    values = (Aresponse, CP + X, V, derivative_V, ring.one, ring.one)
    primitive_a = specialize_sparse(generic["a"], values, ring).exquo(boundary)
    primitive_b = specialize_sparse(generic["b"], values, ring).exquo(boundary)
    common_unit = primitive_a.LC
    a_bivariate = primitive_a / common_unit
    b_bivariate = primitive_b / common_unit
    a = a_bivariate.evaluate(CP, field.zero).monic()
    b0 = b_bivariate.evaluate(CP, field.zero)
    b1 = b_bivariate.diff(CP).evaluate(CP, field.zero)
    require(
        (a.degree(), b0.degree(), b1.degree()) == (36, 44, 36),
        f"{name} inherited primitive row degrees",
    )
    require(a.factor_list()[1][0][0].degree() == 36, f"{name} irreducible a")
    cstar = (-b0 * inverse_mod(b1, a)) % a

    x_ring = a.ring
    XX = x_ring.gens[0]
    Vx = V.evaluate(CP, field.zero)
    Aresponse_x = Aresponse.evaluate(CP, field.zero)
    gx = case["g"].evaluate(CP, field.zero)
    derivative_Vx = Vx.diff(XX)

    a_x = evaluate_c_mod(a_bivariate.diff(X), cstar, a)
    b_x = evaluate_c_mod(b_bivariate.diff(X), cstar, a)
    b_c = evaluate_c_mod(b_bivariate.diff(CP), cstar, a)
    require(
        all(a.gcd(item).degree() == 0 for item in (a_x, b_c, Vx, gx)),
        f"{name} transverse/localization units",
    )
    inverse_a_x = inverse_mod(a_x, a)
    inverse_b_c = inverse_mod(b_c, a)

    tangents = {}
    frozen_failures = {}
    for parameter in ("d", "k"):
        a_parameter_numerator = total_parameter_derivative(
            generic["a"], parameter, generic, case, values
        )
        b_parameter_numerator = total_parameter_derivative(
            generic["b"], parameter, generic, case, values
        )
        require(
            a_parameter_numerator
            == central_parameter_difference(generic["a"], parameter, case)
            and b_parameter_numerator
            == central_parameter_difference(generic["b"], parameter, case),
            f"{name} {parameter} symbolic/finite-difference tangent disagreement",
        )
        a_parameter_raw = a_parameter_numerator.exquo(boundary)
        b_parameter_raw = b_parameter_numerator.exquo(boundary)
        a_parameter = evaluate_c_mod(a_parameter_raw / common_unit, cstar, a)
        b_parameter = evaluate_c_mod(b_parameter_raw / common_unit, cstar, a)
        xdot = (-a_parameter * inverse_a_x) % a
        cdot = (-(b_x * xdot + b_parameter) * inverse_b_c) % a
        require(
            (a_x * xdot + a_parameter) % a == a.ring.zero,
            f"{name} {parameter} first implicit equation",
        )
        require(
            (b_x * xdot + b_c * cdot + b_parameter) % a == a.ring.zero,
            f"{name} {parameter} second implicit equation",
        )
        require(
            a_parameter != a.ring.zero or b_parameter != a.ring.zero,
            f"{name} {parameter} frozen-base hostile must move",
        )
        require(
            xdot != a.ring.zero and cdot != a.ring.zero,
            f"{name} {parameter} tangent must move both base coordinates",
        )
        tangents[parameter] = (xdot, cdot)
        frozen_failures[parameter] = (a_parameter, b_parameter)

    tangent_det = (
        tangents["d"][0] * tangents["k"][1]
        - tangents["k"][0] * tangents["d"][1]
    ) % a
    require(tangent_det != a.ring.zero, f"{name} clutch tangent rank")

    # Redundant consequence path: the resultant also has zero directional
    # derivative along each implicit tangent.  The all-order implication uses
    # the standard subresultant recurrence after the quadratic leading row is
    # made a unit.
    resultant_bivariate = specialize_sparse(generic["resultant"], values, ring)
    resultant_value = evaluate_c_mod(resultant_bivariate, cstar, a)
    require(resultant_value == a.ring.zero, f"{name} base resultant")
    resultant_x = evaluate_c_mod(resultant_bivariate.diff(X), cstar, a)
    resultant_c = evaluate_c_mod(resultant_bivariate.diff(CP), cstar, a)
    for parameter in ("d", "k"):
        resultant_parameter = evaluate_c_mod(
            total_parameter_derivative(
                generic["resultant"], parameter, generic, case, values
            ),
            cstar,
            a,
        )
        xdot, cdot = tangents[parameter]
        require(
            (
                resultant_x * xdot
                + resultant_c * cdot
                + resultant_parameter
            )
            % a
            == a.ring.zero,
            f"{name} {parameter} resultant tangent",
        )

    quadratic_value = tuple(
        evaluate_c_mod(specialize_sparse(item, values, ring), cstar, a)
        for item in generic["quadratic"]
    )
    quadratic_nonzero_pattern = tuple(
        int(item != a.ring.zero) for item in quadratic_value
    )
    require(
        quadratic_nonzero_pattern == (1, 1, 1),
        f"{name} quadratic subresultant must survive",
    )

    Cstar = XX + cstar
    P2 = (2 * derivative_Vx + 8 * Vx**2) % a
    Q2 = (2 * derivative_Vx + 12 * Vx**2) % a
    R2 = (
        Vx
        * (
            4 * Vx**2
            + 8 * Cstar * Vx**2
            + derivative_Vx * (2 * Cstar + Aresponse_x)
        )
    ) % a
    delta = (Q2**2 - 4 * P2 * R2) % a
    finite_clutch = (x_ring.one - Aresponse_x.diff(XX) * Cstar) % a
    require(
        all(a.gcd(item).degree() == 0 for item in (P2, delta, finite_clutch)),
        f"{name} deck and finite-clutch units",
    )
    require(element_digest(delta) == EXPECTED_DELTA_DIGESTS[name], f"{name} deck drift")

    print(
        f"case={name};base_field_degree=36;linear_row=(a,b);"
        "base_Jacobian_det=a_x*b_c_is_unit=1;"
        "V(a,b)_etale_over_(d,k)_at_the_point=1;"
        "completion_parameters=(d-1,k-1);formal_solution=(x(d,k),c(d,k))_unique"
    )
    for parameter in ("d", "k"):
        xdot, cdot = tangents[parameter]
        print(
            f"case={name};parameter={parameter};"
            f"frozen_base_failure_digest={packet_digest(frozen_failures[parameter])};"
            f"xdot_degree={xdot.degree()};cdot_degree={cdot.degree()};"
            f"tangent_digest={packet_digest((xdot, cdot))};"
            f"xdot_nonzero={int(xdot != a.ring.zero)};cdot_nonzero={int(cdot != a.ring.zero)}"
        )
    print(
        f"case={name};two_clutch_tangent_det_nonzero={int(tangent_det != a.ring.zero)};"
        f"two_clutch_tangent_det_degree={tangent_det.degree()};"
        f"two_clutch_tangent_det_digest={element_digest(tangent_det)}"
    )
    print(
        f"case={name};P2_delta_finite_clutch_are_units=1;"
        f"delta_digest={element_digest(delta)};"
        f"quadratic_subresultant_nonzero_pattern={quadratic_nonzero_pattern};"
        "resultant_tangent_vanishes_in_both_clutch_directions=1;"
        "nonsplit_residue_deck_lifts_to_a_connected_finite_etale_C2_cover_on_the_local_germ;"
        "quadratic_subresultant_and_cubic_degrees_persist_formally;"
        "critical_C2_deck_and_gradient_vanishing_persist_formally"
    )
    packet = (
        packet_digest(tangents["d"]),
        packet_digest(tangents["k"]),
        element_digest(tangent_det),
    )
    require(packet == EXPECTED_TANGENT_PACKETS[name], f"{name} tangent packet drift")
    return packet


def source_audit() -> None:
    tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    assertion_nodes = sum(isinstance(node, ast.Assert) for node in ast.walk(tree))
    float_literals = sum(
        isinstance(node, ast.Constant) and isinstance(node.value, float)
        for node in ast.walk(tree)
    )
    require(assertion_nodes == 0 and float_literals == 0, "source exactness audit")
    print(f"source_ast=(assert_nodes={assertion_nodes},float_literals={float_literals})")


def main() -> None:
    print("JC EXCEPTIONAL QUADRATIC TWO-CLUTCH ALGEBRAIC ETALE PERSISTENCE SCOUT")
    print("status=VERIFIED-EXACT_TANGENT_PLUS_PROVED_ETALE_GERM_CONSEQUENCE")
    for relative_path, expected in DEPENDENCIES.items():
        actual = lf_hash(ROOT / relative_path)
        require(actual == expected, f"dependency drift: {relative_path}")
        print(f"dependency_sha256={actual}  {relative_path}")
    generic = derive_generic_rows()
    packets = tuple(audit_case(name, generic) for name in ("4111", "3211"))
    require(packets[0] != packets[1], "accessory cases must remain distinguishable")
    print(
        "mechanism=transverse_(x,c)_base_ideal_absorbs_both_physical_clutch_slopes;"
        "the_relative_Jacobian_criterion_produces_an_algebraic_etale_germ;"
        "completion_gives_the_unique_formal_lift;"
        "finite_etale_idempotents_over_the_local_base_are_detected_on_the_residue_deck"
    )
    print(
        "hostile=frozen_(x,c)_base_fails_in_each_d_and_k_direction;"
        "repair=move_the_point_by_the_exact_implicit_tangent;"
        "gradient_gate=still_proper_on_the_algebraic_C2_deck"
    )
    print(
        "scope=two_THM3212_accessory_fields;etale_germ_over_(d,k)=(1,1);"
        "tangent_paths=symbolic_chain_rule_plus_exact_central_difference;"
        "no_global_component_no_owner_wall_classification_no_mate_no_inverse_no_JC2_no_DC2"
    )
    source_audit()
    print("ALL EXACT TWO-CLUTCH ETALE-PERSISTENCE CHECKS PASSED")


if __name__ == "__main__":
    main()
