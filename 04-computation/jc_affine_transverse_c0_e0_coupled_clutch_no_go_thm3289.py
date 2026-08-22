#!/usr/bin/env python3
"""Exact scout for simultaneous affine C_0,E_0 clutches with B=1.

For either cubic accessory response pair of THM-3212, this companion studies

    P(x,z) = (V(x) z^2 + z + C(x))^2 + A(x) z + E(x),

with C'=d and E'=k constant.  The lanes d=0 and k=0 are inherited from
THM-3212 and THM-3279, respectively.  In the genuinely coupled lane d*k!=0,
the script derives the localized resultant with a literal Sylvester
determinant, proves the finite/T/infinity ledgers, derives the two live S
walls by truncated exact series arithmetic, and proves by exact cubic-field
Euclidean algorithms that the fifth and sixth S jets have no common zero.

This is a critical-point obstruction for one explicit first-coordinate
family.  It constructs no inverse cover and proves no instance of JC(2).
"""

from __future__ import annotations

import ast
from hashlib import sha256
from pathlib import Path

import sympy as sp
from sympy.polys.domains import QQ


ROOT = Path(__file__).resolve().parents[1]
DEPENDENCIES = {
    "04-computation/jc_centered_heptic_source_morse_obstruction_thm3212.py":
        "03d121e57dd2edaece7cd8693d792349f03757c6e781eb5d9d0c897fcc889448",
    "05-knowledge/results/jc_centered_heptic_source_morse_obstruction_thm3212.out":
        "729e0c7b9fa51fa5c4ac5f18f50dc4413c8a8bb7bf5f0ebf2a7709650304bc85",
    "04-computation/jc_heptic_constant_parameter_discriminant_audit_thm3212.py":
        "adeb3c548d5fe3966eefc7ec4eeadfd1410a62356eca8a6c387e39cbe8fc6122",
    "05-knowledge/results/jc_heptic_constant_parameter_discriminant_audit_thm3212.out":
        "d170cf2212848ef76722579a40b65506bedf6e65a031012ca06c27dcd1ef77bb",
    "04-computation/jc_affine_transverse_c0_clutch_no_go_thm3279.py":
        "06820b2476fc0f2cefe3982d054a7db09bb88b4892503580550a1b154564508a",
    "05-knowledge/results/jc_affine_transverse_c0_clutch_no_go_thm3279.out":
        "4a88b5ab31eed4c9a5f90f814a6a24a614db0afecc4cf1cab7fa32dae7c991c4",
}
EXPECTED_CASE_DIGESTS = {
    "4111": (
        "92c368a8101e6ccebc69c4961bda10b41a85dd3f367bdbfeb7ee04747421c319",
        "ae98818930edc795b9a7a73a957ddf2f2e70dee40695bb3ea489c9a8c16780d4",
    ),
    "3211": (
        "f7268035c556bfdc85afd670405986ddbe6839d2ddfb3a13900cba2baabc901a",
        "621f57adfa31900c88fa44d0443d9baf81d199a201a4b9963d99914d134b1d21",
    ),
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def lf_hash(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def literal_cubic_sylvester(first, second, variable):
    """Return the literal 6-by-6 Sylvester determinant of two cubics."""

    first_poly = sp.Poly(first, variable)
    second_poly = sp.Poly(second, variable)
    require(first_poly.degree() == 3, "first Sylvester input is cubic")
    require(second_poly.degree() == 3, "second Sylvester input is cubic")
    first_row = first_poly.all_coeffs()
    second_row = second_poly.all_coeffs()
    zero = sp.S.Zero
    matrix = sp.Matrix(
        [
            first_row + [zero, zero],
            [zero] + first_row + [zero],
            [zero, zero] + first_row,
            second_row + [zero, zero],
            [zero] + second_row + [zero],
            [zero, zero] + second_row,
        ]
    )
    return sp.expand(matrix.det(method="domain-ge"))


def derive_coupled_factor():
    y, A, C, V, derivative_V, d, k = sp.symbols("y A C V Vp d k")
    derivative_A = sp.symbols("Ap")
    L = y**2 + y + C * V
    R1 = 2 * L * (2 * y + 1) + V * A
    R2 = V**3 * k + V**2 * y + L * (-derivative_V * y + 2 * V**2 * d)

    resultant = literal_cubic_sylvester(R1, R2, y)
    K = sp.cancel(resultant / V**3)
    require(sp.expand(resultant - V**3 * K) == 0, "exact V^3 resultant factor")
    K_poly = sp.Poly(K, A, C, V, derivative_V, d, k)
    require(len(K_poly.terms()) == 40, "coupled critical factor has 40 terms")

    raw = (
        2 * L * (derivative_V * y**2 + d * V**2)
        + derivative_A * y * V**2
        + k * V**3
    )
    response_A_prime = (2 * V + A * derivative_V) / (2 * V)
    reduction_defect = sp.cancel(
        raw.subs(derivative_A, response_A_prime)
        - R2
        - derivative_V * y * R1 / 2
    )
    require(reduction_defect == 0, "localized gradient reduction sign")

    weights = (8, 1, 16, 15, 0, 0)
    degree_rows = [
        (sum(exponent * weight for exponent, weight in zip(monomial, weights)),
         monomial, coefficient)
        for monomial, coefficient in K_poly.terms()
    ]
    top_rows = [row for row in degree_rows if row[0] == 97]
    require(len(top_rows) == 1, "unique coupled degree-97 monomial")
    require(
        top_rows[0] == (97, (0, 1, 6, 0, 1, 2), sp.Integer(128)),
        "degree-97 term is 128*C*V^6*d*k^2",
    )
    print(
        "localized=raw-R2=(Vp*y/2)*R1;"
        "resultant=V^3*K40;"
        "unique_degree97=128*C*V^6*d*k^2;"
        "affine_lead=128*d^2*k^2*V_lead^6"
    )
    return K_poly, (A, C, V, derivative_V, d, k)


def t_boundary_checks(K_poly, symbols) -> None:
    A, C, V, derivative_V, d_symbol, k_symbol = symbols
    t, v, c, d, k = sp.symbols("t v c d k", nonzero=True)
    rows = []
    for multiplicity in (3, 4, 5, 6):
        response_slope = sp.Rational(2, 2 - multiplicity)
        specialization = {
            A: response_slope * t,
            C: c + d * t,
            V: v * t**multiplicity,
            derivative_V: multiplicity * v * t ** (multiplicity - 1),
            d_symbol: d,
            k_symbol: k,
        }
        local = sp.Poly(sp.expand(K_poly.as_expr().subs(specialization)), t)
        target = 3 * multiplicity - 1
        expected = (
            sp.Rational(16 * multiplicity * (multiplicity - 1), multiplicity - 2)
            * v**3
            * (k - response_slope * c)
        )
        require(
            all(local.nth(index) == 0 for index in range(target)),
            f"T lower rows vanish for multiplicity {multiplicity}",
        )
        require(
            sp.factor(local.nth(target) - expected) == 0,
            f"T gate row for multiplicity {multiplicity}",
        )
        rows.append((multiplicity, target))
    print(
        f"T_boundary_rows={tuple(rows)};"
        "lead=16*m*(m-1)/(m-2)*v^3*(k-A_prime*c);"
        "gate=finite_clutch"
    )


def series_multiply(first, second, order):
    result = [sp.S.Zero] * (order + 1)
    for left_index, left in enumerate(first):
        if left == 0:
            continue
        for right_index, right in enumerate(second[:order + 1 - left_index]):
            if right != 0:
                result[left_index + right_index] += left * right
    return result


def series_power(series, exponent, order):
    result = [sp.S.One] + [sp.S.Zero] * order
    for _ in range(exponent):
        result = series_multiply(result, series, order)
    return result


def response_coefficients(v_jet, order):
    """Solve 2VA'-AV'=2V recursively at a simple V zero."""

    coefficients = [sp.S.Zero] * (order + 2)
    coefficients[1] = sp.Integer(2)
    for degree in range(2, order + 2):
        subtotal = sp.S.Zero
        for v_index in range(1, degree + 1):
            a_index = degree + 1 - v_index
            if 1 <= a_index <= degree - 1:
                subtotal += (
                    (2 * a_index - v_index)
                    * v_jet[v_index]
                    * coefficients[a_index]
                )
        coefficients[degree] = sp.factor(
            (2 * v_jet[degree] - subtotal)
            / ((2 * degree - 1) * v_jet[1])
        )
    return coefficients


def wall_polynomials(c, v1, v2, v3, v4, v5):
    """Return the live q5 factor and q6 numerator after the q3/q4 walls."""

    slope = 3 * c * v1**2 + v2
    q5_core = (
        -210 * c**2 * v1**4 * v2
        + 630 * c * v1**5
        + 504 * c * v1**3 * v3
        - 518 * c * v1**2 * v2**2
        + 105 * v1**3 * v2
        - 360 * v1**2 * v4
        + 750 * v1 * v2 * v3
        - 400 * v2**3
    )
    q6_core = (
        12600 * c**4 * v1**8 * v2
        + 37800 * c**3 * v1**9
        - 37800 * c**3 * v1**7 * v3
        - 13650 * c**3 * v1**6 * v2**2
        + 189000 * c**2 * v1**7 * v2
        + 48600 * c**2 * v1**6 * v4
        + 32310 * c**2 * v1**5 * v2 * v3
        - 100420 * c**2 * v1**4 * v2**3
        + 18900 * c * v1**8
        - 9450 * c * v1**6 * v3
        + 97125 * c * v1**5 * v2**2
        - 21000 * c * v1**5 * v5
        - 28800 * c * v1**4 * v2 * v4
        - 504 * c * v1**4 * v3**2
        + 172056 * c * v1**3 * v2**2 * v3
        - 124946 * c * v1**2 * v2**4
        + 1575 * v1**6 * v2
        - 16200 * v1**5 * v4
        + 25875 * v1**4 * v2 * v3
        - 1025 * v1**3 * v2**3
        - 7000 * v1**3 * v2 * v5
        + 9720 * v1**3 * v3 * v4
        - 23640 * v1**2 * v2**2 * v4
        - 15882 * v1**2 * v2 * v3**2
        + 73098 * v1 * v2**3 * v3
        - 37168 * v2**5
    )
    return slope, q5_core, slope * q5_core, q6_core


def s_wall_checks(K_poly, symbols) -> None:
    A, C, V, derivative_V, d_symbol, k_symbol = symbols
    order = 6
    v_symbols = [None] + list(sp.symbols("v1:8"))
    v1, v2, v3, v4, v5 = v_symbols[1:6]
    c, d, k = sp.symbols("c d k")
    a_jet = response_coefficients(v_symbols, order)

    base_series = [
        [sp.S.Zero] + a_jet[1:order + 1],
        [c, d] + [sp.S.Zero] * (order - 1),
        [sp.S.Zero] + v_symbols[1:order + 1],
        [index * v_symbols[index] for index in range(1, order + 2)],
        [d] + [sp.S.Zero] * order,
        [k] + [sp.S.Zero] * order,
    ]
    powers = {
        (base_index, exponent): series_power(series, exponent, order)
        for base_index, series in enumerate(base_series)
        for exponent in range(7)
    }
    q = [sp.S.Zero] * (order + 1)
    for monomial, coefficient in K_poly.terms():
        term_series = [sp.S.One] + [sp.S.Zero] * order
        for base_index, exponent in enumerate(monomial):
            term_series = series_multiply(
                term_series, powers[(base_index, exponent)], order
            )
        for degree in range(order + 1):
            q[degree] += coefficient * term_series[degree]

    expected_q3 = (
        sp.Rational(8, 3)
        * v1**2
        * (2 * c - k)
        * (6 * c * v1**2 + 3 * k * v1**2 + 4 * v2)
    )
    require(
        all(sp.factor(q[index]) == 0 for index in range(3)),
        "universal S order at least three",
    )
    require(sp.factor(q[3] - expected_q3) == 0, "factored S q3")

    live_k = -2 * c - sp.Rational(4, 3) * v2 / v1**2
    slope = 3 * c * v1**2 + v2
    d_wall_numerator = (
        15 * v1**3
        - 18 * v1 * v3
        + 16 * v2**2
        - 30 * c**2 * v1**4
        - 10 * c * v1**2 * v2
    )
    live_d = d_wall_numerator / (30 * v1**3)
    expected_q4 = (
        sp.Rational(64, 45)
        * slope
        / v1
        * (
            30 * c**2 * v1**4
            + 10 * c * v1**2 * v2
            + 30 * d * v1**3
            - 15 * v1**3
            + 18 * v1 * v3
            - 16 * v2**2
        )
    )
    q4_live = sp.factor(q[4].subs(k, live_k))
    require(sp.factor(q4_live - expected_q4) == 0, "live S q4")
    require(
        sp.factor(sp.diff(q4_live, d) - sp.Rational(128, 3) * v1**2 * slope)
        == 0,
        "live d-wall slope",
    )
    exceptional_c = -v2 / (3 * v1**2)
    require(
        sp.factor(live_k.subs(c, exceptional_c) - 2 * exceptional_c) == 0,
        "zero d-wall slope is finite-clutch failure",
    )
    require(sp.factor(expected_q4.subs(d, live_d)) == 0, "unique live d wall")

    symbolic_slope, q5_core, q5_full, q6_core = wall_polynomials(
        c, v1, v2, v3, v4, v5
    )
    require(sp.factor(symbolic_slope - slope) == 0, "wall slope consistency")
    q5_expected = -sp.Rational(32, 315) * q5_full / v1**2
    q6_expected = -sp.Rational(32, 4725) * q6_core / v1**3
    wall_substitution = {k: live_k, d: live_d}
    require(
        sp.factor(q[5].subs(wall_substitution) - q5_expected) == 0,
        "coupled S q5",
    )
    require(
        sp.factor(q[6].subs(wall_substitution) - q6_expected) == 0,
        "coupled S q6",
    )
    require(sp.Poly(q5_full, c).degree() == 3, "q5 numerator degree three")
    require(sp.Poly(q6_core, c).degree() == 4, "q6 numerator degree four")
    print(
        "S_q3=(8/3)*v1^2*(2c-k)*(6c*v1^2+3k*v1^2+4v2);"
        "live_k=-2c-4v2/(3v1^2)"
    )
    print(
        "S_q4_slope=(128/3)*v1^2*(3c*v1^2+v2);"
        "slope_zero=>k=2c=>finite_clutch_failure;"
        "live_d=(15v1^3-18v1v3+16v2^2-30c^2v1^4-10cv1^2v2)/(30v1^3)"
    )
    print("S_post_walls=(deg(q5_num),deg(q6_num))=(3,4)")


def shifted_coefficients(poly, root, field, symbol, count):
    shifted_ring = field.poly_ring(symbol)
    t = shifted_ring.gens[0]
    shifted = shifted_ring.zero
    for (exponent,), coefficient in poly.to_dict().items():
        shifted += coefficient * (t + root) ** exponent
    entries = shifted.to_dict()
    return [entries.get((index,), field.zero) for index in range(count)]


def canonical_poly_text(poly) -> str:
    return "|".join(
        f"{degree}:{coefficient}"
        for (degree,), coefficient in sorted(poly.to_dict().items())
    )


def polynomial_pair_digest(first, second) -> str:
    payload = canonical_poly_text(first.monic()) + "\n" + canonical_poly_text(second.monic())
    return sha256(payload.encode("ascii")).hexdigest()


def euclidean_degree_profile(first, second):
    if first.degree() < second.degree():
        first, second = second, first
    profile = [(first.degree(), second.degree())]
    while second:
        remainder = first.rem(second)
        if not remainder:
            break
        first, second = second, remainder.monic()
        profile.append((first.degree(), second.degree()))
    return tuple(profile), second if second else first


def specialize_factor(K_poly, base_values, polynomial_ring):
    result = polynomial_ring.zero
    for monomial, coefficient in K_poly.terms():
        term = polynomial_ring.convert(coefficient)
        for base, exponent in zip(base_values, monomial):
            term *= base**exponent
        result += term
    return result


def build_case(name: str, K_poly) -> tuple[str, str]:
    u, x, t, c = sp.symbols("u x t c")
    if name == "4111":
        accessory = sp.Poly(
            100 * u**3 + 244 * u**2 + 237 * u + 44, u, domain=QQ
        )
        exponent_a, exponent_b = 4, 1
    else:
        accessory = sp.Poly(
            75 * u**3 - 89 * u**2 - 31 * u + 61, u, domain=QQ
        )
        exponent_a, exponent_b = 3, 2

    field = QQ.alg_field_from_poly(accessory, alias="u")
    alpha = field.ext
    x_ring = field.poly_ring(x)
    X = x_ring.gens[0]
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
    D = X**exponent_a * (X - 1) ** exponent_b * quadratic
    T = X * (X - 1) * quadratic
    E = (
        exponent_a * (X - 1) * quadratic
        + exponent_b * X * quadratic
        + X * (X - 1) * (2 * X - alpha)
    ) / 7
    S = X + shift
    V = 4 * S * D * T**2 / gamma**2
    A = 2 * S * E * T / gamma
    g = S * T
    require(V.degree() == 16 and A.degree() == 8, f"{name} response degrees")
    require(
        2 * V * A.diff(X) - A * V.diff(X) == 2 * V,
        f"{name} response identity",
    )

    V_shift = shifted_coefficients(V, -shift, field, t, 6)
    v1, v2, v3, v4, v5 = V_shift[1:6]
    require(v1 != field.zero, f"{name} simple S root")
    c_ring = field.poly_ring(c)
    local_c = c_ring.gens[0]
    slope, q5_core, q5_full, q6_core = wall_polynomials(
        local_c, v1, v2, v3, v4, v5
    )
    require(
        (slope.degree(), q5_core.degree(), q5_full.degree(), q6_core.degree())
        == (1, 2, 3, 4),
        f"{name} S polynomial degrees",
    )
    require(q5_full.gcd(q6_core).degree() == 0, f"{name} full wall gcd")
    require(q5_core.gcd(q6_core).degree() == 0, f"{name} live wall gcd")
    require(slope.gcd(q6_core).degree() == 0, f"{name} hostile slope gcd")
    full_profile, full_gcd = euclidean_degree_profile(q6_core, q5_full)
    live_profile, live_gcd = euclidean_degree_profile(q6_core, q5_core)
    require(
        full_profile == ((4, 3), (3, 2), (2, 1), (1, 0)),
        f"{name} full PRS profile",
    )
    require(
        live_profile == ((4, 2), (2, 1), (1, 0)),
        f"{name} live PRS profile",
    )
    require(
        full_gcd.degree() == 0 and live_gcd.degree() == 0,
        f"{name} PRS terminal units",
    )
    wall_digest = polynomial_pair_digest(q5_full, q6_core)

    control_C = X + 1
    control_d = field.one
    control_k = field.one
    control_K = specialize_factor(
        K_poly,
        (A, control_C, V, V.diff(X), control_d, control_k),
        x_ring,
    )
    finite_gate = control_k - A.diff(X) * control_C
    require(finite_gate.gcd(g).degree() == 0, f"{name} control finite gate")
    require(control_K.degree() == 97, f"{name} coupled degree 97")
    require(
        control_K.LC == 128 * V.LC**6,
        f"{name} coupled top coefficient",
    )
    boundary = S**3 * T**8 * X ** extras[0] * (X - 1) ** extras[1]
    require(boundary.degree() == 44, f"{name} boundary degree")
    control_H = control_K.exquo(boundary)
    require(control_H.degree() == 53, f"{name} residual degree")
    require(control_H.gcd(g).degree() == 0, f"{name} control boundary disjoint")
    require(
        control_H.gcd(control_H.diff(X)).degree() == 0,
        f"{name} control squarefree residual",
    )
    control_digest = sha256(
        canonical_poly_text(control_H.monic()).encode("ascii")
    ).hexdigest()

    expected_wall, expected_control = EXPECTED_CASE_DIGESTS[name]
    if expected_wall != "TBD":
        require(wall_digest == expected_wall, f"{name} wall digest drift")
    if expected_control != "TBD":
        require(control_digest == expected_control, f"{name} control digest drift")
    print(
        f"case={name};S_degrees=(slope,q5core,q5full,q6)=(1,2,3,4);"
        f"PRS_full={full_profile};PRS_live={live_profile};gcds=(0,0,0);"
        f"wall_digest={wall_digest}"
    )
    print(
        f"case={name};control=(C,Eprime)=(1+x,1);finite_gate=unit_mod_ST;"
        "degrees=(K,boundary,H)=(97,44,53);"
        f"H_squarefree_boundary_disjoint=1;control_digest={control_digest}"
    )
    return wall_digest, control_digest


def source_audit() -> None:
    tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    assert_nodes = sum(isinstance(node, ast.Assert) for node in ast.walk(tree))
    float_literals = sum(
        isinstance(node, ast.Constant) and isinstance(node.value, float)
        for node in ast.walk(tree)
    )
    require(assert_nodes == 0 and float_literals == 0, "source AST gate")
    print(f"source_ast=(assert_nodes={assert_nodes},float_literals={float_literals})")


def main() -> None:
    print("affine transverse C0-E0 coupled clutch critical no-go scout")
    for dependency, expected_hash in DEPENDENCIES.items():
        require(
            lf_hash(ROOT / dependency) == expected_hash,
            f"dependency drift: {dependency}",
        )
    print(f"dependency_hash_checks={len(DEPENDENCIES)}")
    K_poly, symbols = derive_coupled_factor()
    t_boundary_checks(K_poly, symbols)
    s_wall_checks(K_poly, symbols)
    case_digests = tuple(build_case(name, K_poly) for name in ("4111", "3211"))
    print(f"case_digests={case_digests}")
    print(
        "consequence=d*k!=0_and_finite_gate_passes=>"
        "degH=53,ord_S(H)<=3,no_T_escape,off_owner_multiplicity>=50"
    )
    print(
        "scope=all_affine_C0_E0_with_B=1;"
        "d=0_inherited_THM3212;k=0_inherited_THM3279;"
        "not_inverse_cover_not_JC2"
    )
    source_audit()


if __name__ == "__main__":
    main()
