#!/usr/bin/env python3
"""Exact exceptional-quadratic scout after THM-3306.

On the two fixed accessory fields and the slice C=c+x, d=k=1, reconstruct
the degree-36 coefficient of the linear subresultant and its transverse base
ideal.  Pull the preceding quadratic subresultant to the exceptional divisor,
certify the squareclass of its discriminant, and compute the first normal
term in the fixed blow-up chart u=a_1, v=b_1=u*t.

The degree-36 linear-row coefficient is named ``linear_a`` throughout.  The
degree-32 leading coefficient of the quadratic row is named ``P2``.  They are
different polynomials.  This is a finite-exact fixed-slice computation, not a
Keller mate or an inverse-cover construction.
"""

from __future__ import annotations

import ast
from hashlib import sha256
from pathlib import Path

import sympy as sp
from sympy.polys.domains import QQ
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.ddm import DDM


ROOT = Path(__file__).resolve().parents[1]
DEPENDENCIES = {
    "01-canon/theorems/THM-3306-affine-c-critical-section-square-discriminant-and-transverse-base-locus.md":
        "514b5c07a70ea0dd0020857af3278008c3df2a03af6826784d96d059a3b26111",
    "04-computation/jc_affine_c_parameter_section_degeneracy_partial_scout_20260803.py":
        "6448f0a8d8238358adca613610cede3fca72dac210ba0487d126b6d466849697",
    "05-knowledge/results/jc_affine_c_parameter_section_degeneracy_partial_scout_20260803.out":
        "28e0395a4e90be88ebd413f27ebdaf82ee4ba51233d72fca42a423448bd6faa5",
    "04-computation/jc_affine_c_parameter_section_degeneracy_independent_audit_20260803.py":
        "4afadde9302723c5af1a9a209525733c453f928332f21d82e1125f3ddb5662f8",
    "05-knowledge/results/jc_affine_c_parameter_section_degeneracy_independent_audit_20260803.out":
        "d72ed2fb78bd3e237ea5893a75573d54f2c58fa9080799a8633e73547aaba81d",
    "04-computation/jc_critical_inverse_cover_cofactor_jacobian_probe_agent.py":
        "a719b2582b93a0a6d110b1f13b65e9d54800e8669914da9f21a9371545bbae31",
    "05-knowledge/results/jc_critical_inverse_cover_cofactor_jacobian_probe_agent.out":
        "67067d9448caa6a809520b190208a561ab4cc14517455d6da0eef9210ccce1ff",
}
CERTIFICATES = {
    "4111": {
        "prime": 13,
        "accessory_root": 6,
        "norm_residue": 7,
        "factor_characters": ((1, -1), (3, 1), (15, -1), (17, -1)),
        "cstar_digest":
            "2cab20bdd322521e6a4309a2c854cdefcb0022767bc449499bf8a4827136d031",
    },
    "3211": {
        "prime": 23,
        "accessory_root": 16,
        "norm_residue": 22,
        "factor_characters": ((3, -1), (13, 1), (20, 1)),
        "cstar_digest":
            "f1e4a546b447a3f84a96dcc8488d4f32583818a54080329d275610d3fc9c7705",
    },
}
EXPECTED_DIGEST_PACKETS = (
    (
        "6b66fda17d31412a565a23649197ca94e9372c7368f01a31f60da61714b328fc",
        "d1143e9c3304e06d55a6e76edf8cf5bf11ac380af6477f3becf022d3f8cbe945",
        "cf9c3718643b56417170affc6b8498c25ba1487c1e4c35fd76600b8dec42900d",
        "e6b73b773c4c1e54d502aa04663b8d29f1125e2ff95063f259dd070b050e2927",
        "5cdc30990101b93ad8f0e1e3f4c66950bd8a4b6df05cfb8a6f701a29d71000ac",
    ),
    (
        "fa7ce86608c083d3e8e722cf5de425a9a3d3ae3ddbdbc5e829beb6e10fb7187a",
        "1c015dff26e8783ffd72a7e163ec28b8e6ac68e03483721e395638b6bc84eb96",
        "b7130c3d5db238b35857e95a79d93833605a194a219e32d5be65fd424fb42209",
        "74432a9c7f78976c9d56791644f6e088f55ea3ee7043e91624872af279155437",
        "3248875cd752b42e2154811bdd2fc1c2d0a1ddff89f62ace4ed703245955bbd3",
    ),
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def lf_hash(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def element_text(element) -> str:
    return "|".join(
        f"{monomial}:{coefficient}"
        for monomial, coefficient in sorted(element.to_dict().items())
    )


def element_digest(element) -> str:
    return sha256(element_text(element).encode("utf-8")).hexdigest()


def packet_digest(elements) -> str:
    text = "\n".join(element_text(element) for element in elements)
    return sha256(text.encode("utf-8")).hexdigest()


def specialize_sparse(poly, values, polynomial_ring):
    result = polynomial_ring.zero
    for monomial, coefficient in poly.terms():
        term = polynomial_ring.convert(coefficient)
        for value, exponent in zip(values, monomial):
            term *= value**exponent
        result += term
    return result


def factor_profile(poly) -> tuple[tuple[int, int], ...]:
    _, factors = poly.factor_list()
    return tuple((factor.degree(), exponent) for factor, exponent in factors)


def derive_generic_rows():
    """Derive the PRS row and the closed-form row by independent formulas."""

    y, A, C, V, derivative_V, d, k = sp.symbols("y A C V Vp d k")
    L = y**2 + y + C * V
    first = sp.expand(2 * L * (2 * y + 1) + V * A)
    second = sp.expand(
        V**3 * k + V**2 * y
        + L * (-derivative_V * y + 2 * V**2 * d)
    )
    sequence = sp.subresultants(first, second, y)
    require(
        tuple(sp.Poly(row, y).degree() for row in sequence) == (3, 3, 2, 1, 0),
        "generic subresultant profile",
    )
    quadratic = sp.Poly(sequence[-3], y)
    linear = sp.Poly(sequence[-2], y)

    # Closed-form cubic cancellation from the inverse-graph sidecar.
    manual_quadratic = sp.Poly(sp.expand(4 * second + derivative_V * first), y)
    P2 = manual_quadratic.nth(2)
    Q2 = manual_quadratic.nth(1)
    R2 = manual_quadratic.nth(0)
    D = sp.expand(6 * P2 - 4 * Q2)
    first_poly = sp.Poly(first, y)
    ell_1 = sp.expand(P2 * (P2 * first_poly.nth(1) - 4 * R2) - D * Q2)
    ell_0 = sp.expand(P2**2 * first_poly.nth(0) - D * R2)

    require(
        all(
            sp.expand(quadratic.nth(index) - manual_quadratic.nth(index)) == 0
            for index in range(3)
        ),
        "closed-form and PRS quadratic rows agree",
    )
    require(
        sp.expand(4 * linear.nth(1) - ell_1) == 0
        and sp.expand(4 * linear.nth(0) - ell_0) == 0,
        "closed-form and PRS linear rows agree",
    )
    chart_identity = sp.expand(
        P2 * ell_0**2 - Q2 * ell_0 * ell_1 + R2 * ell_1**2
        - 16 * P2**2 * sp.Poly(sequence[-1], y).nth(0)
    )
    require(chart_identity == 0, "quadratic exceptional-chart identity")
    require(sp.Poly(first, y).LC() == 4, "no common projective y-root")

    symbols = (A, C, V, derivative_V, d, k)
    print(
        "generic_subresultants=(3,3,2,1,0);"
        "quadratic_row=4*R2+Vprime*R1;linear_row=(ell1*y+ell0)/4;"
        "P2*ell0^2-Q2*ell0*ell1+R2*ell1^2=16*P2^2*resultant;"
        "R1_y_lead=4"
    )
    return {
        "symbols": symbols,
        "linear_a": sp.Poly(linear.nth(1), *symbols),
        "linear_b": sp.Poly(linear.nth(0), *symbols),
        "quadratic": tuple(
            sp.Poly(quadratic.nth(index), *symbols) for index in (2, 1, 0)
        ),
    }


def response_case(name: str):
    u, x = sp.symbols("u x")
    if name == "4111":
        accessory = sp.Poly(
            100 * u**3 + 244 * u**2 + 237 * u + 44, u, domain=QQ
        )
        exponent_a, exponent_b = 4, 1
    elif name == "3211":
        accessory = sp.Poly(
            75 * u**3 - 89 * u**2 - 31 * u + 61, u, domain=QQ
        )
        exponent_a, exponent_b = 3, 2
    else:
        raise RuntimeError(f"unknown accessory case: {name}")

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
    D_source = X**exponent_a * (X - 1) ** exponent_b * quadratic
    T = X * (X - 1) * quadratic
    E_source = (
        exponent_a * (X - 1) * quadratic
        + exponent_b * X * quadratic
        + X * (X - 1) * (2 * X - alpha)
    ) / 7
    S = X + shift
    V = 4 * S * D_source * T**2 / gamma**2
    A = 2 * S * E_source * T / gamma
    boundary = S**3 * T**8 * X**extras[0] * (X - 1) ** extras[1]
    require((V.degree(), A.degree(), boundary.degree()) == (16, 8, 44), name)
    require(2 * V * A.diff(X) - A * V.diff(X) == 2 * V, f"{name} response")
    return {
        "name": name,
        "accessory": accessory,
        "field": field,
        "x_symbol": x,
        "ring": x_ring,
        "X": X,
        "V": V,
        "A": A,
        "boundary": boundary,
    }


def inverse_mod(element, modulus):
    inverse, _, common = element.gcdex(modulus)
    require(common.degree() == 0, "requested quotient element is a unit")
    return (inverse / common.LC) % modulus


def multiplication_charpoly(element, modulus, field, symbol):
    degree = modulus.degree()
    X = modulus.ring.gens[0]
    columns = []
    for exponent in range(degree):
        residue = (element * X**exponent) % modulus
        columns.append(
            [residue.get((row,), field.zero) for row in range(degree)]
        )
    rows = [
        [columns[column][row] for column in range(degree)]
        for row in range(degree)
    ]
    matrix = DomainMatrix.from_rep(DDM(rows, (degree, degree), field))
    coefficients = matrix.charpoly()
    ring = field.poly_ring(symbol)
    Z = ring.gens[0]
    characteristic = ring.zero
    for index, coefficient in enumerate(coefficients):
        characteristic += coefficient * Z ** (degree - index)
    return characteristic


def rational_mod(value, prime: int) -> int:
    numerator = int(value.numerator)
    denominator = int(value.denominator)
    require(denominator % prime != 0, "bad rational denominator")
    return numerator * pow(denominator, -1, prime) % prime


def algebraic_mod(value, prime: int, root: int) -> int:
    result = 0
    for coefficient in value.to_list():
        result = (result * root + rational_mod(coefficient, prime)) % prime
    return result


def coefficients_good(poly, prime: int) -> bool:
    return all(
        int(rational.denominator) % prime != 0
        for coefficient in poly.to_dict().values()
        for rational in coefficient.to_list()
    )


def reduce_to_prime(poly, prime: int, root: int, symbol):
    finite_ring = sp.GF(prime).poly_ring(symbol)
    X = finite_ring.gens[0]
    result = finite_ring.zero
    for (degree,), coefficient in poly.to_dict().items():
        result += finite_ring.convert(algebraic_mod(coefficient, prime, root)) * X**degree
    return result


def quotient_power(element, exponent: int, modulus):
    result = modulus.ring.one
    base = element % modulus
    while exponent:
        if exponent & 1:
            result = (result * base) % modulus
        base = (base * base) % modulus
        exponent //= 2
    return result


def certify_nonsquare(name, accessory, linear_a, discriminant, norm):
    certificate = CERTIFICATES[name]
    prime = certificate["prime"]
    root = certificate["accessory_root"]
    accessory_derivative = accessory.diff()
    accessory_value = sum(
        rational_mod(coefficient, prime) * pow(root, degree, prime)
        for (degree,), coefficient in accessory.as_dict().items()
    ) % prime
    derivative_value = sum(
        rational_mod(coefficient, prime) * pow(root, degree, prime)
        for (degree,), coefficient in accessory_derivative.as_dict().items()
    ) % prime
    require(accessory_value == 0 and derivative_value != 0, f"{name} good K-prime")
    require(
        coefficients_good(linear_a, prime)
        and coefficients_good(discriminant, prime),
        f"{name} good quotient reduction",
    )
    reduced_a = reduce_to_prime(linear_a, prime, root, sp.Symbol("x"))
    reduced_delta = reduce_to_prime(discriminant, prime, root, sp.Symbol("x"))
    finite_X = reduced_a.ring.gens[0]
    require(
        reduced_a.degree() == 36
        and reduced_a.gcd(reduced_a.diff(finite_X)).degree() == 0,
        f"{name} etale quotient reduction",
    )
    characters = []
    for factor, exponent in reduced_a.factor_list()[1]:
        require(exponent == 1, f"{name} squarefree factor")
        residue = reduced_delta % factor
        require(residue != factor.ring.zero, f"{name} discriminant unit mod factor")
        character = quotient_power(
            residue, (prime**factor.degree() - 1) // 2, factor
        )
        if character == factor.ring.one:
            sign = 1
        elif character == -factor.ring.one:
            sign = -1
        else:
            raise RuntimeError(f"{name} invalid quadratic character")
        characters.append((factor.degree(), sign))
    characters = tuple(characters)
    require(characters == certificate["factor_characters"], f"{name} factor characters")

    norm_residue = algebraic_mod(norm, prime, root)
    require(norm_residue == certificate["norm_residue"], f"{name} norm residue")
    require(
        pow(norm_residue, (prime - 1) // 2, prime) == prime - 1,
        f"{name} norm is a residue non-square",
    )
    require(
        sp.prod(sign for _, sign in characters) == -1,
        f"{name} factor characters recover norm character",
    )
    return prime, root, norm_residue, characters


def audit_case(case, rows):
    name = case["name"]
    accessory = case["accessory"]
    field = case["field"]
    ring = case["ring"]
    X = case["X"]
    V = case["V"]
    A = case["A"]
    boundary = case["boundary"]
    derivative_V = V.diff(X)
    second_derivative_V = derivative_V.diff(X)

    def specialize_linear(parameter):
        values = (
            A, X + field.convert(parameter), V, derivative_V,
            field.one, field.one,
        )
        return (
            specialize_sparse(rows["linear_a"], values, ring).exquo(boundary),
            specialize_sparse(rows["linear_b"], values, ring).exquo(boundary),
        )

    raw_a_0, raw_b_0 = specialize_linear(0)
    raw_a_1, raw_b_1 = specialize_linear(1)
    require(raw_a_0 == raw_a_1, f"{name} linear coefficient independent of c")
    common_unit = raw_a_0.LC
    linear_a = raw_a_0 / common_unit
    b0 = raw_b_0 / common_unit
    b1 = (raw_b_1 - raw_b_0) / common_unit
    require(linear_a == linear_a.monic(), f"{name} monic linear coefficient")
    require(
        (linear_a.degree(), b0.degree(), b1.degree()) == (36, 44, 36),
        f"{name} linear row degrees",
    )
    require(factor_profile(linear_a) == ((36, 1),), f"{name} irreducible a")
    require(linear_a.gcd(b1).degree() == 0, f"{name} b1 unit")
    cstar = (-b0 * inverse_mod(b1, linear_a)) % linear_a
    require(cstar.degree() == 24, f"{name} cstar degree")
    require((b0 + cstar * b1) % linear_a == ring.zero, f"{name} base ideal")
    require(
        element_digest(cstar) == CERTIFICATES[name]["cstar_digest"],
        f"{name} cstar digest drift",
    )

    Cstar = X + cstar
    P2_raw = 2 * derivative_V + 8 * V**2
    Q2_raw = 2 * derivative_V + 12 * V**2
    R2_raw = V * (
        4 * V**2 + 8 * Cstar * V**2 + derivative_V * (2 * Cstar + A)
    )
    require(
        (P2_raw.degree(), Q2_raw.degree()) == (32, 32),
        f"{name} degree-32 quadratic pivot",
    )
    P2, Q2, R2 = (P2_raw % linear_a, Q2_raw % linear_a, R2_raw % linear_a)
    values_at_base = (
        A, Cstar, V, derivative_V, field.one, field.one,
    )
    prs_quadratic = tuple(
        specialize_sparse(coefficient, values_at_base, ring) % linear_a
        for coefficient in rows["quadratic"]
    )
    require(prs_quadratic == (P2, Q2, R2), f"{name} independent row agreement")
    require(
        all(linear_a.gcd(coefficient).degree() == 0 for coefficient in (P2, Q2, R2)),
        f"{name} quadratic coefficients are units",
    )
    discriminant = (Q2**2 - 4 * P2 * R2) % linear_a
    require(linear_a.gcd(discriminant).degree() == 0, f"{name} separable quadratic")

    z = sp.Symbol("z")
    discriminant_charpoly = multiplication_charpoly(
        discriminant, linear_a, field, z
    )
    require(
        discriminant_charpoly.degree() == 36
        and factor_profile(discriminant_charpoly) == ((36, 1),),
        f"{name} discriminant is primitive over accessory field",
    )
    norm = linear_a.resultant(discriminant)
    require(
        discriminant_charpoly.to_dict()[(0,)] == norm,
        f"{name} characteristic constant equals norm",
    )
    modular_certificate = certify_nonsquare(
        name, accessory, linear_a, discriminant, norm
    )

    # Blow up the fixed normalized base coordinates
    #   u=linear_a(x), v=b0(x)+c*b1(x)=u*t.
    # The transverse Jacobian makes x,c formal functions of u,t.  The first
    # derivatives at u=0 are dx/du=1/a' and
    # dc/du=(t-b_x/a')/b1, with x-derivatives taken at fixed c.
    a_prime = linear_a.diff(X) % linear_a
    b_x = (b0.diff(X) + cstar * b1.diff(X)) % linear_a
    inverse_a_prime = inverse_mod(a_prime, linear_a)
    inverse_b1 = inverse_mod(b1, linear_a)
    require((a_prime * inverse_a_prime) % linear_a == ring.one, name)
    require((b1 * inverse_b1) % linear_a == ring.one, name)

    P2_x = (2 * second_derivative_V + 16 * V * derivative_V) % linear_a
    Q2_x = (2 * second_derivative_V + 24 * V * derivative_V) % linear_a
    bracket = (
        4 * V**2 + 8 * Cstar * V**2
        + derivative_V * (2 * Cstar + A)
    )
    R2_x = (
        derivative_V * bracket
        + V * (
            8 * V * derivative_V + 8 * V**2
            + 16 * Cstar * V * derivative_V
            + second_derivative_V * (2 * Cstar + A)
            + derivative_V * (2 + A.diff(X))
        )
    ) % linear_a
    R2_c = (V * (8 * V**2 + 2 * derivative_V)) % linear_a

    F1_2 = (P2_x * inverse_a_prime) % linear_a
    F1_1 = (-Q2_x * inverse_a_prime + R2_c * inverse_b1) % linear_a
    F1_0 = (
        R2_x * inverse_a_prime
        - R2_c * b_x * inverse_a_prime * inverse_b1
    ) % linear_a
    require(
        all(coefficient != ring.zero for coefficient in (F1_2, F1_1, F1_0)),
        f"{name} nonzero first normal packet",
    )

    # Reduce F1 modulo F0=P2*t^2-Q2*t+R2.  If the remainder vanished, the
    # exceptional pair would be stationary to first order.
    leading_ratio = (F1_2 * inverse_mod(P2, linear_a)) % linear_a
    remainder_1 = (F1_1 + leading_ratio * Q2) % linear_a
    remainder_0 = (F1_0 - leading_ratio * R2) % linear_a
    common_root_core = (
        P2 * remainder_0**2
        + Q2 * remainder_0 * remainder_1
        + R2 * remainder_1**2
    ) % linear_a
    require(common_root_core != ring.zero, f"{name} first normal is coprime to F0")

    inverse_discriminant = inverse_mod(discriminant, linear_a)
    velocity_1 = (
        -(remainder_1 * Q2 + 2 * P2 * remainder_0)
        * inverse_discriminant
    ) % linear_a
    velocity_0 = (
        (2 * remainder_1 * R2 + remainder_0 * Q2)
        * inverse_discriminant
    ) % linear_a
    require(
        velocity_1 != ring.zero and velocity_0 != ring.zero,
        f"{name} distinct nonzero conjugate velocities",
    )

    print(
        f"case={name};linear_pivot_degree=36;"
        "quadratic_lead_P2_degree_before_mod=32;symbols_distinct=1;"
        "quadratic_row_agrees_by=(generic_PRS,direct_4R2+VprimeR1)"
    )
    print(
        f"case={name};exceptional_F0=P2*t^2-Q2*t+R2;"
        f"coefficient_degrees={tuple(c.degree() for c in (P2, Q2, R2))};"
        "P2_Q2_R2_units=1;t=0_not_root=1;t=infinity_not_root=1;"
        "projective_y_infinity_not_common=1"
    )
    print(
        f"case={name};delta_degree={discriminant.degree()};"
        f"delta_digest={element_digest(discriminant)};"
        "delta_charpoly_factor_profile=((36,1),);delta_primitive_degree=36;"
        f"norm_certificate=(p={modular_certificate[0]},alpha={modular_certificate[1]},"
        f"norm_mod_p={modular_certificate[2]},Legendre=-1);"
        f"a_mod_p_factor_characters={modular_certificate[3]};"
        "delta_nonsquare_in_A=1;F0_irreducible_separable_over_A=1"
    )
    print(
        f"case={name};exceptional_intersection=Spec(A[sqrt(delta)]);"
        "relative_degree_over_A=2;degree_over_K=72;"
        "relative_geometric_fibre_points=2;total_geometric_points_over_K=72;"
        "fibrewise_deck_exchange=C2;"
        "no_A_rational_exceptional_direction=1"
    )
    print(
        f"case={name};first_normal_F1_degrees="
        f"{tuple(c.degree() for c in (F1_2, F1_1, F1_0))};"
        f"F1_digest={packet_digest((F1_2, F1_1, F1_0))};"
        f"F1_mod_F0_degrees={(remainder_1.degree(), remainder_0.degree())};"
        f"remainder_digest={packet_digest((remainder_1, remainder_0))};"
        "gcd(F0,F1)=1;neither_geometric_branch_stationary=1"
    )
    print(
        f"case={name};branch_velocity=-F1/F0prime="
        "velocity_1*t+velocity_0;"
        f"velocity_degrees={(velocity_1.degree(), velocity_0.degree())};"
        f"velocity_digest={packet_digest((velocity_1, velocity_0))};"
        "velocity_1_nonzero=1;conjugate_velocities_distinct=1"
    )
    return (
        element_digest(discriminant),
        element_digest(discriminant_charpoly),
        packet_digest((F1_2, F1_1, F1_0)),
        packet_digest((remainder_1, remainder_0)),
        packet_digest((velocity_1, velocity_0)),
    )


def source_audit() -> None:
    source = Path(__file__).read_text(encoding="utf-8")
    tree = ast.parse(source)
    assert_nodes = sum(isinstance(node, ast.Assert) for node in ast.walk(tree))
    float_literals = sum(
        isinstance(node, ast.Constant) and isinstance(node.value, float)
        for node in ast.walk(tree)
    )
    imported_roots = {
        alias.name.split(".")[0]
        for node in ast.walk(tree)
        if isinstance(node, ast.Import)
        for alias in node.names
    }
    imported_roots.update(
        node.module.split(".")[0]
        for node in ast.walk(tree)
        if isinstance(node, ast.ImportFrom) and node.module is not None
    )
    require(assert_nodes == 0 and float_literals == 0, "source AST gate")
    require(
        not imported_roots.intersection({"importlib", "runpy", "subprocess"}),
        "no imported execution path",
    )
    print(
        f"source_ast=(assert_nodes={assert_nodes},float_literals={float_literals},"
        "imported_execution_paths=0)"
    )


def main() -> None:
    print("finite-exact affine-c exceptional-quadratic blow-up scout")
    for relative_path, expected_hash in DEPENDENCIES.items():
        require(
            lf_hash(ROOT / relative_path) == expected_hash,
            f"dependency drift: {relative_path}",
        )
    print(f"dependency_hash_checks={len(DEPENDENCIES)}")
    rows = derive_generic_rows()
    digests = tuple(
        audit_case(response_case(name), rows) for name in ("4111", "3211")
    )
    require(digests == EXPECTED_DIGEST_PACKETS, "case digest packet drift")
    print(f"case_digest_packets={digests}")
    print(
        "scope=FINITE-EXACT_FIXED_SLICE_EXCEPTIONAL_QUADRATIC_AND_FIRST_NORMAL;"
        "fixed_C=c+x_d=k=1_two_accessory_fields_only;"
        "algebraic_C2_deck_exchange_not_global_family_monodromy;"
        "no_root_choice_no_Keller_mate_no_inverse_cover_not_JC2_not_DC2"
    )
    source_audit()


if __name__ == "__main__":
    main()
