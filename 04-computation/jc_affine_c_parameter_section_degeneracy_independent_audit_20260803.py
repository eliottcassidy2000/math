#!/usr/bin/env python3
"""Independent exact audit of the affine-C critical-section degeneracy.

The audit rebuilds the two THM-3212 response pairs and the localized gradient
subresultants from their defining formulas.  It does not import or execute the
parameter scout being audited.  Its main route is quotient algebra: after the
linear subresultant is written as ``a(x)y+b_0(x)+c b_1(x)``, all calculations
are performed in the degree-36 field K[x]/(a).  This gives the degeneracy
polynomial as a multiplication characteristic polynomial and proves the two
resultant identities by exact norm factorization, without asking a bivariate
degree-36 Sylvester determinant to discover them.

This is a FINITE-EXACT PARTIAL audit of one explicit source family.  It does
not construct a Keller mate or an inverse cover and proves no case of JC(2).
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
    "04-computation/jc_centered_heptic_source_morse_obstruction_thm3212.py":
        "03d121e57dd2edaece7cd8693d792349f03757c6e781eb5d9d0c897fcc889448",
    "05-knowledge/results/jc_centered_heptic_source_morse_obstruction_thm3212.out":
        "729e0c7b9fa51fa5c4ac5f18f50dc4413c8a8bb7bf5f0ebf2a7709650304bc85",
    "04-computation/jc_affine_transverse_c0_e0_coupled_clutch_no_go_thm3289.py":
        "f63ff06e3f5ed30f3f6bc5be99756c347d6af5f8e9b220ce8336abff2cd2ca31",
    "05-knowledge/results/jc_affine_transverse_c0_e0_coupled_clutch_no_go_thm3289.out":
        "1aef4341650cdfaf1043a8699e3a1725a0100af6d9848d99dfa924b6f054dba1",
    "04-computation/jc_affine_transverse_c0_e0_coupled_clutch_no_go_thm3289_independent_audit.py":
        "b2fa8c96854549ccb9e515485214c119b685b31456fb7c53c5e2bd83f7933831",
    "05-knowledge/results/jc_affine_transverse_c0_e0_coupled_clutch_no_go_thm3289_independent_audit.out":
        "48d50289c98c9dd17e099497d21cd9648cac27a097b339a1f75a1e13d8fd8837",
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def lf_hash(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def canonical_poly_text(poly) -> str:
    return "|".join(
        f"{monomial[0]}:{coefficient}"
        for monomial, coefficient in sorted(poly.monic().to_dict().items())
    )


def poly_digest(poly) -> str:
    return sha256(canonical_poly_text(poly).encode("utf-8")).hexdigest()


def specialize_sparse(poly, base_values, polynomial_ring):
    result = polynomial_ring.zero
    for monomial, coefficient in poly.terms():
        term = polynomial_ring.convert(coefficient)
        for value, exponent in zip(base_values, monomial):
            term *= value**exponent
        result += term
    return result


def derive_generic_packet():
    """Derive all elimination rows directly from the two localized cubics."""

    y, A, C, V, derivative_V, d, k = sp.symbols("y A C V Vp d k")
    inner = y**2 + y + C * V
    first = sp.expand(2 * inner * (2 * y + 1) + V * A)
    second = sp.expand(
        V**3 * k + V**2 * y
        + inner * (-derivative_V * y + 2 * V**2 * d)
    )
    sequence = sp.subresultants(first, second, y)
    profile = tuple(sp.Poly(row, y).degree() for row in sequence)
    require(profile == (3, 3, 2, 1, 0), "generic subresultant profile")
    require(sp.Poly(first, y).LC() == 4, "first cubic leading coefficient")

    resultant = sp.Poly(sequence[-1], y).nth(0)
    quotient = sp.cancel(resultant / V**3)
    require(sp.denom(quotient) == 1, "resultant has exact V^3 factor")
    K_poly = sp.Poly(quotient, A, C, V, derivative_V, d, k)
    linear = sp.Poly(sequence[-2], y)
    linear_a = sp.Poly(linear.nth(1), A, C, V, derivative_V, d, k)
    linear_b = sp.Poly(linear.nth(0), A, C, V, derivative_V, d, k)
    quadratic = sp.Poly(sequence[-3], y)
    quadratic_coefficients = tuple(
        sp.Poly(quadratic.nth(index), A, C, V, derivative_V, d, k)
        for index in range(3)
    )
    require(len(K_poly.terms()) == 40, "localized quotient has 40 terms")
    require(K_poly.degree(C) == 2, "localized quotient is quadratic in C")
    require(linear_a.degree(C) == 0, "linear coefficient is independent of C")
    require(linear_b.degree(C) == 1, "constant coefficient is affine in C")
    print(
        "generic_subresultants=(3,3,2,1,0);Res_y=V^3*K40;"
        "degrees_in_C=(K,a,b)=(2,0,1);R1_y_lead=4"
    )
    return {
        "K": K_poly,
        "a": linear_a,
        "b": linear_b,
        "quadratic": quadratic_coefficients,
        "symbols": (A, C, V, derivative_V, d, k),
    }


def multiplication_matrix(element, modulus, field, X):
    degree = modulus.degree()
    columns = []
    for column in range(degree):
        residue = (element * X**column) % modulus
        columns.append(
            [residue.get((row,), field.zero) for row in range(degree)]
        )
    rows = [
        [columns[column][row] for column in range(degree)]
        for row in range(degree)
    ]
    return DomainMatrix.from_rep(DDM(rows, (degree, degree), field))


def krylov_determinant(element, modulus, field, X):
    degree = modulus.degree()
    powers = []
    current = modulus.ring.one
    for _ in range(degree):
        powers.append(
            [current.get((row,), field.zero) for row in range(degree)]
        )
        current = (current * element) % modulus
    rows = [
        [powers[column][row] for column in range(degree)]
        for row in range(degree)
    ]
    matrix = DomainMatrix.from_rep(DDM(rows, (degree, degree), field))
    return matrix.det()


def response_case(name: str):
    """Rebuild one THM-3212 cubic response pair from the theorem formulas."""

    u, x = sp.symbols("u x")
    if name == "4111":
        modulus = sp.Poly(
            100 * u**3 + 244 * u**2 + 237 * u + 44, u, domain=QQ
        )
        exponent_a, exponent_b = 4, 1
    elif name == "3211":
        modulus = sp.Poly(
            75 * u**3 - 89 * u**2 - 31 * u + 61, u, domain=QQ
        )
        exponent_a, exponent_b = 3, 2
    else:
        raise RuntimeError(f"unknown case {name}")

    field = QQ.alg_field_from_poly(modulus, alias="u")
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
    E = (
        exponent_a * (X - 1) * quadratic
        + exponent_b * X * quadratic
        + X * (X - 1) * (2 * X - alpha)
    ) / 7
    S = X + shift
    V = 4 * S * D_source * T**2 / gamma**2
    A = 2 * S * E * T / gamma
    owner = S * T
    boundary = (
        S**3 * T**8
        * X**extras[0] * (X - 1) ** extras[1]
    )
    require((V.degree(), A.degree()) == (16, 8), f"{name} response degrees")
    require(boundary.degree() == 44, f"{name} boundary degree")
    require(
        2 * V * A.diff(X) - A * V.diff(X) == 2 * V,
        f"{name} response identity",
    )
    radical = V.exquo(V.gcd(V.diff(X))).monic()
    require(radical == owner.monic(), f"{name} radical of V")
    return {
        "name": name,
        "field": field,
        "x_symbol": x,
        "x_ring": x_ring,
        "X": X,
        "shift": field.convert(shift),
        "S": S,
        "T": T,
        "owner": owner,
        "boundary": boundary,
        "V": V,
        "A": A,
    }


def factor_profile(poly):
    coefficient, factors = poly.factor_list()
    require(coefficient != 0, "factorization leading unit")
    require(all(exponent == 1 for _, exponent in factors), "squarefree factors")
    return tuple(sorted(factor.degree() for factor, _ in factors))


def audit_case(case, packet):
    name = case["name"]
    field = case["field"]
    x_symbol = case["x_symbol"]
    x_ring = case["x_ring"]
    X = case["X"]
    shift = case["shift"]
    V = case["V"]
    A = case["A"]
    owner = case["owner"]
    boundary = case["boundary"]
    derivative_V = V.diff(X)
    generic_values = packet["symbols"]

    def specialize_at(parameter):
        values = (
            A,
            X + field.convert(parameter),
            V,
            derivative_V,
            field.one,
            field.one,
        )
        K_value = specialize_sparse(packet["K"], values, x_ring)
        raw_a = specialize_sparse(packet["a"], values, x_ring)
        raw_b = specialize_sparse(packet["b"], values, x_ring)
        return (
            K_value.exquo(boundary),
            raw_a.exquo(boundary),
            raw_b.exquo(boundary),
        )

    H_at_0, a_at_0, b_at_0 = specialize_at(0)
    H_at_1, a_at_1, b_at_1 = specialize_at(1)
    H_at_2, a_at_2, b_at_2 = specialize_at(2)
    H_at_3, a_at_3, b_at_3 = specialize_at(3)
    H2 = (H_at_2 - 2 * H_at_1 + H_at_0) / 2
    H1 = H_at_1 - H_at_0 - H2
    H0 = H_at_0
    common_unit = a_at_0.LC
    a = a_at_0 * (field.one / common_unit)
    b0 = b_at_0 * (field.one / common_unit)
    b1 = (b_at_1 - b_at_0) * (field.one / common_unit)
    require(a == a.monic(), f"{name} monic a")
    require(a_at_0 == a_at_1 == a_at_2 == a_at_3, f"{name} a independent c")
    require(b_at_2 == b_at_0 + 2 * (b_at_1 - b_at_0), f"{name} b affine c")
    require(b_at_3 == b_at_0 + 3 * (b_at_1 - b_at_0), f"{name} b affine check")
    require(H_at_3 == H0 + 3 * H1 + 9 * H2, f"{name} H quadratic c")
    require(
        (H0.degree(), H1.degree(), H2.degree(), a.degree(), b0.degree(), b1.degree())
        == (53, 52, 36, 36, 44, 36),
        f"{name} bidegrees",
    )
    require(factor_profile(a) == (36,), f"{name} a irreducible")
    require(a.gcd(a.diff(X)).degree() == 0, f"{name} a squarefree")
    require(a.gcd(b1).degree() == 0, f"{name} b1 unit modulo a")

    inverse_b1, _, gcd_b1 = b1.gcdex(a)
    inverse_b1 *= field.one / gcd_b1.LC
    cstar = (-b0 * inverse_b1) % a
    require(cstar.degree() == 24, f"{name} cstar degree")
    require((b0 + b1 * cstar) % a == x_ring.zero, f"{name} b graph")
    require(H2.gcd(a).degree() == 0, f"{name} H2 unit modulo a")
    require((H1 + 2 * H2 * cstar) % a == x_ring.zero, f"{name} H linear identity")
    require((H0 - H2 * cstar**2) % a == x_ring.zero, f"{name} H constant identity")

    multiplication = multiplication_matrix(cstar, a, field, X)
    characteristic = multiplication.charpoly()
    require(len(characteristic) == 37, f"{name} characteristic degree")
    require(
        krylov_determinant(cstar, a, field, X) != field.zero,
        f"{name} cstar is primitive",
    )
    parameter = sp.Symbol("c")
    parameter_ring = field.poly_ring(parameter)
    parameter_generator = parameter_ring.gens[0]
    D = parameter_ring.zero
    for index, coefficient in enumerate(characteristic):
        D += coefficient * parameter_generator ** (36 - index)
    require(D.degree() == 36 and D == D.monic(), f"{name} monic D")
    require(factor_profile(D) == (36,), f"{name} D irreducible")
    require(D.gcd(D.diff(parameter_generator)).degree() == 0, f"{name} D squarefree")

    # These quotient identities are the norm proof of both resultant claims:
    # b mod a = b1(c-cstar) and H mod a = H2(c-cstar)^2.  Since a is monic
    # irreducible and b1,H2 are units, the omitted norms are nonzero field
    # constants.  The characteristic polynomial is exactly Norm(c-cstar).
    jacobian_residue = (a.diff(X) * b1) % a
    require(jacobian_residue != x_ring.zero, f"{name} base Jacobian nonzero")
    require(a.gcd(a.diff(X) * b1).degree() == 0, f"{name} transverse base ideal")

    quadratic_values = (
        A,
        X + cstar,
        V,
        derivative_V,
        field.one,
        field.one,
    )
    quadratic_residues = tuple(
        specialize_sparse(row, quadratic_values, x_ring) % a
        for row in packet["quadratic"]
    )
    require(
        all(residue != x_ring.zero for residue in quadratic_residues),
        f"{name} quadratic subresultant survives",
    )
    require(a.gcd(derivative_V).degree() == 0, f"{name} second cubic stays cubic")

    H_x_at_incidence = (
        H0.diff(X) + H1.diff(X) * cstar + H2.diff(X) * cstar**2
    ) % a
    H_xx_at_incidence = (
        H0.diff(X).diff(X)
        + H1.diff(X).diff(X) * cstar
        + H2.diff(X).diff(X) * cstar**2
    ) % a
    require(H_x_at_incidence == x_ring.zero, f"{name} first H derivative")
    require(H_xx_at_incidence != x_ring.zero, f"{name} second H derivative")
    require(a.gcd(H_xx_at_incidence).degree() == 0, f"{name} exact H multiplicity two")

    c_over_x_ring = parameter_ring.poly_ring(x_symbol)
    X_parameter = c_over_x_ring.gens[0]
    c_as_x_constant = c_over_x_ring({(0,): parameter_generator})
    H_bivariate = (
        c_over_x_ring.convert(H0)
        + c_over_x_ring.convert(H1) * c_as_x_constant
        + c_over_x_ring.convert(H2) * c_as_x_constant**2
    )
    owner_bivariate = c_over_x_ring.convert(owner)
    finite_gate_bivariate = (
        c_over_x_ring.one
        - c_over_x_ring.convert(A.diff(X))
        * (X_parameter + c_as_x_constant)
    )
    finite_wall = owner_bivariate.resultant(finite_gate_bivariate).monic()
    boundary_wall = owner_bivariate.resultant(H_bivariate).monic()
    require((finite_wall.degree(), boundary_wall.degree()) == (5, 6), f"{name} wall degrees")
    require(factor_profile(finite_wall) == (1, 1, 1, 2), f"{name} finite wall profile")
    require(
        factor_profile(boundary_wall) == (1, 1, 1, 1, 2),
        f"{name} boundary wall profile",
    )
    live_wall = boundary_wall.exquo(finite_wall)
    require(live_wall.degree() == 1, f"{name} live wall degree")
    require(boundary_wall == finite_wall * live_wall, f"{name} wall factorization")

    s_root = -shift
    v1 = V.diff(X).evaluate(X, s_root)
    v2 = V.diff(X).diff(X).evaluate(X, s_root) / 2
    half = field.convert(sp.Rational(1, 2))
    finite_s_parameter = shift + half
    live_s_parameter = shift - half - 2 * v2 / (3 * v1**2)
    exceptional_parameter = shift - v2 / (3 * v1**2)
    require(
        exceptional_parameter == (finite_s_parameter + live_s_parameter) / 2,
        f"{name} exceptional midpoint",
    )
    require(
        finite_wall.evaluate(parameter_generator, finite_s_parameter) == field.zero,
        f"{name} finite S wall",
    )
    require(
        live_wall.evaluate(parameter_generator, live_s_parameter) == field.zero,
        f"{name} live q3 wall",
    )
    require(
        finite_wall.evaluate(parameter_generator, exceptional_parameter) != field.zero
        and boundary_wall.evaluate(parameter_generator, exceptional_parameter) != field.zero,
        f"{name} exceptional point off finite and boundary walls",
    )
    require(D.gcd(boundary_wall * finite_wall).degree() == 0, f"{name} D wall disjointness")
    require(
        all(
            D.evaluate(parameter_generator, value) != field.zero
            for value in (
                finite_s_parameter,
                live_s_parameter,
                exceptional_parameter,
            )
        ),
        f"{name} D avoids three named S walls",
    )

    print(
        f"case={name};H_bidegree=(c:2,x:53);a_degree=36;"
        "a_irreducible_squarefree=1;b=b0+c*b1;b1_unit_mod_a=1;"
        "cstar_degree=24;H_mod_a=H2*(c-cstar)^2;H2_unit_mod_a=1"
    )
    print(
        f"case={name};D_degree=36;D_irreducible_squarefree=1;"
        "Krylov_rank=36;Res_x(a,b)=field_unit*D;"
        "Res_x(a,H)=field_unit*D^2;gcd_a_H_nontrivial_iff_D=0"
    )
    print(
        f"case={name};base_ideal=(a,b)_is_transverse;"
        "Jacobian_det=a_prime*b1_is_a_unit=1;"
        "quadratic_subresultant_nonzero_pattern=(1,1,1);"
        "generic_y_gcd_degree_on_D=2;R2_y_lead_unit_on_D=1"
    )
    print(
        f"case={name};H_x_at_incidence=0;H_xx_at_incidence_is_unit=1;"
        "D_incidence_exact_x_multiplicity=2"
    )
    print(
        f"case={name};finite_wall_degree=5;finite_factor_profile=(1,1,1,2);"
        "boundary_wall_degree=6;boundary_factor_profile=(1,1,1,1,2);"
        "boundary=field_unit*finite*live_q3;live_q3_degree=1;"
        "exceptional_q4_midpoint=1;D_disjoint_all_named_walls=1"
    )
    digest = poly_digest(D)
    print(f"case={name};D_digest={digest}")
    return digest


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
        "no indirect execution mechanism",
    )
    audited_primary = (
        "jc_affine_c_parameter_section_"
        + "degeneracy_partial_scout_20260803.py"
    )
    require(audited_primary not in source, "audited primary filename absent")
    print(
        f"source_ast=(assert_nodes={assert_nodes},float_literals={float_literals},"
        "audited_primary_import_or_execution_paths=0)"
    )


def main() -> None:
    print("independent affine-C critical-section degeneracy audit")
    for relative_path, expected_hash in DEPENDENCIES.items():
        require(
            lf_hash(ROOT / relative_path) == expected_hash,
            f"dependency drift: {relative_path}",
        )
    print(f"dependency_hash_checks={len(DEPENDENCIES)}")
    packet = derive_generic_packet()
    digests = tuple(
        audit_case(response_case(name), packet)
        for name in ("4111", "3211")
    )
    print(f"independent_D_digests={digests}")
    print(
        "scope=FINITE-EXACT_PARTIAL_independent_quotient_norm_audit;"
        "full_repeated_H_parameter_discriminant_not_claimed;"
        "no_inverse_cover_not_JC2"
    )
    source_audit()


if __name__ == "__main__":
    main()
