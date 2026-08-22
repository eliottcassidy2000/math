#!/usr/bin/env python3
"""Finite-exact parameter scout for the THM-3289 critical section.

Fix d=k=1 and C=c+x in each THM-3212 accessory field.  This script
constructs the degree-53 saturated resultant H(c,x) and the primitive linear
subresultant a(x)y+b(c,x).  It replaces a prohibitively large direct
parameter resultant by an exact quotient-algebra norm computation.

The output is a FINITE-EXACT PARTIAL parameter result plus a bounded rational
squarefreeness census.  It is not an inverse cover and proves neither JC(2)
nor any case of the planar Jacobian conjecture.
"""

from __future__ import annotations

import ast
from fractions import Fraction
from hashlib import sha256
from math import gcd
from pathlib import Path

import sympy as sp
from sympy.polys.domains import QQ
from sympy.polys.matrices import DomainMatrix


ROOT = Path(__file__).resolve().parents[1]
DEPENDENCIES = {
    "04-computation/jc_affine_transverse_c0_e0_coupled_clutch_no_go_thm3289.py":
        "f63ff06e3f5ed30f3f6bc5be99756c347d6af5f8e9b220ce8336abff2cd2ca31",
    "05-knowledge/results/jc_affine_transverse_c0_e0_coupled_clutch_no_go_thm3289.out":
        "1aef4341650cdfaf1043a8699e3a1725a0100af6d9848d99dfa924b6f054dba1",
    "04-computation/jc_affine_c0_e0_linear_subresultant_section_partial_scout_20260803.py":
        "eed05abc68745bb300b57c4a669f34edd2c2e959f2ac461bf1f0bd919763330c",
    "05-knowledge/results/jc_affine_c0_e0_linear_subresultant_section_partial_scout_20260803.out":
        "26175ccdac159f5a281c9eaba6e91ac4bd928f50bd7f005b518469deaa0cf600",
}
SCAN_HEIGHT = 12
EXPECTED = {
    "4111": {
        "D": "d2b1a4307d903a4873d86a5d1b4764f0b4056bbe3810dc77cd2338c41a510215",
        "boundary": "2d92a3ab24798afa7c25b21860902b7762cfccabd580802f9414711703fadcfc",
        "finite": "2cdbdfee4220ee6f21761c8bd17ba68e1032e9899037a3bf9b57af458c9872b4",
        "live": "7e30dfaba63a345f7f92109c6a3a0df51764ab6cbdd3bcb921a4b007af937524",
        "cclass": "2cab20bdd322521e6a4309a2c854cdefcb0022767bc449499bf8a4827136d031",
        "embeddings": ((1013, 29), (1019, 22), (1031, 33)),
    },
    "3211": {
        "D": "7d843359be408b76f51d88d808b8ebb6c1a63c8b9f4273dc2d34e078e040b8ad",
        "boundary": "895dd22395834a1a037f7a646abfa4e8c97396b4aba7dfdac33b8974f28ebaec",
        "finite": "fc4c97ddf8522caead685c01a8af6a60cad6e88e41d10bb8662ef3177dbaf37f",
        "live": "62828a1361ac027bec86de80c954d7bde2f3c9ba00a81c6c459e783d5f5949c1",
        "cclass": "f1e4a546b447a3f84a96dcc8488d4f32583818a54080329d275610d3fc9c7705",
        "embeddings": ((1013, 121), (1019, 39), (1021, 239)),
    },
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def lf_hash(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def canonical_poly_text(poly) -> str:
    return "|".join(
        f"{monomial}:{coefficient}"
        for monomial, coefficient in sorted(poly.to_dict().items())
    )


def poly_digest(poly) -> str:
    return sha256(canonical_poly_text(poly).encode("ascii")).hexdigest()


def factor_profile(poly) -> tuple[tuple[int, int], ...]:
    _, factors = poly.factor_list()
    return tuple((factor.degree(), exponent) for factor, exponent in factors)


def specialize_sparse(poly, base_values, polynomial_ring):
    result = polynomial_ring.zero
    for monomial, coefficient in poly.terms():
        term = polynomial_ring.convert(coefficient)
        for base, exponent in zip(base_values, monomial):
            term *= base**exponent
        result += term
    return result


def derive_generic_packet():
    y, A, C, V, derivative_V, d, k = sp.symbols("y A C V Vp d k")
    L = y**2 + y + C * V
    R1 = sp.expand(2 * L * (2 * y + 1) + V * A)
    R2 = sp.expand(
        V**3 * k + V**2 * y + L * (-derivative_V * y + 2 * V**2 * d)
    )
    rows = sp.subresultants(R1, R2, y)
    require(
        tuple(sp.Poly(row, y).degree() for row in rows) == (3, 3, 2, 1, 0),
        "generic subresultant profile",
    )
    resultant = sp.Poly(rows[-1], y).nth(0)
    quotient = sp.cancel(resultant / V**3)
    require(sp.denom(quotient) == 1, "V^3 divides the generic resultant")
    K_poly = sp.Poly(quotient, A, C, V, derivative_V, d, k)
    linear = sp.Poly(rows[-2], y)
    raw_a = sp.Poly(linear.nth(1), A, C, V, derivative_V, d, k)
    raw_b = sp.Poly(linear.nth(0), A, C, V, derivative_V, d, k)
    quadratic = sp.Poly(rows[-3], y)
    quadratic_coefficients = tuple(
        sp.Poly(quadratic.nth(degree), A, C, V, derivative_V, d, k)
        for degree in (2, 1, 0)
    )
    require(len(K_poly.terms()) == 40, "generic K has forty terms")
    require(K_poly.degree(C) == 2, "K is quadratic in C")
    require(raw_a.degree(C) == 0, "linear y coefficient is C-independent")
    require(raw_b.degree(C) == 1, "linear constant coefficient is affine in C")
    require(sp.Poly(R1, y).LC() == 4, "no common projective y-root at infinity")
    print(
        "generic=(deg_C(K),deg_C(a),deg_C(b))=(2,0,1);"
        "subresultants=(3,3,2,1,0);Res_y=V^3*K40;R1_y_lead=4;"
        "rows_used_only_up_to_common_field_units=1"
    )
    return K_poly, raw_a, raw_b, quadratic_coefficients


def quotient_charpoly(element, modulus, parameter_symbol):
    """Characteristic polynomial of multiplication by element mod modulus."""

    modulus = modulus.monic()
    degree = modulus.degree()
    x_generator = modulus.ring.gens[0]
    field = modulus.ring.domain
    columns = []
    for exponent in range(degree):
        column_poly = (element * x_generator**exponent) % modulus
        coefficient_dict = column_poly.to_dict()
        columns.append(
            [coefficient_dict.get((row,), field.zero) for row in range(degree)]
        )
    rows = [
        [columns[column][row] for column in range(degree)]
        for row in range(degree)
    ]
    coefficients = DomainMatrix(rows, (degree, degree), field).charpoly()
    parameter_ring = field.poly_ring(parameter_symbol)
    parameter = parameter_ring.gens[0]
    polynomial = parameter_ring.zero
    for index, coefficient in enumerate(coefficients):
        polynomial += coefficient * parameter ** (degree - index)
    return polynomial


def evaluate_parameter_poly_mod(poly, value, modulus):
    result = modulus.ring.zero
    for (degree,), coefficient in poly.to_dict().items():
        result += modulus.ring.ground_new(coefficient) * value**degree
    return result % modulus


def evaluate_bivariate_at_parameter_mod(poly, value, modulus):
    x_generator = modulus.ring.gens[0]
    result = modulus.ring.zero
    for (x_degree, parameter_degree), coefficient in poly.to_dict().items():
        result += (
            modulus.ring.ground_new(coefficient)
            * x_generator**x_degree
            * value**parameter_degree
        )
    return result % modulus


def rational_universe(height: int) -> tuple[Fraction, ...]:
    return tuple(
        sorted(
            {
                Fraction(numerator, denominator)
                for denominator in range(1, height + 1)
                for numerator in range(-height, height + 1)
                if gcd(abs(numerator), denominator) == 1
            }
        )
    )


def rational_mod(value, prime: int) -> int:
    numerator = int(value.numerator)
    denominator = int(value.denominator)
    require(denominator % prime != 0, "modular denominator is invertible")
    return numerator * pow(denominator, -1, prime) % prime


def algebraic_mod(value, prime: int, root: int) -> int:
    result = 0
    for coefficient in value.to_list():
        result = (result * root + rational_mod(coefficient, prime)) % prime
    return result


def coefficient_denominators_good(poly, prime: int) -> bool:
    return all(
        int(rational.denominator) % prime != 0
        for coefficient in poly.to_dict().values()
        for rational in coefficient.to_list()
    )


def good_residue_embeddings(accessory, coefficient_poly, count: int):
    embeddings = []
    derivative = accessory.diff()
    for prime in sp.primerange(1009, 5000):
        if int(accessory.LC()) % prime == 0:
            continue
        if not coefficient_denominators_good(coefficient_poly, prime):
            continue
        for root in range(prime):
            value = sum(
                rational_mod(coefficient, prime) * pow(root, degree, prime)
                for (degree,), coefficient in accessory.as_dict().items()
            ) % prime
            derivative_value = sum(
                rational_mod(coefficient, prime) * pow(root, degree, prime)
                for (degree,), coefficient in derivative.as_dict().items()
            ) % prime
            if value == 0 and derivative_value != 0:
                embeddings.append((prime, root))
                break
        if len(embeddings) == count:
            break
    require(len(embeddings) == count, "enough good residue embeddings")
    return tuple(embeddings)


def reduce_univariate_to_prime(poly, prime: int, root: int, variable):
    finite_ring = sp.GF(prime).poly_ring(variable)
    finite_variable = finite_ring.gens[0]
    result = finite_ring.zero
    for (degree,), coefficient in poly.to_dict().items():
        result += finite_ring.convert(
            algebraic_mod(coefficient, prime, root)
        ) * finite_variable**degree
    return result


def build_case(name: str, generic_packet) -> None:
    K_poly, raw_a_poly, raw_b_poly, quadratic_coefficients = generic_packet
    u, x, c = sp.symbols("u x c")
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
    derivative_V = V.diff(X)
    g = S * T
    boundary = S**3 * T**8 * X**extras[0] * (X - 1) ** extras[1]
    require(boundary.degree(X) == 44, f"{name} boundary degree")
    require(2 * V * A.diff(X) - A * derivative_V == 2 * V, f"{name} response")

    values = (A, CP + X, V, derivative_V, ring.one, ring.one)
    K_control = specialize_sparse(K_poly, values, ring)
    H = K_control.exquo(boundary)
    primitive_a = specialize_sparse(raw_a_poly, values, ring).exquo(boundary)
    primitive_b = specialize_sparse(raw_b_poly, values, ring).exquo(boundary)
    require(
        (
            H.degree(CP),
            H.degree(X),
            primitive_a.degree(CP),
            primitive_a.degree(X),
            primitive_b.degree(CP),
            primitive_b.degree(X),
        )
        == (2, 53, 0, 36, 1, 44),
        f"{name} bidegrees",
    )

    leading_unit = primitive_a.LC
    a_bivariate = primitive_a * (field.one / leading_unit)
    b_bivariate = primitive_b * (field.one / leading_unit)
    zero = field.zero
    one = field.one
    minus_one = field.convert(-1)
    a = a_bivariate.evaluate(CP, zero).monic()
    b0 = b_bivariate.evaluate(CP, zero)
    b1 = b_bivariate.evaluate(CP, one) - b0
    H0 = H.evaluate(CP, zero)
    H_plus = H.evaluate(CP, one)
    H_minus = H.evaluate(CP, minus_one)
    H1 = (H_plus - H_minus) / 2
    H2 = (H_plus + H_minus) / 2 - H0
    require(
        (H0.degree(), H1.degree(), H2.degree(), b0.degree(), b1.degree())
        == (53, 52, 36, 44, 36),
        f"{name} coefficient degrees",
    )
    require(factor_profile(a) == ((36, 1),), f"{name} a irreducible")
    require(a.gcd(a.diff(a.ring.gens[0])).degree() == 0, f"{name} a squarefree")
    require(a.gcd(b1).degree() == 0, f"{name} b1 is a unit modulo a")
    require(a.gcd(H2).degree() == 0, f"{name} H2 is a unit modulo a")

    inverse_b1, _, gcd_b1 = b1.gcdex(a)
    inverse_b1 *= field.one / gcd_b1.LC
    c_class = (-b0 * inverse_b1) % a
    require(c_class.degree() == 24, f"{name} c-class degree")
    require((b0 + b1 * c_class) % a == a.ring.zero, f"{name} b graph")
    require(
        (H1 + 2 * H2 * c_class) % a == a.ring.zero,
        f"{name} double-root derivative identity",
    )
    require(
        (H0 - H2 * c_class**2) % a == a.ring.zero,
        f"{name} double-root constant identity",
    )

    degeneracy = quotient_charpoly(c_class, a, c).monic()
    parameter = degeneracy.ring.gens[0]
    require(degeneracy.degree() == 36, f"{name} degeneracy degree")
    require(
        degeneracy.gcd(degeneracy.diff(parameter)).degree() == 0,
        f"{name} degeneracy squarefree",
    )
    require(
        factor_profile(degeneracy) == ((36, 1),),
        f"{name} degeneracy irreducible",
    )
    require(
        evaluate_parameter_poly_mod(degeneracy, c_class, a) == a.ring.zero,
        f"{name} Cayley-Hamilton",
    )
    require(
        a.gcd(a.diff(a.ring.gens[0]) * b1).degree() == 0,
        f"{name} transverse coefficient-pair base locus",
    )

    specialized_quadratic = tuple(
        specialize_sparse(coefficient, values, ring)
        for coefficient in quadratic_coefficients
    )
    reduced_quadratic = tuple(
        evaluate_bivariate_at_parameter_mod(coefficient, c_class, a)
        for coefficient in specialized_quadratic
    )
    quadratic_nonzero_pattern = tuple(
        int(coefficient != a.ring.zero) for coefficient in reduced_quadratic
    )
    require(any(quadratic_nonzero_pattern), f"{name} quadratic row survives")

    # These are invariant zero-locus and multiplicity statements.  No raw PRS
    # representative is declared canonical (MISTAKE-360).
    print(
        f"case={name};H_bidegree=(c:2,x:53);a_degree=36;"
        "a_irreducible=1;b=(b0+c*b1);b1_unit_mod_a=1;"
        f"cstar_degree={c_class.degree()};"
        "H_mod_a=H2*(c-cstar)^2;H2_unit_mod_a=1"
    )
    print(
        f"case={name};degeneracy_D_degree={degeneracy.degree()};"
        "D_squarefree_irreducible=1;"
        "Res_x(a,H)=field_unit*D^2;Res_x(a,b)=field_unit*D;"
        "gcd(a,H)!=1_iff_D(c)=0;"
        "on_D=(a,b)=(0,0)_so_projective_selector_has_base_locus;"
        "base_Jacobian_det=a_prime*b1_is_a_unit=1;"
        f"quadratic_subresultant_nonzero_pattern={quadratic_nonzero_pattern};"
        "generic_y_gcd_degree_on_D=2"
    )

    finite_delta = ring.one - A.diff(X) * (CP + X)
    boundary_resultant = g.resultant(H).monic()
    finite_resultant = g.resultant(finite_delta).monic()
    require(boundary_resultant.ring == degeneracy.ring, f"{name} parameter ring")
    require(
        (
            boundary_resultant.degree(),
            finite_resultant.degree(),
            factor_profile(boundary_resultant),
            factor_profile(finite_resultant),
        )
        == (6, 5, ((1, 1), (1, 1), (1, 1), (1, 1), (2, 1)),
            ((1, 1), (1, 1), (1, 1), (2, 1))),
        f"{name} boundary and finite profiles",
    )
    live_S_wall = boundary_resultant.exquo(finite_resultant).monic()
    require(live_S_wall.degree() == 1, f"{name} live S wall")
    require(
        degeneracy.gcd(boundary_resultant).degree() == 0
        and degeneracy.gcd(finite_resultant).degree() == 0,
        f"{name} D disjoint from old walls",
    )

    x_ring = a.ring
    XX = x_ring.gens[0]
    S_univariate = S.evaluate(CP, zero)
    s_root = -S_univariate.to_dict()[(0,)]
    V_bar = V.exquo(S).evaluate(CP, zero)
    v1 = V_bar.evaluate(XX, s_root)
    v2 = V_bar.diff(XX).evaluate(XX, s_root)
    half = field.convert(sp.Rational(1, 2))
    local_finite = half
    local_live = -half - 2 * v2 / (3 * v1**2)
    local_exceptional = -v2 / (3 * v1**2)
    wall_finite = local_finite - s_root
    wall_live = local_live - s_root
    wall_exceptional = local_exceptional - s_root
    require(
        wall_exceptional == (wall_finite + wall_live) / 2,
        f"{name} exceptional midpoint",
    )
    finite_S_factor = parameter - wall_finite
    predicted_live = parameter - wall_live
    exceptional_factor = parameter - wall_exceptional
    require(
        finite_resultant % finite_S_factor == degeneracy.ring.zero,
        f"{name} finite S factor",
    )
    require(live_S_wall == predicted_live.monic(), f"{name} q3 live factor")
    require(
        boundary_resultant.evaluate(parameter, wall_exceptional) != field.zero
        and finite_resultant.evaluate(parameter, wall_exceptional) != field.zero,
        f"{name} exceptional midpoint is not an S-boundary root",
    )
    for wall in (wall_finite, wall_live, wall_exceptional):
        require(
            degeneracy.evaluate(parameter, wall) != field.zero,
            f"{name} D avoids named wall",
        )

    print(
        f"case={name};boundary_parameter_degree=6;boundary_factor_profile="
        "(1,1,1,1,2);finite_gate_degree=5;finite_factor_profile=(1,1,1,2);"
        "boundary=field_unit*finite_gate*live_q3_S_wall;gcd(D,boundary*finite)=1"
    )
    print(
        f"case={name};S_walls=finite(k=2C_S)|live_q3|exceptional_q4;"
        "exceptional_q4_is_exact_midpoint_of_finite_and_live_q3;"
        "exceptional_not_boundary_or_finite_root=1;D_avoids_all_three=1"
    )

    # At the distinguished D-incidence, H has an exact double x-root.
    Hx_on_D = (
        H0.diff(XX) + H1.diff(XX) * c_class + H2.diff(XX) * c_class**2
    ) % a
    Hxx_on_D = (
        H0.diff(XX).diff(XX) + H1.diff(XX).diff(XX) * c_class
        + H2.diff(XX).diff(XX) * c_class**2
    ) % a
    incidence_derivative_degrees = (
        a.gcd(Hx_on_D).degree(), a.gcd(Hxx_on_D).degree()
    )
    require(
        incidence_derivative_degrees == (36, 0),
        f"{name} exact double x-root along D",
    )
    print(
        f"case={name};D_incidence_derivative_gcd_degrees="
        f"{incidence_derivative_degrees}"
    )

    # Eliminating c from H=H_x=0 is cheap because both are quadratic in c.
    # Factoring/projecting the resulting degree-191 x-polynomial is not part
    # of the replay; repeated-H is therefore censused only on a stated finite
    # rational universe below.
    derivative_H0 = H0.diff(XX)
    derivative_H1 = H1.diff(XX)
    derivative_H2 = H2.diff(XX)
    P = H2 * derivative_H0 - derivative_H2 * H0
    Q = H2 * derivative_H1 - derivative_H2 * H1
    R = H1 * derivative_H0 - derivative_H1 * H0
    repeated_x_eliminant = P**2 - Q * R
    require(
        (P.degree(), Q.degree(), R.degree(), repeated_x_eliminant.degree())
        == (88, 87, 104, 191),
        f"{name} repeated-x eliminant degrees",
    )
    repeated_a_order = 0
    repeated_residual = repeated_x_eliminant
    while True:
        quotient, remainder = divmod(repeated_residual, a)
        if remainder != a.ring.zero:
            break
        repeated_a_order += 1
        repeated_residual = quotient
    require(
        (repeated_a_order, repeated_residual.degree()) == (2, 119),
        f"{name} repeated-H distinguished factor",
    )

    universe = rational_universe(SCAN_HEIGHT)
    residue_embeddings = good_residue_embeddings(accessory, H, 3)
    require(
        residue_embeddings == EXPECTED[name]["embeddings"],
        f"{name} residue embedding drift",
    )
    embedding_use = [0 for _ in residue_embeddings]
    rational_failures = []
    rational_gate_hits = []
    for rational_c in universe:
        c_value = field.convert(
            sp.Rational(rational_c.numerator, rational_c.denominator)
        )
        H_c = H.evaluate(CP, c_value)
        symbolic_gate_values = (
            degeneracy.evaluate(parameter, c_value),
            boundary_resultant.evaluate(parameter, c_value),
            finite_resultant.evaluate(parameter, c_value),
        )
        modular_certificate = None
        for embedding_index, (prime, root) in enumerate(residue_embeddings):
            reduced_H = reduce_univariate_to_prime(H_c, prime, root, x)
            finite_x = reduced_H.ring.gens[0]
            if reduced_H.degree() != 53:
                continue
            if reduced_H.gcd(reduced_H.diff(finite_x)).degree() == 0:
                modular_certificate = embedding_index
                embedding_use[embedding_index] += 1
                break
        gate_pattern = tuple(
            int(value == field.zero) for value in symbolic_gate_values
        )
        if gate_pattern[1] or gate_pattern[2]:
            rational_gate_hits.append((str(rational_c), gate_pattern))
        if gate_pattern[0] or modular_certificate is None:
            rational_failures.append((str(rational_c), modular_certificate))
    if rational_failures:
        print(f"case={name};rational_failures={rational_failures}")
    require(not rational_failures, f"{name} rational hostile census")
    require(
        rational_gate_hits
        == [("-2", (0, 1, 1)), ("-3/2", (0, 1, 1))],
        f"{name} rational wall hits",
    )
    require(Fraction(1, 1) in universe and Fraction(1, 2) in universe, "controls")
    print(
        f"case={name};repeated_H_symbolic_c_eliminant_x_degree=191;"
        f"a_factor_order={repeated_a_order};residual_x_degree="
        f"{repeated_residual.degree()};D_incidence_exact_x_multiplicity=2;"
        "full_residual_parameter_discriminant_not_claimed"
    )
    print(
        f"case={name};rational_scan=max(abs(p),q)<={SCAN_HEIGHT};"
        f"reduced_parameters={len(universe)};includes=(1,1/2);"
        f"good_residue_embeddings={residue_embeddings};"
        f"embedding_use={tuple(embedding_use)};"
        f"boundary_or_finite_gate_hits={rational_gate_hits};"
        "failures=(D,modular_H_Hprime_certificate)=0"
    )

    digests = {
        "D": poly_digest(degeneracy),
        "boundary": poly_digest(boundary_resultant),
        "finite": poly_digest(finite_resultant),
        "live": poly_digest(live_S_wall),
        "cclass": poly_digest(c_class),
    }
    for label, digest in digests.items():
        expected = EXPECTED[name][label]
        if expected != "TBD":
            require(digest == expected, f"{name} {label} digest drift")
    print(f"case={name};digests={digests}")


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
    print("finite-exact partial affine-c critical-section degeneracy scout")
    for dependency, expected_hash in DEPENDENCIES.items():
        require(lf_hash(ROOT / dependency) == expected_hash, f"dependency drift: {dependency}")
    print(f"dependency_hash_checks={len(DEPENDENCIES)}")
    generic_packet = derive_generic_packet()
    for name in ("4111", "3211"):
        build_case(name, generic_packet)
    print(
        "scope=FINITE-EXACT_PARTIAL_exact_D_and_wall_loci_plus_bounded_rational_scan;"
        "full_repeated_H_parameter_discriminant_not_claimed;"
        "no_inverse_cover_not_JC2"
    )
    source_audit()


if __name__ == "__main__":
    main()
