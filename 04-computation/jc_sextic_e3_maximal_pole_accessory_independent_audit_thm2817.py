from __future__ import annotations

from itertools import permutations

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


x, lam, mu = sp.symbols("x lambda mu")


ROWS = (
    ((1, 1, 1, 3), 3 * mu - lam - 1, lam**2 - lam + 1),
    ((1, 1, 2, 2), lam + mu - 1, (4 * lam - 3) * (4 * lam - 1)),
    ((1, 1, 3, 1), mu - 3 * lam + 1, 3 * lam**2 - 3 * lam + 1),
    ((1, 2, 1, 2), mu - lam + 1, (lam - 4) * (3 * lam - 4)),
    ((1, 2, 2, 1), mu - lam - 1, (lam - 3) * (3 * lam - 1)),
    ((1, 3, 1, 1), lam + mu - 3, lam**2 - 3 * lam + 3),
    ((2, 1, 1, 2), mu - lam - 1, (lam + 3) * (3 * lam + 1)),
    ((2, 1, 2, 1), mu - lam + 1, (lam + 2) * (3 * lam - 2)),
    ((2, 2, 1, 1), lam + mu - 1, (2 * lam - 3) * (2 * lam + 1)),
    ((3, 1, 1, 1), lam + mu + 1, lam**2 + lam + 1),
)


def monic(expression: sp.Expr, variable: sp.Symbol) -> sp.Expr:
    return sp.Poly(expression, variable, domain=sp.QQ).monic().as_expr()


def ideals_equal(first: sp.GroebnerBasis, second: sp.GroebnerBasis) -> bool:
    return all(first.reduce(poly.as_expr())[1] == 0 for poly in second.polys) and all(
        second.reduce(poly.as_expr())[1] == 0 for poly in first.polys
    )


def polynomial_zero_mod(
    expression: sp.Expr, variable: sp.Symbol, quotient: sp.GroebnerBasis
) -> bool:
    numerator = sp.cancel(expression).as_numer_denom()[0]
    polynomial = sp.Poly(sp.expand(numerator), variable)
    return all(quotient.reduce(coefficient)[1] == 0 for coefficient in polynomial.all_coeffs())


def quotient_unit(
    expression: sp.Expr, quotient: sp.GroebnerBasis, modulus: sp.Expr
) -> bool:
    remainder = sp.factor(quotient.reduce(expression)[1])
    require(not remainder.has(mu), "linear quotient did not eliminate mu")
    remainder_poly = sp.Poly(remainder, lam, domain=sp.QQ)
    modulus_poly = sp.Poly(modulus, lam, domain=sp.QQ)
    if sp.gcd(remainder_poly, modulus_poly).degree() != 0:
        return False
    inverse = sp.invert(remainder_poly.as_expr(), modulus_poly.as_expr(), domain=sp.QQ)
    return sp.rem(
        remainder_poly.as_expr() * inverse - 1,
        modulus_poly.as_expr(),
        domain=sp.QQ,
    ) == 0


def sextic_data(passport: tuple[int, int, int, int]):
    a, b, c, d = passport
    divisor = sp.expand(x**a * (x - 1) ** b * (x - lam) ** c * (x - mu) ** d)
    pole_product = sp.expand(x * (x - 1) * (x - lam) * (x - mu))
    critical = sp.expand(
        a * (x - 1) * (x - lam) * (x - mu)
        + b * x * (x - lam) * (x - mu)
        + c * x * (x - 1) * (x - mu)
        + d * x * (x - 1) * (x - lam)
    )
    _, remainder = sp.div(divisor, critical, x)
    remainder_poly = sp.Poly(remainder, x)
    r2 = remainder_poly.coeff_monomial(x**2)
    r1 = remainder_poly.coeff_monomial(x)
    r0 = remainder_poly.coeff_monomial(1)
    return divisor, pole_product, critical, r2, r1, r0


def factor_signature(factors: tuple[tuple[sp.Expr, int], ...]) -> str:
    return ",".join(str(sp.degree(factor, lam)) for factor, _ in factors)


def audit_algebraic_row(
    passport: tuple[int, int, int, int], linear: sp.Expr, quadratic: sp.Expr
):
    divisor, pole_product, critical, r2, r1, r0 = sextic_data(passport)
    require(
        sp.expand(pole_product * sp.diff(divisor, x) - divisor * critical) == 0,
        f"{passport}: logarithmic derivative identity failed",
    )
    require(sp.degree(divisor, x) == 6 and sp.LC(sp.Poly(divisor, x)) == 1, "D not monic sextic")
    require(sp.degree(critical, x) == 3 and sp.LC(sp.Poly(critical, x)) == 6, "K not 6-monic cubic")

    quotient = sp.groebner([linear, quadratic], mu, lam, order="lex", domain=sp.QQ)
    isolated = sp.groebner([r2, r1, quadratic], mu, lam, order="lex", domain=sp.QQ)
    require(ideals_equal(quotient, isolated), f"{passport}: proposed local ideal is wrong")
    require(
        sp.gcd(sp.Poly(quadratic, lam), sp.Poly(sp.diff(quadratic, lam), lam)).degree()
        == 0,
        f"{passport}: accessory quadratic is not squarefree",
    )

    discriminant = sp.discriminant(critical, x)
    collision = lam * mu * (lam - 1) * (mu - 1) * (lam - mu)
    value_gate = -r0
    gates = (collision, discriminant, value_gate)
    require(
        all(quotient_unit(gate, quotient, quadratic) for gate in gates),
        f"{passport}: a collision/discriminant/value gate is not a unit",
    )
    gamma = sp.expand(collision * discriminant * value_gate)

    resultant = sp.resultant(r2, r1, mu)
    resultant_factors = tuple(sp.factor_list(resultant, lam)[1])
    require(resultant != 0, f"{passport}: projection resultant vanished")
    require(
        all(exponent == 1 for _, exponent in resultant_factors),
        f"{passport}: projection resultant has a repeated branch",
    )

    expected_factors = tuple(sp.factor_list(quadratic, lam)[1])
    expected_monic = {monic(factor, lam) for factor, _ in expected_factors}
    admissible_monic: set[sp.Expr] = set()
    rejected_degrees = []
    for factor, _ in resultant_factors:
        branch = sp.groebner([r2, r1, factor], mu, lam, order="lex", domain=sp.QQ)
        require(len(branch.polys) == 2, f"{passport}: non-triangular resultant branch")
        branch_gamma = sp.factor(branch.reduce(gamma)[1])
        factor_monic = monic(factor, lam)
        if factor_monic in expected_monic:
            require(
                quotient_unit(branch_gamma, branch, factor),
                f"{passport}: alleged admissible branch meets Gamma=0",
            )
            admissible_monic.add(factor_monic)
        else:
            require(
                branch_gamma == 0,
                f"{passport}: resultant branch survives saturation unexpectedly",
            )
            rejected_degrees.append(sp.degree(factor, lam))
    require(
        admissible_monic == expected_monic,
        f"{passport}: admissible resultant factors do not reconstruct the quadratic",
    )

    jacobian = sp.det(
        sp.Matrix(
            (
                (sp.diff(r2, lam), sp.diff(r2, mu)),
                (sp.diff(r1, lam), sp.diff(r1, mu)),
            )
        )
    )
    require(
        quotient_unit(jacobian, quotient, quadratic),
        f"{passport}: accessory points are not simple",
    )

    energy = critical / 6
    accessory = -r0
    numerator = divisor + accessory
    response_constant = -6 * accessory
    response = numerator / divisor
    response_g = response_constant * energy / (2 * divisor * pole_product)
    response_v = 4 * divisor * pole_product**2 / response_constant**2

    require(
        polynomial_zero_mod(numerator - energy**2, x, quotient),
        f"{passport}: B=E^2 reconstruction failed",
    )
    require(
        polynomial_zero_mod(sp.diff(response, x) - 2 * response_g, x, quotient),
        f"{passport}: F'=2G failed",
    )
    require(
        polynomial_zero_mod(response - response_v * response_g**2, x, quotient),
        f"{passport}: F=VG^2 failed",
    )
    require(
        polynomial_zero_mod(
            2 * response_v * sp.diff(response_g, x)
            + sp.diff(response_v, x) * response_g
            - 2,
            x,
            quotient,
        ),
        f"{passport}: response ODE failed",
    )
    require(
        polynomial_zero_mod(response - 1 - accessory / divisor, x, quotient),
        f"{passport}: source normalization identity failed",
    )
    for pole in (sp.Integer(0), sp.Integer(1), lam, mu):
        require(
            polynomial_zero_mod(numerator.subs(x, pole) - accessory, x, quotient),
            f"{passport}: numerator cancels at a pole",
        )

    coordinates = (sp.Integer(0), sp.Integer(1), lam, mu)
    if 3 in passport:
        triple_pole = coordinates[passport.index(3)]
        simple_poles = [
            coordinates[index] for index, multiplicity in enumerate(passport) if multiplicity == 1
        ]
        carrier_coordinate = (x - triple_pole) / (simple_poles[0] - triple_pole)
        inner = 2 * carrier_coordinate**3 - 1
        carrier = "power"
    else:
        double_poles = [
            coordinates[index] for index, multiplicity in enumerate(passport) if multiplicity == 2
        ]
        center = (double_poles[0] + double_poles[1]) / 2
        carrier_coordinate = (x - center) / (double_poles[1] - double_poles[0])
        inner = 4 * carrier_coordinate**3 - 3 * carrier_coordinate
        carrier = "Chebyshev"
    canonical_response = inner**2 / (inner**2 - 1)
    require(
        polynomial_zero_mod(response - canonical_response, x, quotient),
        f"{passport}: unmarked carrier normalization failed",
    )

    return {
        "passport": passport,
        "resultant_factors": factor_signature(resultant_factors),
        "admissible_factors": factor_signature(expected_factors),
        "rejected_degrees": ",".join(str(degree) for degree in rejected_degrees),
        "carrier": carrier,
    }


def all_matchings(points: tuple[int, ...]):
    if not points:
        yield ()
        return
    first = points[0]
    for index in range(1, len(points)):
        second = points[index]
        remainder = points[1:index] + points[index + 1 :]
        for tail in all_matchings(remainder):
            yield ((first, second),) + tail


def compose(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(left[right[index]] for index in range(len(left)))


def inverse(permutation: tuple[int, ...]) -> tuple[int, ...]:
    answer = [0] * len(permutation)
    for source, target in enumerate(permutation):
        answer[target] = source
    return tuple(answer)


def cycles(permutation: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    seen: set[int] = set()
    answer = []
    for start in range(len(permutation)):
        if start in seen:
            continue
        cycle = []
        point = start
        while point not in seen:
            seen.add(point)
            cycle.append(point)
            point = permutation[point]
        answer.append(tuple(cycle))
    return tuple(answer)


def crossing(first: tuple[int, int], second: tuple[int, int]) -> bool:
    a, b = sorted(first)
    c, d = sorted(second)
    return (a < c < b < d) or (c < a < d < b)


def rotate_set(points: frozenset[int], shift: int) -> tuple[int, ...]:
    return tuple(sorted((point + shift) % 6 for point in points))


def marked_rotation_key(
    involution: tuple[int, ...], labelled_cycles: tuple[frozenset[int], ...]
):
    candidates = []
    for shift in range(6):
        rotated_involution = tuple(
            (involution[(point - shift) % 6] + shift) % 6 for point in range(6)
        )
        rotated_cycles = tuple(rotate_set(cycle, shift) for cycle in labelled_cycles)
        candidates.append((rotated_involution, rotated_cycles))
    return min(candidates)


def audit_nielsen():
    total_cycle = tuple((point + 1) % 6 for point in range(6))
    identity = tuple(range(6))
    all_records = []
    planar_records = []
    noncrossing_records = []
    for matching in all_matchings(tuple(range(6))):
        involution_list = list(range(6))
        for left, right in matching:
            involution_list[left], involution_list[right] = right, left
        involution = tuple(involution_list)
        pole_permutation = inverse(compose(involution, total_cycle))
        require(
            compose(compose(involution, total_cycle), pole_permutation) == identity,
            "branch-cycle product is not the identity",
        )
        pole_cycles = cycles(pole_permutation)
        record = (involution, tuple(frozenset(cycle) for cycle in pole_cycles))
        all_records.append(record)
        if len(pole_cycles) == 4:
            planar_records.append(record)
        if not any(
            crossing(matching[left], matching[right])
            for left in range(len(matching))
            for right in range(left + 1, len(matching))
        ):
            noncrossing_records.append(record)

    require(len(all_records) == 15, "there are not 15 perfect matchings")
    require(len(planar_records) == 5, "there are not five genus-zero matchings")
    require(
        {record[0] for record in planar_records}
        == {record[0] for record in noncrossing_records},
        "genus zero and noncrossing conditions disagree",
    )

    counts = {}
    all_marked = set()
    for passport, _, _ in ROWS:
        marked_classes = set()
        for involution, pole_cycles in planar_records:
            for labelled in permutations(pole_cycles):
                if tuple(len(cycle) for cycle in labelled) != passport:
                    continue
                key = marked_rotation_key(involution, labelled)
                marked_classes.add(key)
                all_marked.add((passport, key))
        require(len(marked_classes) == 2, f"{passport}: Nielsen count is not two")
        counts[passport] = len(marked_classes)
    require(len(all_marked) == 20, "marked Nielsen atlas does not have 20 charts")
    return counts


def main() -> None:
    algebraic_records = [
        audit_algebraic_row(passport, linear, quadratic)
        for passport, linear, quadratic in ROWS
    ]
    nielsen_counts = audit_nielsen()
    require(sum(record["carrier"] == "power" for record in algebraic_records) == 4, "power orbit row count")
    require(
        sum(record["carrier"] == "Chebyshev" for record in algebraic_records) == 6,
        "Chebyshev orbit row count",
    )

    print("THM-2817 INDEPENDENT SEXTIC e=3 ACCESSORY AUDIT")
    print("method=Sylvester projection branches + exact quotient gates; no primary companion import")
    print("passport | resultant factor degrees | admissible degrees | rejected degrees | carrier | Nielsen")
    for record in algebraic_records:
        passport = record["passport"]
        print(
            f"{passport} | {record['resultant_factors']} | "
            f"{record['admissible_factors']} | {record['rejected_degrees']} | "
            f"{record['carrier']} | {nielsen_counts[passport]}"
        )
    print("collision_discriminant_value_gates=units on all 20 accessory points")
    print("accessory_ideals=ten radical linear-plus-squarefree-quadratic ideals")
    print("reconstruction=B=E^2; F'=2G, F=VG^2, 2VG'+V'G=2 verified coefficientwise")
    print("source_normalization=monic D degree 6; monic E degree 3; F(infinity)=1 with order 6")
    print("Nielsen=15 matchings; 5 planar/noncrossing; 2 per passport; 20 marked charts")
    print("unmarked_carriers=power,Chebyshev")
    print("ALL INDEPENDENT EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
