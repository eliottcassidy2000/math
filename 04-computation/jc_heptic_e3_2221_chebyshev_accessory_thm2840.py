from __future__ import annotations

from itertools import permutations

import sympy as sp


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


x, lam, mu = sp.symbols("x lambda mu")


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
    require(not remainder.has(mu), "triangular quotient did not eliminate mu")
    remainder_poly = sp.Poly(remainder, lam, domain=sp.QQ)
    modulus_poly = sp.Poly(modulus, lam, domain=sp.QQ)
    if sp.gcd(remainder_poly, modulus_poly).degree() != 0:
        return False
    inverse = sp.invert(remainder_poly.as_expr(), modulus_poly.as_expr(), domain=sp.QQ)
    return (
        sp.rem(
            remainder_poly.as_expr() * inverse - 1,
            modulus_poly.as_expr(),
            domain=sp.QQ,
        )
        == 0
    )


def all_three_edge_matchings(points: tuple[int, ...], edge_count: int):
    if edge_count == 0:
        yield ()
        return
    if len(points) < 2 * edge_count:
        return
    first = points[0]
    yield from all_three_edge_matchings(points[1:], edge_count)
    for index in range(1, len(points)):
        second = points[index]
        remainder = points[1:index] + points[index + 1 :]
        for tail in all_three_edge_matchings(remainder, edge_count - 1):
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


def rotate_set(points: frozenset[int], shift: int) -> tuple[int, ...]:
    return tuple(sorted((point + shift) % 7 for point in points))


def marked_rotation_key(
    involution: tuple[int, ...], labelled_cycles: tuple[frozenset[int], ...]
):
    candidates = []
    for shift in range(7):
        rotated_involution = tuple(
            (involution[(point - shift) % 7] + shift) % 7 for point in range(7)
        )
        rotated_cycles = tuple(rotate_set(cycle, shift) for cycle in labelled_cycles)
        candidates.append((rotated_involution, rotated_cycles))
    return min(candidates)


def unmarked_rotation_key(involution: tuple[int, ...]) -> tuple[int, ...]:
    return min(
        tuple((involution[(point - shift) % 7] + shift) % 7 for point in range(7))
        for shift in range(7)
    )


def nielsen_census():
    total_cycle = tuple((point + 1) % 7 for point in range(7))
    identity = tuple(range(7))
    unique_matchings = {
        tuple(sorted(tuple(sorted(edge)) for edge in matching))
        for matching in all_three_edge_matchings(tuple(range(7)), 3)
    }
    require(len(unique_matchings) == 105, "wrong number of three-edge matchings")

    planar = []
    passport_histogram: dict[tuple[int, ...], int] = {}
    for matching in unique_matchings:
        involution_list = list(range(7))
        for left, right in matching:
            involution_list[left], involution_list[right] = right, left
        involution = tuple(involution_list)
        pole_permutation = inverse(compose(involution, total_cycle))
        require(
            compose(compose(involution, total_cycle), pole_permutation) == identity,
            "branch product failed",
        )
        pole_cycles = cycles(pole_permutation)
        if len(pole_cycles) != 4:
            continue
        passport = tuple(sorted((len(cycle) for cycle in pole_cycles), reverse=True))
        passport_histogram[passport] = passport_histogram.get(passport, 0) + 1
        planar.append((involution, tuple(frozenset(cycle) for cycle in pole_cycles), passport))

    require(len(planar) == 35, "wrong genus-zero matching count")
    require(
        passport_histogram
        == {(4, 1, 1, 1): 7, (3, 2, 1, 1): 21, (2, 2, 2, 1): 7},
        "wrong planar passport histogram",
    )

    target = (2, 2, 2, 1)
    marked = set()
    unmarked = set()
    for involution, pole_cycles, passport in planar:
        if passport != target:
            continue
        unmarked.add(unmarked_rotation_key(involution))
        for labelled_cycles in permutations(pole_cycles):
            if tuple(len(cycle) for cycle in labelled_cycles) == target:
                marked.add(marked_rotation_key(involution, labelled_cycles))
    require(len(marked) == 6, "ordered passport does not have six Nielsen charts")
    require(len(unmarked) == 1, "passport does not have one unmarked Nielsen class")
    return passport_histogram, len(marked), len(unmarked)


def main() -> None:
    divisor = sp.expand(x**2 * (x - 1) ** 2 * (x - lam) ** 2 * (x - mu))
    pole_product = sp.expand(x * (x - 1) * (x - lam) * (x - mu))
    forced_pole_factor = sp.cancel(divisor / pole_product)
    critical = sp.expand(
        2 * (x - 1) * (x - lam) * (x - mu)
        + 2 * x * (x - lam) * (x - mu)
        + 2 * x * (x - 1) * (x - mu)
        + x * (x - 1) * (x - lam)
    )
    energy = critical / 7
    require(
        sp.expand(sp.diff(divisor, x) - 7 * forced_pole_factor * energy) == 0,
        "logarithmic derivative reconstruction failed",
    )
    require(sp.Poly(energy, x).degree() == 3 and sp.Poly(energy, x).LC() == 1, "E not monic")

    simple_factor, remainder = sp.div(divisor, energy**2, x)
    derivative_quotient, derivative_remainder = sp.div(sp.diff(remainder, x), energy, x)
    require(derivative_remainder == 0, "R' is not exactly divisible by E")
    require(sp.Poly(simple_factor, x).degree() == 1, "simple-zero quotient is not linear")
    require(sp.Poly(simple_factor, x).LC() == 1, "simple-zero quotient is not monic")
    derivative_poly = sp.Poly(derivative_quotient, x)

    equation_one = sp.factor(
        -sp.Rational(49, 40) * derivative_poly.coeff_monomial(x)
    )
    equation_zero = sp.factor(
        sp.Rational(49, 8) * derivative_poly.coeff_monomial(1)
    )
    expected_one = 2 * lam**2 + 2 * lam * mu - 3 * lam - 3 * mu**2 + 2 * mu + 2
    expected_zero = (
        4 * lam**2 * mu
        + 3 * lam**2
        - 5 * lam * mu**2
        + 6 * lam * mu
        + 3 * lam
        - 5 * mu**2
        + 4 * mu
    )
    require(
        sp.expand(equation_one - expected_one) == 0
        and sp.expand(equation_zero - expected_zero) == 0,
        "accessory equations differ",
    )

    quotient = sp.groebner(
        [equation_one, equation_zero], mu, lam, order="lex", domain=sp.QQ
    )
    p = lam**3 - 2 * lam**2 - lam + 1
    q_cubic = lam**3 - lam**2 - 2 * lam + 1
    modulus = sp.expand(p * q_cubic)
    expected_relation = (
        2 * lam**5
        - 5 * lam**4
        - 6 * lam**3
        + 14 * lam**2
        + 6 * lam
        + mu
        - 6
    )
    require(
        len(quotient.polys) == 2
        and quotient.reduce(expected_relation)[1] == 0
        and quotient.reduce(modulus)[1] == 0,
        "six-point triangular basis differs",
    )
    require(
        sp.discriminant(p, lam) == 49
        and sp.discriminant(q_cubic, lam) == 49
        and sp.resultant(p, q_cubic, lam) == -1
        and sp.discriminant(modulus, lam) == 2401,
        "sextic accessory modulus is not squarefree",
    )

    first_branch = sp.groebner(
        [equation_one, equation_zero, p], mu, lam, order="lex", domain=sp.QQ
    )
    second_branch = sp.groebner(
        [equation_one, equation_zero, q_cubic],
        mu,
        lam,
        order="lex",
        domain=sp.QQ,
    )
    require(
        first_branch.reduce(mu - lam**2 + lam)[1] == 0
        and second_branch.reduce(mu + lam**2 - lam - 1)[1] == 0,
        "cubic branch relations differ",
    )

    constant_remainder = sp.factor(remainder.subs(x, 0))
    accessory = -constant_remainder
    numerator = divisor + accessory
    alpha = sp.factor(-simple_factor.subs(x, 0))
    require(
        sp.expand(simple_factor - (x - alpha)) == 0,
        "simple zero is not x-alpha",
    )

    collision = lam * mu * (lam - 1) * (mu - 1) * (lam - mu)
    discriminant_gate = sp.discriminant(energy, x)
    zero_separation = sp.resultant(simple_factor, energy, x)
    affine_scale_denominator = 6 * mu - 2 * lam - 2
    gates = {
        "collision": collision,
        "critical_discriminant": discriminant_gate,
        "nonzero_accessory": accessory,
        "simple_double_separation": zero_separation,
        "affine_scale": affine_scale_denominator,
    }
    for name, gate in gates.items():
        require(
            quotient_unit(gate, quotient, modulus),
            f"{name} is not a unit on the accessory algebra",
        )

    jacobian = sp.det(
        sp.Matrix(
            (
                (sp.diff(equation_one, lam), sp.diff(equation_one, mu)),
                (sp.diff(equation_zero, lam), sp.diff(equation_zero, mu)),
            )
        )
    )
    require(
        quotient_unit(jacobian, quotient, modulus),
        "the six accessory points are not simple",
    )
    require(
        polynomial_zero_mod(numerator - simple_factor * energy**2, x, quotient),
        "D+A=S E^2 failed",
    )

    response_constant = -7 * accessory
    derivative_cross_product = (
        (sp.diff(numerator, x) * divisor - numerator * sp.diff(divisor, x))
        * pole_product
        - response_constant * energy * divisor
    )
    require(
        polynomial_zero_mod(derivative_cross_product, x, quotient),
        "F'=2G failed",
    )
    # With G=C E/(2DT) and V=4 S D T^2/C^2, the already checked
    # B=S E^2 is exactly F=VG^2.  Differentiating it and using F'=2G
    # gives 2VG'+V'G=2 without another large rational simplification.

    affine_scale = 7 / affine_scale_denominator
    chebyshev_coordinate = 1 + affine_scale * (x - mu)
    chebyshev_7 = sp.chebyshevt(7, chebyshev_coordinate)
    require(
        polynomial_zero_mod(
            accessory * chebyshev_7 - (2 * divisor + accessory), x, quotient
        ),
        "degree-seven Chebyshev identification failed",
    )
    # A*T7(Y)=2D+A is equivalent by one cross multiplication to
    # (D+A)/D=(T7(Y)+1)/(T7(Y)-1).
    require(
        polynomial_zero_mod(
            chebyshev_coordinate.subs(x, alpha) + 1, x, quotient
        )
        and polynomial_zero_mod(
            chebyshev_coordinate.subs(x, mu) - 1, x, quotient
        ),
        "simple zero/pole do not normalize to -1/+1",
    )

    histogram, marked_count, unmarked_count = nielsen_census()
    require(marked_count == sp.degree(modulus, lam), "algebraic and Nielsen counts differ")

    print("THM-2840 HEPTIC e=3 PASSPORT (2,2,2,1) CHEBYSHEV ACCESSORY AUDIT")
    print("universe=N=7,e=3,s=1,h=4; poles (0,1,lambda,mu) with (2,2,2,1)")
    print(f"equations={sp.factor(equation_one)} ; {sp.factor(equation_zero)}")
    print(f"accessory_modulus=({p})*({q_cubic})")
    print("radical_saturated_length=6; all collision/value/discriminant/separation gates are units")
    print("cubic_branches: mu=lambda^2-lambda ; mu=1+lambda-lambda^2")
    print("response=D+A=S*E^2 with F'=2G, F=VG^2, 2VG'+V'G=2")
    print("source_normalization=F(infinity)=1 with total order 7")
    print("carrier=F=(T7(Y)+1)/(T7(Y)-1), Y=1+7(x-mu)/(6mu-2lambda-2)")
    print(f"Nielsen_histogram={histogram}")
    print("Nielsen_target=6 marked charts, 1 unmarked class")
    print("ALL EXACT PROBE CHECKS PASSED")


if __name__ == "__main__":
    main()
