#!/usr/bin/env python3
"""Exact controls for the balanced THM-2796 Stieltjes--Pade theorem.

The script independently:

* enumerates clean genus-zero permutation triples through degree eight;
* isolates the genuine nonsplit squareclass packets;
* proves the e=0 cyclic and e=h=1 full-symmetric subchambers in bounded
  controls;
* verifies the two explicit all-degree rational families;
* checks the factorized Stieltjes identity, moment vanishing, interpolation
  constraints, and exact first square-prefix defect.
"""

from itertools import permutations
from math import factorial, gcd, lcm
from pathlib import Path
import ast

import sympy as sp


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


GATES = 0


def gate(condition, message):
    global GATES
    require(bool(condition), message)
    GATES += 1


def compose(left, right):
    return tuple(left[right[i]] for i in range(len(left)))


def inverse(perm):
    result = [0] * len(perm)
    for i, value in enumerate(perm):
        result[value] = i
    return tuple(result)


def conjugate(conjugator, perm):
    return compose(compose(conjugator, perm), inverse(conjugator))


def cycle_type(perm):
    unseen = set(range(len(perm)))
    lengths = []
    while unseen:
        start = min(unseen)
        current = start
        length = 0
        while current in unseen:
            unseen.remove(current)
            length += 1
            current = perm[current]
        lengths.append(length)
    return tuple(sorted(lengths, reverse=True))


def involutions(size, transpositions):
    def recurse(remaining, pairs):
        if not remaining:
            if len(pairs) == transpositions:
                perm = list(range(size))
                for left, right in pairs:
                    perm[left] = right
                    perm[right] = left
                yield tuple(perm)
            return
        if len(pairs) > transpositions:
            return
        if len(pairs) + len(remaining) // 2 < transpositions:
            return
        first = min(remaining)
        rest = remaining - {first}
        yield from recurse(rest, pairs)
        if len(pairs) < transpositions:
            for second in sorted(rest):
                yield from recurse(rest - {second}, pairs + ((first, second),))

    yield from recurse(set(range(size)), ())


def is_transitive(generators):
    reached = {0}
    frontier = [0]
    while frontier:
        current = frontier.pop()
        for generator in generators:
            image = generator[current]
            if image not in reached:
                reached.add(image)
                frontier.append(image)
    return len(reached) == len(generators[0])


def centralizer_of_standard_cycle(size, cycle_length):
    fixed = tuple(range(cycle_length, size))
    for shift in range(cycle_length):
        cycle_part = tuple((i + shift) % cycle_length for i in range(cycle_length))
        for tail in permutations(fixed):
            yield cycle_part + tail


def generated_group(generators):
    identity = tuple(range(len(generators[0])))
    group = {identity}
    frontier = [identity]
    while frontier:
        current = frontier.pop()
        for generator in generators:
            product_perm = compose(current, generator)
            if product_perm not in group:
                group.add(product_perm)
                frontier.append(product_perm)
    return group


def permutation_order(perm):
    unseen = set(range(len(perm)))
    order = 1
    while unseen:
        start = min(unseen)
        current = start
        length = 0
        while current in unseen:
            unseen.remove(current)
            length += 1
            current = perm[current]
        order = lcm(order, length)
    return order


def is_cyclic_group(group):
    # This is a genuine cyclicity test, not merely an abelianness proxy:
    # a finite group is cyclic exactly when it contains an element whose
    # order is the cardinality of the group.
    return any(permutation_order(element) == len(group) for element in group)


def census_degree(size):
    records = []
    if size == 1:
        return [
            {
                "e": 0,
                "h": 1,
                "s": 1,
                "r": 1,
                "p": (1,),
                "genuine": True,
                "cyclic": True,
                "group_order": 1,
            }
        ]

    for e in range(size // 2 + 1):
        s = size - 2 * e
        for h in range(1, e + 2):
            r = size - e + h - 1
            if r > size:
                continue
            sigma_one = (
                tuple(list(range(1, r)) + [0] + list(range(r, size)))
                if r > 1
                else tuple(range(size))
            )
            candidates = []
            for sigma_zero in involutions(size, e):
                if not is_transitive((sigma_zero, sigma_one)):
                    continue
                sigma_infinity = inverse(compose(sigma_zero, sigma_one))
                if len(cycle_type(sigma_infinity)) == h:
                    candidates.append(sigma_zero)

            remaining = set(candidates)
            centralizer = tuple(centralizer_of_standard_cycle(size, r))
            while remaining:
                representative = min(remaining)
                orbit = {
                    conjugate(element, representative)
                    for element in centralizer
                }
                remaining -= orbit
                sigma_infinity = inverse(compose(representative, sigma_one))
                pole_partition = cycle_type(sigma_infinity)
                group = generated_group((representative, sigma_one))
                genuine = s > 0 or any(part % 2 for part in pole_partition)
                records.append(
                    {
                        "e": e,
                        "h": h,
                        "s": s,
                        "r": r,
                        "p": pole_partition,
                        "genuine": genuine,
                        "cyclic": is_cyclic_group(group),
                        "group_order": len(group),
                    }
                )
    return records


x = sp.symbols("x")


def cancel(expr):
    return sp.cancel(sp.together(expr))


def exact_degree(expr):
    numerator, denominator = sp.fraction(cancel(expr))
    return max(sp.degree(numerator, x), sp.degree(denominator, x))


def root_power_sums(poly, maximum_power):
    """Newton sums of a monic polynomial, including p_0=degree.

    This stays exact for irreducible polynomials of arbitrarily high degree;
    unlike ``sympy.roots``, it never silently drops roots lacking a radical
    representation.
    """

    monic = sp.Poly(poly, x, domain=sp.QQ).monic()
    degree = monic.degree()
    coefficients = monic.all_coeffs()
    sums = [sp.Integer(degree)]
    for power in range(1, maximum_power + 1):
        if power <= degree:
            relation_tail = sum(
                coefficients[index] * sums[power - index]
                for index in range(1, power)
            )
            relation_tail += power * coefficients[power]
        else:
            relation_tail = sum(
                coefficients[index] * sums[power - index]
                for index in range(1, degree + 1)
            )
        sums.append(sp.expand(-relation_tail))
    return sums


def verify_factor_packet(
    simple_poly,
    double_poly,
    high_points,
    pole_parts,
    v_scale,
    g_scale,
    kappa,
    label,
):
    simple_poly = sp.Poly(simple_poly, x, domain=sp.QQ).monic().as_expr()
    double_poly = sp.Poly(double_poly, x, domain=sp.QQ).monic().as_expr()
    high_poly = sp.prod(x - point for point in high_points)
    denominator_p = sp.prod(
        (x - point) ** part
        for point, part in zip(high_points, pole_parts)
    )
    denominator_g = sp.prod(
        (x - point) ** (part + 1)
        for point, part in zip(high_points, pole_parts)
    )
    potential = sp.expand(v_scale * simple_poly * denominator_p * high_poly**2)
    base_g = cancel(g_scale * double_poly / denominator_g)
    function = cancel(potential * base_g**2)
    marked_poly = sp.expand(simple_poly * double_poly * high_poly)
    lam = sp.cancel(v_scale * g_scale)

    t_derivative_sum = sum(
        part * sp.cancel(high_poly / (x - point))
        for point, part in zip(high_points, pole_parts)
    )
    stieltjes = sp.expand(
        2 * simple_poly * high_poly * sp.diff(double_poly, x)
        + double_poly * high_poly * sp.diff(simple_poly, x)
        - double_poly * simple_poly * t_derivative_sum
    )
    first_order_coefficient = sp.expand(
        high_poly * sp.diff(simple_poly, x)
        - simple_poly * t_derivative_sum
    )
    heine_stieltjes_operator = sp.expand(
        2 * simple_poly * high_poly * sp.diff(double_poly, x, 2)
        + (
            2 * sp.diff(simple_poly * high_poly, x)
            + first_order_coefficient
        )
        * sp.diff(double_poly, x)
        + sp.diff(first_order_coefficient, x) * double_poly
    )
    gate(sp.degree(stieltjes, x) <= 0, f"{label}: Stieltjes polynomial constant")
    gate(
        heine_stieltjes_operator == 0,
        f"{label}: differentiated Heine-Stieltjes operator",
    )
    gate(
        cancel(lam * stieltjes - 2 * kappa) == 0,
        f"{label}: linear ODE scale",
    )
    gate(
        cancel(2 * potential * sp.diff(base_g, x) + sp.diff(potential, x) * base_g - 2 * kappa)
        == 0,
        f"{label}: 2VG'+V'G",
    )
    gate(
        cancel(sp.diff(function, x) - 2 * kappa * base_g) == 0,
        f"{label}: F'=2kappa G",
    )
    gate(
        cancel(potential * sp.diff(function, x) ** 2 - 4 * kappa**2 * function)
        == 0,
        f"{label}: square-potential equation",
    )
    gate(
        cancel(marked_poly * sp.diff(function, x) / function - stieltjes)
        == 0,
        f"{label}: exact constant-numerator Pade identity",
    )

    s = sp.degree(simple_poly, x)
    e = sp.degree(double_poly, x)
    h = len(high_points)
    total_degree = sum(pole_parts)
    r = s + e + h - 1
    gate(total_degree == s + 2 * e, f"{label}: N balance")
    gate(sp.degree(marked_poly, x) == r + 1, f"{label}: marked degree")

    simple_factor = sp.Poly(simple_poly, x, domain=sp.QQ)
    double_factor = sp.Poly(double_poly, x, domain=sp.QQ)
    gate(
        sp.gcd(simple_factor, simple_factor.diff()).degree() == 0,
        f"{label}: simple-root polynomial squarefree",
    )
    gate(
        sp.gcd(double_factor, double_factor.diff()).degree() == 0,
        f"{label}: double-zero location polynomial squarefree",
    )
    gate(
        sp.gcd(simple_factor, double_factor).degree() == 0,
        f"{label}: simple and double loci disjoint",
    )
    for point in high_points:
        gate(simple_poly.subs(x, point) != 0, f"{label}: simple/high loci disjoint")
        gate(double_poly.subs(x, point) != 0, f"{label}: double/high loci disjoint")
    gate(
        len(set(high_points)) == len(high_points),
        f"{label}: distinct high-root locations",
    )

    # Moment equations through the first active moment.
    simple_sums = root_power_sums(simple_poly, r)
    double_sums = root_power_sums(double_poly, r)
    for power_index in range(0, r + 1):
        moment = (
            simple_sums[power_index]
            + 2 * double_sums[power_index]
            - sum(
                part * point**power_index
                for point, part in zip(high_points, pole_parts)
            )
        )
        if power_index < r:
            gate(sp.simplify(moment) == 0, f"{label}: vanished moment {power_index}")
        else:
            gate(
                sp.simplify(lam * moment - 2 * kappa) == 0,
                f"{label}: first active moment",
            )

    # Exact half-square prefix: degree r+2, with a nonzero prescribed first
    # defect coefficient.
    square_prefix_defect = sp.expand(potential - v_scale * marked_poly**2)
    gate(
        sp.degree(square_prefix_defect, x) == r + 2,
        f"{label}: exact first square-prefix defect degree",
    )
    expected_lead = sp.cancel(2 * kappa * v_scale / (r * lam))
    actual_lead = sp.LC(sp.Poly(square_prefix_defect, x))
    gate(
        sp.simplify(actual_lead - expected_lead) == 0,
        f"{label}: exact first square-prefix coefficient",
    )

    # Pointwise interpolation/principal-part constraints.
    stieltjes_constant = sp.simplify(stieltjes)
    simple_interpolant = sp.expand(
        double_poly * high_poly * sp.diff(simple_poly, x)
        - stieltjes_constant
    )
    double_interpolant = sp.expand(
        2 * simple_poly * high_poly * sp.diff(double_poly, x)
        - stieltjes_constant
    )
    gate(
        sp.Poly(simple_interpolant, x, domain=sp.QQ)
        .rem(sp.Poly(simple_poly, x, domain=sp.QQ))
        .is_zero,
        f"{label}: simple-root interpolation",
    )
    gate(
        sp.Poly(double_interpolant, x, domain=sp.QQ)
        .rem(sp.Poly(double_poly, x, domain=sp.QQ))
        .is_zero,
        f"{label}: double-zero interpolation",
    )
    for point, part in zip(high_points, pole_parts):
        high_without = sp.cancel(high_poly / (x - point))
        gate(
            sp.simplify(
                -part
                * double_poly.subs(x, point)
                * simple_poly.subs(x, point)
                * high_without.subs(x, point)
                - stieltjes_constant
            )
            == 0,
            f"{label}: high-root principal part",
        )

    return function, potential, base_g


def family_controls():
    cyclic_rows = 0
    symmetric_rows = 0
    chord_rows = 0
    smallest_noncyclic = None

    for size in range(1, 11):
        # Cyclic e=0 family.
        simple_poly = x**size - 1
        v_scale = sp.Rational(1, size**2)
        g_scale = 2 * size
        function, potential, _ = verify_factor_packet(
            simple_poly,
            1,
            (0,),
            (size,),
            v_scale,
            g_scale,
            1,
            f"cyclic-{size}",
        )
        gate(
            cancel(function - 4 * (1 - x ** (-size))) == 0,
            f"cyclic-{size}: installed hostile formula",
        )
        cyclic_rows += 1

        if size >= 2:
            # Unique e=h=1 one-pole family.
            numerator = sp.expand(x**size - size * x + size - 1)
            quotient, remainder = sp.div(numerator, (x - 1) ** 2)
            gate(remainder == 0, f"symmetric-{size}: double-root quotient")
            simple_poly = quotient
            v_scale = sp.Rational(4, size**2 * (size - 1) ** 2)
            g_scale = sp.Rational(size * (size - 1), 2)
            function, potential, _ = verify_factor_packet(
                simple_poly,
                x - 1,
                (0,),
                (size,),
                v_scale,
                g_scale,
                1,
                f"symmetric-{size}",
            )
            expected_function = cancel(numerator / x**size)
            gate(
                cancel(function - expected_function) == 0,
                f"symmetric-{size}: explicit map",
            )
            scalar, factors = sp.factor_list(potential)
            square = (
                sp.sqrt(scalar).is_Rational
                and all(exponent % 2 == 0 for _, exponent in factors)
            )
            gate(
                square == (size == 2),
                f"symmetric-{size}: genuine nonsplit boundary",
            )
            if size == 3:
                smallest_noncyclic = (function, potential)
            symmetric_rows += 1

        if size >= 2:
            # Complete e=1,h=2 two-pole chord family.  Put b=N-d and
            # D=x^d(x-1)^b.  The unique off-pole critical point of D is
            # gamma=d/N, and c=D(gamma).  Hence F=1-c/D has exactly the
            # desired double zero and no other finite critical point.
            for distance in range(1, size // 2 + 1):
                complement = size - distance
                gamma = sp.Rational(distance, size)
                denominator = x**distance * (x - 1) ** complement
                critical_value = sp.expand(denominator.subs(x, gamma))
                numerator = sp.expand(denominator - critical_value)
                quotient, remainder = sp.div(numerator, (x - gamma) ** 2)
                gate(remainder == 0, f"chord-{size}-{distance}: double quotient")
                v_scale = sp.cancel(4 / (critical_value**2 * size**2))
                g_scale = sp.cancel(critical_value * size / 2)
                function, potential, _ = verify_factor_packet(
                    quotient,
                    x - gamma,
                    (0, 1),
                    (distance, complement),
                    v_scale,
                    g_scale,
                    1,
                    f"chord-{size}-{distance}",
                )
                expected_function = cancel(1 - critical_value / denominator)
                gate(
                    cancel(function - expected_function) == 0,
                    f"chord-{size}-{distance}: reciprocal-monomial map",
                )
                gate(
                    cancel(
                        sp.diff(function, x)
                        - critical_value
                        * size
                        * (x - gamma)
                        / (
                            x ** (distance + 1)
                            * (x - 1) ** (complement + 1)
                        )
                    )
                    == 0,
                    f"chord-{size}-{distance}: forced critical point",
                )
                scalar, factors = sp.factor_list(potential)
                square = (
                    sp.sqrt(scalar).is_Rational
                    and all(exponent % 2 == 0 for _, exponent in factors)
                )
                gate(not square, f"chord-{size}-{distance}: genuine nonsplit")
                chord_rows += 1

    gate(smallest_noncyclic is not None, "smallest noncyclic map missing")
    function, potential = smallest_noncyclic
    gate(
        cancel(function - (x - 1) ** 2 * (x + 2) / x**3) == 0,
        "smallest noncyclic F3",
    )
    gate(
        cancel(potential - x**5 * (x + 2) / 9) == 0,
        "smallest noncyclic V3",
    )
    return cyclic_rows, symmetric_rows, chord_rows


def main():
    census = {}
    e0_rows = 0
    eh1_rows = 0
    e1h2_rows = 0
    for size in range(1, 9):
        records = census_degree(size)
        genuine = [record for record in records if record["genuine"]]
        census[size] = {
            "all": len(records),
            "genuine": len(genuine),
            "passports": len(
                {
                    (record["e"], record["h"], record["p"])
                    for record in genuine
                }
            ),
            "noncyclic": sum(not record["cyclic"] for record in genuine),
        }

        e0 = [record for record in genuine if record["e"] == 0]
        gate(len(e0) == 1, f"N={size}: e=0 uniqueness")
        gate(e0[0]["h"] == 1 and e0[0]["p"] == (size,), f"N={size}: e=0 passport")
        gate(e0[0]["cyclic"] and e0[0]["group_order"] == size, f"N={size}: cyclic group")
        e0_rows += 1

        if size >= 3:
            eh1 = [
                record
                for record in genuine
                if record["e"] == 1 and record["h"] == 1
            ]
            gate(len(eh1) == 1, f"N={size}: e=h=1 uniqueness")
            gate(
                eh1[0]["group_order"] == factorial(size),
                f"N={size}: e=h=1 full symmetric group",
            )
            eh1_rows += 1

            e1h2 = [
                record
                for record in genuine
                if record["e"] == 1 and record["h"] == 2
            ]
            gate(len(e1h2) == size // 2, f"N={size}: e=1,h=2 chord census")
            expected_partitions = {
                tuple(sorted((distance, size - distance), reverse=True))
                for distance in range(1, size // 2 + 1)
            }
            gate(
                {record["p"] for record in e1h2} == expected_partitions,
                f"N={size}: e=1,h=2 pole partitions",
            )
            for record in e1h2:
                distance = min(record["p"])
                divisor = gcd(size, distance)
                expected_order = divisor * factorial(size // divisor) ** divisor
                gate(
                    record["group_order"] == expected_order,
                    f"N={size}: chord monodromy order",
                )
            e1h2_rows += len(e1h2)

    gate(census[2]["noncyclic"] == 0, "degree two must be cyclic")
    gate(census[3]["noncyclic"] == 2, "degree three first noncyclic census")

    (
        cyclic_family_rows,
        symmetric_family_rows,
        chord_family_rows,
    ) = family_controls()

    source = Path(__file__).read_text(encoding="utf-8")
    tree = ast.parse(source)
    assert_nodes = sum(isinstance(node, ast.Assert) for node in ast.walk(tree))
    gate(assert_nodes == 0, "truth-bearing Python assert")

    print("THM-2796 BALANCED RESPONSE STIELTJES-PADE CLASSIFICATION")
    for size in range(1, 9):
        row = census[size]
        print(
            f"N={size}:all_dessins={row['all']},genuine={row['genuine']},"
            f"passports={row['passports']},noncyclic={row['noncyclic']}"
        )
    print(f"e0_cyclic_census_rows={e0_rows}")
    print(f"eh1_full_symmetric_census_rows={eh1_rows}")
    print(f"e1h2_chord_census_rows={e1h2_rows}")
    print(f"cyclic_symbolic_family_rows={cyclic_family_rows}")
    print(f"symmetric_symbolic_family_rows={symmetric_family_rows}")
    print(f"chord_symbolic_family_rows={chord_family_rows}")
    print("smallest_genuine_noncyclic_degree=3")
    print("smallest_noncyclic_F=(x-1)^2(x+2)/x^3")
    print("smallest_noncyclic_V=x^5(x+2)/9")
    print("balanced_factor_Stieltjes_identity=PASS")
    print("differentiated_Heine_Stieltjes_operator=PASS")
    print("constant_numerator_Pade_identity=PASS")
    print("weighted_moments_through_rminus1=ZERO")
    print("first_active_moment=lambda^-1*2kappa")
    print("exact_square_prefix_defect_degree=r+2")
    print(f"assert_nodes={assert_nodes}")
    print(f"exact_gates={GATES}")
    print("FAILED CHECKS: NONE")


if __name__ == "__main__":
    main()
