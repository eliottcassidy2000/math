#!/usr/bin/env python3
"""Exact finite companion for provisional THM-3492.

The universal exact-sequence arguments are in the theorem.  This companion
checks bit-packed restriction maps, explicit kernel and target bases, direct
cyclic-arc images, every pointed Pascal carrier in the declared universes,
and an independently reconstructed physical depth-six hostile.  Every gate
is an explicit ``check`` call and remains active under ``python -O``.
"""

from __future__ import annotations

from functools import lru_cache


PERIOD_EXPONENT_MAX = 6
SLACK_DEGREE_MAX = 4


def check(condition: bool, label: str) -> None:
    if not condition:
        raise RuntimeError(f"CHECK FAILED: {label}")


def valuation_two(value: int) -> int:
    check(value > 0, "two-adic valuation has positive input")
    return (value & -value).bit_length() - 1


def submasks(mask: int) -> list[int]:
    values = []
    current = mask
    while True:
        values.append(current)
        if current == 0:
            return sorted(values)
        current = (current - 1) & mask


def binomial_mod_two(top: int, bottom: int) -> int:
    if top < 0 or bottom < 0 or bottom > top:
        return 0
    return int((bottom & ~top) == 0)


def y_power_profile(power: int) -> int:
    """X-coefficient mask of (X+1)^power over F_2."""

    return sum(1 << phase for phase in submasks(power))


def profile_evaluation(profile: int, phase: int) -> int:
    return (profile >> phase) & 1


def hasse_moment(profile: int, order: int, period: int) -> int:
    value = 0
    remaining = profile & ((1 << period) - 1)
    while remaining:
        bit = remaining & -remaining
        phase = bit.bit_length() - 1
        value ^= binomial_mod_two(phase, order)
        remaining ^= bit
    return value


def hasse_vector(profile: int, period: int) -> list[int]:
    return [hasse_moment(profile, order, period) for order in range(period)]


def polynomial_evaluation_one(polynomial: int) -> int:
    return polynomial.bit_count() & 1


def gf2_span_basis(vectors: list[int]) -> dict[int, int]:
    basis: dict[int, int] = {}
    for vector in vectors:
        current = vector
        while current:
            pivot = current.bit_length() - 1
            if pivot in basis:
                current ^= basis[pivot]
            else:
                basis[pivot] = current
                break
    return basis


def gf2_rank(vectors: list[int]) -> int:
    return len(gf2_span_basis(vectors))


def gf2_in_span(vector: int, basis: dict[int, int]) -> bool:
    current = vector
    while current:
        pivot = current.bit_length() - 1
        if pivot not in basis:
            return False
        current ^= basis[pivot]
    return True


def encode_target(profile: int, polynomial: int, period: int) -> int:
    return profile | (polynomial << period)


def encode_bivariate(slices: list[int], period: int) -> int:
    value = 0
    for exponent, profile in enumerate(slices):
        value ^= profile << (exponent * period)
    return value


def bivariate_restriction(
    slices: list[int], period: int, marked_phase: int
) -> tuple[int, int]:
    phase_axis = 0
    slack_axis = 0
    for exponent, profile in enumerate(slices):
        phase_axis ^= profile
        if profile_evaluation(profile, marked_phase):
            slack_axis ^= 1 << exponent
    return phase_axis, slack_axis


def compatible(profile: int, polynomial: int, marked_phase: int) -> bool:
    return profile_evaluation(profile, marked_phase) == polynomial_evaluation_one(
        polynomial
    )


def separated_lift(
    profile: int,
    polynomial: int,
    corner: int,
    carrier: int,
    degree: int,
) -> list[int]:
    """Return F+(D+c)g as q-coefficient profile slices."""

    check(polynomial.bit_length() <= degree + 1, "lift polynomial degree bound")
    slices = [0] * (degree + 1)
    slices[0] = profile
    reduced = polynomial ^ corner
    for exponent in range(degree + 1):
        if (reduced >> exponent) & 1:
            slices[exponent] ^= carrier
    return slices


def q_plus_one_tensor(profile: int, shift: int, degree: int) -> list[int]:
    check(0 <= shift < degree, "kernel shift in range")
    slices = [0] * (degree + 1)
    slices[shift] = profile
    slices[shift + 1] = profile
    return slices


def bivariate_hasse_moment(
    slices: list[int], order: int, period: int
) -> int:
    polynomial = 0
    for exponent, profile in enumerate(slices):
        if hasse_moment(profile, order, period):
            polynomial ^= 1 << exponent
    return polynomial


def profile_face(period: int, length: int) -> list[int]:
    check(0 < length < period, "proper profile face")
    endpoint = period - length
    return sorted(endpoint + item for item in submasks(length - 1))


def direct_arc_columns(period: int, length: int) -> list[int]:
    """Columns of q(h) -> xor_(0<=j<length) q(h+j)."""

    columns = []
    for input_phase in range(period):
        column = 0
        for offset in range(length):
            column ^= 1 << ((input_phase - offset) % period)
        columns.append(column)
    return columns


def audit_abstract_fiber() -> dict[str, int]:
    parameter_triples = 0
    restriction_columns = 0
    kernel_basis_vectors = 0
    target_basis_lifts = 0
    maximum_domain_dimension = 0
    maximum_kernel_dimension = 0

    for exponent in range(1, PERIOD_EXPONENT_MAX + 1):
        period = 1 << exponent
        constant_profile = (1 << period) - 1
        for marked_phase in range(period):
            check(
                profile_evaluation(constant_profile, marked_phase) == 1,
                f"unconstrained carrier p={period} h={marked_phase}",
            )
            for degree in range(1, SLACK_DEGREE_MAX + 1):
                domain_images = []
                for q_exponent in range(degree + 1):
                    for phase in range(period):
                        phase_axis = 1 << phase
                        slack_axis = (
                            1 << q_exponent if phase == marked_phase else 0
                        )
                        check(
                            compatible(phase_axis, slack_axis, marked_phase),
                            f"domain restriction compatibility p={period} h={marked_phase} N={degree}",
                        )
                        domain_images.append(
                            encode_target(phase_axis, slack_axis, period)
                        )
                expected_target_dimension = period + degree
                map_rank = gf2_rank(domain_images)
                check(
                    map_rank == expected_target_dimension,
                    f"abstract restriction rank p={period} h={marked_phase} N={degree}",
                )

                kernel_vectors = []
                for shift in range(degree):
                    for phase in range(period):
                        if phase == marked_phase:
                            continue
                        slices = q_plus_one_tensor(1 << phase, shift, degree)
                        check(
                            bivariate_restriction(slices, period, marked_phase)
                            == (0, 0),
                            f"abstract kernel restriction p={period} h={marked_phase} N={degree}",
                        )
                        kernel_vectors.append(encode_bivariate(slices, period))
                expected_kernel_dimension = degree * (period - 1)
                check(
                    gf2_rank(kernel_vectors) == expected_kernel_dimension,
                    f"abstract kernel basis p={period} h={marked_phase} N={degree}",
                )
                domain_dimension = (degree + 1) * period
                check(
                    domain_dimension - map_rank == expected_kernel_dimension,
                    f"abstract rank-nullity p={period} h={marked_phase} N={degree}",
                )

                fiber_basis: list[tuple[int, int]] = []
                for phase in range(period):
                    if phase != marked_phase:
                        fiber_basis.append((1 << phase, 0))
                fiber_basis.append((1 << marked_phase, 1))
                for q_exponent in range(1, degree + 1):
                    fiber_basis.append((0, 1 | (1 << q_exponent)))
                check(
                    len(fiber_basis) == expected_target_dimension,
                    f"abstract target basis count p={period} h={marked_phase} N={degree}",
                )
                check(
                    gf2_rank(
                        [
                            encode_target(profile, polynomial, period)
                            for profile, polynomial in fiber_basis
                        ]
                    )
                    == expected_target_dimension,
                    f"abstract target basis rank p={period} h={marked_phase} N={degree}",
                )
                for profile, polynomial in fiber_basis:
                    corner = profile_evaluation(profile, marked_phase)
                    check(
                        compatible(profile, polynomial, marked_phase),
                        f"abstract target compatibility p={period} h={marked_phase} N={degree}",
                    )
                    lift = separated_lift(
                        profile,
                        polynomial,
                        corner,
                        constant_profile,
                        degree,
                    )
                    check(
                        bivariate_restriction(lift, period, marked_phase)
                        == (profile, polynomial),
                        f"abstract separated lift p={period} h={marked_phase} N={degree}",
                    )

                parameter_triples += 1
                restriction_columns += len(domain_images)
                kernel_basis_vectors += len(kernel_vectors)
                target_basis_lifts += len(fiber_basis)
                maximum_domain_dimension = max(
                    maximum_domain_dimension, domain_dimension
                )
                maximum_kernel_dimension = max(
                    maximum_kernel_dimension, expected_kernel_dimension
                )

    return {
        "parameter_triples": parameter_triples,
        "restriction_columns": restriction_columns,
        "kernel_basis_vectors": kernel_basis_vectors,
        "target_basis_lifts": target_basis_lifts,
        "maximum_domain_dimension": maximum_domain_dimension,
        "maximum_kernel_dimension": maximum_kernel_dimension,
    }


def audit_arc_fiber() -> dict[str, int]:
    parameter_triples = 0
    direct_arc_matrices = 0
    pointed_carriers = 0
    restricted_kernel_vectors = 0
    restricted_target_lifts = 0
    two_carrier_hostiles = 0
    singleton_face_positive_kernels = 0
    maximum_constrained_kernel_dimension = 0

    for exponent in range(1, PERIOD_EXPONENT_MAX + 1):
        period = 1 << exponent
        for length in range(1, period):
            marked_phase = period - length
            radical_order = (1 << valuation_two(length)) - 1
            image_basis = [
                y_power_profile(order)
                for order in range(radical_order, period)
            ]
            image_span = gf2_span_basis(image_basis)
            image_dimension = period - radical_order
            check(
                len(image_span) == image_dimension,
                f"arc image ideal dimension p={period} s={length}",
            )

            arc_columns = direct_arc_columns(period, length)
            arc_span = gf2_span_basis(arc_columns)
            check(
                len(arc_span) == image_dimension,
                f"direct arc rank p={period} s={length}",
            )
            for column in arc_columns:
                check(
                    gf2_in_span(column, image_span),
                    f"direct arc column in ideal p={period} s={length}",
                )
            for vector in image_basis:
                check(
                    gf2_in_span(vector, arc_span),
                    f"ideal basis in direct arc image p={period} s={length}",
                )
            direct_arc_matrices += 1

            carrier = y_power_profile(period - 1)
            check(
                profile_evaluation(carrier, marked_phase) == 1,
                f"constant image carrier p={period} s={length}",
            )
            check(
                gf2_in_span(carrier, image_span),
                f"constant carrier in image p={period} s={length}",
            )
            image_kernel_basis = []
            for order in range(radical_order, period - 1):
                vector = y_power_profile(order)
                if profile_evaluation(vector, marked_phase):
                    vector ^= carrier
                check(
                    profile_evaluation(vector, marked_phase) == 0,
                    f"image mark-kernel vector p={period} s={length} j={order}",
                )
                check(
                    gf2_in_span(vector, image_span),
                    f"image mark-kernel membership p={period} s={length} j={order}",
                )
                image_kernel_basis.append(vector)
            check(
                gf2_rank(image_kernel_basis) == image_dimension - 1,
                f"image mark-kernel basis p={period} s={length}",
            )

            face = profile_face(period, length)
            check(
                len(face) == 1 << (length - 1).bit_count(),
                f"pointed face size p={period} s={length}",
            )
            check(
                min(face) >= radical_order + 1,
                f"pointed face transverse p={period} s={length}",
            )

            for degree in range(1, SLACK_DEGREE_MAX + 1):
                domain_images = []
                for q_exponent in range(degree + 1):
                    for vector in image_basis:
                        phase_axis = vector
                        slack_axis = (
                            1 << q_exponent
                            if profile_evaluation(vector, marked_phase)
                            else 0
                        )
                        check(
                            compatible(phase_axis, slack_axis, marked_phase),
                            f"constrained restriction compatibility p={period} s={length} N={degree}",
                        )
                        domain_images.append(
                            encode_target(phase_axis, slack_axis, period)
                        )
                expected_target_dimension = image_dimension + degree
                map_rank = gf2_rank(domain_images)
                check(
                    map_rank == expected_target_dimension,
                    f"constrained restriction rank p={period} s={length} N={degree}",
                )

                kernel_vectors = []
                for shift in range(degree):
                    for vector in image_kernel_basis:
                        slices = q_plus_one_tensor(vector, shift, degree)
                        check(
                            bivariate_restriction(
                                slices, period, marked_phase
                            )
                            == (0, 0),
                            f"constrained kernel restriction p={period} s={length} N={degree}",
                        )
                        kernel_vectors.append(encode_bivariate(slices, period))
                expected_kernel_dimension = degree * (image_dimension - 1)
                check(
                    gf2_rank(kernel_vectors) == expected_kernel_dimension,
                    f"constrained kernel basis p={period} s={length} N={degree}",
                )
                check(
                    (degree + 1) * image_dimension - map_rank
                    == expected_kernel_dimension,
                    f"constrained rank-nullity p={period} s={length} N={degree}",
                )

                fiber_basis: list[tuple[int, int]] = [
                    (vector, 0) for vector in image_kernel_basis
                ]
                fiber_basis.append((carrier, 1))
                for q_exponent in range(1, degree + 1):
                    fiber_basis.append((0, 1 | (1 << q_exponent)))
                check(
                    len(fiber_basis) == expected_target_dimension,
                    f"constrained target basis count p={period} s={length} N={degree}",
                )
                check(
                    gf2_rank(
                        [
                            encode_target(profile, polynomial, period)
                            for profile, polynomial in fiber_basis
                        ]
                    )
                    == expected_target_dimension,
                    f"constrained target basis rank p={period} s={length} N={degree}",
                )
                for profile, polynomial in fiber_basis:
                    corner = profile_evaluation(profile, marked_phase)
                    lift = separated_lift(
                        profile, polynomial, corner, carrier, degree
                    )
                    check(
                        bivariate_restriction(lift, period, marked_phase)
                        == (profile, polynomial),
                        f"constrained separated lift p={period} s={length} N={degree}",
                    )
                    for coefficient_profile in lift:
                        check(
                            gf2_in_span(coefficient_profile, image_span),
                            f"constrained lift coefficient in image p={period} s={length} N={degree}",
                        )

                test_profile = y_power_profile(radical_order)
                corner = profile_evaluation(test_profile, marked_phase)
                test_polynomial = (
                    corner ^ (1 << (degree - 1)) ^ (1 << degree)
                )
                check(
                    polynomial_evaluation_one(test_polynomial) == corner,
                    f"pointed-route test compatibility p={period} s={length} N={degree}",
                )
                routes = []
                for face_order in face:
                    pointed_carrier = y_power_profile(face_order)
                    check(
                        gf2_in_span(pointed_carrier, image_span),
                        f"pointed carrier in image p={period} s={length} j={face_order}",
                    )
                    check(
                        profile_evaluation(pointed_carrier, marked_phase) == 1,
                        f"pointed carrier mark p={period} s={length} j={face_order}",
                    )
                    moments = hasse_vector(pointed_carrier, period)
                    check(
                        moments
                        == [
                            int(order == face_order)
                            for order in range(period)
                        ],
                        f"pointed carrier singleton Hasse p={period} s={length} j={face_order}",
                    )
                    route = separated_lift(
                        test_profile,
                        test_polynomial,
                        corner,
                        pointed_carrier,
                        degree,
                    )
                    check(
                        bivariate_restriction(route, period, marked_phase)
                        == (test_profile, test_polynomial),
                        f"pointed route axes p={period} s={length} N={degree} j={face_order}",
                    )
                    for coefficient_profile in route:
                        check(
                            gf2_in_span(coefficient_profile, image_span),
                            f"pointed route coefficient in image p={period} s={length} N={degree} j={face_order}",
                        )
                    reduced = test_polynomial ^ corner
                    for order in range(period):
                        expected = hasse_moment(
                            test_profile, order, period
                        ) ^ 0
                        expected_polynomial = expected
                        if order == face_order:
                            expected_polynomial ^= reduced
                        check(
                            bivariate_hasse_moment(route, order, period)
                            == expected_polynomial,
                            f"pointed route Hasse table p={period} s={length} N={degree} j={face_order} r={order}",
                        )
                    check(
                        (
                            bivariate_hasse_moment(
                                route, face_order, period
                            )
                            >> degree
                        )
                        & 1
                        == 1,
                        f"pointed route monic top p={period} s={length} N={degree} j={face_order}",
                    )
                    routes.append(route)
                    pointed_carriers += 1

                if len(routes) >= 2:
                    check(
                        routes[0] != routes[1],
                        f"two pointed routes distinct p={period} s={length} N={degree}",
                    )
                    difference = [
                        left ^ right
                        for left, right in zip(routes[0], routes[1])
                    ]
                    check(
                        bivariate_restriction(
                            difference, period, marked_phase
                        )
                        == (0, 0),
                        f"two pointed routes kernel p={period} s={length} N={degree}",
                    )
                    two_carrier_hostiles += 1
                else:
                    check(
                        expected_kernel_dimension > 0,
                        f"singleton face retains mixed kernel p={period} s={length} N={degree}",
                    )
                    singleton_face_positive_kernels += 1

                parameter_triples += 1
                restricted_kernel_vectors += len(kernel_vectors)
                restricted_target_lifts += len(fiber_basis)
                maximum_constrained_kernel_dimension = max(
                    maximum_constrained_kernel_dimension,
                    expected_kernel_dimension,
                )

    return {
        "parameter_triples": parameter_triples,
        "direct_arc_matrices": direct_arc_matrices,
        "pointed_carriers": pointed_carriers,
        "restricted_kernel_vectors": restricted_kernel_vectors,
        "restricted_target_lifts": restricted_target_lifts,
        "two_carrier_hostiles": two_carrier_hostiles,
        "singleton_face_positive_kernels": singleton_face_positive_kernels,
        "maximum_constrained_kernel_dimension": maximum_constrained_kernel_dimension,
    }


def rule30(left: int, center: int, right: int) -> int:
    return left ^ center ^ right ^ (center & right)


def centered_rows(horizon: int) -> list[frozenset[int]]:
    rows = []
    row = frozenset({0})
    for time in range(horizon + 1):
        rows.append(row)
        row = frozenset(
            site
            for site in range(-time - 1, time + 2)
            if rule30(
                int(site - 1 in row),
                int(site in row),
                int(site + 1 in row),
            )
        )
    return rows


def folded_source(
    rows: list[frozenset[int]], source_depth: int, distance: int
) -> int:
    source_time = source_depth + distance
    if source_time < 2:
        return 0
    row = rows[source_time]

    def edge(site: int) -> int:
        return int(site in row and site + 1 in row)

    if distance == 0:
        return edge(0)
    return edge(distance) ^ edge(-distance)


@lru_cache(maxsize=None)
def ternary_coefficient(exponent: int, degree: int) -> int:
    """[x^degree](1+x+x^2)^exponent modulo two by binary recursion."""

    if exponent < 0 or degree < 0 or degree > 2 * exponent:
        return 0
    if exponent == 0:
        return int(degree == 0)
    half_exponent = exponent // 2
    half_degree = degree // 2
    if exponent % 2 == 0:
        if degree % 2:
            return 0
        return ternary_coefficient(half_exponent, half_degree)
    if degree % 2:
        return ternary_coefficient(half_exponent, half_degree)
    return ternary_coefficient(
        half_exponent, half_degree
    ) ^ ternary_coefficient(half_exponent, half_degree - 1)


def inward_target_polynomial(
    rows: list[frozenset[int]], target: int
) -> int:
    polynomial = 0
    for source_depth in range(3, target):
        remainder = target - source_depth - 1
        for distance in range(remainder // 2 + 1):
            slack = remainder - 2 * distance
            source = folded_source(rows, source_depth, distance)
            transport = ternary_coefficient(distance + slack, slack)
            if source & transport:
                polynomial ^= 1 << slack
    return polynomial


def phi_width(row: int, width: int) -> int:
    return (row ^ ((row << 1) | (row << 2))) & ((1 << width) - 1)


def seed_cycle(width: int) -> list[int]:
    rows = []
    seen: set[int] = set()
    row = 1
    while row not in seen:
        seen.add(row)
        rows.append(row)
        row = phi_width(row, width)
    check(row == 1, f"seed cycle closes width={width}")
    return rows


def orbit_prefix(width: int, horizon: int) -> list[int]:
    rows = []
    row = 1
    for _ in range(horizon + 1):
        rows.append(row)
        row = phi_width(row, width)
    return rows


def arc_profile(current: list[int], length: int) -> int:
    period = len(current)
    profile = 0
    for phase in range(period):
        value = 0
        for offset in range(length):
            value ^= current[(phase + offset) % period]
        profile ^= value << phase
    return profile


def tuple_bits(value: int, width: int) -> tuple[int, ...]:
    return tuple((value >> index) & 1 for index in range(width))


def audit_physical_depth_six() -> dict[str, object]:
    target = 6
    rows = centered_rows(target)
    packed_direct = 0
    for inward in range(2 * target + 1):
        if target - inward in rows[target]:
            packed_direct ^= 1 << inward

    period = len(seed_cycle(target))
    check(period == 8, "physical P_6=8")
    orbit = orbit_prefix(target + 1, target + period)
    check(
        orbit[target] == packed_direct & ((1 << (target + 1)) - 1),
        "packed/local Rule 30 row at depth six",
    )
    center = (orbit[target] >> target) & 1
    length = target % period
    marked_phase = period - length
    radical_order = (1 << valuation_two(length)) - 1
    degree = target - 4
    check(
        (length, marked_phase, radical_order, degree) == (6, 2, 1, 2),
        "physical depth-six parameters",
    )

    current = []
    for phase in range(period):
        row = orbit[phase + target]
        current.append(
            ((row >> (target - 1)) & 1)
            | ((row >> (target - 2)) & 1)
        )
    terminal_profile = arc_profile(current, target)
    check(
        profile_evaluation(terminal_profile, marked_phase) == center,
        "physical depth-six marked phase",
    )
    terminal_moments = hasse_vector(terminal_profile, period)

    inward_polynomial = inward_target_polynomial(rows, target)
    backbone = ternary_coefficient(target - 2, target - 1)
    completed_polynomial = inward_polynomial ^ backbone
    check(backbone == 0, "physical B_6=0")
    check(
        inward_polynomial == ((1 << 0) | (1 << 2)),
        "physical C_6=1+q^2",
    )
    check(
        completed_polynomial == inward_polynomial,
        "physical D_6=C_6",
    )
    check(
        polynomial_evaluation_one(completed_polynomial) == center == 0,
        "physical depth-six common corner",
    )
    check(
        completed_polynomial.bit_length() - 1 == degree
        and ((completed_polynomial >> degree) & 1) == 1,
        "physical depth-six monic slack",
    )

    face = profile_face(period, length)
    check(face == [2, 3, 6, 7], "physical depth-six pointed face")
    check(
        sum(terminal_moments[order] for order in face) & 1 == center,
        "physical depth-six Pascal point",
    )
    check(
        [(order, terminal_moments[order]) for order in face]
        == [(2, 0), (3, 1), (6, 0), (7, 1)],
        "physical depth-six face values",
    )

    image_basis = [
        y_power_profile(order) for order in range(radical_order, period)
    ]
    image_span = gf2_span_basis(image_basis)
    routes = []
    for face_order in (2, 3):
        carrier = y_power_profile(face_order)
        check(
            gf2_in_span(carrier, image_span),
            f"physical carrier g_{face_order} in YV_8",
        )
        check(
            profile_evaluation(carrier, marked_phase) == 1,
            f"physical carrier g_{face_order} mark",
        )
        check(
            hasse_vector(carrier, period)
            == [int(order == face_order) for order in range(period)],
            f"physical carrier g_{face_order} singleton Hasse",
        )
        route = separated_lift(
            terminal_profile,
            completed_polynomial,
            center,
            carrier,
            degree,
        )
        check(
            bivariate_restriction(route, period, marked_phase)
            == (terminal_profile, completed_polynomial),
            f"physical Z_{face_order} axes",
        )
        for coefficient_profile in route:
            check(
                gf2_in_span(coefficient_profile, image_span),
                f"physical Z_{face_order} coefficient in YV_8",
            )
        for order in range(period):
            expected = terminal_moments[order]
            if order == face_order:
                expected ^= completed_polynomial
            check(
                bivariate_hasse_moment(route, order, period) == expected,
                f"physical Z_{face_order} Hasse routing order={order}",
            )
        routes.append(route)

    check(routes[0] != routes[1], "physical Z_2 and Z_3 distinct")
    route_difference = [
        left ^ right for left, right in zip(routes[0], routes[1])
    ]
    check(
        bivariate_restriction(route_difference, period, marked_phase)
        == (0, 0),
        "physical Z_2+Z_3 is mixed kernel",
    )
    constrained_kernel_dimension = degree * (
        period - radical_order - 1
    )
    check(constrained_kernel_dimension == 12, "physical mixed kernel dimension")

    return {
        "period": period,
        "length": length,
        "marked_phase": marked_phase,
        "radical_order": radical_order,
        "degree": degree,
        "center": center,
        "backbone": backbone,
        "inward_polynomial": inward_polynomial,
        "terminal_profile": tuple_bits(terminal_profile, period),
        "terminal_moments": tuple(terminal_moments),
        "face": tuple(face),
        "face_values": tuple(
            (order, terminal_moments[order]) for order in face
        ),
        "kernel_dimension": constrained_kernel_dimension,
    }


def main() -> None:
    abstract = audit_abstract_fiber()
    arc = audit_arc_fiber()
    physical = audit_physical_depth_six()

    check(abstract["parameter_triples"] == 504, "declared abstract universe")
    check(arc["parameter_triples"] == 480, "declared arc universe")
    check(arc["direct_arc_matrices"] == 120, "declared direct arc matrices")

    print("THM-3492 RULE 30 SLACK-PHASE FIBER EXACT AUDIT")
    print(
        "abstract_universe="
        "p=2^d,d=1..6;all_marks;N=1..4;"
        f"parameter_triples={abstract['parameter_triples']}"
    )
    print(
        "abstract_linear_algebra="
        f"restriction_columns={abstract['restriction_columns']};"
        f"kernel_basis_vectors={abstract['kernel_basis_vectors']};"
        f"target_basis_lifts={abstract['target_basis_lifts']};"
        f"max_domain_dimension={abstract['maximum_domain_dimension']};"
        f"max_kernel_dimension={abstract['maximum_kernel_dimension']}"
    )
    print(
        "arc_universe="
        "p=2^d,d=1..6;1<=s<p;N=1..4;"
        f"parameter_triples={arc['parameter_triples']};"
        f"direct_arc_matrices={arc['direct_arc_matrices']}"
    )
    print(
        "arc_linear_algebra="
        f"restricted_kernel_vectors={arc['restricted_kernel_vectors']};"
        f"restricted_target_lifts={arc['restricted_target_lifts']};"
        f"max_constrained_kernel_dimension={arc['maximum_constrained_kernel_dimension']}"
    )
    print(
        "pointed_routing="
        f"carriers={arc['pointed_carriers']};"
        f"two_carrier_hostiles={arc['two_carrier_hostiles']};"
        f"singleton_face_positive_kernels={arc['singleton_face_positive_kernels']}"
    )
    print(
        "physical_depth6="
        f"p={physical['period']};s={physical['length']};"
        f"hstar={physical['marked_phase']};ell={physical['radical_order']};"
        f"N={physical['degree']};c={physical['center']};B={physical['backbone']};"
        f"D_exponents={(0, 2)};kernel_dimension={physical['kernel_dimension']}"
    )
    print(f"physical_F6={physical['terminal_profile']}")
    print(f"physical_M6={physical['terminal_moments']}")
    print(
        f"physical_point_face={physical['face']};"
        f"values={physical['face_values']};routes=(2,3)"
    )
    print("ordinary_mode_and_optimized_mode_use_identical_explicit_gates")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
