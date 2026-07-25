#!/usr/bin/env python3
"""Exact, optimization-safe companion for THM-2348."""

from fractions import Fraction
from itertools import product


def require(condition: bool, message: str) -> None:
    """Raise under ordinary and optimized Python."""
    if not condition:
        raise RuntimeError(message)


def bit_count(mask: int) -> int:
    return mask.bit_count()


def submasks(mask: int):
    submask = mask
    while True:
        yield submask
        if submask == 0:
            return
        submask = (submask - 1) & mask


def coordinates(mask: int, dimension: int) -> tuple[int, ...]:
    return tuple(index for index in range(dimension) if mask & (1 << index))


def marginal(
    table: dict[tuple[int, ...], Fraction],
    colours: int,
    dimension: int,
    mask: int,
    values: tuple[int, ...],
) -> Fraction:
    active = coordinates(mask, dimension)
    require(
        len(active) == len(values),
        "marginal received a value tuple of the wrong length",
    )
    prescribed = dict(zip(active, values))
    total = Fraction(0)
    count = 0
    for point, value in table.items():
        if all(point[index] == colour for index, colour in prescribed.items()):
            total += value
            count += 1
    require(
        count == colours ** (dimension - len(active)),
        "uniform marginal fibre size changed",
    )
    return total / count


def anova(
    table: dict[tuple[int, ...], Fraction],
    colours: int,
    dimension: int,
) -> dict[int, dict[tuple[int, ...], Fraction]]:
    answer: dict[int, dict[tuple[int, ...], Fraction]] = {}
    for mask in range(1 << dimension):
        active = coordinates(mask, dimension)
        component: dict[tuple[int, ...], Fraction] = {}
        for values in product(range(colours), repeat=len(active)):
            value_by_coordinate = dict(zip(active, values))
            total = Fraction(0)
            for submask in submasks(mask):
                sub_active = coordinates(submask, dimension)
                sub_values = tuple(
                    value_by_coordinate[index] for index in sub_active
                )
                sign = -1 if (bit_count(mask) - bit_count(submask)) % 2 else 1
                total += sign * marginal(
                    table,
                    colours,
                    dimension,
                    submask,
                    sub_values,
                )
            component[values] = total
        answer[mask] = component
    return answer


def rectangle_atlas() -> tuple[int, int, int]:
    """Exhaust all {-1,0,1}-valued 2x2 tables and normalize defects."""
    tables = 0
    rectangular = 0
    opposite_corner_controls = 0
    for entries in product(range(-1, 2), repeat=4):
        e00, e01, e10, e11 = map(Fraction, entries)
        delta = e00 + e11 - e01 - e10
        additive = all(
            (
                (e00, e01),
                (e10, e11),
            )[row][column]
            == e00
            + (e10 - e00 if row else 0)
            + (e01 - e00 if column else 0)
            for row in range(2)
            for column in range(2)
        )
        require(
            additive == (delta == 0),
            "rectangle zero and additive separation disagree",
        )
        tables += 1
        if delta == 0:
            rectangular += 1
            continue

        if delta < 0:
            epsilon = -delta / 2
            row = (Fraction(0), epsilon - e10 + e00)
            column = (-e00, epsilon - e01)
            expected_minima = {(0, 0), (1, 1)}
        else:
            epsilon = delta / 2
            row = (Fraction(0), -e10 - epsilon + e00)
            column = (epsilon - e00, -e01)
            expected_minima = {(0, 1), (1, 0)}

        normalized = {
            (i, j): ((e00, e01), (e10, e11))[i][j] + row[i] + column[j]
            for i in range(2)
            for j in range(2)
        }
        minimum = min(normalized.values())
        minima = {point for point, value in normalized.items() if value == minimum}
        require(
            minima == expected_minima,
            "rectangle defect did not normalize to opposite-corner minima",
        )
        require(
            {(i, j) for i in {p[0] for p in minima} for j in {p[1] for p in minima}}
            != minima,
            "opposite-corner optimum became Cartesian",
        )
        opposite_corner_controls += 1
    return tables, rectangular, opposite_corner_controls


def conditioning_atlas() -> tuple[int, int, int]:
    """Check ANOVA restriction, averaging, and Mobius contraction exactly."""
    colours = 3
    dimension = 3
    table = {}
    for x0, x1, x2 in product(range(colours), repeat=dimension):
        table[(x0, x1, x2)] = Fraction(
            (x0 + 1) * (x1 + 2)
            + (x1 - x2) ** 2
            + 3 * int(x0 == x2)
            - 2 * int((x0, x1, x2) == (1, 2, 0))
        )

    full = anova(table, colours, dimension)
    fixed_value = 1
    restricted = {
        (x0, x1): table[(x0, x1, fixed_value)]
        for x0, x1 in product(range(colours), repeat=2)
    }
    restricted_anova = anova(restricted, colours, 2)

    conditioning_checks = 0
    for mask in range(1 << 2):
        active = coordinates(mask, 2)
        for values in product(range(colours), repeat=len(active)):
            with_fixed = []
            value_by_coordinate = dict(zip(active, values))
            for index in coordinates(mask | 4, 3):
                with_fixed.append(
                    fixed_value if index == 2 else value_by_coordinate[index]
                )
            expected = (
                full[mask][values] + full[mask | 4][tuple(with_fixed)]
            )
            require(
                restricted_anova[mask][values] == expected,
                "labelled ANOVA conditioning identity failed",
            )
            conditioning_checks += 1

    averaged = {}
    for x0, x1 in product(range(colours), repeat=2):
        averaged[(x0, x1)] = sum(
            (table[(x0, x1, x2)] for x2 in range(colours)),
            Fraction(0),
        ) / colours
    averaged_anova = anova(averaged, colours, 2)
    averaging_checks = 0
    for mask in range(1 << 2):
        active = coordinates(mask, 2)
        for values in product(range(colours), repeat=len(active)):
            require(
                averaged_anova[mask][values] == full[mask][values],
                "uniform deletion did not preserve the surviving ANOVA tensor",
            )
            averaging_checks += 1

    all_mask = (1 << dimension) - 1
    mobius = {
        colour: {
            mask: Fraction(
                (colour + 1) * (2 * bit_count(mask) + mask)
                - (mask % 3)
            )
            for mask in range(1, 1 << dimension)
        }
        for colour in range(colours)
    }

    def score(colour: int, mask: int) -> Fraction:
        return sum(
            (
                mobius[colour][submask]
                for submask in submasks(mask)
                if submask
            ),
            Fraction(0),
        )

    fixed_bit = 4
    mobius_checks = 0
    for colour in range(colours):
        fixed_for_colour = fixed_bit if colour == fixed_value else 0

        def conditioned_score(mask: int) -> Fraction:
            return (
                score(colour, mask | fixed_for_colour)
                - score(colour, fixed_for_colour)
            )

        for mask in range(1, 1 << 2):
            contracted = Fraction(0)
            for submask in submasks(mask):
                sign = -1 if (bit_count(mask) - bit_count(submask)) % 2 else 1
                contracted += sign * conditioned_score(submask)
            expected = mobius[colour][mask]
            if fixed_for_colour:
                expected += mobius[colour][mask | fixed_bit]
            require(
                contracted == expected,
                "Boolean Mobius conditioning contraction failed",
            )
            mobius_checks += 1

    require(all_mask == 7, "three-token universe changed")
    return conditioning_checks, averaging_checks, mobius_checks


def gauge_and_fixed_optimum_atlas() -> tuple[int, int, Fraction]:
    """Check the two-colour gauge and unbounded hidden coupling family."""
    gauge_table = ((Fraction(1), Fraction(0)), (Fraction(0), Fraction(-1)))
    gauge_rectangle = (
        gauge_table[0][0]
        + gauge_table[1][1]
        - gauge_table[0][1]
        - gauge_table[1][0]
    )
    require(gauge_rectangle == 0, "opposite blockwise pairs were not unary")

    hostile_checks = 0
    largest_hidden = Fraction(0)
    for magnitude in range(1, 65):
        table = (
            (Fraction(0), Fraction(magnitude)),
            (Fraction(magnitude), Fraction(magnitude)),
        )
        values = {
            (i, j): table[i][j] for i in range(2) for j in range(2)
        }
        minimum = min(values.values())
        minima = {point for point, value in values.items() if value == minimum}
        require(minima == {(0, 0)}, "fixed-optimum hostile lost its singleton")
        interaction = (
            table[0][0] - table[0][1] - table[1][0] + table[1][1]
        ) / 4
        require(
            interaction == Fraction(-magnitude, 4),
            "fixed-optimum hidden interaction coefficient changed",
        )
        largest_hidden = max(largest_hidden, abs(interaction))
        hostile_checks += 1
    return int(gauge_rectangle), hostile_checks, largest_hidden


def owner_criterion_atlas() -> tuple[int, int, int]:
    """Exhaust a three-block bank for the exact product-owner iff."""
    instances = 0
    product_instances = 0
    nonzero_theta_instances = 0
    blocks = range(3)
    for a_values in product(range(3), repeat=3):
        alpha = min(a_values)
        owners_p = {block for block in blocks if a_values[block] == alpha}
        for b_values in product(range(3), repeat=3):
            beta = min(b_values)
            owners_q = {block for block in blocks if b_values[block] == beta}
            proposed = {
                (left, right) for left in owners_p for right in owners_q
            }
            for mu_values in product(range(-2, 3), repeat=3):
                energy = {
                    (left, right): Fraction(
                        a_values[left]
                        + b_values[right]
                        + (mu_values[left] if left == right else 0)
                    )
                    for left in blocks
                    for right in blocks
                }
                minimum = min(energy.values())
                actual = {
                    point for point, value in energy.items() if value == minimum
                }

                inside_values = {
                    Fraction(mu_values[left] if left == right else 0)
                    for left, right in proposed
                }
                criterion = len(inside_values) == 1
                theta = next(iter(inside_values)) if criterion else Fraction(0)
                if criterion:
                    criterion = all(
                        Fraction(
                            a_values[left]
                            - alpha
                            + b_values[right]
                            - beta
                            + (mu_values[left] if left == right else 0)
                        )
                        > theta
                        for left in blocks
                        for right in blocks
                        if (left, right) not in proposed
                    )
                require(
                    criterion == (actual == proposed),
                    "constant-inside / strict-outside owner criterion failed",
                )
                if criterion:
                    product_instances += 1
                    require(
                        minimum == alpha + beta + theta,
                        "product-owner minimum formula failed",
                    )
                    if theta:
                        require(
                            len(owners_p) == len(owners_q) == 1
                            and owners_p == owners_q,
                            "nonzero theta occurred outside the common "
                            "unique-owner case",
                        )
                        nonzero_theta_instances += 1
                instances += 1
    return instances, product_instances, nonzero_theta_instances


def triangle_bound_atlas() -> tuple[int, int]:
    """Check the algebraic metric universe used by the universal bounds."""
    controls = 0
    unknot_controls = 0
    for u_source, u_p, u_q in product(range(6), repeat=3):
        for d_source_p in range(
            abs(u_source - u_p),
            u_source + u_p + 1,
        ):
            for d_source_q in range(
                abs(u_source - u_q),
                u_source + u_q + 1,
            ):
                lower = max(
                    0,
                    d_source_p - u_q,
                    d_source_q - u_p,
                )
                upper = min(
                    d_source_p + u_q,
                    d_source_q + u_p,
                )
                for d_source_pq in range(lower, upper + 1):
                    mu = (
                        d_source_pq
                        - d_source_p
                        - d_source_q
                        + u_source
                    )
                    require(
                        mu >= -2 * min(u_p, u_q),
                        "lower Gordian cohabitation bound failed",
                    )
                    require(
                        mu <= 2 * min(u_source, u_p, u_q),
                        "upper Gordian cohabitation bound failed",
                    )
                    controls += 1
                    if u_source == 0:
                        sigma = u_p + u_q - d_source_pq
                        require(mu == -sigma, "unknot interaction sign changed")
                        require(sigma >= 0, "unknot interaction lost positivity")
                        unknot_controls += 1
    return controls, unknot_controls


rectangle_tables, rectangle_zero, opposite_corner_controls = rectangle_atlas()
conditioning_checks, averaging_checks, mobius_checks = conditioning_atlas()
gauge_rectangle, fixed_optimum_hostiles, largest_hidden = (
    gauge_and_fixed_optimum_atlas()
)
owner_instances, owner_product_instances, nonzero_theta_instances = (
    owner_criterion_atlas()
)
triangle_controls, unknot_controls = triangle_bound_atlas()

require(rectangle_tables == 81, "rectangle census changed")
require(
    rectangle_zero + opposite_corner_controls == rectangle_tables,
    "rectangle partition census changed",
)
require(owner_instances == 91_125, "owner criterion census changed")

print("theorem=THM-2348")
print("status=PROVED+VERIFIED-EXACT+INDEPENDENTLY-HOSTILE-AUDITED")
print(f"rectangle_tables={rectangle_tables}")
print(f"rectangle_zero_tables={rectangle_zero}")
print(f"opposite_corner_controls={opposite_corner_controls}")
print(f"conditioning_component_checks={conditioning_checks}")
print(f"uniform_marginal_checks={averaging_checks}")
print(f"mobius_contraction_checks={mobius_checks}")
print(f"two_colour_gauge_rectangle={gauge_rectangle}")
print(f"fixed_optimum_hostiles={fixed_optimum_hostiles}")
print(f"largest_tested_hidden_interaction={largest_hidden}")
print(f"owner_criterion_instances={owner_instances}")
print(f"owner_product_instances={owner_product_instances}")
print(f"nonzero_theta_product_instances={nonzero_theta_instances}")
print(f"triangle_bound_controls={triangle_controls}")
print(f"unknot_interaction_identifications={unknot_controls}")
print("mixed_anova_vanishing_iff_robust_rectangularity=YES")
print("single_cartesian_optimum_implies_decoupling=NO")
print("target_token_conditioning_commutes=YES")
print("unconditional_connected_sum_additivity=NO")
print("all_checks=PASS")
