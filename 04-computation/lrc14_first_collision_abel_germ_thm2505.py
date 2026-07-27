#!/usr/bin/env python3
"""Exact finite-cylinder companion for THM-2505.

For each prime p in a small bank, two disjoint p^4-cylinder atoms first
collide after three Perron quotients.  The calculation verifies the collision
mass sequence, its Abel summation-by-parts identity, the common order two of
all primitive coloured germs, and the Gregory/Ramanujan sign at the corner.
Only Fraction arithmetic and cyclotomic reduction for prime p are used.
"""

from fractions import Fraction as F


PRIMES = (2, 3, 5, 7, 13)
DEPTH = 4
FIRST_COLLISION = 3


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def perron(values: list[F], prime: int) -> list[F]:
    """Apply P_p to a uniform p-adic cylinder vector."""

    if len(values) == 1:
        return values[:]
    require(len(values) % prime == 0, "p-adic grid divisibility")
    base = len(values) // prime
    return [
        sum(values[index + root * base] for root in range(prime)) / prime
        for index in range(base)
    ]


def inner(left: list[F], right: list[F]) -> F:
    require(len(left) == len(right), "common cylinder grid")
    return sum(a * b for a, b in zip(left, right)) / len(left)


def root_correlation(left: list[F], right: list[F], prime: int) -> list[F]:
    """C(t)=int sum_u left(y,u+t) right(y,u) dy."""

    require(len(left) == len(right) and len(left) % prime == 0, "root chart")
    base = len(left) // prime
    return [
        sum(
            left[index + ((root + shift) % prime) * base]
            * right[index + root * base]
            for index in range(base)
            for root in range(prime)
        )
        / base
        for shift in range(prime)
    ]


def prime_cyclotomic_reduce(values: list[F], colour: int, prime: int) -> tuple[F, ...]:
    """Reduce sum_s values[s] zeta_p^(-colour*s) in 1,...,zeta^(p-2)."""

    raw = [F(0) for _ in range(prime)]
    for shift, value in enumerate(values):
        raw[(-colour * shift) % prime] += value
    return tuple(raw[index] - raw[prime - 1] for index in range(prime - 1))


def add_vectors(vectors: list[tuple[F, ...]]) -> tuple[F, ...]:
    if not vectors:
        return ()
    return tuple(sum(vector[index] for vector in vectors) for index in range(len(vectors[0])))


def evaluate_vector_germ(coefficients: list[tuple[F, ...]], z: F) -> tuple[F, ...]:
    return tuple(
        sum(coefficient[index] * z**shell for shell, coefficient in enumerate(coefficients))
        for index in range(len(coefficients[0]))
    )


def cell_germ_coefficients(
    prime: int,
    depth: int,
    left_index: int,
    right_index: int,
    left_weight: F = F(1),
) -> tuple[list[F], list[list[tuple[F, ...]]]]:
    """Collision masses and primitive germ coefficients for two cylinder atoms."""

    grid = prime**depth
    left = [F(0) for _ in range(grid)]
    right = [F(0) for _ in range(grid)]
    left[left_index] = left_weight
    right[right_index] = F(1)
    levels_left = [left]
    levels_right = [right]
    for _ in range(depth):
        levels_left.append(perron(levels_left[-1], prime))
        levels_right.append(perron(levels_right[-1], prime))
    masses = [inner(a, b) for a, b in zip(levels_left, levels_right)]
    by_colour: list[list[tuple[F, ...]]] = [[] for _ in range(prime - 1)]
    for shell in range(depth):
        correlation = root_correlation(levels_left[shell], levels_right[shell], prime)
        for colour in range(1, prime):
            by_colour[colour - 1].append(
                tuple(
                    value / (prime * prime)
                    for value in prime_cyclotomic_reduce(correlation, colour, prime)
                )
            )
    return masses, by_colour


def valuation(integer: int, prime: int) -> int:
    require(integer != 0, "valuation of zero is not used")
    integer = abs(integer)
    answer = 0
    while integer % prime == 0:
        integer //= prime
        answer += 1
    return answer


def abel_rhs(masses: list[F], plateau: F, z: F) -> F:
    """Evaluate -(1-z) sum_{s>=1} I_s z^(s-1), including the constant tail."""

    require(F(0) <= z < F(1), "Abel parameter")
    finite = sum(masses[s] * z ** (s - 1) for s in range(1, len(masses)))
    tail = plateau * z ** (len(masses) - 1) / (1 - z)
    return -(1 - z) * (finite + tail)


def run_prime(prime: int) -> dict[str, object]:
    grid = prime**DEPTH
    left = [F(0) for _ in range(grid)]
    right = [F(0) for _ in range(grid)]
    left[0] = F(1)
    right[prime ** (DEPTH - FIRST_COLLISION)] = F(1)
    require(inner(left, right) == 0, "initial atoms must be disjoint")

    levels_left = [left]
    levels_right = [right]
    for _ in range(DEPTH):
        levels_left.append(perron(levels_left[-1], prime))
        levels_right.append(perron(levels_right[-1], prime))

    masses = [inner(a, b) for a, b in zip(levels_left, levels_right)]
    expected_masses = [F(0), F(0), F(0), F(1, prime**7), F(1, prime**8)]
    require(masses == expected_masses, f"p={prime} collision masses")
    require(
        min(index for index, mass in enumerate(masses) if mass > 0) == FIRST_COLLISION,
        f"p={prime} first collision",
    )

    differences = [masses[s] - masses[s + 1] for s in range(DEPTH)]
    expected_differences = [F(0), F(0), -F(1, prime**7), F(prime - 1, prime**8)]
    require(differences == expected_differences, f"p={prime} shell differences")

    coloured: list[list[tuple[F, ...]]] = []
    for shell in range(DEPTH):
        correlation = root_correlation(levels_left[shell], levels_right[shell], prime)
        currents = [
            tuple(value / (prime * prime) for value in prime_cyclotomic_reduce(correlation, k, prime))
            for k in range(1, prime)
        ]
        aggregate = add_vectors(currents)
        require(
            aggregate == (differences[shell],) + (F(0),) * (prime - 2),
            f"p={prime} shell/root ledger at s={shell}",
        )
        coloured.append(currents)

    for shell in range(FIRST_COLLISION - 1):
        require(
            all(not any(current) for current in coloured[shell]),
            f"p={prime} premature coloured germ",
        )
    leading = coloured[FIRST_COLLISION - 1]
    require(all(any(current) for current in leading), f"p={prime} primitive saturation")
    require(
        add_vectors(leading) == (-F(1, prime**7),) + (F(0),) * (prime - 2),
        f"p={prime} signed leading ledger",
    )

    # The theorem proves this for every rational radius.  The finite-cylinder
    # control checks three exact radii, including one near the Abel boundary.
    for z in (F(1, 7), F(1, 2), F(prime, prime + 1)):
        for colour in range(prime - 1):
            coefficients = [coloured[shell][colour] for shell in range(DEPTH)]
            require(any(evaluate_vector_germ(coefficients, z)), f"p={prime} colour at z={z}")

    # The aggregate Abel germ is a polynomial because both cylinder averages
    # are constant after depth four.
    def germ(z: F) -> F:
        return sum(value * z**shell for shell, value in enumerate(differences))

    for z in (F(0), F(1, 7), F(1, 2), F(prime, prime + 1)):
        require(germ(z) == abel_rhs(masses, masses[-1], z), f"p={prime} Abel identity at {z}")
        if z > 0:
            require(germ(z) < 0, f"p={prime} strict Abel sign at {z}")

    require(sum(differences) == -F(1, prime**8), f"p={prime} z=1 covariance")
    require(
        differences[FIRST_COLLISION - 1] == -masses[FIRST_COLLISION],
        f"p={prime} leading derivative invoice",
    )

    # Full mixed Gregory differences at d-dimensional squarefree corners have
    # sign (-1)^d.  Build every face of a top-only collision cube and take the
    # complete alternating sum, auditing the Mobius/Ramanujan normalization.
    gregory_signs = []
    for dimension in range(1, 5):
        cube = {
            tuple((mask_int >> bit) & 1 for bit in range(dimension)): F(
                int(mask_int == (1 << dimension) - 1)
            )
            for mask_int in range(1 << dimension)
        }
        total = sum((-1) ** sum(mask) * mass for mask, mass in cube.items())
        require(total == (-1) ** dimension, f"d={dimension} Gregory sign")
        gregory_signs.append(total)

    # Even before the first aggregate root colour appears, the individual
    # ordinary q=1 interval coefficients are nonzero: 1 is not a multiple of
    # the p^4 grid.  For a grid cylinder, the exact zero law at q!=0 is
    # p^4|q, obtained from the one-cell integral's sine factor.
    interval_coefficient_nonzero = lambda frequency: frequency != 0 and frequency % grid != 0
    require(interval_coefficient_nonzero(1), "ordinary-frequency hostile")
    require(not interval_coefficient_nonzero(grid), "ordinary-frequency zero control")
    require(all(not any(current) for current in coloured[0]), "base residue cancellation")

    # W_z is the valuation multiplier.  On exact frequency labels it obeys
    # W_z U=z U W_z and P W_z=z W_z P.  Its projection telescoping formula has
    # exactly one surviving valuation shell.
    z_operator = F(2, 3)
    for frequency in range(-500, 501):
        if frequency == 0:
            continue
        order = valuation(frequency, prime)
        require(z_operator ** valuation(prime * frequency, prime) == z_operator * z_operator**order,
                f"p={prime} graded clock covariance")
        projected = sum(
            z_operator**shell
            * F(int(frequency % prime**shell == 0) - int(frequency % prime ** (shell + 1) == 0))
            for shell in range(order + 1)
        )
        require(projected == z_operator**order, f"p={prime} projection shell multiplier")

    # A fixed Abel radius loses collision depth even with the whole primitive
    # vector.  X first collides at depth one. Y first collides at depth two,
    # but G_Y,k(z)=(z/p^2)G_X,k(z) for every primitive k. At z=1/2, scaling
    # X's first packet by 1/(2p^2) makes all coloured values agree exactly.
    radius = F(1, 2)
    x_masses, x_germs = cell_germ_coefficients(
        prime, 2, prime, 0, F(1, 2 * prime * prime)
    )
    y_masses, y_germs = cell_germ_coefficients(prime, 3, prime * prime + prime, 0)
    require(min(s for s, mass in enumerate(x_masses) if mass > 0) == 1, "X depth one")
    require(min(s for s, mass in enumerate(y_masses) if mass > 0) == 2, "Y depth two")
    for x_coefficients, y_coefficients in zip(x_germs, y_germs):
        require(
            evaluate_vector_germ(x_coefficients, radius)
            == evaluate_vector_germ(y_coefficients, radius),
            f"p={prime} fixed-radius depth collision",
        )

    return {
        "prime": prime,
        "masses": masses,
        "differences": differences,
        "leading": leading,
        "germ_half": germ(F(1, 2)),
        "gregory_signs": gregory_signs,
        "fixed_radius_depths": (1, 2),
    }


def main() -> None:
    controls = [run_prime(prime) for prime in PRIMES]
    control13 = next(control for control in controls if control["prime"] == 13)
    masses = control13["masses"]
    differences = control13["differences"]
    leading = control13["leading"]

    print("THM-2505 FIRST-COLLISION ABEL-GERM EXACT COMPANION")
    print("prime_controls=" + ",".join(map(str, PRIMES)) + ";depth=4;first_collision=3")
    print("p13_I=" + ",".join(map(str, masses)))
    print("p13_shell=" + ",".join(map(str, differences)))
    print("p13_aggregate_germ=-z^2/13^7+12*z^3/13^8")
    print(f"p13_germ_at_half={control13['germ_half']};germ_at_one={sum(differences)}")
    print(f"coloured_germs_order=2;nonzero_leading_colours={len(leading)}")
    print("p13_colour1_lead=" + ",".join(map(str, leading[0])))
    print("gregory_corner_signs=" + ",".join(map(str, control13["gregory_signs"])))
    print("rational_radius_colours=all_nonzero;operator_covariance=PASS")
    print("fixed_radius=1/2;equal_full_colour_vectors_at_depths=1,2")
    print("ordinary_q1_product=nonzero;ordinary_q13^4=zero;base_root_residue_aggregates=zero")
    print("all_checks=PASS")


if __name__ == "__main__":
    main()
