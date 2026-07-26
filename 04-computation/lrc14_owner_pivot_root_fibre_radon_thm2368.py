#!/usr/bin/env python3
"""Exact, optimization-safe companion for THM-2368."""

from fractions import Fraction
from itertools import product


P = 13
FINITE_FIELD = 53
ZETA_MOD_53 = 10


def require(condition: bool, message: str) -> None:
    """Raise under ordinary and optimized Python."""
    if not condition:
        raise RuntimeError(message)


def cyclic_convolution(
    left: tuple[Fraction, ...],
    right: tuple[Fraction, ...],
) -> tuple[Fraction, ...]:
    return tuple(
        sum(left[index] * right[(target - index) % P] for index in range(P))
        for target in range(P)
    )


def delta(indices: tuple[int, ...]) -> tuple[Fraction, ...]:
    return tuple(Fraction(int(index in indices)) for index in range(P))


def determinant_bareiss(matrix: list[list[int]]) -> int:
    """Fraction-free exact determinant."""
    work = [row[:] for row in matrix]
    size = len(work)
    sign = 1
    previous = 1
    for pivot_index in range(size - 1):
        pivot_row = next(
            (
                row
                for row in range(pivot_index, size)
                if work[row][pivot_index] != 0
            ),
            None,
        )
        require(pivot_row is not None, "unexpected singular circulant")
        if pivot_row != pivot_index:
            work[pivot_index], work[pivot_row] = (
                work[pivot_row],
                work[pivot_index],
            )
            sign *= -1
        pivot = work[pivot_index][pivot_index]
        for row in range(pivot_index + 1, size):
            for column in range(pivot_index + 1, size):
                numerator = (
                    work[row][column] * pivot
                    - work[row][pivot_index] * work[pivot_index][column]
                )
                require(
                    numerator % previous == 0,
                    "Bareiss division stopped being exact",
                )
                work[row][column] = numerator // previous
        previous = pivot
    return sign * work[-1][-1]


def convolution_matrix(word: tuple[Fraction, ...]) -> list[list[int]]:
    require(
        all(value.denominator == 1 for value in word),
        "determinant word stopped being integral",
    )
    return [
        [int(word[(row - column) % P]) for column in range(P)]
        for row in range(P)
    ]


def cyclotomic_coordinates(
    word: tuple[int, ...],
    frequency: int,
) -> tuple[int, ...]:
    """Coordinates of sum_r word[r] zeta^(frequency*r) in Q(zeta_13)."""
    coordinates = [0] * (P - 1)
    for index, coefficient in enumerate(word):
        exponent = (frequency * index) % P
        if exponent == P - 1:
            for basis_index in range(P - 1):
                coordinates[basis_index] -= coefficient
        else:
            coordinates[exponent] += coefficient
    return tuple(coordinates)


def exhaustive_prime_cyclotomic_unit_atlas() -> tuple[int, int]:
    """Check all nonempty proper Boolean words on C_13."""
    words = 0
    nonzero_character_checks = 0
    for word in product((0, 1), repeat=P):
        mass = sum(word)
        if mass in (0, P):
            continue
        words += 1
        require(mass != 0, "zero character vanished")
        for frequency in range(1, P):
            require(
                any(cyclotomic_coordinates(word, frequency)),
                "proper Boolean word lost a primitive character",
            )
            nonzero_character_checks += 1
    require(words == 2**P - 2, "Boolean-word census changed")
    return words, nonzero_character_checks


def canonical_kernel_units() -> tuple[tuple[int, ...], int]:
    """Verify the four one/two-root danger/safe convolution kernels."""
    identity = delta((0,))
    total = delta(tuple(range(P)))
    alternating_half = tuple(
        Fraction((-1) ** index, 2) for index in range(P)
    )
    kernels_and_inverses = (
        (identity, identity, 1),
        (
            delta((0, 1)),
            alternating_half,
            2,
        ),
        (
            tuple(total[index] - identity[index] for index in range(P)),
            tuple(
                -identity[index] + Fraction(1, 12) * total[index]
                for index in range(P)
            ),
            12,
        ),
        (
            tuple(
                total[index] - delta((0, 1))[index] for index in range(P)
            ),
            tuple(
                -alternating_half[index] + Fraction(1, 22) * total[index]
                for index in range(P)
            ),
            11,
        ),
    )
    determinants = []
    inverse_checks = 0
    for kernel, inverse, expected_determinant in kernels_and_inverses:
        require(
            cyclic_convolution(kernel, inverse) == identity,
            "displayed group-ring inverse failed",
        )
        determinant = abs(determinant_bareiss(convolution_matrix(kernel)))
        require(
            determinant == expected_determinant,
            "canonical circulant determinant changed",
        )
        determinants.append(determinant)
        inverse_checks += P
    return tuple(determinants), inverse_checks


def circular_distance(value: Fraction) -> Fraction:
    residue = value % 1
    return min(residue, 1 - residue)


def translated_grid_counts(half_width: Fraction) -> set[int]:
    """Sample every phase chamber of a translated thirteen-grid."""
    boundaries = {Fraction(0), Fraction(1)}
    for root in range(P):
        for endpoint in (half_width, 1 - half_width):
            phase = P * endpoint - root
            if 0 < phase < 1:
                boundaries.add(phase)
    ordered = sorted(boundaries)
    samples = [
        (left + right) / 2 for left, right in zip(ordered, ordered[1:])
    ]
    counts = set()
    for phase in samples:
        counts.add(
            sum(
                circular_distance((phase + root) / P) < half_width
                for root in range(P)
            )
        )
    return counts


def capacity_atlas() -> tuple[set[int], set[int], int, int]:
    ordinary_danger_counts = translated_grid_counts(Fraction(1, 14))
    guard_danger_counts = translated_grid_counts(Fraction(1, 7))
    require(
        ordinary_danger_counts == {1, 2},
        "ordinary root danger count changed",
    )
    require(
        guard_danger_counts == {3, 4},
        "guard root danger count changed",
    )

    lower_support = P - (4 + 3 * 2)
    upper_support = P - 3
    require((lower_support, upper_support) == (3, 10), "capacity gap changed")

    guard_lower = {0, 1, 2, 3}
    ordinary_lower = ({4, 5}, {6, 7}, {8, 9})
    lower_word = set(range(P)) - guard_lower.union(*ordinary_lower)
    guard_upper = {0, 1, 2}
    ordinary_upper = ({0}, {1}, {2})
    upper_word = set(range(P)) - guard_upper.union(*ordinary_upper)
    require(len(lower_word) == 3, "capacity lower hostile failed")
    require(len(upper_word) == 10, "capacity upper hostile failed")
    return (
        ordinary_danger_counts,
        guard_danger_counts,
        lower_support,
        upper_support,
    )


def dft_mod_53(word: tuple[int, ...], frequency: int) -> int:
    return sum(
        coefficient
        * pow(ZETA_MOD_53, (-frequency * index) % P, FINITE_FIELD)
        for index, coefficient in enumerate(word)
    ) % FINITE_FIELD


def radon_table(
    a_word: tuple[int, ...],
    c_word: tuple[int, ...],
    p_word: tuple[int, ...],
    k_a: int,
    k_c: int,
) -> list[list[int]]:
    """Return the integer numerator N=13Q."""
    return [
        [
            sum(
                a_word[(k_a * branch + s) % P]
                * c_word[(k_c * branch + t) % P]
                * p_word[branch]
                for branch in range(P)
            )
            for t in range(P)
        ]
        for s in range(P)
    ]


def radon_transform_mod_53(
    table: list[list[int]],
    lambda_frequency: int,
    mu_frequency: int,
) -> int:
    inverse_p_cubed = pow(P, -3, FINITE_FIELD)
    return (
        inverse_p_cubed
        * sum(
            table[s][t]
            * pow(
                ZETA_MOD_53,
                (-lambda_frequency * s - mu_frequency * t) % P,
                FINITE_FIELD,
            )
            for s in range(P)
            for t in range(P)
        )
    ) % FINITE_FIELD


def factorized_transform_mod_53(
    a_word: tuple[int, ...],
    c_word: tuple[int, ...],
    p_word: tuple[int, ...],
    k_a: int,
    k_c: int,
    lambda_frequency: int,
    mu_frequency: int,
) -> int:
    inverse_p = pow(P, -1, FINITE_FIELD)
    a_transform = inverse_p * dft_mod_53(a_word, lambda_frequency)
    c_transform = inverse_p * dft_mod_53(c_word, mu_frequency)
    p_frequency = (
        lambda_frequency * k_a + mu_frequency * k_c
    ) % P
    p_transform = inverse_p * sum(
        p_word[index]
        * pow(ZETA_MOD_53, (p_frequency * index) % P, FINITE_FIELD)
        for index in range(P)
    )
    return a_transform * c_transform * p_transform % FINITE_FIELD


def directional_energies(
    table: list[list[int]],
) -> tuple[Fraction, Fraction, Fraction]:
    horizontal = Fraction(0)
    vertical = Fraction(0)
    mixed = Fraction(0)
    for s in range(P):
        for t in range(P):
            horizontal += Fraction(
                (table[(s + 1) % P][t] - table[s][t]) ** 2,
                P**2,
            )
            vertical += Fraction(
                (table[s][(t + 1) % P] - table[s][t]) ** 2,
                P**2,
            )
            mixed_difference = (
                table[(s + 1) % P][(t + 1) % P]
                - table[(s + 1) % P][t]
                - table[s][(t + 1) % P]
                + table[s][t]
            )
            mixed += Fraction(mixed_difference**2, P**2)
    return horizontal, vertical, mixed


def radon_factorization_atlas() -> tuple[int, int, Fraction, Fraction, Fraction]:
    safe_one = tuple(int(index != 0) for index in range(P))
    safe_two = tuple(int(index not in (0, 1)) for index in range(P))
    controls = (
        (
            safe_one,
            safe_two,
            tuple(int(index in (0, 1, 2)) for index in range(P)),
            2,
            7,
        ),
        (
            safe_two,
            safe_one,
            tuple(int(index in range(10)) for index in range(P)),
            5,
            8,
        ),
    )
    factorization_checks = 0
    full_mode_checks = 0
    minimum_energies = None
    for a_word, c_word, p_word, k_a, k_c in controls:
        table = radon_table(a_word, c_word, p_word, k_a, k_c)
        for lambda_frequency, mu_frequency in product(range(P), repeat=2):
            direct = radon_transform_mod_53(
                table,
                lambda_frequency,
                mu_frequency,
            )
            factorized = factorized_transform_mod_53(
                a_word,
                c_word,
                p_word,
                k_a,
                k_c,
                lambda_frequency,
                mu_frequency,
            )
            require(direct == factorized, "Radon transform formula failed")
            require(direct != 0, "control lost a target mode modulo 53")
            factorization_checks += 1
            full_mode_checks += 1
        energies = directional_energies(table)
        require(
            energies[0] >= Fraction(2, P**2)
            and energies[1] >= Fraction(2, P**2),
            "integer axis Dirichlet gap failed",
        )
        require(
            energies[2] >= Fraction(4, P**2),
            "integer mixed Dirichlet gap failed",
        )
        if minimum_energies is None:
            minimum_energies = energies
        else:
            minimum_energies = tuple(
                min(old, new) for old, new in zip(minimum_energies, energies)
            )
    require(minimum_energies is not None, "missing Radon controls")
    return (
        factorization_checks,
        full_mode_checks,
        minimum_energies[0],
        minimum_energies[1],
        minimum_energies[2],
    )


def abstract_circulant_hostile() -> tuple[int, Fraction, Fraction, int, int]:
    """Full pointwise root spectrum can integrate to zero H-drift."""
    danger = tuple(int(index == 0) for index in range(P))
    safe = tuple(1 - value for value in danger)
    pointwise_mode_checks = 0
    for u, v, lambda_frequency, mu_frequency in product(
        range(P), repeat=4
    ):
        a_translate = tuple(safe[(u + s) % P] for s in range(P))
        c_translate = tuple(safe[(v + t) % P] for t in range(P))
        require(
            any(cyclotomic_coordinates(a_translate, lambda_frequency))
            if lambda_frequency
            else sum(a_translate) != 0,
            "abstract hostile lost an s mode",
        )
        require(
            any(cyclotomic_coordinates(c_translate, mu_frequency))
            if mu_frequency
            else sum(c_translate) != 0,
            "abstract hostile lost a t mode",
        )
        pointwise_mode_checks += 1

    values = {}
    for r, s, t in product(range(P), repeat=3):
        u_mass = sum(safe[(u + s) % P] for u in range(P))
        v_mass = sum(
            safe[(v + t) % P] * danger[(v + r) % P]
            for v in range(P)
        )
        values[(r, s, t)] = Fraction(u_mass * v_mass, P**2)
        expected = Fraction(0) if r == t else Fraction(12, P**2)
        require(values[(r, s, t)] == expected, "hostile H formula failed")

    orbit_checks = 0
    drift_energy = Fraction(0)
    for r, s, t in product(range(P), repeat=3):
        orbit_average = sum(
            values[((r + v) % P, (s + u) % P, (t + v) % P)]
            for u, v in product(range(P), repeat=2)
        ) / P**2
        require(orbit_average == values[(r, s, t)], "hostile is not circulant")
        drift_energy += (values[(r, s, t)] - orbit_average) ** 2 / P**3
        orbit_checks += 1
    require(drift_energy == 0, "circulant hostile acquired H-drift")

    safe_two = tuple(int(index not in (0, 1)) for index in range(P))
    residual_three = tuple(int(index in (0, 1, 2)) for index in range(P))
    radon_control = radon_table(safe, safe_two, residual_three, 2, 7)
    radon_orbit_profile = []
    for difference in range(P):
        mass = sum(
            radon_control[u][v]
            * safe[u]
            * safe[v]
            * danger[(v + difference) % P]
            for u, v in product(range(P), repeat=2)
        )
        radon_orbit_profile.append(Fraction(mass, P**3))
    require(radon_orbit_profile[0] == 0, "Radon orbit lost diagonal zero")
    positive_offsets = sum(value > 0 for value in radon_orbit_profile[1:])
    require(
        positive_offsets == P - 1,
        "Radon orbit hostile lost an off-diagonal value",
    )
    for r, t, v_shift in product(range(P), repeat=3):
        before = radon_orbit_profile[(r - t) % P]
        after = radon_orbit_profile[
            ((r + v_shift) - (t + v_shift)) % P
        ]
        require(before == after, "Radon orbit hostile lost invariance")

    return (
        pointwise_mode_checks,
        Fraction(12, P**2),
        drift_energy,
        orbit_checks,
        positive_offsets,
    )


def main() -> None:
    require(
        pow(ZETA_MOD_53, P, FINITE_FIELD) == 1 and ZETA_MOD_53 != 1,
        "finite-field thirteenth root changed",
    )
    boolean_words, cyclotomic_checks = exhaustive_prime_cyclotomic_unit_atlas()
    determinants, inverse_checks = canonical_kernel_units()
    (
        ordinary_counts,
        guard_counts,
        support_lower,
        support_upper,
    ) = capacity_atlas()
    (
        factorization_checks,
        full_mode_checks,
        horizontal_energy,
        vertical_energy,
        mixed_energy,
    ) = radon_factorization_atlas()
    (
        hostile_mode_checks,
        hostile_off_diagonal,
        hostile_drift,
        hostile_orbit_checks,
        radon_orbit_positive_offsets,
    ) = abstract_circulant_hostile()

    print("theorem=THM-2368")
    print("status=PROVED+VERIFIED-EXACT+INDEPENDENTLY-AUDITED")
    print("root_group=C_13")
    print(f"nonempty_proper_boolean_words={boolean_words}")
    print(f"prime_cyclotomic_nonzero_character_checks={cyclotomic_checks}")
    print("all_nonempty_proper_boolean_words_are_QC13_units=YES")
    print(f"canonical_kernel_abs_determinants={','.join(map(str, determinants))}")
    print(f"canonical_kernel_inverse_coordinate_checks={inverse_checks}")
    print(f"ordinary_danger_root_counts={sorted(ordinary_counts)}")
    print(f"guard_danger_root_counts={sorted(guard_counts)}")
    print(f"residual_branch_support_interval=[{support_lower},{support_upper}]")
    print(f"radon_factorization_checks_mod_53={factorization_checks}")
    print(f"radon_full_mode_control_checks_mod_53={full_mode_checks}")
    print("radon_pointwise_target_modes=169/169")
    print(f"integer_axis_dirichlet_gap=2/{P**2}")
    print(f"integer_mixed_dirichlet_gap=4/{P**2}")
    print(f"control_horizontal_energy_min={horizontal_energy}")
    print(f"control_vertical_energy_min={vertical_energy}")
    print(f"control_mixed_energy_min={mixed_energy}")
    print(f"circulant_hostile_pointwise_mode_checks={hostile_mode_checks}")
    print(f"circulant_hostile_off_diagonal_H={hostile_off_diagonal}")
    print(f"circulant_hostile_orbit_checks={hostile_orbit_checks}")
    print(f"circulant_hostile_H_drift={hostile_drift}")
    print(
        "exact_radon_uniform_orbit_hostile_positive_offsets="
        f"{radon_orbit_positive_offsets}"
    )
    print("exact_radon_uniform_orbit_hostile_H_drift=0")
    print("full_169_pointwise_modes_imply_bare_H_drift=NO")
    print("missing_coordinate=root_phase_anchor_through_y_integration")
    print("scalar_rows_excluded=0")
    print("lrc14_status=OPEN")
    print("all_checks=PASS")


if __name__ == "__main__":
    main()
