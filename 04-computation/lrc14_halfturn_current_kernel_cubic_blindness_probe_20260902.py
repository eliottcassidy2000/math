#!/usr/bin/env python3
"""Exact audits for the half-turn current kernel and cubic-moment blindness.

For an odd speed w, put

    sigma_w(t) = 1_{||wt||<1/14} - 1_{||w(t+1/2)||<1/14}.

The first part of this script compares an exact wall sweep for
integral sigma_u sigma_v with its closed mod-14 formula.  It then audits the
global Dirichlet-convolution Fourier coefficient formula, the exact 7-adic
shell/reverse-martingale decomposition, and a geometric-chain hostile to a
quadratic-current closure.  The final part checks two exchangeable labelled
occupation laws which agree on every factorial moment through total degree
three but have different zero-minimum mass.  The latter is an abstract moment
obstruction, not a realizable LRC counterexample.
"""

from __future__ import annotations

from fractions import Fraction
from math import ceil, comb, floor, gcd


def distance_mod_14(value: int) -> int:
    residue = value % 14
    return min(residue, 14 - residue)


def predicted_correlation(u: int, v: int) -> Fraction:
    common = gcd(u, v)
    reduced_u = u // common
    reduced_v = v // common
    assert reduced_u % 2 == reduced_v % 2 == 1
    return Fraction(
        distance_mod_14(reduced_u + reduced_v)
        - distance_mod_14(reduced_u - reduced_v),
        7 * reduced_u * reduced_v,
    )


def danger_intervals(w: int, left: Fraction, right: Fraction) -> list[tuple[Fraction, Fraction]]:
    answer: list[tuple[Fraction, Fraction]] = []
    low_index = floor(w * left - Fraction(1, 14)) - 1
    high_index = ceil(w * right + Fraction(1, 14)) + 1
    for index in range(low_index, high_index + 1):
        tooth_left = Fraction(14 * index - 1, 14 * w)
        tooth_right = Fraction(14 * index + 1, 14 * w)
        clipped_left = max(left, tooth_left)
        clipped_right = min(right, tooth_right)
        if clipped_left < clipped_right:
            answer.append((clipped_left, clipped_right))
    return answer


def direct_correlation(u: int, v: int) -> Fraction:
    """Integrate on [0,1/2) with measure 2dt by an exact wall sweep."""
    left = Fraction()
    right = Fraction(1, 2)
    shift = Fraction(1, 2)
    data: list[tuple[list[tuple[Fraction, Fraction]], list[tuple[Fraction, Fraction]]]] = []
    walls = {left, right}

    for speed in (u, v):
        low = danger_intervals(speed, left, right)
        high = [
            (interval_left - shift, interval_right - shift)
            for interval_left, interval_right in danger_intervals(
                speed, left + shift, right + shift
            )
        ]
        data.append((low, high))
        for interval_left, interval_right in low + high:
            walls.add(interval_left)
            walls.add(interval_right)

    correlation = Fraction()
    ordered = sorted(walls)
    for cell_left, cell_right in zip(ordered, ordered[1:]):
        midpoint = (cell_left + cell_right) / 2
        values = []
        for low, high in data:
            low_value = sum(a < midpoint < b for a, b in low)
            high_value = sum(a < midpoint < b for a, b in high)
            values.append(low_value - high_value)
        correlation += 2 * (cell_right - cell_left) * values[0] * values[1]
    return correlation


def rational_rank(matrix: list[list[Fraction]]) -> int:
    work = [row[:] for row in matrix]
    rows = len(work)
    columns = len(work[0])
    pivot_row = 0
    for column in range(columns):
        pivot = next(
            (row for row in range(pivot_row, rows) if work[row][column]),
            None,
        )
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        scale = work[pivot_row][column]
        work[pivot_row] = [entry / scale for entry in work[pivot_row]]
        for row in range(rows):
            if row == pivot_row or not work[row][column]:
                continue
            multiple = work[row][column]
            work[row] = [
                entry - multiple * pivot_entry
                for entry, pivot_entry in zip(work[row], work[pivot_row])
            ]
        pivot_row += 1
        if pivot_row == rows:
            break
    return pivot_row


def brownian_coordinate(residue: int) -> tuple[int, int]:
    """Return (orientation, time) for the three-step Brownian factorization."""
    reduced = residue % 14
    coordinates = {
        1: (1, 1),
        5: (1, 2),
        3: (1, 3),
        13: (-1, 1),
        9: (-1, 2),
        11: (-1, 3),
        7: (0, 0),
    }
    return coordinates[reduced]


SineVector = tuple[Fraction, Fraction, Fraction]
ZERO_SINE_VECTOR: SineVector = (Fraction(), Fraction(), Fraction())


def valuation(number: int, prime: int) -> int:
    exponent = 0
    while number % prime == 0:
        number //= prime
        exponent += 1
    return exponent


def add_vectors(left: SineVector, right: SineVector) -> SineVector:
    return tuple(a + b for a, b in zip(left, right))  # type: ignore[return-value]


def scale_vector(scale: Fraction, vector: SineVector) -> SineVector:
    return tuple(scale * entry for entry in vector)  # type: ignore[return-value]


def sine_vector(k: int) -> SineVector:
    """Represent sin(pi*k/7) in the basis sin(pi/7),sin(2pi/7),sin(3pi/7)."""
    assert k > 0 and k % 2 == 1
    table: dict[int, SineVector] = {
        1: (Fraction(1), Fraction(), Fraction()),
        3: (Fraction(), Fraction(), Fraction(1)),
        5: (Fraction(), Fraction(1), Fraction()),
        7: ZERO_SINE_VECTOR,
        9: (Fraction(), Fraction(-1), Fraction()),
        11: (Fraction(), Fraction(), Fraction(-1)),
        13: (Fraction(-1), Fraction(), Fraction()),
    }
    return table[k % 14]


def direct_normalized_fourier_coefficient(speeds: tuple[int, ...], n: int) -> SineVector:
    """Return (pi/2)*C_hat(n) as an exact formal sine-basis vector."""
    answer = ZERO_SINE_VECTOR
    for speed in speeds:
        if n % speed:
            continue
        quotient = n // speed
        if quotient % 2 == 0:
            continue
        answer = add_vectors(
            answer,
            scale_vector(Fraction(1, quotient), sine_vector(quotient)),
        )
    return answer


def convolution_normalized_fourier_coefficient(
    speeds: tuple[int, ...], n: int
) -> SineVector:
    """The equivalent divisor convolution: n^-1 sum_(w|n) w sin(pi*n/(7w))."""
    answer = ZERO_SINE_VECTOR
    speed_set = set(speeds)
    for divisor in range(1, n + 1):
        if n % divisor or divisor not in speed_set:
            continue
        quotient = n // divisor
        if quotient % 2 == 0:
            continue
        answer = add_vectors(
            answer,
            scale_vector(Fraction(divisor, n), sine_vector(quotient)),
        )
    return answer


def current_energy(speeds: tuple[int, ...]) -> Fraction:
    return sum(
        (predicted_correlation(left, right) for left in speeds for right in speeds),
        start=Fraction(),
    )


def shell_minimal_speeds(speeds: tuple[int, ...]) -> tuple[int, ...]:
    """Divisibility-minimal speeds separately inside each exact 7-adic shell."""
    return tuple(
        speed
        for speed in speeds
        if not any(
            other != speed
            and valuation(other, 7) == valuation(speed, 7)
            and speed % other == 0
            for other in speeds
        )
    )


def audit_dirichlet_square_and_filtration() -> None:
    # Every subset of this bounded odd universe is checked.  Formal sine
    # vectors make the coefficient comparison exact rather than floating.
    universe = tuple(range(1, 16, 2))
    frequency_cap = 315
    coefficient_rows = 0
    projector_rows = 0
    shell_energy_rows = 0
    private_frequency_rows = 0

    for mask in range(1 << len(universe)):
        speeds = tuple(
            speed for index, speed in enumerate(universe) if mask & (1 << index)
        )
        shells: dict[int, tuple[int, ...]] = {}
        for exponent in sorted({valuation(speed, 7) for speed in speeds}):
            shells[exponent] = tuple(
                speed for speed in speeds if valuation(speed, 7) == exponent
            )

        # Pairwise kernel orthogonality agrees exactly with the shell split.
        assert current_energy(speeds) == sum(
            (current_energy(shell) for shell in shells.values()),
            start=Fraction(),
        )
        shell_energy_rows += 1

        for n in range(1, frequency_cap + 1, 2):
            direct = direct_normalized_fourier_coefficient(speeds, n)
            convolution = convolution_normalized_fourier_coefficient(speeds, n)
            assert direct == convolution
            coefficient_rows += 1

            # Q_e keeps precisely the modes divisible by 7^e.  Therefore
            # Q_e-Q_(e+1) selects the exact valuation shell e.  At e=0,
            # Q_0 is the identity.  The current has zero constant mode.
            for exponent in range(4):
                projected = (
                    direct
                    if n % 7**exponent == 0 and n % 7 ** (exponent + 1) != 0
                    else ZERO_SINE_VECTOR
                )
                shell_coefficient = direct_normalized_fourier_coefficient(
                    shells.get(exponent, ()), n
                )
                assert projected == shell_coefficient
                projector_rows += 1

        # At n=w, every shell-minimal w has its uncancelled fundamental
        # coefficient 2*sin(pi/7)/pi.  Different w give different modes.
        for speed in shell_minimal_speeds(speeds):
            assert direct_normalized_fourier_coefficient(speeds, speed) == (
                Fraction(1),
                Fraction(),
                Fraction(),
            )
            private_frequency_rows += 1

    # A one-shell geometric chain is an exact hostile to any attempted
    # universal EC^2>4 closure.  Its reduced pair residues have period three:
    # 9^d mod 14 = 9,11,1 and signs -,-,+.
    chain = tuple(9**exponent for exponent in range(12))
    chain_energy = current_energy(chain)
    assert chain_energy == Fraction(585613119200, 219667417263)
    assert chain_energy < 4
    assert {valuation(speed, 7) for speed in chain} == {0}
    for left_index, left in enumerate(chain):
        for right_index, right in enumerate(chain):
            distance = abs(left_index - right_index)
            if distance == 0:
                expected = Fraction(2, 7)
            else:
                sign = 1 if distance % 3 == 0 else -1
                expected = Fraction(2 * sign, 7 * 9**distance)
            assert predicted_correlation(left, right) == expected

    signed_geometric_tail = (
        -Fraction(1, 9) - Fraction(1, 81) + Fraction(1, 729)
    ) / (1 - Fraction(1, 729))
    asymptotic_energy_per_speed = Fraction(2, 7) + Fraction(4, 7) * signed_geometric_tail
    assert signed_geometric_tail == Fraction(-89, 728)
    assert asymptotic_energy_per_speed == Fraction(275, 1274)

    print("dirichlet_square_and_septimal_filtration")
    print(
        " exact_universe=all_256_subsets_of_odd_1_to_15"
        f" odd_frequency_cap={frequency_cap}"
    )
    print(f" exact_formal_coefficient_rows={coefficient_rows}")
    print(
        " coefficient=(1_W*_D beta)(n),"
        " beta(k)=2sin(pi*k/7)/(pi*k)_for_odd_k"
    )
    print(
        " energy=2*sum_n>=1_odd|(1_W*_D_beta)(n)|^2"
        "=8/pi^2*sum_n>=1_odd n^-2*(sum_w|n w*sin(pi*n/(7w)))^2"
    )
    print(
        f" exact_shell_energy_rows={shell_energy_rows}"
        f" exact_projector_rows={projector_rows}"
    )
    print(
        " Q_e=fourier_projection_to_7^e_dividing_n;"
        " (Q_e-Q_(e+1))C=C_{nu7(speed)=e}"
    )
    print(f" exact_private_fundamental_rows={private_frequency_rows}")
    print(
        " private_floor_per_shell_minimal_speed="
        "8*sin(pi/7)^2/pi^2"
    )
    print(" chain_9_speeds=1,9,...,9^11")
    print(f" chain_9_energy={chain_energy}<4")
    print(
        " chain_9_pair_sign_period=-,-,+"
        f" asymptotic_energy_per_speed={asymptotic_energy_per_speed}"
    )


def audit_current_kernel() -> None:
    checked = 0
    for u in range(1, 100, 2):
        for v in range(1, 100, 2):
            direct = direct_correlation(u, v)
            predicted = predicted_correlation(u, v)
            assert direct == predicted, (u, v, direct, predicted)
            checked += 1

    residues = (1, 3, 5, 7, 9, 11, 13)
    numerator_matrix = [
        [
            Fraction(
                distance_mod_14(row + column)
                - distance_mod_14(row - column)
            )
            for column in residues
        ]
        for row in residues
    ]
    assert rational_rank(numerator_matrix) == 3
    positive_block = [row[:3] for row in numerator_matrix[:3]]
    determinant = (
        positive_block[0][0]
        * (
            positive_block[1][1] * positive_block[2][2]
            - positive_block[1][2] * positive_block[2][1]
        )
        - positive_block[0][1]
        * (
            positive_block[1][0] * positive_block[2][2]
            - positive_block[1][2] * positive_block[2][0]
        )
        + positive_block[0][2]
        * (
            positive_block[1][0] * positive_block[2][1]
            - positive_block[1][1] * positive_block[2][0]
        )
    )
    assert determinant == 8

    # In the non-geometric order 1,5,3 the positive block is twice the
    # Brownian covariance min(i,j).  Opposite residues negate the feature
    # vector, while residue 7 is null.
    for i, row_residue in enumerate(residues):
        row_sign, row_time = brownian_coordinate(row_residue)
        for j, column_residue in enumerate(residues):
            column_sign, column_time = brownian_coordinate(column_residue)
            brownian_value = (
                2 * row_sign * column_sign * min(row_time, column_time)
            )
            assert numerator_matrix[i][j] == brownian_value

    # The sign partition is the additive half-circle partition of odd
    # residues.  Residue 7 is null; opposite representatives negate rows.
    positive = {1, 3, 5}
    negative = {9, 11, 13}
    for i, row_residue in enumerate(residues):
        for j, column_residue in enumerate(residues):
            entry = numerator_matrix[i][j]
            if 7 in (row_residue, column_residue):
                assert entry == 0
            elif (row_residue in positive) == (column_residue in positive):
                assert entry > 0
            else:
                assert entry < 0

    print("halfturn_current_kernel")
    print(f" exact_wall_pairs_checked={checked}")
    print(" formula=(d14(U+V)-d14(U-V))/(7UV), U=u/gcd, V=v/gcd")
    print(" residues=" + ",".join(map(str, residues)))
    print(
        " numerator_matrix="
        + ";".join(",".join(str(int(entry)) for entry in row) for row in numerator_matrix)
    )
    print(" numerator_rank=3 positive_block_determinant=8")
    print(" sign=positive_same_half_negative_opposite_half_zero_at_residue7")
    print(" brownian_order=1,5,3 numerator=2*sign(r)*sign(s)*min(time(r),time(s))")
    print(" brownian_inverse_path_precision=1,-1/2,0;-1/2,1,-1/2;0,-1/2,1/2")


Distribution = dict[tuple[int, int], Fraction]


def factorial_moment(distribution: Distribution, left_order: int, right_order: int) -> Fraction:
    return sum(
        probability * comb(left_depth, left_order) * comb(right_depth, right_order)
        for (left_depth, right_depth), probability in distribution.items()
    )


def zero_minimum_mass(distribution: Distribution) -> Fraction:
    return sum(
        probability
        for (left_depth, right_depth), probability in distribution.items()
        if min(left_depth, right_depth) == 0
    )


def free_sheet_mass(distribution: Distribution) -> Fraction:
    return sum(
        probability * ((left_depth == 0) + (right_depth == 0))
        for (left_depth, right_depth), probability in distribution.items()
    )


def bonferroni_three(depth: int) -> int:
    return 1 - depth + comb(depth, 2) - comb(depth, 3)


def audit_cubic_blindness() -> None:
    # Law A has a positive zero-minimum atom.  Law B has both sheets occupied
    # everywhere.  Conditional on a depth pair, labels are chosen uniformly
    # as two disjoint subsets of a twelve-label set.
    law_a: Distribution = {
        (0, 4): Fraction(5, 56),
        (4, 0): Fraction(5, 56),
        (2, 2): Fraction(15, 28),
        (1, 1): Fraction(2, 7),
    }
    law_b: Distribution = {
        (1, 3): Fraction(5, 14),
        (3, 1): Fraction(5, 14),
        (1, 1): Fraction(2, 7),
    }

    assert sum(law_a.values(), start=Fraction()) == 1
    assert sum(law_b.values(), start=Fraction()) == 1
    assert all(a + b <= 12 for a, b in law_a | law_b)
    assert factorial_moment(law_a, 1, 0) == Fraction(12, 7)
    assert factorial_moment(law_b, 1, 0) == Fraction(12, 7)
    assert factorial_moment(law_a, 0, 1) == Fraction(12, 7)
    assert factorial_moment(law_b, 0, 1) == Fraction(12, 7)
    # Uniform labelling therefore gives every label marginal 1/7 on each
    # sheet, and a label never occupies both sheets.
    assert factorial_moment(law_a, 1, 0) / 12 == Fraction(1, 7)

    cubic_rows: list[str] = []
    for total_degree in range(4):
        for left_order in range(total_degree + 1):
            right_order = total_degree - left_order
            moment_a = factorial_moment(law_a, left_order, right_order)
            moment_b = factorial_moment(law_b, left_order, right_order)
            assert moment_a == moment_b
            cubic_rows.append(f"{left_order},{right_order}:{moment_a}")

    assert zero_minimum_mass(law_a) == Fraction(5, 28)
    assert free_sheet_mass(law_a) == Fraction(5, 28)
    assert zero_minimum_mass(law_b) == 0
    assert free_sheet_mass(law_b) == 0

    quartic_differences = []
    for left_order in range(5):
        right_order = 4 - left_order
        difference = factorial_moment(law_a, left_order, right_order) - factorial_moment(
            law_b, left_order, right_order
        )
        quartic_differences.append(difference)
    assert quartic_differences == (
        [Fraction(5, 56), Fraction(-5, 14), Fraction(15, 28), Fraction(-5, 14), Fraction(5, 56)]
    )

    # Along a+b=4, the signed difference is (5/56) times the fourth
    # finite-difference stencil.  It annihilates every polynomial of degree
    # at most three and fires on degree four.
    stencil = (1, -4, 6, -4, 1)
    for degree in range(4):
        assert sum(weight * point**degree for point, weight in enumerate(stencil)) == 0
    assert sum(weight * point**4 for point, weight in enumerate(stencil)) == 24

    # The nonlinear min-cubic is the maximum of two ordinary third
    # Bonferroni polynomials.  Its extra signal is an absolute-current term,
    # with a useful variance lower bound on the twelve-label state simplex.
    values = tuple(bonferroni_three(depth) for depth in range(13))
    assert all(values[depth + 1] <= values[depth] for depth in range(12))
    for left_depth in range(13):
        for right_depth in range(13 - left_depth):
            total = left_depth + right_depth
            current = left_depth - right_depth
            left_value = bonferroni_three(left_depth)
            right_value = bonferroni_three(right_depth)
            assert bonferroni_three(min(left_depth, right_depth)) == max(
                left_value, right_value
            )
            exact_rebate = Fraction(
                abs(current)
                * (current**2 + 3 * (total - 4) ** 2 - 4),
                24,
            )
            assert abs(left_value - right_value) == exact_rebate
            assert exact_rebate >= Fraction(max(current**2 - 4, 0), 6)

    law_a_q_cubic = sum(
        probability * bonferroni_three(min(left_depth, right_depth))
        for (left_depth, right_depth), probability in law_a.items()
    )
    law_b_q_cubic = sum(
        probability * bonferroni_three(min(left_depth, right_depth))
        for (left_depth, right_depth), probability in law_b.items()
    )
    law_a_one_sheet = sum(
        probability * bonferroni_three(left_depth)
        for (left_depth, _right_depth), probability in law_a.items()
    )
    law_b_one_sheet = sum(
        probability * bonferroni_three(left_depth)
        for (left_depth, _right_depth), probability in law_b.items()
    )
    assert law_a_one_sheet == law_b_one_sheet == 0
    assert law_a_q_cubic == Fraction(5, 28)
    assert law_b_q_cubic == 0

    # Uniform labelling makes every distinct-label signed-current covariance
    # equal.  Its value is incompatible with the exact odd-speed arithmetic
    # kernel:  M/(7UV)=-1/231 would force 33M=-UV, even=odd.
    same_sheet_pair = factorial_moment(law_a, 2, 0) / comb(12, 2)
    mixed_sheet_pair = factorial_moment(law_a, 1, 1) / (12 * 11)
    labelled_current_covariance = 2 * same_sheet_pair - 2 * mixed_sheet_pair
    assert labelled_current_covariance == Fraction(-1, 231)
    assert all(
        predicted_correlation(left, right) != labelled_current_covariance
        for left in range(1, 100, 2)
        for right in range(1, 100, 2)
        if left != right
    )

    print("cubic_moment_blindness")
    print(" labels=12 conditional_labels=uniform_disjoint sheet_swap_symmetric")
    print(" law_A=(0,4):5/56,(4,0):5/56,(2,2):15/28,(1,1):2/7")
    print(" law_B=(1,3):5/14,(3,1):5/14,(1,1):2/7")
    print(" each_label_each_sheet_marginal=1/7 same_label_mixed_intersection=0")
    print(" all_labelled_factorial_tensors_total_degree_le_3=MATCH")
    print(" common_moments=" + ";".join(cubic_rows))
    print(" H0_A=5/28 free_A=5/28 H0_B=0 free_B=0")
    print(" quartic_differences=" + ",".join(map(str, quartic_differences)))
    print(" difference_mechanism=(5/56)*[1,-4,6,-4,1]_along_a_plus_b_eq_4")
    print(" first_mixed_separator=C(a,2)C(b,2), difference=15/28")
    print(" q_cubic=max(bonf3(a),bonf3(b))")
    print(" nonlinear_current=abs(a-b)*((a-b)^2+3*(a+b-4)^2-4)/24")
    print(" current_rebate_lower=abs_bonf3_difference>=max((a-b)^2-4,0)/6")
    print(" law_A_one_sheet_bonf3=0 law_A_q_cubic=5/28")
    print(" law_B_one_sheet_bonf3=0 law_B_q_cubic=0")
    print(" exchangeable_distinct_label_current_covariance=-1/231")
    print(" physical_kernel_exclusion=33*even_numerator_cannot_equal_negative_odd_UV")


def main() -> None:
    audit_current_kernel()
    audit_dirichlet_square_and_filtration()
    audit_cubic_blindness()
    print("SCOPE abstract_moment_obstruction_not_physical_LRC_counterexample")
    print("VERDICT PASS")


if __name__ == "__main__":
    main()
