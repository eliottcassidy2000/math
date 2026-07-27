#!/usr/bin/env python3
"""Dependency-free exact referee for THM-2526.

The script separates three different notions which are easy to conflate:

* a nonzero skew operator after choosing an affine origin;
* a nonzero skew circulant after choosing an oriented slope; and
* a skew operator invariant under the full AGL(1,F_13) gauge.

The first two exist.  The third is exactly zero.  All calculations use
integer/rational arithmetic only.
"""

from fractions import Fraction


P = 13
CHI7 = {1: 1, 2: 1, 3: -1, 4: 1, 5: -1, 6: -1}


def require(condition, label):
    if not condition:
        raise AssertionError(label)


def zero():
    return [[0 for _ in range(P)] for _ in range(P)]


def identity():
    return [[int(row == col) for col in range(P)] for row in range(P)]


def transpose(matrix):
    return [list(row) for row in zip(*matrix)]


def add(first, second):
    return [
        [first[row][col] + second[row][col] for col in range(P)]
        for row in range(P)
    ]


def subtract(first, second):
    return [
        [first[row][col] - second[row][col] for col in range(P)]
        for row in range(P)
    ]


def scale(value, matrix):
    return [[value * entry for entry in row] for row in matrix]


def multiply(first, second):
    return [
        [
            sum(first[row][middle] * second[middle][col] for middle in range(P))
            for col in range(P)
        ]
        for row in range(P)
    ]


def product(*matrices):
    result = identity()
    for matrix in matrices:
        result = multiply(result, matrix)
    return result


def commutator(first, second):
    return subtract(multiply(first, second), multiply(second, first))


def conjugate(pullback, matrix):
    return product(transpose(pullback), matrix, pullback)


def matrix_sum(matrices):
    result = zero()
    for matrix in matrices:
        result = add(result, matrix)
    return result


def rational_rank(matrix):
    work = [[Fraction(entry) for entry in row] for row in matrix]
    pivot_row = 0
    for col in range(P):
        pivot = next(
            (row for row in range(pivot_row, P) if work[row][col]), None
        )
        if pivot is None:
            continue
        work[pivot_row], work[pivot] = work[pivot], work[pivot_row]
        value = work[pivot_row][col]
        work[pivot_row] = [entry / value for entry in work[pivot_row]]
        for row in range(P):
            if row == pivot_row or not work[row][col]:
                continue
            factor = work[row][col]
            work[row] = [
                work[row][j] - factor * work[pivot_row][j] for j in range(P)
            ]
        pivot_row += 1
    return pivot_row


def augmentation_determinant(matrix):
    """Determinant on basis e_1-e_0,...,e_12-e_0 of the augmentation."""
    size = P - 1
    work = [
        [Fraction(matrix[row][col] - matrix[row][0]) for col in range(1, P)]
        for row in range(1, P)
    ]
    determinant = Fraction(1)
    for col in range(size):
        pivot = next((row for row in range(col, size) if work[row][col]), None)
        require(pivot is not None, "augmentation determinant pivot")
        if pivot != col:
            work[col], work[pivot] = work[pivot], work[col]
            determinant = -determinant
        value = work[col][col]
        determinant *= value
        for row in range(col + 1, size):
            if not work[row][col]:
                continue
            factor = work[row][col] / value
            for index in range(col, size):
                work[row][index] -= factor * work[col][index]
    require(determinant.denominator == 1, "integral augmentation determinant")
    return determinant.numerator


def matrix_vector(matrix, vector):
    return [
        sum(matrix[row][col] * vector[col] for col in range(P))
        for row in range(P)
    ]


def is_symmetric(matrix):
    return matrix == transpose(matrix)


def is_skew(matrix):
    return matrix == scale(-1, transpose(matrix))


def translation(amount):
    """(T_a p)_x=p_(x+a)."""
    return [
        [int(col == (row + amount) % P) for col in range(P)]
        for row in range(P)
    ]


def dilation(unit):
    """(D_u p)_x=p_(u x)."""
    return [
        [int(col == unit * row % P) for col in range(P)]
        for row in range(P)
    ]


def hamilton(tau):
    """The symmetric THM-2523 operator A_tau."""
    matrix = zero()
    for row in range(P):
        for s in range(1, 7):
            matrix[row][(row + 2 * tau * s) % P] -= CHI7[s]
            matrix[row][(row - 2 * tau * s) % P] -= CHI7[s]
    return matrix


def oriented_radon(tau):
    """Twice the skew part of the aligned THM-2521 Radon marginal.

    H_tau=sum_s(T_(-2 tau s)-T_(2 tau s)).
    """
    matrix = zero()
    for row in range(P):
        for s in range(1, 7):
            matrix[row][(row - 2 * tau * s) % P] += 1
            matrix[row][(row + 2 * tau * s) % P] -= 1
    return matrix


def paley13():
    residues = {1, 3, 4, 9, 10, 12}
    row = [0] + [1 if value in residues else -1 for value in range(1, P)]
    return [[row[(col - index) % P] for col in range(P)] for index in range(P)]


def audit_affine_commutant_and_halfsets():
    # AGL(1,13) has exactly the diagonal and off-diagonal orbits on X x X.
    maps = [
        tuple((unit * value + shift) % P for value in range(P))
        for unit in range(1, P)
        for shift in range(P)
    ]
    require(len(set(maps)) == P * (P - 1), "affine group size")

    unseen = {(x, y) for x in range(P) for y in range(P)}
    orbits = []
    while unseen:
        seed = next(iter(unseen))
        orbit = {(mapping[seed[0]], mapping[seed[1]]) for mapping in maps}
        unseen -= orbit
        orbits.append(orbit)
    orbit_sizes = sorted(len(orbit) for orbit in orbits)
    require(orbit_sizes == [P, P * (P - 1)], "two pair orbits")

    # Every ordered off-diagonal pair shares its orbit with its transpose.
    require(
        all(
            (right, left)
            in next(orbit for orbit in orbits if (left, right) in orbit)
            for left in range(P)
            for right in range(P)
            if left != right
        ),
        "affine swap lies in off-diagonal orbit",
    )

    # There are 2^6 cyclic half-systems.  Inversion sends every one to its
    # complement, so none is reflection-invariant.
    pairs = tuple((value, P - value) for value in range(1, 7))
    halfsets = []
    for mask in range(1 << 6):
        chosen = {
            pair[(mask >> index) & 1] for index, pair in enumerate(pairs)
        }
        halfsets.append(chosen)
    require(len({tuple(sorted(choice)) for choice in halfsets}) == 64, "halfsets")
    reflection_fixed = sum(choice == {-value % P for value in choice} for choice in halfsets)
    require(reflection_fixed == 0, "no reflection-fixed halfset")
    return orbit_sizes, len(halfsets), reflection_fixed


def audit_origin_skew():
    dilate = dilation(5)
    reflection = dilation(-1 % P)
    require(multiply(dilate, dilate) == reflection, "D5 square is reflection")
    require(transpose(dilate) == product(reflection, dilate), "D5 adjoint")

    origin_skew = subtract(dilate, transpose(dilate))
    require(is_skew(origin_skew), "origin operator is skew")
    require(rational_rank(origin_skew) == 6, "origin operator rank")
    require(
        conjugate(dilation(-1 % P), origin_skew) == origin_skew,
        "origin operator commutes with reflection",
    )
    require(
        all(conjugate(dilation(unit), origin_skew) == origin_skew for unit in range(1, P)),
        "origin operator commutes with all linear gauges",
    )
    translation_failures = sum(
        conjugate(translation(shift), origin_skew) != origin_skew
        for shift in range(1, P)
    )
    require(translation_failures == 12, "every nonzero translation moves origin skew")
    orbit_sum = matrix_sum(
        conjugate(translation(center), origin_skew) for center in range(P)
    )
    require(orbit_sum == zero(), "origin Reynolds average vanishes")

    # Three directed four-cycles away from the fixed origin.
    nonzero_entries = sum(entry != 0 for row in origin_skew for entry in row)
    require(nonzero_entries == 24, "three oriented four-cycles")
    return origin_skew, translation_failures, nonzero_entries


def audit_oriented_slope_skew():
    ranks = []
    determinants = []
    tournament_entries = []
    for tau in range(1, P):
        skew = oriented_radon(tau)
        ranks.append(rational_rank(skew))
        determinants.append(augmentation_determinant(skew))
        require(is_skew(skew), "oriented Radon is skew")
        require(all(sum(row) == 0 for row in skew), "oriented Radon row sum")
        off_diagonal = [
            skew[row][col]
            for row in range(P)
            for col in range(P)
            if row != col
        ]
        require(set(off_diagonal) == {-1, 1}, "oriented Radon is a tournament")
        tournament_entries.append(len(off_diagonal))

        # The exact Cayley-transform identity proves rank 12 without a
        # floating-point Fourier calculation.
        require(
            multiply(add(identity(), translation(tau)), skew)
            == subtract(translation(tau), identity()),
            "Cayley transform identity",
        )
        require(
            all(
                conjugate(translation(shift), skew) == skew
                for shift in range(P)
            ),
            "slope skew is translation-invariant",
        )
        require(oriented_radon(-tau % P) == scale(-1, skew), "slope converse")
        require(
            conjugate(dilation(-1 % P), skew) == scale(-1, skew),
            "reflection flips slope skew",
        )
        require(
            all(
                conjugate(dilation(unit), skew) == oriented_radon(unit * tau % P)
                for unit in range(1, P)
            ),
            "multiplicative slope covariance",
        )
    require(set(ranks) == {12}, "all oriented-slope ranks")
    require(set(determinants) == {13}, "all oriented-slope determinants")
    require(set(tournament_entries) == {P * (P - 1)}, "complete oriented support")
    require(matrix_sum(oriented_radon(tau) for tau in range(1, P)) == zero(), "slope average")
    return ranks[0], determinants[0], oriented_radon(1)[0]


def audit_hamilton_dilation_and_fano_products(origin_skew):
    dilate = dilation(5)
    dilate_star = transpose(dilate)
    reflection = dilation(-1 % P)
    plus = add(dilate, dilate_star)
    minus = origin_skew
    paley = paley13()

    require(is_symmetric(plus), "D+D* symmetric")
    require(is_skew(minus), "D-D* skew")
    require(is_symmetric(paley), "C13 symmetric")
    require(rational_rank(paley) == 12, "C13 augmentation rank")
    require(
        product(hamilton(1), hamilton(2), hamilton(4)) == scale(5, paley),
        "Fano slope product is 5 C13",
    )
    all_ones = [[1 for _ in range(P)] for _ in range(P)]
    require(
        multiply(paley, paley) == subtract(scale(P, identity()), all_ones),
        "C13 square",
    )
    require(conjugate(dilate, paley) == scale(-1, paley), "D5 flips C13")
    require(
        multiply(paley, reflection) == multiply(reflection, paley),
        "C13 commutes with reflection",
    )

    energy_commutators = 0
    for tau in range(1, P):
        energy = hamilton(tau)
        require(is_symmetric(energy), "A_tau symmetric")
        require(rational_rank(energy) == 12, "A_tau augmentation rank")
        require(conjugate(dilate, energy) == scale(-1, energy), "D5 anti-isometry")
        require(multiply(energy, reflection) == multiply(reflection, energy), "A reflection")
        require(commutator(energy, paley) == zero(), "A and C commute")
        for sigma in range(1, P):
            require(commutator(energy, hamilton(sigma)) == zero(), "energy algebra commutes")
            energy_commutators += 1

        even_skew = multiply(energy, plus)
        odd_symmetric = multiply(energy, minus)
        require(is_skew(even_skew), "A(D+D*) skew")
        require(is_symmetric(odd_symmetric), "A(D-D*) symmetric")
        require(rational_rank(even_skew) == 6, "even-sector skew rank")
        require(rational_rank(odd_symmetric) == 6, "odd-sector symmetric rank")

        energy_dilate = multiply(energy, dilate)
        skew_part_twice = subtract(energy_dilate, transpose(energy_dilate))
        symmetric_part_twice = add(energy_dilate, transpose(energy_dilate))
        require(skew_part_twice == even_skew, "skew part of A D5")
        require(symmetric_part_twice == odd_symmetric, "symmetric part of A D5")

        translated_average = matrix_sum(
            conjugate(translation(center), even_skew) for center in range(P)
        )
        require(translated_average == zero(), "origin-skew Hamilton average")

    # The Fano-weighted oriented product is a full-rank skew operator, but
    # it has exactly the same converse sheet as the oriented Radon operator.
    fano_skews = []
    for tau in range(1, P):
        slope_skew = oriented_radon(tau)
        fano_skew = multiply(hamilton(tau), slope_skew)
        fano_skews.append(fano_skew)
        require(is_skew(fano_skew), "A H skew")
        require(rational_rank(fano_skew) == 12, "A H rank")
        require(
            all(
                fano_skew[0][difference] * slope_skew[0][difference] < 0
                for difference in range(1, P)
            ),
            "Fano weights retain only the converse halfset",
        )
    energy = hamilton(1)
    slope_skew = oriented_radon(1)
    fano_skew = fano_skews[0]

    paley_skew = multiply(paley, slope_skew)
    require(is_skew(paley_skew), "C13 H skew")
    require(rational_rank(paley_skew) == 12, "C13 H rank")
    require(
        all(paley_skew[0][difference] != 0 for difference in range(1, P)),
        "C13 H has no off-diagonal ties",
    )

    # A commutator of the two complementary skew carriers exists but spends
    # both their gauges.  Its centre-average vanishes.
    mixed_ranks = []
    for tau in range(1, P):
        current_slope = oriented_radon(tau)
        mixed_commutator = commutator(origin_skew, current_slope)
        require(is_skew(mixed_commutator), "skew-skew commutator")
        mixed_ranks.append(rational_rank(mixed_commutator))
        require(
            conjugate(reflection, mixed_commutator)
            == scale(-1, mixed_commutator),
            "mixed commutator still has converse",
        )
        mixed_average = matrix_sum(
            commutator(
                conjugate(translation(center), origin_skew), current_slope
            )
            for center in range(P)
        )
        require(mixed_average == zero(), "mixed centre average")
    require(set(mixed_ranks) == {12}, "all mixed commutator ranks")

    return (
        energy_commutators,
        fano_skew[0],
        paley_skew[0],
        mixed_ranks[0],
    )


def audit_self_pairing_boundary():
    controls = (
        [1] + [0] * 12,
        list(range(P)),
        [index * index - 3 * index + 1 for index in range(P)],
    )
    operators = (
        subtract(dilation(5), transpose(dilation(5))),
        oriented_radon(1),
        multiply(hamilton(1), oriented_radon(1)),
    )
    checks = 0
    cross_nonzero = 0
    for operator in operators:
        for vector in controls:
            image = [sum(operator[row][col] * vector[col] for col in range(P)) for row in range(P)]
            require(sum(vector[row] * image[row] for row in range(P)) == 0, "skew diagonal")
            checks += 1
        for left in controls:
            for right in controls:
                value = sum(
                    left[row]
                    * operator[row][col]
                    * right[col]
                    for row in range(P)
                    for col in range(P)
                )
                cross_nonzero += value != 0
    require(cross_nonzero > 0, "ordered cross-pairing survives")
    return checks, cross_nonzero


def audit_live_guard_sheet_odd_bank():
    # THM-2198 uses the live sheet coordinate s=-Hk.  If u=-H, the
    # relabelling V_k -> V_s is L=M_(u^-1), so pulling H_1 back to k gives
    # H_(u^-1)=H_(-H^-1).
    gauge_checks = 0
    baseline = oriented_radon(1)
    for guard in range(1, P):
        u = -guard % P
        tau = pow(u, -1, P)
        relabel = dilation(tau)
        require(
            conjugate(relabel, baseline) == oriented_radon(tau),
            "live reversed-guard sheet gauge",
        )
        gauge_checks += 1

    # Exhaust all nonconstant Boolean predecessor profiles as hostile finite
    # controls.  Scale the centred autocorrelation by 13^2:
    #   b_scaled(t)=13 sum_r q_r q_(r+t)-(sum_r q_r)^2.
    # In the s-coordinate tau=1.  The live odd bank H_1 A_1 b_scaled is
    # nonzero, reflection-odd, and has all twelve primitive modes.
    energy = hamilton(1)
    skew = baseline
    profile_count = 0
    primitive_mode_checks = 0
    positive_pair_controls = 0
    for mask in range(1, (1 << P) - 1):
        profile = [(mask >> index) & 1 for index in range(P)]
        collision = [
            sum(profile[index] * profile[(index + shift) % P] for index in range(P))
            for shift in range(P)
        ]
        total = sum(profile)
        centred_scaled = [P * value - total * total for value in collision]
        require(sum(centred_scaled) == 0, "centred Boolean collision")
        require(
            all(centred_scaled[-index % P] == centred_scaled[index] for index in range(P)),
            "Boolean collision is even",
        )
        cayley_bank = matrix_vector(skew, centred_scaled)
        require(any(cayley_bank), "live Cayley bank nonzero")
        require(
            all(cayley_bank[-index % P] == -cayley_bank[index] for index in range(P)),
            "live Cayley bank is reflection-odd",
        )
        translated_fano = matrix_vector(energy, centred_scaled)
        odd_bank = matrix_vector(skew, translated_fano)
        require(any(odd_bank), "live odd bank nonzero")
        require(sum(odd_bank) == 0, "live odd bank centred")
        require(
            all(odd_bank[-index % P] == -odd_bank[index] for index in range(P)),
            "live bank is reflection-odd",
        )
        require(
            any(odd_bank[index] > 0 and odd_bank[-index % P] < 0 for index in range(1, P)),
            "live odd bank has an antipodal sign pair",
        )
        positive_pair_controls += 1
        # A rational polynomial of degree at most 12 which vanishes at one
        # primitive thirteenth root is a scalar multiple of Phi_13, hence
        # has all coefficients equal.  The nonzero centred odd vector above
        # is not constant, so all twelve conjugate primitive modes survive.
        require(len(set(odd_bank)) > 1, "cyclotomic primitive-mode certificate")
        primitive_mode_checks += P - 1
        profile_count += 1
    require(profile_count == (1 << P) - 2, "all nonconstant Boolean profiles")
    return gauge_checks, profile_count, primitive_mode_checks, positive_pair_controls


def main():
    orbit_sizes, halfsets, reflection_fixed = audit_affine_commutant_and_halfsets()
    origin_skew, translation_failures, origin_support = audit_origin_skew()
    slope_rank, slope_determinant, slope_row = audit_oriented_slope_skew()
    energy_commutators, fano_row, paley_row, mixed_rank = (
        audit_hamilton_dilation_and_fano_products(origin_skew)
    )
    diagonal_checks, cross_nonzero = audit_self_pairing_boundary()
    guard_gauges, boolean_profiles, odd_modes, positive_pairs = (
        audit_live_guard_sheet_odd_bank()
    )

    print("THM-2526 affine skew-orientation boundary: PASS")
    print(
        "affine_pair_orbits=" + ",".join(map(str, orbit_sizes))
        + f"; cyclic_halfsets={halfsets}; reflection_fixed={reflection_fixed}"
    )
    print(
        f"origin_skew_rank=6; support={origin_support}; "
        f"nonzero_translation_failures={translation_failures}; centre_average=0"
    )
    print(
        f"oriented_slope_rank={slope_rank}; determinant={slope_determinant}; "
        "cayley_identity=(I+T_tau)H_tau=T_tau-I"
    )
    print("oriented_slope_first_row=" + ",".join(map(str, slope_row)))
    print(
        "A_tau(D5+D5*)=skew_rank_6; "
        "A_tau(D5-D5*)=symmetric_rank_6"
    )
    print(
        f"commuting_energy_checks={energy_commutators}; "
        f"mixed_skew_commutator_rank={mixed_rank}; all_Reynolds_averages=0"
    )
    print("Fano_oriented_first_row=" + ",".join(map(str, fano_row)))
    print("Paley_oriented_first_row=" + ",".join(map(str, paley_row)))
    print(
        f"skew_diagonal_checks={diagonal_checks}; "
        f"nonzero_ordered_cross_controls={cross_nonzero}"
    )
    print(
        f"live_guard_gauges={guard_gauges}; Boolean_profiles={boolean_profiles}; "
        f"odd_primitive_modes={odd_modes}; antipodal_sign_pairs={positive_pairs}"
    )
    print(
        "full_AGL_invariant_skew_dimension=0; live sheet orientation retained; "
        "Boolean owner emission not inferred"
    )


if __name__ == "__main__":
    main()
