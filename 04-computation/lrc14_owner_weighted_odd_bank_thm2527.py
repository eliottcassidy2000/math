#!/usr/bin/env python3
"""Dependency-free exact referee for THM-2527.

The finite part of the theorem is a Boolean cut inequality on one
thirteen-point root fibre.  All calculations below use integers or
``Fraction``.  In particular, no floating Fourier approximation is used.
"""

from collections import Counter
from fractions import Fraction


P = 13
CHI7 = {1: 1, 2: 1, 3: -1, 4: 1, 5: -1, 6: -1}


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def translation(a):
    """Matrix of (T_a x)_r=x_(r+a)."""
    return [[int(c == (r + a) % P) for c in range(P)] for r in range(P)]


def zero_matrix():
    return [[0 for _ in range(P)] for _ in range(P)]


def add_scaled(target, scalar, source):
    for r in range(P):
        for c in range(P):
            target[r][c] += scalar * source[r][c]


def matmul(left, right):
    return [
        [sum(left[r][j] * right[j][c] for j in range(P)) for c in range(P)]
        for r in range(P)
    ]


def matvec(matrix, vector):
    return [sum(matrix[r][c] * vector[c] for c in range(P)) for r in range(P)]


def transpose(matrix):
    return [list(row) for row in zip(*matrix)]


def operator_a(tau):
    matrix = zero_matrix()
    for s in range(1, 7):
        add_scaled(matrix, -CHI7[s], translation(2 * tau * s))
        add_scaled(matrix, -CHI7[s], translation(-2 * tau * s))
    return matrix


def operator_h(tau):
    matrix = zero_matrix()
    for s in range(1, 7):
        add_scaled(matrix, 1, translation(-2 * tau * s))
        add_scaled(matrix, -1, translation(2 * tau * s))
    return matrix


def correlation(mask):
    bits = [(mask >> r) & 1 for r in range(P)]
    return [
        sum(bits[r] * bits[(r + u) % P] for r in range(P))
        for u in range(P)
    ]


def rotate_mask(mask, shift):
    return sum(
        ((mask >> r) & 1) << ((r + shift) % P)
        for r in range(P)
    )


def cyclic_representative(mask):
    return min(rotate_mask(mask, shift) for shift in range(P))


def has_cyclic_zero_run(mask, length=3):
    return any(
        all(((mask >> ((start + j) % P)) & 1) == 0 for j in range(length))
        for start in range(P)
    )


def rational_rank(matrix):
    work = [[Fraction(value) for value in row] for row in matrix]
    rows = len(work)
    cols = len(work[0])
    rank = 0
    for col in range(cols):
        pivot = next((r for r in range(rank, rows) if work[r][col]), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        scale = work[rank][col]
        work[rank] = [value / scale for value in work[rank]]
        for r in range(rows):
            if r == rank or not work[r][col]:
                continue
            scale = work[r][col]
            work[r] = [
                work[r][c] - scale * work[rank][c] for c in range(cols)
            ]
        rank += 1
    return rank


def main():
    a1 = operator_a(1)
    h1 = operator_h(1)
    f1 = matmul(a1, h1)

    require(a1 == transpose(a1), "A_1 must be symmetric")
    require(h1 == [[-x for x in row] for row in transpose(h1)], "H_1 skew")
    require(matmul(a1, h1) == matmul(h1, a1), "A_1 and H_1 commute")
    require(rational_rank(f1) == 12, "A_1 H_1 augmentation rank")
    expected_row = [0, -1, 3, -5, 7, -7, 5, -5, 7, -7, 5, -3, 1]
    require(f1[0] == expected_row, "A_1 H_1 first row drift")
    require(sum(value * value for value in expected_row) == 316, "row norm")

    # Exact inverse-coordinate certificate.  If O=13 A_1 H_1 b and sum b=0,
    # this row recovers b_0.  The harmless constant -65 disappears on U.
    inverse_numerator = [
        0,
        -17,
        -37,
        -37,
        -33,
        -15,
        1,
        -1,
        15,
        33,
        37,
        37,
        17,
    ]
    certificate = [
        13 * sum(inverse_numerator[r] * f1[r][c] for r in range(P))
        for c in range(P)
    ]
    require(certificate == [780] + [-65] * 12, "inverse row certificate")
    require(sum(abs(value) for value in inverse_numerator) == 280, "inverse l1")
    require(sum(value * value for value in inverse_numerator) == 8684, "inverse l2")

    distribution = Counter()
    equality_masks = []
    zero_masks = []
    guard_masks = 0
    guard_equality_masks = 0

    for mask in range(1 << P):
        c = correlation(mask)
        fc = matvec(f1, c)
        n = mask.bit_count()

        # This is the fixed Boolean positive coordinate.  The displayed
        # correlation formula is independently checked against A_1 H_1.
        psi = (
            7 * c[0]
            - 12 * c[1]
            + 8 * c[2]
            - 6 * c[3]
            + 7 * c[4]
            - 6 * c[5]
            + 2 * c[6]
        )
        require(psi == -fc[4] == fc[-4], "fixed-coordinate identity")
        require(all(fc[-t] == -fc[t] for t in range(P)), "odd fibre bank")
        cut_lengths = [
            sum(
                (((mask >> r) & 1) - ((mask >> ((r + d) % P)) & 1)) ** 2
                for r in range(P)
            )
            for d in range(1, 7)
        ]
        require(
            2 * psi
            == 12 * cut_lengths[0]
            - 8 * cut_lengths[1]
            + 6 * cut_lengths[2]
            - 7 * cut_lengths[3]
            + 6 * cut_lengths[4]
            - 2 * cut_lengths[5],
            "cyclic cut-polytope form",
        )
        require(psi >= 0, "Boolean cut inequality")
        require(42 * psi >= n * (P - n), "relative Boolean cut floor")
        require(sum(int(psi >= j) for j in range(1, 99)) == psi, "layer factor")

        distribution[psi] += 1
        if psi == 0:
            zero_masks.append(mask)
        if 42 * psi == n * (P - n) and 0 < n < P:
            equality_masks.append(mask)

        if mask and has_cyclic_zero_run(mask):
            guard_masks += 1
            if 42 * psi == n * (P - n):
                guard_equality_masks += 1

    require(zero_masks == [0, (1 << P) - 1], "only constant masks vanish")
    require(min(value for value in distribution if value) == 1, "positive minimum")
    require(max(distribution) == 98, "positive maximum")
    require(len(equality_masks) == 52, "sharp-floor equality count")
    require(
        {mask.bit_count() for mask in equality_masks} == {6, 7},
        "sharp floor occurs at balanced masks",
    )
    require(guard_masks == 5434, "three-zero-run mask count")
    require(guard_equality_masks == 52, "guard class retains all equality masks")
    equality_necklaces = sorted({cyclic_representative(mask) for mask in equality_masks})
    require(equality_necklaces == [399, 463, 483, 487], "four equality necklaces")

    # Booleanity is load-bearing.  This rational [0,1]-valued root vector
    # is a minimal hostile control for the fixed-sign assertion.
    non_boolean = [Fraction(1, 2)] + [Fraction(0)] * 4 + [Fraction(1)] * 4 + [Fraction(0)] * 4
    non_boolean_correlation = [
        sum(non_boolean[r] * non_boolean[(r + u) % P] for r in range(P))
        for u in range(P)
    ]
    non_boolean_psi = (
        7 * non_boolean_correlation[0]
        - 12 * non_boolean_correlation[1]
        + 8 * non_boolean_correlation[2]
        - 6 * non_boolean_correlation[3]
        + 7 * non_boolean_correlation[4]
        - 6 * non_boolean_correlation[5]
        + 2 * non_boolean_correlation[6]
    )
    require(non_boolean_psi == Fraction(-1, 4), "non-Boolean sign boundary")

    # Gauge covariance: for every slope, -4 tau is the same positive
    # coordinate after multiplicative relabelling.  Exhausting all masks a
    # second time would be redundant; checking the twelve matrix conjugates
    # on the 13 delta correlations fixes the full linear identity.
    for tau in range(1, P):
        ftau = matmul(operator_a(tau), operator_h(tau))
        for d in range(7):
            even_delta = [0] * P
            even_delta[d] = 1
            if d:
                even_delta[-d] = 1
            relabelled = [0] * P
            # c'_(u)=c_(tau u) puts F_tau at chart slope one.
            for u in range(P):
                relabelled[u] = even_delta[(tau * u) % P]
            lhs = matvec(ftau, even_delta)[(-4 * tau) % P]
            rhs = matvec(f1, relabelled)[-4]
            require(lhs == rhs, "slope covariance at positive coordinate")

    # A small exact weighted-disintegration control.  The 1/13 branch
    # normalization cancels the leading 13 in O=13 A H b.
    sample_masks = [1, 0b11, 0b110001111, 0b1010101010101]
    sample_weights = [Fraction(2, 7), Fraction(3, 11), Fraction(5, 13), Fraction(7, 17)]
    b_profile = [
        sum(weight * correlation(mask)[u] for mask, weight in zip(sample_masks, sample_weights))
        / P
        for u in range(P)
    ]
    mean_b = sum(b_profile) / P
    centred = [value - mean_b for value in b_profile]
    odd_bank = [P * value for value in matvec(f1, centred)]
    weighted_drift = sum(
        weight * mask.bit_count() * (P - mask.bit_count()) / (P * P)
        for mask, weight in zip(sample_masks, sample_weights)
    )
    require(centred[0] == weighted_drift, "weighted total-defect disintegration")
    integrated_psi = sum(
        weight * (-matvec(f1, correlation(mask))[4])
        for mask, weight in zip(sample_masks, sample_weights)
    )
    require(odd_bank[-4] == integrated_psi, "same-integrand owner coupling")
    require(42 * odd_bank[-4] >= 169 * centred[0], "weighted 169/42 floor")

    print("THM-2527 owner-weighted odd bank: PASS")
    print("operator_rank=12; AH_row_l2_squared=316")
    print("inverse_coordinate_l1=56/169; inverse_coordinate_l2_squared=8684/845^2")
    print("boolean_masks=8192; zero_masks=2; min_positive_psi=1; max_psi=98")
    print("sharp_cut_floor=13/42; equality_masks=52; equality_sizes=6,7")
    print("equality_C13_necklaces=399,463,483,487; nonboolean_control=-1/4")
    print("three_zero_run_masks=5434; sharp_guard_equalities=52")
    print("fixed_positive_coordinate=-4*tau; owner_floor=(169/42)*weighted_drift")
    print("rational total drift plus Phi_13 irreducibility gives all twelve modes")


if __name__ == "__main__":
    main()
