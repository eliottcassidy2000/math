#!/usr/bin/env python3
"""Exact algebra audit for the proposed THM-2792 cyclic intertwiner."""

from fractions import Fraction


P = 13
ENDPOINT_DIRECTIONS = P * P - 1
ENDPOINT_CYCLES = ENDPOINT_DIRECTIONS * P
NONZERO_ENDPOINT_MODES = ENDPOINT_CYCLES * (P - 1)
ACTUAL_CHORD_EXPONENT = 53028
ACTUAL_TRANSgression = 4079


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def convolution(left, right):
    return tuple(
        sum(
            left[index] * right[(degree - index) % P]
            for index in range(P)
        )
        for degree in range(P)
    )


def delta(exponent):
    return tuple(1 if index == exponent % P else 0 for index in range(P))


def scalar(value, row):
    return tuple(value * entry for entry in row)


def add(*rows):
    return tuple(sum(row[index] for row in rows) for index in range(P))


def automorphism(row, exponent):
    """The group-algebra map z -> z^exponent."""
    output = [0] * P
    for index, value in enumerate(row):
        output[index * exponent % P] += value
    return tuple(output)


def bareiss_determinant(matrix):
    work = [list(row) for row in matrix]
    size = len(work)
    sign = 1
    previous = 1
    for pivot_index in range(size - 1):
        if work[pivot_index][pivot_index] == 0:
            swap = next(
                (
                    row for row in range(pivot_index + 1, size)
                    if work[row][pivot_index] != 0
                ),
                None,
            )
            require(swap is not None, "singular Bareiss pivot")
            work[pivot_index], work[swap] = work[swap], work[pivot_index]
            sign *= -1
        pivot = work[pivot_index][pivot_index]
        for row in range(pivot_index + 1, size):
            for column in range(pivot_index + 1, size):
                numerator = (
                    work[row][column] * pivot
                    - work[row][pivot_index]
                    * work[pivot_index][column]
                )
                require(
                    numerator % previous == 0,
                    "Bareiss division stopped being exact",
                )
                work[row][column] = numerator // previous
        previous = pivot
    return sign * work[-1][-1]


def convolution_matrix(kernel):
    return tuple(
        tuple(kernel[(row - column) % P] for column in range(P))
        for row in range(P)
    )


def remove_row_column(matrix, removed_row, removed_column):
    return tuple(
        tuple(
            value
            for column, value in enumerate(row)
            if column != removed_column
        )
        for row_index, row in enumerate(matrix)
        if row_index != removed_row
    )


def main():
    identity = delta(0)
    z = delta(1)
    one_plus_z = add(identity, z)
    alternating = tuple(
        Fraction(1 if index % 2 == 0 else -1, 2)
        for index in range(P)
    )
    require(
        convolution(one_plus_z, alternating) == identity
        and convolution(alternating, one_plus_z) == identity,
        "alternating inverse of 1+z failed",
    )

    a_normalized = convolution(delta(6), one_plus_z)
    a_inverse = convolution(delta(-6), alternating)
    require(
        convolution(a_normalized, a_inverse) == identity
        and convolution(a_inverse, a_normalized) == identity,
        "inverse of z^6(1+z) failed",
    )

    matrix = convolution_matrix(one_plus_z)
    determinant = bareiss_determinant(matrix)
    cofactor = bareiss_determinant(remove_row_column(matrix, 0, 0))
    require(
        (determinant, cofactor) == (2, 1),
        "integral determinant/cofactor boundary changed",
    )

    # Exact normalized integrality criterion:
    # J*A_0^{-1} is integral iff J(1) is even.  Exhaust the Boolean cube.
    integral_boolean_intertwiners = 0
    half_integral_boolean_intertwiners = 0
    for mask in range(1 << P):
        boolean_j = tuple(
            (mask >> index) & 1
            for index in range(P)
        )
        boolean_t = convolution(boolean_j, a_inverse)
        is_integral = all(value.denominator == 1 for value in boolean_t)
        even_zero_mode = sum(boolean_j) % 2 == 0
        require(
            is_integral == even_zero_mode,
            "normalized integral-intertwiner criterion failed",
        )
        integral_boolean_intertwiners += is_integral
        half_integral_boolean_intertwiners += not is_integral
    require(
        (
            integral_boolean_intertwiners,
            half_integral_boolean_intertwiners,
        )
        == (4096, 4096),
        "Boolean normalized-integrality census changed",
    )

    # A nontrivial rational unit control J=2+z.
    j_control = add(scalar(2, identity), z)
    j_inverse_numerator = tuple(
        (-1) ** index * 2 ** (P - 1 - index)
        for index in range(P)
    )
    j_inverse = scalar(Fraction(1, 2**P + 1), j_inverse_numerator)
    require(
        convolution(j_control, j_inverse) == identity,
        "2+z control inverse failed",
    )
    t_control = convolution(j_control, a_inverse)
    t_control_inverse = convolution(a_normalized, j_inverse)
    require(
        convolution(t_control, a_normalized) == j_control
        and convolution(t_control, t_control_inverse) == identity,
        "unit intertwiner control failed",
    )

    # Independent source/target origin changes are coordinate changes of the
    # same abstract module map.  The kernel acquires only their relative
    # monomial z^(h_A-h_E); common re-origining cancels.
    origin_pair_checks = 0
    for endpoint_shift in range(P):
        endpoint_origin = delta(-endpoint_shift)
        j_shifted = convolution(endpoint_origin, j_control)
        for semantic_shift in range(P):
            semantic_origin = delta(-semantic_shift)
            a_shifted = convolution(semantic_origin, a_normalized)
            a_shifted_inverse = convolution(
                a_inverse,
                delta(semantic_shift),
            )
            coordinate_t = convolution(j_shifted, a_shifted_inverse)
            expected_t = convolution(
                delta(semantic_shift - endpoint_shift),
                t_control,
            )
            require(
                coordinate_t == expected_t
                and convolution(coordinate_t, a_shifted) == j_shifted,
                "independent-origin coordinate covariance failed",
            )
            if semantic_shift == endpoint_shift:
                require(
                    coordinate_t == t_control,
                    "common-origin covariance failed",
                )
            if semantic_shift == 0 and endpoint_shift:
                require(
                    coordinate_t != t_control,
                    "endpoint basis change left one coordinate row fixed",
                )
            origin_pair_checks += 1
    require(origin_pair_checks == P * P, "origin-pair census changed")

    # A pointwise fixed-action torsor identification linearizes to a
    # monomial.  Even after a scalar weight it preserves the two-point
    # support of A, while the actual endpoint J has inherited support 13.
    monomial_support_checks = 0
    for shift in range(P):
        transported = convolution(delta(shift), a_normalized)
        require(
            sum(value != 0 for value in transported) == 2,
            "monomial transport changed the semantic support size",
        )
        monomial_support_checks += 1
    require(monomial_support_checks == P, "monomial support census changed")

    # Ordered-cone sharpening.  If nonnegative kernels B,C satisfy
    # B*C=delta_0, then every x in supp(B), y in supp(C) must have x+y=0:
    # there is no cancellation away from zero.  Exhaust the possible
    # nonempty support of B; only a singleton leaves any possible support
    # for C.  Positive values are then reciprocal scalar monomials.
    bipositive_support_patterns = 0
    for mask in range(1, 1 << P):
        support = tuple(
            index
            for index in range(P)
            if (mask >> index) & 1
        )
        compatible_inverse_support = set(range(P))
        for index in support:
            compatible_inverse_support.intersection_update(
                {(-index) % P}
            )
        if compatible_inverse_support:
            require(
                len(support) == 1
                and len(compatible_inverse_support) == 1,
                "a nonsingleton support admitted a nonnegative inverse",
            )
            bipositive_support_patterns += 1
    require(
        bipositive_support_patterns == P,
        "bi-positive support-pattern census changed",
    )

    # A formal common coefficient reindexing is covariant.  It must not be
    # identified with changing the physical endpoint increment s to u*s.
    for unit in range(1, P):
        inverse = pow(unit, -1, P)
        a_reparametrized = automorphism(a_normalized, inverse)
        j_reparametrized = automorphism(j_control, inverse)
        t_reparametrized = convolution(
            j_reparametrized,
            automorphism(a_inverse, inverse),
        )
        require(
            t_reparametrized == automorphism(t_control, inverse),
            "common cycle reparametrization failed",
        )

    # Sharp formal frame hostile: independently tag the two endpoints by
    # unique prime powers.  For every nonidentity scaling u, the physical
    # u-step edge differs from every formal u-reindexing even after an
    # arbitrary cycle-origin shift.  This rules out a universal algebraic
    # identity; it is not a computation on the canonical endpoint current.
    endpoint_left = tuple(2**index for index in range(P))
    endpoint_right = tuple(3**index for index in range(P))
    step_one = tuple(
        endpoint_left[(index + 1) % P] * endpoint_right[index]
        for index in range(P)
    )
    frame_hostile_checks = 0
    for unit in range(2, P):
        step_unit = tuple(
            endpoint_left[(unit * index + unit) % P]
            * endpoint_right[unit * index % P]
            for index in range(P)
        )
        formal_reindex = tuple(
            step_one[unit * index % P]
            for index in range(P)
        )
        for shift in range(P):
            shifted_reindex = tuple(
                formal_reindex[(index + shift) % P]
                for index in range(P)
            )
            require(
                step_unit != shifted_reindex,
                "physical frame scaling became a formal reindexing "
                "after an origin shift",
            )
            frame_hostile_checks += 1
    require(
        frame_hostile_checks == (P - 2) * P,
        "formal frame-hostile census changed",
    )

    # The unit hypotheses alone cannot force positivity.
    positive_hostile_j = identity
    positive_hostile_t = a_inverse
    positive_count = sum(value > 0 for value in positive_hostile_t)
    negative_count = sum(value < 0 for value in positive_hostile_t)
    require(
        (positive_count, negative_count) == (7, 6)
        and convolution(positive_hostile_t, a_normalized)
        == positive_hostile_j,
        "sharp alternating-sign hostile changed",
    )

    # Quotient transfer forgets the first lower-central transgression.
    require(
        ACTUAL_CHORD_EXPONENT
        == 1 + P * ACTUAL_TRANSgression
        and ACTUAL_TRANSgression % P == 10,
        "THM-2791 chord arithmetic changed",
    )
    near_support_mod_p = tuple(sorted((0, 1)))
    actual_support_mod_p = tuple(
        sorted((0, ACTUAL_CHORD_EXPONENT % P))
    )
    require(
        near_support_mod_p == actual_support_mod_p == (0, 1),
        "coarse transfer stopped identifying the two graded lifts",
    )

    total_endpoint_modes = ENDPOINT_CYCLES * P
    require(
        ENDPOINT_CYCLES == 2184
        and NONZERO_ENDPOINT_MODES == 26208
        and total_endpoint_modes == 28392,
        "endpoint-cycle Fourier census changed",
    )

    print("THM-2792 CYCLIC UNIT INTERTWINER SCRATCH AUDIT")
    print(
        "dependency_invoice="
        "THM2625 zero_modes=2184; "
        "THM2790 nonzero_modes=26208; total_modes=28392"
    )
    print(
        "split_group_algebra="
        "K[C13] -> K^13 by evaluation at all 13 characters; "
        "unit iff every Fourier coordinate is nonzero"
    )
    print(
        "A_normalized=z^6(1+z); "
        "A_inverse=(1/2)z^7(1-z+z^2-...+z^12)"
    )
    print(
        f"integral_boundary=det(I+shift)={determinant}; "
        f"one_12x12_minor={cofactor}; Smith=(1^12,2)"
    )
    print(
        "integrality_criterion=for integral J, "
        "J*A_normalized^(-1) is integral iff J(1) is divisible by 2; "
        "Boolean census integral/half=4096/4096"
    )
    print(
        "intertwiner=the cyclic generators define the coordinate-free "
        "module isomorphism F(r*A)=r*J; in chosen origins its unique "
        "circulant kernel is T=J*A_inverse"
    )
    print(
        f"origin_law={origin_pair_checks} source/target origin pairs obey "
        "T' = z^(h_A-h_E)T; endpoint-only row change is covariance of "
        "the same abstract F"
    )
    print(
        "set_level_boundary=all "
        f"{monomial_support_checks} scalar-monomial transports preserve "
        "support(A)=2, while THM2790 gives inherited support(J)=13; "
        "F is not a pointwise torsor identification"
    )
    print(
        "bipositive_boundary=only "
        f"{bipositive_support_patterns} singleton support patterns admit "
        "a coefficientwise-nonnegative inverse; bi-positive circulant "
        "isomorphisms are positive scalar monomials, while A_normalized "
        "is the sharp one-way-positive/signed-inverse control"
    )
    print(
        "coordinate_law=a formal common polynomial reindexing "
        "z->z^(u^-1) carries T covariantly"
    )
    print(
        "frame_hostile=tagged P_j=2^j,Q_j=3^j reject all "
        f"{frame_hostile_checks} nonidentity-scaling/origin pairs; "
        "there is no formal frame-reindexing identity, but this is not "
        "a computation on the canonical current"
    )
    print(
        "positive_hostile=J=delta0 has all 13 modes but "
        f"T=A_inverse has signs +:{positive_count} -:{negative_count}"
    )
    print(
        "lower_central_hostile="
        "1+u and 1+u^53028 have the same transfer 1+z, "
        "but the latter has first transgression 4079 mod13=10"
    )
    print(
        "scope=origin-free characteristic-zero module isomorphism after "
        "choosing the oriented group action; no integral/bi-positive/"
        "set-level realization, one-way positivity decision, global "
        "frame-natural family, physical descent, or ancestry intertwiner"
    )
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
