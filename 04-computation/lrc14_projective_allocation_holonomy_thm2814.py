#!/usr/bin/env python3
"""Exact companion for THM-2814.

The script separates two non-idempotent allocation branches:

* a physically normalized common-line Segre square, whose additive face can
  be nonzero even though its projective cross-ratio is one;
* a four-line square modulo row/column rephasing, whose complete invariant is
  projective square holonomy.

It also checks idempotent provenance, the THM-2593 charged coefficient atlas,
the THM-2716 quarter-turn control, and the exact THM-2779/2806/2807 constants.
It is dependency-free and uses no Python ``assert`` statement.
"""

from collections import Counter
from itertools import product
from math import comb


P = 13
Q = 7
FIELD = 53


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def mobius(vector, modulus=None):
    value = vector[0] - vector[1] - vector[2] + vector[3]
    return value if modulus is None else value % modulus


def cross_ratio(vector, modulus):
    denominator = vector[1] * vector[2] % modulus
    require(denominator != 0, "cross-ratio denominator vanished")
    return (
        vector[0] * vector[3] * pow(denominator, -1, modulus)
    ) % modulus


def normalize_square(vector, modulus):
    """Normalize (v00,v10,v01,v11) by source-row/target-column gauge."""
    v00, v10, v01, v11 = vector
    require(all(value % modulus for value in vector), "zero square entry")
    row0 = pow(v00, -1, modulus)
    row1 = pow(v10, -1, modulus)
    col0 = 1
    col1 = v00 * pow(v01, -1, modulus) % modulus
    normalized = (
        row0 * col0 * v00 % modulus,
        row1 * col0 * v10 % modulus,
        row0 * col1 * v01 % modulus,
        row1 * col1 * v11 % modulus,
    )
    return normalized, (row0, row1, col0, col1)


def allocation_spectrum(vector, modulus):
    out = []
    for source_character, target_character in product(range(2), repeat=2):
        value = 0
        for source_flag, target_flag in product(range(2), repeat=2):
            index = 2 * source_flag + target_flag
            exponent = (
                source_character * source_flag
                + target_character * target_flag
            )
            value += (-1 if exponent % 2 else 1) * vector[index]
        out.append(value % modulus)
    return tuple(out)


def idempotent_and_provenance_probe():
    idempotents = tuple(
        value for value in range(P) if value * value % P == value
    )
    require(idempotents == (0, 1), "field idempotents changed")

    common_faces = []
    for alpha, beta in product(idempotents, repeat=2):
        vector = (1, alpha, beta, alpha * beta % P)
        if all(vector):
            common_faces.append(mobius(vector, P))
    require(common_faces == [0], "idempotent common atom became curved")

    contrast_pairs = []
    for alpha, beta in product(range(1, P), repeat=2):
        vector = (1, alpha, beta, alpha * beta % P)
        require(cross_ratio(vector, P) == 1, "Segre square gained holonomy")
        require(
            mobius(vector, P) == (1 - alpha) * (1 - beta) % P,
            "Segre mixed-face factorization changed",
        )
        if mobius(vector, P):
            contrast_pairs.append((alpha, beta))
    require(
        len(contrast_pairs) == (P - 2) ** 2 == 121,
        "normalized common-line contrast census changed",
    )

    checks = 0
    face_histogram = Counter()
    for left_bits in range(1 << 4):
        left = tuple((left_bits >> index) & 1 for index in range(4))
        for right_bits in range(1 << 4):
            right = tuple((right_bits >> index) & 1 for index in range(4))
            for weights in product((-1, 0, 1), repeat=4):
                vector = (
                    sum(weights),
                    sum(weights[i] * left[i] for i in range(4)),
                    sum(weights[i] * right[i] for i in range(4)),
                    sum(
                        weights[i] * left[i] * right[i]
                        for i in range(4)
                    ),
                )
                face = mobius(vector)
                complement = sum(
                    weights[i] * (1 - left[i]) * (1 - right[i])
                    for i in range(4)
                )
                require(face == complement, "idempotent provenance changed")
                face_histogram[face] += 1
                checks += 1
    require(checks == 16 * 16 * 3**4 == 20736, "probe universe changed")

    quotient_hostile = (4, 2, 2, 1)
    require(mobius(quotient_hostile) == 1, "coarse hostile changed")
    return idempotents, len(contrast_pairs), checks, quotient_hostile


def projective_square_probe():
    counts = Counter()
    for vector in product(range(1, P), repeat=4):
        kappa = cross_ratio(vector, P)
        normalized, _gauge = normalize_square(vector, P)
        require(
            normalized == (1, 1, 1, kappa),
            "row/column normal form changed",
        )
        counts[kappa] += 1
    require(
        counts == Counter({value: (P - 1) ** 3 for value in range(1, P)}),
        "cross-ratio fibre census changed",
    )

    gauge_checks = 0
    for kappa in range(1, P):
        representative = (1, 1, 1, kappa)
        for row0, row1, col0, col1 in product(range(1, P), repeat=4):
            transformed = (
                row0 * col0 % P,
                row1 * col0 % P,
                row0 * col1 % P,
                row1 * col1 * kappa % P,
            )
            require(
                cross_ratio(transformed, P) == kappa,
                "cross-ratio lost gauge invariance",
            )
            gauge_checks += 1
        require(cross_ratio(representative, P) == kappa, "bad representative")
    require(
        gauge_checks == (P - 1) ** 5 == 248832,
        "row/column gauge universe changed",
    )

    rank_one_chart = (1, 2, 3, 6)
    normalized, _gauge = normalize_square(rank_one_chart, P)
    require(cross_ratio(rank_one_chart, P) == 1, "rank-one hostile changed")
    require(mobius(rank_one_chart, P) == 2, "rank-one chart face changed")
    require(normalized == (1, 1, 1, 1), "rank-one chart did not flatten")

    flat = (1, 1, 1, 1)
    signed = (1, -1, -1, 1)
    require(
        allocation_spectrum(flat, P) == (4, 0, 0, 0)
        and allocation_spectrum(signed, P) == (0, 0, 0, 4),
        "allocation Fourier sign control changed",
    )
    require(cross_ratio(tuple(value % P for value in signed), P) == 1,
            "external sign acquired joint holonomy")
    return tuple(sorted(counts.items())), gauge_checks, rank_one_chart, normalized


def primitive_thirteenth_root():
    for candidate in range(2, FIELD):
        if pow(candidate, P, FIELD) == 1 and candidate != 1:
            require(
                all(
                    pow(candidate, exponent, FIELD) != 1
                    for exponent in range(1, P)
                ),
                "root was not primitive",
            )
            return candidate
    raise RuntimeError("no primitive thirteenth root")


ZETA = primitive_thirteenth_root()


def character_and_quarter_turn_probe():
    faces = []
    for origin in range(P):
        alpha = pow(ZETA, origin, FIELD)
        vector = (1, alpha, alpha, alpha * alpha % FIELD)
        require(cross_ratio(vector, FIELD) == 1, "edge character gained cocycle")
        require(
            mobius(vector, FIELD) == (1 - alpha) ** 2 % FIELD,
            "edge-character face changed",
        )
        faces.append(mobius(vector, FIELD))
    require(faces[0] == 0 and all(faces[1:]), "character origin census changed")

    representative_factors = (
        1,
        ZETA,
        pow(ZETA, -1, FIELD),
        1,
    )
    require(
        len(set(representative_factors)) > 1
        and cross_ratio(representative_factors, FIELD) == 1,
        "representative edge gauge changed",
    )
    normalized, _gauge = normalize_square(representative_factors, FIELD)
    require(normalized == (1, 1, 1, 1), "character gauge did not flatten")

    vector = (1, 0)
    turn = lambda pair: (-pair[1], pair[0])
    once = turn(vector)
    twice = turn(once)
    mixed = (
        vector[0] - 2 * once[0] + twice[0],
        vector[1] - 2 * once[1] + twice[1],
    )
    scalar_readout = (sum(vector), sum(once), sum(once), sum(twice))
    require(
        twice == (-1, 0) and mixed == (0, -2),
        "quarter-turn contrast changed",
    )
    require(all(scalar_readout) and mobius(scalar_readout) == -2,
            "quarter-turn readout failed")
    return tuple(faces), representative_factors, normalized, mixed, scalar_readout


# THM-2593/2585 charged target rows before reduction modulo Phi_7.
Y_RAW = (
    (0, 9, 9, 0, 0, 4, 4),
    (5, 9, 7, 11, 11, 4, 9),
    (3, 11, 7, 10, 10, 11, 1),
    (9, 11, 4, 2, 10, 5, 10),
    (11, 5, 4, 10, 12, 11, 8),
    (11, 1, 11, 2, 8, 1, 1),
    (6, 9, 12, 8, 4, 4, 7),
    (7, 6, 9, 9, 5, 1, 4),
    (2, 12, 12, 5, 11, 2, 12),
    (2, 5, 2, 1, 3, 9, 8),
    (4, 3, 8, 3, 11, 9, 2),
    (10, 12, 2, 3, 3, 6, 2),
    (8, 4, 9, 2, 2, 6, 4),
)


def r7_reduce(raw):
    return tuple((raw[index] - raw[Q - 1]) % P for index in range(Q - 1))


def r7_mul(left, right):
    coefficients = [0] * Q
    for i, x in enumerate(left):
        for j, y in enumerate(right):
            coefficients[(i + j) % Q] += x * y
    top = coefficients[Q - 1]
    return tuple((coefficients[index] - top) % P for index in range(Q - 1))


R7_ONE = (1, 0, 0, 0, 0, 0)


def r7_inverse(value):
    columns = [
        r7_mul(value, tuple(int(index == power) for index in range(Q - 1)))
        for power in range(Q - 1)
    ]
    matrix = [
        [columns[column][row] for column in range(Q - 1)]
        + [R7_ONE[row]]
        for row in range(Q - 1)
    ]
    for column in range(Q - 1):
        pivot = next(
            (row for row in range(column, Q - 1) if matrix[row][column] % P),
            None,
        )
        require(pivot is not None, "attempted to invert an R7 nonunit")
        if pivot != column:
            matrix[column], matrix[pivot] = matrix[pivot], matrix[column]
        inverse = pow(matrix[column][column] % P, -1, P)
        matrix[column] = [entry * inverse % P for entry in matrix[column]]
        for row in range(Q - 1):
            if row == column:
                continue
            factor = matrix[row][column] % P
            matrix[row] = [
                (matrix[row][j] - factor * matrix[column][j]) % P
                for j in range(Q)
            ]
    answer = tuple(matrix[row][-1] for row in range(Q - 1))
    require(r7_mul(value, answer) == R7_ONE, "bad R7 inverse")
    return answer


def quadratic_remainder(value, linear_coefficient):
    work = list(value)
    for degree in range(len(work) - 1, 1, -1):
        coefficient = work[degree] % P
        work[degree] = 0
        work[degree - 1] -= coefficient * linear_coefficient
        work[degree - 2] -= coefficient
    return work[0] % P, work[1] % P


def polynomial_product(left, right):
    out = [0] * (len(left) + len(right) - 1)
    for i, x in enumerate(left):
        for j, y in enumerate(right):
            out[i + j] = (out[i + j] + x * y) % P
    return tuple(out)


def charged_atlas_probe():
    factors = ((1, 3, 1), (1, 5, 1), (1, 6, 1))
    factor_product = (1,)
    for factor in factors:
        factor_product = polynomial_product(factor_product, factor)
    require(factor_product == (1,) * Q, "Phi7 factorization changed")

    rows = tuple(r7_reduce(row) for row in Y_RAW)
    require(len(set(rows)) == P, "charged atlas rows collided")
    inverses = tuple(r7_inverse(row) for row in rows)
    contrasts = []
    nonunits = []
    zero_component_histogram = Counter()
    for q in range(P):
        for step in range(1, P):
            endpoint = (q + step) % P
            transport = r7_mul(rows[endpoint], inverses[q])
            require(r7_mul(transport, rows[q]) == rows[endpoint],
                    "charged transport changed")
            normalized = r7_mul(inverses[endpoint], r7_mul(transport, rows[q]))
            require(normalized == R7_ONE, "vertex gauge did not trivialize edge")
            difference = tuple(
                (R7_ONE[index] - transport[index]) % P
                for index in range(Q - 1)
            )
            require(any(difference), "nontrivial atlas edge became identity")
            require(any(r7_mul(difference, difference)),
                    "atlas contrast became square-zero")
            zero_components = tuple(
                linear
                for linear in (3, 5, 6)
                if quadratic_remainder(difference, linear) == (0, 0)
            )
            zero_component_histogram[zero_components] += 1
            if zero_components:
                nonunits.append((q, step, endpoint))
            contrasts.append((difference, not zero_components))
    expected_nonunits = (
        (9, 5, 1),
        (12, 5, 4),
        (1, 8, 9),
        (4, 8, 12),
    )
    require(
        tuple(sorted(nonunits)) == tuple(sorted(expected_nonunits)),
        "nonunit contrasts changed",
    )
    require(
        zero_component_histogram == Counter({(): 152, (6,): 4}),
        "CRT contrast census changed",
    )

    nonzero_products = 0
    unit_products = 0
    for (left, left_unit), (right, right_unit) in product(contrasts, repeat=2):
        value = r7_mul(left, right)
        require(any(value), "two charged contrasts annihilated")
        nonzero_products += 1
        unit_products += int(left_unit and right_unit)
    require(nonzero_products == 156**2, "mixed product census changed")
    require(unit_products == 152**2, "unit product census changed")
    return (
        len(contrasts),
        sum(is_unit for _value, is_unit in contrasts),
        expected_nonunits,
        nonzero_products,
        unit_products,
        dict(zero_component_histogram),
    )


def rees_heisenberg_and_address_probe():
    central = (P**2, P, P, 1)
    valuations = (2, 1, 1, 0)
    leading_units = tuple(
        central[index] // P ** valuations[index] % P
        for index in range(4)
    )
    valuation_curvature = valuations[0] - valuations[1] - valuations[2] + valuations[3]
    common_grade = (central[0], P * central[1], P * central[2], P**2 * central[3])
    require(
        mobius(central) == 144
        and leading_units == (1, 1, 1, 1)
        and valuation_curvature == 0
        and len(set(common_grade)) == 1,
        "THM-2806 Rees profile changed",
    )

    weyl_checks = 0
    for root in range(P):
        tm = ((root + 1) % P, pow(ZETA, -root, FIELD))
        mt = ((root + 1) % P, pow(ZETA, -(root + 1), FIELD))
        require(tm[0] == mt[0] and tm[1] == ZETA * mt[1] % FIELD,
                "Weyl relation changed")
        weyl_checks += 1
    cocycle_square = (1, 1, 1, ZETA)
    require(
        cross_ratio(cocycle_square, FIELD) == ZETA
        and mobius(cocycle_square, FIELD) == ZETA - 1,
        "Heisenberg square changed",
    )

    determinant_counts = Counter()
    swap_checks = 0
    for sx, sy, tx, ty in product(range(P), repeat=4):
        determinant = (sx * ty - sy * tx) % P
        determinant_counts[determinant] += 1
        if determinant:
            phase = pow(ZETA, determinant, FIELD)
            reverse = pow(ZETA, -determinant, FIELD)
            require(phase * reverse % FIELD == 1, "orientation inversion changed")
            swap_checks += 1
    require(
        determinant_counts[0] == 2353
        and all(determinant_counts[value] == 2184 for value in range(1, P)),
        "symplectic frame census changed",
    )
    require(swap_checks == 26208, "nondegenerate frame census changed")

    pure = 13
    vertical = 169 * 4079
    diagonal = 689364
    require(pure + vertical == diagonal and 4079 % P == 10,
            "THM-2807 triangle changed")

    address_shifts = (-13, 1248)
    residues = tuple((shift // P) % P for shift in address_shifts)
    epsilon = tuple(
        (
            (comb(12, degree) if degree <= 12 else 0)
            + (comb(5, degree) if degree <= 5 else 0)
        ) % P
        for degree in range(P)
    )
    weighted = tuple(9 * coefficient % P for coefficient in epsilon)
    require(residues == (12, 5) and epsilon[:2] == (2, 4)
            and weighted[:2] == (5, 10), "cofiber phase germ changed")
    return (
        central,
        valuations,
        leading_units,
        valuation_curvature,
        common_grade,
        weyl_checks,
        cocycle_square,
        tuple(sorted(determinant_counts.items())),
        swap_checks,
        (pure, vertical, diagonal),
        address_shifts,
        residues,
        epsilon,
        weighted,
    )


def main():
    idempotent = idempotent_and_provenance_probe()
    projective = projective_square_probe()
    character = character_and_quarter_turn_probe()
    atlas = charged_atlas_probe()
    geometry = rees_heisenberg_and_address_probe()

    print("THM-2814 PROJECTIVE ALLOCATION HOLONOMY EXACT COMPANION")
    print("status=VERIFIED-EXACT; no row exclusion or LRC14 conclusion")
    print(
        f"branch_A=idempotents={idempotent[0]};"
        f"normalized_nonzero_Segre_contrasts={idempotent[1]};"
        f"linear_provenance_checks={idempotent[2]}"
    )
    print(
        f"coarse_hostile={idempotent[3]};mu={mobius(idempotent[3])};"
        "source=joint-absent sector"
    )
    print(
        f"branch_B=cross_ratio_fibres={projective[0]};"
        f"gauge_checks={projective[1]}"
    )
    print(
        f"kappa_one_chart_hostile={projective[2]};"
        f"mu={mobius(projective[2], P)};normalized={projective[3]}"
    )
    print(
        f"character_faces_by_origin={character[0]};"
        f"representative_factors={character[1]};"
        f"normalized={character[2]}"
    )
    print(
        f"quarter_turn=mixed_vector={character[3]};"
        f"scalar_readout={character[4]}"
    )
    print(
        f"thm2593_atlas=edges={atlas[0]};unit_contrasts={atlas[1]};"
        f"nonunits={atlas[2]};nonzero_ordered_products={atlas[3]};"
        f"unit_products={atlas[4]};crt={atlas[5]}"
    )
    print(
        f"thm2806_rees=central={geometry[0]};valuations={geometry[1]};"
        f"leading_units={geometry[2]};curvature={geometry[3]};"
        f"common_grade={geometry[4]}"
    )
    print(
        f"heisenberg=weyl_checks={geometry[5]};"
        f"cocycle_square={geometry[6]};kappa={ZETA}"
    )
    print(
        f"symplectic_frames=determinant_counts={geometry[7]};"
        f"oriented_nondegenerate={geometry[8]};swap_inverts_holonomy=yes"
    )
    print(
        f"thm2807=edges={geometry[9]};holonomy=1;vertical_residue=10"
    )
    print(
        f"thm2771_copy_germ=shifts={geometry[10]};residues={geometry[11]};"
        f"epsilon={geometry[12]};weighted={geometry[13]}"
    )
    print(
        "boundary=Branch A needs independently normalized non-idempotent "
        "toggles on one raw physical atom. Branch B needs nontrivial square "
        "holonomy when only vertex lines are supplied. Existing fixed, "
        "charged-atlas, address-simplex, and cofiber objects pay neither "
        "physical invoice; THM-2779 pays Branch B only abstractly."
    )
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
