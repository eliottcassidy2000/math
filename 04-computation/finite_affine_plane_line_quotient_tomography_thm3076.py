#!/usr/bin/env python3
"""Exact controls for finite affine-plane line-quotient tomography.

The proof is symbolic.  This companion exhausts the prime fields 2,3,5,7,
checks every projector and every labelled view subbank over Q, constructs the
sharp missing-direction hostiles, verifies both the good-characteristic and
characteristic-p boundaries, and identifies generatorwise the common A4
permutation model arising from F2^2 semidirect C3 and PSL2(F3).
"""

from __future__ import annotations

from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations


PRIMES = (2, 3, 5, 7)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def add_vectors(x, y, p):
    return ((x[0] + y[0]) % p, (x[1] + y[1]) % p)


def scale_vector(a, x, p):
    return (a * x[0] % p, a * x[1] % p)


def points(p):
    return tuple((a, b) for a in range(p) for b in range(p))


def lines(p):
    generators = tuple((1, slope) for slope in range(p)) + ((0, 1),)
    return tuple(
        tuple(sorted({scale_vector(a, generator, p) for a in range(p)}))
        for generator in generators
    )


def zero_matrix(rows, columns):
    return [[F(0) for _ in range(columns)] for _ in range(rows)]


def identity(size):
    result = zero_matrix(size, size)
    for index in range(size):
        result[index][index] = F(1)
    return result


def matrix_add(left, right):
    require(len(left) == len(right), "row mismatch")
    require(not left or len(left[0]) == len(right[0]), "column mismatch")
    return [
        [a + b for a, b in zip(left_row, right_row)]
        for left_row, right_row in zip(left, right)
    ]


def matrix_scale(value, matrix):
    return [[value * entry for entry in row] for row in matrix]


def matrix_multiply(left, right):
    require(not left or len(left[0]) == len(right), "product shape")
    columns = len(right[0]) if right else 0
    return [
        [
            sum(left[i][k] * right[k][j] for k in range(len(right)))
            for j in range(columns)
        ]
        for i in range(len(left))
    ]


def matrix_vector(matrix, vector):
    return [sum(entry * value for entry, value in zip(row, vector)) for row in matrix]


def matrix_equal(left, right):
    return left == right


def rank(matrix):
    rows = [list(row) for row in matrix]
    if not rows:
        return 0
    row_count = len(rows)
    column_count = len(rows[0])
    pivot_row = 0
    for column in range(column_count):
        pivot = next(
            (row for row in range(pivot_row, row_count) if rows[row][column]),
            None,
        )
        if pivot is None:
            continue
        rows[pivot_row], rows[pivot] = rows[pivot], rows[pivot_row]
        pivot_value = rows[pivot_row][column]
        rows[pivot_row] = [value / pivot_value for value in rows[pivot_row]]
        for row in range(row_count):
            if row == pivot_row or not rows[row][column]:
                continue
            factor = rows[row][column]
            rows[row] = [
                value - factor * pivot_entry
                for value, pivot_entry in zip(rows[row], rows[pivot_row])
            ]
        pivot_row += 1
        if pivot_row == row_count:
            break
    return pivot_row


def rank_mod_prime(matrix, p):
    rows = [
        [
            int(entry.numerator) * pow(int(entry.denominator), -1, p) % p
            for entry in row
        ]
        for row in matrix
    ]
    if not rows:
        return 0
    row_count = len(rows)
    column_count = len(rows[0])
    pivot_row = 0
    for column in range(column_count):
        pivot = next(
            (row for row in range(pivot_row, row_count) if rows[row][column]),
            None,
        )
        if pivot is None:
            continue
        rows[pivot_row], rows[pivot] = rows[pivot], rows[pivot_row]
        inverse = pow(rows[pivot_row][column], -1, p)
        rows[pivot_row] = [value * inverse % p for value in rows[pivot_row]]
        for row in range(row_count):
            if row == pivot_row or not rows[row][column]:
                continue
            factor = rows[row][column]
            rows[row] = [
                (value - factor * pivot_entry) % p
                for value, pivot_entry in zip(rows[row], rows[pivot_row])
            ]
        pivot_row += 1
        if pivot_row == row_count:
            break
    return pivot_row


def matrix_equal_mod_prime(left, right, p):
    require(len(left) == len(right), "mod-p row mismatch")
    require(not left or len(left[0]) == len(right[0]), "mod-p column mismatch")
    for left_row, right_row in zip(left, right):
        for a, b in zip(left_row, right_row):
            difference = a - b
            require(difference.denominator % p, (p, "noninvertible denominator"))
            if difference.numerator % p:
                return False
    return True


def flatten_matrix(matrix):
    return [entry for row in matrix for entry in row]


def translation_matrix(p, shift):
    universe = points(p)
    position = {point: index for index, point in enumerate(universe)}
    matrix = zero_matrix(len(universe), len(universe))
    for row, point in enumerate(universe):
        matrix[row][position[add_vectors(point, shift, p)]] = F(1)
    return matrix


def line_sum_matrix(p, line):
    result = zero_matrix(p * p, p * p)
    for shift in line:
        result = matrix_add(result, translation_matrix(p, shift))
    return result


def quotient_label(p, line_index, point):
    if line_index < p:
        return (point[1] - line_index * point[0]) % p
    return point[0]


def coset_hostile(p, missing_line_index):
    # Difference of two quotient fibres: mean zero, fixed by the missing-line
    # average, and killed by every transverse line average.
    return tuple(
        F(int(quotient_label(p, missing_line_index, point) == 0))
        - F(int(quotient_label(p, missing_line_index, point) == 1))
        for point in points(p)
    )


def finite_plane_record(p):
    universe = points(p)
    direction_lines = lines(p)
    require(len(direction_lines) == p + 1, (p, "line count"))
    require(all(len(line) == p for line in direction_lines), (p, "line size"))
    nonzero_memberships = {
        point: sum(point in line for line in direction_lines)
        for point in universe
        if point != (0, 0)
    }
    require(set(nonzero_memberships.values()) == {1}, (p, "unique direction"))

    line_sums = tuple(line_sum_matrix(p, line) for line in direction_lines)
    projectors = tuple(matrix_scale(F(1, p), matrix) for matrix in line_sums)
    all_sum = zero_matrix(p * p, p * p)
    for shift in universe:
        all_sum = matrix_add(all_sum, translation_matrix(p, shift))
    global_projector = matrix_scale(F(1, p * p), all_sum)
    ident = identity(p * p)

    for index, projector in enumerate(projectors):
        require(matrix_equal(matrix_multiply(projector, projector), projector), (p, index, "idempotent"))
        require(rank(projector) == p, (p, index, "line rank"))

    centered_projectors = tuple(
        matrix_add(projector, matrix_scale(F(-1), global_projector))
        for projector in projectors
    )
    centered_sum = zero_matrix(p * p, p * p)
    for index, centered in enumerate(centered_projectors):
        centered_sum = matrix_add(centered_sum, centered)
        require(matrix_multiply(centered, centered) == centered, (p, index, "centered idempotent"))
        require(rank(centered) == p - 1, (p, index, "centered rank"))
    for i, j in combinations(range(p + 1), 2):
        require(
            matrix_multiply(centered_projectors[i], centered_projectors[j])
            == zero_matrix(p * p, p * p),
            (p, i, j, "centered annihilation"),
        )
    for i, j in combinations(range(p + 1), 2):
        require(
            matrix_equal(matrix_multiply(projectors[i], projectors[j]), global_projector),
            (p, i, j, "pair product"),
        )
        joint_labels = {
            (quotient_label(p, i, point), quotient_label(p, j, point))
            for point in universe
        }
        require(len(joint_labels) == p * p, (p, i, j, "joint labels"))

    projector_sum = zero_matrix(p * p, p * p)
    for projector in projectors:
        projector_sum = matrix_add(projector_sum, projector)
    require(
        matrix_equal(projector_sum, matrix_add(ident, matrix_scale(F(p), global_projector))),
        (p, "reconstruction identity"),
    )
    require(
        centered_sum == matrix_add(ident, matrix_scale(F(-1), global_projector)),
        (p, "centered resolution"),
    )
    operator_span_rank = rank(
        [flatten_matrix(projector) for projector in projectors]
        + [flatten_matrix(global_projector)]
    )
    require(operator_span_rank == p + 2, (p, "operator-span rank", operator_span_rank))

    rank_census = []
    for size in range(1, p + 2):
        ranks = set()
        for chosen in combinations(range(p + 1), size):
            observation = []
            for index in chosen:
                observation.extend(projectors[index])
            observed_rank = rank(observation)
            require(observed_rank == 1 + size * (p - 1), (p, chosen, observed_rank))
            observable_projector = matrix_scale(F(1 - size), global_projector)
            for index in chosen:
                observable_projector = matrix_add(observable_projector, projectors[index])
            require(
                matrix_multiply(observable_projector, observable_projector)
                == observable_projector,
                (p, chosen, "observable projector"),
            )
            require(rank(observable_projector) == observed_rank, (p, chosen, "observable rank"))
            ranks.add(observed_rank)
        require(len(ranks) == 1, (p, size, ranks))
        rank_census.append(next(iter(ranks)))

    hostile_checks = 0
    for missing in range(p + 1):
        hostile = coset_hostile(p, missing)
        require(any(hostile), (p, missing, "zero hostile"))
        require(matrix_vector(projectors[missing], hostile) == list(hostile), (p, missing, "missing view"))
        require(matrix_vector(global_projector, hostile) == [F(0)] * (p * p), (p, missing, "hostile mean"))
        for observed in range(p + 1):
            if observed == missing:
                continue
            require(
                matrix_vector(projectors[observed], hostile) == [F(0)] * (p * p),
                (p, missing, observed, "hostile leakage"),
            )
            hostile_checks += 1

    # Good-characteristic control at a prime ell dividing p+1.  Nothing fails:
    # p remains invertible, all line averages retain their ranks, and the full
    # bank still reconstructs.  Only the *raw sum* kills constants because
    # p+1=0 in the coefficient field.  This rules out the tempting but false
    # extra exclusion char(k) | p+1.
    good_prime = {2: 3, 3: 2, 5: 2, 7: 2}[p]
    require((p + 1) % good_prime == 0 and p % good_prime, (p, good_prime, "good-prime choice"))
    for index, projector in enumerate(projectors):
        require(rank_mod_prime(projector, good_prime) == p, (p, good_prime, index, "good rank"))
        require(
            matrix_equal_mod_prime(matrix_multiply(projector, projector), projector, good_prime),
            (p, good_prime, index, "good idempotent"),
        )
    require(rank_mod_prime([row for projector in projectors for row in projector], good_prime) == p * p,
            (p, good_prime, "good full-bank rank"))
    require(
        matrix_equal_mod_prime(
            projector_sum,
            matrix_add(ident, matrix_scale(F(p), global_projector)),
            good_prime,
        ),
        (p, good_prime, "good reconstruction"),
    )
    constant = [F(1)] * (p * p)
    raw_constant = matrix_vector(projector_sum, constant)
    require(
        all(value.numerator % good_prime == 0 for value in raw_constant),
        (p, good_prime, "raw sum does not kill constants"),
    )

    # Integral identity and the exact characteristic-p failure: each line norm
    # squares to p times itself, hence becomes square-zero modulo p.  Distinct
    # norms multiply to the global norm, and their sum is the same global norm
    # modulo p, but they no longer define normalized Reynolds projectors.
    integral_sum = zero_matrix(p * p, p * p)
    for matrix in line_sums:
        integral_sum = matrix_add(integral_sum, matrix)
        square = matrix_multiply(matrix, matrix)
        require(square == matrix_scale(F(p), matrix), (p, "integral square"))
        require(
            all(int(entry) % p == 0 for row in square for entry in row),
            (p, "mod-p square-zero"),
        )
    for i, j in combinations(range(p + 1), 2):
        require(
            matrix_multiply(line_sums[i], line_sums[j]) == all_sum,
            (p, i, j, "integral transverse product"),
        )
    require(
        integral_sum == matrix_add(matrix_scale(F(p), ident), all_sum),
        (p, "integral reconstruction"),
    )
    require(
        matrix_equal_mod_prime(integral_sum, all_sum, p),
        (p, "mod-p norm sum"),
    )
    cleared_observation = []
    for matrix in line_sums:
        cleared_observation.extend(matrix)
    cleared_rank = rank_mod_prime(cleared_observation, p)
    require(cleared_rank == p * (p + 1) // 2, (p, "cleared characteristic-p rank", cleared_rank))

    return (
        p,
        len(universe),
        len(direction_lines),
        tuple(rank_census),
        hostile_checks,
        cleared_rank,
        good_prime,
        operator_span_rank,
    )


def permutation_compose(left, right):
    return tuple(left[right[index]] for index in range(len(left)))


def permutation_power(permutation, exponent):
    result = tuple(range(len(permutation)))
    for _ in range(exponent):
        result = permutation_compose(permutation, result)
    return result


def permutation_inverse(permutation):
    result = [None] * len(permutation)
    for index, image in enumerate(permutation):
        result[image] = index
    return tuple(result)


def permutation_conjugate(conjugator, permutation):
    return permutation_compose(
        conjugator,
        permutation_compose(permutation, permutation_inverse(conjugator)),
    )


def permutation_sign(permutation):
    inversions = sum(
        permutation[i] > permutation[j]
        for i in range(len(permutation))
        for j in range(i + 1, len(permutation))
    )
    return -1 if inversions % 2 else 1


def generated_group(generators):
    identity_permutation = tuple(range(len(generators[0])))
    group = {identity_permutation}
    frontier = [identity_permutation]
    while frontier:
        element = frontier.pop()
        for generator in generators:
            product = permutation_compose(generator, element)
            if product not in group:
                group.add(product)
                frontier.append(product)
    return frozenset(group)


def linear_image(matrix, vector, p):
    return (
        (matrix[0][0] * vector[0] + matrix[0][1] * vector[1]) % p,
        (matrix[1][0] * vector[0] + matrix[1][1] * vector[1]) % p,
    )


def affine_a4_over_f2():
    p = 2
    universe = points(p)
    position = {point: index for index, point in enumerate(universe)}
    rho = ((0, 1), (1, 1))
    translations = []
    for shift in universe:
        translations.append(
            tuple(position[add_vectors(point, shift, p)] for point in universe)
        )
    rho_permutation = tuple(position[linear_image(rho, point, p)] for point in universe)
    group = generated_group((translations[1], translations[2], rho_permutation))
    require(len(group) == 12, "affine A4 size")
    require(all(permutation_sign(element) == 1 for element in group), "affine A4 parity")
    require(frozenset(translations).issubset(group), "V4 translations absent")
    require(permutation_power(rho_permutation, 3) == tuple(range(4)), "rho order")
    line_sets = lines(2)
    direction_permutation = tuple(
        line_sets.index(tuple(sorted(linear_image(rho, point, 2) for point in line)))
        for line in line_sets
    )
    require(len(set(direction_permutation)) == 3, "rho direction cycle")
    require(permutation_power(direction_permutation, 3) == tuple(range(3)), "direction order")
    translation_x = tuple(position[add_vectors(point, (1, 0), p)] for point in universe)
    return group, direction_permutation, translation_x, rho_permutation


def determinant_two(matrix, p):
    return (matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]) % p


def matrix_two_multiply(left, right, p):
    return tuple(
        tuple(
            sum(left[i][k] * right[k][j] for k in range(2)) % p
            for j in range(2)
        )
        for i in range(2)
    )


def projective_permutation(matrix, p):
    line_sets = lines(p)
    return tuple(
        line_sets.index(tuple(sorted(linear_image(matrix, point, p) for point in line)))
        for line in line_sets
    )


def projective_a4_over_f3():
    p = 3
    matrices = tuple(
        ((a, b), (c, d))
        for a in range(p)
        for b in range(p)
        for c in range(p)
        for d in range(p)
        if determinant_two(((a, b), (c, d)), p) == 1
    )
    require(len(matrices) == 24, "SL2(F3) size")
    projective_group = frozenset(projective_permutation(matrix, p) for matrix in matrices)
    require(len(projective_group) == 12, "PSL2(F3) size")
    require(all(permutation_sign(element) == 1 for element in projective_group), "projective A4 parity")

    modular_s = ((0, -1), (1, 0))
    modular_t = ((1, 1), (0, 1))
    modular_r = matrix_two_multiply(modular_s, modular_t, p)
    s_permutation = projective_permutation(modular_s, p)
    r_permutation = projective_permutation(modular_r, p)
    require(permutation_power(s_permutation, 2) == tuple(range(4)), "modular S order")
    require(permutation_power(r_permutation, 3) == tuple(range(4)), "modular ST order")
    require(
        permutation_power(permutation_compose(s_permutation, r_permutation), 3)
        == tuple(range(4)),
        "modular S(ST) order",
    )
    require(generated_group((s_permutation, r_permutation)) == projective_group, "modular image")
    return projective_group, s_permutation, r_permutation


records = tuple(finite_plane_record(p) for p in PRIMES)
affine_a4, direction_cycle, affine_translation, affine_rho = affine_a4_over_f2()
projective_a4, modular_s, modular_r = projective_a4_over_f3()
require(affine_a4 == projective_a4, "two A4 permutation models differ")
conjugator = (1, 0, 2, 3)
require(
    permutation_conjugate(conjugator, modular_s) == affine_translation,
    "modular S does not conjugate to the affine translation",
)
require(
    permutation_conjugate(conjugator, modular_r) == affine_rho,
    "modular ST does not conjugate to the affine order-three map",
)

semantic = repr(
    (
        records,
        tuple(sorted(affine_a4)),
        direction_cycle,
        modular_s,
        modular_r,
        conjugator,
    )
).encode()
semantic_sha256 = sha256(semantic).hexdigest()

print("THM3076_FINITE_AFFINE_PLANE_LINE_QUOTIENT_TOMOGRAPHY")
print("scope=prime_fields:2,3,5,7;coefficient_field=Q;symbolic_proof=all_prime_p_with_p_invertible")
for p, point_count, line_count, ranks, hostile_checks, cleared_rank, good_prime, operator_span_rank in records:
    print(
        f"FIELD;p={p};points={point_count};directions={line_count};"
        f"subbank_ranks={','.join(map(str, ranks))};hostile_checks={hostile_checks};"
        f"char_p_cleared_rank={cleared_rank};good_char={good_prime};"
        f"operator_span_rank={operator_span_rank}"
    )
print("operator_law=sum_L_P_L=I+pP_V;distinct_P_L_P_M=P_V")
print("centered_law=E_L=P_L-P_V;E_L^2=E_L;E_LE_M=0;sum_E_L=I-P_V")
print("rank_law=rank(s_labelled_views)=1+s(p-1);kernel=(p+1-s)(p-1)")
print("label_table_boundary=two_joint_quotient_labels_are_bijective;two_separate_averages_have_kernel_(p-1)^2")
print("good_char_boundary=char_dividing_(p+1)_is_allowed;raw_projector_sum_kills_constants_but_corrected_formula_reconstructs")
print("char_p_boundary=A_L^2=pA_L;A_LA_M=A_V;sum_A_L=A_V_mod_p;full_cleared_rank=p(p+1)/2")
print("small_prime_bridge=F2^2_semidirect_C3=A4=PSL2(F3);generatorwise_S_to_translation,ST_to_rho;orders=2,3,3")
print(f"p2_direction_cycle={direction_cycle};A4_permutations={len(affine_a4)}")
print("scope_boundary=no_binary_or_ternary_tree_intertwiner;no_quartic,Farey,Keller,LRC_realization")
print(f"semantic_sha256={semantic_sha256}")
print("all_exact_controls=PASS")
