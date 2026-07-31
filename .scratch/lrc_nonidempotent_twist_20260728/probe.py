#!/usr/bin/env python3
"""Exact hostile probes for non-idempotent same-atom allocation amplitudes.

This is a dependency-free scratch companion.  It tests:

* the universal idempotent common-atom no-go;
* the support provenance of every linear coarse-graining;
* allocation-sign twisting versus genuine mixed interaction;
* a C13 local-system phase and its representative-gauge obstruction;
* the THM-2806 Rees/leading-unit square;
* the minimal Heisenberg cocycle escape;
* the THM-2807 translation-triangle holonomy;
* the two equal THM-2771 cofiber-copy address phases.
"""

from itertools import product
from math import comb


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


P = 13
FIELD = 53


def primitive_thirteenth_root():
    for candidate in range(2, FIELD):
        if (
            pow(candidate, P, FIELD) == 1
            and candidate != 1
        ):
            return candidate
    raise RuntimeError("no thirteenth root found")


ZETA = primitive_thirteenth_root()
require(pow(ZETA, P, FIELD) == 1, "root order does not divide 13")
require(
    all(pow(ZETA, exponent, FIELD) != 1 for exponent in range(1, P)),
    "root does not have exact order 13",
)


def mobius(vector, modulus=None):
    value = vector[0] - vector[1] - vector[2] + vector[3]
    return value if modulus is None else value % modulus


def allocation_spectrum(vector, modulus):
    """C2^2 Fourier coefficients in order 00,10,01,11."""
    out = []
    for source_character, target_character in product(range(2), repeat=2):
        value = 0
        for source_flag, target_flag in product(range(2), repeat=2):
            index = 2 * source_flag + target_flag
            sign = (
                -1
                if (
                    source_character * source_flag
                    + target_character * target_flag
                ) % 2
                else 1
            )
            value += sign * vector[index]
        out.append(value % modulus)
    return tuple(out)


def exhaustive_idempotent_and_coarse_graining():
    idempotents = tuple(
        value for value in range(P) if value * value % P == value
    )
    require(idempotents == (0, 1), "field idempotents changed")
    common_idempotent_faces = []
    for left, right in product(idempotents, repeat=2):
        vector = (1, left, right, left * right % P)
        if all(vector):
            common_idempotent_faces.append(mobius(vector, P))
    require(
        common_idempotent_faces == [0],
        "idempotent fourfold co-support acquired curvature",
    )

    successful_units = []
    for left, right in product(range(1, P), repeat=2):
        vector = (1, left, right, left * right % P)
        if mobius(vector, P):
            successful_units.append((left, right))
    require(
        len(successful_units) == (P - 2) ** 2,
        "non-idempotent unit census changed",
    )

    checks = 0
    nonzero_examples = []
    for left_mask in range(1 << 4):
        left = tuple((left_mask >> index) & 1 for index in range(4))
        for right_mask in range(1 << 4):
            right = tuple((right_mask >> index) & 1 for index in range(4))
            for weights in product((-1, 0, 1), repeat=4):
                bare = sum(weights)
                source = sum(
                    weights[index] * left[index] for index in range(4)
                )
                target = sum(
                    weights[index] * right[index] for index in range(4)
                )
                both = sum(
                    weights[index] * left[index] * right[index]
                    for index in range(4)
                )
                face = bare - source - target + both
                complement = sum(
                    weights[index]
                    * (1 - left[index])
                    * (1 - right[index])
                    for index in range(4)
                )
                common_part = sum(
                    weights[index]
                    * left[index]
                    * right[index]
                    * (1 - left[index])
                    * (1 - right[index])
                    for index in range(4)
                )
                require(face == complement, "coarse support identity failed")
                require(
                    common_part == 0,
                    "common-support atom contributed to Boolean face",
                )
                if face and len(nonzero_examples) < 3:
                    nonzero_examples.append(
                        (left, right, weights, face)
                    )
                checks += 1
    require(checks == 16 * 16 * 3**4, "exhaustive universe changed")

    quotient_example = (4, 2, 2, 1)
    require(
        mobius(quotient_example) == 1,
        "conditional-expectation example changed",
    )
    return (
        idempotents,
        len(successful_units),
        checks,
        tuple(nonzero_examples),
        quotient_example,
    )


def sign_and_character_probes():
    flat = (1, 1, 1, 1)
    alternating_twist = (1, -1, -1, 1)
    flat_spectrum = allocation_spectrum(flat, P)
    twisted_spectrum = allocation_spectrum(alternating_twist, P)
    require(flat_spectrum == (4, 0, 0, 0), "flat spectrum changed")
    require(
        twisted_spectrum == (0, 0, 0, 4),
        "allocation sign did not shift trivial mass to mixed character",
    )

    phases = []
    for origin in range(P):
        alpha = pow(ZETA, origin, FIELD)
        beta = pow(ZETA, 2 * origin, FIELD)
        vector = (
            1,
            alpha,
            beta,
            alpha * beta % FIELD,
        )
        phase = mobius(vector, FIELD)
        require(
            phase
            == (1 - alpha) * (1 - beta) % FIELD,
            "rank-one local-system face stopped factoring",
        )
        phases.append(phase)
    require(phases[0] == 0, "normalized marked origin acquired phase")
    require(
        all(phases[index] for index in range(1, P)),
        "nontrivial character origin lost its abstract face",
    )

    source_gauge = ZETA
    target_gauge = pow(ZETA, -1, FIELD)
    vertex_gauge = (
        1,
        source_gauge,
        target_gauge,
        source_gauge * target_gauge % FIELD,
    )
    require(
        len(set(vertex_gauge)) > 1,
        "nontrivial allocation local system unexpectedly scalar-descended",
    )

    # The canonical allocation-graded parallel transport divides out the
    # edge characters.  It returns the square to the flat vector.
    origin = 4
    alpha = pow(ZETA, origin, FIELD)
    beta = pow(ZETA, 2 * origin, FIELD)
    raw = (1, alpha, beta, alpha * beta % FIELD)
    transported = (
        raw[0],
        raw[1] * pow(alpha, -1, FIELD) % FIELD,
        raw[2] * pow(beta, -1, FIELD) % FIELD,
        raw[3] * pow(alpha * beta % FIELD, -1, FIELD) % FIELD,
    )
    require(transported == flat, "parallel transport did not flatten square")
    require(mobius(transported, FIELD) == 0, "covariant face became nonzero")
    return (
        flat_spectrum,
        twisted_spectrum,
        tuple(phases),
        vertex_gauge,
        raw,
        transported,
    )


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


def r7_reduce(row):
    return tuple((row[index] - row[6]) % P for index in range(6))


def r7_mul(left, right):
    coefficients = [0] * 7
    for i, x in enumerate(left):
        for j, y in enumerate(right):
            coefficients[(i + j) % 7] += x * y
    top = coefficients[6]
    return tuple((coefficients[index] - top) % P for index in range(6))


R7_ONE = (1, 0, 0, 0, 0, 0)


def r7_matrix(value):
    basis = tuple(
        tuple(int(index == column) for index in range(6))
        for column in range(6)
    )
    columns = tuple(r7_mul(value, vector) for vector in basis)
    return [
        [columns[column][row] for column in range(6)]
        for row in range(6)
    ]


def mod_determinant(matrix):
    work = [[value % P for value in row] for row in matrix]
    answer = 1
    for column in range(len(work)):
        pivot = next(
            (
                row for row in range(column, len(work))
                if work[row][column]
            ),
            None,
        )
        if pivot is None:
            return 0
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            answer = -answer
        diagonal = work[column][column] % P
        answer = answer * diagonal % P
        inverse = pow(diagonal, -1, P)
        for index in range(column, len(work)):
            work[column][index] = work[column][index] * inverse % P
        for row in range(column + 1, len(work)):
            factor = work[row][column]
            for index in range(column, len(work)):
                work[row][index] = (
                    work[row][index]
                    - factor * work[column][index]
                ) % P
    return answer % P


def r7_inverse(value):
    matrix = r7_matrix(value)
    augmented = [
        matrix[row][:]
        + [R7_ONE[row]]
        for row in range(6)
    ]
    for column in range(6):
        pivot = next(
            (
                row for row in range(column, 6)
                if augmented[row][column]
            ),
            None,
        )
        require(pivot is not None, "attempted to invert an R7 nonunit")
        if pivot != column:
            augmented[column], augmented[pivot] = (
                augmented[pivot], augmented[column]
            )
        inverse = pow(augmented[column][column], -1, P)
        for index in range(column, 7):
            augmented[column][index] = (
                augmented[column][index] * inverse % P
            )
        for row in range(6):
            if row == column:
                continue
            factor = augmented[row][column]
            for index in range(column, 7):
                augmented[row][index] = (
                    augmented[row][index]
                    - factor * augmented[column][index]
                ) % P
    answer = tuple(augmented[row][6] for row in range(6))
    require(r7_mul(value, answer) == R7_ONE, "R7 inverse failed")
    return answer


def quadratic_mul(left, right, trace_parameter):
    """Multiply modulo z^2+trace_parameter*z+1 over F_13."""
    return (
        (
            left[0] * right[0]
            - left[1] * right[1]
        ) % P,
        (
            left[0] * right[1]
            + left[1] * right[0]
            - trace_parameter * left[1] * right[1]
        ) % P,
    )


def quadratic_evaluate(polynomial, trace_parameter):
    answer = (0, 0)
    power = (1, 0)
    z = (0, 1)
    for coefficient in polynomial:
        answer = (
            (answer[0] + coefficient * power[0]) % P,
            (answer[1] + coefficient * power[1]) % P,
        )
        power = quadratic_mul(power, z, trace_parameter)
    return answer


def charged_atlas_probe():
    """Independent census of THM-2593 edge-unit mixed amplitudes."""
    rows = tuple(r7_reduce(row) for row in Y_RAW)
    require(len(set(rows)) == P, "charged atlas rows stopped being distinct")
    inverses = tuple(r7_inverse(row) for row in rows)
    differences = []
    exceptional = []
    component_zero_histogram = {}
    for step in range(1, P):
        for sheet in range(P):
            transport = r7_mul(rows[(sheet + step) % P], inverses[sheet])
            difference = tuple(
                (R7_ONE[index] - transport[index]) % P
                for index in range(6)
            )
            square = r7_mul(difference, difference)
            determinant = mod_determinant(r7_matrix(difference))
            zero_components = tuple(
                parameter
                for parameter in (3, 5, 6)
                if quadratic_evaluate(difference, parameter) == (0, 0)
            )
            component_zero_histogram[zero_components] = (
                component_zero_histogram.get(zero_components, 0) + 1
            )
            require(
                (determinant == 0) == bool(zero_components),
                "matrix and CRT unit tests disagree",
            )
            require(any(difference), "distinct atlas edge became identity")
            require(any(square), "atlas edge contrast became square-zero")
            differences.append((sheet, step, difference, determinant))
            if determinant == 0:
                exceptional.append((sheet, step, (sheet + step) % P))
    require(len(differences) == 156, "atlas edge universe changed")
    require(
        len(exceptional) == 4,
        "atlas edge-difference unit census changed",
    )
    require(
        component_zero_histogram == {(): 152, (6,): 4},
        "atlas CRT zero-component census changed",
    )

    product_nonzero = 0
    product_unit = 0
    for left in differences:
        for right in differences:
            mixed = r7_mul(left[2], right[2])
            require(any(mixed), "two charged atlas contrasts annihilated")
            product_nonzero += 1
            product_unit += bool(mod_determinant(r7_matrix(mixed)))
    require(product_nonzero == 156**2, "atlas pair universe changed")
    require(product_unit == 23104, "atlas pair unit census changed")
    require(
        product_unit == 152**2,
        "atlas product units stopped being exactly unit-by-unit pairs",
    )
    return (
        len(differences),
        len(differences) - len(exceptional),
        tuple(exceptional),
        product_nonzero,
        product_unit,
        component_zero_histogram,
    )


def quarter_turn_probe():
    """The minimal integral non-idempotent phase local system."""
    vector = (1, 0)

    def turn(value):
        return (-value[1], value[0])

    once = turn(vector)
    twice = turn(once)
    require(twice == (-1, 0), "quarter turn stopped squaring to minus one")
    mixed_vector = (
        vector[0] - once[0] - once[0] + twice[0],
        vector[1] - once[1] - once[1] + twice[1],
    )
    require(mixed_vector == (0, -2), "quarter-turn mixed face changed")
    scalar_readout = (
        sum(vector),
        sum(once),
        sum(once),
        sum(twice),
    )
    require(all(scalar_readout), "quarter-turn readout lost co-support")
    require(mobius(scalar_readout) == -2, "quarter-turn readout flattened")
    return vector, once, twice, mixed_vector, scalar_readout


def rees_probe():
    central = (P**2, P, P, 1)
    valuations = (2, 1, 1, 0)
    leading_units = tuple(
        central[index] // P ** valuations[index] % P
        for index in range(4)
    )
    valuation_curvature = (
        valuations[0] - valuations[1]
        - valuations[2] + valuations[3]
    )
    require(leading_units == (1, 1, 1, 1), "leading square changed")
    require(valuation_curvature == 0, "valuation square gained curvature")
    require(mobius(central) == 144, "central face changed")
    require(mobius(central) % P == 1, "central face residue changed")
    require(mobius(leading_units, P) == 0, "graded common square not flat")

    lifted_to_grade_two = (
        central[0],
        P * central[1],
        P * central[2],
        P**2 * central[3],
    )
    require(
        len(set(lifted_to_grade_two)) == 1,
        "canonical common-grade transport did not flatten square",
    )
    return (
        central,
        valuations,
        leading_units,
        valuation_curvature,
        lifted_to_grade_two,
    )


def projective_square_gauge_probe():
    """Classify nonzero F_13 squares modulo row/column rephasing."""
    units = tuple(range(1, P))
    normalized_histogram = {kappa: 0 for kappa in units}
    mobius_histogram_at_one = {value: 0 for value in range(P)}
    normalization_checks = 0

    for v00, v10, v01, v11 in product(units, repeat=4):
        square = (v00, v10, v01, v11)
        kappa = (
            v00
            * v11
            * pow(v10 * v01 % P, -1, P)
        ) % P

        # Fix the one-dimensional kernel of row/column rephasing by r0=1.
        # The remaining uniquely determined factors send the first three
        # vertices to one.
        r0 = 1
        r1 = v00 * pow(v10, -1, P) % P
        c0 = pow(v00, -1, P)
        c1 = pow(v01, -1, P)
        normalized = (
            r0 * c0 * v00 % P,
            r1 * c0 * v10 % P,
            r0 * c1 * v01 % P,
            r1 * c1 * v11 % P,
        )
        require(
            normalized == (1, 1, 1, kappa),
            "anchored row/column normalization failed",
        )
        normalized_histogram[kappa] += 1
        if kappa == 1:
            mobius_histogram_at_one[mobius(square, P)] += 1
        normalization_checks += 1

    require(
        normalization_checks == (P - 1) ** 4,
        "nonzero square universe changed",
    )
    require(
        set(normalized_histogram.values()) == {(P - 1) ** 3},
        "cross-ratio fibres stopped having uniform size",
    )

    orbit_parametrizations = {}
    for kappa in units:
        orbit = set()
        for r1, c0, c1 in product(units, repeat=3):
            orbit.add(
                (
                    c0,
                    c1,
                    r1 * c0 % P,
                    r1 * c1 * kappa % P,
                )
            )
        require(
            len(orbit) == (P - 1) ** 3,
            "anchored gauge parameters ceased to be injective",
        )
        require(
            all(
                (
                    square[0]
                    * square[3]
                    * pow(square[1] * square[2] % P, -1, P)
                ) % P
                == kappa
                for square in orbit
            ),
            "row/column orbit changed cross-ratio",
        )
        orbit_parametrizations[kappa] = len(orbit)

    flat = (1, 1, 1, 1)
    segre_rephasing = (1, 2, 3, 6)
    require(
        mobius(flat, P) == 0
        and mobius(segre_rephasing, P) == 2,
        "kappa-one additive hostile examples changed",
    )
    require(
        (
            flat[0] * flat[3]
            * pow(flat[1] * flat[2] % P, -1, P)
        ) % P
        == 1
        and (
            segre_rephasing[0] * segre_rephasing[3]
            * pow(
                segre_rephasing[1] * segre_rephasing[2] % P,
                -1,
                P,
            )
        ) % P
        == 1,
        "hostile examples left the kappa-one orbit",
    )
    require(
        sum(bool(count) for count in mobius_histogram_at_one.values()) > 1,
        "additive Mobius value became projectively intrinsic at kappa one",
    )
    return (
        normalization_checks,
        normalized_histogram,
        orbit_parametrizations,
        mobius_histogram_at_one,
        (flat, mobius(flat, P)),
        (segre_rephasing, mobius(segre_rephasing, P)),
    )


def symplectic_heisenberg_frame_probe():
    """Enumerate nondegenerate frames and their Weyl commutator phases."""
    vectors = tuple(product(range(P), repeat=2))
    nonzero = tuple(vector for vector in vectors if vector != (0, 0))
    phase_histogram = {}
    determinant_histogram = {value: 0 for value in range(1, P)}
    frames = 0
    swap_checks = 0

    for left in nonzero:
        for right in nonzero:
            determinant = (
                left[0] * right[1] - left[1] * right[0]
            ) % P
            if determinant == 0:
                continue
            phase = pow(ZETA, determinant, FIELD)
            reverse_determinant = (
                right[0] * left[1] - right[1] * left[0]
            ) % P
            reverse_phase = pow(ZETA, reverse_determinant, FIELD)
            require(
                reverse_phase == pow(phase, -1, FIELD),
                "frame swap stopped inverting Heisenberg holonomy",
            )
            determinant_histogram[determinant] += 1
            phase_histogram[phase] = phase_histogram.get(phase, 0) + 1
            frames += 1
            swap_checks += 1

    require(
        frames == (P**2 - 1) * (P**2 - P),
        "ordered nondegenerate frame census changed",
    )
    require(
        set(determinant_histogram.values()) == {2184},
        "nonzero determinant fibres stopped being uniform",
    )
    require(
        len(phase_histogram) == P - 1
        and set(phase_histogram.values()) == {2184},
        "nontrivial Heisenberg phase fibres stopped being uniform",
    )
    require(swap_checks == frames, "frame swap audit was incomplete")
    return (
        frames,
        determinant_histogram,
        phase_histogram,
        swap_checks,
    )


def heisenberg_probe():
    # T e_r=e_(r+1), M e_r=zeta^(-r)e_r.  Compare TM and zeta MT.
    checks = 0
    for root in range(P):
        tm_target = (root + 1) % P
        tm_coefficient = pow(ZETA, -root, FIELD)
        mt_target = (root + 1) % P
        mt_coefficient = pow(ZETA, -(root + 1), FIELD)
        require(tm_target == mt_target, "Weyl supports diverged")
        require(
            tm_coefficient == ZETA * mt_coefficient % FIELD,
            "Weyl commutator phase changed",
        )
        checks += 1

    cocycle_square = (1, 1, 1, ZETA)
    require(
        mobius(cocycle_square, FIELD) == (ZETA - 1) % FIELD,
        "cocycle escape changed",
    )
    cross_ratio = (
        cocycle_square[0]
        * cocycle_square[3]
        * pow(cocycle_square[1] * cocycle_square[2], -1, FIELD)
    ) % FIELD
    require(cross_ratio == ZETA, "cocycle cross-ratio changed")
    require(
        all(
            (left * right - ZETA * right * left) % FIELD != 0
            for left, right in product(range(1, FIELD), repeat=2)
        ),
        "a one-dimensional scalar Weyl pair appeared",
    )
    return checks, cocycle_square, cross_ratio


def address_and_cofiber_probes():
    # THM-2807: the three translation exponents form an honest boundary.
    pure = 13
    vertical = 169 * 4079
    diagonal = 689364
    require(pure + vertical == diagonal, "address triangle stopped closing")
    require(4079 % P == 10, "vertical transgression changed")
    exponent_boundary = (1, 0), (0, 4079), (1, 4079)
    require(
        (
            exponent_boundary[0][0] + exponent_boundary[1][0],
            exponent_boundary[0][1] + exponent_boundary[1][1],
        ) == exponent_boundary[2],
        "formal translation character acquired holonomy",
    )

    interval = (142004992589460, 142005019034340)
    minus = (142004190428100, 142004216872980)
    plus = (142082000080020, 142082026524900)
    require(
        minus[1] < interval[0] < interval[1] < plus[0],
        "cofiber copies stopped being disjoint",
    )
    address_shifts = (-13, 1248)
    depth_one_residues = tuple((shift // 13) % P for shift in address_shifts)
    require(
        depth_one_residues == (12, 5),
        "cofiber address residues changed",
    )

    epsilon_coefficients = tuple(
        (
            (comb(12, degree) if degree <= 12 else 0)
            + (comb(5, degree) if degree <= 5 else 0)
        ) % P
        for degree in range(P)
    )
    require(
        epsilon_coefficients[:2] == (2, 4),
        "cofiber copy augmentation/slope changed",
    )
    copy_bockstein = 9
    weighted_coefficients = tuple(
        copy_bockstein * value % P for value in epsilon_coefficients
    )
    require(
        weighted_coefficients[:2] == (5, 10),
        "weighted cofiber copy germ changed",
    )
    require(
        (copy_bockstein + copy_bockstein) % P == 5,
        "cofiber Bockstein sum changed",
    )
    return (
        (pure, vertical, diagonal),
        exponent_boundary,
        address_shifts,
        depth_one_residues,
        epsilon_coefficients,
        weighted_coefficients,
    )


def main():
    idempotent = exhaustive_idempotent_and_coarse_graining()
    sign_character = sign_and_character_probes()
    charged_atlas = charged_atlas_probe()
    quarter_turn = quarter_turn_probe()
    rees = rees_probe()
    projective_gauge = projective_square_gauge_probe()
    symplectic_frames = symplectic_heisenberg_frame_probe()
    heisenberg = heisenberg_probe()
    address = address_and_cofiber_probes()

    print("NON-IDEMPOTENT SAME-ATOM ALLOCATION PROBE")
    print("status=FINITE-EXACT scratch; abstract lemmas proved in REPORT.md")
    print(f"field={FIELD}; primitive_13th_root={ZETA}")
    print(
        "idempotent="
        f"field_idempotents={idempotent[0]};"
        f"fourfold_nonidempotent_unit_successes={idempotent[1]};"
        f"coarse_graining_checks={idempotent[2]}"
    )
    print(f"coarse_nonzero_examples={idempotent[3]}")
    print(
        "conditional_expectation_hostile="
        f"(B,P,Q,H)={idempotent[4]};mu={mobius(idempotent[4])};"
        "source_of_mu=both-absent atom, not fourfold-common atom"
    )
    print(
        "allocation_sign="
        f"flat_spectrum={sign_character[0]};"
        f"alternating_twist_spectrum={sign_character[1]};"
        "verdict=twist shifts trivial mass into mixed character"
    )
    print(f"character_faces_by_marked_origin={sign_character[2]}")
    print(
        "representative_gauge="
        f"vertex_factors={sign_character[3]};"
        "verdict=no common scalar descent"
    )
    print(
        "canonical_local_system_transport="
        f"raw={sign_character[4]};transported={sign_character[5]};"
        "covariant_mu=0"
    )
    print(
        "thm2593_charged_atlas="
        f"nontrivial_edges={charged_atlas[0]};"
        f"unit_edge_contrasts={charged_atlas[1]};"
        f"nonunit_edge_contrasts={charged_atlas[2]};"
        f"ordered_nonzero_mixed_products={charged_atlas[3]};"
        f"unit_mixed_products={charged_atlas[4]};"
        f"crt_zero_component_histogram={charged_atlas[5]};"
        "verdict=coefficient-positive-control but coboundary between "
        "different sheets"
    )
    print(
        "quarter_turn="
        f"v={quarter_turn[0]};Jv={quarter_turn[1]};"
        f"J2v={quarter_turn[2]};mixed_vector={quarter_turn[3]};"
        f"scalar_readout={quarter_turn[4]};"
        "verdict=minimal integral signed local-system escape"
    )
    print(
        "rees="
        f"central={rees[0]};valuations={rees[1]};"
        f"leading_units={rees[2]};valuation_curvature={rees[3]};"
        f"common_grade={rees[4]};"
        "verdict=ungraded_mu_is_orbit-cardinality_contrast_not_2-curvature"
    )
    print(
        "projective_square_gauge="
        f"nonzero_squares={projective_gauge[0]};"
        f"kappa_fibre_histogram={projective_gauge[1]};"
        f"anchored_orbit_sizes={projective_gauge[2]};"
        f"kappa_one_mu_histogram={projective_gauge[3]};"
        f"flat_example={projective_gauge[4]};"
        f"rephased_segre_example={projective_gauge[5]};"
        "verdict=kappa classifies row/column gauge orbits, while additive "
        "mu is not projectively intrinsic"
    )
    print(
        "symplectic_heisenberg_frames="
        f"ordered_nondegenerate_frames={symplectic_frames[0]};"
        f"determinant_histogram={symplectic_frames[1]};"
        f"phase_histogram={symplectic_frames[2]};"
        f"swap_inverse_checks={symplectic_frames[3]};"
        "verdict=each nontrivial central phase occurs 2184 times and "
        "frame swap inverts it"
    )
    print(
        "heisenberg_escape="
        f"weyl_checks={heisenberg[0]};"
        f"cocycle_square={heisenberg[1]};"
        f"cross_ratio={heisenberg[2]};"
        "verdict=nontrivial_gauge-invariant_joint_phase_requires_internal_fibre"
    )
    print(
        "thm2807_triangle="
        f"integer_edges={address[0]};formal_edges={address[1]};"
        "holonomy=1;vertical_residue=10"
    )
    print(
        "thm2771_cofiber_copies="
        f"address_shifts={address[2]};depth_one_residues={address[3]};"
        f"epsilon_coefficients={address[4]};"
        f"bockstein_weighted_coefficients={address[5]};"
        "verdict=two disjoint atoms give a unit copy-germ, not same-atom face"
    )
    print(
        "frontier=a fixed nontrivial edge phase is the cheapest additive "
        "escape, and a joint cocycle kappa!=1 is the cheapest "
        "edge-rephasing-invariant escape.  Existing phases are signed or "
        "are coboundaries between distinct sheets; positive averaging "
        "preserves bare-only provenance.  THM-2791/2807/cofiber sheets do "
        "not realize either escape on one full physical atom."
    )


if __name__ == "__main__":
    main()
