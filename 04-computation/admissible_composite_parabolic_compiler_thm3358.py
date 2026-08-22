#!/usr/bin/env python3
"""Exact companion for the THM-3358 composite parabolic compiler.

The two independent construction paths are:

1. Gaussian/Berggren path: enumerate every primitive representation
   N=a^2+b^2, compile its two determinant-matched cusp gates, and attach
   every unit normal prefix P_(N,alpha).
2. Arithmetic path: brute-force the roots modulo N, Hensel-lift them by
   Newton's formula modulo N^2, and scan all N lifts above every root.

Universe
--------
* every admissible odd N <= 1,200 (all prime divisors are 1 mod 4);
* the least rank-four grade N=32,045=5*13*17*29;
* every primitive representation, both A/D gates, every
  alpha in (Z/NZ)^*, and two lawful rows on every compiled unary ray;
* the entire normal fibre above every root, including nonunit and Hensel
  strata;
* explicit even, inert, ramified, partial-power, determinant-charge, and
  cusp-parity hostiles.

All checks use exact integer arithmetic and explicit exceptions, so normal
and optimized Python execute the same validation path.
"""

from collections import Counter
from hashlib import sha256
from itertools import product
from math import gcd, isqrt


U2 = ((2, -1), (1, 0))
A2 = ((2, 1), (1, 0))
D2 = ((1, 2), (0, 1))
PARAMETER_MATRICES = {"U": U2, "A": A2, "D": D2}

U3 = ((1, -2, 2), (2, -1, 2), (2, -2, 3))
A3 = ((1, 2, 2), (2, 1, 2), (2, 2, 3))
D3 = ((-1, 2, 2), (-2, 1, 2), (-2, 2, 3))
BERGGREN_MATRICES = {"U": U3, "A": A3, "D": D3}

I2 = ((1, 0), (0, 1))
SWAP = ((0, 1), (1, 0))
ROOT_PARAMETERS = (2, 1)
ROOT_TRIPLE = (3, 4, 5)
GRADE_LIMIT = 1_200
RANK_FOUR_GRADE = 32_045
CONTROL_3D_GRADES = {5, 25, 65, 169, 1_105, RANK_FOUR_GRADE}

EXPECTED_GRADE_COUNT = 146
EXPECTED_RANK_DISTRIBUTION = ((1, 101), (2, 43), (3, 1), (4, 1))
EXPECTED_REPRESENTATIONS = 199
EXPECTED_ORIENTED_CENTERS = 398
EXPECTED_UNIT_PREFIXES = 100_384
EXPECTED_NORMAL_POINTS = 743_214
EXPECTED_EXACT_BRANCHES = 553_216
EXPECTED_QUOTIENT_BRANCHES = 276_608
EXPECTED_WORD_ROWS = 1_106_432
EXPECTED_AFFINE_NORMAL = (
    44,
    10_814,
    9_024,
    4_512,
    "9a5f3058e13140091429c7deff7e17b2093b8ab9d12f7d6058d7fd3ab4431cd6",
)
EXPECTED_BOOLEAN_FACES = (
    584,
    1_592,
    3_620,
    "f170bc47620385e61628d5db387867f3a29bfe041010b7b8810ae74487cfd82a",
)
EXPECTED_EXACT_DIGEST = "385b513612da847bfed44fcd49fe26d0215c12be35353dab2aa77189a803e333"
EXPECTED_SEMANTIC_DIGEST = "a9cff64857122f215f3a7c38fcc34853246041973d98164ff714475035777711"


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def matrix_vector(matrix, vector):
    return tuple(
        sum(matrix[i][j] * vector[j] for j in range(len(vector)))
        for i in range(len(matrix))
    )


def matrix_multiply(left, right):
    rows = len(left)
    middle = len(right)
    columns = len(right[0])
    return tuple(
        tuple(
            sum(left[i][k] * right[k][j] for k in range(middle))
            for j in range(columns)
        )
        for i in range(rows)
    )


def transpose(matrix):
    return tuple(
        tuple(matrix[j][i] for j in range(len(matrix)))
        for i in range(len(matrix[0]))
    )


def determinant2(matrix):
    return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]


def word_matrix(word):
    matrix = I2
    for letter in word:
        matrix = matrix_multiply(PARAMETER_MATRICES[letter], matrix)
    return matrix


def apply_parameter_word(word, vector=ROOT_PARAMETERS):
    for letter in word:
        vector = matrix_vector(PARAMETER_MATRICES[letter], vector)
    return vector


def apply_berggren_word(word, vector=ROOT_TRIPLE):
    for letter in word:
        vector = matrix_vector(BERGGREN_MATRICES[letter], vector)
    return vector


def euclid_triple(parameters):
    m, n = parameters
    return m * m - n * n, 2 * m * n, m * m + n * n


def norm2(vector):
    return vector[0] * vector[0] + vector[1] * vector[1]


def dot(left, right):
    return left[0] * right[0] + left[1] * right[1]


def det(left, right):
    return left[0] * right[1] - left[1] * right[0]


def content(vector):
    return gcd(abs(vector[0]), abs(vector[1]))


def gauss_multiply(left, right):
    a, b = left
    c, d = right
    return a * c - b * d, a * d + b * c


def gauss_conjugate(value):
    return value[0], -value[1]


def C(index):
    return 2 * index * index + 2 * index + 1


def ceil_div(numerator, denominator):
    return -((-numerator) // denominator)


def valuation(number, prime):
    exponent = 0
    while number % prime == 0:
        exponent += 1
        number //= prime
    return exponent


def factorint(number):
    factors = []
    divisor = 2
    while divisor * divisor <= number:
        if number % divisor == 0:
            exponent = 0
            while number % divisor == 0:
                number //= divisor
                exponent += 1
            factors.append((divisor, exponent))
        divisor = 3 if divisor == 2 else divisor + 2
    if number > 1:
        factors.append((number, 1))
    return tuple(factors)


def euler_phi(number):
    result = number
    for prime, _ in factorint(number):
        result = result // prime * (prime - 1)
    return result


def admissible(number):
    return number > 1 and number % 2 == 1 and all(
        prime % 4 == 1 for prime, _ in factorint(number)
    )


def primitive_representations(number):
    representations = []
    for b in range(1, isqrt(number) + 1):
        a = isqrt(number - b * b)
        if (
            a > b
            and a * a + b * b == number
            and gcd(a, b) == 1
            and (a - b) % 2 == 1
        ):
            representations.append((a, b))
    return tuple(representations)


def odd_cusp_tail(target):
    """Unique R with B(R)(3,1)=target for a coprime positive odd pair."""
    x, y = target
    require(
        x > y > 0 and gcd(x, y) == 1 and x % 2 == y % 2 == 1,
        ("bad odd cusp", target),
    )
    inverse_letters = []
    while (x, y) != (3, 1):
        old_sum = x + y
        if x < 2 * y:
            inverse_letters.append("U")
            x, y = y, 2 * y - x
        elif x < 3 * y:
            inverse_letters.append("A")
            x, y = y, x - 2 * y
        elif x > 3 * y:
            inverse_letters.append("D")
            x, y = x - 2 * y, y
        else:
            require(y == 1, ("nonterminal ratio-three cusp", x, y))
            x, y = 3, 1
            continue
        require(
            x > y > 0 and x + y < old_sum,
            ("cusp descent did not decrease", target, x, y),
        )
        require(
            gcd(x, y) == 1 and x % 2 == y % 2 == 1,
            ("cusp descent lost type", target, x, y),
        )
    return "".join(reversed(inverse_letters))


def berggren_address(target):
    """Unique root-to-node word for a primitive opposite-parity pair."""
    x, y = target
    require(
        x > y > 0 and gcd(x, y) == 1 and (x - y) % 2 == 1,
        ("bad Berggren target", target),
    )
    inverse_letters = []
    while (x, y) != ROOT_PARAMETERS:
        old_sum = x + y
        if x < 2 * y:
            inverse_letters.append("U")
            x, y = y, 2 * y - x
        elif x == 2 * y:
            raise RuntimeError(("unexpected ratio-two boundary", target, x, y))
        elif x < 3 * y:
            inverse_letters.append("A")
            x, y = y, x - 2 * y
        elif x == 3 * y:
            raise RuntimeError(("unexpected ratio-three boundary", target, x, y))
        else:
            inverse_letters.append("D")
            x, y = x - 2 * y, y
        require(
            x > y > 0 and x + y < old_sum,
            ("Berggren descent did not decrease", target, x, y),
        )
        require(
            gcd(x, y) == 1 and (x - y) % 2 == 1,
            ("Berggren descent lost type", target, x, y),
        )
    word = "".join(reversed(inverse_letters))
    require(apply_parameter_word(word) == target, ("address mismatch", target, word))
    return word


def leading_u_length(word):
    length = 0
    for letter in word:
        if letter != "U":
            break
        length += 1
    return length


def roots_mod_grade(number):
    return tuple(index for index in range(number) if C(index) % number == 0)


def hensel_lift(number, root):
    """Newton lift of a C-root from N to N^2."""
    require(0 <= root < number and C(root) % number == 0, (number, root))
    derivative = 4 * root + 2
    require(gcd(derivative, number) == 1, (number, root, derivative))
    correction = (-(C(root) // number) * pow(derivative, -1, number)) % number
    lifted = (root + number * correction) % (number * number)
    require(C(lifted) % (number * number) == 0, (number, root, lifted))
    require(lifted % number == root, (number, root, lifted))
    return lifted


def compile_representation(number, representation):
    """Gaussian/cusp construction, independent of the root scan."""
    a, b = representation
    d = a * a - b * b
    e = 2 * a * b
    orientation = 1 if d > e else -1
    require(d != e, (number, representation, d, e))
    base_h = ((d, e), (-orientation * e, orientation * d))
    cusp = matrix_vector(base_h, (1, 1))
    require(cusp == (d + e, abs(d - e)), (number, representation, cusp))
    require(
        matrix_multiply(transpose(base_h), base_h)
        == ((number * number, 0), (0, number * number)),
        ("orthogonality", number, representation, base_h),
    )
    require(abs(determinant2(base_h)) == number * number, (number, base_h))
    require(
        cusp[0] > cusp[1] > 0
        and gcd(*cusp) == 1
        and cusp[0] % 2 == cusp[1] % 2 == 1,
        ("nonprimitive cusp", number, representation, cusp),
    )
    tail = odd_cusp_tail(cusp)
    require(apply_parameter_word(tail, (3, 1)) == cusp, (number, tail, cusp))

    records = []
    for gate in "AD":
        suffix = gate + tail
        suffix_matrix = word_matrix(suffix)
        eta = determinant2(suffix_matrix)
        require(abs(eta) == 1, (number, representation, gate, eta))
        require(
            matrix_vector(suffix_matrix, (1, 1)) == cusp,
            ("cusp collision", number, representation, gate),
        )
        toggle = base_h
        if determinant2(toggle) != eta * number * number:
            toggle = matrix_multiply(toggle, SWAP)
        require(
            determinant2(toggle) == eta * number * number,
            ("determinant gauge", number, representation, gate),
        )
        require(matrix_vector(toggle, (1, 1)) == cusp, (number, gate, toggle))
        require(
            matrix_multiply(transpose(toggle), toggle)
            == ((number * number, 0), (0, number * number)),
            ("gauged orthogonality", number, representation, gate),
        )

        second_column = suffix_matrix[0][1], suffix_matrix[1][1]
        transposed_column = matrix_vector(transpose(toggle), second_column)
        u = number - transposed_column[0], number - transposed_column[1]
        require(u[0] == u[1] + 1, ("nonconsecutive intercept", number, gate, u))
        raw_offset = u[1]
        center = raw_offset - number
        center_vector = center + 1, center
        require(
            center_vector
            == tuple(-value for value in transposed_column),
            ("center vector", number, representation, gate),
        )
        require(
            C(center) == number * number * norm2(second_column),
            ("Hensel center norm", number, representation, gate, center),
        )
        records.append(
            {
                "representation": representation,
                "gate": gate,
                "cusp": cusp,
                "tail": tail,
                "suffix": suffix,
                "suffix_matrix": suffix_matrix,
                "toggle": toggle,
                "second_column": second_column,
                "raw_offset": raw_offset,
                "center": center,
            }
        )

    left, right = records
    require(
        left["raw_offset"] + right["raw_offset"] + 1
        == 2 * number - number * number,
        ("raw gate shear", number, representation, left, right),
    )
    require(
        left["center"] + right["center"] + 1 == -(number * number),
        ("center antipode", number, representation, left, right),
    )
    return tuple(records)


def compile_grade(number):
    factors = factorint(number)
    rank = len(factors)
    representations = primitive_representations(number)
    require(
        len(representations) == 2 ** (rank - 1),
        ("representation count", number, factors, representations),
    )
    records = []
    for representation in representations:
        records.extend(compile_representation(number, representation))
    return representations, tuple(records)


def arithmetic_normal_atlas(number, factors, exact_digest):
    """Independent brute root/lift scan, with no Gaussian compiler input."""
    rank = len(factors)
    roots = roots_mod_grade(number)
    require(len(roots) == 2 ** rank, ("root count", number, factors, roots))
    centers = tuple(hensel_lift(number, root) for root in roots)
    require(len(set(centers)) == len(centers), (number, roots, centers))

    exact_offsets = set()
    quotient_coordinates = set()
    normal_points = 0
    hensel_points = 0
    for root, center in zip(roots, centers):
        local_hensel = []
        for lift_index in range(number):
            offset = root + number * lift_index
            value = C(offset)
            require(value % number == 0, (number, root, offset, value))
            difference = (offset - center) % (number * number)
            require(difference % number == 0, (number, root, center, offset))
            alpha = (difference // number) % number
            quotient_gcd = gcd(number, value // number)
            require(
                quotient_gcd == gcd(number, alpha),
                ("normal gcd law", number, root, center, offset, alpha, quotient_gcd),
            )
            if value % (number * number) == 0:
                local_hensel.append(offset)
            for prime, exponent in factors:
                alpha_depth = exponent if alpha == 0 else min(valuation(alpha, prime), exponent)
                require(
                    min(valuation(value, prime), 2 * exponent)
                    == exponent + alpha_depth,
                    ("truncated valuation", number, root, offset, alpha, prime, exponent),
                )
            if quotient_gcd == 1:
                nu = (value // number) % number
                require(gcd(nu, number) == 1, ("nonunit nu", number, offset, nu))
                require(
                    nu == (2 * alpha * (2 * center + 1)) % number,
                    ("normal nu formula", number, root, center, offset, alpha, nu),
                )
                derivative = 4 * offset + 2
                retracted = (
                    offset - value * pow(derivative, -1, number * number)
                ) % (number * number)
                require(retracted == center, ("Newton retraction", number, offset, center, retracted))
                reconstructed = (
                    center
                    + number * nu * pow(4 * center + 2, -1, number)
                ) % (number * number)
                require(
                    reconstructed == offset,
                    ("nu inverse", number, root, center, offset, nu, reconstructed),
                )
                reflected = (-1 - offset) % (number * number)
                require(C(reflected) % number == 0, (number, offset, reflected))
                reflected_nu = (C(reflected) // number) % number
                require(reflected_nu == nu, ("nu conjugation", number, offset, reflected, nu))
                reflected_root = (number - 1 - root) % number
                require(reflected % number == reflected_root, (number, root, reflected))
                parent_root = min(root, reflected_root)
                quotient_coordinates.add((parent_root, nu))
                exact_offsets.add(offset)
            normal_points += 1
        require(local_hensel == [center], ("nonunique Hensel lift", number, root, local_hensel, center))
        hensel_points += 1

    expected_exact = (2 ** rank) * euler_phi(number)
    require(len(exact_offsets) == expected_exact, (number, len(exact_offsets), expected_exact))
    require(
        len(quotient_coordinates) == expected_exact // 2,
        ("quotient coordinate count", number, len(quotient_coordinates), expected_exact // 2),
    )
    parent_roots = {
        min(root, (number - 1 - root) % number) for root in roots
    }
    expected_quotient_coordinates = {
        (parent_root, unit)
        for parent_root in parent_roots
        for unit in range(1, number)
        if gcd(unit, number) == 1
    }
    require(
        quotient_coordinates == expected_quotient_coordinates,
        ("parent-root times unit atlas", number),
    )
    require(
        {(-1 - offset) % (number * number) for offset in exact_offsets} == exact_offsets,
        ("exact shell not conjugation-stable", number),
    )
    require(hensel_points == 2 ** rank, (number, hensel_points))
    for offset in sorted(exact_offsets):
        exact_digest.update(f"{number}:{offset}\n".encode())
    return roots, centers, exact_offsets, normal_points, len(quotient_coordinates)


def prefix_bank(number, semantic_digest):
    prefixes = []
    q = (number - 1) // 2
    for alpha in range(1, number):
        if gcd(alpha, number) != 1:
            continue
        seed = number + alpha, alpha
        word = berggren_address(seed)
        require(apply_parameter_word(word) == seed, (number, alpha, seed, word))
        if alpha == 1:
            require(word == "D" * q, ("pure-D prefix", number, word, q))
        else:
            require(word != "D" * q, ("nonunique pure-D prefix", number, alpha, word))
        require(
            leading_u_length(word) < len(word),
            ("prefix is an all-U word", number, alpha, word),
        )
        semantic_digest.update((repr(("prefix", number, alpha, word)) + "\n").encode())
        prefixes.append((alpha, word, leading_u_length(word)))
    require(len(prefixes) == euler_phi(number), (number, len(prefixes), euler_phi(number)))
    return tuple(prefixes)


def audit_compiled_shell(
    number,
    factors,
    records,
    prefixes,
    arithmetic_offsets,
    semantic_digest,
):
    modulus = number * number
    compiler_offsets = set()
    word_rows = 0
    all_start_one = True
    max_suffix = (0, None)
    max_prefix = (0, None)
    three_dimensional_checks = 0

    for alpha, prefix, _ in prefixes:
        if len(prefix) > max_prefix[0]:
            max_prefix = len(prefix), alpha

    for record_index, record in enumerate(records):
        center = record["center"]
        center_residue = center % modulus
        coarse_root = center_residue % number
        require(C(center_residue) % modulus == 0, (number, record, center_residue))
        semantic_digest.update(
            (
                repr(
                    (
                        "compiler",
                        number,
                        record["representation"],
                        record["gate"],
                        record["cusp"],
                        record["tail"],
                        record["suffix"],
                        record["raw_offset"],
                        center,
                        center_residue,
                        coarse_root,
                    )
                )
                + "\n"
            ).encode()
        )
        if len(record["suffix"]) > max_suffix[0]:
            max_suffix = len(record["suffix"]), (
                record["representation"],
                record["gate"],
            )

        for prefix_index, (alpha, prefix, prefix_u_run) in enumerate(prefixes):
            raw_offset = center + alpha * number
            residue = raw_offset % modulus
            require(residue not in compiler_offsets, ("duplicate compiler residue", number, residue))
            compiler_offsets.add(residue)
            require(residue in arithmetic_offsets, ("compiler outside exact shell", number, residue))
            require(residue % number == coarse_root, (number, residue, coarse_root))

            start = max(1, ceil_div(2 - raw_offset, modulus))
            all_start_one = all_start_one and start == 1
            for s in (start, start + 1):
                t = modulus * s + raw_offset
                require(t >= 2, (number, record, alpha, s, t))
                source = t + 1, t
                numerator = matrix_vector(record["toggle"], source)
                require(
                    numerator[0] % number == numerator[1] % number == 0,
                    ("nonintegral toggle", number, record, alpha, s, source, numerator),
                )
                toggled = numerator[0] // number, numerator[1] // number

                seed = number + alpha, alpha
                intermediate = number * s + alpha, number * (s - 1) + alpha
                require(
                    intermediate
                    == (
                        seed[0] + (s - 1) * number,
                        seed[1] + (s - 1) * number,
                    ),
                    (number, alpha, s, seed, intermediate),
                )
                require(
                    gcd(*intermediate) == 1 and (intermediate[0] - intermediate[1]) % 2 == 1,
                    ("nonprimitive intermediate", number, alpha, s, intermediate),
                )
                target = matrix_vector(record["suffix_matrix"], intermediate)
                require(target == toggled, ("word/toggle mismatch", number, record, alpha, s))
                require(
                    target[0] > target[1] > 0
                    and gcd(*target) == 1
                    and (target[0] - target[1]) % 2 == 1,
                    ("target left primitive chamber", number, record, alpha, s, target),
                )
                hypotenuse = C(t)
                require(norm2(source) == hypotenuse == norm2(target), (number, alpha, s))
                require(
                    hypotenuse % number == 0 and gcd(hypotenuse // number, number) == 1,
                    ("exact grade", number, alpha, s, t, hypotenuse),
                )
                nu = (hypotenuse // number) % number
                require(
                    nu == (2 * alpha * (2 * center + 1)) % number,
                    ("compiler nu", number, record, alpha, s, nu),
                )
                retracted = (
                    t
                    - hypotenuse
                    * pow(4 * t + 2, -1, number * number)
                ) % (number * number)
                require(
                    retracted == center_residue,
                    ("compiler Newton retraction", number, record, alpha, s, retracted),
                )
                reconstructed = (
                    center_residue
                    + number * nu * pow(4 * center_residue + 2, -1, number)
                ) % (number * number)
                require(
                    reconstructed == t % (number * number),
                    ("compiler nu inverse", number, record, alpha, s, reconstructed),
                )
                for prime, exponent in factors:
                    require(
                        valuation(hypotenuse, prime) == exponent,
                        ("compiler valuation", number, alpha, s, prime, exponent, hypotenuse),
                    )

                product = gauss_multiply(source, target)
                conjugate_product = gauss_multiply(source, gauss_conjugate(target))
                actual_contents = sorted((content(product), content(conjugate_product)))
                expected_contents = sorted((number, hypotenuse // number))
                require(
                    actual_contents == expected_contents,
                    (
                        "folded content",
                        number,
                        record["representation"],
                        record["gate"],
                        alpha,
                        s,
                        actual_contents,
                        expected_contents,
                    ),
                )

                source_depth = t - 1
                target_depth = len(prefix) + (s - 1) + len(record["suffix"])
                path_cost = source_depth + target_depth - 2 * prefix_u_run
                require(path_cost >= 0, (number, alpha, s, path_cost))
                require(
                    abs(source_depth - target_depth)
                    == abs(
                        modulus * s
                        + raw_offset
                        - len(prefix)
                        - s
                        - len(record["suffix"])
                    ),
                    ("depth jump", number, alpha, s),
                )
                word_rows += 1

            if (
                number in CONTROL_3D_GRADES
                and record_index in {0, len(records) - 1}
                and prefix_index in {0, len(prefixes) - 1}
            ):
                for controlled_s in (start, start + 1):
                    full_word = (
                        prefix
                        + "U" * (controlled_s - 1)
                        + record["suffix"]
                    )
                    parameters = apply_parameter_word(full_word)
                    triple_from_parameters = euclid_triple(parameters)
                    triple_from_tree = apply_berggren_word(full_word)
                    require(
                        parameters
                        == matrix_vector(
                            record["suffix_matrix"],
                            (
                                number * controlled_s + alpha,
                                number * (controlled_s - 1) + alpha,
                            ),
                        ),
                        ("literal parameter word", number, alpha, record),
                    )
                    require(
                        triple_from_parameters == triple_from_tree,
                        ("2D/3D mismatch", number, alpha, record, full_word),
                    )
                    three_dimensional_checks += 1

    require(compiler_offsets == arithmetic_offsets, ("two-path mismatch", number))
    expected = (2 ** len(factors)) * euler_phi(number)
    require(len(compiler_offsets) == expected, (number, len(compiler_offsets), expected))
    return {
        "exact_branches": len(compiler_offsets),
        "word_rows": word_rows,
        "all_start_one": all_start_one,
        "max_prefix": max_prefix,
        "max_suffix": max_suffix,
        "three_dimensional_checks": three_dimensional_checks,
    }


def general_affine_normal_audit():
    """Cross-check the normal/nu law on tame non-U-spine affine quadratics."""
    controls = (
        ((1, 0), (1, 1)),
        ((1, 0), (2, 1)),
        ((1, 0), (3, 2)),
        ((1, -1), (4, 1)),
        ((0, 1), (2, 1)),
    )
    moduli = (5, 13, 17, 25, 29, 37, 41, 65, 85, 125, 145, 169, 185, 221, 325)
    digest = sha256()
    atlas_count = 0
    normal_points = 0
    exact_points = 0
    quotient_points = 0

    for x0, direction in controls:
        require(gcd(*direction) == 1, ("nonprimitive affine direction", direction))
        A = norm2(direction)
        h = dot(x0, direction)
        C0 = norm2(x0)
        charge = det(x0, direction)
        require(charge != 0, (x0, direction, charge))
        require(A * C0 - h * h == charge * charge, (x0, direction, A, h, C0))

        def quadratic(index):
            return A * index * index + 2 * h * index + C0

        for modulus in moduli:
            if not admissible(modulus) or gcd(modulus, 2 * A * charge) != 1:
                continue
            factors = factorint(modulus)
            rank = len(factors)
            roots = tuple(index for index in range(modulus) if quadratic(index) % modulus == 0)
            require(len(roots) == 2 ** rank, (x0, direction, modulus, roots))
            axis = (-2 * h * pow(A, -1, modulus * modulus)) % (modulus * modulus)
            exact = set()
            quotient_coordinates = set()

            for root in roots:
                derivative = 2 * (A * root + h)
                require(gcd(derivative, modulus) == 1, (x0, direction, modulus, root))
                correction = (
                    -(quadratic(root) // modulus) * pow(derivative, -1, modulus)
                ) % modulus
                center = (root + modulus * correction) % (modulus * modulus)
                require(quadratic(center) % (modulus * modulus) == 0, (modulus, root, center))
                require(center % modulus == root, (modulus, root, center))

                for lift_index in range(modulus):
                    offset = root + modulus * lift_index
                    value = quadratic(offset)
                    require(value % modulus == 0, (modulus, root, offset, value))
                    difference = (offset - center) % (modulus * modulus)
                    require(difference % modulus == 0, (modulus, root, center, offset))
                    alpha = (difference // modulus) % modulus
                    require(
                        gcd(modulus, value // modulus) == gcd(modulus, alpha),
                        ("general affine normal gcd", x0, direction, modulus, offset, alpha),
                    )
                    if gcd(alpha, modulus) == 1:
                        nu = (value // modulus) % modulus
                        require(
                            nu == (2 * alpha * (A * center + h)) % modulus,
                            ("general affine nu", x0, direction, modulus, offset, alpha, nu),
                        )
                        derivative_at_offset = 2 * (A * offset + h)
                        retracted = (
                            offset
                            - value
                            * pow(derivative_at_offset, -1, modulus * modulus)
                        ) % (modulus * modulus)
                        require(
                            retracted == center,
                            ("general affine Newton retraction", x0, direction, modulus, offset),
                        )
                        reconstructed = (
                            center
                            + modulus
                            * nu
                            * pow(2 * (A * center + h), -1, modulus)
                        ) % (modulus * modulus)
                        require(
                            reconstructed == offset,
                            ("general affine nu inverse", x0, direction, modulus, offset),
                        )
                        reflected = (axis - offset) % (modulus * modulus)
                        require(quadratic(reflected) % modulus == 0, (modulus, offset, reflected))
                        reflected_nu = (quadratic(reflected) // modulus) % modulus
                        require(
                            reflected_nu == nu,
                            ("general affine nu reflection", x0, direction, modulus, offset, reflected),
                        )
                        reflected_root = reflected % modulus
                        quotient_coordinates.add((min(root, reflected_root), nu))
                        exact.add(offset)
                    normal_points += 1

            expected_exact = (2 ** rank) * euler_phi(modulus)
            require(len(exact) == expected_exact, (x0, direction, modulus, len(exact)))
            require(
                {(axis - offset) % (modulus * modulus) for offset in exact} == exact,
                ("general affine reflection closure", x0, direction, modulus),
            )
            require(
                len(quotient_coordinates) == expected_exact // 2,
                ("general affine quotient", x0, direction, modulus, quotient_coordinates),
            )
            parent_roots = {
                min(root, (axis - root) % modulus) for root in roots
            }
            expected_coordinates = {
                (parent_root, unit)
                for parent_root in parent_roots
                for unit in range(1, modulus)
                if gcd(unit, modulus) == 1
            }
            require(
                quotient_coordinates == expected_coordinates,
                ("general affine parent-root times unit", x0, direction, modulus),
            )
            digest.update(
                (
                    repr(
                        (
                            x0,
                            direction,
                            A,
                            h,
                            charge,
                            modulus,
                            roots,
                            tuple(sorted(exact)),
                            tuple(sorted(quotient_coordinates)),
                        )
                    )
                    + "\n"
                ).encode()
            )
            atlas_count += 1
            exact_points += len(exact)
            quotient_points += len(quotient_coordinates)

    return atlas_count, normal_points, exact_points, quotient_points, digest.hexdigest()


def unary_word_invariant_audit():
    """Exhaust the universal fixed-prefix/unary-U/fixed-suffix identity."""
    words = ("",) + tuple(
        "".join(letters)
        for length in range(1, 5)
        for letters in product("UAD", repeat=length)
    )
    require(len(words) == 121 and len(set(words)) == len(words), len(words))

    prefix_data = []
    suffix_data = []
    for word in words:
        prefix_point = apply_parameter_word(word)
        delta = prefix_point[0] - prefix_point[1]
        require(delta != 0, ("zero unary coefficient", word, prefix_point))
        prefix_data.append((word, prefix_point, delta))

        direction = apply_parameter_word(word, (1, 1))
        require(
            gcd(*direction) == 1
            and direction[0] % 2 == direction[1] % 2 == 1,
            ("non-odd primitive suffix direction", word, direction),
        )
        suffix_data.append((word, direction))

    digest = sha256()
    row_count = 0
    for prefix, prefix_point, delta in prefix_data:
        for suffix, direction in suffix_data:
            initial = apply_parameter_word(suffix, prefix_point)
            for unary_steps in range(6):
                unary_point = apply_parameter_word("U" * unary_steps, prefix_point)
                require(
                    unary_point
                    == (
                        prefix_point[0] + unary_steps * delta,
                        prefix_point[1] + unary_steps * delta,
                    ),
                    ("unary transvection", prefix, unary_steps, unary_point),
                )
                point = apply_parameter_word(suffix, unary_point)
                expected = (
                    initial[0] + unary_steps * delta * direction[0],
                    initial[1] + unary_steps * delta * direction[1],
                )
                require(
                    point == expected,
                    ("universal unary word identity", prefix, suffix, unary_steps),
                )
                require(
                    abs(det(point, direction)) == abs(delta),
                    ("unary determinant charge", prefix, suffix, unary_steps),
                )
                digest.update(
                    (
                        repr(
                            (
                                prefix,
                                suffix,
                                unary_steps,
                                prefix_point,
                                delta,
                                direction,
                                point,
                            )
                        )
                        + "\n"
                    ).encode()
                )
                row_count += 1

    require(row_count == len(words) * len(words) * 6, row_count)
    return len(words), len(words) * len(words), row_count, digest.hexdigest()


def unitary_boolean_face_audit(grades):
    """Compile every unitary-divisor vertex on selected full-grade sources."""
    compiler_cache = {}
    digest = sha256()
    source_count = 0
    vertex_count = 0
    pair_count = 0

    def root_compilers(grade):
        if grade not in compiler_cache:
            _, grade_records = compile_grade(grade)
            compiler_cache[grade] = {
                record["center"] % grade: record for record in grade_records
            }
            require(
                len(compiler_cache[grade]) == 2 ** len(factorint(grade)),
                ("restricted compiler root count", grade, compiler_cache[grade]),
            )
        return compiler_cache[grade]

    for number in grades:
        factors = factorint(number)
        prime_power_blocks = tuple(prime ** exponent for prime, exponent in factors)
        unitary_divisors = []
        for mask in range(1 << len(prime_power_blocks)):
            divisor = 1
            for index, block in enumerate(prime_power_blocks):
                if (mask >> index) & 1:
                    divisor *= block
            unitary_divisors.append((mask, divisor))

        _, records = compile_grade(number)
        selected_records = (records[0], records[-1])
        for record in selected_records:
            for alpha in (1, number - 1):
                raw_offset = record["center"] + alpha * number
                start = max(1, ceil_div(2 - raw_offset, number * number))
                source_index = number * number * start + raw_offset
                source = source_index + 1, source_index
                hypotenuse = C(source_index)
                require(
                    hypotenuse % number == 0
                    and gcd(hypotenuse // number, number) == 1
                    and hypotenuse // number > 1,
                    ("Boolean-face source grade", number, record, alpha, source_index),
                )
                full_retraction = (
                    source_index
                    - hypotenuse
                    * pow(4 * source_index + 2, -1, number * number)
                ) % (number * number)
                require(
                    full_retraction == record["center"] % (number * number),
                    ("full-grade face retraction", number, source_index, full_retraction),
                )
                full_nu = (hypotenuse // number) % number

                targets = {}
                restricted_metadata = []
                for mask, divisor in unitary_divisors:
                    if divisor == 1:
                        target = source
                        restricted_metadata.append((mask, divisor, 0, "", 0))
                    else:
                        restricted_record = root_compilers(divisor)[source_index % divisor]
                        center = restricted_record["center"]
                        restricted_retraction = (
                            source_index
                            - hypotenuse
                            * pow(4 * source_index + 2, -1, divisor * divisor)
                        ) % (divisor * divisor)
                        require(
                            restricted_retraction == center % (divisor * divisor),
                            ("unitary retraction", number, divisor, source_index),
                        )
                        require(
                            restricted_retraction == full_retraction % (divisor * divisor),
                            ("unitary retraction restriction", number, divisor, source_index),
                        )
                        restricted_nu = (hypotenuse // divisor) % divisor
                        require(
                            restricted_nu
                            == ((number // divisor) * full_nu) % divisor,
                            ("unitary nu restriction", number, divisor, source_index),
                        )
                        require((source_index - center) % divisor == 0, (number, divisor, source_index))
                        restricted_alpha = ((source_index - center) // divisor) % divisor
                        require(
                            1 <= restricted_alpha < divisor
                            and gcd(restricted_alpha, divisor) == 1,
                            ("restricted normal unit", number, divisor, source_index, restricted_alpha),
                        )
                        restricted_raw = center + restricted_alpha * divisor
                        require(
                            (source_index - restricted_raw) % (divisor * divisor) == 0,
                            ("restricted ray residue", number, divisor, source_index, restricted_raw),
                        )
                        restricted_s = (source_index - restricted_raw) // (divisor * divisor)
                        require(
                            restricted_s >= 1,
                            ("restricted compiler before lawful tail", number, divisor, source_index, restricted_s),
                        )
                        prefix = berggren_address((divisor + restricted_alpha, restricted_alpha))
                        intermediate = (
                            divisor * restricted_s + restricted_alpha,
                            divisor * (restricted_s - 1) + restricted_alpha,
                        )
                        target = matrix_vector(restricted_record["suffix_matrix"], intermediate)
                        numerator = matrix_vector(restricted_record["toggle"], source)
                        require(
                            numerator[0] % divisor == numerator[1] % divisor == 0,
                            ("restricted toggle integrality", number, divisor, source_index),
                        )
                        require(
                            target == (numerator[0] // divisor, numerator[1] // divisor),
                            ("restricted word/toggle", number, divisor, source_index),
                        )
                        require(
                            apply_parameter_word(prefix) == (divisor + restricted_alpha, restricted_alpha),
                            ("restricted prefix", number, divisor, restricted_alpha, prefix),
                        )
                        restricted_metadata.append(
                            (mask, divisor, restricted_alpha, prefix, restricted_s)
                        )
                    require(norm2(target) == hypotenuse, (number, divisor, target, hypotenuse))
                    targets[mask] = target

                require(len(targets) == 2 ** len(factors), (number, targets))
                masks = tuple(sorted(targets))
                for left_index, left_mask in enumerate(masks):
                    for right_mask in masks[left_index:]:
                        symmetric_divisor = 1
                        symmetric_mask = left_mask ^ right_mask
                        for index, block in enumerate(prime_power_blocks):
                            if (symmetric_mask >> index) & 1:
                                symmetric_divisor *= block
                        left_target = targets[left_mask]
                        right_target = targets[right_mask]
                        pair_contents = sorted(
                            (
                                content(gauss_multiply(left_target, right_target)),
                                content(
                                    gauss_multiply(
                                        left_target,
                                        gauss_conjugate(right_target),
                                    )
                                ),
                            )
                        )
                        require(
                            pair_contents
                            == sorted((symmetric_divisor, hypotenuse // symmetric_divisor)),
                            (
                                "Boolean-face pair",
                                number,
                                source_index,
                                left_mask,
                                right_mask,
                                pair_contents,
                                symmetric_divisor,
                            ),
                        )
                        if left_mask != right_mask:
                            require(1 not in pair_contents, ("parent-face collision", number, pair_contents))
                        pair_count += 1

                digest.update(
                    (
                        repr(
                            (
                                number,
                                record["representation"],
                                record["gate"],
                                alpha,
                                source_index,
                                hypotenuse,
                                tuple(restricted_metadata),
                                tuple(sorted(targets.items())),
                            )
                        )
                        + "\n"
                    ).encode()
                )
                source_count += 1
                vertex_count += len(targets)

    return source_count, vertex_count, pair_count, digest.hexdigest()


def hostile_audit():
    require(roots_mod_grade(3) == (), ("inert root hostile", roots_mod_grade(3)))
    require(roots_mod_grade(45) == (), ("inert composite roots", roots_mod_grade(45)))

    even_representation = (3, 1)
    d = even_representation[0] ** 2 - even_representation[1] ** 2
    e = 2 * even_representation[0] * even_representation[1]
    even_cusp = d + e, abs(d - e)
    require(norm2(even_representation) == 10, even_representation)
    require(gcd(*even_cusp) == 2 and (10 - 1) % 2 == 1, ("even hostile", even_cusp))

    inert_representation = (6, 3)
    d = inert_representation[0] ** 2 - inert_representation[1] ** 2
    e = 2 * inert_representation[0] * inert_representation[1]
    inert_cusp = d + e, abs(d - e)
    require(norm2(inert_representation) == 45, inert_representation)
    require(gcd(*inert_cusp) == 9, ("nonprimitive inert cusp", inert_cusp))

    bad_split_representation = (10, 5)
    d = bad_split_representation[0] ** 2 - bad_split_representation[1] ** 2
    e = 2 * bad_split_representation[0] * bad_split_representation[1]
    bad_split_cusp = d + e, abs(d - e)
    require(norm2(bad_split_representation) == 125, bad_split_representation)
    require(gcd(*bad_split_cusp) == 25, ("nonprimitive split cusp", bad_split_cusp))
    require(primitive_representations(125) == ((11, 2),), primitive_representations(125))

    _, five_records = compile_grade(5)
    raw = tuple(sorted(record["raw_offset"] % 25 for record in five_records))
    centered = tuple(sorted(record["center"] % 25 for record in five_records))
    require(raw == (1, 8), ("raw p=5 offsets", raw))
    require(centered == (3, 21), ("centered p=5 offsets", centered))
    require((25 - 1 - raw[0]) % 25 != raw[1], ("raw offsets accidentally antipodal", raw))
    require((25 - 1 - centered[0]) % 25 == centered[1], ("centered antipode", centered))
    exact_five = tuple(
        offset for offset in range(25)
        if C(offset) % 5 == 0 and gcd(C(offset) // 5, 5) == 1
    )
    require(exact_five == (1, 6, 8, 11, 13, 16, 18, 23), exact_five)

    ramified_c_direction = (1, 0)
    ramified_c_origin = (0, 5)
    ramified_c_modulus = 5
    ramified_c_A = norm2(ramified_c_direction)
    ramified_c_h = dot(ramified_c_origin, ramified_c_direction)
    ramified_c_C0 = norm2(ramified_c_origin)
    ramified_c_charge = det(ramified_c_origin, ramified_c_direction)

    def ramified_c_quadratic(index):
        return (
            ramified_c_A * index * index
            + 2 * ramified_c_h * index
            + ramified_c_C0
        )

    ramified_c_roots = tuple(
        index
        for index in range(ramified_c_modulus)
        if ramified_c_quadratic(index) % ramified_c_modulus == 0
    )
    ramified_c_lifts = tuple(
        index
        for index in range(ramified_c_modulus * ramified_c_modulus)
        if ramified_c_quadratic(index)
        % (ramified_c_modulus * ramified_c_modulus)
        == 0
    )
    require(ramified_c_charge == -5, ramified_c_charge)
    require(ramified_c_roots == (0,), ramified_c_roots)
    require(ramified_c_lifts == (0, 5, 10, 15, 20), ramified_c_lifts)

    ramified_A_direction = (1, 2)
    ramified_A_origin = (1, 0)
    ramified_A_modulus = 5
    ramified_A_A = norm2(ramified_A_direction)
    ramified_A_h = dot(ramified_A_origin, ramified_A_direction)
    ramified_A_C0 = norm2(ramified_A_origin)
    ramified_A_charge = det(ramified_A_origin, ramified_A_direction)

    def ramified_A_quadratic(index):
        return (
            ramified_A_A * index * index
            + 2 * ramified_A_h * index
            + ramified_A_C0
        )

    ramified_A_roots = tuple(
        index
        for index in range(ramified_A_modulus)
        if ramified_A_quadratic(index) % ramified_A_modulus == 0
    )
    ramified_A_lifts = tuple(
        index
        for index in range(ramified_A_modulus * ramified_A_modulus)
        if ramified_A_quadratic(index)
        % (ramified_A_modulus * ramified_A_modulus)
        == 0
    )
    require(ramified_A_A == 5 and ramified_A_charge == 2, (ramified_A_A, ramified_A_charge))
    require(ramified_A_roots == (2,), ramified_A_roots)
    require(ramified_A_lifts == (2,), ramified_A_lifts)

    partial_number = 25
    partial_index = 428
    partial_divisor = 5
    require(
        C(partial_index) % partial_number == 0
        and gcd(C(partial_index) // partial_number, partial_number) == 1,
        ("source is not exact 25-grade", partial_index, C(partial_index)),
    )
    partial_record = next(
        record
        for record in five_records
        if record["center"] % partial_divisor == partial_index % partial_divisor
    )
    partial_numerator = matrix_vector(
        partial_record["toggle"],
        (partial_index + 1, partial_index),
    )
    require(
        partial_numerator[0] % partial_divisor
        == partial_numerator[1] % partial_divisor
        == 0,
        ("partial-power nonintegrality", partial_numerator),
    )
    partial_target = (
        partial_numerator[0] // partial_divisor,
        partial_numerator[1] // partial_divisor,
    )
    require(partial_target == (600, 85), partial_target)
    require(content(partial_target) == 5, partial_target)

    charge_direction = (1, 1)
    charge_origin = (2, -1)
    charge_modulus = 5
    charge_allocation = ((3, -4), (4, 3))
    charge_root = 0
    require(
        norm2(charge_direction) == 2
        and dot(charge_origin, charge_direction) == 1
        and norm2(charge_origin) == 5
        and det(charge_origin, charge_direction) == 3,
        ("charge hostile quadratic", charge_origin, charge_direction),
    )
    require(
        gcd(charge_modulus, 2 * norm2(charge_direction) * det(charge_origin, charge_direction))
        == 1,
        "charge hostile is not tame",
    )
    require(
        matrix_multiply(transpose(charge_allocation), charge_allocation)
        == ((charge_modulus * charge_modulus, 0), (0, charge_modulus * charge_modulus)),
        charge_allocation,
    )
    charge_rotated_direction = matrix_vector(charge_allocation, charge_direction)

    def charge_rotated_point(step_index):
        affine_index = charge_root + charge_modulus * charge_modulus * step_index
        affine_point = (
            charge_origin[0] + affine_index * charge_direction[0],
            charge_origin[1] + affine_index * charge_direction[1],
        )
        numerator = matrix_vector(charge_allocation, affine_point)
        require(
            numerator[0] % charge_modulus == numerator[1] % charge_modulus == 0,
            ("charge hostile allocation", step_index, numerator),
        )
        return numerator[0] // charge_modulus, numerator[1] // charge_modulus

    charge_initial = charge_rotated_point(0)
    charge_next = charge_rotated_point(1)
    charge_step = (
        charge_next[0] - charge_initial[0],
        charge_next[1] - charge_initial[1],
    )
    require(gcd(*charge_rotated_direction) == 1, charge_rotated_direction)
    require(
        charge_step
        == (
            charge_modulus * charge_rotated_direction[0],
            charge_modulus * charge_rotated_direction[1],
        ),
        charge_step,
    )
    require(
        abs(det(charge_initial, charge_rotated_direction)) == 15,
        (charge_initial, charge_rotated_direction),
    )

    parity_direction = (2, 1)
    parity_origin = (1, 0)
    parity_modulus = 13
    parity_allocation = ((5, -12), (12, 5))
    parity_root = 9
    require(
        norm2(parity_direction) == 5
        and det(parity_origin, parity_direction) == 1
        and gcd(
            parity_modulus,
            2 * norm2(parity_direction) * det(parity_origin, parity_direction),
        )
        == 1,
        ("parity hostile is not tame", parity_origin, parity_direction),
    )
    parity_source = (
        parity_origin[0] + parity_root * parity_direction[0],
        parity_origin[1] + parity_root * parity_direction[1],
    )
    parity_numerator = matrix_vector(parity_allocation, parity_source)
    require(
        parity_numerator[0] % parity_modulus
        == parity_numerator[1] % parity_modulus
        == 0,
        ("parity hostile allocation", parity_numerator),
    )
    parity_rotated_direction = matrix_vector(parity_allocation, parity_direction)
    require(
        (parity_rotated_direction[0] - parity_rotated_direction[1]) % 2 == 1,
        ("allocation did not preserve opposite parity", parity_rotated_direction),
    )
    primitive_parity_allocations = tuple(
        ((real, -imaginary), (imaginary, real))
        for real in range(-parity_modulus, parity_modulus + 1)
        for imaginary in range(-parity_modulus, parity_modulus + 1)
        if real * real + imaginary * imaginary == parity_modulus * parity_modulus
        and gcd(abs(real), abs(imaginary)) == 1
    )
    require(len(primitive_parity_allocations) == 8, primitive_parity_allocations)
    parity_directions = tuple(
        matrix_vector(allocation, parity_direction)
        for allocation in primitive_parity_allocations
    )
    require(
        all(
            gcd(*rotated_direction) == 1
            and (rotated_direction[0] - rotated_direction[1]) % 2 == 1
            for rotated_direction in parity_directions
        ),
        ("primitive allocation parity bank", parity_directions),
    )
    require(
        all(
            (apply_parameter_word(word, (1, 1))[0] % 2,
             apply_parameter_word(word, (1, 1))[1] % 2)
            == (1, 1)
            for length in range(5)
            for word in ("".join(letters) for letters in product("UAD", repeat=length))
        ),
        "bounded suffix bank left the odd cusp",
    )
    return {
        "even": (10, even_cusp),
        "inert": (3, 45, inert_cusp),
        "bad_split": (125, bad_split_cusp),
        "p5_raw": raw,
        "p5_centered": centered,
        "p5_exact": exact_five,
        "ramified_c": (ramified_c_roots, ramified_c_lifts),
        "ramified_A": (ramified_A_roots, ramified_A_lifts),
        "partial_power": (partial_number, partial_index, partial_divisor, partial_target, content(partial_target)),
        "charge": (
            charge_root,
            charge_initial,
            charge_rotated_direction,
            charge_step,
            abs(det(charge_initial, charge_rotated_direction)),
        ),
        "parity": (
            parity_root,
            parity_source,
            parity_rotated_direction,
            tuple(sorted(parity_directions)),
        ),
    }


def main():
    grades = tuple(
        number
        for number in range(3, GRADE_LIMIT + 1, 2)
        if admissible(number)
    ) + (RANK_FOUR_GRADE,)
    require(len(grades) == EXPECTED_GRADE_COUNT, len(grades))
    require(len(set(grades)) == len(grades), grades)

    rank_distribution = Counter()
    exact_digest = sha256()
    semantic_digest = sha256()
    total_representations = 0
    total_centers = 0
    total_prefixes = 0
    total_normal_points = 0
    total_exact_branches = 0
    total_quotient_branches = 0
    total_word_rows = 0
    total_3d_checks = 0
    all_start_one = True
    global_max_prefix = (0, None)
    global_max_suffix = (0, None)
    summaries = {}

    for number in grades:
        factors = factorint(number)
        rank = len(factors)
        rank_distribution[rank] += 1
        require(all(prime % 4 == 1 for prime, _ in factors), (number, factors))

        (
            roots,
            arithmetic_centers,
            arithmetic_offsets,
            normal_points,
            quotient_branches,
        ) = arithmetic_normal_atlas(number, factors, exact_digest)
        representations, records = compile_grade(number)
        compiler_centers = tuple(sorted(record["center"] % (number * number) for record in records))
        require(
            compiler_centers == tuple(sorted(arithmetic_centers)),
            ("compiler/Hensel center mismatch", number, compiler_centers, arithmetic_centers),
        )
        compiler_roots = tuple(sorted(record["center"] % number for record in records))
        require(compiler_roots == roots, ("compiler coarse roots", number, compiler_roots, roots))

        prefixes = prefix_bank(number, semantic_digest)
        compiled = audit_compiled_shell(
            number,
            factors,
            records,
            prefixes,
            arithmetic_offsets,
            semantic_digest,
        )

        total_representations += len(representations)
        total_centers += len(records)
        total_prefixes += len(prefixes)
        total_normal_points += normal_points
        total_exact_branches += compiled["exact_branches"]
        total_quotient_branches += quotient_branches
        total_word_rows += compiled["word_rows"]
        total_3d_checks += compiled["three_dimensional_checks"]
        all_start_one = all_start_one and compiled["all_start_one"]
        if compiled["max_prefix"][0] > global_max_prefix[0]:
            global_max_prefix = compiled["max_prefix"][0], (
                number,
                compiled["max_prefix"][1],
            )
        if compiled["max_suffix"][0] > global_max_suffix[0]:
            global_max_suffix = compiled["max_suffix"][0], (
                number,
                compiled["max_suffix"][1],
            )

        summary = (
            rank,
            len(representations),
            len(roots),
            normal_points,
            compiled["exact_branches"],
            quotient_branches,
            len(prefixes),
            compiled["max_prefix"],
            compiled["max_suffix"],
        )
        summaries[number] = summary
        semantic_digest.update((repr(("grade", number, factors, summary)) + "\n").encode())

    require(tuple(sorted(rank_distribution.items())) == EXPECTED_RANK_DISTRIBUTION, rank_distribution)
    require(total_representations == EXPECTED_REPRESENTATIONS, total_representations)
    require(total_centers == EXPECTED_ORIENTED_CENTERS, total_centers)
    require(total_prefixes == EXPECTED_UNIT_PREFIXES, total_prefixes)
    require(total_normal_points == EXPECTED_NORMAL_POINTS, total_normal_points)
    require(total_exact_branches == EXPECTED_EXACT_BRANCHES, total_exact_branches)
    require(total_quotient_branches == EXPECTED_QUOTIENT_BRANCHES, total_quotient_branches)
    require(total_word_rows == EXPECTED_WORD_ROWS, total_word_rows)

    rank_four_roots = roots_mod_grade(RANK_FOUR_GRADE)
    require(
        rank_four_roots
        == (
            1081, 3546, 4851, 5501, 7316, 7966, 9271, 11736,
            20308, 22773, 24078, 24728, 26543, 27193, 28498, 30963,
        ),
        rank_four_roots,
    )
    require(summaries[25][:7] == (1, 1, 2, 50, 40, 20, 20), summaries[25])
    require(summaries[125][:7] == (1, 1, 2, 250, 200, 100, 100), summaries[125])
    require(summaries[169][:7] == (1, 1, 2, 338, 312, 156, 156), summaries[169])
    require(
        summaries[RANK_FOUR_GRADE][:7]
        == (4, 8, 16, 512_720, 344_064, 172_032, 21_504),
        summaries[RANK_FOUR_GRADE],
    )

    affine_normal = general_affine_normal_audit()
    semantic_digest.update((repr(("general_affine_normal", affine_normal)) + "\n").encode())
    unary_words = unary_word_invariant_audit()
    semantic_digest.update((repr(("unary_word_invariant", unary_words)) + "\n").encode())
    boolean_faces = unitary_boolean_face_audit(grades)
    semantic_digest.update((repr(("unitary_boolean_faces", boolean_faces)) + "\n").encode())
    hostiles = hostile_audit()
    semantic_digest.update((repr(("hostiles", hostiles)) + "\n").encode())
    semantic_digest.update(
        (
            repr(
                (
                    "totals",
                    tuple(sorted(rank_distribution.items())),
                    total_representations,
                    total_centers,
                    total_prefixes,
                    total_normal_points,
                    total_exact_branches,
                    total_quotient_branches,
                    total_word_rows,
                    total_3d_checks,
                    all_start_one,
                    global_max_prefix,
                    global_max_suffix,
                )
            )
            + "\n"
        ).encode()
    )

    exact_hash = exact_digest.hexdigest()
    semantic_hash = semantic_digest.hexdigest()
    if EXPECTED_EXACT_DIGEST:
        require(exact_hash == EXPECTED_EXACT_DIGEST, (exact_hash, EXPECTED_EXACT_DIGEST))
    if EXPECTED_SEMANTIC_DIGEST:
        require(
            semantic_hash == EXPECTED_SEMANTIC_DIGEST,
            (semantic_hash, EXPECTED_SEMANTIC_DIGEST),
        )

    print("THM-3358 ADMISSIBLE COMPOSITE PARABOLIC COMPILER AUDIT")
    print(
        "grades", len(grades), "rank_distribution", tuple(sorted(rank_distribution.items())),
        "max_grade", max(grades),
    )
    print(
        "primitive_representations", total_representations,
        "oriented_hensel_centers", total_centers,
        "unit_prefixes", total_prefixes,
    )
    print(
        "normal_fibre_points", total_normal_points,
        "exact_grade_branches", total_exact_branches,
        "conjugation_quotient_branches", total_quotient_branches,
        "compiled_word_rows", total_word_rows,
    )
    print(
        "two_path_center_and_exact_shell_match=True",
        "all_sampled_branch_starts_equal_one", all_start_one,
        "independent_2d_3d_controls", total_3d_checks,
    )
    print("prime_power_controls")
    for number in (25, 125, 169):
        print(" ", number, summaries[number][:7])
    print("rank_four_control", RANK_FOUR_GRADE, summaries[RANK_FOUR_GRADE][:7])
    print("rank_four_roots", rank_four_roots)
    print("max_prefix", global_max_prefix, "max_suffix", global_max_suffix)
    print("general_affine_normal", affine_normal)
    print("unary_word_invariant", unary_words)
    print("unitary_boolean_faces", boolean_faces)
    print(
        "hostiles",
        "even", hostiles["even"],
        "inert", hostiles["inert"],
        "bad_split", hostiles["bad_split"],
    )
    print(
        "p5_raw", hostiles["p5_raw"],
        "p5_centered", hostiles["p5_centered"],
        "p5_exact", hostiles["p5_exact"],
    )
    print(
        "affine_hostiles",
        "ramified_c", hostiles["ramified_c"],
        "ramified_A", hostiles["ramified_A"],
        "partial_power", hostiles["partial_power"],
    )
    print(
        "compiler_no_go_hostiles",
        "charge", hostiles["charge"],
        "parity", hostiles["parity"],
    )
    print("exact_shell_sha256", exact_hash)
    print("semantic_sha256", semantic_hash)
    print("normalization=all checks exact; compare stdout after LF normalization")
    print("ALL CHECKS PASS")


if __name__ == "__main__":
    main()
