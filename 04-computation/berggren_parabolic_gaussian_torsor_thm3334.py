#!/usr/bin/env python3
"""Exact consequence audit for THM-3334.

Universe
--------
* the first 12,500 points of the consecutive-parameter Berggren U-spine;
* every primitive Euclid parameter pair 1 <= n < m <= 80;
* every odd fixed-hypotenuse fibre c <= 10,000;
* the complete Berggren tree below hypotenuse 1,105;
* all 64 tournaments on the first four-parent collision fibre.

All truth-bearing arithmetic is integral.  The checks use explicit exceptions,
so ``python`` and ``python -O`` execute the same validation path.
"""

from itertools import combinations, product
from math import gcd, isqrt


U = ((1, -2, 2), (2, -1, 2), (2, -2, 3))
A = ((1, 2, 2), (2, 1, 2), (2, 2, 3))
D = ((-1, 2, 2), (-2, 1, 2), (-2, 2, 3))
PARAMETER_U = ((2, -1), (1, 0))
IDENTITY3 = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
LORENTZ = ((1, 0, 0), (0, 1, 0), (0, 0, -1))
SPINE_RECORD_BOUND = 12_500
FIXED_HYPOTENUSE_BOUND = 10_000
EUCLID_BOUND = 80


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def add_matrix(left, right):
    return tuple(
        tuple(left[i][j] + right[i][j] for j in range(len(left[0])))
        for i in range(len(left))
    )


def sub_matrix(left, right):
    return tuple(
        tuple(left[i][j] - right[i][j] for j in range(len(left[0])))
        for i in range(len(left))
    )


def multiply_matrix(left, right):
    return tuple(
        tuple(
            sum(left[i][k] * right[k][j] for k in range(len(right)))
            for j in range(len(right[0]))
        )
        for i in range(len(left))
    )


def transpose(matrix):
    return tuple(tuple(matrix[j][i] for j in range(len(matrix)))
                 for i in range(len(matrix[0])))


def matrix_vector(matrix, vector):
    return tuple(
        sum(matrix[i][j] * vector[j] for j in range(len(vector)))
        for i in range(len(matrix))
    )


def matrix_power(matrix, exponent):
    require(exponent >= 0, "matrix exponent must be nonnegative")
    result = IDENTITY3
    base = matrix
    remaining = exponent
    while remaining:
        if remaining & 1:
            result = multiply_matrix(result, base)
        base = multiply_matrix(base, base)
        remaining //= 2
    return result


def dot(left, right):
    return sum(x * y for x, y in zip(left, right))


def det2(left, right):
    return left[0] * right[1] - left[1] * right[0]


def sub_vector(left, right):
    return tuple(x - y for x, y in zip(left, right))


def cross(left, right):
    return (
        left[1] * right[2] - left[2] * right[1],
        left[2] * right[0] - left[0] * right[2],
        left[0] * right[1] - left[1] * right[0],
    )


def euclid_triple(m, n):
    return (m * m - n * n, 2 * m * n, m * m + n * n)


def spinor_scalar(u):
    return (u[0] * u[0] - u[1] * u[1]) ** 2 + 2


def spine_parent(t):
    require(t >= 1, "the positive U-spine starts at t=1")
    return (2 * t + 1, 2 * t * (t + 1), 2 * t * (t + 1) + 1)


def spine_q(t):
    return 4 * t * (t + 1) + 3


def factor_integer(value):
    require(value >= 1, "factorization input must be positive")
    factors = {}
    remaining = value
    while remaining % 2 == 0:
        factors[2] = factors.get(2, 0) + 1
        remaining //= 2
    divisor = 3
    while divisor * divisor <= remaining:
        while remaining % divisor == 0:
            factors[divisor] = factors.get(divisor, 0) + 1
            remaining //= divisor
        divisor += 2
    if remaining > 1:
        factors[remaining] = factors.get(remaining, 0) + 1
    return factors


def primitive_parameter_representations(c):
    """Return canonical m>n>0 with m^2+n^2=c and primitive parity."""
    rows = []
    for n in range(1, isqrt(c // 2) + 1):
        m_squared = c - n * n
        m = isqrt(m_squared)
        if m > n and m * m == m_squared and gcd(m, n) == 1 \
                and (m - n) % 2 == 1:
            rows.append((m, n))
    return tuple(rows)


def expected_primitive_count(c):
    require(c % 2 == 1, "fixed-hypotenuse count is scoped to odd c")
    if c == 1:
        return 0
    factors = factor_integer(c)
    if any(prime % 4 == 3 for prime in factors):
        return 0
    return 1 << (len(factors) - 1)


def gaussian_multiply(left, right):
    return (left[0] * right[0] - left[1] * right[1],
            left[0] * right[1] + left[1] * right[0])


def gaussian_power(value, exponent):
    result = (1, 0)
    base = value
    remaining = exponent
    while remaining:
        if remaining & 1:
            result = gaussian_multiply(result, base)
        base = gaussian_multiply(base, base)
        remaining //= 2
    return result


def split_prime_factor(prime):
    require(prime % 4 == 1, "only split Gaussian primes are accepted")
    for y in range(1, isqrt(prime) + 1):
        x_squared = prime - y * y
        x = isqrt(x_squared)
        if x > y and x * x == x_squared:
            return (x, y)
        if x == y and x * x == x_squared:
            return (x, y)
    raise RuntimeError(f"failed to split prime {prime}")


def canonical_parameters(gaussian):
    entries = sorted((abs(gaussian[0]), abs(gaussian[1])), reverse=True)
    require(entries[1] > 0, "positive primitive representation expected")
    return tuple(entries)


def complement_bits(bits):
    return tuple(1 - bit for bit in bits)


def bit_class(bits):
    opposite = complement_bits(bits)
    return min(tuple(bits), opposite)


def xor_bits(left, right):
    return tuple(x ^ y for x, y in zip(left, right))


def quotient_add(left, right):
    return bit_class(xor_bits(left, right))


def gaussian_choice_map(c):
    """Map F_2^r/<1> classes to canonical primitive representations."""
    factors = factor_integer(c)
    require(c > 1 and all(p % 4 == 1 for p in factors),
            "Gaussian choice map needs an admissible odd hypotenuse")
    primes = tuple(sorted(factors))
    prime_gaussians = tuple(split_prime_factor(p) for p in primes)
    raw = {}
    for bits in product((0, 1), repeat=len(primes)):
        value = (1, 0)
        for bit, prime, gaussian in zip(bits, primes, prime_gaussians):
            chosen = gaussian if bit == 0 else (gaussian[0], -gaussian[1])
            value = gaussian_multiply(
                value, gaussian_power(chosen, factors[prime])
            )
        raw[bits] = canonical_parameters(value)

    classes = {}
    for bits, parameters in raw.items():
        key = bit_class(bits)
        if key in classes:
            require(classes[key] == parameters,
                    f"global conjugation mismatch at c={c}, bits={bits}")
        else:
            classes[key] = parameters
        require(raw[complement_bits(bits)] == parameters,
                f"complement did not realize conjugation at c={c}")

    require(len(classes) == 1 << (len(primes) - 1),
            f"wrong Boolean quotient size at c={c}")
    require(set(classes.values()) == set(primitive_parameter_representations(c)),
            f"Gaussian choices missed a primitive representation at c={c}")
    return primes, classes


def roots_of_spine_polynomial(prime):
    return tuple(t for t in range(prime)
                 if (2 * t * t + 2 * t + 1) % prime == 0)


def extended_gcd(a, b):
    if b == 0:
        return (a, 1, 0)
    common, x1, y1 = extended_gcd(b, a % b)
    return (common, y1, x1 - (a // b) * y1)


def crt_pair(a, modulus_a, b, modulus_b):
    common, inverse_a, _ = extended_gcd(modulus_a, modulus_b)
    require(common == 1, "CRT moduli must be coprime")
    step = ((b - a) * inverse_a) % modulus_b
    modulus = modulus_a * modulus_b
    return ((a + modulus_a * step) % modulus, modulus)


def crt_spine_index(primes):
    residue = 0
    modulus = 1
    for prime in primes:
        roots = roots_of_spine_polynomial(prime)
        require(len(roots) == 2,
                f"split prime {prime} did not have two spine roots")
        residue, modulus = crt_pair(residue, modulus, roots[0], prime)
    if residue == 0:
        residue += modulus
    return residue, modulus


def build_bounded_berggren_tree(limit):
    root = (3, 4, 5)
    addresses = {root: ""}
    frontier = [root]
    cursor = 0
    while cursor < len(frontier):
        parent = frontier[cursor]
        cursor += 1
        for letter, matrix in (("U", U), ("A", A), ("D", D)):
            child = matrix_vector(matrix, parent)
            require(child[2] > parent[2],
                    f"hypotenuse failed to increase at {parent}->{child}")
            if child[2] <= limit:
                require(child not in addresses,
                        f"bounded Berggren tree duplicate at {child}")
                addresses[child] = addresses[parent] + letter
                frontier.append(child)
    return addresses


def tournament_orientation(points, edge_bits):
    edges = tuple(combinations(sorted(points), 2))
    oriented = {}
    for edge, bit in zip(edges, edge_bits):
        oriented[edge] = edge if bit == 0 else (edge[1], edge[0])
    return oriented


def orientation_value(oriented, left, right):
    edge = tuple(sorted((left, right)))
    return oriented[edge] == (left, right)


def translation_invariant_tournament_count(points):
    edges = tuple(combinations(sorted(points), 2))
    invariant = 0
    for edge_bits in product((0, 1), repeat=len(edges)):
        oriented = tournament_orientation(points, edge_bits)
        lawful = True
        for shift in points:
            for left, right in edges:
                translated_left = quotient_add(left, shift)
                translated_right = quotient_add(right, shift)
                if orientation_value(oriented, left, right) != orientation_value(
                        oriented, translated_left, translated_right):
                    lawful = False
                    break
            if not lawful:
                break
        if lawful:
            invariant += 1
    return invariant


def fibonacci_numbers(count):
    values = [0, 1]
    while len(values) < count:
        values.append(values[-1] + values[-2])
    return tuple(values[:count])


def main():
    # Matrix mechanism and exact U-spine.
    nilpotent = sub_matrix(U, IDENTITY3)
    nilpotent_squared = multiply_matrix(nilpotent, nilpotent)
    nilpotent_cubed = multiply_matrix(nilpotent_squared, nilpotent)
    require(nilpotent_squared != ((0, 0, 0),) * 3,
            "U-I unexpectedly had nilpotence index below three")
    require(nilpotent_cubed == ((0, 0, 0),) * 3,
            "(U-I)^3 did not vanish")
    for matrix in (U, A, D):
        require(multiply_matrix(transpose(matrix),
                                multiply_matrix(LORENTZ, matrix)) == LORENTZ,
                "Berggren matrix failed the Lorentz-form check")

    root = (3, 4, 5)
    cusp_spinor = (1, 1)
    cusp_null = euclid_triple(*cusp_spinor)
    require(cusp_null == (0, 2, 2), "raw fixed cusp had wrong content")
    require(matrix_vector(PARAMETER_U, cusp_spinor) == cusp_spinor,
            "parameter U-step did not fix the cusp spinor")
    require(matrix_vector(U, cusp_null) == cusp_null,
            "triple U-step did not fix the raw cusp")
    require(spinor_scalar(cusp_spinor) == 2,
            "fixed-cusp scalar did not produce the exceptional two")
    for depth in range(0, 501):
        t = depth + 1
        expected = spine_parent(t)
        spinor = (depth + 2, depth + 1)
        next_spinor = (depth + 3, depth + 2)
        require(matrix_vector(matrix_power(U, depth), root) == expected,
                f"U-spine power formula failed at depth {depth}")
        require(euclid_triple(*spinor) == expected,
                f"spinor lift failed at depth {depth}")
        require(matrix_vector(PARAMETER_U, spinor) == next_spinor,
                f"parameter U-step failed at depth {depth}")
        require(matrix_vector(U, euclid_triple(*spinor))
                == euclid_triple(*next_spinor),
                f"symmetric-square intertwining failed at depth {depth}")
        require(tuple(spinor[i] + cusp_spinor[i] for i in range(2))
                == next_spinor, f"cusp-fan addition failed at depth {depth}")
        require(det2(cusp_spinor, spinor) == -1
                and det2(spinor, next_spinor) == 1
                and det2(next_spinor, cusp_spinor) == 1,
                f"Farey-face determinants failed at depth {depth}")
        a, b, c = expected
        q = spine_q(t)
        require((a * a + 2, 2 * c + 1, 2 * b + 3)
                == (q, q, q), f"Q identities failed at t={t}")
        require(c - b == 1 and a * a == b + c,
                f"invariant parabola failed at t={t}")
        require(2 * b + 3 == spine_parent(b + 1)[0],
                f"self-index identity failed at t={t}")
        require(spinor_scalar(spinor) == q,
                f"fixed-cusp scalar scheme failed at depth {depth}")

    boundary = (1, 0, 1)
    fixed_null = (0, 1, 1)
    require(matrix_vector(U, boundary) == root,
            "degenerate boundary predecessor failed")
    require(matrix_vector(U, fixed_null) == fixed_null,
            "primitive null direction was not fixed")
    require(dot(fixed_null, fixed_null) == 2,
            "fixed primitive generator did not have Euclidean norm squared two")
    require(dot(matrix_vector(U, (1, 0, 0)),
                matrix_vector(U, (1, 0, 0))) != 1,
            "Euclidean norm hostile unexpectedly became U-invariant")
    require(isqrt(13 - 2) ** 2 != 13 - 2,
            "untyped odd Q hostile unexpectedly reconstructed a spine parent")
    image_depths = (spine_parent(1)[1], spine_parent(2)[1])
    image_spinors = tuple((depth + 2, depth + 1) for depth in image_depths)
    require(det2(*image_spinors) == 8,
            "self-index hostile did not destroy Farey adjacency")

    # Full descendant-angle classification and plane/area identities.
    angle_counts = {"acute": 0, "U-obtuse": 0, "D-obtuse": 0}
    euclid_rows = 0
    for m in range(2, EUCLID_BOUND + 1):
        for n in range(1, m):
            if gcd(m, n) != 1 or (m - n) % 2 == 0:
                continue
            euclid_rows += 1
            parent = euclid_triple(m, n)
            a, b, c = parent
            child_u = matrix_vector(U, parent)
            child_a = matrix_vector(A, parent)
            child_d = matrix_vector(D, parent)
            for child in (child_u, child_a, child_d):
                require(2 * child[0] + 2 * child[1] - 3 * child[2] + c == 0,
                        f"descendant plane failed at parent {parent}")

            vector_u = sub_vector(child_a, child_u)
            vector_v = sub_vector(child_d, child_u)
            vector_w = sub_vector(child_d, child_a)
            require(vector_u == (4 * b, 2 * b, 4 * b),
                    "published U-to-A vector formula failed")
            require(vector_w == (-2 * a, -4 * a, -4 * a),
                    "published A-to-D vector formula failed")
            require(dot(vector_u, vector_u) == 36 * b * b
                    and dot(vector_w, vector_w) == 36 * a * a,
                    "descendant side lengths failed")
            area_cross = cross(vector_u, vector_v)
            require(dot(area_cross, area_cross) == 272 * a * a * b * b,
                    "descendant area formula failed")

            dot_at_a = dot(tuple(-x for x in vector_u), vector_w)
            dot_at_u = dot(vector_u, vector_v)
            dot_at_d = dot(vector_v, vector_w)
            require(dot_at_a == 32 * a * b,
                    "constant-angle numerator failed")
            require(dot_at_u == 4 * b * (9 * b - 8 * a),
                    "U-angle sign formula failed")
            require(dot_at_d == 4 * a * (9 * a - 8 * b),
                    "D-angle sign formula failed")
            require(dot_at_u != 0 and dot_at_d != 0,
                    "right descendant triangle entered primitive universe")
            if dot_at_u < 0:
                require(dot_at_d > 0, "two obtuse angles appeared")
                angle_counts["U-obtuse"] += 1
            elif dot_at_d < 0:
                angle_counts["D-obtuse"] += 1
            else:
                angle_counts["acute"] += 1

    for t in range(1, 2_001):
        a, b, _ = spine_parent(t)
        require(9 * a - 8 * b < 0 and 9 * b - 8 * a > 0,
                f"U-spine obtuse classification failed at t={t}")

    # Fixed-hypotenuse count, Gaussian torsor, and split-prime boundary.
    fixed_rows = 0
    gaussian_rows = 0
    for c in range(1, FIXED_HYPOTENUSE_BOUND + 1, 2):
        representations = primitive_parameter_representations(c)
        require(len(representations) == expected_primitive_count(c),
                f"fixed-hypotenuse count failed at c={c}")
        fixed_rows += 1
        if representations:
            primes, classes = gaussian_choice_map(c)
            require(len(classes) == len(representations),
                    f"Gaussian torsor cardinality failed at c={c}")
            require(all(prime % 4 == 1 for prime in primes),
                    f"nonsplit prime entered primitive fibre at c={c}")
            gaussian_rows += 1

    # Record collision grades on the spine and audit absence of 3 mod 4 primes.
    record_grades = []
    maximum_omega = 0
    first_two_parent = None
    first_four_parent = None
    for t in range(1, SPINE_RECORD_BOUND + 1):
        c = t * t + (t + 1) * (t + 1)
        factors = factor_integer(c)
        require(all(prime % 4 == 1 for prime in factors),
                f"3 mod 4 prime divided the U-spine at t={t}")
        omega = len(factors)
        if omega > maximum_omega:
            maximum_omega = omega
            record_grades.append((omega, t, c, tuple(sorted(factors))))
        if omega >= 2 and first_two_parent is None:
            first_two_parent = (t, c)
        if omega >= 3 and first_four_parent is None:
            first_four_parent = (t, c)
    require(first_two_parent == (6, 85),
            f"wrong first two-parent spine fibre: {first_two_parent}")
    require(first_four_parent == (23, 1105),
            f"wrong first four-parent spine fibre: {first_four_parent}")
    require(record_grades[:6] == [
        (1, 1, 5, (5,)),
        (2, 6, 85, (5, 17)),
        (3, 23, 1105, (5, 13, 17)),
        (4, 223, 99905, (5, 13, 29, 53)),
        (5, 1081, 2339285, (5, 13, 17, 29, 73)),
        (6, 12131, 294346585, (5, 13, 17, 41, 73, 89)),
    ], f"unexpected Boolean record ladder: {record_grades[:6]}")
    require(set(euclid_triple(*row)
                for row in primitive_parameter_representations(65)) == {
                    (63, 16, 65), (33, 56, 65)},
            "minimal global plane-label collision at c=65 failed")

    # CRT controls for the discriminant -4 split-prime theorem.
    crt_controls = (
        (5,),
        (5, 13),
        (5, 13, 17),
        (5, 13, 17, 29),
        (5, 13, 17, 29, 37),
    )
    for primes in crt_controls:
        for prime in primes:
            require(len(roots_of_spine_polynomial(prime)) == 2,
                    f"prime {prime} failed the discriminant -4 root count")
        residue, modulus = crt_spine_index(primes)
        for lift in range(4):
            t = residue + lift * modulus
            c = t * t + (t + 1) * (t + 1)
            require(all(c % prime == 0 for prime in primes),
                    f"CRT spine construction failed for {primes}")

    for prime in (3, 7, 11, 19, 23, 31, 43, 47):
        require(prime % 4 == 3 and not roots_of_spine_polynomial(prime),
                f"inert-prime hostile failed at p={prime}")

    # Explicit ancestry/Boolean branch transplant at c=85 and c=1105.
    addresses = build_bounded_berggren_tree(1105)
    fibres = {}
    for c in (85, 1105):
        fibres[c] = tuple(
            sorted((triple, addresses[triple])
                   for triple in addresses if triple[2] == c)
        )
        require(len(fibres[c]) == expected_primitive_count(c),
                f"bounded tree missed a parent at c={c}")
    require(dict(fibres[85]) == {
        (13, 84, 85): "UUUUU",
        (77, 36, 85): "AD",
    }, f"unexpected c=85 ancestry fibre: {fibres[85]}")
    require(dict(fibres[1105]) == {
        (47, 1104, 1105): "U" * 22,
        (817, 744, 1105): "UDUA",
        (943, 576, 1105): "DAUD",
        (1073, 264, 1105): "DADDD",
    }, f"unexpected c=1105 ancestry fibre: {fibres[1105]}")

    primes_1105, choice_1105 = gaussian_choice_map(1105)
    require(primes_1105 == (5, 13, 17), "wrong split primes at 1105")
    point_to_parent = {
        point: euclid_triple(*parameters)
        for point, parameters in choice_1105.items()
    }
    points = tuple(sorted(point_to_parent))
    require(len(points) == 4, "1105 fibre was not an affine V4 torsor")
    direction_matchings = {}
    for coordinate, prime in enumerate(primes_1105):
        basis = tuple(1 if index == coordinate else 0
                      for index in range(len(primes_1105)))
        direction = bit_class(basis)
        matching = set()
        for point in points:
            partner = quotient_add(point, direction)
            matching.add(tuple(sorted((point_to_parent[point],
                                       point_to_parent[partner]))))
        require(len(matching) == 2,
                f"prime direction {prime} did not form a perfect matching")
        direction_matchings[prime] = tuple(sorted(matching))
    all_matching_edges = set().union(
        *(set(matching) for matching in direction_matchings.values())
    )
    require(len(all_matching_edges) == 6,
            "three prime matchings did not partition K4")
    require(translation_invariant_tournament_count(points) == 0,
            "an affine-V4-invariant tournament unexpectedly existed")

    # Internal generalized-Fibonacci box versus the distinct golden path.
    for t in range(1, 2_001):
        box = (1, t, t + 1, 2 * t + 1)
        a, b, c = spine_parent(t)
        require(box[0] + box[1] == box[2]
                and box[1] + box[2] == box[3],
                f"generalized Fibonacci box failed at t={t}")
        require(box[0] * box[3] == a
                and 2 * box[1] * box[2] == b,
                f"first K4 matching failed at t={t}")
        require(box[0] * box[2] + box[1] * box[3] == c,
                f"crossing K4 matching failed at t={t}")
        require(box[2] * box[3] - box[0] * box[1] == c,
                f"adjacent K4 matching failed at t={t}")

    fibonacci = fibonacci_numbers(20)
    unit_gap_indices = tuple(
        index for index in range(2, 15)
        if fibonacci[index + 1] - fibonacci[index] == 1
    )
    require(unit_gap_indices == (2, 3),
            f"U/Fibonacci path intersection changed: {unit_gap_indices}")
    parabolic_step = PARAMETER_U
    golden_step = ((1, 1), (1, 0))
    parabolic_nilpotent = sub_matrix(parabolic_step,
                                     ((1, 0), (0, 1)))
    require(multiply_matrix(parabolic_nilpotent, parabolic_nilpotent)
            == ((0, 0), (0, 0)),
            "parameter U-step was not parabolic")
    state = (0, 1)
    parabolic_mod2 = []
    golden_mod2 = []
    for _ in range(6):
        parabolic_mod2.append(state)
        state = tuple(value % 2 for value in matrix_vector(parabolic_step, state))
    state = (0, 1)
    for _ in range(6):
        golden_mod2.append(state)
        state = tuple(value % 2 for value in matrix_vector(golden_step, state))
    require(tuple(parabolic_mod2) == ((0, 1), (1, 0)) * 3,
            f"wrong parabolic parity orbit: {parabolic_mod2}")
    require(tuple(golden_mod2) == ((0, 1), (1, 0), (1, 1)) * 2,
            f"wrong golden parity orbit: {golden_mod2}")

    print("THM-3334 exact audit")
    print("matrix: (U-I)^3=0, (U-I)^2!=0, Lorentz U/A/D=PASS")
    print("cusp fan: s=(1,1) fixed; u_(n+1)=u_n+s; every {s,u_n,u_(n+1)} is Farey")
    print("spine: U^n formula and Q=a^2+2=2c+1=2b+3 checked n=0..500")
    print("boundary: U(1,0,1)=(3,4,5); primitive fixed null=(0,1,1), Euclidean norm^2=2")
    print(f"angles: {euclid_rows} primitive rows through m={EUCLID_BOUND}; "
          f"acute/U-obtuse/D-obtuse={angle_counts['acute']}/"
          f"{angle_counts['U-obtuse']}/{angle_counts['D-obtuse']}")
    print("angles: cos(A)=8/9 universally; every audited U-spine row is D-obtuse")
    print(f"fixed hypotenuse: {fixed_rows} odd c rows through "
          f"{FIXED_HYPOTENUSE_BOUND}; {gaussian_rows} nonempty Gaussian torsors")
    print("Boolean record ladder (omega,t,c,primes):")
    for record in record_grades[:6]:
        print(f"  {record}")
    print("first spine collision: t=6, c=85, fibre=C2")
    print("first spine affine-V4 fibre: t=23, c=1105=5*13*17")
    print("global scalar hostile: c=65 (Q=131) already has two represented parents")
    print("ancestry c=85: (13,84,85)->UUUUU; (77,36,85)->AD")
    print("ancestry c=1105: (47,1104,1105)->U^22; (817,744,1105)->UDUA; "
          "(943,576,1105)->DAUD; (1073,264,1105)->DADDD")
    for prime in primes_1105:
        readable = tuple(
            tuple((left[0], right[0]) for left, right in direction_matchings[prime])
        )
        print(f"prime-XOR matching p={prime}: odd-leg pairs={readable}")
    print("affine-V4 tournament control: 0/64 translation-invariant orientations")
    print(f"CRT: {len(crt_controls)} split-prime sets and four lifts each PASS; "
          "3 mod 4 hostiles PASS")
    print("Fibonacci: U boxes=(1,t,t+1,2t+1); true golden path intersects at t=1,2 only")
    print("Fibonacci parity: U-step has 2-cycle; golden step has 3-cycle including odd/odd")
    print("scope hostiles: Q=13 does not reconstruct; Q=a_b sends one Farey edge to determinant 8")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
