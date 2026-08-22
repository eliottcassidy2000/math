#!/usr/bin/env python3
"""Exact companion audit for THM-3336.

The script checks the primitive-content cocycle for Gaussian multiplication,
its Smith/root-of-minus-one normal form, the complete weighted law on bounded
Farey faces, the opposite-vertex (diamond) exchange, radius-matrix covariance,
and the parity correction required by primitive Pythagorean composition.

Universe
--------
* primitive Gaussian multipliers in ``[-8,8]^2`` and primitive vectors in
  ``[-7,7]^2`` for raw products, radii, and primitive triples;
* all triples of primitive vectors in ``[-4,4]^2`` for the content cocycle;
* primitive multipliers in ``[-10,10]^2`` and primitive vectors in
  ``[-8,8]^2`` for the Smith/content formula;
* every oriented Farey edge with endpoints in ``[-6,6]^2``, for every
  primitive multiplier in ``[-8,8]^2``;
* exact divisor-pattern censuses at eleven admissible Gaussian norms;
* eight explicit hostile controls and four lawful saturated 13-column LRC
  decks, including the two norm-2 flips and their norm-101 amplifications.

All truth-bearing arithmetic is integral (apart from one exact ``Fraction``
check of the general non-Gaussian radius action).  Validation uses explicit
exceptions, never ``assert``, so normal and optimized Python follow the same
path and print byte-identical output.
"""

from fractions import Fraction
from itertools import combinations
from math import gcd


RAW_MULTIPLIER_BOUND = 8
RAW_VECTOR_BOUND = 7
COCYCLE_BOUND = 4
SMITH_MULTIPLIER_BOUND = 10
SMITH_VECTOR_BOUND = 8
FACE_MULTIPLIER_BOUND = 8
FACE_VECTOR_BOUND = 6
PATTERN_NORMS = (1, 2, 5, 10, 13, 17, 25, 50, 65, 125, 169)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def gcd_many(values):
    answer = 0
    for value in values:
        answer = gcd(answer, abs(value))
    return answer


def lcm(left, right):
    if left == 0 or right == 0:
        return 0
    return abs(left // gcd(left, right) * right)


def lcm_many(values):
    answer = 1
    for value in values:
        answer = lcm(answer, value)
    return answer


def divisors(value):
    return tuple(d for d in range(1, value + 1) if value % d == 0)


def primitive(vector):
    return vector != (0, 0) and gcd(abs(vector[0]), abs(vector[1])) == 1


def primitive_vectors(bound):
    return tuple(
        (x, y)
        for x in range(-bound, bound + 1)
        for y in range(-bound, bound + 1)
        if primitive((x, y))
    )


def add2(left, right):
    return (left[0] + right[0], left[1] + right[1])


def sub2(left, right):
    return (left[0] - right[0], left[1] - right[1])


def scale2(scalar, vector):
    return (scalar * vector[0], scalar * vector[1])


def dot2(left, right):
    return left[0] * right[0] + left[1] * right[1]


def det2(left, right):
    return left[0] * right[1] - left[1] * right[0]


def norm2(vector):
    return dot2(vector, vector)


def gaussian_multiply(left, right):
    return (
        left[0] * right[0] - left[1] * right[1],
        left[0] * right[1] + left[1] * right[0],
    )


def gaussian_matrix(multiplier):
    a, b = multiplier
    return ((a, -b), (b, a))


def matrix_vector(matrix, vector):
    return tuple(
        sum(matrix[i][j] * vector[j] for j in range(len(vector)))
        for i in range(len(matrix))
    )


def matrix_multiply(left, right):
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


def matrix_add(left, right):
    return tuple(
        tuple(left[i][j] + right[i][j] for j in range(len(left[0])))
        for i in range(len(left))
    )


def matrix_scale(scalar, matrix):
    return tuple(tuple(scalar * entry for entry in row) for row in matrix)


def det_matrix2(matrix):
    return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]


def content2(vector):
    return gcd(abs(vector[0]), abs(vector[1]))


def multiplier_content(multiplier, vector):
    return content2(gaussian_multiply(multiplier, vector))


def primitive_image(multiplier, vector):
    product = gaussian_multiply(multiplier, vector)
    common = content2(product)
    require(common > 0, (multiplier, vector, "zero product"))
    return common, (product[0] // common, product[1] // common)


def phi(vector):
    m, n = vector
    return (m * m - n * n, 2 * m * n, m * m + n * n)


def add3(left, right):
    return tuple(x + y for x, y in zip(left, right))


def sub3(left, right):
    return tuple(x - y for x, y in zip(left, right))


def scale3(scalar, vector):
    return tuple(scalar * entry for entry in vector)


def content3(vector):
    return gcd_many(vector)


def primitive_triple(vector):
    raw = phi(vector)
    epsilon = content3(raw)
    require(epsilon in (1, 2), (vector, "unexpected Phi content", epsilon))
    return epsilon, tuple(entry // epsilon for entry in raw)


def brahmagupta(left, right):
    a, b, c = left
    x, y, z = right
    return (a * x - b * y, a * y + b * x, c * z)


def lorentz(left, right):
    return left[2] * right[2] - left[0] * right[0] - left[1] * right[1]


def det3_columns(x, y, z):
    return (
        x[0] * (y[1] * z[2] - y[2] * z[1])
        - y[0] * (x[1] * z[2] - x[2] * z[1])
        + z[0] * (x[1] * y[2] - x[2] * y[1])
    )


def radius_matrix(triple):
    a, b, c = triple
    numerators = (a + b - c, a - b + c, -a + b + c, a + b + c)
    require(all(value % 2 == 0 for value in numerators),
            (triple, "nonintegral radius matrix"))
    r, r_a, r_b, r_c = (value // 2 for value in numerators)
    return ((r, r_a), (r_b, r_c))


def signed_radius_matrix(triple):
    radius = radius_matrix(triple)
    return ((radius[0][0], -radius[0][1]),
            (-radius[1][0], radius[1][1]))


def bezout_unit(a, b):
    """Return r,t with a*r+b*t=1, for coprime signed integers a,b."""
    old_r, new_r = 1, 0
    old_t, new_t = 0, 1
    old_value, new_value = abs(a), abs(b)
    while new_value:
        quotient = old_value // new_value
        old_value, new_value = new_value, old_value - quotient * new_value
        old_r, new_r = new_r, old_r - quotient * new_r
        old_t, new_t = new_t, old_t - quotient * new_t
    require(old_value == 1, (a, b, "Bezout requested for nonprimitive pair"))
    if a < 0:
        old_r = -old_r
    if b < 0:
        old_t = -old_t
    require(a * old_r + b * old_t == 1, (a, b, old_r, old_t))
    return old_r, old_t


def smith_data(multiplier):
    a, b = multiplier
    require(primitive(multiplier), (multiplier, "Smith data needs primitive s"))
    r, t = bezout_unit(a, b)
    h = a * t - b * r
    norm = a * a + b * b
    left = ((r, t), (-b, a))
    right = ((1, -h), (0, 1))
    diagonal = matrix_multiply(matrix_multiply(left, gaussian_matrix(multiplier)),
                               right)
    require(diagonal == ((1, 0), (0, norm)),
            (multiplier, "Smith diagonal", diagonal))
    require((h * h + 1) % norm == 0, (multiplier, h, "not root of -1"))
    require((a * h + b) % norm == 0, (multiplier, h, "a*h != -b"))
    require((b * h - a) % norm == 0, (multiplier, h, "b*h != a"))
    h_shift = a * (t - a) - b * (r + b)
    require(h_shift == h - norm, (multiplier, "Bezout gauge shift"))
    return norm, h, r, t


def lattice_index(columns):
    minors = (det2(left, right) for left, right in combinations(columns, 2))
    return gcd_many(minors)


def farey_edges(bound):
    vectors = primitive_vectors(bound)
    return tuple(
        (left, right)
        for left in vectors
        for right in vectors
        if abs(det2(left, right)) == 1
    )


def audit_raw_products():
    multipliers = primitive_vectors(RAW_MULTIPLIER_BOUND)
    vectors = primitive_vectors(RAW_VECTOR_BOUND)
    raw_count = 0
    radius_count = 0
    primitive_triple_count = 0
    sign = ((1, 0), (0, -1))
    for multiplier in multipliers:
        matrix = gaussian_matrix(multiplier)
        matrix_t = transpose(matrix)
        for vector in vectors:
            product = gaussian_multiply(multiplier, vector)
            require(matrix_vector(matrix, vector) == product,
                    (multiplier, vector, "matrix/product mismatch"))
            common, image = primitive_image(multiplier, vector)
            require(primitive(image), (multiplier, vector, "image not primitive"))
            require(product == scale2(common, image),
                    (multiplier, vector, "primitive factorization"))
            require(norm2(product) == norm2(multiplier) * norm2(vector),
                    (multiplier, vector, "norm product"))
            require(phi(product) == scale3(common * common, phi(image)),
                    (multiplier, vector, "quadratic content"))
            require(brahmagupta(phi(multiplier), phi(vector)) == phi(product),
                    (multiplier, vector, "raw Brahmagupta"))
            raw_count += 1

            radius = radius_matrix(phi(vector))
            raw_radius = radius_matrix(phi(product))
            transformed = matrix_multiply(matrix_multiply(matrix, radius), matrix)
            require(raw_radius == transformed,
                    (multiplier, vector, "H covariance", raw_radius, transformed))
            require(matrix_scale(common * common, radius_matrix(phi(image)))
                    == transformed,
                    (multiplier, vector, "normalized H covariance"))
            require(det_matrix2(radius) == 0,
                    (multiplier, vector, "null radius determinant"))
            signed = signed_radius_matrix(phi(vector))
            signed_product = signed_radius_matrix(phi(product))
            require(signed_product
                    == matrix_multiply(matrix_multiply(matrix_t, signed), matrix_t),
                    (multiplier, vector, "signed V4 radius covariance"))
            # The sign matrix relation itself is independently checked here.
            require(signed == matrix_multiply(matrix_multiply(sign, radius), sign),
                    (multiplier, vector, "signed/all-positive radius relation"))
            radius_count += 1

            epsilon_s, triple_s = primitive_triple(multiplier)
            epsilon_u, triple_u = primitive_triple(vector)
            epsilon_image, triple_image = primitive_triple(image)
            numerator = common * common * epsilon_image
            denominator = epsilon_s * epsilon_u
            require(numerator % denominator == 0,
                    (multiplier, vector, "nonintegral parity factor"))
            factor = numerator // denominator
            primitive_product = brahmagupta(triple_s, triple_u)
            require(primitive_product == scale3(factor, triple_image),
                    (multiplier, vector, "primitive-triple factor"))
            require(content3(primitive_product) == factor,
                    (multiplier, vector, "wrong primitive-product content"))
            primitive_triple_count += 1
    return len(multipliers), len(vectors), raw_count, radius_count, \
        primitive_triple_count


def audit_cocycle():
    values = primitive_vectors(COCYCLE_BOUND)
    count = 0
    crt_root_count = 0
    for first in values:
        for second in values:
            product_multiplier = gaussian_multiply(first, second)
            for vector in values:
                d_second, after_second = primitive_image(second, vector)
                d_first, after_both = primitive_image(first, after_second)
                d_product, direct = primitive_image(product_multiplier, vector)
                require(d_product == d_second * d_first,
                        (first, second, vector, "content cocycle"))
                require(direct == after_both,
                        (first, second, vector, "projective monoid action"))
                count += 1
            if gcd(norm2(first), norm2(second)) == 1:
                require(primitive(product_multiplier),
                        (first, second, "coprime-norm product not primitive"))
                norm_first, h_first, _, _ = smith_data(first)
                norm_second, h_second, _, _ = smith_data(second)
                _, h_product, _, _ = smith_data(product_multiplier)
                require(h_product % norm_first == h_first % norm_first
                        and h_product % norm_second == h_second % norm_second,
                        (first, second, "CRT root gluing"))
                crt_root_count += 1
    return len(values), count, crt_root_count


def audit_smith_content():
    multipliers = primitive_vectors(SMITH_MULTIPLIER_BOUND)
    vectors = primitive_vectors(SMITH_VECTOR_BOUND)
    smith_cache = {multiplier: smith_data(multiplier)
                   for multiplier in multipliers}
    count = 0
    divisor_witness_count = 0
    for multiplier in multipliers:
        norm, h, _, _ = smith_cache[multiplier]
        for vector in vectors:
            common = multiplier_content(multiplier, vector)
            linear_content = gcd(abs(vector[0] + h * vector[1]), norm)
            require(common == linear_content,
                    (multiplier, vector, h, common, linear_content))
            require(norm % common == 0,
                    (multiplier, vector, "content does not divide norm"))
            count += 1
        for divisor in divisors(norm):
            witness = (divisor - h, 1)
            require(primitive(witness)
                    and multiplier_content(multiplier, witness) == divisor,
                    (multiplier, divisor, "divisor occurrence witness"))
            divisor_witness_count += 1
    return len(multipliers), len(vectors), count, divisor_witness_count


def audit_face(multiplier, u, v):
    delta = det2(u, v)
    require(abs(delta) == 1, (u, v, "not a Farey edge"))
    w = add2(u, v)
    d_u, image_u = primitive_image(multiplier, u)
    d_v, image_v = primitive_image(multiplier, v)
    d_w, image_w = primitive_image(multiplier, w)
    norm = norm2(multiplier)
    require(gcd(d_u, d_v) == gcd(d_v, d_w) == gcd(d_w, d_u) == 1,
            (multiplier, u, v, "face contents not pairwise coprime"))
    product = d_u * d_v * d_w
    require(norm % product == 0, (multiplier, u, v, "face product !| N"))
    kappa = norm // product

    minors = (
        det2(image_v, image_w),
        det2(image_w, image_u),
        det2(image_u, image_v),
    )
    expected = (-delta * kappa * d_u,
                -delta * kappa * d_v,
                delta * kappa * d_w)
    require(minors == expected,
            (multiplier, u, v, "signed opposite labels", minors, expected))
    labels = tuple(abs(value) for value in minors)
    require(labels == (kappa * d_u, kappa * d_v, kappa * d_w),
            (multiplier, u, v, "absolute opposite labels"))
    require(gcd_many(labels) == kappa,
            (multiplier, u, v, "kappa/gcd reconstruction"))
    require(lcm_many(labels) == norm,
            (multiplier, u, v, "norm/lcm reconstruction"))
    require(labels[0] * labels[1] * labels[2] == kappa * kappa * norm,
            (multiplier, u, v, "label product"))
    require(lattice_index((image_u, image_v, image_w)) == kappa,
            (multiplier, u, v, "image lattice index"))
    require(scale2(d_w, image_w)
            == add2(scale2(d_u, image_u), scale2(d_v, image_v)),
            (multiplier, u, v, "weighted mediant"))

    a, b = multiplier
    adjugate = ((a, b), (-b, a))
    for source, common, image in (
            (u, d_u, image_u), (v, d_v, image_v), (w, d_w, image_w)):
        numerator = matrix_vector(adjugate, scale2(common, image))
        require(all(entry % norm == 0 for entry in numerator),
                (multiplier, source, "inverse divisibility"))
        recovered = tuple(entry // norm for entry in numerator)
        require(recovered == source,
                (multiplier, source, "fixed-multiplier reconstruction"))

    raw_volume = det3_columns(phi(u), phi(v), phi(w))
    image_volume = det3_columns(phi(image_u), phi(image_v), phi(image_w))
    require(raw_volume == -4 * delta,
            (multiplier, u, v, "raw Phi volume", raw_volume))
    require(image_volume == -4 * delta * kappa * kappa * norm,
            (multiplier, u, v, "content-curved Phi volume", image_volume))

    if norm % 2 == 0:
        require(sum(common % 2 == 0 for common in (d_u, d_v, d_w)) == 1,
                (multiplier, u, v, "even-norm face parity"))
        require(kappa % 2 == 1,
                (multiplier, u, v, "even norm left 2 in kappa"))

    for left, right, label in (
            (image_v, image_w, labels[0]),
            (image_w, image_u, labels[1]),
            (image_u, image_v, labels[2])):
        require(lorentz(phi(left), phi(right)) == 2 * label * label,
                (multiplier, u, v, "Lorentz edge label"))
        h_left = radius_matrix(phi(left))
        h_right = radius_matrix(phi(right))
        cross_determinant = (det_matrix2(matrix_add(h_left, h_right))
                             - det_matrix2(h_left) - det_matrix2(h_right))
        require(cross_determinant == -2 * label * label,
                (multiplier, u, v, "radius edge label"))

    # Opposite vertices of the Farey diamond.
    minus = sub2(u, v)
    d_plus, image_plus = d_w, image_w
    d_minus, image_minus = primitive_image(multiplier, minus)
    require(add2(scale2(d_plus, image_plus), scale2(d_minus, image_minus))
            == scale2(2 * d_u, image_u),
            (multiplier, u, v, "diamond sum"))
    require(sub2(scale2(d_plus, image_plus), scale2(d_minus, image_minus))
            == scale2(2 * d_v, image_v),
            (multiplier, u, v, "diamond difference"))
    phi_left = add3(scale3(d_plus * d_plus, phi(image_plus)),
                    scale3(d_minus * d_minus, phi(image_minus)))
    phi_right = add3(scale3(2 * d_u * d_u, phi(image_u)),
                     scale3(2 * d_v * d_v, phi(image_v)))
    require(phi_left == phi_right,
            (multiplier, u, v, "symmetric-square diamond"))
    return 3


def audit_faces():
    multipliers = primitive_vectors(FACE_MULTIPLIER_BOUND)
    edges = farey_edges(FACE_VECTOR_BOUND)
    face_count = 0
    edge_label_count = 0
    for multiplier in multipliers:
        for u, v in edges:
            edge_label_count += audit_face(multiplier, u, v)
            face_count += 1
    return len(multipliers), len(edges), face_count, edge_label_count


def expected_content_patterns(norm):
    rows = set()
    divs = divisors(norm)
    for d_u in divs:
        for d_v in divs:
            for d_w in divs:
                if gcd(d_u, d_v) != 1 or gcd(d_v, d_w) != 1 \
                        or gcd(d_w, d_u) != 1:
                    continue
                if norm % 2 == 0 \
                        and sum(value % 2 == 0 for value in (d_u, d_v, d_w)) != 1:
                    continue
                rows.add((d_u, d_v, d_w))
    return rows


def audit_pattern_converse():
    total_patterns = 0
    residue_rows = 0
    for norm in PATTERN_NORMS:
        observed = set()
        for alpha in range(norm):
            for beta in range(norm):
                if gcd(gcd(alpha, beta), norm) != 1:
                    continue
                observed.add((gcd(alpha, norm), gcd(beta, norm),
                              gcd(alpha + beta, norm)))
                residue_rows += 1
        expected = expected_content_patterns(norm)
        require(observed == expected,
                (norm, "content-pattern converse", observed ^ expected))
        total_patterns += len(observed)
    return len(PATTERN_NORMS), residue_rows, total_patterns


def add_fraction_candidate(candidates, numerator, denominator):
    require(denominator > 0, "candidate denominator must be positive")
    numerator %= denominator
    if 2 * numerator > denominator:
        numerator = denominator - numerator
    common = gcd(numerator, denominator)
    candidates.add((numerator // common, denominator // common))


def exact_lonely_margin(speeds):
    """Exact maximin from all triangle-wave kinks and pair intersections."""
    speeds = tuple(sorted(set(abs(speed) for speed in speeds)))
    require(speeds and all(speed > 0 for speed in speeds), "positive speeds needed")
    candidates = {(0, 1), (1, 2)}
    for speed in speeds:
        for index in range(speed):
            add_fraction_candidate(candidates, 2 * index + 1, 2 * speed)
    for left, right in combinations(speeds, 2):
        for denominator in (left + right, abs(left - right)):
            if denominator == 0:
                continue
            for numerator in range(denominator // 2 + 1):
                add_fraction_candidate(candidates, numerator, denominator)

    best_numerator, best_denominator = -1, 1
    best_time = (0, 1)
    for numerator, denominator in candidates:
        minimum = denominator
        for speed in speeds:
            residue = (speed * numerator) % denominator
            distance = min(residue, denominator - residue)
            if distance < minimum:
                minimum = distance
        if minimum * best_denominator > best_numerator * denominator:
            best_numerator, best_denominator = minimum, denominator
            best_time = (numerator, denominator)
    return Fraction(best_numerator, best_denominator), Fraction(*best_time), \
        len(candidates)


def make_affine_bank(direction, base, step_multiplier=1):
    step = scale2(step_multiplier, direction)
    return (direction,) + tuple(
        add2(base, scale2(index, step)) for index in range(12)
    )


def audit_lrc_deck(name, multiplier, direction, columns, expected_contents,
                   expected_q, expected_q_image, expected_gate,
                   expected_gate_image, expected_margin, expected_margin_image,
                   witness_covector=None, image_witness_covector=None):
    require(columns[0] == direction and len(columns) == 13,
            (name, "malformed bank"))
    require(all(primitive(column) for column in columns),
            (name, "nonprimitive source column"))
    contents_images = tuple(primitive_image(multiplier, column)
                            for column in columns)
    contents = tuple(row[0] for row in contents_images)
    images = tuple(row[1] for row in contents_images)
    image_direction = images[0]
    require(contents == expected_contents,
            (name, "content pattern", contents, expected_contents))
    require(all(primitive(column) for column in images),
            (name, "nonprimitive image column"))
    require(lattice_index(columns) == lattice_index(images) == 1,
            (name, "source/image bank not saturated"))
    require(len(set(columns)) == len(set(images)) == 13,
            (name, "bank columns not distinct"))

    speeds = tuple(dot2(direction, column) for column in columns)
    image_speeds = tuple(dot2(image_direction, column) for column in images)
    require(all(speed > 0 for speed in speeds + image_speeds),
            (name, "target row is not positive", speeds, image_speeds))
    require(len(set(speeds)) == len(set(image_speeds)) == 13,
            (name, "target speeds not distinct"))
    for label, covector, bank in (
            ("source", witness_covector, columns),
            ("image", image_witness_covector, images)):
        if covector is None:
            continue
        witness_speeds = tuple(dot2(covector, column) for column in bank)
        require(primitive(covector) and all(speed > 0 for speed in witness_speeds)
                and len(set(witness_speeds)) == 13,
                (name, label, "claimed positive covector", witness_speeds))
    defect = max(abs(det2(direction, column)) for column in columns)
    image_defect = max(abs(det2(image_direction, column)) for column in images)
    q = norm2(direction)
    q_image = norm2(image_direction)
    require((q, q_image, defect, image_defect)
            == (expected_q, expected_q_image, 1, 1),
            (name, "q/D data", q, q_image, defect, image_defect))
    require((q >= 91 * defect) == expected_gate,
            (name, "source gate orientation"))
    require((q_image >= 91 * image_defect) == expected_gate_image,
            (name, "image gate orientation"))

    margin, time, candidates = exact_lonely_margin(speeds)
    image_margin, image_time, image_candidates = exact_lonely_margin(image_speeds)
    require(margin == expected_margin,
            (name, "source exact margin", margin, time))
    require(image_margin == expected_margin_image,
            (name, "image exact margin", image_margin, image_time))
    return {
        "name": name,
        "q": q,
        "q_image": q_image,
        "margin": margin,
        "margin_image": image_margin,
        "candidates": candidates,
        "image_candidates": image_candidates,
    }


def audit_lrc_decks():
    norm2_multiplier = (1, 1)
    safe_direction = (9, 7)
    safe_bank = make_affine_bank(safe_direction, (5, 4))
    fail_direction = (8, 1)
    fail_bank = make_affine_bank(fail_direction, (7, 1), step_multiplier=2)

    norm101_multiplier = (10, 1)
    amplified_safe_direction = (10, -1)
    amplified_safe_bank = make_affine_bank(amplified_safe_direction, (1, 0))
    amplified_fail_direction = (1, 0)
    amplified_fail_bank = make_affine_bank(
        amplified_fail_direction, (10, -1), step_multiplier=101
    )

    return (
        audit_lrc_deck(
            "N2-pass-to-fail", norm2_multiplier, safe_direction, safe_bank,
            (2,) + (1,) * 12, 130, 65, True, False,
            Fraction(1, 4), Fraction(1, 2),
            witness_covector=(-38, 49), image_witness_covector=(-50, 7),
        ),
        audit_lrc_deck(
            "N2-fail-to-pass", norm2_multiplier, fail_direction, fail_bank,
            (1,) + (2,) * 12, 65, 130, False, True,
            Fraction(1, 2), Fraction(1, 4),
            witness_covector=(-6, 49), image_witness_covector=(-50, 39),
        ),
        audit_lrc_deck(
            "N101-pass-to-fail", norm101_multiplier,
            amplified_safe_direction, amplified_safe_bank,
            (101,) + (1,) * 12, 101, 1, True, False,
            Fraction(280, 1131), Fraction(505, 1131),
        ),
        audit_lrc_deck(
            "N101-fail-to-pass", norm101_multiplier,
            amplified_fail_direction, amplified_fail_bank,
            (1,) + (101,) * 12, 1, 101, False, True,
            Fraction(505, 1131), Fraction(280, 1131),
        ),
    )


def audit_hostiles():
    count = 0

    # Primitive Gaussian factors need not have a primitive Gaussian product.
    product = gaussian_multiply((2, 1), (4, 3))
    require(product == (5, 10) and content2(product) == 5,
            "primitive Gaussian product-collapse hostile")
    require(brahmagupta((3, 4, 5), (7, 24, 25))
            == (-75, 100, 125), "Brahmagupta collapse hostile")
    _, first_root, _, _ = smith_data((2, 1))
    _, second_root, _, _ = smith_data((4, 3))
    require(first_root % 5 == 2 and second_root % 5 == 3
            and (first_root + second_root) % 5 == 0,
            "opposite-root content hostile")
    count += 1

    # The primitive-triple collapse factor is not always the spinor d^2.
    d, image = primitive_image((3, 1), (3, 1))
    epsilon_s, triple_s = primitive_triple((3, 1))
    epsilon_image, triple_image = primitive_triple(image)
    factor = d * d * epsilon_image // (epsilon_s * epsilon_s)
    require((d, factor, triple_s, triple_image) == (2, 1, (4, 3, 5), (7, 24, 25)),
            "primitive-triple parity hostile data")
    require(brahmagupta(triple_s, triple_s) == triple_image,
            "primitive-triple parity hostile product")
    count += 1

    # Multiplier primitivity is indispensable for pairwise face contents.
    bad_multiplier = (2, 0)
    bad_face = ((1, 0), (0, 1), (1, 1))
    bad_contents = tuple(multiplier_content(bad_multiplier, vertex)
                         for vertex in bad_face)
    require(bad_contents == (2, 2, 2)
            and norm2(bad_multiplier) % (2 * 2 * 2) != 0,
            "nonprimitive multiplier hostile")
    count += 1

    # Source primitivity is separately needed for d|N.
    require(content2(matrix_vector(((1, 0), (0, 1)), (2, 0))) == 2,
            "nonprimitive source hostile")
    count += 1

    # Farey adjacency is separately needed for pairwise content coprimality.
    general_matrix = ((1, 0), (0, 2))
    nonfarey_u, nonfarey_v = (0, 1), (2, 1)
    require(det2(nonfarey_u, nonfarey_v) == -2
            and content2(matrix_vector(general_matrix, nonfarey_u)) == 2
            and content2(matrix_vector(general_matrix, nonfarey_v)) == 2,
            "non-Farey hostile")
    count += 1

    # Content can be a nonunitary divisor: a prime exponent may split with kappa.
    hostile_s, hostile_u = (100, 71), (100, 59)
    hostile_norm = norm2(hostile_s)
    hostile_d = multiplier_content(hostile_s, hostile_u)
    require(hostile_norm == 15041 == 13 * 13 * 89 and hostile_d == 13
            and gcd(hostile_d, hostile_norm // hostile_d) == 13,
            "nonunitary-divisor hostile")
    count += 1

    # Omitting kappa loses a whole odd prime from an axis face.
    omission_s = (1, 2)
    omission_face = ((1, 0), (0, 1), (1, 1))
    omission_contents = tuple(multiplier_content(omission_s, vertex)
                              for vertex in omission_face)
    require(omission_contents == (1, 1, 1) and norm2(omission_s) == 5,
            "kappa-omission hostile")
    count += 1

    # G H G is special to Gaussian matrices; a shear needs two factors.
    shear = ((1, 1), (0, 1))
    vector = (2, 1)
    h_vector = radius_matrix(phi(vector))
    h_shear = radius_matrix(phi(matrix_vector(shear, vector)))
    naive = matrix_multiply(matrix_multiply(shear, h_vector), shear)
    require(naive != h_shear, "non-Gaussian radius hostile did not fire")
    l_matrix = ((Fraction(1), Fraction(-1)),
                (Fraction(1), Fraction(1)))
    l_inverse = ((Fraction(1, 2), Fraction(1, 2)),
                 (Fraction(-1, 2), Fraction(1, 2)))
    p_matrix = ((Fraction(0), Fraction(1)),
                (Fraction(1), Fraction(0)))
    shear_q = tuple(tuple(Fraction(entry) for entry in row) for row in shear)
    left_factor = matrix_multiply(matrix_multiply(l_matrix, shear_q), l_inverse)
    right_factor = transpose(matrix_multiply(matrix_multiply(p_matrix, shear_q),
                                             p_matrix))
    corrected = matrix_multiply(matrix_multiply(left_factor, h_vector), right_factor)
    require(corrected == h_shear, "general two-factor radius covariance")
    count += 1

    require(count == 8, (count, "hostile count drift"))
    return count


def main():
    multiplier_count, vector_count, raw_count, radius_count, triple_count = \
        audit_raw_products()
    cocycle_pool, cocycle_count, crt_root_count = audit_cocycle()
    smith_multipliers, smith_vectors, smith_count, divisor_witness_count = \
        audit_smith_content()
    face_multipliers, farey_edge_count, face_count, edge_label_count = \
        audit_faces()
    pattern_norms, pattern_rows, pattern_count = audit_pattern_converse()
    hostile_count = audit_hostiles()
    deck_rows = audit_lrc_decks()

    print("THM-3336 primitive Gaussian content-curvature exact audit")
    print(f"raw multipliers/vectors/checks: {multiplier_count}/{vector_count}/{raw_count}")
    print(f"content cocycle pool/checks: {cocycle_pool}/{cocycle_count}")
    print(f"coprime-norm CRT root checks: {crt_root_count}")
    print(f"Smith multipliers/vectors/content checks: "
          f"{smith_multipliers}/{smith_vectors}/{smith_count}")
    print(f"Smith divisor-occurrence witnesses: {divisor_witness_count}")
    print(f"Farey multipliers/oriented edges/faces: "
          f"{face_multipliers}/{farey_edge_count}/{face_count}")
    print(f"opposite edge-label and radius-label checks: {edge_label_count}")
    print(f"diamond and reconstruction checks: {face_count}/{3 * face_count}")
    print(f"radius covariance checks: {radius_count}")
    print(f"primitive-triple parity-factor checks: {triple_count}")
    print(f"content-pattern norms/residue rows/patterns: "
          f"{pattern_norms}/{pattern_rows}/{pattern_count}")
    print(f"hostile controls: {hostile_count}/8")
    for row in deck_rows:
        print(f"{row['name']}: q {row['q']}->{row['q_image']}; "
              f"M {row['margin']}->{row['margin_image']}; "
              f"breakpoints {row['candidates']}/{row['image_candidates']}")
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
