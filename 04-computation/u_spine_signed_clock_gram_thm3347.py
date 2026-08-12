#!/usr/bin/env python3
"""Exact audit for the U-spine signed prime-clock Gram/projective boundary.

Universe
--------
* every pair 0 <= r,s <= 300 for the unrestricted Hensel-layer Gram law;
* all positive integer weight tuples through 8 in ranks 1,2,3 for the
  low-rank conditional-negative-type boundary;
* balanced equal-weight ranks 4..10, direct Catalan level-four controls
  through rank 30, and dominant-weight controls in ranks 4..7 for the
  antipodal folded metric;
* exact projector and cosh kernels on selected prime-power cubes through
  rank seven;
* the literal arithmetic hostile N=32045=5*13*17*29, using formal prime-log
  coefficient vectors rather than floating-point logarithms.

All decisions are integral or rational.  The universal claims in THM-3347
rest on its displayed factorizations; this companion supplies exhaustive
bounded controls and sharp hostile witnesses.
"""

from fractions import Fraction
from itertools import product
from math import comb, factorial, gcd


PAIR_LIMIT = 300
COHOM_CUBES = (
    (5,),
    (5, 13),
    (5, 13, 17),
    (5, 13, 17, 29),
    (5, 13, 17, 29, 37),
    (5, 13, 17, 29, 37, 41),
    (5, 13, 17, 29, 37, 41, 53),
    (25, 13, 17, 29),
)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def factorization(n):
    result = {}
    p = 2
    while p * p <= n:
        while n % p == 0:
            result[p] = result.get(p, 0) + 1
            n //= p
        p = 3 if p == 2 else p + 2
    if n > 1:
        result[n] = result.get(n, 0) + 1
    return result


def C(t):
    return 2 * t * t + 2 * t + 1


def content_channels(r, s):
    return gcd(C(r), abs(s - r)), gcd(C(r), r + s + 1)


def clean_form(form):
    return {p: exponent for p, exponent in form.items() if exponent}


def add_log(form, n, coefficient=1):
    for p, exponent in factorization(n).items():
        form[p] = form.get(p, 0) + coefficient * exponent


def subtract_forms(left, right):
    result = dict(left)
    for p, exponent in right.items():
        result[p] = result.get(p, 0) - exponent
    return clean_form(result)


def signed_clock_features(t):
    """Return {(p,e): sign} using the smaller reflected t-root as +."""
    result = {}
    for p, valuation in factorization(C(t)).items():
        require(p % 4 == 1, "U-spine norm acquired a nonsplit odd prime")
        for exponent in range(1, valuation + 1):
            modulus = p ** exponent
            residue = t % modulus
            reflected = modulus - 1 - residue
            require(C(residue) % modulus == 0
                    and C(reflected) % modulus == 0,
                    "Hensel feature did not lie on the two root branches")
            result[(p, exponent)] = 1 if residue < reflected else -1
    return result


def feature_dot_form(features_r, features_s):
    result = {}
    for coordinate in features_r.keys() & features_s.keys():
        p, _ = coordinate
        result[p] = result.get(p, 0) \
                    + features_r[coordinate] * features_s[coordinate]
    return clean_form(result)


def expected_lambda_form(r, s):
    d_minus, d_plus = content_channels(r, s)
    result = {}
    add_log(result, d_minus, 1)
    add_log(result, d_plus, -1)
    return clean_form(result)


def popcount(n):
    return n.bit_count()


def quotient_representatives(dimension):
    """Chart the antipodal quotient by fixing coordinate zero to zero."""
    return tuple(value << 1 for value in range(1 << (dimension - 1)))


def weighted_hamming(weights, difference):
    return sum(weight for index, weight in enumerate(weights)
               if difference & (1 << index))


def folded_distance(weights, difference):
    h = weighted_hamming(weights, difference)
    return min(h, sum(weights) - h)


def quotient_walsh(values, character):
    """Unnormalized Fourier coefficient on the bit-zero quotient chart."""
    return sum(value * (-1 if popcount(character & index) % 2 else 1)
               for index, value in enumerate(values))


def even_full_character(quotient_character):
    full = quotient_character << 1
    if popcount(quotient_character) % 2:
        full |= 1
    require(popcount(full) % 2 == 0, "quotient character was not even")
    return full


def cosh_log_rational(q):
    return Fraction(q * q + 1, 2 * q)


def sinh_log_rational(q):
    return Fraction(q * q - 1, 2 * q)


def exp_lambda_ratio(prime_powers, difference):
    numerator = 1
    denominator = 1
    for index, q in enumerate(prime_powers):
        if difference & (1 << index):
            denominator *= q
        else:
            numerator *= q
    return Fraction(numerator, denominator)


def cosh_from_ratio(ratio):
    return (ratio + 1 / ratio) / 2


def catalan(index):
    return comb(2 * index, index) // (index + 1)


def equal_weight_level_four_coefficient(dimension):
    """Exact quotient coefficient via the level-four Krawtchouk sum."""
    total = 0
    for weight in range(dimension + 1):
        krawtchouk = sum(
            (-1) ** overlap * comb(4, overlap)
            * comb(dimension - 4, weight - overlap)
            for overlap in range(5)
            if 0 <= weight - overlap <= dimension - 4
        )
        total += min(weight, dimension - weight) * krawtchouk
    require(total % 2 == 0, "level-four raw coefficient was not even")
    return total // 2


def even_power_predicted_eigenvalue(weights, full_character, half_degree):
    """Return the quotient eigenvalue from the parity-constrained multinomial."""
    dimension = len(weights)
    target_parities = tuple(
        1 if full_character & (1 << index) else 0
        for index in range(dimension)
    )
    total_degree = 2 * half_degree
    coefficient = 0

    def visit(index, remaining, exponents):
        nonlocal coefficient
        if index == dimension - 1:
            exponent = remaining
            if exponent % 2 != target_parities[index]:
                return
            row = exponents + [exponent]
            term = factorial(total_degree)
            for value in row:
                term //= factorial(value)
            for weight, value in zip(weights, row):
                term *= weight ** value
            coefficient += term
            return
        for exponent in range(target_parities[index], remaining + 1, 2):
            visit(index + 1, remaining - exponent, exponents + [exponent])

    visit(0, total_degree, [])
    return (1 << (dimension - 1)) * coefficient


def roots(n):
    return tuple(t for t in range(n) if C(t) % n == 0)


def root_bits(n, root, base, prime_powers):
    x = 2 * root + 1
    base_x = 2 * base + 1
    bits = 0
    for index, q in enumerate(prime_powers):
        if x % q != base_x % q:
            require(x % q == (-base_x) % q,
                    "root was neither the base nor reflected CRT branch")
            bits |= 1 << index
    return bits


def main():
    unrestricted_pairs = 0
    distance_invoices = 0
    feature_cache = {t: signed_clock_features(t)
                     for t in range(PAIR_LIMIT + 1)}
    for t, features in feature_cache.items():
        norm_form = {}
        add_log(norm_form, C(t))
        diagonal = feature_dot_form(features, features)
        require(diagonal == norm_form,
                "signed clock feature norm did not equal log C_t")
        require(sum(1 for _ in features) == sum(factorization(C(t)).values()),
                "clock feature support did not equal total valuation")

    for r in range(PAIR_LIMIT + 1):
        for s in range(PAIR_LIMIT + 1):
            observed = feature_dot_form(feature_cache[r], feature_cache[s])
            expected = expected_lambda_form(r, s)
            require(observed == expected,
                    "Hensel-layer dot product did not equal signed content")

            d_minus, d_plus = content_channels(r, s)
            common = gcd(C(r), C(s))
            require(d_minus * d_plus == common and gcd(d_minus, d_plus) == 1,
                    "inherited two-channel splitter failed")

            # ||Psi(r)-Psi(s)||^2 = log(A B d_+^4), where
            # A=C_r/common and B=C_s/common.  Check as formal log vectors.
            distance_form = {}
            add_log(distance_form, C(r))
            add_log(distance_form, C(s))
            for p, exponent in observed.items():
                distance_form[p] = distance_form.get(p, 0) - 2 * exponent
            target_minus = {}
            add_log(target_minus, C(r) // common)
            add_log(target_minus, C(s) // common)
            add_log(target_minus, d_plus, 4)
            require(clean_form(distance_form) == clean_form(target_minus),
                    "unrestricted minus-distance invoice failed")

            plus_form = {}
            add_log(plus_form, C(r))
            add_log(plus_form, C(s))
            for p, exponent in observed.items():
                plus_form[p] = plus_form.get(p, 0) + 2 * exponent
            target_plus = {}
            add_log(target_plus, C(r) // common)
            add_log(target_plus, C(s) // common)
            add_log(target_plus, d_minus, 4)
            require(clean_form(plus_form) == clean_form(target_plus),
                    "unrestricted plus-distance invoice failed")
            unrestricted_pairs += 1
            distance_invoices += 2

    # In ranks <=3 every weighted antipodal folded metric is CND.  Exhaust
    # all small positive integer weights and verify the Fourier criterion.
    low_rank_weight_rows = 0
    low_rank_fourier_checks = 0
    for dimension in range(1, 4):
        for weights in product(range(1, 9), repeat=dimension):
            reps = quotient_representatives(dimension)
            values = [folded_distance(weights, representative)
                      for representative in reps]
            absolute_values = [sum(weights) - 2 * value for value in values]
            require(values[0] == 0, "folded quotient distance lost origin")
            require(all(value >= 0 for value in absolute_values),
                    "absolute folded Gram acquired a negative entry")
            require(quotient_walsh(absolute_values, 0) >= 0,
                    "absolute folded Gram lost its constant mode")
            for character in range(1, 1 << (dimension - 1)):
                coefficient = quotient_walsh(values, character)
                require(coefficient <= 0,
                        "rank-at-most-three folded metric was not CND")
                require(quotient_walsh(absolute_values, character)
                        == -2 * coefficient >= 0,
                        "rank-at-most-three absolute kernel was not PSD")
                low_rank_fourier_checks += 1
            low_rank_weight_rows += 1

    # Equal weights expose the first positive level-four Fourier mode in every
    # rank >=4.  The exact values are integers.
    balanced_obstructions = []
    for dimension in range(4, 11):
        weights = (1,) * dimension
        reps = quotient_representatives(dimension)
        values = [folded_distance(weights, representative)
                  for representative in reps]
        positives = []
        for character in range(1, 1 << (dimension - 1)):
            coefficient = quotient_walsh(values, character)
            full = even_full_character(character)
            if coefficient > 0:
                positives.append((popcount(full), coefficient))
        require(positives and min(level for level, _ in positives) == 4,
                "balanced fold did not first fail at Fourier level four")
        balanced_obstructions.append(
            (dimension, min(value for level, value in positives if level == 4)))

    catalan_controls = 0
    for dimension in range(4, 31):
        half = dimension // 2
        expected = (2 * catalan(half - 2) if dimension % 2 == 0
                    else catalan(half - 1))
        require(equal_weight_level_four_coefficient(dimension) == expected,
                "all-rank Catalan level-four formula failed")
        catalan_controls += 1

    # A dominant coordinate supplies a canonical hemisphere section.  In the
    # chart fixing that bit to zero, folded distance is ordinary Hamming on
    # all remaining coordinates and is therefore CND.
    dominant_rows = 0
    dominant_fourier_checks = 0
    for dimension in range(4, 8):
        for tail in product(range(1, 5), repeat=dimension - 1):
            weights = (sum(tail),) + tail
            reps = quotient_representatives(dimension)
            values = []
            for representative in reps:
                folded = folded_distance(weights, representative)
                ordinary = sum(weights[index]
                               for index in range(1, dimension)
                               if representative & (1 << index))
                require(folded == ordinary,
                        "dominant hemisphere did not unfold the metric")
                values.append(folded)
            for character in range(1, 1 << (dimension - 1)):
                require(quotient_walsh(values, character) <= 0,
                        "dominant-coordinate folded metric was not CND")
                dominant_fourier_checks += 1
            dominant_rows += 1

    projector_rows = 0
    projector_fourier_checks = 0
    raw_gram_fourier_checks = 0
    cosh_rows = 0
    cosh_fourier_checks = 0
    tensor_rows = 0
    tensor_fourier_checks = 0
    for prime_powers in COHOM_CUBES:
        dimension = len(prime_powers)
        # Use integer weights here.  The identities are polynomial in the
        # weights; fixed arithmetic logs are inserted by the theorem.
        weights = tuple(index + 2 for index in range(dimension))
        total = sum(weights)
        reps = quotient_representatives(dimension)

        raw_signed_values = [total - 2 * weighted_hamming(weights, difference)
                             for difference in range(1 << dimension)]
        for character in range(1 << dimension):
            observed = sum(
                value * (-1 if popcount(character & difference) % 2 else 1)
                for difference, value in enumerate(raw_signed_values)
            )
            expected = 0
            if popcount(character) == 1:
                index = next(index for index in range(dimension)
                             if character & (1 << index))
                expected = (1 << dimension) * weights[index]
            require(observed == expected,
                    "raw fixed-grade Gram spectrum failed")
            raw_gram_fourier_checks += 1

        squared_lambda = []
        for representative in reps:
            h = weighted_hamming(weights, representative)
            signed = total - 2 * h
            d = min(h, total - h)
            require(total * total - signed * signed == 4 * d * (total - d),
                    "projector chord transform failed")
            require(2 * total - 2 * abs(signed) == 4 * d,
                    "quotient-sphere chord transform failed")
            if representative:
                require(d > 0 and total * total - signed * signed > 0,
                        "projector failed to separate parent vertices")
            squared_lambda.append(signed * signed)
            projector_rows += 1

        quotient_size = len(reps)
        nonzero_projector_modes = 0
        for character in range(1 << (dimension - 1)):
            full = even_full_character(character)
            coefficient = quotient_walsh(squared_lambda, character)
            level = popcount(full)
            if level == 0:
                expected = quotient_size * sum(weight * weight
                                               for weight in weights)
            elif level == 2:
                indices = [index for index in range(dimension)
                           if full & (1 << index)]
                expected = quotient_size * 2 * weights[indices[0]] \
                           * weights[indices[1]]
            else:
                expected = 0
            require(coefficient == expected,
                    "squared Gram kernel left levels zero and two")
            if coefficient:
                nonzero_projector_modes += 1
            projector_fourier_checks += 1
        require(nonzero_projector_modes == 1 + dimension * (dimension - 1) // 2,
                "projector kernel rank was not 1+binomial(k,2)")

        # At tau=1, cosh(log q) and sinh(log q) are rational.  Check the exact
        # strictly-positive quotient spectrum against the product formula.
        cosh_values = []
        for representative in reps:
            ratio = exp_lambda_ratio(prime_powers, representative)
            cosh_values.append(cosh_from_ratio(ratio))
            cosh_rows += 1
        for character in range(1 << (dimension - 1)):
            full = even_full_character(character)
            expected = Fraction(quotient_size)
            for index, q in enumerate(prime_powers):
                expected *= (sinh_log_rational(q)
                             if full & (1 << index)
                             else cosh_log_rational(q))
            observed = quotient_walsh(cosh_values, character)
            require(observed == expected and observed > 0,
                    "cosh kernel spectrum was not strictly positive")
            cosh_fourier_checks += 1

        # Even tensor powers expose exactly even characters through level 2m.
        for half_degree in range(1, 5):
            tensor_values = []
            for representative in reps:
                signed = total - 2 * weighted_hamming(weights, representative)
                tensor_values.append(signed ** (2 * half_degree))
                tensor_rows += 1
            observed_rank = 0
            for character in range(1 << (dimension - 1)):
                full = even_full_character(character)
                observed = quotient_walsh(tensor_values, character)
                expected = even_power_predicted_eigenvalue(
                    weights, full, half_degree
                )
                require(observed == expected,
                        "even tensor multinomial spectrum failed")
                should_be_positive = popcount(full) <= 2 * half_degree
                require((observed > 0) == should_be_positive,
                        "even tensor support boundary failed")
                observed_rank += int(observed > 0)
                tensor_fourier_checks += 1
            expected_rank = sum(
                comb(dimension, 2 * level)
                for level in range(min(half_degree, dimension // 2) + 1)
            )
            require(observed_rank == expected_rank,
                    "even tensor kernel rank failed")

    # Literal arithmetic rank-four hostile.  Fix the 5-adic bit to zero and
    # use the all-three character on the quotient chart.
    hostile_n = 5 * 13 * 17 * 29
    hostile_prime_powers = (5, 13, 17, 29)
    hostile_roots = roots(hostile_n)
    require(len(hostile_roots) == 16 and hostile_roots[0] == 1081,
            "rank-four hostile root atlas changed")
    base = hostile_roots[0]
    bits_to_root = {
        root_bits(hostile_n, root, base, hostile_prime_powers): root
        for root in hostile_roots
    }
    require(len(bits_to_root) == 16, "rank-four CRT root bits collided")
    hostile_reps = [bits_to_root[index << 1] for index in range(8)]
    hostile_coefficients = [(-1) ** popcount(index) for index in range(8)]
    require(sum(hostile_coefficients) == 0,
            "rank-four hostile coefficient vector was not centered")

    distance_quadratic = {}
    absolute_quadratic = {}
    for i, r in enumerate(hostile_reps):
        for j, s in enumerate(hostile_reps):
            coefficient = hostile_coefficients[i] * hostile_coefficients[j]
            delta_minus = gcd(hostile_n, abs(s - r))
            delta_plus = gcd(hostile_n, r + s + 1)
            require(delta_minus * delta_plus == hostile_n,
                    "rank-four hostile channels did not complement")
            smaller = min(delta_minus, delta_plus)
            larger = max(delta_minus, delta_plus)
            add_log(distance_quadratic, smaller, coefficient)
            add_log(absolute_quadratic, larger, coefficient)
            add_log(absolute_quadratic, smaller, -coefficient)
    require(clean_form(distance_quadratic) == {5: 16},
            "rank-four folded distance hostile changed")
    require(clean_form(absolute_quadratic) == {5: -32},
            "rank-four absolute Gram hostile changed")

    print("U-SPINE SIGNED PRIME-CLOCK GRAM -- EXACT AUDIT")
    print(f"unrestricted signed-content Gram pairs: {unrestricted_pairs}")
    print(f"unrestricted Hilbert distance invoices: {distance_invoices}")
    print(f"low-rank CND weight rows/Fourier checks: "
          f"{low_rank_weight_rows}/{low_rank_fourier_checks}")
    print("balanced first-positive level-4 rows (rank,coefficient):",
          balanced_obstructions)
    print(f"Catalan level-4 formula controls, ranks 4..30: {catalan_controls}")
    print(f"dominant-section rows/Fourier checks: "
          f"{dominant_rows}/{dominant_fourier_checks}")
    print(f"raw fixed-grade Gram Fourier checks: {raw_gram_fourier_checks}")
    print(f"projector rows/Fourier checks: "
          f"{projector_rows}/{projector_fourier_checks}")
    print(f"cosh rows/strict-positive Fourier checks: "
          f"{cosh_rows}/{cosh_fourier_checks}")
    print(f"even-tensor rows/Fourier checks: "
          f"{tensor_rows}/{tensor_fourier_checks}")
    print("rank-four arithmetic hostile N=32045 roots:", hostile_reps)
    print("rank-four hostile forms: d -> 16 log(5), |Lambda| -> -32 log(5)")
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
