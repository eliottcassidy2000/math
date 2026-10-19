#!/usr/bin/env python3
"""Stdlib exact controls for the elliptic path through 482, 3362, 6242.

Checks remain active under python -O. Finite evidence does not replace the
non-torsion, positivity, distinctness, or asymptotic arguments in the note.
"""

import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
import json
from math import gcd, isqrt, lcm


MODEL = (F(-482), F(6195836), F(-299007352))
P = (F(482), F(51840))
O = None
checks = Counter()


def check(condition, name):
    checks[name] += 1
    if not condition:
        raise RuntimeError(name)


def cubic(x, model=MODEL):
    a2, a4, a6 = model
    return x ** 3 + a2 * x * x + a4 * x + a6


def curve_discriminant(model):
    a2, a4, a6 = model
    return 16 * (a2 * a2 * a4 * a4 - 4 * a4 ** 3 - 4 * a2 ** 3 * a6
                 - 27 * a6 * a6 + 18 * a2 * a4 * a6)


def on_curve(point, model=MODEL):
    return point is O or point[1] * point[1] == cubic(point[0], model)


def slope(left, right, model=MODEL):
    x, y = left
    u, v = right
    a2, a4, _ = model
    return ((3 * x * x + 2 * a2 * x + a4) / (2 * y)
            if left == right else (v - y) / (u - x))


def add(left, right, model=MODEL):
    if left is O:
        return right
    if right is O:
        return left
    x, y = left
    u, v = right
    if x == u and y == -v:
        return O
    m = slope(left, right, model)
    xx = m * m - model[0] - x - u
    return xx, m * (x - xx) - y


def multiply(n, point, model=MODEL):
    result = O
    while n:
        if n & 1:
            result = add(result, point, model)
        point = add(point, point, model)
        n >>= 1
    return result


def rational_bits(value):
    return max(abs(value.numerator).bit_length(), value.denominator.bit_length())


def poly_mul(left, right):
    result = [F(0)] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            result[i + j] += a * b
    return result


def digest_integers(values):
    digest = sha256()
    for value in values:
        data = value.to_bytes(max(1, (value.bit_length() + 7) // 8), "big")
        digest.update(len(data).to_bytes(8, "big"))
        digest.update(data)
    return digest.hexdigest()


def interpolate(A, B, C, s1, s2):
    A, B, C, s1, s2 = map(F, (A, B, C, s1, s2))
    check(A + B == s1 * s1 and B + C == s2 * s2, "interpolation_input_edges")
    check(A != B, "interpolation_distinct_first_two")
    Y = -(s1 + s2) * (B - A) / 2
    check(Y != 0, "interpolation_nonzero_y")
    a2 = -A
    a4 = 2 * Y * s1 - A * A
    a6 = Y * Y - A * a4
    model = (a2, a4, a6)
    point = (A, Y)
    check(on_curve(point, model), "interpolated_point_on_cubic")
    check(slope(point, point, model) == s1, "interpolated_tangent_slope")
    twice = add(point, point, model)
    check(twice[0] == B and on_curve(twice, model), "interpolated_second_abscissa")
    check(slope(point, twice, model) == s2, "interpolated_second_chord_slope")
    thrice = add(point, twice, model)
    check(thrice[0] == C and on_curve(thrice, model), "interpolated_third_abscissa")
    return model, point


def sp_add(left, right):
    result = dict(left)
    for degree, coefficient in right.items():
        result[degree] = result.get(degree, 0) + coefficient
    return {degree: coefficient for degree, coefficient in result.items() if coefficient}


def sp_scale(poly, factor):
    return {degree: coefficient * factor for degree, coefficient in poly.items()
            if coefficient * factor}


def sp_mul(left, right):
    result = {}
    for i, a in left.items():
        for j, b in right.items():
            result[i + j] = result.get(i + j, 0) + a * b
    return {degree: coefficient for degree, coefficient in result.items() if coefficient}


def sp_sub(left, right):
    return sp_add(left, sp_scale(right, -1))


def sp_square(poly):
    return sp_mul(poly, poly)


def sp_coefficients(poly):
    return [poly.get(i, 0) for i in range(max(poly, default=0) + 1)]


def sp_evaluate(poly, value):
    return sum(coefficient * value ** degree for degree, coefficient in poly.items())


def bp_add(left, right):
    result = dict(left)
    for powers, coefficient in right.items():
        result[powers] = result.get(powers, 0) + coefficient
    return {powers: coefficient for powers, coefficient in result.items() if coefficient}


def bp_scale(poly, scalar):
    return {powers: coefficient * scalar for powers, coefficient in poly.items()
            if coefficient * scalar}


def bp_sub(left, right):
    return bp_add(left, bp_scale(right, -1))


def bp_mul(left, right):
    result = {}
    for (i, j), a in left.items():
        for (k, ell), b in right.items():
            powers = (i + k, j + ell)
            result[powers] = result.get(powers, 0) + a * b
    return {powers: coefficient for powers, coefficient in result.items() if coefficient}


def bp_pow(poly, exponent):
    result = {(0, 0): 1}
    for _ in range(exponent):
        result = bp_mul(result, poly)
    return result


def bp_evaluate(poly, a, b):
    return sum(coefficient * a ** i * b ** j for (i, j), coefficient in poly.items())


def valuation(value, prime):
    value = F(value)
    result = 0
    numerator, denominator = abs(value.numerator), value.denominator
    if numerator == 0:
        raise ValueError("valuation of zero is not used")
    while numerator % prime == 0:
        numerator //= prime
        result += 1
    while denominator % prime == 0:
        denominator //= prime
        result -= 1
    return result


def primitive_ppt_certificate():
    # Independent exact arithmetic in the sparse bivariate ring Z[a,b].
    a, b = {(1, 0): 1}, {(0, 1): 1}
    a2, b2, ab = bp_pow(a, 2), bp_pow(b, 2), bp_mul(a, b)
    A = bp_scale(bp_sub(bp_add(a2, b2), bp_scale(ab, 4)), 2)
    Y = bp_scale(bp_mul(a2, b), 16)
    U = bp_sub(bp_scale(bp_mul(Y, bp_sub(b, a)), 4), bp_pow(A, 2))
    V = bp_sub(bp_pow(Y, 2), bp_mul(A, U))
    disc = bp_add(bp_sub(bp_mul(bp_pow(A, 2), bp_pow(U, 2)), bp_scale(bp_pow(U, 3), 4)),
                  bp_scale(bp_mul(bp_pow(A, 3), V), 4))
    disc = bp_sub(bp_sub(disc, bp_scale(bp_pow(V, 2), 27)), bp_scale(bp_mul(bp_mul(A, U), V), 18))
    factor = {(3, 0): 1, (2, 1): 6, (1, 2): 14, (0, 3): -2}
    expected_disc = bp_scale(bp_mul({(5, 4): 1}, factor), -65536)
    check(disc == expected_disc, "ppt_exact_bivariate_cubic_discriminant")
    rational_root_values = {str(t): t ** 3 + 6 * t * t + 14 * t - 2 for t in (-2, -1, 1, 2)}
    check(all(rational_root_values.values()), "ppt_cubic_factor_has_no_rational_root_candidates")
    xs = [A, bp_scale(bp_add(a2, b2), 2),
          bp_scale(bp_add(bp_add(a2, b2), bp_scale(ab, 4)), 2),
          bp_scale(bp_sub(bp_sub(b2, a2), bp_scale(ab, 4)), 2),
          bp_scale(bp_add(bp_add(a2, bp_scale(b2, 17)), bp_scale(ab, 4)), 2)]
    ys = [Y, bp_scale(bp_mul(a, b2), -16),
          bp_scale(bp_mul(ab, bp_add(a, bp_scale(b, 2))), 16),
          bp_scale(bp_mul(a2, b), -8),
          bp_scale(bp_mul(b, bp_add(bp_add(a2, bp_scale(ab, 6)), bp_scale(b2, 12))), -16)]
    ms = [bp_scale(bp_sub(b, a), 2), bp_scale(bp_add(a, b), -2),
          bp_scale(b, 2), bp_scale(b, 6)]
    for i, (x, y) in enumerate(zip(xs, ys)):
        rhs = bp_add(bp_add(bp_sub(bp_pow(x, 3), bp_mul(A, bp_pow(x, 2))), bp_mul(U, x)), V)
        check(bp_pow(y, 2) == rhs, "ppt_first_five_symbolic_curve_points")
        if i < 4:
            check(bp_add(x, xs[i + 1]) == bp_pow(ms[i], 2), "ppt_first_four_symbolic_edges")
            check(ys[i + 1] == bp_sub(bp_mul(ms[i], bp_sub(x, xs[i + 1])), y),
                  "ppt_first_four_symbolic_y_updates")
            if i == 0:
                check(bp_scale(bp_mul(Y, ms[i]), 2) == bp_add(bp_pow(A, 2), U),
                      "ppt_symbolic_tangent")
            else:
                check(bp_sub(y, Y) == bp_mul(ms[i], bp_sub(x, A)), "ppt_symbolic_chord_incidence")
    d = bp_add(a, bp_scale(b, 2))
    N = {(4, 0): 2, (3, 1): 8, (2, 2): 10, (1, 3): -24, (0, 4): 8}
    M = bp_scale(bp_mul({(2, 2): 1}, {(2, 0): 1, (1, 1): 6, (0, 2): 16}), 16)
    K = {(2, 0): -2, (1, 1): -6, (0, 2): -12}
    check(bp_mul(bp_sub(ys[4], Y), d) == bp_mul(K, bp_sub(xs[4], A)),
          "ppt_fifth_symbolic_chord_incidence")
    check(bp_add(bp_mul(xs[4], bp_pow(d, 2)), N) == bp_pow(K, 2), "ppt_sixth_symbolic_abscissa")
    check(M == bp_sub(bp_mul(K, bp_sub(bp_mul(xs[4], bp_pow(d, 2)), N)), bp_mul(ys[4], bp_pow(d, 3))),
          "ppt_sixth_symbolic_ordinate")
    rhs6 = bp_sub(bp_pow(N, 3), bp_mul(bp_mul(A, bp_pow(N, 2)), bp_pow(d, 2)))
    rhs6 = bp_add(bp_add(rhs6, bp_mul(bp_mul(U, N), bp_pow(d, 4))), bp_mul(V, bp_pow(d, 6)))
    check(bp_pow(M, 2) == rhs6, "ppt_sixth_symbolic_curve_point")
    remainder = {}
    for (i, j), coefficient in N.items():
        remainder[i + j] = remainder.get(i + j, 0) + coefficient * (-2) ** i
    remainder = {power: coefficient for power, coefficient in remainder.items() if coefficient}
    check(remainder == {4: 64}, "ppt_numerator_remainder_mod_a_plus_two_b")
    short_a4 = bp_sub(bp_scale(U, 81), bp_scale(bp_pow(A, 2), 27))
    short_a6 = bp_sub(bp_add(bp_scale(V, 729), bp_scale(bp_mul(A, U), 243)), bp_scale(bp_pow(A, 3), 54))
    short_N = bp_sub(bp_scale(N, 9), bp_scale(bp_mul(A, bp_pow(d, 2)), 3))
    short_M = bp_scale(M, 27)
    short_rhs = bp_add(bp_add(bp_pow(short_N, 3), bp_mul(bp_mul(short_a4, short_N), bp_pow(d, 4))),
                       bp_mul(short_a6, bp_pow(d, 6)))
    check(bp_pow(short_M, 2) == short_rhs, "ppt_sixth_integral_short_model_identity")
    # Prove the finite residue implication behind 3 not dividing a+2b.
    for ar in range(3):
        for br in range(3):
            for cr in range(3):
                if (ar, br) != (0, 0) and (ar * ar + br * br - cr * cr) % 3 == 0:
                    check((ar + 2 * br) % 3 != 0, "ppt_mod_three_denominator_exclusion")
    examples = []
    count = 0
    for m in range(2, 31):
        for n in range(1, m):
            if gcd(m, n) != 1 or (m - n) % 2 == 0:
                continue
            leg1, leg2 = m * m - n * n, 2 * m * n
            for aa, bb in ((leg1, leg2), (leg2, leg1)):
                count += 1
                dd = aa + 2 * bb
                check(gcd(aa, bb) == 1 and gcd(dd, bb) == 1, "ppt_sample_coprimality")
                check(dd > 6 and valuation(dd, 2) <= 1 and dd % 3 != 0,
                      "ppt_sample_denominator_prime_conditions")
                av, yv, uv, vv = [bp_evaluate(poly, aa, bb) for poly in (A, Y, U, V)]
                model = tuple(map(F, (-av, uv, vv)))
                point = (F(av), F(yv))
                Q = O
                for j in range(1, 7):
                    Q = add(Q, point, model)
                    expected = ((F(bp_evaluate(xs[j - 1], aa, bb)), F(bp_evaluate(ys[j - 1], aa, bb)))
                                if j <= 5 else
                                (F(bp_evaluate(N, aa, bb), dd * dd), F(bp_evaluate(M, aa, bb), dd ** 3)))
                    check(Q == expected and on_curve(Q, model), "ppt_sample_first_six_points")
                exact_denominator = (dd // gcd(dd, 2)) ** 2
                X6 = 9 * Q[0] - 3 * av
                check(Q[0].denominator == exact_denominator and X6.denominator == exact_denominator,
                      "ppt_sample_exact_sixth_denominator")
                check(exact_denominator > 1 and Q[1] > 0, "ppt_sample_lutz_nontorsion_certificate")
                if len(examples) < 8:
                    examples.append({"a": aa, "b": bb, "d": dd, "x6": str(Q[0]),
                                     "short_X6": str(X6), "exact_denominator": exact_denominator})
    auxiliary_three_adic = []
    # This branch is unnecessary for PPTs, but verifies the broader coprime-pair argument.
    for aa, bb in ((7, 1), (4, 7), (25, 1)):
        dd = aa + 2 * bb
        e = valuation(dd, 3)
        av = bp_evaluate(A, aa, bb)
        x6 = F(bp_evaluate(N, aa, bb), dd * dd)
        X6 = 9 * x6 - 3 * av
        check(e >= 2 and valuation(X6, 3) == 2 - 2 * e, "auxiliary_non_ppt_three_adic_branch")
        auxiliary_three_adic.append({"a": aa, "b": bb, "d": dd, "e": e, "v3_short_X6": valuation(X6, 3)})
    return {
        "status": "Exact identities in Z[a,b], finite modular implication, and finite PPT samples",
        "cubic_discriminant_factorization": "-65536*a^5*b^4*(a^3+6*a^2*b+14*a*b^2-2*b^3)",
        "elliptic_discriminant_note": "The Weierstrass curve discriminant is 16 times the cubic discriminant.",
        "rational_root_candidate_values": rational_root_values,
        "x6_numerator_terms": [{"a_power": i, "b_power": j, "coefficient": coefficient}
                               for (i, j), coefficient in sorted(N.items())],
        "x6_denominator": "(a+2*b)^2",
        "short_model": {"X": "9*x-3*A", "Z": "27*y", "a4": "81*U-27*A^2",
                        "a6": "729*V+243*A*U-54*A^3"},
        "exact_reduced_denominator_for_primitive_ppt": "((a+2*b)/gcd(a+2*b,2))^2",
        "ordered_ppt_samples": count, "sample_examples": examples,
        "auxiliary_non_ppt_three_adic_samples": auxiliary_three_adic
    }


def family_certificate():
    # Sparse polynomial ring Z[T], represented canonically by degree -> coefficient.
    m, n = {0: 5, 1: 6}, {0: 4, 1: 5}
    a, b, c = {0: 9, 1: 20, 2: 11}, {0: 40, 1: 98, 2: 60}, {0: 41, 1: 100, 2: 61}
    check(sp_sub(sp_square(m), sp_square(n)) == a, "family_ppt_first_leg_identity")
    check(sp_scale(sp_mul(m, n), 2) == b, "family_ppt_second_leg_identity")
    check(sp_add(sp_square(m), sp_square(n)) == c, "family_ppt_hypotenuse_identity")
    check(sp_add(sp_square(a), sp_square(b)) == sp_square(c), "family_pythagorean_identity")
    check(sp_sub(sp_scale(n, 6), sp_scale(m, 5)) == {0: -1}, "family_bezout_identity")
    # L^4 is computed by integer matrix multiplication, not just assumed.
    L = ((2, -1), (1, 0))
    L4 = ((1, 0), (0, 1))
    for _ in range(4):
        L4 = tuple(tuple(sum(L4[i][r] * L[r][j] for r in range(2))
                         for j in range(2)) for i in range(2))
    check(L4 == ((5, -4), (4, -3)), "family_fourth_parabolic_matrix")
    source_m, source_n = {0: 1, 1: 2}, {1: 1}
    check(sp_add(sp_scale(source_m, L4[0][0]), sp_scale(source_n, L4[0][1])) == m
          and sp_add(sp_scale(source_m, L4[1][0]), sp_scale(source_n, L4[1][1])) == n,
          "family_outer_parameter_bridge")
    A = sp_scale(sp_sub(sp_square(c), sp_scale(sp_mul(a, b), 4)), 2)
    B = sp_scale(sp_square(c), 2)
    C = sp_scale(sp_add(sp_square(c), sp_scale(sp_mul(a, b), 4)), 2)
    Y = sp_scale(sp_mul(sp_square(a), b), 16)
    s1, s2 = sp_scale(sp_sub(b, a), 2), sp_scale(sp_add(a, b), -2)
    U = sp_sub(sp_scale(sp_mul(Y, s1), 2), sp_square(A))
    V = sp_sub(sp_square(Y), sp_mul(A, U))
    check(sp_add(A, B) == sp_square(s1), "family_first_square_edge_identity")
    check(sp_add(B, C) == sp_square(s2), "family_second_square_edge_identity")
    check(sp_add(A, C) == sp_square(sp_scale(c, 2)), "family_triangle_closing_edge_identity")
    check(sp_add(A, C) == sp_scale(B, 2), "family_arithmetic_progression_identity")
    check(sp_add(sp_scale(Y, 2), sp_mul(sp_add(s1, s2), sp_sub(B, A))) == {},
          "family_interpolation_y_identity")
    certificates = {
        "A": A,
        "3U_minus_A_squared": sp_sub(sp_scale(U, 3), sp_square(A)),
        "minus_V": sp_scale(V, -1)
    }
    expected = {
        "A": [482, 2944, 6484, 6176, 2162],
        "3U_minus_A_squared": [18355184, 170125952, 686630080, 1575122048,
                              2244442016, 2032158592, 1140210880, 361836928, 49605104],
        "minus_V": [299007352, 8974023936, 88616121872, 465807661952,
                    1540102485064, 3461402023936, 5494932857312, 6256099417344,
                    5094172876424, 2902177296128, 1100255644432, 249542740864, 25624171192]
    }
    for name, poly in certificates.items():
        coefficients = sp_coefficients(poly)
        check(coefficients == expected[name], "family_exact_positive_coefficient_array")
        check(all(coefficient > 0 for coefficient in coefficients), "family_all_coefficients_strictly_positive")
    # Substitute T=154k symbolically. Every nonconstant coefficient vanishes mod 7 and 11.
    inherited_polys = {"a2": sp_scale(A, -1), "a4": U, "a6": V, "Px": A, "Py": Y}
    expected_constants = {"a2": -482, "a4": 6195836, "a6": -299007352, "Px": 482, "Py": 51840}
    residues = {}
    for prime in (7, 11):
        residue = {}
        for name, poly in inherited_polys.items():
            reduced = {degree: coefficient * 154 ** degree % prime
                       for degree, coefficient in poly.items()
                       if coefficient * 154 ** degree % prime}
            expected_constant = expected_constants[name] % prime
            check(reduced == ({0: expected_constant} if expected_constant else {}),
                  "family_symbolic_modular_inheritance")
            residue[name] = expected_constant
        residues[str(prime)] = residue
    samples = []
    for T_step, family_index in [(step, index) for step in (2, 154) for index in (0, 1, 2, 5)]:
        T = T_step * family_index
        mm, nn, aa, bb, cc = [sp_evaluate(poly, T) for poly in (m, n, a, b, c)]
        check(mm > nn > 0 and gcd(mm, nn) == 1 and (mm - nn) % 2 == 1,
              "family_sample_primitive_parameters")
        check(aa > 0 and bb > 0 and gcd(aa, bb) == 1 and aa * aa + bb * bb == cc * cc,
              "family_sample_primitive_triangle")
        labelA, labelB, labelC, yy = [sp_evaluate(poly, T) for poly in (A, B, C, Y)]
        model = tuple(map(F, (-labelA, sp_evaluate(U, T), sp_evaluate(V, T))))
        point = (F(labelA), F(yy))
        check(curve_discriminant(model) != 0 and on_curve(point, model), "family_sample_nonsingular_point")
        check(interpolate(labelA, labelB, labelC, sp_evaluate(s1, T), sp_evaluate(s2, T))
              == (model, point), "family_sample_interpolation")
        orbit, Q = [], O
        for j in range(1, 9):
            Q = add(Q, point, model)
            check(Q is not O and Q[0] > 0 and on_curve(Q, model), "family_sample_positive_orbit")
            check(Q == multiply(j, point, model), "family_sample_binary_agreement")
            orbit.append(Q)
        check([Q[0] for Q in orbit[:3]] == [labelA, labelB, labelC], "family_sample_first_triple")
        check(len({Q[0] for Q in orbit}) == len(orbit), "family_sample_distinct_orbit")
        denominators = []
        for Q in orbit:
            d = isqrt(Q[0].denominator)
            check(d * d == Q[0].denominator, "family_sample_square_denominator")
            denominators.append(d)
        H = lcm(*denominators)
        labels = [H * H * Q[0] for Q in orbit]
        check(all(label.denominator == 1 for label in labels), "family_sample_scaled_integer_labels")
        for j in range(7):
            edge_slope = slope(orbit[j], point, model)
            check(orbit[j][0] + orbit[j + 1][0] == edge_slope * edge_slope,
                  "family_sample_rational_edge")
            scaled_root = H * edge_slope
            check(scaled_root.denominator == 1
                  and labels[j] + labels[j + 1] == scaled_root * scaled_root,
                  "family_sample_integer_edge")
        samples.append({"T_step": T_step, "family_index": family_index, "T": T, "m": mm, "n": nn,
                        "ppt": [aa, bb, cc], "seed_labels": [labelA, labelB, labelC],
                        "tested_vertices": len(orbit), "scale_root_bits": H.bit_length(),
                        "maximum_scaled_label_bits": max(v.numerator.bit_length() for v in labels)})
    return {
        "coefficient_order": "ascending degree in T",
        "positive_coefficient_certificates": expected,
        "all_polynomials": {name: sp_coefficients(poly) for name, poly in
                            {"a": a, "b": b, "hypotenuse": c, "A": A, "B": B,
                             "C_label": C, "Y": Y, "s1": s1, "s2": s2, "U": U, "V": V}.items()},
        "L_fourth_power": [list(row) for row in L4],
        "symbolic_reductions_after_T_equals_154k": residues,
        "finite_sample_paths": samples
    }


def run(length):
    check(length >= 6, "length_at_least_six")
    a2, a4, a6 = MODEL
    check(a2 + P[0] == 0, "offset_cancellation")
    check(on_curve(P), "base_point_on_curve")
    D = curve_discriminant(MODEL)
    check(D == -14863179074355855360000, "exact_curve_discriminant")
    derivative_discriminant = (2 * a2) ** 2 - 12 * a4
    check(derivative_discriminant == -73420736, "exact_derivative_discriminant")
    check(derivative_discriminant < 0 and a6 < 0, "strict_monotonicity_positive_root")
    reductions = {}
    for prime in (7, 11):
        check(D.denominator == 1 and D.numerator % prime != 0, "good_reduction")
        points_mod_p = [(x, y) for x in range(prime) for y in range(prime)
                        if (y * y - cubic(F(x))).numerator % prime == 0]
        reductions[str(prime)] = {
            "order": 1 + len(points_mod_p), "discriminant_residue": int(D % prime),
            "affine_points": [list(Q) for Q in points_mod_p]}
    check(reductions["7"]["order"] == 10, "order_mod_seven")
    check(reductions["11"]["order"] == 12, "order_mod_eleven")
    check(gcd(reductions["7"]["order"], reductions["11"]["order"]) == 2,
          "torsion_order_divides_two")
    check(P[1] != 0 and add(P, P) is not O, "base_point_not_two_torsion")

    interpolated = interpolate(482, 3362, 6242, 62, -98)
    check(interpolated == (MODEL, P), "recovers_seeded_curve_and_point")
    for A, B, C, s1, s2 in ((1, 3, 6, 2, 3), (14, 11, 110, 5, 11),
                           (482, 3362, 6242, 62, 98)):
        for sign1, sign2 in ((1, 1), (1, -1), (-1, 1), (-1, -1)):
            interpolate(A, B, C, sign1 * s1, sign2 * s2)
    singular_model, singular_point = interpolate(14, 11, 110, 5, 11)
    check(singular_model == (F(-14), F(44), F(-40))
          and singular_point == (F(14), F(24)), "singular_control_coefficients")
    check(curve_discriminant(singular_model) == 0, "singularity_is_not_automatic")
    check(poly_mul(poly_mul([F(-2), F(1)], [F(-2), F(1)]), [F(-10), F(1)])
          == [F(-40), F(44), F(-14), F(1)], "singular_control_factorization")

    points = []
    Q = O
    for n in range(1, length + 1):
        Q = add(Q, P)
        check(Q is not O and on_curve(Q), "multiple_on_curve")
        check(Q == multiply(n, P), "sequential_vs_binary_multiplication")
        check(Q[0] > 0, "positive_abscissa")
        points.append(Q)
    values = [Q[0] for Q in points]
    check(len(set(values)) == length, "distinct_path_labels")
    check(values[:6] == list(map(F, (482, 3362, 6242, 158, 57442)))
          + [F(8198402, 7921)], "first_six_abscissas")
    check(values[0] + values[2] == 82 ** 2, "retained_first_third_chord")

    slopes = []
    for n in range(1, length):
        x, y = points[n - 1]
        xx, yy = points[n]
        m = slope(points[n - 1], P)
        b = P[1] - m * P[0]
        slopes.append(m)
        check(m * x + b == y and m * xx + b == -yy, "chord_line_incidence")
        polynomial = [a6 - b * b, a4 - 2 * m * b, a2 - m * m, F(1)]
        factored = poly_mul(poly_mul([-x, F(1)], [-P[0], F(1)]), [-xx, F(1)])
        check(polynomial == factored, "full_chord_polynomial_identity")
        check(x + xx == m * m, "adjacent_square_sum")
    check(slopes[:5] == [F(62), F(-98), F(80), F(240), F(-21522, 89)],
          "first_five_slopes")

    denominators = []
    observations = []
    for n, Q in enumerate(points, 1):
        d = isqrt(Q[0].denominator)
        check(d * d == Q[0].denominator, "square_reduced_x_denominator")
        denominators.append(d)
        observations.append({"n": n, "x_bits": rational_bits(Q[0]),
                             "y_bits": rational_bits(Q[1]), "denominator_root_bits": d.bit_length()})
    scalings = []
    for k in sorted(set((5, 6, min(8, length), min(16, length), length))):
        H = lcm(*denominators[:k])
        scaled = [H * H * x for x in values[:k]]
        check(all(x.denominator == 1 for x in scaled), "scaled_labels_integral")
        labels = [x.numerator for x in scaled]
        check(min(labels) > 0 and len(set(labels)) == k, "scaled_labels_positive_distinct")
        roots = []
        for i, m in enumerate(slopes[:k - 1]):
            root = abs(H * m)
            check(root.denominator == 1, "scaled_edge_root_integral")
            r = root.numerator
            check(labels[i] + labels[i + 1] == r * r, "scaled_edge_square_sum")
            roots.append(r)
        check(labels[0] + labels[2] == (82 * H) ** 2, "scaled_triangle_extra_edge")
        item = {"vertices": k, "H_bits": H.bit_length(),
                "maximum_label_bits": max(x.bit_length() for x in labels),
                "labels_sha256": digest_integers(labels)}
        if k <= 6:
            item.update({"H": H, "labels": labels, "edge_roots": roots})
        if k == 5:
            check(H == 1 and labels == [482, 3362, 6242, 158, 57442],
                  "unscaled_integer_five_path")
        if k == 6:
            check(H == 89, "six_prefix_requires_scaling")
            check(labels[:3] != [482, 3362, 6242], "scaling_changes_initial_labels")
        scalings.append(item)

    family = family_certificate()
    primitive_ppt = primitive_ppt_certificate()
    return {
        "status": "FINITE-EXACT; standard library; checks active under python -O",
        "length": length, "checks": dict(sorted(checks.items())),
        "check_count": sum(checks.values()),
        "curve": {"a2": int(a2), "a4": int(a4), "a6": int(a6),
                  "P": list(map(int, P)), "discriminant": int(D),
                  "derivative_discriminant": int(derivative_discriminant)},
        "good_reduction_controls": reductions,
        "first_six_points": [{"n": n, "x": str(Q[0]), "y": str(Q[1])}
                             for n, Q in enumerate(points[:6], 1)],
        "first_five_slopes": list(map(str, slopes[:5])),
        "common_square_scalings": scalings,
        "height_observations_not_asymptotic_proof": observations,
        "singular_interpolation_counterexample": {
            "A_B_C": [14, 11, 110], "slopes": [5, 11], "Y": 24,
            "a2_a4_a6": [-14, 44, -40], "factorization": "(x-2)^2(x-10)"},
        "infinite_ppt_family_certificate": family,
        "general_primitive_ppt_certificate": primitive_ppt,
        "scope": "The rational path preserves the exact triple. Integer scaling preserves a square-scaled copy; it changes the initial labels for prefixes of length at least six."
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--length", type=int, default=32)
    args = parser.parse_args()
    print(json.dumps(run(args.length), indent=2, sort_keys=True))
