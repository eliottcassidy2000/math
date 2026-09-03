#!/usr/bin/env python3
"""Independent exact audit of the proposed THM-4392 Poisson duality.

This program imports no producer or repository mathematics.  Exact finite
checks use integers/Fractions and explicit exceptions active under python -O.
No numerical Fejer truncation is used as theorem evidence.
"""

from fractions import Fraction
from hashlib import sha256
from itertools import permutations, product
from math import gcd, isqrt
import json


LAMBDA = Fraction(1, 14)
R = Fraction(3, 14)
CHECKS = 0


def check(condition, label):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(f"CHECK FAILED: {label}")


def ff(value):
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def gcd_many(values):
    g = 0
    for value in values:
        g = gcd(g, abs(value))
    return g


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def cross(a, b):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def extended_gcd(a, b):
    old_r, r = abs(a), abs(b)
    old_s, s = 1, 0
    old_t, t = 0, 1
    while r:
        q = old_r // r
        old_r, r = r, old_r - q * r
        old_s, s = s, old_s - q * s
        old_t, t = t, old_t - q * t
    if a < 0:
        old_s = -old_s
    if b < 0:
        old_t = -old_t
    return old_r, old_s, old_t


def bezout_vector(w):
    g12, x, y = extended_gcd(w[0], w[1])
    g, s, t = extended_gcd(g12, w[2])
    check(g == 1, f"primitive Bezout input {w}")
    u = (s * x, s * y, t)
    check(dot(w, u) == 1, f"Bezout identity {w}, {u}")
    return u


def kernel_basis(w):
    """A primitive oriented basis of Lambda_w with cross product -w."""
    g, x, y = extended_gcd(w[0], w[1])
    b1 = (w[1] // g, -w[0] // g, 0)
    b2 = (-w[2] * x, -w[2] * y, g)
    check(dot(w, b1) == dot(w, b2) == 0, f"kernel basis {w}")
    check(cross(b1, b2) == tuple(-x for x in w), f"primitive kernel area {w}")
    return b1, b2


def dual_basis(b1, b2):
    g11, g12, g22 = dot(b1, b1), dot(b1, b2), dot(b2, b2)
    determinant = g11 * g22 - g12 * g12
    d1 = tuple(Fraction(g22 * b1[i] - g12 * b2[i], determinant) for i in range(3))
    d2 = tuple(Fraction(-g12 * b1[i] + g11 * b2[i], determinant) for i in range(3))
    check(dot(d1, b1) == 1 and dot(d1, b2) == 0, "first dual vector")
    check(dot(d2, b1) == 0 and dot(d2, b2) == 1, "second dual vector")
    return d1, d2, determinant


def integer_tuple(v):
    check(all(isinstance(x, int) or x.denominator == 1 for x in v), f"integral tuple {v}")
    return tuple(int(x) for x in v)


def audit_lattice_normalizations(w):
    norm2 = dot(w, w)
    b1, b2 = kernel_basis(w)
    d1, d2, gram_det = dual_basis(b1, b2)
    check(gram_det == norm2, f"covolume squared {w}, {gram_det}")
    j1 = integer_tuple(cross(w, d1))
    j2 = integer_tuple(cross(w, d2))
    check(dot(w, j1) == dot(w, j2) == 0, f"J image in Lambda {w}")
    check(cross(j1, j2) in (w, tuple(-x for x in w)), f"J dual basis is lattice basis {w}")

    # T(e)=-w cross e has T T^t=||w||^2 I-w w^t.
    standard = ((1, 0, 0), (0, 1, 0), (0, 0, 1))
    columns = [tuple(-x for x in cross(w, e)) for e in standard]
    for i in range(3):
        for j in range(3):
            actual = sum(columns[k][i] * columns[k][j] for k in range(3))
            expected = norm2 * (1 if i == j else 0) - w[i] * w[j]
            check(actual == expected, f"cross-map singular normalization {w}, {i}, {j}")

    # Lift the two raw mod-three classes to honest elements of Lambda_w.
    residue_u = tuple(1 if x % 3 == 1 else -1 for x in w)
    check(all((residue_u[i] * w[i] - 1) % 3 == 0 for i in range(3)),
          f"inverse speed residue {w}")
    quotient = dot(w, residue_u) // 3
    check(dot(w, residue_u) == 3 * quotient, f"raw residue lies in kernel mod three {w}")
    bezout = bezout_vector(w)
    c0 = tuple(residue_u[i] - 3 * quotient * bezout[i] for i in range(3))
    check(dot(w, c0) == 0, f"raw coset lift {w}, {c0}")
    check(all((c0[i] - residue_u[i]) % 3 == 0 for i in range(3)), f"raw coset residue {w}")

    # Exact constants: rank-two dilation 3 gives 9, three transforms give 27.
    check(Fraction(1, 9) * 27 == 3, "Poisson residual 27/9 factor")
    return {
        "w": w,
        "norm2": norm2,
        "kernel_basis": (b1, b2),
        "dual_J_basis": (j1, j2),
        "c0": c0,
        "rank2_dilation": 9,
        "transform_product": 27,
        "residual": 3,
    }


def root_pair(exponent):
    """zeta^exponent in the integral basis (1,zeta), zeta^2=-1-zeta."""
    return ((1, 0), (0, 1), (-1, -1))[exponent % 3]


def sheet_character(w_residue, n_residue):
    a = b = 0
    for assignment in permutations((0, 1, 2)):
        exponent = sum(w_residue[i] * n_residue[i] * assignment[i] for i in range(3))
        x, y = root_pair(exponent)
        a += x
        b += y
    check(b == 0, f"real sheet character {w_residue}, {n_residue}, {(a,b)}")
    return a


def audit_mod_three_character():
    weight_counts = {6: 0, -3: 0}
    for w in product((1, 2), repeat=3):
        inverse = tuple(pow(x, -1, 3) for x in w)
        kernel = {
            C for C in product(range(3), repeat=3)
            if dot(w, C) % 3 == 0
        }
        raw_allowed = {C for C in kernel if all(x % 3 for x in C)}
        expected_raw = {inverse, tuple((-x) % 3 for x in inverse)}
        check(len(kernel) == 9, f"rank-two residue kernel {w}")
        check(raw_allowed == expected_raw, f"two raw residue cosets {w}, {raw_allowed}")
        for z in product(range(3), repeat=3):
            n = tuple(x % 3 for x in cross(w, z))
            weighted = tuple((w[i] * n[i]) % 3 for i in range(3))
            check(sum(weighted) % 3 == 0, f"weighted resonance residues {w}, {z}")
            equal = len(set(weighted)) == 1
            distinct = len(set(weighted)) == 3
            check(equal or distinct, f"zero-sum ternary dichotomy {w}, {z}, {weighted}")
            phase = dot(w, z) % 3
            check(equal == (phase == 0), f"phase/weight correspondence {w}, {z}")
            character = sheet_character(w, n)
            expected = 6 if equal else -3
            check(character == expected, f"sheet weight {w}, {z}, {character}")
            check(3 * (2 if phase == 0 else -1) == character,
                  f"two-coset transform times 27/9 {w}, {z}")
            weight_counts[character] += 1
    check(weight_counts == {6: 72, -3: 144}, f"finite character counts {weight_counts}")
    return weight_counts


def strict_bound(value):
    return (value.numerator - 1) // value.denominator


def component_length(w, C):
    w1, w2, w3 = w
    c1, c2, c3 = C
    return max(Fraction(0), min(
        2 * R / w1,
        2 * R / w2,
        2 * R / w3,
        R / w1 + R / w2 - Fraction(abs(c3), w1 * w2),
        R / w1 + R / w3 - Fraction(abs(c2), w1 * w3),
        R / w2 + R / w3 - Fraction(abs(c1), w2 * w3),
    ))


def raw_carriers(w):
    bounds = (
        strict_bound(R * (w[1] + w[2])),
        strict_bound(R * (w[0] + w[2])),
        strict_bound(R * (w[0] + w[1])),
    )
    result = {}
    for c1 in range(-bounds[0], bounds[0] + 1):
        for c2 in range(-bounds[1], bounds[1] + 1):
            numerator = -(w[0] * c1 + w[1] * c2)
            if numerator % w[2]:
                continue
            c3 = numerator // w[2]
            if abs(c3) > bounds[2]:
                continue
            C = (c1, c2, c3)
            if any(x % 3 == 0 for x in C):
                continue
            length = component_length(w, C)
            if length:
                result[C] = length
    return result


def eligible_intervals(speed, radius):
    inverse = pow(speed, -1, 3)
    intervals = []
    for n in range(speed + 1):
        left = max(Fraction(0), (Fraction(n) - radius) / speed)
        right = min(Fraction(1), (Fraction(n) + radius) / speed)
        if left < right:
            intervals.append((left, right, n, (-inverse * n) % 3))
    return intervals


def physical_y_components(w):
    lists = tuple(eligible_intervals(speed, R) for speed in w)
    indices = [0, 0, 0]
    result = {}
    while all(indices[i] < len(lists[i]) for i in range(3)):
        current = tuple(lists[i][indices[i]] for i in range(3))
        left = max(row[0] for row in current)
        right = min(row[1] for row in current)
        owners = tuple(row[3] for row in current)
        if left < right and len(set(owners)) == 3:
            n = tuple(row[2] for row in current)
            C = cross(w, n)
            result[C] = result.get(C, Fraction(0)) + right - left
        first_right = min(row[1] for row in current)
        for i in range(3):
            if current[i][1] == first_right:
                indices[i] += 1
    return result


def sheet_intervals(speed, sheet):
    radius = LAMBDA / speed
    intervals = []
    for k in range(speed):
        center = Fraction(k, speed) - Fraction(sheet, 3)
        center -= center.numerator // center.denominator
        left, right = center - radius, center + radius
        if left < 0:
            intervals.append((Fraction(0), right))
            intervals.append((left + 1, Fraction(1)))
        elif right > 1:
            intervals.append((left, Fraction(1)))
            intervals.append((Fraction(0), right - 1))
        else:
            intervals.append((left, right))
    intervals.sort()
    return intervals


def intersect_three_interval_lists(lists):
    indices = [0, 0, 0]
    total = Fraction(0)
    pieces = 0
    while all(indices[i] < len(lists[i]) for i in range(3)):
        current = tuple(lists[i][indices[i]] for i in range(3))
        left = max(row[0] for row in current)
        right = min(row[1] for row in current)
        if left < right:
            total += right - left
            pieces += 1
        first_right = min(row[1] for row in current)
        for i in range(3):
            if current[i][1] == first_right:
                indices[i] += 1
    return total, pieces


def shifted_six_sheet_measure(w):
    total = Fraction(0)
    pieces = 0
    sheet_cache = {(i, j): sheet_intervals(w[i], j) for i in range(3) for j in range(3)}
    for assignment in permutations((0, 1, 2)):
        subtotal, subpieces = intersect_three_interval_lists(
            tuple(sheet_cache[(i, assignment[i])] for i in range(3))
        )
        total += subtotal
        pieces += subpieces
    return total, pieces


CONTROL_CLAIMS = {
    (1, 5, 11): (Fraction(6, 77), 2, 6),
    (5, 19, 23): (Fraction(2, 665), 2, 6),
    (11, 13, 17): (Fraction(20, 1547), 2, 6),
    (101, 103, 107): (Fraction(36488, 7791847), 10, 30),
    (1, 11, 175): (Fraction(36, 1225), 12, 36),
}

EXPECTED_RAW = {
    (11, 13, 17): {
        (-5, -1, 4): Fraction(10, 1547),
        (5, 1, -4): Fraction(10, 1547),
    },
    (101, 103, 107): {
        (-41, 8, 31): Fraction(4, 11021),
        (-35, -1, 34): Fraction(10, 11021),
        (-29, -10, 37): Fraction(47, 72821),
        (-23, -19, 40): Fraction(26, 72821),
        (-17, -28, 43): Fraction(5, 72821),
        (17, 28, -43): Fraction(5, 72821),
        (23, 19, -40): Fraction(26, 72821),
        (29, 10, -37): Fraction(47, 72821),
        (35, 1, -34): Fraction(10, 11021),
        (41, -8, -31): Fraction(4, 11021),
    },
}


def audit_exact_physical_controls():
    summaries = []
    for w, (expected_measure, expected_raw_count, expected_sheet_count) in CONTROL_CLAIMS.items():
        raw = raw_carriers(w)
        physical = physical_y_components(w)
        shifted_measure, shifted_count = shifted_six_sheet_measure(w)
        measure = sum(raw.values(), Fraction(0))
        check(raw == physical, f"raw/literal-y carrier dictionary {w}")
        check(measure == expected_measure, f"exact raw measure {w}, {measure}")
        check(shifted_measure == expected_measure, f"six shifted-sheet measure {w}, {shifted_measure}")
        check(len(raw) == expected_raw_count, f"raw component count {w}")
        check(shifted_count == expected_sheet_count, f"shifted component count {w}, {shifted_count}")
        if w in EXPECTED_RAW:
            check(raw == EXPECTED_RAW[w], f"complete hostile carrier dictionary {w}")
        summaries.append((w, measure, len(raw), shifted_count))
    return tuple(summaries)


def resonance_vectors_euclidean(w, bound_squared):
    radius = isqrt(bound_squared)
    result = []
    for n1 in range(-radius, radius + 1):
        for n2 in range(-radius, radius + 1):
            for n3 in range(-radius, radius + 1):
                n = (n1, n2, n3)
                norm2 = dot(n, n)
                if n != (0, 0, 0) and norm2 <= bound_squared and dot(w, n) == 0:
                    result.append(n)
    return tuple(result)


def resonance_vectors_product(w, product_bound):
    result = set()
    for n1 in range(-product_bound, product_bound + 1):
        if n1 == 0 or n1 % 7 == 0:
            continue
        for n2 in range(-product_bound, product_bound + 1):
            if n2 == 0 or n2 % 7 == 0:
                continue
            numerator = -(w[0] * n1 + w[1] * n2)
            if numerator % w[2]:
                continue
            n3 = numerator // w[2]
            if n3 == 0 or n3 % 7 == 0:
                continue
            value = abs(n1 * n2 * n3)
            if value <= product_bound:
                result.add((n1, n2, n3))
    return tuple(sorted(result))


def concrete_sheet_weight(w, n):
    residues_w = tuple(x % 3 for x in w)
    residues_n = tuple(x % 3 for x in n)
    return sheet_character(residues_w, residues_n)


def fourier_factor_nonzero(n):
    """At lambda=1/14: zero only at a nonzero coordinate divisible by seven."""
    return all(coordinate == 0 or coordinate % 7 != 0 for coordinate in n)


def audit_shortest_orbit_and_termwise_scope():
    w1 = (11, 13, 17)
    w2 = (101, 103, 107)
    minimizer = (-2, 3, -1)
    expected_orbit = (minimizer, tuple(-x for x in minimizer))
    summaries = []
    populations = {}
    for w in (w1, w2):
        vectors = resonance_vectors_euclidean(w, 14)
        minimum_norm = min(dot(n, n) for n in vectors)
        minimum_orbit = tuple(sorted(n for n in vectors if dot(n, n) == minimum_norm))
        check(minimum_norm == 14, f"Euclidean minimum {w}")
        check(minimum_orbit == tuple(sorted(expected_orbit)), f"complete Euclidean orbit {w}")

        product_vectors = resonance_vectors_product(w, 6)
        minimum_product = min(abs(n[0] * n[1] * n[2]) for n in product_vectors)
        product_orbit = tuple(n for n in product_vectors
                              if abs(n[0] * n[1] * n[2]) == minimum_product)
        check(minimum_product == 6, f"m7 minimum {w}")
        check(product_orbit == tuple(sorted(expected_orbit)), f"complete m7 orbit {w}")
        check(concrete_sheet_weight(w, minimizer) == -3, f"minimizer sheet weight {w}")
        counts = tuple((bound, len(resonance_vectors_product(w, bound)))
                       for bound in (6, 12, 24, 48, 96))
        populations[w] = counts
        summaries.append((w, minimum_norm, minimum_orbit, minimum_product, product_orbit))

    check(populations[w1] == ((6, 2), (12, 2), (24, 4), (48, 6), (96, 12)),
          f"near-minimal population first {populations[w1]}")
    check(populations[w2] == ((6, 2), (12, 2), (24, 2), (48, 4), (96, 4)),
          f"near-minimal population second {populations[w2]}")

    mu1 = CONTROL_CLAIMS[w1][0]
    mu2 = CONTROL_CLAIMS[w2][0]
    check(mu1 - mu2 == Fraction(14198572, 1721998187), "hostile measure difference")
    check(mu1 / mu2 == Fraction(5565605, 2015962), "hostile measure ratio")

    # These witnesses establish incomparability only under the natural
    # identity identification of the same embedded lattice vector.
    raw_deleted_fourier_live = minimizer
    check(dot(w2, raw_deleted_fourier_live) == 0, "first termwise witness in Lambda")
    check(any(x % 3 == 0 for x in raw_deleted_fourier_live), "first witness raw-deleted")
    check(fourier_factor_nonzero(raw_deleted_fourier_live), "first witness Fourier-live")
    check(concrete_sheet_weight(w2, raw_deleted_fourier_live) == -3,
          "first witness nonzero sheet character")

    raw_live_fourier_killed = (-35, -1, 34)
    raw2 = raw_carriers(w2)
    check(raw_live_fourier_killed in raw2, "second termwise witness raw-live")
    check(raw2[raw_live_fourier_killed] == Fraction(10, 11021), "second witness exact length")
    check(all(x % 3 for x in raw_live_fourier_killed), "second witness raw mod-three gate")
    check(not fourier_factor_nonzero(raw_live_fourier_killed), "second witness Fourier-killed")
    check(raw_live_fourier_killed[0] != 0 and raw_live_fourier_killed[0] % 7 == 0,
          "second witness nonzero mod-seven sinc zero")
    check(fourier_factor_nonzero((0, 1, -1)), "zero coordinate is not a sinc zero")

    return {
        "summaries": tuple(summaries),
        "populations": tuple((w, populations[w]) for w in (w1, w2)),
        "difference": mu1 - mu2,
        "ratio": mu1 / mu2,
        "identity_witnesses": (
            ("raw_deleted_fourier_live", raw_deleted_fourier_live),
            ("raw_live_fourier_killed", raw_live_fourier_killed, raw2[raw_live_fourier_killed]),
        ),
    }


def json_ready(value):
    if isinstance(value, Fraction):
        return ff(value)
    if isinstance(value, dict):
        return {str(key): json_ready(item) for key, item in value.items()}
    if isinstance(value, (tuple, list, set)):
        return [json_ready(item) for item in value]
    return value


def main():
    normalizations = tuple(audit_lattice_normalizations(w) for w in CONTROL_CLAIMS)
    weight_counts = audit_mod_three_character()
    physical = audit_exact_physical_controls()
    hostile = audit_shortest_orbit_and_termwise_scope()
    semantic = {
        "normalizations": normalizations,
        "weights": weight_counts,
        "physical": physical,
        "hostile": hostile,
        "convergence_dependency": "THM-1092-fejer-regularized-kfold-resonance-identity",
        "fejer_numeric_status": "diagnostic_not_used",
        "hybrid_status": "conjectural_not_used",
        "termwise_scope": "natural_identity_identification_only",
    }
    digest = sha256(json.dumps(json_ready(semantic), sort_keys=True).encode()).hexdigest()

    print("THM4392 BOXSPLINE-FOURIER CLEAN-ROOM AUDIT")
    print("status=PASS exact_local_identity; LRC14=OPEN")
    print("normalization=B=L/||w||; covol(Lambda)=||w||; covol(3Lambda)=9||w||")
    print("poisson_constants=one_ninth_times_twenty_seven_equals_three")
    print("finite_character_counts=" + str(tuple(sorted(weight_counts.items()))))
    print("mod3_role=two_primal_cosets_to_character_weights_6_or_minus3")
    print("mod7_role=interval_sinc_zero_only_for_nonzero_multiples; hhat(0)=1/7")
    print("physical_controls=" + str(physical))
    print("shortest_orbit=norm2_14;m7_6;orbit=+/-(-2,3,-1);same_residue_profile=(2,1,2)")
    print("hostile_measures=20/1547,36488/7791847 difference=" + ff(hostile["difference"])
          + " ratio=" + ff(hostile["ratio"]))
    print("near_minimal_counts=" + str(hostile["populations"]))
    print("termwise_verdict=identity_on_shared_embedded_lattice_does_not_intertwine_gates")
    print("termwise_nonclaim=arbitrary_bijections_or_blockwise_correspondences_not_excluded")
    print("convergence=THM1092_absolute_full_support_plus_direct_rank1_zero_strata")
    print("fejer_rows=NUMERIC_DIAGNOSTIC_NOT_USED")
    print("low_dual_high_primal_hybrid=CONJECTURAL_NOT_USED")
    print(f"checks={CHECKS}")
    print(f"semantic_sha256={digest}")
    print("verdict=PASS")


if __name__ == "__main__":
    main()
