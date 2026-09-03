#!/usr/bin/env python3
"""Exact/controlled-numeric verifier for the THM-4392 LRC Fourier bridge.

No repository implementation is imported.  Exact rational arithmetic checks
the raw carrier sum and the original six shifted combs.  Floating point is
used only for explicitly labelled Fejer convergence diagnostics.
"""

from fractions import Fraction
from itertools import permutations, product
from math import fsum, gcd, isqrt, pi, sin
import sys


Q = Fraction
LAMBDA = Q(1, 14)
RAW_RADIUS = Q(3, 14)
CHECKS = 0


class VerificationError(RuntimeError):
    pass


def require(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise VerificationError(message)


def gcd3(values):
    return gcd(gcd(abs(values[0]), abs(values[1])), abs(values[2]))


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def cross(a, b):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def ceil_fraction(value):
    return -((-value.numerator) // value.denominator)


def carrier_length(w, carrier):
    w1, w2, w3 = w
    c1, c2, c3 = carrier
    six = (
        2 * RAW_RADIUS / w1,
        2 * RAW_RADIUS / w2,
        2 * RAW_RADIUS / w3,
        RAW_RADIUS / w1 + RAW_RADIUS / w2 - Q(abs(c3), w1 * w2),
        RAW_RADIUS / w1 + RAW_RADIUS / w3 - Q(abs(c2), w1 * w3),
        RAW_RADIUS / w2 + RAW_RADIUS / w3 - Q(abs(c1), w2 * w3),
    )
    return max(Q(0), min(six))


def raw_carrier_measure(w):
    """Finite exact sum over ker_Z(w), retaining all-nonzero mod-three C."""
    w1, w2, w3 = w
    bounds = (
        ceil_fraction(RAW_RADIUS * (w2 + w3)),
        ceil_fraction(RAW_RADIUS * (w1 + w3)),
        ceil_fraction(RAW_RADIUS * (w1 + w2)),
    )
    rows = []
    for c1 in range(-bounds[0], bounds[0] + 1):
        for c2 in range(-bounds[1], bounds[1] + 1):
            numerator = -(w1 * c1 + w2 * c2)
            if numerator % w3:
                continue
            c3 = numerator // w3
            if abs(c3) > bounds[2]:
                continue
            carrier = (c1, c2, c3)
            require(dot(w, carrier) == 0, ("carrier-relation", w, carrier))
            if any(entry % 3 == 0 for entry in carrier):
                continue
            length = carrier_length(w, carrier)
            if length:
                rows.append((carrier, length))

    # The displayed strict pair bounds make every omitted carrier have length 0.
    require(all(abs(c[i]) < bounds[i] for c, _ in rows for i in range(3)),
            ("strict-carrier-support", w, bounds, rows))
    return sum((length for _, length in rows), Q(0)), tuple(rows)


def merge_intervals(intervals):
    merged = []
    for left, right in sorted(intervals):
        if not merged or left > merged[-1][1]:
            merged.append([left, right])
        elif right > merged[-1][1]:
            merged[-1][1] = right
    return tuple((left, right) for left, right in merged if left < right)


SHIFT_CACHE = {}


def shifted_danger_intervals(speed, sheet):
    """D_(speed,sheet) on x in [0,1], directly from its definition."""
    key = speed, sheet
    if key in SHIFT_CACHE:
        return SHIFT_CACHE[key]
    intervals = []
    radius = LAMBDA / speed
    for nearest in range(-speed, 2 * speed + 1):
        center = Q(nearest, speed) - Q(sheet, 3)
        left = max(Q(0), center - radius)
        right = min(Q(1), center + radius)
        if left < right:
            intervals.append((left, right))
    result = merge_intervals(intervals)
    SHIFT_CACHE[key] = result
    return result


def intersect_lists(left_intervals, right_intervals):
    result = []
    i = j = 0
    while i < len(left_intervals) and j < len(right_intervals):
        left = max(left_intervals[i][0], right_intervals[j][0])
        right = min(left_intervals[i][1], right_intervals[j][1])
        if left < right:
            result.append((left, right))
        left_end = left_intervals[i][1]
        right_end = right_intervals[j][1]
        if left_end <= right_end:
            i += 1
        if right_end <= left_end:
            j += 1
    return tuple(result)


def original_six_sheet_measure(w):
    total = Q(0)
    component_count = 0
    for sheets in permutations((0, 1, 2)):
        common = shifted_danger_intervals(w[0], sheets[0])
        common = intersect_lists(common, shifted_danger_intervals(w[1], sheets[1]))
        common = intersect_lists(common, shifted_danger_intervals(w[2], sheets[2]))
        total += sum((right - left for left, right in common), Q(0))
        component_count += len(common)
    return total, component_count


def resonance_vectors_in_cube(w, bound):
    w1, w2, w3 = w
    for n1 in range(-bound, bound + 1):
        for n2 in range(-bound, bound + 1):
            numerator = -(w1 * n1 + w2 * n2)
            if numerator % w3:
                continue
            n3 = numerator // w3
            if abs(n3) <= bound:
                yield n1, n2, n3


def sheet_weight(w, frequency):
    residues = tuple((wi * ni) % 3 for wi, ni in zip(w, frequency))
    if residues[0] == residues[1] == residues[2]:
        return 6
    if set(residues) == {0, 1, 2}:
        return -3
    raise VerificationError(("bad-resonance-residues", w, frequency, residues))


def sheet_weight_by_permutations(w, frequency):
    counts = [0, 0, 0]
    residues = tuple((wi * ni) % 3 for wi, ni in zip(w, frequency))
    for sheets in permutations((0, 1, 2)):
        exponent = sum(residue * sheet for residue, sheet in zip(residues, sheets)) % 3
        counts[exponent] += 1
    require(counts[1] == counts[2], ("nonreal-sheet-sum", w, frequency, counts))
    # 1 + omega + omega^2 = 0, so c0+c1*omega+c1*omega^2=c0-c1.
    return counts[0] - counts[1], tuple(counts)


def finite_character_audit():
    """Check that the primal two-coset mask transforms to the sheet weight."""
    counts = {-3: 0, 6: 0}
    for w in product((1, 2), repeat=3):
        for z in product(range(3), repeat=3):
            frequency = cross(w, z)
            require(dot(w, frequency) == 0,
                    ("finite-cross-relation", w, z, frequency))
            direct, root_counts = sheet_weight_by_permutations(w, frequency)

            # In F_3 every unit is its own inverse, so the retained raw-carrier
            # cosets are +/-u with u_i=w_i.  If eta=P_w z, its character on u
            # has phase z.u.  The Fourier transform of {+u,-u} is 2 or -1.
            phase = dot(w, z) % 3
            mask_transform = 2 if phase == 0 else -1
            require(direct == 3 * mask_transform,
                    ("finite-mask-transform", w, z, frequency, phase,
                     root_counts, direct, mask_transform))
            require(direct == sheet_weight(w, frequency),
                    ("sheet-classification", w, z, frequency))
            counts[direct] += 1
    require(counts == {-3: 144, 6: 72}, ("finite-character-counts", counts))
    return tuple(sorted(counts.items()))


def fourier_coefficient(n):
    if n == 0:
        return 1.0 / 7.0
    if n % 7 == 0:
        return 0.0
    return sin(pi * n / 7.0) / (pi * n)


def fejer_weight(n, bound):
    return max(1.0 - abs(n) / (bound + 1.0), 0.0)


def fejer_resonance_sum(w, bound):
    terms = []
    killed = 0
    for frequency in resonance_vectors_in_cube(w, bound):
        factors = []
        zero = False
        for n in frequency:
            coefficient = fourier_coefficient(n)
            if coefficient == 0.0:
                zero = True
                break
            factors.append(fejer_weight(n, bound) * coefficient)
        if zero:
            killed += 1
            continue
        terms.append(sheet_weight(w, frequency) * factors[0] * factors[1] * factors[2])
    return fsum(terms), killed, len(terms)


def full_support_m7(w):
    discovery_bound = 2
    candidate = None
    while candidate is None:
        for vector in resonance_vectors_in_cube(w, discovery_bound):
            if 0 in vector or any(entry % 7 == 0 for entry in vector):
                continue
            value = abs(vector[0] * vector[1] * vector[2])
            candidate = value if candidate is None else min(candidate, value)
        discovery_bound *= 2
        require(discovery_bound <= 4096, ("m7-discovery", w))

    # Any product at most candidate has every absolute coordinate at most it.
    minimizers = []
    best = candidate
    w1, w2, w3 = w
    for n1 in range(-candidate, candidate + 1):
        for n2 in range(-candidate, candidate + 1):
            if n1 == 0 or n2 == 0:
                continue
            numerator = -(w1 * n1 + w2 * n2)
            if numerator % w3:
                continue
            n3 = numerator // w3
            if n3 == 0 or any(entry % 7 == 0 for entry in (n1, n2, n3)):
                continue
            value = abs(n1 * n2 * n3)
            if value < best:
                best = value
                minimizers = [(n1, n2, n3)]
            elif value == best:
                minimizers.append((n1, n2, n3))
    require(best == candidate, ("m7-completeness", w, candidate, best))
    return best, tuple(sorted(set(minimizers)))


def euclidean_shortest(w):
    discovery_bound = 1
    candidate = None
    while candidate is None:
        for vector in resonance_vectors_in_cube(w, discovery_bound):
            if vector == (0, 0, 0):
                continue
            norm2 = dot(vector, vector)
            candidate = norm2 if candidate is None else min(candidate, norm2)
        discovery_bound *= 2
        require(discovery_bound <= 4096, ("euclidean-discovery", w))

    coordinate_bound = isqrt(candidate)
    vectors = []
    best = candidate
    w1, w2, w3 = w
    for n1 in range(-coordinate_bound, coordinate_bound + 1):
        for n2 in range(-coordinate_bound, coordinate_bound + 1):
            numerator = -(w1 * n1 + w2 * n2)
            if numerator % w3:
                continue
            n3 = numerator // w3
            vector = n1, n2, n3
            if vector == (0, 0, 0):
                continue
            norm2 = dot(vector, vector)
            if norm2 < best:
                best = norm2
                vectors = [vector]
            elif norm2 == best:
                vectors.append(vector)
    require(best == candidate, ("euclidean-completeness", w, candidate, best))
    return best, tuple(sorted(set(vectors)))


def near_multiplicative_count(w, ceiling):
    vectors = []
    w1, w2, w3 = w
    for n1 in range(-ceiling, ceiling + 1):
        for n2 in range(-ceiling, ceiling + 1):
            if n1 == 0 or n2 == 0:
                continue
            numerator = -(w1 * n1 + w2 * n2)
            if numerator % w3:
                continue
            n3 = numerator // w3
            vector = n1, n2, n3
            if n3 == 0 or any(entry % 7 == 0 for entry in vector):
                continue
            if abs(n1 * n2 * n3) <= ceiling:
                vectors.append(vector)
    return tuple(sorted(vectors))


def shortest_vector_obstruction():
    first = (11, 13, 17)
    second = (101, 103, 107)
    require(tuple(entry % 3 for entry in first) == tuple(entry % 3 for entry in second),
            "obstruction residue profiles differ")
    euclidean_first = euclidean_shortest(first)
    euclidean_second = euclidean_shortest(second)
    m7_first = full_support_m7(first)
    m7_second = full_support_m7(second)
    require(euclidean_first == euclidean_second,
            ("euclidean-shortest-not-equal", euclidean_first, euclidean_second))
    require(m7_first == m7_second,
            ("m7-not-equal", m7_first, m7_second))

    expected_orbit = ((-2, 3, -1), (2, -3, 1))
    require(euclidean_first == (14, expected_orbit),
            ("euclidean-shortest-orbit", euclidean_first))
    require(m7_first == (6, expected_orbit), ("m7-shortest-orbit", m7_first))
    require(sheet_weight(first, expected_orbit[0]) == -3,
            "first shortest sheet weight")
    require(sheet_weight(second, expected_orbit[0]) == -3,
            "second shortest sheet weight")

    first_measure, first_carriers = raw_carrier_measure(first)
    second_measure, second_carriers = raw_carrier_measure(second)
    require(first_measure != second_measure,
            ("shortest-data-falsely-determined-measure", first_measure, second_measure))
    populations = tuple(
        (ceiling,
         len(near_multiplicative_count(first, ceiling)),
         len(near_multiplicative_count(second, ceiling)))
        for ceiling in (6, 12, 24, 48, 96)
    )
    # The two coordinate deletions are dual, not termwise-identical.
    fourier_only = expected_orbit[0]
    raw_only = (-35, -1, 34)
    require(dot(second, raw_only) == 0, ("raw-only-relation", second, raw_only))
    require(all(entry % 3 for entry in raw_only), ("raw-only-mod3", raw_only))
    require(carrier_length(second, raw_only) == Q(10, 11021),
            ("raw-only-positive-length", raw_only))
    require(any(entry % 7 == 0 for entry in raw_only), ("raw-only-mod7", raw_only))
    require(any(entry % 3 == 0 for entry in fourier_only),
            ("fourier-only-is-raw-deleted", fourier_only))
    require(all(entry and entry % 7 for entry in fourier_only),
            ("fourier-only-is-mod7-live", fourier_only))
    termwise_hostiles = (
        ("raw_deleted_but_fourier_live", fourier_only),
        ("raw_live_but_fourier_killed", raw_only, Q(10, 11021)),
    )
    return ((first, first_measure, first_carriers),
            (second, second_measure, second_carriers),
            euclidean_first, m7_first, populations, termwise_hostiles)


def main():
    global CHECKS
    if "--tripwire" in sys.argv:
        try:
            require(False, "deliberate optimized-mode tripwire")
        except VerificationError:
            print("optimization_tripwire=LIVE")
            print(f"explicit_checks={CHECKS}")
            return
        raise VerificationError("tripwire vanished")

    require(RAW_RADIUS == 3 * LAMBDA, "radius scaling")
    require(Q(27, 9) == 3, "Poisson normalization factor")
    character_counts = finite_character_audit()

    samples = (
        (1, 5, 11),
        (5, 19, 23),
        (11, 13, 17),
        (101, 103, 107),
        (1, 11, 175),
    )
    exact_rows = []
    for w in samples:
        require(gcd3(w) == 1, ("sample-not-primitive", w))
        require(all(entry % 2 and entry % 3 for entry in w),
                ("sample-not-odd-three-units", w))
        raw_measure, carriers = raw_carrier_measure(w)
        direct_measure, components = original_six_sheet_measure(w)
        require(raw_measure == direct_measure,
                ("raw-vs-six-sheet", w, raw_measure, direct_measure))
        exact_rows.append((w, raw_measure, len(carriers), components))

    obstruction = shortest_vector_obstruction()

    # Product-Fejer values are diagnostics only; exact proof is analytic.
    fejer_rows = []
    for w in ((1, 5, 11), (11, 13, 17), (101, 103, 107)):
        exact = dict((row[0], row[1]) for row in exact_rows)[w]
        approximations = []
        for bound in (24, 48, 96, 192):
            value, killed, retained = fejer_resonance_sum(w, bound)
            approximations.append((bound, value, abs(value - float(exact)), killed, retained))
        require(approximations[-1][2] < 0.01,
                ("Fejer-diagnostic-not-close", w, exact, approximations[-1]))
        fejer_rows.append((w, exact, tuple(approximations)))

    print("status=PASS")
    print("scope=local_three_tail_scale_three_comb_only; LRC(14)=OPEN")
    print("identity=raw_mod3_coset_boxspline_sum_equals_six_sheet_mod7_sinc_resonance_sum")
    print(f"finite_character_weight_counts={character_counts}")
    print(f"exact_raw_vs_original_six_sheet={tuple(exact_rows)}")
    print(f"shortest_vector_obstruction={obstruction}")
    for w, exact, approximations in fejer_rows:
        formatted = tuple((bound, f"{value:.12f}", f"{error:.12f}", killed, retained)
                          for bound, value, error, killed, retained in approximations)
        print(f"fejer_diagnostic w={w} exact={exact} rows={formatted}")
    print("mod7_note=only_nonzero_multiples_of_7_are_killed; hhat(0)=1/7_survives")
    print("shortest_vector_verdict=bare_lambda1_and_m7_and_complete_minimizer_orbit_do_not_determine_measure")
    print("proof_status=boxspline_poisson_identity_PROVED_IN_REPORT; finite_checks_EXACT; Fejer_rows_NUMERIC_DIAGNOSTIC")
    print("optimization_safe_checks=yes")
    print(f"explicit_checks={CHECKS}")


if __name__ == "__main__":
    main()
