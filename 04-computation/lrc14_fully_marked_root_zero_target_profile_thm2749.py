#!/usr/bin/env python3
"""Exact two-sided fully marked root-zero clutch and target profile (THM-2749).

The source E3/lawful-target mask and its translated-target pullback, the
actual D^6 terminal fork, relative-present complements, rail, and half-tooth
seam are all inserted before the delayed-carry functional is evaluated.  It
is a two-sided common-section fibre product, not the natural one-sheet
coefficient.  The target is recomputed both in the pulled-back source
coordinate and directly after translation.

All arithmetic is integral/rational.  Dependency pins use LF-normalized
bytes so normal and optimized replay is portable across CRLF checkouts.
There are no truth-bearing Python assertions.
"""

from fractions import Fraction
from hashlib import sha256
from math import gcd
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
COMP = ROOT / "04-computation"
sys.path.insert(0, str(COMP))

DEPENDENCIES = {
    "lrc14_relative_present_semantic_lift_probe_20260728.py":
        "f16754bd38ae0dfa0d7d91cc404b4447dbf359635101aa7b4223363f8064352f",
    "lrc14_central_half_odometer_full_local_cycle_thm2698.py":
        "45cc393a856c00342fdf84875a0bc5a6d4c3df196ab35bb9ac2aad3cfc966c25",
    "lrc14_full_physical_lift_fibre_thm2707.py":
        "f05a07b2fb22cb5b39ed7d14e66d26154ecc50fc214861dc6576c3bcfaed2412",
    "lrc14_half_c221_semantic_source_fibre_census_20260728.py":
        "0fbeb041007fea1b9e14f0ff6e82fc97ebf724b26c2c10ef85732b4c994b94cd",
    "lrc14_predecessor_carry_private_root_atlas_thm2640.py":
        "a28b03a5903256c1c1c294ea5af389c7991fc0a5ad6908f0f25a5b0cc6e71abf",
    "lrc14_replica_dichotomy_typed_row_opus_20260727.py":
        "6ba64a68a9fd008d2e06949b1f1cf75012f1f4e734f75f55ce0af58ae20ad7b9",
    "lrc14_two_target_present_semantic_attachment_probe_20260728.py":
        "062b352f4db12a5f01822b293cdbb10629632dacc5fa27b406d8dd321e550709",
}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def lf_bytes(path):
    return path.read_bytes().replace(bytes((13, 10)), bytes((10,)))


for dependency, expected in DEPENDENCIES.items():
    actual = sha256(lf_bytes(COMP / dependency)).hexdigest()
    require(actual == expected,
            f"audited dependency changed: {dependency}: {actual} != {expected}")


import lrc14_full_physical_lift_fibre_thm2707 as lift
import lrc14_predecessor_carry_private_root_atlas_thm2640 as private
import lrc14_two_target_present_semantic_attachment_probe_20260728 as two


P = 13
R = P**6
T = private.T
SHIFT = 7 * T // R
CONTENT = 26
RAIL = 8
AMPLITUDE = 339633525654239542165440


def merge(intervals):
    result = []
    for left, right in sorted(intervals):
        if left >= right:
            continue
        if result and left <= result[-1][1]:
            result[-1] = (result[-1][0], max(result[-1][1], right))
        else:
            result.append((left, right))
    return tuple(result)


def intersect(first, second):
    i = j = 0
    result = []
    while i < len(first) and j < len(second):
        left = max(first[i][0], second[j][0])
        right = min(first[i][1], second[j][1])
        if left < right:
            result.append((left, right))
        if first[i][1] <= second[j][1]:
            i += 1
        else:
            j += 1
    return tuple(result)


def shift_union(intervals, shift, modulus=T):
    """Pull a circular union back under x -> x+shift."""
    result = []
    for left, right in intervals:
        start = (left - shift) % modulus
        stop = start + right - left
        if stop <= modulus:
            result.append((start, stop))
        else:
            result.extend(((start, modulus), (0, stop - modulus)))
    return merge(result)


def shift_weighted(pieces, shift, modulus=T):
    """Pull a weighted circular profile back under x -> x+shift."""
    result = []
    for left, right, weight in pieces:
        start = (left - shift) % modulus
        stop = start + right - left
        if stop <= modulus:
            result.append((start, stop, weight))
        else:
            result.extend(((start, modulus, weight),
                           (0, stop - modulus, weight)))
    return tuple(sorted(result))


def complement(intervals, modulus=T):
    result = []
    cursor = 0
    for left, right in merge(intervals):
        if cursor < left:
            result.append((cursor, left))
        cursor = max(cursor, right)
    if cursor < modulus:
        result.append((cursor, modulus))
    return tuple(result)


def restrict_weighted(pieces, intervals):
    if not intervals:
        return ()
    return tuple(private.old.intersect_weighted_union(
        pieces, intervals, tuple(left for left, _right in intervals)
    ))


def prefix_intervals(prefix):
    starts, lengths, _cumulative = prefix
    return tuple((start, start + length)
                 for start, length in zip(starts, lengths))


def marked_prefixes(module, delayed, fork):
    """Insert the actual terminal D^6 fork before delayed integration."""
    result = []
    for ell5 in range(7):
        by_kappa = []
        for kappa in range(2):
            intervals = intersect(
                prefix_intervals(delayed[0][ell5][6][kappa]), fork
            )
            prefix = module.make_prefix(intervals)
            require(prefix[2][-1] % P == 0,
                    "marked terminal prefix lost one-cell descent")
            by_kappa.append(prefix)
        result.append(tuple(by_kappa))
    return tuple(result)


def unit_reduction(values, root):
    require(root % P, "root normalization requested at zero")
    if not all(value % CONTENT == 0 for value in values):
        return None
    inverse = pow(root, -1, P)
    normalized = tuple(value // CONTENT * inverse % P for value in values)
    reduced = tuple((normalized[index] - normalized[-1]) % P
                    for index in range(6))
    determinant = private.old.sat.multiplication_determinant_7(reduced)
    return reduced, determinant


def static_semantic_common(module, source, clock_comb, rail_common):
    """Fixed (ell,s)=(1,0) semantic factors on both seam endpoints."""
    static_source = two.intersect_sorted(source, clock_comb[1])
    static_source = module.subtract_comb(
        static_source, module.W[1], 182, -13, 13
    )
    static_source = module.subtract_comb(
        static_source, module.C2, 182, -13, 13
    )
    static_target = shift_union(tuple(static_source), SHIFT)
    return intersect(
        rail_common,
        intersect(tuple(static_source), static_target),
    )


def target_label_common(module, static_common, t):
    """Insert both moving q2/c3 factors at source and target."""
    common = tuple(module.subtract_comb(
        static_common, module.W[2], 182,
        -14 * t - 13, -14 * t + 13,
    ))
    common = tuple(module.subtract_comb(
        common, module.C3, 182,
        14 * t - 13, 14 * t + 13,
    ))
    target_w2_danger = shift_union(
        tuple(module.make_comb(
            module.W[2], 182, -14 * t - 13, -14 * t + 13
        )),
        SHIFT,
    )
    common = intersect(common, complement(target_w2_danger))
    # c3*(7/R)=14/13, hence a 14-unit pullback on the 182-grid.
    return tuple(module.subtract_comb(
        common, module.C3, 182,
        14 * t - 13 - 14, 14 * t + 13 - 14,
    ))


def rail_data(rails, rail_index):
    source_weight = tuple(rails[rail_index][3])
    target_weight = shift_weighted(source_weight, SHIFT)
    source_support = merge(
        (left, right) for left, right, weight in source_weight if weight
    )
    target_support = merge(
        (left, right) for left, right, weight in target_weight if weight
    )
    return (source_weight, target_weight,
            intersect(source_support, target_support))


def clutch_vector(module, delayed, present, source_weight, target_weight,
                  common):
    source_vector = []
    target_vector = []
    for ell5 in range(7):
        source_safe = complement(present[ell5, 7])
        target_safe = complement(
            shift_union(present[ell5, 7], SHIFT)
        )
        row_common = intersect(common, source_safe)
        row_common = intersect(row_common, target_safe)
        source_row = restrict_weighted(source_weight, row_common)
        target_row = restrict_weighted(target_weight, row_common)
        source_row = private.old.intersect_weighted_comb(
            source_row, module.C3, 182, 169, 181
        )
        target_row = private.old.intersect_weighted_comb(
            target_row, module.C3, 182, 169, 181
        )
        source_value = private.delayed_carry_pair(
            source_row, delayed[ell5], {}
        )[12][1]
        target_value = private.delayed_carry_pair(
            target_row, delayed[ell5], {}
        )[12][1]

        # Independent forward-coordinate target recomputation on Sigma_+.
        direct_common = shift_union(row_common, -SHIFT)
        direct_row = restrict_weighted(source_weight, direct_common)
        direct_row = private.old.intersect_weighted_comb(
            direct_row, module.C3, 182, 1, 13
        )
        direct_value = private.delayed_carry_pair(
            direct_row, delayed[ell5], {}
        )[6][1]
        require(target_value == direct_value,
                "pulled and direct target coefficients differ")
        source_vector.append(source_value)
        target_vector.append(target_value)
    return tuple(source_vector), tuple(target_vector)


def multiply_polynomials(left, right):
    result = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            result[i + j] += a * b
    return tuple(result)


def sylvester_resultant(first, second):
    """Exact resultant via a fraction-free Sylvester determinant."""
    while len(first) > 1 and first[-1] == 0:
        first = first[:-1]
    while len(second) > 1 and second[-1] == 0:
        second = second[:-1]
    m = len(first) - 1
    n = len(second) - 1
    f = tuple(reversed(first))
    g = tuple(reversed(second))
    size = m + n
    matrix = [[0] * size for _ in range(size)]
    for row in range(n):
        matrix[row][row:row + m + 1] = f
    for row in range(m):
        matrix[n + row][row:row + n + 1] = g

    sign = 1
    previous = 1
    for pivot_index in range(size - 1):
        pivot_row = next((row for row in range(pivot_index, size)
                          if matrix[row][pivot_index]), None)
        require(pivot_row is not None, "zero Sylvester determinant")
        if pivot_row != pivot_index:
            matrix[pivot_index], matrix[pivot_row] = (
                matrix[pivot_row], matrix[pivot_index]
            )
            sign *= -1
        pivot = matrix[pivot_index][pivot_index]
        for row in range(pivot_index + 1, size):
            for column in range(pivot_index + 1, size):
                numerator = (
                    matrix[row][column] * pivot
                    - matrix[row][pivot_index]
                    * matrix[pivot_index][column]
                )
                require(numerator % previous == 0,
                        "Bareiss division lost integrality")
                matrix[row][column] = numerator // previous
            matrix[row][pivot_index] = 0
        previous = pivot
    return sign * matrix[-1][-1]


def main():
    require((R, T, SHIFT) == (4826809, 297836897838480, 431933040),
            "canonical scales changed")
    module, _prefixes, _a, _b, rails, present, _starts = (
        lift.m.core.build_carrier_data()
    )
    require(module.C3 == 742586 and len(rails) == 162,
            "canonical carrier changed")
    deep_numerator = module.C3 * SHIFT * 182
    require(deep_numerator % T == 0
            and deep_numerator // T % 182 == 14,
            "minimal translation lost its fourteen-unit deep shift")
    require(R * SHIFT == 7 * T,
            "minimal translation lost delayed-phase invariance")

    q_minus = Fraction(47850889647341, 100360982066072)
    q_plus = Fraction(47851035194197, 100360982066072)
    require(q_plus - q_minus == Fraction(7, R),
            "strict adjacent witness translation changed")
    deep_source = Fraction(125553481, 742586)
    deep_target = Fraction(799033, 742586)
    require(169 < deep_source < 181 and 1 < deep_target < 13
            and deep_source - deep_target == 168,
            "strict witness left the two seam charts")

    delayed = marked_prefixes(
        module,
        private.build_pair_prefixes(module),
        two.deepest_fork(module),
    )
    source = two.exclusive_source(module, 3)
    clock_comb = tuple(
        module.make_comb(module.C1, 182, 26 * ell - 13, 26 * ell + 13)
        for ell in range(7)
    )

    source_weight, target_weight, rail_common = rail_data(rails, RAIL)
    static_common = static_semantic_common(
        module, source, clock_comb, rail_common
    )
    rows = []
    for t in range(P):
        common = target_label_common(module, static_common, t)
        if t == 3:
            # Slow direct construction checks the factored scan.
            direct = tuple(two.source_present_section(
                module, source, 1, 0, t, clock_comb
            ))
            direct_target = shift_union(direct, SHIFT)
            require(common == intersect(
                rail_common, intersect(direct, direct_target)
            ), "factored and direct semantic seams differ")
        source_vector, target_vector = clutch_vector(
            module, delayed, present, source_weight, target_weight, common
        )
        source_unit = unit_reduction(source_vector, 12)
        target_unit = unit_reduction(target_vector, 1)
        rows.append((
            t,
            sum(right - left for left, right in common),
            source_vector,
            target_vector,
            source_unit,
            target_unit,
        ))

    expected_vector = (0,) + (AMPLITUDE,) * 6
    zero_vector = (0,) * 7
    require(tuple(row[0] for row in rows if any(row[2]))
            == tuple(range(3, 12)),
            "target-label support changed")
    for t, mass, source_vector, target_vector, source_unit, target_unit in rows:
        expected = expected_vector if 3 <= t <= 11 else zero_vector
        expected_mass = 6320326320 if 3 <= t <= 11 else 0
        require(source_vector == target_vector == expected,
                f"source/target vector changed at t={t}")
        require(mass == expected_mass,
                f"semantic seam mass changed at t={t}")
        if expected == expected_vector:
            require(source_unit == ((9, 0, 0, 0, 0, 0), 1)
                    and target_unit == ((4, 0, 0, 0, 0, 0), 1),
                    f"root unit changed at t={t}")
        else:
            require(source_unit == ((0, 0, 0, 0, 0, 0), 0)
                    and target_unit == ((0, 0, 0, 0, 0, 0), 0),
                    f"zero target label became a unit at t={t}")

    value = AMPLITUDE
    valuation = 0
    while value % P == 0:
        valuation += 1
        value //= P
    require(valuation == 1 and (AMPLITUDE // CONTENT) % P == 9,
            "amplitude content/valuation changed")
    gain = 4 * pow(9, -1, P) % P
    require(gain == 12, "clutch gain changed")

    scalar_profile = tuple(row[2][1] for row in rows)
    target_fourier = tuple(
        two.primitive_fourier_coordinates(scalar_profile, frequency)
        for frequency in range(1, P)
    )
    require(all(any(coordinate for coordinate in row)
                for row in target_fourier),
            "a primitive target character vanished")

    # W=z^3+...+z^11 is a cyclotomic unit with a three-term integral inverse.
    W_profile = tuple(int(3 <= exponent <= 11) for exponent in range(13))
    W = W_profile[:-1]
    V = tuple(int(exponent in (2, 6, 10)) for exponent in range(11))
    phi13 = (1,) * 13
    quotient = (-1, 1, 0, 0, 0, 1, 0, 0, 0, 1)
    left = list(multiply_polynomials(W, V))
    left[0] -= 1
    right = multiply_polynomials(phi13, quotient)
    require(tuple(left) == right,
            "three-section Bezout identity changed")
    resultant = sylvester_resultant(phi13, W)
    require(resultant == 1, "target-window cyclotomic norm changed")
    circular = tuple(
        sum(W_profile[index]
            * (1 if (output - index) % P in (2, 6, 10) else 0)
            for index in range(P))
        for output in range(P)
    )
    require(circular == (3,) + (2,) * 12,
            "three-section circular recombination changed")
    target_det_mod13 = pow(sum(W) % P, 12, P)
    inverse_det_mod13 = pow(sum(V) % P, 12, P)
    require((sum(W) % P, sum(V) % P,
             target_det_mod13, inverse_det_mod13)
            == (9, 3, 1, 1),
            "characteristic-thirteen target unit changed")
    bigraded_determinants = tuple(
        pow(pow(root_scalar, 6, P), 12, P)
        * pow(target_det_mod13, 6, P) % P
        for root_scalar in (9, 4)
    )
    require(bigraded_determinants == (1, 1),
            "formal 72-dimensional tensor unit changed")

    print("LRC14 FULLY MARKED ROOT-ZERO CLUTCH AND TARGET PROFILE")
    print("status=PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT")
    print(f"dependency_pins=LF-normalized count={len(DEPENDENCIES)}")
    print(f"R={R} T={T} tau=7/{R} shift_grid={SHIFT} deep_shift=14")
    print("source_seam=H^R_12_intersect_H^L_0=(169,181)/182")
    print("target_seam=H^R_0_intersect_H^L_1=(1,13)/182")
    print(f"witness={q_minus}->{q_plus} deep={deep_source}->{deep_target}")
    print("retained=two-sided common E3/lawful-(ell,s,t) section + D^6 Q_(3,{1,2}) + relative present + rail8 + h6/kappa1")
    print("t_support=(3,4,5,6,7,8,9,10,11); zeros=(0,1,2,12); common_grid_mass=6320326320")
    print(f"raw_source_target_vector={expected_vector}")
    print(f"C={AMPLITUDE} v13(C)={valuation} (C/26)_mod13={(AMPLITUDE // CONTENT) % P}")
    print("root_profiles=source12:(9,0,0,0,0,0),det1 target1:(4,0,0,0,0,0),det1 gain=12=-1")
    print("target_window=W=z^3+...+z^11 all_primitive_characters=12/12")
    print(f"bezout_inverse=V=z^2+z^6+z^10 resultant_norm={resultant}")
    print(f"positive_three_section_circular_recombination={circular}")
    print("recombination_scope=coefficient-derived lawful t-tables; not three physical translates of one Boolean packet")
    print("char13_target_unit=W(1)=9 det=1 inverse_V(1)=3 det=1")
    print(f"formal_root_times_target_72d_determinants={bigraded_determinants}")
    print("inherited_unmarked_boundary=THM2744 target nonunits rails1,3; eleven gain-1 rails; rail13 gain8")
    print("single_sheet_boundary=the natural U_(0,3) coefficients are unequal (profiles 5 versus 9); no single-sheet naturality")
    print("SCOPE: one two-sided fully marked partial Cech clutch and coefficient target decoder; no global transporter, endpoint current, row exclusion, or LRC14")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
