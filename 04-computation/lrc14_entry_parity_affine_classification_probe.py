#!/usr/bin/env python3
"""Exact probes for the affine 2-adic representation obstruction in LRC(14).

The load-bearing argument is elementary.  This script freezes its boundary
cases, the minority-anchor half-turn clock, and the fixed-pool tests needed by
THM-4326 + THM-4150. Re-referencing changes the distinguished runner;
existence of a convenient reference is a diagnostic and is not physical entry
for a different anchored runner.
"""

from fractions import Fraction
from itertools import combinations
from math import gcd


POOL = {
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63,
    80, 84, 85, 88, 95, 120, 126, 132, 143, 145,
    168, 170, 176, 190, 193, 240, 252, 264, 286, 290,
}
EVEN_DENOMINATOR_BANK = (4, 6, 8, 10, 12, 14)

CHECKS = 0


def check(condition):
    global CHECKS
    CHECKS += 1
    assert condition


def content(velocities):
    values = sorted(velocities)
    check(len(values) == 14 and len(set(values)) == 14)
    g = 0
    for value in values[1:]:
        g = gcd(g, value - values[0])
    check(g > 0)
    return g


def v2(value):
    value = abs(value)
    check(value > 0)
    exponent = 0
    while value % 2 == 0:
        value //= 2
        exponent += 1
    return exponent


def primitive_parity_sizes(velocities):
    values = sorted(velocities)
    g = content(values)
    base = values[0]
    counts = [0, 0]
    for value in values:
        counts[((value - base) // g) % 2] += 1
    check(sum(counts) == 14 and min(counts) > 0)
    return tuple(sorted(counts))


def minimum_valuation_edge_count(velocities):
    values = sorted(velocities)
    s = min(v2(b - a) for a, b in combinations(values, 2))
    edges = sum(v2(b - a) == s for a, b in combinations(values, 2))
    return s, edges


def odd_relative_count(velocities, reference):
    g = content(velocities)
    return sum(
        ((value - reference) // g) % 2 != 0
        for value in velocities
        if value != reference
    )


def distinct_positive_row(velocities, reference):
    g = content(velocities)
    row = [abs(value - reference) // g for value in velocities if value != reference]
    return len(row) == 13 and len(set(row)) == 13


def valid_two_tail_references(velocities):
    return tuple(
        reference
        for reference in sorted(velocities)
        if odd_relative_count(velocities, reference) == 2
        and distinct_positive_row(velocities, reference)
    )


def divisors(value):
    return tuple(d for d in range(1, value + 1) if value % d == 0)


def circle_distance(value, time):
    phase = value * time
    residue = Fraction(phase.numerator % phase.denominator, phase.denominator)
    return min(residue, 1 - residue)


def minority_clock(h, sign):
    check(h >= 1 and sign in (-1, 1))
    return Fraction(1, 2) + sign * Fraction(1, 28 * h)


def minority_even_denominator_escape(h):
    """Return a THM-366 clock denominator, or None on its exact 420-wall."""
    check(h >= 1)
    return next((m for m in EVEN_DENOMINATOR_BANK if (2 * h) % m), None)


def missing_small_denominators(speeds):
    return tuple(
        m for m in range(2, 15) if not any(speed % m == 0 for speed in speeds)
    )


def audit_minority_denominator_sieve():
    audited_odd_residues = 0
    for h in range(1, 1681):
        denominator = minority_even_denominator_escape(h)
        check((denominator is None) == (h % 420 == 0))
        if denominator is None:
            for m in EVEN_DENOMINATOR_BANK:
                check((2 * h) % m == 0)
                check(circle_distance(2 * h, Fraction(1, m)) == 0)
            continue
        time = Fraction(1, denominator)
        check(denominator in EVEN_DENOMINATOR_BANK)
        check((2 * h) % denominator != 0)
        check(circle_distance(2 * h, time) >= Fraction(1, 14))
        for speed in range(1, 2 * denominator, 2):
            check(speed % denominator != 0)
            check(circle_distance(speed, time) >= Fraction(1, 14))
            audited_odd_residues += 1

    # At the exact even-denominator wall, 9, 11, and 13 are the only
    # remaining THM-366 obligations not automatically paid by 2h=840.
    h = 420
    sieve_complete_odds = tuple(range(1, 24, 2))
    missing_thirteen_odds = (1, 3, 5, 7, 9, 11, 15, 17, 19, 21, 23, 25)
    check(missing_small_denominators((2 * h,) + sieve_complete_odds) == ())
    check(
        missing_small_denominators((2 * h,) + missing_thirteen_odds) == (13,)
    )
    check(
        min(
            circle_distance(speed, Fraction(1, 13))
            for speed in (2 * h,) + missing_thirteen_odds
        )
        >= Fraction(1, 14)
    )
    return audited_odd_residues


def audit_minority_clock_bands():
    """Check the exact periodic residue criterion over three full periods."""
    audited_speeds = 0
    for h in range(1, 41):
        plus = minority_clock(h, 1)
        minus = minority_clock(h, -1)
        check(minus == 1 - plus)
        check(circle_distance(2 * h, plus) == Fraction(1, 14))
        check(circle_distance(2 * h, minus) == Fraction(1, 14))
        for speed in range(1, 84 * h, 2):
            residue = speed % (28 * h)
            expected = Fraction(abs(residue - 14 * h), 28 * h)
            safe_residue = residue <= 12 * h or residue >= 16 * h
            plus_distance = circle_distance(speed, plus)
            minus_distance = circle_distance(speed, minus)
            check(residue > 0)
            check(plus_distance == expected)
            check(minus_distance == expected)
            check((plus_distance >= Fraction(1, 14)) == safe_residue)
            check((minus_distance >= Fraction(1, 14)) == safe_residue)
            audited_speeds += 1

    # The first nonvacuous twelve-odd row lies at h=2.
    h = 2
    initial_band = tuple(range(1, 24, 2))
    boundary_fail = initial_band[:-1] + (25,)
    for sign in (-1, 1):
        time = minority_clock(h, sign)
        check(len(initial_band) == len(set(initial_band)) == 12)
        check(max(initial_band) < 12 * h)
        check(circle_distance(2 * h, time) == Fraction(1, 14))
        check(
            min(circle_distance(speed, time) for speed in initial_band)
            == Fraction(1, 14) + Fraction(1, 28 * h)
        )
        check(
            min(circle_distance(speed, time) for speed in boundary_fail)
            == Fraction(1, 14) - Fraction(1, 28 * h)
        )
        check(circle_distance(12 * h + 1, time) < Fraction(1, 14))
        check(circle_distance(16 * h + 1, time) > Fraction(1, 14))
    return audited_speeds


def adaptive_divisor_capacity(speeds, scale):
    capacity = Fraction(0, 1)
    for speed in speeds:
        if speed % scale == 0:
            continue
        order = scale // gcd(scale, speed)
        capacity += Fraction((order + 6) // 7, order)
    return capacity


def audit_minority_residual_atlas():
    unit_checks = 0
    for p in range(8, 15):
        units = tuple(a for a in range(1, 2 * p) if gcd(a, 2 * p) == 1)
        sign_representatives = tuple(a for a in units if a < p)
        for h_residue in range(1, p):
            for a in units:
                check(
                    circle_distance(2 * h_residue, Fraction(a, 2 * p))
                    >= Fraction(1, p)
                )
                for speed in range(1, 2 * p, 2):
                    bad = (
                        circle_distance(speed, Fraction(a, 2 * p))
                        < Fraction(1, 14)
                    )
                    check(
                        bad
                        == (
                            gcd(speed, 2 * p) == 1
                            and (a * speed) % (2 * p) in (1, 2 * p - 1)
                        )
                    )
                    unit_checks += 1
        classes = tuple(
            min(a, 2 * p - a) for a in units if a < 2 * p - a
        )
        check(len(classes) == len(sign_representatives))
        for mask in range(1 << len(classes)):
            represented = {
                classes[index]
                for index in range(len(classes))
                if mask & (1 << index)
            }
            bank_has_witness = any(
                all(
                    (a * speed) % (2 * p) not in (1, 2 * p - 1)
                    for speed in represented
                )
                for a in sign_representatives
            )
            check(bank_has_witness == (len(represented) < len(classes)))

    h = 420
    combined_hostile_odds = (
        1287, 1555, 3235, 4915, 6595, 6825,
        7035, 8275, 8925, 9955, 11025, 11635,
    )
    combined_hostile = (2 * h,) + combined_hostile_odds
    check(len(set(combined_hostile_odds)) == 12)
    check(all(speed % 2 == 1 for speed in combined_hostile_odds))
    check(gcd(*combined_hostile) == 1)
    check(missing_small_denominators(combined_hostile) == ())

    modulus = 28 * h
    grid_anchor_safe = 0
    for j in range(modulus):
        anchor_safe = j % 14 != 0
        if anchor_safe:
            grid_anchor_safe += 1
            check(
                any(
                    12 * h < (j * speed) % modulus < 16 * h
                    for speed in combined_hostile_odds
                )
            )
        else:
            check(
                circle_distance(
                    2 * h, Fraction(1, 2) + Fraction(j, modulus)
                )
                == 0
            )
    check(grid_anchor_safe == 10920)

    base_modulus = 2352
    base_bad_low, base_bad_high = 1008, 1344
    block_d = {311 + 336 * k for k in range(7)}
    block_c = {7 * value for value in (195, 201, 255, 315)}
    base_blockers = block_d | block_c
    base_anchor_safe = tuple(j for j in range(base_modulus) if j % 14 != 0)
    check(len(base_anchor_safe) == 2184)
    check(
        all(
            any(
                base_bad_low < (j * speed) % base_modulus < base_bad_high
                for speed in base_blockers
            )
            for j in base_anchor_safe
        )
    )
    for speed in block_d:
        reopened = sum(
            not any(
                base_bad_low < (j * other) % base_modulus < base_bad_high
                for other in base_blockers - {speed}
            )
            for j in base_anchor_safe
        )
        check(reopened == 156)
    for speed in block_c:
        reopened = sum(
            not any(
                base_bad_low < (j * other) % base_modulus < base_bad_high
                for other in base_blockers - {speed}
            )
            for j in base_anchor_safe
        )
        check(reopened == 36)
    check(
        {5 * speed for speed in base_blockers}
        == set(combined_hostile_odds) - {1287}
    )

    eligible_unit_clocks = 0
    for p in (8, 9, 11, 13):
        check(h % p != 0)
        for a in range(1, 2 * p):
            if gcd(a, 2 * p) != 1:
                continue
            eligible_unit_clocks += 1
            check(
                any(
                    circle_distance(speed, Fraction(a, 2 * p))
                    < Fraction(1, 14)
                    for speed in combined_hostile_odds
                )
            )
    check(eligible_unit_clocks == 36)
    check(
        min(circle_distance(speed, Fraction(7, 19)) for speed in combined_hostile)
        == Fraction(2, 19)
    )

    deep_h = 180180
    deep_odds = (5, 7, 11, 13, 17, 23, 29, 31, 49, 75, 111, 135)
    deep_row = (2 * deep_h,) + deep_odds
    deep_killers = {5, 13, 23, 29, 31, 49, 75, 111, 135}
    check(gcd(*deep_row) == 1)
    check(missing_small_denominators(deep_row) == ())
    retained_by_denominator = []
    retained_clocks = 0
    for denominator in range(2, 41):
        retained_here = 0
        for numerator in range(1, (denominator - 1) // 2 + 1):
            if gcd(numerator, denominator) != 1:
                continue
            time = Fraction(numerator, denominator)
            if circle_distance(2 * deep_h, time) < Fraction(1, 14):
                continue
            retained_here += 1
            retained_clocks += 1
            check(
                any(
                    circle_distance(speed, time) < Fraction(1, 14)
                    for speed in deep_killers
                )
            )
        if retained_here:
            retained_by_denominator.append((denominator, retained_here))
    check(retained_clocks == 112)
    check(
        tuple(retained_by_denominator)
        == (
            (16, 4), (17, 7), (19, 8), (23, 10), (25, 10), (27, 9),
            (29, 12), (31, 13), (32, 8), (34, 7), (37, 16), (38, 8),
        )
    )
    check(
        min(circle_distance(speed, Fraction(1, 41)) for speed in deep_row)
        == Fraction(5, 41)
    )
    represented_scales = sorted(
        set().union(*(set(divisors(speed)) for speed in deep_row)) - {1}
    )
    capacities = tuple(
        (adaptive_divisor_capacity(deep_row, scale), scale)
        for scale in represented_scales
    )
    check(len(represented_scales) == 202)
    check(min(capacities) == (Fraction(10, 7), 7))

    return unit_checks, grid_anchor_safe, eligible_unit_clocks, retained_clocks


def fixed_pool_entries(velocities):
    """Return primitive-literal THM-4326 entries (reference, common scale)."""
    g = content(velocities)
    entries = []
    for reference in valid_two_tail_references(velocities):
        row = [abs(value - reference) // g for value in velocities if value != reference]
        half_body = tuple(sorted(value // 2 for value in row if value % 2 == 0))
        tails = tuple(sorted(value for value in row if value % 2 == 1))
        check(len(half_body) == 11 and len(tails) == 2 and tails[0] != tails[1])
        body_content = 0
        for value in half_body:
            body_content = gcd(body_content, value)
        for scale in divisors(body_content):
            quotient = {value // scale for value in half_body}
            if len(quotient & POOL) == 9 and len(quotient - POOL) == 2:
                entries.append((reference, scale))
    return tuple(entries)


def projective_fixed_pool_entries_at_reference(velocities, reference):
    """Return mass-equivalent rational rescalings that enter THM-4326."""
    if reference not in valid_two_tail_references(velocities):
        return ()
    g = content(velocities)
    row = [abs(value - reference) // g for value in velocities if value != reference]
    half_body = tuple(sorted(value // 2 for value in row if value % 2 == 0))
    candidate_scales = {
        Fraction(pool_label, body_label)
        for pool_label in POOL
        for body_label in half_body
    }
    entries = []
    for scale in sorted(candidate_scales):
        # If lambda=p/q and H=qK, then lambda H=pK. Multiplication by p and
        # q both preserve Haar measure, so G_H and G_(lambda H) have equal
        # mass.  The physical tails are not rescaled.
        image = tuple(scale * value for value in half_body)
        if not all(value.denominator == 1 for value in image):
            continue
        quotient = {int(value) for value in image}
        if len(quotient & POOL) == 9 and len(quotient - POOL) == 2:
            entries.append((scale.numerator, scale.denominator))
    return tuple(entries)


def projective_fixed_pool_entries(velocities):
    return tuple(
        (reference, numerator, denominator)
        for reference in valid_two_tail_references(velocities)
        for numerator, denominator in projective_fixed_pool_entries_at_reference(
            velocities, reference
        )
    )


def audit_configuration(velocities):
    sizes = primitive_parity_sizes(velocities)
    expected_entry = sizes == (2, 12)
    references = valid_two_tail_references(velocities)
    s, edge_count = minimum_valuation_edge_count(velocities)
    check((edge_count == 24) == expected_entry)
    check(bool(references) == expected_entry)
    if expected_entry:
        # All twelve majority-class references have two signed odd relatives;
        # the endpoint argument guarantees at least one collision-free row.
        signed_count = sum(
            odd_relative_count(velocities, reference) == 2
            for reference in velocities
        )
        check(signed_count == 12)
        for reference in references:
            g = content(velocities)
            row = [abs(value - reference) // g for value in velocities if value != reference]
            check(len(set(row)) == 13)
            check(sum(value % 2 == 0 for value in row) == 11)
            check(sum(value % 2 == 1 for value in row) == 2)
    return s, edge_count, references


def main():
    minority_denominator_residues = audit_minority_denominator_sieve()
    minority_clock_speeds = audit_minority_clock_bands()
    (
        minority_unit_checks,
        minority_grid_anchor_safe,
        minority_eligible_unit_clocks,
        minority_deep_retained_clocks,
    ) = audit_minority_residual_atlas()
    profile_table = []
    for minority_size in range(1, 8):
        majority_size = 14 - minority_size
        velocities = tuple(2 * j for j in range(majority_size)) + tuple(
            2 * j + 1 for j in range(minority_size)
        )
        sizes = primitive_parity_sizes(velocities)
        _, edge_count = minimum_valuation_edge_count(velocities)
        reference_counts = sorted(
            {odd_relative_count(velocities, reference) for reference in velocities}
        )
        check(sizes == (minority_size, majority_size))
        check(edge_count == minority_size * majority_size)
        check(reference_counts == sorted({minority_size, majority_size}))
        profile_table.append((minority_size, edge_count))

    # A bounded nonsplit universe checks the negative direction and all
    # reference choices, translations implicit in the difference formulas.
    small_count = 0
    for velocities in combinations(range(-2, 16), 14):
        audit_configuration(velocities)
        small_count += 1
    check(small_count == 3060)

    # Complete structured 12+2 universe: omit one of thirteen even sites and
    # choose every odd pair from twelve sites.
    split_count = 0
    even_sites = tuple(range(-12, 14, 2))
    odd_sites = tuple(range(-11, 12, 2))
    for majority in combinations(even_sites, 12):
        for minority in combinations(odd_sites, 2):
            audit_configuration(majority + minority)
            split_count += 1
    check(split_count == 858)

    ap_hostile = tuple(range(14))
    sharp_one_reference = tuple(range(0, 23, 2)) + (-1, 1)
    parity_pass_pool_fail = tuple(range(0, 23, 2)) + (1, 3)

    positive_body = tuple(sorted(POOL))[:9] + (1, 2)
    pool_positive = (0,) + tuple(2 * value for value in positive_body) + (1, 3)

    even_numerator_body = (1, 2, 4, 5, 8, 10, 15, 20, 21, 30, 40)
    even_numerator_positive = (
        (0,) + tuple(2 * value for value in even_numerator_body) + (1, 3)
    )

    check(audit_configuration(ap_hostile)[1:] == (49, ()))
    check(valid_two_tail_references(sharp_one_reference) == (22,))
    check(fixed_pool_entries(parity_pass_pool_fail) == ())
    check(projective_fixed_pool_entries(parity_pass_pool_fail) == ())
    pool_hit_counts = {
        scale: len({scale * value for value in range(1, 12)} & POOL)
        for scale in range(1, max(POOL) + 1)
    }
    max_pool_hits = max(pool_hit_counts.values())
    max_pool_hit_scales = tuple(
        scale for scale, hits in pool_hit_counts.items() if hits == max_pool_hits
    )
    check((max_pool_hits, max_pool_hit_scales) == (6, (10,)))
    check((0, 1) in fixed_pool_entries(pool_positive))
    check((1, 1) in projective_fixed_pool_entries_at_reference(pool_positive, 0))
    check(
        not any(
            reference == 0
            for reference, _ in fixed_pool_entries(even_numerator_positive)
        )
    )
    check(
        (2, 1)
        in projective_fixed_pool_entries_at_reference(even_numerator_positive, 0)
    )

    # Translation and every nonzero common dilation, including even ones,
    # leave the classification unchanged after affine-content normalization.
    examples = (
        ap_hostile,
        sharp_one_reference,
        parity_pass_pool_fail,
        pool_positive,
        even_numerator_positive,
    )
    affine_checks = 0
    for velocities in examples:
        baseline = primitive_parity_sizes(velocities)
        baseline_edge_count = minimum_valuation_edge_count(velocities)[1]
        baseline_pool_entry = bool(fixed_pool_entries(velocities))
        baseline_projective_entry = bool(projective_fixed_pool_entries(velocities))
        for scale in (-6, -3, -1, 1, 2, 5):
            for shift in (-17, 0, 8):
                image = tuple(scale * value + shift for value in velocities)
                check(primitive_parity_sizes(image) == baseline)
                check(minimum_valuation_edge_count(image)[1] == baseline_edge_count)
                check(bool(fixed_pool_entries(image)) == baseline_pool_entry)
                check(
                    bool(projective_fixed_pool_entries(image))
                    == baseline_projective_entry
                )
                affine_checks += 1

    print(
        "ROOT_PROFILE_MINIMAL_VALUATION_EDGES="
        + ",".join(f"{minority}:{edges}" for minority, edges in profile_table)
    )
    print("SMALL_NONSPLIT_CONFIGURATIONS=3060")
    print("STRUCTURED_SPLIT_12_2_CONFIGURATIONS=858")
    print(f"AFFINE_NORMALIZATION_CHECKS={affine_checks}")
    print("AP_HOSTILE_SPLIT=7+7 MINIMAL_VALUATION_EDGES=49")
    print("ANCHORED_GUARDRAIL=REFERENCE_DEGREE_2_REQUIRED;REREFERENCE_CHANGES_RUNNER")
    print("SHARP_COLLISION_FREE_REFERENCE_LOWER_BOUND=1 REFERENCE=22")
    print(
        "MINORITY_ANCHOR_EVEN_DENOMINATOR_SIEVE=SAFE_UNLESS_420_DIVIDES_H "
        f"ODD_RESIDUES_AUDITED={minority_denominator_residues}"
    )
    print("MINORITY_ANCHOR_420_WALL_REMAINING_SMALL_DENOMINATORS=9,11,13")
    print(
        "MINORITY_ANCHOR_CLOCK=T_PLUS_MINUS_1_OVER_28H "
        f"ODD_SPEEDS_AUDITED={minority_clock_speeds}"
    )
    print("MINORITY_ANCHOR_SAFE_RESIDUES=R_LE_12H_OR_R_GE_16H_MOD_28H")
    print("MINORITY_ANCHOR_INITIAL_BAND_POSITIVE=H_2_W_1_TO_23_ODD")
    print("MINORITY_ANCHOR_FIRST_CLOCK_FAILURE=H_2_W_25_CLEARANCE_3_OVER_56")
    print(
        "MINORITY_ANCHOR_UNIT_BANK_IFF_CHECKS="
        f"{minority_unit_checks} P_RANGE=8_TO_14"
    )
    print(
        "MINORITY_ANCHOR_COMBINED_BANK_HOSTILE=H_420 "
        f"GRID_ANCHOR_SAFE={minority_grid_anchor_safe} GRID_WITNESSES=0 "
        f"ELIGIBLE_UNIT_CLOCKS={minority_eligible_unit_clocks} "
        "UNIT_WITNESSES=0 SAFE_CONTROL=7_OVER_19 MIN_CLEARANCE=2_OVER_19"
    )
    print(
        "MINORITY_ANCHOR_DEEP_HOSTILE=H_180180 "
        f"DENOM_LE_40_ANCHOR_SAFE={minority_deep_retained_clocks} "
        "WITNESSES=0 ADAPTIVE_SCALES=202 MIN_CAPACITY=10_OVER_7_AT_7 "
        "FIRST_WITNESS=1_OVER_41 MIN_CLEARANCE=5_OVER_41"
    )
    print("PARITY_PASS_POOL_FAIL_VALID_REFS=0,22 FIXED_POOL_ENTRIES=0")
    print("PARITY_PASS_PROJECTIVE_POOL_MAX_HITS=6 AT_SCALE=10")
    print("FIXED_POOL_POSITIVE_ENTRY=REFERENCE_0_SCALE_1")
    print("EVEN_NUMERATOR_PROJECTIVE_ENTRY=REFERENCE_0_SCALE_2_OVER_1")
    print(f"CHECKS={CHECKS}")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
