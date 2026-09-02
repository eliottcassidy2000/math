#!/usr/bin/env python3
"""Exact probes for the affine 2-adic representation obstruction in LRC(14).

The load-bearing argument is elementary.  This script freezes its boundary
cases and the fixed-pool tests needed by THM-4326 + THM-4150. Re-referencing
changes the distinguished runner; existence of a convenient reference is a
diagnostic and is not physical entry for a different anchored runner.
"""

from fractions import Fraction
from itertools import combinations
from math import gcd


POOL = {
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63,
    80, 84, 85, 88, 95, 120, 126, 132, 143, 145,
    168, 170, 176, 190, 193, 240, 252, 264, 286, 290,
}

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
    print("PARITY_PASS_POOL_FAIL_VALID_REFS=0,22 FIXED_POOL_ENTRIES=0")
    print("PARITY_PASS_PROJECTIVE_POOL_MAX_HITS=6 AT_SCALE=10")
    print("FIXED_POOL_POSITIVE_ENTRY=REFERENCE_0_SCALE_1")
    print("EVEN_NUMERATOR_PROJECTIVE_ENTRY=REFERENCE_0_SCALE_2_OVER_1")
    print(f"CHECKS={CHECKS}")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
