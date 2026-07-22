#!/usr/bin/env python3
"""Exact unit-annulus certificate for the LRC(14) scalar 6+1 tail.

Suppose H and q_1,...,q_6 are units modulo 13 and the seventh tooth is
13v.  If s=v_13(v) and N=13^(s+2), sampling t=z/N and multiplying z by
H modulo N reduces a putative cover to six multiplicative translates of

    B_N={z in (Z/NZ)^*: ||z/N|| < 1/14}

covering

    A_N={z in (Z/NZ)^*: ||z/N|| > 1/7}.

This script excludes s=0,1,2,3 exactly.  At N=169 it supplies both a
zero-graph/component certificate and an independent marginal-gain branch
and bound.  At N=13^3,13^4,13^5 it scans every unit sign class using an
integer Euclidean floor-sum evaluator.  No floating point or assert is used.
"""

from collections import Counter
from itertools import combinations


def require(condition, message):
    """Runtime check that remains active under python -O."""
    if not condition:
        raise RuntimeError(message)


def norm_numerator(x, modulus):
    """Return modulus * ||x/modulus|| for an integer residue x."""
    residue = x % modulus
    return min(residue, modulus - residue)


def unit_guard_universe(modulus):
    """The exact finite unit annulus A_N, in standard representatives."""
    return tuple(
        z
        for z in range(modulus)
        if z % 13 != 0 and 7 * norm_numerator(z, modulus) > modulus
    )


def sign_class_labels(modulus):
    """Represent (Z/NZ)^*/{+/-1} by 1 <= a < N/2."""
    return tuple(a for a in range(1, (modulus + 1) // 2) if a % 13 != 0)


def direct_mask(modulus, universe, multiplier):
    """Return the exact terminal mask as a Python-integer bitset."""
    mask = 0
    for index, z in enumerate(universe):
        if 14 * norm_numerator(multiplier * z, modulus) < modulus:
            mask |= 1 << index
    return mask


def floor_sum(n, modulus, multiplier, offset):
    """Sum floor((multiplier*x+offset)/modulus), 0 <= x < n.

    This is the Euclidean algorithm usually called ``floor_sum``.  All
    inputs here are nonnegative.  Its logarithmic recurrence is exact.
    """
    require(n >= 0 and modulus > 0 and multiplier >= 0 and offset >= 0,
            "floor_sum received an invalid input")
    answer = 0
    while True:
        if multiplier >= modulus:
            answer += (n - 1) * n * (multiplier // modulus) // 2
            multiplier %= modulus
        if offset >= modulus:
            answer += n * (offset // modulus)
            offset %= modulus
        height = multiplier * n + offset
        if height < modulus:
            return answer
        n = height // modulus
        offset = height % modulus
        modulus, multiplier = multiplier, modulus


def floor_sum_range(modulus, multiplier, first, last, offset):
    """Sum floor((multiplier*z+offset)/modulus), first <= z <= last."""
    if last < first:
        return 0
    a = multiplier % modulus
    return (
        floor_sum(last + 1, modulus, a, offset)
        - floor_sum(first, modulus, a, offset)
    )


def full_grid_guard_safe_tooth_count(modulus, multiplier):
    """Count guard-safe z with multiplier*z terminal-dangerous on Z/NZ.

    Write L=floor(N/7), R=floor((N-1)/14).  Guard-safe representatives are
    L+1,...,N-L-1.  If r=multiplier*z mod N, terminal danger is the disjoint
    alternative r<R+1 or r>=N-R.  Each alternative is a difference of two
    floor sums, giving the formula below.
    """
    left = modulus // 7 + 1
    right = modulus - modulus // 7 - 1
    radius = (modulus - 1) // 14
    threshold = radius + 1
    length = right - left + 1

    base = floor_sum_range(modulus, multiplier, left, right, 0)
    at_least_threshold = (
        floor_sum_range(
            modulus, multiplier, left, right, modulus - threshold
        )
        - base
    )
    below_threshold = length - at_least_threshold
    at_least_modulus_minus_radius = (
        floor_sum_range(modulus, multiplier, left, right, radius) - base
    )
    return below_threshold + at_least_modulus_minus_radius


def annular_mask_count(depth, multiplier):
    """Count A_N intersect multiplier^(-1) B_N for N=13^depth."""
    modulus = 13 ** depth
    quotient = modulus // 13
    # Multiples z=13w contribute the same full-grid count at modulus N/13.
    return (
        full_grid_guard_safe_tooth_count(modulus, multiplier)
        - full_grid_guard_safe_tooth_count(quotient, multiplier % quotient)
    )


def direct_annular_mask_count(depth, multiplier):
    """Independent definition-level count, used as a hostile cross-check."""
    modulus = 13 ** depth
    total = 0
    for z in range(modulus):
        if z % 13 == 0:
            continue
        if 7 * norm_numerator(z, modulus) <= modulus:
            continue
        if 14 * norm_numerator(multiplier * z, modulus) < modulus:
            total += 1
    return total


def expected_universe_size(depth):
    """Closed count of A_(13^depth); parity uses 13=-1 modulo 14."""
    modulus = 13 ** depth
    parity = 1 if depth % 2 == 0 else -1
    numerator = 60 * modulus - 130 * parity
    require(numerator % 91 == 0, "universe formula is not integral")
    return numerator // 91


def validate_floor_sum():
    """Compare the Euclidean evaluator with its literal definition."""
    state = 20260722
    for case in range(512):
        state = (1103515245 * state + 12345) % (2 ** 31)
        modulus = 2 + state % 96
        state = (1103515245 * state + 12345) % (2 ** 31)
        n = state % 101
        state = (1103515245 * state + 12345) % (2 ** 31)
        multiplier = state % 501
        state = (1103515245 * state + 12345) % (2 ** 31)
        offset = state % 501
        exact = sum(
            (multiplier * x + offset) // modulus for x in range(n)
        )
        got = floor_sum(n, modulus, multiplier, offset)
        require(got == exact, f"floor_sum mismatch in hostile case {case}")


def connected(vertices, masks):
    """Connectivity in the graph of positive pair intersections."""
    if not vertices:
        return False
    seen = {vertices[0]}
    stack = [vertices[0]]
    while stack:
        i = stack.pop()
        for j in vertices:
            if j not in seen and masks[i] & masks[j]:
                seen.add(j)
                stack.append(j)
    return len(seen) == len(vertices)


def mod_169_certificate():
    """Return the compact component certificate and independent optimum."""
    modulus = 169
    universe = unit_guard_universe(modulus)
    labels = sign_class_labels(modulus)
    masks = tuple(direct_mask(modulus, universe, a) for a in labels)
    sizes = tuple(mask.bit_count() for mask in masks)

    require(len(universe) == 110, "bad mod-169 universe size")
    require(len(labels) == 78, "bad mod-169 sign-class count")
    size_histogram = Counter(sizes)
    expected_histogram = Counter({0: 2, 8: 1, 12: 3, 16: 21, 18: 40, 20: 11})
    require(size_histogram == expected_histogram, "bad mod-169 mask histogram")

    pair_histogram = Counter(
        (masks[i] & masks[j]).bit_count()
        for i, j in combinations(range(len(masks)), 2)
    )
    require(
        pair_histogram == Counter({0: 632, 2: 1508, 4: 613, 6: 148, 8: 58, 10: 44}),
        "bad mod-169 pair histogram",
    )

    # A six-cover must have total incidence at least 110.  Since five masks
    # contribute at most 100, no mask of size <=8 can occur.
    relevant = tuple(i for i, size in enumerate(sizes) if size >= 12)
    require(len(relevant) == 75, "bad relevant-mask count")
    zero_neighbors = {
        i: {
            j
            for j in relevant
            if j != i and masks[i] & masks[j] == 0
        }
        for i in relevant
    }

    component_rows = []
    expected_rows = {
        1: (72, 2172, 96, (16, 17, 30, 31, 46, 61)),
        2: (345, 35, 92, (31, 46, 16, 30, 47, 61)),
        3: (608, 0, 0, ()),
    }
    for component_size in (1, 2, 3):
        seed_count = 0
        completion_count = 0
        best_union = 0
        best_labels = ()
        for component in combinations(relevant, component_size):
            if component_size >= 2 and not connected(component, masks):
                continue
            common_zero = set(relevant) - set(component)
            for i in component:
                common_zero &= zero_neighbors[i]
            needed = 6 - component_size
            if len(common_zero) < needed:
                continue
            seed_count += 1
            component_sum = sum(sizes[i] for i in component)
            component_union = 0
            for i in component:
                component_union |= masks[i]
            for outside in combinations(sorted(common_zero), needed):
                if component_sum + sum(sizes[j] for j in outside) < 110:
                    continue
                completion_count += 1
                union = component_union
                for j in outside:
                    union |= masks[j]
                union_size = union.bit_count()
                if union_size > best_union:
                    best_union = union_size
                    best_labels = tuple(labels[j] for j in component + outside)
        row = (seed_count, completion_count, best_union, best_labels)
        require(row == expected_rows[component_size],
                f"bad component row {component_size}: {row}")
        component_rows.append((component_size,) + row)

    # If the positive graph of a hypothetical cover is connected, Hunter's
    # tree bound and even intersections leave only six size-20 masks with all
    # positive intersections equal to two.  Check that exceptional residue.
    size_20 = tuple(i for i in relevant if sizes[i] == 20)
    require(len(size_20) == 11, "bad size-20 mask count")
    exceptional_connected = 0
    for packet in combinations(size_20, 6):
        if any(
            (masks[i] & masks[j]).bit_count() > 2
            for i, j in combinations(packet, 2)
        ):
            continue
        if connected(packet, masks):
            exceptional_connected += 1
    require(exceptional_connected == 0,
            "connected Hunter equality packet unexpectedly exists")

    # Independent branch and bound.  The sum of the largest remaining
    # individual marginal gains is a rigorous (overlap-ignoring) upper bound.
    best_size = 0
    best_packet = ()
    nodes = 0
    leaves = 0
    prunes = 0

    def search(start, chosen, union):
        nonlocal best_size, best_packet, nodes, leaves, prunes
        nodes += 1
        needed = 6 - len(chosen)
        if needed == 0:
            leaves += 1
            union_size = union.bit_count()
            if union_size > best_size:
                best_size = union_size
                best_packet = chosen
            return
        gains = [
            ((masks[j] & ~union).bit_count(), j)
            for j in range(start, len(masks))
        ]
        if len(gains) < needed:
            return
        optimistic = union.bit_count() + sum(
            sorted((gain for gain, _ in gains), reverse=True)[:needed]
        )
        if optimistic <= best_size:
            prunes += 1
            return
        for _, j in gains:
            search(j + 1, chosen + (j,), union | masks[j])

    search(0, (), 0)
    best_labels = tuple(labels[j] for j in best_packet)
    best_union_mask = 0
    for j in best_packet:
        best_union_mask |= masks[j]
    uncovered = tuple(
        z
        for index, z in enumerate(universe)
        if (best_union_mask >> index) & 1 == 0
    )
    require(best_size == 96, "bad independent six-union optimum")
    require(best_labels == (14, 19, 27, 29, 63, 77),
            "bad independent maximizing packet")
    require(
        uncovered == (28, 30, 34, 40, 42, 74, 77, 92, 95, 127, 129, 135, 139, 141),
        "bad independent uncovered witness",
    )
    require((nodes, leaves, prunes) == (95680, 407, 87118),
            "bad branch-and-bound telemetry")

    positive_control = (3, 4, 5, 7, 8, 11, 18, 75, 76)
    mask_by_label = {label: mask for label, mask in zip(labels, masks)}
    positive_union = 0
    for label in positive_control:
        positive_union |= mask_by_label[label]
    require(positive_union.bit_count() == 110,
            "nine-mask positive cover control failed")

    return {
        "universe": len(universe),
        "classes": len(labels),
        "size_histogram": tuple(sorted(size_histogram.items())),
        "pair_histogram": tuple(sorted(pair_histogram.items())),
        "component_rows": tuple(component_rows),
        "connected_candidates": exceptional_connected,
        "best_size": best_size,
        "best_labels": best_labels,
        "uncovered": uncovered,
        "nodes": nodes,
        "leaves": leaves,
        "prunes": prunes,
        "positive_control": positive_control,
    }


def scan_depth(depth):
    """Scan every sign class at N=13^depth by exact floor sums."""
    modulus = 13 ** depth
    labels = sign_class_labels(modulus)
    expected_classes = 6 * 13 ** (depth - 1)
    require(len(labels) == expected_classes, f"bad class count at depth {depth}")

    direct_universe = len(unit_guard_universe(modulus))
    formula_universe = expected_universe_size(depth)
    require(direct_universe == formula_universe,
            f"bad universe formula at depth {depth}")

    counts = []
    for label in labels:
        counts.append((annular_mask_count(depth, label), label))
    counts.sort(reverse=True)
    maximum = counts[0][0]
    maximizers = tuple(sorted(label for count, label in counts if count == maximum))
    second = max(count for count, _ in counts if count < maximum)
    runners_up = tuple(sorted(label for count, label in counts if count == second))

    expected = {
        2: (110, 20, (5, 6, 14, 17, 19, 27, 41, 46, 58, 75, 82), 18),
        3: (1450, 240, (6, 183), 236),
        4: (18830, 3140, (6, 2380), 3084),
        5: (244810, 40800, (6, 30941), 40058),
    }
    expected_universe, expected_maximum, expected_maximizers, expected_second = expected[depth]
    require(formula_universe == expected_universe,
            f"unexpected universe at depth {depth}")
    require(maximum == expected_maximum,
            f"unexpected maximum at depth {depth}")
    require(maximizers == expected_maximizers,
            f"unexpected maximizers at depth {depth}")
    require(second == expected_second,
            f"unexpected runner-up at depth {depth}")

    # Definition-level cross-checks.  Depths two and three are cheap enough
    # to check exhaustively.  At depths four and five check every extremizer,
    # every runner-up, and a deterministic hostile sample of other labels.
    if depth <= 3:
        check_labels = labels
    else:
        sample = {5, 6, (modulus - 1) // 12}
        sample.update(maximizers)
        sample.update(runners_up)
        state = 13 ** depth + 20260722
        while len(sample) < 20:
            state = (1664525 * state + 1013904223) % modulus
            candidate = min(state, modulus - state)
            if 0 < candidate < (modulus + 1) // 2 and candidate % 13 != 0:
                sample.add(candidate)
        check_labels = tuple(sorted(sample))
    for label in check_labels:
        direct = direct_annular_mask_count(depth, label)
        fast = annular_mask_count(depth, label)
        require(direct == fast,
                f"direct/floor-sum mismatch at depth {depth}, label {label}")

    return {
        "depth": depth,
        "modulus": modulus,
        "universe": formula_universe,
        "classes": len(labels),
        "maximum": maximum,
        "maximizers": maximizers,
        "second": second,
        "runners_up": runners_up,
        "direct_checks": len(check_labels),
    }


def main():
    validate_floor_sum()
    certificate = mod_169_certificate()
    scans = tuple(scan_depth(depth) for depth in (2, 3, 4, 5))

    # The depth-two scan must agree with the stronger union certificate.
    require(scans[0]["maximum"] == 20, "depth-two scan mismatch")
    require(certificate["best_size"] < certificate["universe"],
            "depth-two obstruction failed")

    # Depth three and five close by raw capacity.
    require(6 * scans[1]["maximum"] == scans[1]["universe"] - 10,
            "depth-three capacity gap failed")
    require(6 * scans[3]["maximum"] == scans[3]["universe"] - 10,
            "depth-five capacity gap failed")

    # At depth four one nonmaximal mask already loses too much capacity.  If
    # all six are maximal, only two distinct masks occur and their union is
    # tiny; repetitions cannot enlarge it.
    depth_four = scans[2]
    require(
        5 * depth_four["maximum"] + depth_four["second"]
        == 18784 < depth_four["universe"],
        "depth-four runner-up gap failed",
    )
    modulus_four = depth_four["modulus"]
    universe_four = unit_guard_universe(modulus_four)
    maximum_union = 0
    for label in depth_four["maximizers"]:
        maximum_union |= direct_mask(modulus_four, universe_four, label)
    require(maximum_union.bit_count() == 5758,
            "depth-four extremizer union changed")

    print("LRC14 SIX-PLUS-ONE UNIT-ANNULUS DEPTH-4 CERTIFICATE")
    print("arithmetic: exact integers only; no floats; no assert dependence")
    print("floor_sum hostile definition checks: 512/512")
    print()
    print("mod 169 base certificate")
    print(f"  universe={certificate['universe']} sign_classes={certificate['classes']}")
    print(f"  mask_size_histogram={certificate['size_histogram']}")
    print(f"  pair_intersection_histogram={certificate['pair_histogram']}")
    print("  smallest_component_size seeds admissible_completions max_union witness")
    for row in certificate["component_rows"]:
        size, seeds, completions, best, witness = row
        print(f"    {size}: {seeds} {completions} {best} {witness}")
    print("  connected size-20/all-edge-2 exceptional packets: 0 of 462")
    print(
        "  independent branch_bound: "
        f"nodes={certificate['nodes']} leaves={certificate['leaves']} "
        f"prunes={certificate['prunes']} max_union={certificate['best_size']}"
    )
    print(f"  maximizing labels={certificate['best_labels']}")
    print(f"  uncovered residues={certificate['uncovered']}")
    print(f"  positive nine-mask cover={certificate['positive_control']}")
    print()
    print("unit-annulus exhaustive floor-sum census")
    print("  m N |U_N| sign_classes max maximizers second runners_up direct_checks")
    for scan in scans:
        print(
            f"  {scan['depth']} {scan['modulus']} {scan['universe']} "
            f"{scan['classes']} {scan['maximum']} {scan['maximizers']} "
            f"{scan['second']} {scan['runners_up']} {scan['direct_checks']}"
        )
    print()
    print("capacity consequences")
    print("  m=2: exact maximum union 96 < 110")
    print("  m=3: 6*240 = 1440 < 1450")
    print("  m=4: nonmax sum <= 18784 < 18830; max-mask union = 5758")
    print("  m=5: 6*40800 = 244800 < 244810")
    print("CONCLUSION: v_13(v) in {0,1,2,3} is impossible; hence 13^4 divides v.")
    print()
    print("reproduce:")
    print("  python3 04-computation/lrc14_six_plus_one_unit_annulus_depth4_codex_20260722.py")
    print("  python3 -O 04-computation/lrc14_six_plus_one_unit_annulus_depth4_codex_20260722.py")


if __name__ == "__main__":
    main()
