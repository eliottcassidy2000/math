#!/usr/bin/env python3
"""Exact referee for the r=6 harmonic tail and finite-box reduction.

The order-sensitive max-T extrapolation is false (THM-1134).  This verifier
instead keeps the seven-speed core safe set E intact and charges each of the
six killer danger sets by its exact one-interval discrepancy.  For

    D_k = {t in R/Z : ||k t|| < 1/14},

the centered primitive of 1_{D_k}-1/7 has oscillation 6/(49k).  Hence, if E
has measure mu and C circular interval components,

    measure(E intersect D_k) <= mu/7 + 6C/(49k).

The resulting sufficient condition for six killers is

    sum_i 1/k_i < 7 mu/(6C).

All 792 seven-subsets of {1,...,12} are checked exactly below.  A second,
paper-level interval argument uses settled LRC for every proper prefix to
show that a jump

    k_{s+1}/k_s > 6(8+s)/(s+1)

absorbs all 6-s remaining killers at once.  Together these statements turn
the unbounded branch into a finite (but still unverified) box.

Tournament Analysis is deliberately diagnostic.  On killer danger-set
vertices, orient by harmonic loss 1/k and break ties by the integer label.
This is transitive and retains the scalar used by the proof, but destroys
which core intervals each danger set hits and all higher intersections.  The
proof-facing object is the component--danger-set loss ledger.
"""

from fractions import Fraction as F
from itertools import combinations


CORE_UNIVERSE = tuple(range(1, 13))
EXPLICIT_CORE = (1, 2, 4, 7, 9, 11, 12)
EXPLICIT_FIRST_FIVE = (290, 292, 294, 296, 298)
EXPLICIT_SIXTH = 338
EXPLICIT_SIMPLE_WITNESS = F(5, 16)
EXPLICIT_SUPPLIED_WITNESS = F(106, 303)


def require(condition: bool, message: str) -> None:
    """Proof-bearing check that remains active under ``python -O``."""
    if not condition:
        raise RuntimeError(message)


def norm_distance(x: F) -> F:
    residue = x % 1
    return min(residue, 1 - residue)


def safe_set(speeds: tuple[int, ...]) -> list[tuple[F, F]]:
    """Positive-length components of {t in [0,1] : ||v t|| >= 1/14}."""
    breakpoints = {F(0), F(1)}
    for speed in speeds:
        radius = F(1, 14 * speed)
        for j in range(speed + 1):
            for endpoint in (F(j, speed) - radius, F(j, speed) + radius):
                if 0 <= endpoint <= 1:
                    breakpoints.add(endpoint)

    ordered = sorted(breakpoints)
    result: list[tuple[F, F]] = []
    for left, right in zip(ordered, ordered[1:]):
        midpoint = (left + right) / 2
        if all(norm_distance(speed * midpoint) >= F(1, 14) for speed in speeds):
            if result and result[-1][1] == left:
                result[-1] = (result[-1][0], right)
            else:
                result.append((left, right))
    return result


def safe_set_by_bad_union(speeds: tuple[int, ...]) -> list[tuple[F, F]]:
    """Independent construction: union all clipped danger teeth, then complement."""
    bad: list[tuple[F, F]] = []
    for speed in speeds:
        radius = F(1, 14 * speed)
        for j in range(speed + 1):
            left = max(F(0), F(j, speed) - radius)
            right = min(F(1), F(j, speed) + radius)
            if left < right:
                bad.append((left, right))
    bad.sort()

    merged: list[tuple[F, F]] = []
    for left, right in bad:
        if merged and left <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], right))
        else:
            merged.append((left, right))

    result: list[tuple[F, F]] = []
    cursor = F(0)
    for left, right in merged:
        if cursor < left:
            result.append((cursor, left))
        cursor = max(cursor, right)
    if cursor < 1:
        result.append((cursor, F(1)))
    return result


def remove_bad(intervals: list[tuple[F, F]], speed: int) -> list[tuple[F, F]]:
    """Subtract {t : ||speed*t|| < 1/14} with exact rational endpoints."""
    result: list[tuple[F, F]] = []
    for left, right in intervals:
        first = (left * speed).__floor__() - 1
        last = (right * speed).__floor__() + 1
        cursor = left
        for j in range(first, last + 1):
            bad_left = F(14 * j - 1, 14 * speed)
            bad_right = F(14 * j + 1, 14 * speed)
            if bad_right <= left or bad_left >= right:
                continue
            bad_left = max(left, bad_left)
            bad_right = min(right, bad_right)
            if cursor < bad_left:
                result.append((cursor, bad_left))
            cursor = max(cursor, bad_right)
        if cursor < right:
            result.append((cursor, right))
    return result


def interval_metrics(intervals: list[tuple[F, F]]) -> tuple[int, F, F]:
    require(bool(intervals), "the safe interval list must be nonempty")
    measure = sum((right - left for left, right in intervals), F(0))
    largest = max(right - left for left, right in intervals)
    return len(intervals), measure, largest


def measure_threshold(core: tuple[int, ...], killers: tuple[int, ...]) -> tuple[F, str, int, F, F]:
    intervals = safe_set(core)
    for killer in killers:
        intervals = remove_bad(intervals, killer)
    components, measure, largest = interval_metrics(intervals)
    count_threshold = F(components, 6) / measure
    component_threshold = F(1, 3) / largest
    if count_threshold <= component_threshold:
        return count_threshold, "count", components, measure, largest
    return component_threshold, "component", components, measure, largest


def covers_two_through_fourteen(speeds: tuple[int, ...]) -> bool:
    return all(any(speed % q == 0 for speed in speeds) for q in range(2, 15))


def main() -> None:
    print("### exact seven-speed core ledger ###")
    rows = []
    for core in combinations(CORE_UNIVERSE, 7):
        intervals = safe_set(core)
        require(
            intervals == safe_set_by_bad_union(core),
            f"independent safe-set reconstructions disagree for core {core}",
        )
        components, measure, largest = interval_metrics(intervals)
        budget = F(7) * measure / (6 * components)
        first_killer = 13 * max(core) + 1
        require(
            budget > F(1, first_killer),
            f"nonpositive second-coordinate denominator for core {core}",
        )
        rows.append((budget, core, components, measure, largest, first_killer))

    require(len(rows) == 792, "the seven-speed core census must contain 792 rows")
    worst_budget = min(rows, key=lambda row: row[0])
    budget, core, components, measure, largest, first_killer = worst_budget
    require(core == (1, 2, 7, 9, 10, 11, 12), "unexpected worst harmonic core")
    require(components == 26, "unexpected worst-core component count")
    require(measure == F(12559, 48510), "unexpected worst-core measure")
    require(budget == F(12559, 1081080), "unexpected worst harmonic budget")
    # Ordered killers are distinct integers.  If k1>=K then
    # sum_i 1/k_i <= sum_{j=K}^{K+5} 1/j, not merely 6/K.
    first_tail_514 = sum((F(1, j) for j in range(514, 520)), F(0))
    first_tail_513 = sum((F(1, j) for j in range(513, 519)), F(0))
    require(
        budget - first_tail_514
        == F(114421278347, 370205035514594280) > 0,
        "the distinct k1=514 harmonic comparison must be strictly positive",
    )
    require(
        budget - first_tail_513 <= 0,
        "k1=513 should be the first value not dispatched by this uniform estimate",
    )

    # If k2>=951, then k2,...,k6 are at least 951,...,955.  The remaining
    # k1 is at least 13*max(P)+1.  Check the exact margin on every core.
    second_tail_951 = sum((F(1, j) for j in range(951, 956)), F(0))
    second_tail_950 = sum((F(1, j) for j in range(950, 955)), F(0))
    second_rows = [
        (row[0] - F(1, row[5]) - second_tail_951, row) for row in rows
    ]
    second_margin, worst_second = min(second_rows, key=lambda item: item[0])
    require(worst_second[1] == core, "unexpected worst core for the k2 escape")
    require(worst_second[5] == 157, "unexpected minimum first killer at the k2 escape")
    require(
        second_margin == F(4670546220583, 4412023437154312980) > 0,
        "the distinct k2=951 harmonic comparison must be strictly positive",
    )
    require(
        budget - F(1, 157) - second_tail_950 <= 0,
        "k2=950 should be the first value not dispatched by this uniform estimate",
    )

    print(f"cores checked exactly: {len(rows)}")
    print("worst core: " + " ".join(map(str, core)))
    print(f"components C: {components}")
    print(f"measure mu: {measure}")
    print(f"harmonic budget 7mu/(6C): {budget}")
    print(f"strict distinct all-large check at k1=514: {budget - first_tail_514}")
    print(f"k1=513 estimate margin (nonpositive): {budget - first_tail_513}")
    print("CONCLUSION 1: k1 >= 514 implies a surviving core-safe time")
    print(f"strict distinct k2=951 margin: {second_margin}")
    print(f"k2=950 estimate margin (nonpositive): {budget - F(1, 157) - second_tail_950}")
    print("CONCLUSION 2: every harmonic-uncertified tuple has k1 <= 513 and k2 <= 950")
    print()

    print("### proper-prefix interval ratios and finite box ###")
    ratios = []
    for s in range(1, 6):
        # A prefix has 7+s speeds.  Settled LRC gives clearance 1/(8+s).
        # Its margin above 1/14 and max-speed Lipschitzness give a safe
        # interval of length 2*(1/(8+s)-1/14)/k_s.  Applying the same exact
        # discrepancy to all 6-s remaining killers gives the ratio below.
        ratio = F(6 * (8 + s), s + 1)
        ratios.append(ratio)
        print(f"s={s}: k_(s+1)/k_s > {ratio} dispatches all {6-s} remaining killers")
    require(
        ratios == [F(27), F(20), F(33, 2), F(72, 5), F(13)],
        "unexpected proper-prefix ratio ladder",
    )

    # The s=1 ratio is superseded by the stronger all-core harmonic k2 cap.
    finite_box = [513, 950]
    for s in range(2, 6):
        finite_box.append((F(6 * (8 + s), s + 1) * finite_box[-1]).__floor__())
    require(
        finite_box == [513, 950, 19000, 313500, 4514400, 58687200],
        "unexpected recursive finite box",
    )
    print("finite residual box (not an exhaustive verification):")
    for index, bound in enumerate(finite_box, 1):
        print(f"  k{index} <= {bound}")
    print("CONCLUSION 3: outside this box, the harmonic or proper-prefix interval horn proves LRC")
    print("GUARDRAIL: no claim is made that the tuples inside this finite box have been enumerated")
    print()

    print("### the explicit covering max-T failure is arithmetically easy ###")
    speeds = EXPLICIT_CORE + EXPLICIT_FIRST_FIVE + (EXPLICIT_SIXTH,)
    require(
        len(speeds) == 13 and len(set(speeds)) == 13,
        "the displayed family must have thirteen distinct speeds",
    )
    require(covers_two_through_fourteen(speeds), "the displayed family must cover 2..14")
    threshold, branch, residual_components, residual_measure, residual_largest = measure_threshold(
        EXPLICIT_CORE, EXPLICIT_FIRST_FIVE
    )
    require(
        threshold == F(1043, 3) > EXPLICIT_SIXTH,
        "the exact max-T certificate must fail at sixth killer 338",
    )

    simple_distances = [
        (speed, norm_distance(speed * EXPLICIT_SIMPLE_WITNESS)) for speed in speeds
    ]
    simple_minimum = min(distance for _, distance in simple_distances)
    simple_owners = [speed for speed, distance in simple_distances if distance == simple_minimum]
    require(
        simple_minimum == F(1, 8) > F(1, 14) and simple_owners == [294, 298],
        "the simple 5/16 witness certificate changed",
    )

    supplied_distances = [
        (speed, norm_distance(speed * EXPLICIT_SUPPLIED_WITNESS)) for speed in speeds
    ]
    supplied_minimum = min(distance for _, distance in supplied_distances)
    supplied_owners = [
        speed for speed, distance in supplied_distances if distance == supplied_minimum
    ]
    require(
        supplied_minimum == F(15, 101) > F(1, 14) and supplied_owners == [9, 294],
        "the supplied 106/303 witness certificate changed",
    )
    print("speeds: " + " ".join(map(str, speeds)))
    print("covering 2..14: yes")
    print(
        f"first-five residual: N={residual_components} mu={residual_measure} "
        f"L={residual_largest} T={threshold} ({branch})"
    )
    print(f"sixth killer {EXPLICIT_SIXTH} < T: the max-T measure horn fails")
    print(
        f"simple exact witness t={EXPLICIT_SIMPLE_WITNESS}: "
        f"minimum distance {simple_minimum}, owners {simple_owners}"
    )
    print(
        f"supplied exact witness t={EXPLICIT_SUPPLIED_WITNESS}: "
        f"minimum distance {supplied_minimum}, owners {supplied_owners}"
    )
    print("CONCLUSION 4: this is a certificate failure, not a hard LRC family")
    print()

    print("### Tournament Analysis (lossy diagnostic) ###")
    killers = EXPLICIT_FIRST_FIVE + (EXPLICIT_SIXTH,)
    path = sorted(killers, key=lambda killer: (F(1, killer), killer), reverse=True)
    print("vertices: the six killer danger sets")
    print("observable: harmonic loss 1/k; tie gauge: increasing integer k")
    print("Hamiltonian path: " + " -> ".join(map(str, path)))
    print("score histogram: 0 1 2 3 4 5")
    print("directed cycles: 0; SCCs: 6 singleton components; Hamiltonian paths: 1")
    print("preserved: the scalar harmonic loss used in the sufficient inequality")
    print("destroyed: interval owners, phases, and every >=2-way danger-set intersection")
    print("proof-facing object: component--danger-set discrepancy ledger")
    print("challenged assumption: proof vertices need not be runners; components carry the loss")
    print("DONE")


if __name__ == "__main__":
    main()
