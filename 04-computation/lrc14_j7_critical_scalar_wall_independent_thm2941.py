#!/usr/bin/env python3
"""Independent exact audit of the critical seven-slot Hunter scalar.

This file deliberately does not import the discovery probe in
``.scratch/lrc_j7_critical_hunter_degeneracy_20260729`` or any of its
aggregates.  It reconstructs every six-subset of ``{1,...,14}`` directly
from the defining danger combs

    D_w = {t in R/Z : ||w t|| < 1/14}

and the carrier

    G_E = (R/Z) \\ union_(e in E) D_e.

All carrier endpoints are represented on the common integer ruler
``14*lcm(1,...,14)``.  Singleton coverages use a separate periodic
antiderivative.  Restricted pair intersections use an exact reduced-period
integer sweep and are checked against the corrected LEM-042 capped-trapezoid
formula.  Thus the implementation is independent of the interval engine and
tail enumeration used by the discovery probe.

The only analytic input is the proved strict THM-735 tail

    c_E(w) < h_E/7 + gamma_E/w,  gamma_E = 99 r_E / 490,

which dynamically seals both the global top seven and the global pair-union
maximum without an arbitrary label cap.
"""

from __future__ import annotations

import argparse
import hashlib
import math
import multiprocessing as mp
from bisect import bisect_right
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations


BASE_LABEL = 15
INITIAL_HEAD = 256
RULER = 14 * math.lcm(*range(1, 15))
ONE_FOURTEENTH = RULER // 14
THIRTEEN_FOURTEENTHS = 13 * ONE_FOURTEENTH
ONE_SEVENTH = RULER // 7


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def ceil_fraction(value: F) -> int:
    return -(-value.numerator // value.denominator)


def ftext(value: F) -> str:
    return f"{value.numerator}/{value.denominator}"


def merge_intervals(
    intervals: list[tuple[int, int]],
) -> tuple[tuple[int, int], ...]:
    """Merge positive-length closed interval carriers; endpoints have zero mass."""
    rows = sorted((left, right) for left, right in intervals if left < right)
    if not rows:
        return ()
    merged: list[list[int]] = [[rows[0][0], rows[0][1]]]
    for left, right in rows[1:]:
        if left <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], right)
        else:
            merged.append([left, right])
    return tuple((left, right) for left, right in merged)


def danger_intervals(speed: int, ruler: int) -> tuple[tuple[int, int], ...]:
    """D_speed on [0,1], with all endpoints integral on ``ruler``."""
    require(ruler % (14 * speed) == 0, "danger ruler is not exact")
    step = ruler // speed
    radius = ruler // (14 * speed)
    rows = [(0, radius)]
    rows.extend(
        (index * step - radius, index * step + radius)
        for index in range(1, speed)
    )
    rows.append((ruler - radius, ruler))
    return tuple(rows)


def carrier_for(body: tuple[int, ...]) -> tuple[tuple[int, int], ...]:
    danger = merge_intervals(
        [
            interval
            for speed in body
            for interval in danger_intervals(speed, RULER)
        ]
    )
    carrier: list[tuple[int, int]] = []
    cursor = 0
    for left, right in danger:
        if cursor < left:
            carrier.append((cursor, left))
        cursor = max(cursor, right)
    if cursor < RULER:
        carrier.append((cursor, RULER))
    # Every danger comb contains zero, so a positive carrier never wraps.
    require(
        not carrier
        or (carrier[0][0] > 0 and carrier[-1][1] < RULER),
        "carrier unexpectedly wraps around zero",
    )
    return tuple(carrier)


def danger_primitive(numerator: int) -> int:
    """RULER times integral_0^(numerator/RULER) 1_{D_1}(s) ds."""
    cycles, remainder = divmod(numerator, RULER)
    return (
        cycles * ONE_SEVENTH
        + min(remainder, ONE_FOURTEENTH)
        + max(0, remainder - THIRTEEN_FOURTEENTHS)
    )


def singleton_coverage(
    carrier: tuple[tuple[int, int], ...],
    label: int,
) -> F:
    numerator = sum(
        danger_primitive(label * right) - danger_primitive(label * left)
        for left, right in carrier
    )
    return F(numerator, RULER * label)


def intersect_intervals(
    first: tuple[tuple[int, int], ...],
    second: tuple[tuple[int, int], ...],
) -> tuple[tuple[int, int], ...]:
    rows: list[tuple[int, int]] = []
    i = 0
    j = 0
    while i < len(first) and j < len(second):
        left = max(first[i][0], second[j][0])
        right = min(first[i][1], second[j][1])
        if left < right:
            rows.append((left, right))
        if first[i][1] <= second[j][1]:
            i += 1
        else:
            j += 1
    return merge_intervals(rows)


def lem042_pair_mass(first: int, second: int) -> F:
    """Corrected full-period capped-trapezoid law at radius 1/14."""
    a, b = sorted((first, second))
    g = math.gcd(a, b)
    width = F(a + b, 14 * a * b)
    delta = F(g, a * b)
    cap = F(1, 7 * b)
    limit = (a + b) // (14 * g) + 1
    total = F(0)
    for index in range(-limit, limit + 1):
        height = width - abs(index) * delta
        if height > 0:
            total += min(height, cap)
    return g * total


@lru_cache(maxsize=512)
def pair_pattern(
    first: int,
    second: int,
) -> tuple[
    int,
    int,
    tuple[int, ...],
    tuple[int, ...],
    tuple[int, ...],
    int,
]:
    """Reduced-period D_first cap D_second pattern.

    Returns ``(g,L,starts,ends,prefix,total)``.  ``prefix[k]`` is the total
    length of intervals with index below ``k``.
    """
    require(first != second, "pair labels must be distinct")
    g = math.gcd(first, second)
    a = first // g
    b = second // g
    ruler = math.lcm(RULER, 14 * a, 14 * b)
    overlap = intersect_intervals(
        danger_intervals(a, ruler),
        danger_intervals(b, ruler),
    )
    starts = tuple(left for left, _ in overlap)
    ends = tuple(right for _, right in overlap)
    prefix_rows = [0]
    for left, right in overlap:
        prefix_rows.append(prefix_rows[-1] + right - left)
    prefix = tuple(prefix_rows)
    total = prefix[-1]
    require(
        F(total, ruler) == lem042_pair_mass(first, second),
        "pair pattern disagrees with corrected LEM-042",
    )
    return g, ruler, starts, ends, prefix, total


def pattern_prefix(
    position: int,
    starts: tuple[int, ...],
    ends: tuple[int, ...],
    prefix: tuple[int, ...],
) -> int:
    complete = bisect_right(ends, position)
    total = prefix[complete]
    if complete < len(starts) and starts[complete] < position:
        total += min(position, ends[complete]) - starts[complete]
    return total


def restricted_pair_intersection(
    carrier: tuple[tuple[int, int], ...],
    first: int,
    second: int,
) -> F:
    g, ruler, starts, ends, prefix, period_mass = pair_pattern(first, second)
    ruler_ratio = ruler // RULER

    def cumulative(numerator: int) -> int:
        cycles, remainder = divmod(numerator, RULER)
        position = remainder * ruler_ratio
        return cycles * period_mass + pattern_prefix(
            position,
            starts,
            ends,
            prefix,
        )

    numerator = sum(
        cumulative(g * right) - cumulative(g * left)
        for left, right in carrier
    )
    return F(numerator, ruler * g)


def full_grid_pair_intersection(
    carrier: tuple[tuple[int, int], ...],
    first: int,
    second: int,
) -> F:
    """Independent literal common-ruler control for the reduced-period engine."""
    ruler = math.lcm(RULER, 14 * first, 14 * second)
    scale = ruler // RULER
    carrier_scaled = tuple(
        (left * scale, right * scale) for left, right in carrier
    )
    overlap = intersect_intervals(
        intersect_intervals(
            carrier_scaled,
            danger_intervals(first, ruler),
        ),
        danger_intervals(second, ruler),
    )
    return F(sum(right - left for left, right in overlap), ruler)


def pair_union(
    carrier: tuple[tuple[int, int], ...],
    coverages: dict[int, F],
    first: int,
    second: int,
) -> F:
    overlap = restricted_pair_intersection(carrier, first, second)
    require(
        0 <= overlap <= min(coverages[first], coverages[second]),
        "restricted pair intersection left its singleton bounds",
    )
    return coverages[first] + coverages[second] - overlap


def finite_pair_cap(
    carrier: tuple[tuple[int, int], ...],
    coverages: dict[int, F],
) -> tuple[F, tuple[int, int], int]:
    labels = tuple(
        sorted(coverages, key=lambda label: (-coverages[label], label))
    )
    require(len(labels) >= 2, "pair head has fewer than two labels")
    warm = tuple(sorted(labels[:2]))
    best = pair_union(carrier, coverages, *warm)
    witness = warm
    paid = 1
    for first_index in range(len(labels) - 1):
        first = labels[first_index]
        if coverages[first] + coverages[labels[first_index + 1]] < best:
            break
        for second_index in range(first_index + 1, len(labels)):
            second = labels[second_index]
            if coverages[first] + coverages[second] < best:
                break
            pair = tuple(sorted((first, second)))
            if pair == warm:
                continue
            value = pair_union(carrier, coverages, *pair)
            paid += 1
            if value > best or (value == best and pair < witness):
                best = value
                witness = pair
    return best, witness, paid


def extend_coverages(
    carrier: tuple[tuple[int, int], ...],
    coverages: dict[int, F],
    target: int,
) -> None:
    current = max(coverages, default=BASE_LABEL - 1)
    for label in range(current + 1, target + 1):
        coverages[label] = singleton_coverage(carrier, label)


def profile_root(body: tuple[int, ...]) -> dict[str, object]:
    carrier = carrier_for(body)
    require(carrier, "six-body carrier is empty")
    mass = F(sum(right - left for left, right in carrier), RULER)
    components = len(carrier)
    gamma = F(99 * components, 490)
    for label in body:
        require(
            singleton_coverage(carrier, label) == 0,
            "body danger comb meets its own good carrier",
        )

    coverages: dict[int, F] = {}
    head = INITIAL_HEAD
    extend_coverages(carrier, coverages, head)
    while True:
        ranked = tuple(
            sorted(coverages, key=lambda label: (-coverages[label], label))
        )
        top7 = tuple((label, coverages[label]) for label in ranked[:7])
        q1 = top7[0][1]
        q7 = top7[-1][1]
        q7_gap = q7 - mass / 7
        require(q7_gap > 0, "q7 did not clear the critical density")
        required = max(INITIAL_HEAD, ceil_fraction(gamma / q7_gap) - 1)
        if required <= head:
            require(
                mass / 7 + gamma / (head + 1) <= q7,
                "top-seven tail did not seal",
            )
            break
        head = required
        extend_coverages(carrier, coverages, head)
    q_head = head

    while True:
        pair_cap, pair_witness, paid = finite_pair_cap(carrier, coverages)
        ranked = tuple(
            sorted(coverages, key=lambda label: (-coverages[label], label))
        )
        q1 = coverages[ranked[0]]
        pair_tail_gap = pair_cap - q1 - mass / 7
        require(pair_tail_gap > 0, "global pair cap lacks a finite tail gap")
        required = max(q_head, ceil_fraction(gamma / pair_tail_gap) - 1)
        if required <= head:
            tau = mass / 7 + gamma / (head + 1)
            require(
                pair_cap <= 2 * q1,
                "finite pair cap exceeded twice the leading singleton",
            )
            require(
                q1 + tau <= pair_cap,
                "one-head/one-tail pair did not seal",
            )
            # The previous two inequalities imply tau<=q1.  This is the
            # missing two-omitted-label case: each omitted singleton is
            # bounded by tau, so their union is at most 2*tau, which is in
            # turn at most q1+tau and hence at most the finite pair winner.
            require(
                tau <= q1
                and 2 * tau <= q1 + tau
                and q1 + tau <= pair_cap,
                "two-tail pair did not seal",
            )
            break
        head = required
        extend_coverages(carrier, coverages, head)

    # Re-rank after any pair-driven extension.
    ranked = tuple(
        sorted(coverages, key=lambda label: (-coverages[label], label))
    )
    top7 = tuple((label, coverages[label]) for label in ranked[:7])
    q1 = top7[0][1]
    q7 = top7[-1][1]
    q7_gap = q7 - mass / 7
    require(
        mass / 7 + gamma / (head + 1) <= q7,
        "final scan lost the top-seven seal",
    )
    pair_gap = pair_cap - 2 * mass / 7
    require(pair_gap > 0, "B2 did not clear twice the critical density")

    critical = mass / 7

    def hunter(center: F) -> F:
        return center + sum(
            min(center, value, pair_cap - center)
            for _, value in top7[1:]
        )

    upper = min(q1, pair_cap)

    def clipped(value: F) -> F:
        return max(F(0), min(upper, value))

    candidates = {F(0), upper, clipped(pair_cap / 2)}
    for _, value in top7[1:]:
        candidates.add(clipped(value))
        candidates.add(clipped(pair_cap - value))
    ordered_candidates = tuple(sorted(candidates))
    for left, right in zip(ordered_candidates, ordered_candidates[1:]):
        midpoint = (left + right) / 2
        require(
            2 * hunter(midpoint) == hunter(left) + hunter(right),
            "Hunter breakpoint list is incomplete",
        )
    star_center = min(
        ordered_candidates,
        key=lambda center: (-hunter(center), center),
    )
    star_envelope = hunter(star_center)
    envelope_gap = star_envelope - mass

    require(
        hunter(critical) == mass
        and all(value > critical for _, value in top7)
        and pair_cap - critical > critical,
        "critical Hunter identity failed",
    )
    # For every a<critical, each of the seven summands is at most a, hence
    # G7(a)<=7a<h.  The exact first hostile level is therefore critical.
    threshold = critical

    return {
        "body": body,
        "mass": mass,
        "components": components,
        "gamma": gamma,
        "top7": top7,
        "q7": q7,
        "q7_gap": q7_gap,
        "q_head": q_head,
        "pair_cap": pair_cap,
        "pair_witness": pair_witness,
        "pair_gap": pair_gap,
        "pair_tail_gap": pair_tail_gap,
        "pair_head": head,
        "paid_pairs": paid,
        "threshold": threshold,
        "hunter_at_threshold": hunter(threshold),
        "star_center": star_center,
        "star_envelope": star_envelope,
        "envelope_gap": envelope_gap,
    }


def row_digest(rows: list[dict[str, object]]) -> str:
    payload = tuple(
        (
            row["body"],
            row["mass"],
            row["components"],
            row["top7"],
            row["q_head"],
            row["pair_cap"],
            row["pair_witness"],
            row["pair_head"],
            row["threshold"],
            row["star_center"],
            row["star_envelope"],
        )
        for row in rows
    )
    return hashlib.sha256(repr(payload).encode()).hexdigest()


def extremum(
    rows: list[dict[str, object]],
    field: str,
    largest: bool = False,
) -> dict[str, object]:
    if largest:
        return min(
            rows,
            key=lambda row: (-row[field], row["body"]),
        )
    return min(rows, key=lambda row: (row[field], row["body"]))


def print_extremum(
    name: str,
    row: dict[str, object],
    field: str,
) -> None:
    value = row[field]
    rendered = ftext(value) if isinstance(value, F) else str(value)
    print(
        f"{name}=(value={rendered},body={row['body']},"
        f"h={ftext(row['mass'])},q7={ftext(row['q7'])},"
        f"top7={tuple(label for label, _ in row['top7'])},"
        f"B2={ftext(row['pair_cap'])},pair={row['pair_witness']},"
        f"q_head={row['q_head']},pair_head={row['pair_head']},"
        f"paid_pairs={row['paid_pairs']})"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--workers",
        type=int,
        default=min(8, mp.cpu_count() or 1),
    )
    args = parser.parse_args()
    require(args.workers >= 1, "worker count must be positive")
    require(RULER == 5_045_040, "base ruler changed")

    roots = tuple(combinations(range(1, 15), 6))
    require(len(roots) == math.comb(14, 6) == 3_003, "root universe changed")
    if args.workers == 1:
        rows = [profile_root(body) for body in roots]
    else:
        with mp.get_context("spawn").Pool(args.workers) as pool:
            rows = list(pool.imap(profile_root, roots, chunksize=1))
    rows.sort(key=lambda row: row["body"])

    # Three literal common-grid controls do not share the reduced-period path.
    control_indices = (0, len(rows) // 2, len(rows) - 1)
    control_pairs = ((15, 16), (23, 39), (50, 51))
    for index, pair in zip(control_indices, control_pairs):
        carrier = carrier_for(rows[index]["body"])
        require(
            restricted_pair_intersection(carrier, *pair)
            == full_grid_pair_intersection(carrier, *pair),
            "literal pair-grid control failed",
        )

    require(
        len(rows) == 3_003
        and all(row["q7_gap"] > 0 for row in rows)
        and all(row["pair_gap"] > 0 for row in rows)
        and all(row["threshold"] == row["mass"] / 7 for row in rows)
        and all(row["hunter_at_threshold"] == row["mass"] for row in rows),
        "all-root critical audit failed",
    )

    print("LRC14 independent critical seven-slot Hunter audit")
    print(
        "universe=(six_subsets=3003,base_labels=1..14,"
        "external_labels=integers>=15)"
    )
    print(
        "counts=(q7_strict=3003,B2_strict=3003,"
        "lambda_equal_h_over_7=3003,direct_scalar_closures=0)"
    )
    for field in (
        "q7_gap",
        "pair_gap",
        "pair_tail_gap",
        "envelope_gap",
        "mass",
        "components",
        "q_head",
        "pair_head",
        "paid_pairs",
    ):
        print_extremum(f"minimum_{field}", extremum(rows, field), field)
        print_extremum(
            f"maximum_{field}",
            extremum(rows, field, largest=True),
            field,
        )
    print(f"row_digest={row_digest(rows)}")
    print(
        "mechanism=if q7>=h/7 and B2>=2h/7 then "
        "G7(h/7)=h; for a<h/7, G7(a)<=7a<h"
    )
    print(
        "tail=THM735 strict singleton discrepancy dynamically seals "
        "top seven and pair cap; no arbitrary label ceiling"
    )
    print(
        "two_tail_seal=H2<=2q1 and H2>=q1+tau imply tau<=q1; "
        "therefore two omitted labels have union <=2tau<=q1+tau<=H2"
    )
    print(
        "pair_controls=corrected LEM042 trapezoid plus three independent "
        "literal common-grid restricted-overlap checks"
    )
    print("scope=FINITE-EXACT audit of scalar G7 only; not LRC14")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    main()
