#!/usr/bin/env python3
"""Exact bounded verifier for the THM-4423 height-499 projection atlas.

This program never uses the physical comb to select a coordinate pair.  It
enumerates the complete raw-carrier dictionary and sums the edgewise-minimum
envelope of each of the three component contact networks. THM-4414 identifies
these sums with the exact degree-zero raw projection capacities. Thus a
minimum at most 6/77 is the scoped network certificate; no converse to
physical entry is asserted.

All load-bearing comparisons are integer comparisons with common denominator
14*a*b*c.  A second relation-lattice enumeration audits every row in a small
box and selected rows outside that box.
"""

from fractions import Fraction as Q
from functools import lru_cache
from hashlib import sha256
from itertools import combinations
from json import dumps
from math import gcd
import argparse
import heapq
import sys
import time


TARGET = Q(6, 77)
PAIR_NAMES = ("01|2", "02|1", "12|0")


class VerificationError(RuntimeError):
    pass


CHECKS = 0


def require(condition, payload):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise VerificationError(payload)


def strict_floor(numerator, denominator):
    """Largest integer strictly smaller than numerator/denominator."""
    return (numerator - 1) // denominator


@lru_cache(maxsize=None)
def danger_components(speed):
    """Closed rational endpoint data for |speed*y-owner|<3/14 on [0,1]."""
    return tuple(
        (max(0, 14 * owner - 3), min(14 * speed, 14 * owner + 3), owner)
        for owner in range(speed + 1)
    )


def cross(w, n):
    return (
        w[1] * n[2] - w[2] * n[1],
        w[2] * n[0] - w[0] * n[2],
        w[0] * n[1] - w[1] * n[0],
    )


def carrier_sweep(w, full_audit=False):
    """Exact simultaneous interval sweep in y=3x coordinates."""
    rows = tuple(danger_components(speed) for speed in w)
    cursor = [0, 0, 0]
    carriers = []
    seen = set()
    while all(cursor[index] < len(rows[index]) for index in range(3)):
        current = tuple(rows[index][cursor[index]] for index in range(3))

        # Endpoints are numerator/(14*w_i); the common factor 14 cancels.
        lower_num, lower_den = current[0][0], w[0]
        upper_num, upper_den = current[0][1], w[0]
        for index in (1, 2):
            if current[index][0] * lower_den > lower_num * w[index]:
                lower_num, lower_den = current[index][0], w[index]
            if current[index][1] * upper_den < upper_num * w[index]:
                upper_num, upper_den = current[index][1], w[index]

        if lower_num * upper_den < upper_num * lower_den:
            owners = tuple(current[index][2] for index in range(3))
            sheet_labels = tuple((-w[index] * owners[index]) % 3 for index in range(3))
            if set(sheet_labels) == {0, 1, 2}:
                carrier = cross(w, owners)
                require(carrier not in seen, ("duplicate-carrier", w, owners, carrier))
                seen.add(carrier)
                carriers.append(carrier)
                if full_audit:
                    require(sum(w[i] * carrier[i] for i in range(3)) == 0,
                            ("relation", w, carrier))
                    require(all(value % 3 for value in carrier),
                            ("live-mod-three", w, carrier))
                    pair_numerators = (
                        3 * (w[0] + w[1]) - 14 * abs(carrier[2]),
                        3 * (w[0] + w[2]) - 14 * abs(carrier[1]),
                        3 * (w[1] + w[2]) - 14 * abs(carrier[0]),
                    )
                    require(all(value > 0 for value in pair_numerators),
                            ("positive-roof", w, carrier, pair_numerators))

        # Advance every interval attaining the earliest right endpoint.
        earliest_num, earliest_den = current[0][1], w[0]
        for index in (1, 2):
            if current[index][1] * earliest_den < earliest_num * w[index]:
                earliest_num, earliest_den = current[index][1], w[index]
        for index in range(3):
            if current[index][1] * earliest_den == earliest_num * w[index]:
                cursor[index] += 1

    return tuple(carriers)


def carrier_lattice_box(w):
    """Independent exact enumeration in ker(w), with no interval sweep."""
    bounds = (
        strict_floor(3 * (w[1] + w[2]), 14),
        strict_floor(3 * (w[0] + w[2]), 14),
        strict_floor(3 * (w[0] + w[1]), 14),
    )
    result = set()
    for first in range(-bounds[0], bounds[0] + 1):
        for second in range(-bounds[1], bounds[1] + 1):
            numerator = -(w[0] * first + w[1] * second)
            if numerator % w[2]:
                continue
            third = numerator // w[2]
            carrier = (first, second, third)
            if abs(third) > bounds[2] or not all(value % 3 for value in carrier):
                continue
            pair_numerators = (
                3 * (w[0] + w[1]) - 14 * abs(third),
                3 * (w[0] + w[2]) - 14 * abs(second),
                3 * (w[1] + w[2]) - 14 * abs(first),
            )
            if all(value > 0 for value in pair_numerators):
                result.add(carrier)
    return result


def edge_envelopes(w, carriers):
    """Return three envelope numerators over the common denominator 14abc."""
    a, b, c = w
    sums = [0, 0, 0]
    physical_numerator = 0
    for first, second, third in carriers:
        pair_numerators = (
            3 * (a + b) - 14 * abs(third),
            3 * (a + c) - 14 * abs(second),
            3 * (b + c) - 14 * abs(first),
        )
        require(all(value > 0 for value in pair_numerators),
                ("envelope-positive", w, (first, second, third), pair_numerators))
        # A pair component is itself capped by the shorter of its two sheet
        # intervals, while the opposite endpoint is capped by the third-sheet
        # interval.  Since a<b<c, the minimum of those three caps is always
        # the global shortest cap 3/(7c), whose numerator over 14abc is 6ab.
        global_cap_numerator = 6 * a * b
        edge_numerators = (
            min(pair_numerators[0] * c, global_cap_numerator),
            min(pair_numerators[1] * b, global_cap_numerator),
            min(pair_numerators[2] * a, global_cap_numerator),
        )
        for index in range(3):
            sums[index] += edge_numerators[index]
        physical_numerator += min(edge_numerators)
    return tuple(sums), physical_numerator, 14 * a * b * c


def eligible_values(height):
    return tuple(value for value in range(1, height + 1, 2) if value % 3)


def selected_for_independent_audit(w, height, audit_height):
    if w[2] <= audit_height:
        return True
    # Deterministic outer controls: the strongest observed rays, a terminal
    # shell control, and a very thin arithmetic sample.  The complete inner
    # H=79 universe is already audited above, so avoid turning the independent
    # O(height^2) lattice enumeration into the main scan.
    fixed = {
        (1, 53, 107), (1, 67, 133), (5, 71, 137), (1, 95, 191),
        (1, 221, 443), (1, 235, 469), (5, 251, 497), (5, 247, 499),
    }
    return w in fixed or (13 * w[0] + 17 * w[1] + 19 * w[2]) % 100003 == 0


def scan(height, audit_height, top_count):
    start = time.time()
    values = eligible_values(height)
    triples = 0
    audited = 0
    target_rows = []
    hostile_rows = []
    pair_success = [0, 0, 0]
    chosen_pairs = [0, 0, 0]
    carrier_min = None
    carrier_max = 0
    physical_above = 0
    best_heap = []
    outer_best = None
    shell_best = {}
    average_counterexample = None
    semantic_rows = []

    for w in combinations(values, 3):
        if gcd(gcd(w[0], w[1]), w[2]) != 1:
            continue
        triples += 1
        audit = selected_for_independent_audit(w, height, audit_height)
        carriers = carrier_sweep(w, full_audit=audit)
        if audit:
            independently = carrier_lattice_box(w)
            require(set(carriers) == independently,
                    ("sweep-vs-lattice", w, set(carriers) ^ independently))
            audited += 1
        envelope_numerators, physical_numerator, denominator = edge_envelopes(w, carriers)
        best_numerator = min(envelope_numerators)
        chosen = envelope_numerators.index(best_numerator)
        chosen_pairs[chosen] += 1
        for index, numerator in enumerate(envelope_numerators):
            if numerator * 77 <= 6 * denominator:
                pair_success[index] += 1

        require(physical_numerator <= best_numerator,
                ("physical-below-envelope", w, physical_numerator, best_numerator))
        if physical_numerator * 77 > 6 * denominator:
            physical_above += 1
        if best_numerator * 77 > 6 * denominator:
            hostile_rows.append((w, best_numerator, denominator, envelope_numerators))
        elif best_numerator * 77 == 6 * denominator:
            target_rows.append((w, envelope_numerators, denominator))

        carrier_count = len(carriers)
        carrier_min = carrier_count if carrier_min is None else min(carrier_min, carrier_count)
        carrier_max = max(carrier_max, carrier_count)

        rank_row = (Q(best_numerator, denominator), w, envelope_numerators,
                    denominator, carrier_count)
        if len(best_heap) < top_count:
            heapq.heappush(best_heap, rank_row)
        elif rank_row > best_heap[0]:
            heapq.heapreplace(best_heap, rank_row)
        if w[2] > 79 and (outer_best is None or rank_row > outer_best):
            outer_best = rank_row
        previous = shell_best.get(w[2])
        if previous is None or rank_row > previous:
            shell_best[w[2]] = rank_row

        if (average_counterexample is None
                and sum(envelope_numerators) * 77 > 18 * denominator):
            average_counterexample = (w, envelope_numerators, denominator)

        semantic_rows.append((w, envelope_numerators, denominator, carrier_count))

    require(triples > 0, "nonempty-universe")
    require(not hostile_rows, ("bounded-hostiles", hostile_rows[:10]))
    require(physical_above == 0, ("physical-controls-above", physical_above))
    require(target_rows == [((1, 5, 11), (60, 60, 60), 770)],
            ("unique-envelope-equality", target_rows))
    require(average_counterexample is not None, "average-counterexample-present")

    best_rows = sorted(best_heap, reverse=True)
    shell_tail = tuple(
        (shell, row[0], row[1], row[4])
        for shell, row in sorted(shell_best.items())[-20:]
    )
    semantic_digest = sha256(dumps(semantic_rows, sort_keys=True).encode()).hexdigest()
    result = {
        "height": height,
        "audit_height": audit_height,
        "eligible_values": values,
        "triple_count": triples,
        "independently_audited_rows": audited,
        "target_rows": target_rows,
        "hostile_count": len(hostile_rows),
        "physical_above_count": physical_above,
        "pair_success": tuple(pair_success),
        "chosen_pairs": tuple(chosen_pairs),
        "carrier_count_range": (carrier_min, carrier_max),
        "best_rows": best_rows,
        "outer_best": outer_best,
        "shell_tail": shell_tail,
        "average_counterexample": average_counterexample,
        "rows_sha256": semantic_digest,
        "checks": CHECKS,
        "seconds": time.time() - start,
    }
    return result


def fraction_tuple(numerators, denominator):
    return tuple(str(Q(value, denominator)) for value in numerators)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--height", type=int, default=499)
    parser.add_argument("--audit-height", type=int, default=79)
    parser.add_argument("--top", type=int, default=20)
    args = parser.parse_args()
    require(args.height >= 79 and args.height % 2 == 1, "height-odd-at-least-79")
    require(1 <= args.audit_height <= args.height, "audit-height-range")
    sys.stdout.reconfigure(newline="\n")
    result = scan(args.height, args.audit_height, args.top)
    print("LRC14 HEIGHT-499 PROJECTION ATLAS THM-4423 VERIFIER")
    print("status=PASS FINITE_EXACT; low_carrier_splice=PROVED; universal_statement=OPEN; LRC14=OPEN")
    print("universe=sorted primitive distinct positive odd ternary-unit triples; max_speed<=%d"
          % result["height"])
    print("eligible_speed_count=%d; triples=%d; independently_audited_rows=%d"
          % (len(result["eligible_values"]), result["triple_count"],
             result["independently_audited_rows"]))
    print("certificate=network_capacity<=edgewise_min_envelope; selection_uses_no_physical_mass")
    print("hostile_count=%d; target_rows=%s; physical_above_count=%d"
          % (result["hostile_count"], result["target_rows"], result["physical_above_count"]))
    print("pair_success_counts=%s; selected_pair_counts=%s"
          % (result["pair_success"], result["chosen_pairs"]))
    print("carrier_count_range=%s" % (result["carrier_count_range"],))
    outer = result["outer_best"]
    if outer is None:
        print("strongest_beyond_79=none_in_universe")
    else:
        print("strongest_beyond_79=triple:%s,best:%s,envelopes:%s,carriers:%d"
              % (outer[1], outer[0], fraction_tuple(outer[2], outer[3]), outer[4]))
    print("top_rows=")
    for row in result["best_rows"]:
        print("  triple=%s best=%s envelopes=%s carriers=%d"
              % (row[1], row[0], fraction_tuple(row[2], row[3]), row[4]))
    avg_w, avg_nums, avg_den = result["average_counterexample"]
    print("unweighted_average_claim_refuted=triple:%s,sum:%s,target_sum:%s,envelopes:%s"
          % (avg_w, Q(sum(avg_nums), avg_den), Q(18, 77),
             fraction_tuple(avg_nums, avg_den)))
    print("last_twenty_shell_maxima=")
    for shell, value, w, count in result["shell_tail"]:
        print("  shell=%d triple=%s best=%s carriers=%d" % (shell, w, value, count))
    print("rows_sha256=" + result["rows_sha256"])
    print("checks=%d" % result["checks"])
    print("verdict=PASS")


if __name__ == "__main__":
    main()
