#!/usr/bin/env python3
"""Finite-exact low-core closure for THM-769/774's two-sheet branch.

For every ten-set U in [1,N], the widest strict component of G_U bounds any
odd speed whose closed doubled-danger teeth could cover all of G_U.  The
remaining finite odd speeds are tested exactly.  At N=19 there is no universal
odd runner, hence no two-runner opposite-colour cover, even though the odd
speeds were unbounded before the intrinsic component-width cap.
"""

import argparse
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
from math import comb


DELTA = Fraction(1, 13)
FOLDED_RADIUS = Fraction(2, 13)
EXPECTED_DIGEST_N19 = (
    "ec206bf06eda11b5f8ee5318b2bdbc97d61ae63c78e508c83610ce3a8a2dcf83"
)


def intersect(left, right):
    answer = []
    i = j = 0
    while i < len(left) and j < len(right):
        a, b = left[i]
        c, d = right[j]
        lo, hi = max(a, c), min(b, d)
        if hi > lo:
            answer.append((lo, hi))
        if b <= d:
            i += 1
        else:
            j += 1
    return answer


def strict_safe_components(speeds):
    """Components of {tau:min_v ||v*tau||>1/13}, with exact endpoints."""
    current = [(Fraction(0), Fraction(1))]
    for speed in speeds:
        safe = [
            (
                (Fraction(k) + DELTA) / speed,
                (Fraction(k + 1) - DELTA) / speed,
            )
            for k in range(speed)
        ]
        current = intersect(current, safe)
        if not current:
            break
    return tuple(current)


def nearest_integer(z):
    q, r = divmod(z.numerator, z.denominator)
    if 2 * r < z.denominator:
        return q
    if 2 * r > z.denominator:
        return q + 1
    return None


def universal_folded_signature(components, speed):
    """Parity/slack word iff one odd speed covers every open component."""
    word = []
    slacks = []
    for left, right in components:
        centre = (left + right) / 2
        half_width = (right - left) / 2
        owner = nearest_integer(speed * centre)
        if owner is None:
            return None
        # Closed eligibility accepts equality at each open component's limit.
        slack = (
            FOLDED_RADIUS
            - abs(speed * centre - owner)
            - speed * half_width
        )
        if slack < 0:
            return None
        word.append(owner & 1)
        slacks.append(slack)
    return tuple(word), tuple(slacks)


def total_order_tournament(vertices, key):
    """Scalar comparison gauge with the core tuple as fixed tie path."""
    order = sorted(vertices, key=lambda row: (key(row), row[0]))
    rank = {row[0]: i for i, row in enumerate(order)}
    edges = set()
    scores = {row[0]: 0 for row in vertices}
    for i, left in enumerate(vertices):
        for right in vertices[i + 1 :]:
            u, v = left[0], right[0]
            edge = (u, v) if rank[u] < rank[v] else (v, u)
            edges.add(edge)
            scores[edge[0]] += 1
    return {
        "edges": edges,
        "score_histogram": tuple(sorted(scores.values())),
        "directed_cycles": 0,
        "scc_sizes": (1,) * len(vertices),
        "hamiltonian_paths": 1,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--N", type=int, default=19)
    args = parser.parse_args()
    if args.N < 10:
        raise SystemExit("N must be at least 10")

    core_count = 0
    divisor_complete = 0
    odd_tests = 0
    pair_capacity = 0
    universal = 0
    summaries = []
    digest = sha256()

    for core in combinations(range(1, args.N + 1), 10):
        core_count += 1
        divisor_complete += all(
            any(speed % modulus == 0 for speed in core)
            for modulus in range(2, 13)
        )
        components = strict_safe_components(core)
        assert components
        widest = max(right - left for left, right in components)
        ratio = Fraction(4, 13) / widest
        cap = ratio.numerator // ratio.denominator
        odd_count = (cap + 1) // 2
        odd_tests += odd_count
        pair_capacity += comb(odd_count, 2)
        summaries.append((core, cap, widest, len(components)))

        flags = []
        for speed in range(1, cap + 1, 2):
            good = universal_folded_signature(components, speed) is not None
            flags.append("1" if good else "0")
            universal += int(good)

        record = (
            str(core)
            + "|"
            + ",".join(
                f"{left.numerator}/{left.denominator}:"
                f"{right.numerator}/{right.denominator}"
                for left, right in components
            )
            + f"|{cap}|"
            + "".join(flags)
            + "\n"
        )
        digest.update(record.encode())

    hardest = sorted(
        summaries, key=lambda row: (-row[1], row[2], row[0])
    )[:13]
    cap_tournament = total_order_tournament(
        hardest, lambda row: (row[1], Fraction(1, 1) / row[2])
    )
    topology_tournament = total_order_tournament(
        hardest, lambda row: (row[3], row[1])
    )
    flips = len(cap_tournament["edges"] ^ topology_tournament["edges"]) // 2

    max_cap = max(row[1] for row in summaries)
    max_record = min(
        (row for row in summaries if row[1] == max_cap),
        key=lambda row: (row[2], row[0]),
    )
    core, cap, widest, component_count = max_record
    hexdigest = digest.hexdigest()

    print("THM-774 finite-exact unbounded-odd low-core closure")
    print("strict G_U / closed odd eligibility endpoint convention")
    print()
    print(
        f"N={args.N} all_10_subsets={core_count} "
        f"expected={comb(args.N,10)} divisor_complete={divisor_complete}"
    )
    print(
        f"odd_(U,w)_tests={odd_tests} implied_pair_capacity={pair_capacity} "
        f"universal_incidences={universal}"
    )
    print(
        f"max_intrinsic_w_cap={cap} at U={core} "
        f"ell_max={widest} components={component_count}"
    )
    print(f"sha256={hexdigest}")
    print()
    print("Tournament Analysis (vertices = 13 hardest core obligations)")
    print("  pairwise observable: relative apparent obstruction difficulty")
    print("  gauge A: intrinsic odd cap, then reciprocal widest component")
    print("  gauge B: component count, then intrinsic odd cap")
    print("  tie Hamiltonian path: lexicographic core tuple")
    print("  score histogram:", cap_tournament["score_histogram"])
    print(
        "  directed cycles:",
        cap_tournament["directed_cycles"],
        topology_tournament["directed_cycles"],
    )
    print("  SCC sizes:", cap_tournament["scc_sizes"])
    print(
        "  Hamiltonian-path counts:",
        cap_tournament["hamiltonian_paths"],
        topology_tournament["hamiltonian_paths"],
    )
    print("  edge flips between gauges:", flips)
    print("  challenged vertex set: runners; selected vertices are core obligations")
    print("  preserved predicate: difficulty ranking only")
    print("  destroyed data: exact endpoints, parity words, and simultaneous ownership")
    print()

    assert core_count == comb(args.N, 10)
    if args.N == 19:
        assert divisor_complete == 3400
        assert odd_tests == 767700
        assert pair_capacity == 3179784
        assert universal == 0
        assert cap == 60
        assert core == (1, 2, 3, 5, 7, 8, 11, 13, 18, 19)
        assert widest == Fraction(6, 1183)
        assert component_count == 18
        assert hexdigest == EXPECTED_DIGEST_N19
        print("FINAL: PASS — no tight s=2 packet with max(U)<=19")
    else:
        print("FINAL: general-N audit complete")


if __name__ == "__main__":
    main()
