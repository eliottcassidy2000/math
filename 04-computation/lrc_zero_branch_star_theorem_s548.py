#!/usr/bin/env python3
"""
lrc_zero_branch_star_theorem_s548.py

Exact verifier for THM-391, the q-grid zero-branch star peeling theorem.

For n >= 2 and 2 <= q <= n, choose centers C among nonzero q-grid points u/q
and speeds S all divisible by q.  The local danger intervals

    (u/q - 1/(n s), u/q + 1/(n s)),    s in S

are separated between different q-grid centers and nested at each center.  The
endpoint-protection core is therefore empty.  More sharply, the peeling layers
are explicit: speeds are removed in increasing speed order, with
|C| * multiplicity(s) intervals removed at the layer for speed s.

The theorem does not use q being a prime power.  Prime-power zero branches are
the LRC application because they are the p-adic tree branches in HYP-2036.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations_with_replacement
from math import gcd


@dataclass(frozen=True)
class StarCase:
    n: int
    q: int
    centers: tuple[int, ...]
    speeds: tuple[int, ...]


def interval_for(n: int, q: int, u: int, speed: int) -> tuple[Fraction, Fraction, str]:
    center = Fraction(u, q)
    radius = Fraction(1, n * speed)
    return (center - radius, center + radius, f"u={u},s={speed}")


def star_intervals(case: StarCase) -> list[tuple[Fraction, Fraction, str]]:
    return [
        interval_for(case.n, case.q, u, speed)
        for u in case.centers
        for speed in case.speeds
    ]


def inside_open(point: Fraction, interval: tuple[Fraction, Fraction, str]) -> bool:
    left, right, _ = interval
    return left < point < right


def endpoint_core(intervals: list[tuple[Fraction, Fraction, str]]) -> tuple[int, tuple[int, ...]]:
    active = set(range(len(intervals)))
    layers: list[int] = []
    while True:
        remove: set[int] = set()
        for idx in active:
            left, right, _ = intervals[idx]
            left_protected = any(
                j != idx and j in active and inside_open(left, intervals[j])
                for j in active
            )
            right_protected = any(
                j != idx and j in active and inside_open(right, intervals[j])
                for j in active
            )
            if not (left_protected and right_protected):
                remove.add(idx)
        if not remove:
            break
        layers.append(len(remove))
        active -= remove
    return len(active), tuple(layers)


def predicted_layers(case: StarCase) -> tuple[int, ...]:
    counts = Counter(case.speeds)
    center_count = len(case.centers)
    return tuple(center_count * counts[speed] for speed in sorted(counts))


def center_families(q: int) -> dict[str, tuple[int, ...]]:
    units = tuple(u for u in range(1, q) if gcd(u, q) == 1)
    nonunits = tuple(u for u in range(1, q) if gcd(u, q) != 1)
    families = {
        "units": units,
        "all_nonzero": tuple(range(1, q)),
        "first": (1,),
    }
    if nonunits:
        families["nonunits"] = nonunits
    return families


def verify_case(case: StarCase) -> tuple[bool, str]:
    intervals = star_intervals(case)
    core_size, layers = endpoint_core(intervals)
    expected = predicted_layers(case)
    if core_size != 0:
        return False, f"nonempty core={core_size}: {case}"
    if layers != expected:
        return False, f"layers {layers} != expected {expected}: {case}"
    return True, "ok"


def bounded_verification() -> tuple[int, int, int]:
    total = 0
    max_intervals = 0
    max_layers = 0
    for n in range(2, 9):
        for q in range(2, n + 1):
            for centers in center_families(q).values():
                for length in range(0, 4):
                    for multipliers in combinations_with_replacement(range(1, 5), length):
                        speeds = tuple(q * m for m in multipliers)
                        case = StarCase(n=n, q=q, centers=centers, speeds=speeds)
                        ok, message = verify_case(case)
                        if not ok:
                            raise AssertionError(message)
                        total += 1
                        max_intervals = max(max_intervals, len(star_intervals(case)))
                        max_layers = max(max_layers, len(predicted_layers(case)))
    return total, max_intervals, max_layers


def selected_examples() -> list[StarCase]:
    return [
        StarCase(n=18, q=9, centers=tuple(u for u in range(1, 9) if gcd(u, 9) == 1), speeds=(9,)),
        StarCase(n=18, q=9, centers=tuple(u for u in range(1, 9) if gcd(u, 9) == 1), speeds=(9, 18, 36)),
        StarCase(n=18, q=18, centers=tuple(u for u in range(1, 18) if gcd(u, 18) == 1), speeds=(18, 36)),
        StarCase(n=14, q=7, centers=tuple(u for u in range(1, 7)), speeds=(7, 14, 21)),
        # Non-prime-power example: the proof still works because only the q-grid
        # spacing and the radius bound s >= q are used.
        StarCase(n=12, q=6, centers=tuple(range(1, 6)), speeds=(6, 12, 18)),
    ]


def main() -> None:
    print("THM-391 q-grid zero-branch star peeling verifier")
    print("=" * 72)
    total, max_intervals, max_layers = bounded_verification()
    print(f"bounded exact cases checked: {total}")
    print(f"max intervals in bounded check: {max_intervals}")
    print(f"max predicted peel layers: {max_layers}")
    print("all bounded cores empty and all peel layers match the formula")
    print()

    print("selected LRC-style stars")
    for case in selected_examples():
        intervals = star_intervals(case)
        core_size, layers = endpoint_core(intervals)
        expected = predicted_layers(case)
        center_label = "units" if all(gcd(u, case.q) == 1 for u in case.centers) else "nonzero"
        print(
            f"  n={case.n:2d} q={case.q:2d} centers={center_label:7s} "
            f"|C|={len(case.centers):2d} speeds={case.speeds}: "
            f"intervals={len(intervals):2d} core={core_size} "
            f"layers={layers} expected={expected}"
        )

    print()
    print("synthesis")
    print("  The local star theorem is q-agnostic: primality is not used.")
    print("  Prime powers matter because they are p-adic tree branches, not because")
    print("  their local interval geometry is harder.")
    print("  A nonempty LRC endpoint core must mix descendant/off-grid/event labels;")
    print("  one zero branch cannot hold the core by itself.")


if __name__ == "__main__":
    main()
