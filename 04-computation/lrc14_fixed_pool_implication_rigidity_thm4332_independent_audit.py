#!/usr/bin/env python3
"""Independent interval audit of the fixed pool's implication closure.

Unlike the discovery program, this constructs G_P once from only pool walls,
then intersects its safe open cells and safe wall points directly with each
test label's strict danger intervals.  No hitting-set or combined arrangement
is used.
"""

from fractions import Fraction
from math import floor
import sys


P = (8, 10, 15, 16, 20, 30, 40, 42, 60, 63, 80, 84, 85, 88, 95,
     120, 126, 132, 143, 145, 168, 170, 176, 190, 193, 240, 252,
     264, 286, 290)


def distance_numerator(v: int, x: Fraction) -> tuple[int, int]:
    residue = (v * x.numerator) % x.denominator
    return min(residue, x.denominator - residue), x.denominator


def is_safe(v: int, x: Fraction) -> bool:
    numerator, denominator = distance_numerator(v, x)
    return 14 * numerator >= denominator


def pool_wall_set() -> list[Fraction]:
    answer = {Fraction(0), Fraction(1)}
    for v in P:
        for k in range(v):
            answer.add(Fraction(14 * k + 1, 14 * v))
            answer.add(Fraction(14 * k + 13, 14 * v))
    return sorted(answer)


def danger_intervals(v: int):
    # Open real lifts, clipped to [0,1].  Only intersections with open pool
    # cells are tested here; pool wall points are checked separately.
    for k in range(v + 1):
        left = max(Fraction(0), Fraction(14 * k - 1, 14 * v))
        right = min(Fraction(1), Fraction(14 * k + 1, 14 * v))
        if left < right:
            yield left, right


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(newline="\n")
    hi = int(sys.argv[1]) if len(sys.argv) > 1 else 400
    event = pool_wall_set()
    safe_cells = []
    for left, right in zip(event, event[1:]):
        mid = (left + right) / 2
        if all(is_safe(v, mid) for v in P):
            safe_cells.append((left, right))
    safe_points = [x for x in event if all(is_safe(v, x) for v in P)]
    cell_endpoints = {x for left, right in safe_cells for x in (left, right)}
    if set(safe_points) != cell_endpoints:
        raise RuntimeError("safe wall points are not exactly the safe-cell endpoints")
    components = []
    for left, right in safe_cells:
        if components and components[-1][1] == left:
            components[-1] = (components[-1][0], right)
        else:
            components.append((left, right))
    largest = max(right - left for left, right in components)
    global_scan_bound = floor(Fraction(6, 7) / largest)

    implied = []
    hostile = {}
    for h in range(1, hi + 1):
        witness = None
        intervals = list(danger_intervals(h))
        for a, b in safe_cells:
            for left, right in intervals:
                lo, hi2 = max(a, left), min(b, right)
                if lo < hi2:
                    witness = (lo + hi2) / 2
                    break
            if witness is not None:
                break
        if witness is None:
            for x in safe_points:
                if not is_safe(h, x):
                    witness = x
                    break
        if witness is None:
            implied.append(h)
        else:
            if not all(is_safe(v, witness) for v in P) or is_safe(h, witness):
                raise RuntimeError(f"bad hostile witness h={h}, x={witness}")
            hostile[h] = witness

    expected = [v for v in P if v <= hi]
    print(f"POOL_INTERVAL_AUDIT range=1..{hi} walls={len(event)} "
          f"safe_cells={len(safe_cells)} safe_points={len(safe_points)} "
          f"components={len(components)} largest={largest} "
          f"global_scan_bound={global_scan_bound}")
    print("IMPLIED", ','.join(map(str, implied)))
    print("EXPECTED", ','.join(map(str, expected)))
    for h in (1, 2, 3, 7, 9, 11, 14, 17, 50, 291, 400):
        if h <= hi and h in hostile:
            print(f"HOSTILE h={h} x={hostile[h]}")
    if implied != expected:
        raise RuntimeError("closure mismatch")
    if hi >= global_scan_bound:
        print("GLOBAL_CLOSURE_CERTIFIED",
              "G_P_subset_G_h iff h_in_P for every positive integer h")
    print("PASS")


if __name__ == "__main__":
    main()
