#!/usr/bin/env python3
"""Exact probe for joint constraint implications from the THM-4326 pool.

For a speed v, D_v={x: ||vx||<1/14}.  A pool subset B implies a label h,
G_B subset G_h, exactly when D_h subset union_{b in B} D_b.  We construct
the complete rational wall arrangement for P union {h}, retain both open
cells and wall points, and solve the resulting 30-variable hitting set.

This is the primary combined-arrangement implementation.  The companion
fixed-pool interval program uses a disjoint construction.
"""

from __future__ import annotations

from fractions import Fraction
from functools import lru_cache
import sys


P = (8, 10, 15, 16, 20, 30, 40, 42, 60, 63, 80, 84, 85, 88, 95,
     120, 126, 132, 143, 145, 168, 170, 176, 190, 193, 240, 252,
     264, 286, 290)


def danger(v: int, x: Fraction) -> bool:
    # Exact strict test 14 ||v x|| < 1 on a rational representative.
    numerator = (v * x.numerator) % x.denominator
    return 14 * min(numerator, x.denominator - numerator) < x.denominator


def walls(v: int) -> set[Fraction]:
    answer: set[Fraction] = set()
    for k in range(v):
        answer.add(Fraction(14 * k + 1, 14 * v))
        answer.add(Fraction(14 * k + 13, 14 * v))
    return answer


POOL_WALLS = set().union(*(walls(v) for v in P), {Fraction(0), Fraction(1)})


def minimal_requirements(h: int) -> tuple[int, ...] | None:
    """Return inclusion-minimal nonzero pool-danger masks on D_h atoms."""
    event = sorted(POOL_WALLS | walls(h))
    probes: list[Fraction] = []
    # Wall points and open cells both matter because all dangers are strict.
    probes.extend(event[:-1])
    probes.extend((a + b) / 2 for a, b in zip(event, event[1:]))
    requirements: set[int] = set()
    for x in probes:
        if not danger(h, x):
            continue
        mask = 0
        for index, speed in enumerate(P):
            if danger(speed, x):
                mask |= 1 << index
        if mask == 0:
            return None
        requirements.add(mask)
    # If A subset B, hitting A automatically hits B.  Keep only minimal A.
    ordered = sorted(requirements, key=lambda m: (m.bit_count(), m))
    minimal: list[int] = []
    for mask in ordered:
        if not any((old & mask) == old for old in minimal):
            minimal.append(mask)
    return tuple(minimal)


def minimum_hitting_set(requirements: tuple[int, ...], cap: int = 9) -> tuple[int, ...] | None:
    full_labels = (1 << len(P)) - 1

    @lru_cache(maxsize=None)
    def search(unhit: tuple[int, ...], chosen: int, room: int) -> int | None:
        if not unhit:
            return chosen
        if room == 0:
            return None
        # Branch on the tightest currently unhit requirement.
        target = min(unhit, key=lambda m: (m & ~chosen).bit_count())
        options = target & ~chosen & full_labels
        while options:
            bit = options & -options
            options -= bit
            nxt = tuple(mask for mask in unhit if not (mask & bit))
            got = search(nxt, chosen | bit, room - 1)
            if got is not None:
                return got
        return None

    for size in range(cap + 1):
        got = search(requirements, 0, size)
        if got is not None:
            return tuple(P[i] for i in range(len(P)) if got & (1 << i))
    return None


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(newline="\n")
    hi = int(sys.argv[1]) if len(sys.argv) > 1 else 400
    print(f"POOL_IMPLICATION_CLOSURE range=1..{hi} pool_size={len(P)} cap=9")
    implied: list[tuple[int, tuple[int, ...], int]] = []
    union_fail: list[int] = []
    over_cap: list[tuple[int, int]] = []
    for h in range(1, hi + 1):
        req = minimal_requirements(h)
        if req is None:
            union_fail.append(h)
            continue
        witness = minimum_hitting_set(req, 9)
        if witness is None:
            over_cap.append((h, len(req)))
            continue
        implied.append((h, witness, len(req)))
        if h not in P:
            print(f"NEW h={h} cover_size={len(witness)} B={witness} minimal_masks={len(req)}")
    print("IMPLIED", ','.join(str(h) for h, _, _ in implied))
    print("NEW_COUNT", sum(h not in P for h, _, _ in implied))
    print("UNION_FAIL_COUNT", len(union_fail))
    print("OVER_CAP_COUNT", len(over_cap))


if __name__ == "__main__":
    main()
