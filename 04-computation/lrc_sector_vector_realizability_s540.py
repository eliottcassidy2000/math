#!/usr/bin/env python3
"""
Realizable sector-vectors for LRC.

codex-2026-06-01-S540

Sector vector:
    c(t) = (c_0,...,c_{n-1}), c_k = # runners in [k/n,(k+1)/n).

This script separates three meanings of "realizable":

  1. existentially realizable:
       some primitive speed set and some open time has sector vector c;

  2. bounded-clock reachable:
       c appears among speed sets with max speed <= B;

  3. forced:
       c appears for every primitive speed set in a bounded box.

The main punchline is deliberately deflationary: every composition of n-1 into
n sectors is existentially realizable.  Thus the restriction in HYP-2022 cannot
mean global existence of sector-vectors.  It lives in fixed-clock menus, low
complexity, boundary compactification, target anchoring, and forced hitting of
the observer-empty face c_0=c_{n-1}=0.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
from functools import reduce
from itertools import combinations
from math import comb, gcd


ZERO = Fraction(0)
ONE = Fraction(1)


def frac(x: Fraction) -> Fraction:
    return x - Fraction(x.numerator // x.denominator)


def compositions(total: int, parts: int):
    if parts == 1:
        yield (total,)
        return
    for x in range(total + 1):
        for rest in compositions(total - x, parts - 1):
            yield (x,) + rest


def clasp_reflect(c: tuple[int, ...]) -> tuple[int, ...]:
    """Reflection fixing the observer clasp between sectors n-1 and 0."""
    n = len(c)
    return tuple(c[n - 1 - k] for k in range(n))


def canon_clasp(c: tuple[int, ...]) -> tuple[int, ...]:
    r = clasp_reflect(c)
    return min(c, r)


def canon_dihedral(c: tuple[int, ...]) -> tuple[int, ...]:
    n = len(c)
    rotations = []
    for s in range(n):
        r = c[s:] + c[:s]
        rotations.append(r)
        rotations.append(tuple(reversed(r)))
    return min(rotations)


def is_good(c: tuple[int, ...]) -> bool:
    return c[0] == 0 and c[-1] == 0


def occupancy(speeds: tuple[int, ...], n: int, t: Fraction) -> tuple[int, ...]:
    counts = [0] * n
    for v in speeds:
        counts[int(n * frac(Fraction(v) * t)) % n] += 1
    return tuple(counts)


def sector_walls(speeds: tuple[int, ...], n: int) -> list[Fraction]:
    walls = {ZERO, ONE}
    for v in speeds:
        for k in range(n * v):
            walls.add(Fraction(k, n * v))
    return sorted(walls)


def open_cell_midpoints(speeds: tuple[int, ...], n: int):
    walls = sector_walls(speeds, n)
    for a, b in zip(walls, walls[1:]):
        if a < b:
            yield (a + b) / 2


def primitive_speed_sets_with_max(n: int, max_speed: int):
    """Yield distinct primitive speed sets of length n-1 with maximum max_speed."""
    if max_speed < n - 1:
        return
    for rest in combinations(range(1, max_speed), n - 2):
        speeds = tuple(rest) + (max_speed,)
        if reduce(gcd, speeds) == 1:
            yield speeds


def choose_witness_speeds(c: tuple[int, ...]):
    """Construct distinct primitive speeds and t=1/(nL) realizing c in open sectors."""
    n = len(c)
    total = sum(c)
    assert total == n - 1

    # q=nL.  Sector k contains open residues kL+1,...,(k+1)L-1.
    # We backtrack only until gcd becomes 1; n is small in this audit.
    for L in range(max(c) + 3, max(c) + 80):
        pools = []
        for k, count in enumerate(c):
            candidates = list(range(k * L + 1, (k + 1) * L))
            if len(candidates) < count:
                break
            for _ in range(count):
                pools.append((k, candidates))
        else:
            used_by_sector = {k: set() for k in range(n)}

            def search(pos: int, chosen: list[int], current_gcd: int):
                if pos == len(pools):
                    return tuple(sorted(chosen)) if current_gcd == 1 else None
                k, candidates = pools[pos]
                ordered = sorted(
                    candidates,
                    key=lambda r: (gcd(current_gcd, r) != 1 if current_gcd else True, r),
                )
                for r in ordered:
                    if r in used_by_sector[k]:
                        continue
                    used_by_sector[k].add(r)
                    ng = r if current_gcd == 0 else gcd(current_gcd, r)
                    ans = search(pos + 1, chosen + [r], ng)
                    if ans is not None:
                        return ans
                    used_by_sector[k].remove(r)
                return None

            speeds = search(0, [], 0)
            if speeds is not None:
                t = Fraction(1, n * L)
                assert reduce(gcd, speeds) == 1
                assert occupancy(speeds, n, t) == c
                return speeds, t, L
    raise RuntimeError(f"failed to construct witness for {c}")


def bounded_reachable(n: int, bound: int):
    seen_first: dict[tuple[int, ...], tuple[int, tuple[int, ...], Fraction]] = {}
    speed_sets = 0
    certified_sets = 0
    universal = None

    for max_speed in range(1, bound + 1):
        for speeds in primitive_speed_sets_with_max(n, max_speed):
            speed_sets += 1
            local = set()
            local_good = False
            for t in open_cell_midpoints(speeds, n):
                c = occupancy(speeds, n, t)
                local.add(c)
                local_good |= is_good(c)
                if c not in seen_first:
                    seen_first[c] = (max_speed, speeds, t)
            if local_good:
                certified_sets += 1
            universal = set(local) if universal is None else universal & local
    return seen_first, speed_sets, certified_sets, (universal or set())


def summarize_n(n: int, bound: int):
    all_vectors = list(compositions(n - 1, n))
    good_vectors = [c for c in all_vectors if is_good(c)]
    clasp_orbits = {canon_clasp(c) for c in all_vectors}
    dihedral_orbits = {canon_dihedral(c) for c in all_vectors}
    good_clasp_orbits = {canon_clasp(c) for c in good_vectors}

    witnesses = {c: choose_witness_speeds(c) for c in all_vectors}
    max_constructed = max(max(speeds) for speeds, _t, _L in witnesses.values())
    worst = sorted(
        ((max(speeds), c, speeds, t) for c, (speeds, t, _L) in witnesses.items()),
        reverse=True,
    )[:3]

    seen_first, speed_sets, certified_sets, universal = bounded_reachable(n, bound)
    seen = set(seen_first)
    seen_good = {c for c in seen if is_good(c)}
    first_hist = Counter(v[0] for v in seen_first.values())
    forced_good = sum(1 for c in universal if is_good(c))

    print(f"n={n}")
    print(f"  all raw sector-vectors:        {len(all_vectors)} = C({2*n-2},{n-1})")
    print(f"  observer-clasp reflection orbits: {len(clasp_orbits)}")
    print(f"  free dihedral orbits:          {len(dihedral_orbits)}")
    print(f"  good raw vectors c0=c[n-1]=0:  {len(good_vectors)}")
    print(f"  good clasp-reflection orbits:  {len(good_clasp_orbits)}")
    print(
        f"  constructive theorem:          {len(witnesses)}/{len(all_vectors)} "
        f"vectors realized, max witness speed <= {max_constructed}"
    )
    print(f"  bounded exact search B<={bound}:")
    print(
        f"    speed sets={speed_sets}, certified by open good vector="
        f"{certified_sets}/{speed_sets}"
    )
    print(
        f"    vectors seen={len(seen)}/{len(all_vectors)}, "
        f"good seen={len(seen_good)}/{len(good_vectors)}, "
        f"forced intersection={len(universal)} (good={forced_good})"
    )
    print(f"    first-seen max-speed histogram: {dict(sorted(first_hist.items()))}")
    if len(universal) <= 10:
        print(f"    forced vectors: {sorted(universal)}")
    print("  worst constructive witnesses:")
    for m, c, speeds, t in worst:
        print(f"    max={m:3d}, c={c}, speeds={speeds}, t={t}")
    print()


def main():
    print("Sector-vector realizability -- codex S540")
    print("=" * 78)
    print("A sector-vector is c=(c_0,...,c_{n-1}), sum c_i=n-1.")
    print("Good means c_0=c_{n-1}=0, i.e. the observer's two sectors are empty.")
    print()

    bounds = {3: 10, 4: 12, 5: 14, 6: 16, 7: 14}
    for n, bound in bounds.items():
        summarize_n(n, bound)

    print("Interpretation")
    print("-" * 78)
    print("Existential sector-vector realizability is trivial in the strong sense:")
    print("  every composition of n-1 into n sectors is realized by an explicit")
    print("  primitive speed set at an open time t=1/(nL).")
    print("Therefore HYP-2022's restriction is not global vector existence.")
    print("The real questions are fixed-clock menus, low-complexity realization,")
    print("observer anchoring, boundary compactification, and forced hitting of")
    print("the good face c_0=c_{n-1}=0 by every primitive clock.")
    print("The bounded forced intersections above are bad boundary fans, not LRC")
    print("witnesses: they come from times near 0 or 1 where almost all runners sit")
    print("next to the observer.  The good face is not a single forced vector.")


if __name__ == "__main__":
    main()
