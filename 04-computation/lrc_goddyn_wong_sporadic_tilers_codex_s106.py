#!/usr/bin/env python3
"""Goddyn-Wong sporadic tiler audit.

Convention: n moving speeds have LRC threshold 1/(n+1).  The AP tiler is
{1,...,n}.  A Goddyn-Wong acceleration replaces selected speeds v by 2v.

The external criterion quoted by Tao: replacing one velocity v in {1,...,n}
by 2v preserves tightness when v shares a nontrivial common factor with every
integer in [n-v+1, 2n-2v+1].

This script treats that criterion as the structural lens, then exact-checks
representative one- and multi-accelerated rows by envelope maximization.
"""

from __future__ import annotations

from fractions import Fraction
from math import gcd


F = Fraction


def norm_frac(x: F) -> F:
    r = x % 1
    return min(r, 1 - r)


def envelope_candidates(speeds: list[int]) -> set[F]:
    speeds = sorted(set(speeds))
    out: set[F] = {F(1, 2)}

    for v in speeds:
        k = 0
        while True:
            t = F(2 * k + 1, 2 * v)
            if t > F(1, 2):
                break
            out.add(t)
            k += 1

    for i, a in enumerate(speeds):
        for b in speeds[i + 1 :]:
            for d in (a + b, b - a):
                if d <= 0:
                    continue
                k = 1
                while True:
                    t = F(k, d)
                    if t > F(1, 2):
                        break
                    out.add(t)
                    k += 1
    return out


def lonely_constant(speeds: list[int]) -> tuple[F, F]:
    best = F(0)
    arg = F(0)
    for t in envelope_candidates(speeds):
        val = min(norm_frac(F(v) * t) for v in speeds)
        if val > best:
            best = val
            arg = t
    return best, arg


def safe_measure(speeds: list[int], q: int) -> F:
    """Exact measure of {t: ||v t|| >= 1/q for all v}.

    Used only on moderate rows.  Breakpoints are the endpoints of the open
    danger arcs ||v t|| < 1/q.
    """
    bps: set[F] = {F(0), F(1)}
    for v in speeds:
        a = abs(v)
        for m in range(a + 1):
            for sign in (-1, 1):
                x = F(q * m + sign, q * a)
                if 0 <= x <= 1:
                    bps.add(x)

    ordered = sorted(bps)
    total = F(0)
    for lo, hi in zip(ordered, ordered[1:]):
        mid = (lo + hi) / 2
        if all(norm_frac(F(v) * mid) >= F(1, q) for v in speeds):
            total += hi - lo
    return total


def gw_accelerable_velocities(n: int) -> list[int]:
    """Velocities v whose Tao/Goddyn-Wong non-coprime interval test holds."""
    good: list[int] = []
    for v in range(n // 2 + 1, n + 1):
        lo = n - v + 1
        hi = 2 * n - 2 * v + 1
        if lo <= hi and all(gcd(v, j) > 1 for j in range(lo, hi + 1)):
            good.append(v)
    return good


def accelerated_tuple(n: int, accelerations: list[int] | None = None) -> list[int]:
    if accelerations is None:
        accelerations = gw_accelerable_velocities(n)
    return sorted([x for x in range(1, n + 1) if x not in accelerations] + [2 * v for v in accelerations])


def noncoprime_window(n: int, v: int) -> tuple[int, int, list[int]]:
    lo = n - v + 1
    hi = 2 * n - 2 * v + 1
    return lo, hi, [gcd(v, j) for j in range(lo, hi + 1)]


def boundary_profile(speeds: list[int], q: int) -> list[tuple[int, list[int], list[int]]]:
    """At q-grid points a/q, record immediate left/right blockers.

    A speed with a*s == 1 mod q blocks the left side of a/q; a speed with
    a*s == -1 mod q blocks the right side.  Rows are restricted to grid points
    at which no speed is exactly at the origin.
    """
    rows = []
    for a in range(1, q):
        residues = [(a * s) % q for s in speeds]
        if 0 in residues:
            continue
        left = [s for s in speeds if (a * s) % q == 1]
        right = [s for s in speeds if (a * s) % q == q - 1]
        rows.append((a, left, right))
    return rows


def print_selected_verifications() -> None:
    selected = [7, 13, 19, 32, 73]
    print("SELECTED EXACT TIGHTNESS CHECKS")
    for n in selected:
        q = n + 1
        acc = gw_accelerable_velocities(n)
        speeds = accelerated_tuple(n, acc)
        m, arg = lonely_constant(speeds)
        target = F(1, q)
        safe = safe_measure(speeds, q) if n <= 32 else None
        safe_text = str(safe) if safe is not None else "skipped-large"
        print(
            f"n={n:3d} q={q:3d} acc={acc!s:10s} max={max(speeds):3d} "
            f"M={m} target={target} arg={arg} safe={safe_text} tight={m == target}"
        )


def print_acceleration_census(limit: int = 150) -> None:
    print("\nGODDYN-WONG ACCELERATION CENSUS")
    for n in range(2, limit + 1):
        acc = gw_accelerable_velocities(n)
        if not acc:
            continue
        windows = []
        for v in acc:
            lo, hi, g = noncoprime_window(n, v)
            windows.append(f"{v}:[{lo},{hi}] gcds={g}")
        speeds = accelerated_tuple(n, acc)
        tail = [x for x in speeds if x >= max(1, n - 8)]
        print(f"n={n:3d} q={n+1:3d} acc={acc!s:12s} tail={tail!s:38s} windows={' ; '.join(windows)}")


def print_lrc14_profile() -> None:
    n = 13
    q = 14
    speeds = accelerated_tuple(n)
    print("\nLRC14 GODDYN-WONG ROW")
    print(f"speeds={speeds}")
    print("missing from AP:", sorted(set(range(1, n + 1)) - set(speeds)))
    print("new speeds:", sorted(set(speeds) - set(range(1, n + 1))))
    lo, hi, g = noncoprime_window(n, 12)
    print(f"accelerated v=12 because every j in [{lo},{hi}] has gcd(12,j)>1: {g}")
    print("q-grid immediate blocker profile (a/q, left blockers, right blockers):")
    for a, left, right in boundary_profile(speeds, q):
        print(f"  a={a:2d}: left={left!s} right={right!s}")


def print_tournament_analysis() -> None:
    carriers = [
        ("jacobsthal_window", 5),
        ("single_acceleration", 4),
        ("multi_acceleration", 3),
        ("difference_closed_AP", 2),
        ("arbitrary_sporadic_search", 1),
    ]
    print("\nTOURNAMENT ANALYSIS")
    print("vertices are proof carriers, not runners.")
    print("observable: retained structure for classifying exact tilers.")
    print("Hamiltonian path:", " > ".join(name for name, _ in sorted(carriers, key=lambda x: -x[1])))


if __name__ == "__main__":
    print_selected_verifications()
    print_acceleration_census()
    print_lrc14_profile()
    print_tournament_analysis()
