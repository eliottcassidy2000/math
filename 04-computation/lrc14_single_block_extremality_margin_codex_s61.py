#!/usr/bin/env python3
"""HYP-2694 / S61: exact single-block wide-cover extremality ledger.

This script isolates the quotient used in HYP-2694:

    E_M = {0} union {M, M+1, ..., M+m-1},   m = k-1.

In the decorrelated model the fixed speed 0 is its own anchored cluster, while
the m nonzero speeds are split into one or more far coherent blocks.  All blocks
share the same slow coordinate x and have independent carrier phases.  The
finite question is therefore:

    among integer partitions of m, does the one-part partition [m] maximize
    the exact shared-x decorrelated cover?

The script proves this finite quotient exactly for m=7..11 by enumerating all
integer partitions and integrating over the exact carrier breakpoints.  It also
records a single-block finite-scale error bound:

    |p0(E_M) - p0_decorr([m])| <= J_m / M,
    J_m = 7 * binom(m, 2),

which follows by freezing x on each 1/M cell while the carrier phase Mx runs
once around the circle; for fixed carrier phase, the block cover indicator has
at most J_m jumps as a function of x.

This closes the large-M single-block branch once J_m/M is below the cap margin.
The multi-block analogue is still the live HYP-2684/HYP-2694 discrepancy task.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from functools import lru_cache
from fractions import Fraction as F
from itertools import combinations
from math import ceil, gcd


CAP = {
    8: F(2243, 5880),
    9: F(1979, 4004),
    10: F(55, 91),
    11: F(66, 91),
    12: F(6, 7),
}

INNER = frozenset(range(1, 7))


def fmt(q: F) -> str:
    return f"{q} ({float(q):.9f})"


@lru_cache(maxsize=None)
def breakpoints(row: tuple[int, ...]) -> tuple[F, ...]:
    bps = {F(0), F(1)}
    for e in row:
        if e == 0:
            continue
        for a in range(7 * e + 1):
            bps.add(F(a, 7 * e))
    return tuple(sorted(b for b in bps if 0 <= b <= 1))


@lru_cache(maxsize=None)
def far_dist_at_x(block: tuple[int, ...], x: F) -> tuple[tuple[frozenset[int], F], ...]:
    """Law of covered inner-sector set for one far block at fixed slow x."""
    base = tuple((c * x) % 1 for c in block)
    tb = {F(0), F(1)}
    for b in base:
        for s in range(7):
            tb.add((F(s, 7) - b) % 1)
    cuts = sorted(tb)
    dist: dict[frozenset[int], F] = defaultdict(F)
    for lo, hi in zip(cuts, cuts[1:]):
        if hi <= lo:
            continue
        mid = (lo + hi) / 2
        hit = frozenset(
            s for s in (int(((b + mid) % 1) * 7) for b in base) if 1 <= s <= 6
        )
        dist[hit] += hi - lo
    return tuple(sorted(dist.items(), key=lambda kv: (len(kv[0]), tuple(kv[0]))))


def anchored_set_at_x(block: tuple[int, ...], x: F) -> frozenset[int]:
    return frozenset(s for s in (int((c * x) % 1 * 7) for c in block) if 1 <= s <= 6)


@lru_cache(maxsize=None)
def p0_decorr(clusters: tuple[tuple[int, ...], ...]) -> F:
    """Exact shared-x decorrelated cover for anchored cluster plus far clusters."""
    xcuts = set()
    for block in clusters:
        xcuts.update(breakpoints(block))
    xs = sorted(xcuts)
    total = F(0)
    for lo, hi in zip(xs, xs[1:]):
        if hi <= lo:
            continue
        mid = (lo + hi) / 2
        cur: dict[frozenset[int], F] = {anchored_set_at_x(clusters[0], mid): F(1)}
        for block in clusters[1:]:
            nxt: dict[frozenset[int], F] = defaultdict(F)
            for have, mass in cur.items():
                for hit, weight in far_dist_at_x(block, mid):
                    nxt[have | hit] += mass * weight
            cur = dict(nxt)
        total += (hi - lo) * sum(mass for have, mass in cur.items() if INNER <= have)
    return total


def p0_exact(row: tuple[int, ...]) -> F:
    """Exact p0 for an actual speed row."""
    xs = breakpoints(tuple(sorted(set(row))))
    total = F(0)
    for lo, hi in zip(xs, xs[1:]):
        if hi <= lo:
            continue
        mid = (lo + hi) / 2
        hit = frozenset(s for s in (int((e * mid) % 1 * 7) for e in row) if 1 <= s <= 6)
        if INNER <= hit:
            total += hi - lo
    return total


def partitions(n: int, max_part: int | None = None) -> list[tuple[int, ...]]:
    """Integer partitions of n in nonincreasing order."""
    if max_part is None:
        max_part = n
    if n == 0:
        return [()]
    out: list[tuple[int, ...]] = []
    for first in range(min(n, max_part), 0, -1):
        for rest in partitions(n - first, first):
            out.append((first,) + rest)
    return out


def clusters_for_part(part: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    """Anchor 0 is separate; each nonzero part is a coherent block {0,...,s-1}."""
    return ((0,),) + tuple(tuple(range(s)) for s in part)


def primitive(row: tuple[int, ...]) -> bool:
    g = 0
    for e in row:
        g = gcd(g, e)
    return g == 1


def finite_single_block_row(m: int, M: int) -> tuple[int, ...]:
    return (0,) + tuple(range(M, M + m))


def part_extremality() -> dict[int, tuple[F, F, tuple[int, ...]]]:
    print("PART A -- exact decorrelated partition-function extremality")
    print("vertices are integer partitions of m=k-1 nonzero runners; anchor 0 is fixed")
    print(
        f"{'k':>3} {'m':>3} {'#parts':>7} {'best':>12} {'p0_decorr(best)':>28} "
        f"{'cap-margin':>24} {'closest split':>18} {'split gap':>24}"
    )
    out: dict[int, tuple[F, F, tuple[int, ...]]] = {}
    for k in range(8, 13):
        m = k - 1
        rows = [(p0_decorr(clusters_for_part(part)), part) for part in partitions(m)]
        rows.sort(key=lambda item: (item[0], item[1]), reverse=True)
        best_val, best_part = rows[0]
        assert best_part == (m,)
        split_val, split_part = rows[1]
        margin = CAP[k] - best_val
        gap = best_val - split_val
        out[k] = (best_val, margin, split_part)
        print(
            f"{k:>3} {m:>3} {len(rows):>7} {str(best_part):>12} {fmt(best_val):>28} "
            f"{fmt(margin):>24} {str(split_part):>18} {fmt(gap):>24}"
        )
    print()
    print("TOURNAMENT ANALYSIS")
    print("  vertices: integer partitions of m, not runners or arcs")
    print("  pairwise observable: exact p0_decorr(partition); larger value wins")
    print("  switch/gauge: scalar comparison in the shared-x carrier quotient")
    print("  tie Hamiltonian path: partitions sorted by exact p0_decorr")
    for k in range(8, 13):
        m = k - 1
        vals = [(p0_decorr(clusters_for_part(part)), part) for part in partitions(m)]
        scores = {}
        ties = 0
        for val_a, part_a in vals:
            wins = 0
            for val_b, part_b in vals:
                if part_a == part_b:
                    continue
                if val_a > val_b:
                    wins += 1
                elif val_a == val_b:
                    ties += 1
            scores[part_a] = wins
        hist = Counter(scores.values())
        directed_cycles = 0
        hpaths = 1 if ties == 0 else "tie-family"
        print(
            f"  k={k}: score_hist={dict(sorted(hist.items()))} "
            f"directed_3cycles={directed_cycles} Hamiltonian_path_count={hpaths} ties={ties // 2}"
        )
    print()
    return out


def part_single_block_error(ext: dict[int, tuple[F, F, tuple[int, ...]]]) -> None:
    print("PART B -- single-block finite-scale error budget")
    print("For E_M={0} union {M,...,M+m-1}, J_m=7*binom(m,2) and")
    print("|p0(E_M)-p0_decorr([m])| <= J_m/M by the cell-freeze BV argument.")
    print(
        f"{'k':>3} {'m':>3} {'J_m':>5} {'margin':>24} {'M for J/M<margin':>20} "
        f"{'observed M':>10} {'p0(E_M)-limit':>24} {'within margin?':>15}"
    )
    for k in range(8, 13):
        m = k - 1
        limit, margin, _split = ext[k]
        jumps = 7 * m * (m - 1) // 2
        cutoff = jumps * margin.denominator // margin.numerator + 1
        # The named HYP-2694 row uses M=19 at k=8.  Use max(19,m+8) uniformly.
        M = max(19, m + 8)
        row = finite_single_block_row(m, M)
        value = p0_exact(row)
        err = value - limit
        print(
            f"{k:>3} {m:>3} {jumps:>5} {fmt(margin):>24} {cutoff:>20} "
            f"{M:>10} {fmt(err):>24} {str(abs(err) < margin):>15}"
        )
        assert primitive(row)
        assert value < CAP[k]
    print()


def part_finite_samples(ext: dict[int, tuple[F, F, tuple[int, ...]]]) -> None:
    print("PART C -- exact finite shifted-block samples")
    print("These are not the proof; they show the BV cutoff is very conservative.")
    for k in range(8, 13):
        m = k - 1
        limit, margin, _split = ext[k]
        print(f"k={k}, m={m}, limit={fmt(limit)}, cap-margin={fmt(margin)}")
        for M in (15, 19, 30, 100):
            if M <= m:
                continue
            row = finite_single_block_row(m, M)
            value = p0_exact(row)
            print(
                f"  M={M:>3}: p0={fmt(value)} err={fmt(value - limit)} "
                f"cap-p0={fmt(CAP[k] - value)}"
            )
        print()


def main() -> None:
    print("HYP-2694 / S61 -- single-block extremality, cap margin, and error budget")
    print("All displayed partition-function values are exact Fractions.")
    print()
    ext = part_extremality()
    part_single_block_error(ext)
    part_finite_samples(ext)
    print("SYNTHESIS")
    print("  1. In the HYP-2694 quotient, the one-part block [m] is the exact")
    print("     decorrelated maximizer for every m=7..11; every split loses.")
    print("  2. The exact single-block cover is below cap_k with margin at least")
    print("     the k=9 margin shown above, about 0.18862.")
    print("  3. For actual shifted single blocks, the cell-freeze BV lemma gives")
    print("     a rigorous large-M closure once 7*binom(m,2)/M is below that margin.")
    print("  4. Remaining HYP-2694 work is the multi-block/multi-carrier version:")
    print("     prove the analogous joint decorrelation error is bounded by the")
    print("     available split gap plus cap margin, then finite-check the small gaps.")


if __name__ == "__main__":
    main()
