#!/usr/bin/env python3
"""
tournament_clock_s24.py

oracle-2026-06-01-S24

TOURNAMENT ANALYSIS, continued (builds on S22/S23).

The "tournament clock".  Take n runners with integer speeds s_0<...<s_{n-1} on a
unit circle.  At time t each runner is at x_i(t)=frac(s_i t).  Lift the
configuration to a tournament by the PHASE comparator:

    i -> j   iff   frac(x_i - x_j) = frac((s_i - s_j) t) in (0, 1/2)

("i leads j by less than half a lap").  Ties (frac in {0,1/2}) are broken by the
fixed Hamiltonian label path 0 -> 1 -> ... -> n-1 (the repo's base path; the
speed-rank labels).  As t runs over [0,1) the tournament T(t) is piecewise
constant and traces a CLOSED WALK through tournament space (the metagraph G_n).

Edge (i,j) flips exactly when frac((s_i-s_j)t) hits 0 or 1/2, i.e. at
    t = m / (2|s_i - s_j|),   m in Z.
So the wall-crossing times are exact rationals.  Between consecutive walls the
tournament is constant ("a cell").  We enumerate cells exactly and, per cell,
compute the tournament's H (number of directed Hamiltonian paths), score
sequence, and isomorphism class.  This turns a runner system into a *tournament
signature* and lets us hunt for patterns as the speeds change.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations, permutations
from math import gcd
from functools import lru_cache


# ----- tournament invariants -------------------------------------------------

def hamiltonian_path_count(adj: tuple[tuple[int, ...], ...]) -> int:
    """# of directed Hamiltonian paths in tournament given by adj[i][j]=1 if i->j."""
    n = len(adj)
    full = (1 << n) - 1
    from functools import lru_cache

    @lru_cache(maxsize=None)
    def dp(mask: int, last: int) -> int:
        if mask == full:
            return 1
        total = 0
        for nxt in range(n):
            if not (mask >> nxt) & 1 and adj[last][nxt]:
                total += dp(mask | (1 << nxt), nxt)
        return total

    out = 0
    for start in range(n):
        out += dp(1 << start, start)
    dp.cache_clear()
    return out


def score_sequence(adj) -> tuple[int, ...]:
    n = len(adj)
    return tuple(sorted(sum(adj[i]) for i in range(n)))


def canonical_form(adj) -> tuple:
    """Canonical (iso-invariant) key: lexicographically minimal over relabelings."""
    n = len(adj)
    best = None
    for perm in permutations(range(n)):
        flat = tuple(
            adj[perm[i]][perm[j]]
            for i in range(n) for j in range(n) if i != j
        )
        if best is None or flat < best:
            best = flat
    return best


# ----- the tournament clock --------------------------------------------------

def frac(x: Fraction) -> Fraction:
    return x - (x.numerator // x.denominator)


def tournament_at(speeds, t: Fraction):
    """Phase-comparator tournament at time t (t not on a wall)."""
    n = len(speeds)
    adj = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            f = frac(Fraction(speeds[i] - speeds[j]) * t)
            if 0 < f < Fraction(1, 2):
                adj[i][j] = 1
    return tuple(tuple(r) for r in adj)


def wall_times(speeds) -> list[Fraction]:
    walls = set()
    n = len(speeds)
    for i, j in combinations(range(n), 2):
        d = abs(speeds[i] - speeds[j])
        if d == 0:
            continue
        for m in range(0, 2 * d):
            walls.add(Fraction(m, 2 * d))
    return sorted(walls)


def clock_cells(speeds):
    """Return list of (t_mid, adj) for each cell of [0,1)."""
    walls = wall_times(speeds)
    walls.append(Fraction(1))
    cells = []
    for a, b in zip(walls, walls[1:]):
        mid = (a + b) / 2
        cells.append((mid, tournament_at(speeds, mid)))
    return cells


def clock_signature(speeds):
    cells = clock_cells(speeds)
    Hs = [hamiltonian_path_count(adj) for _, adj in cells]
    canons = [canonical_form(adj) for _, adj in cells]
    scores = [score_sequence(adj) for _, adj in cells]
    distinct_iso = len(set(canons))
    # iso-class walk (sequence of class ids around the clock)
    ids = {}
    walk = []
    for c in canons:
        ids.setdefault(c, len(ids))
        walk.append(ids[c])
    transitive_cells = [k for k, s in enumerate(scores)
                        if s == tuple(range(len(speeds)))]  # 0,1,..,n-1 = transitive
    return {
        "n": len(speeds),
        "n_cells": len(cells),
        "distinct_iso": distinct_iso,
        "H_multiset": tuple(sorted(set(Hs))),
        "H_per_cell": Hs,
        "n_transitive_cells": len(transitive_cells),
        "walk": walk,
        "cells": cells,
    }


def fmt_speeds(s):
    return "(" + ",".join(map(str, s)) + ")"


def report(label, speeds):
    sig = clock_signature(speeds)
    print(f"[{label}] speeds={fmt_speeds(speeds)}  n={sig['n']}")
    print(f"   cells={sig['n_cells']}  distinct_iso_classes={sig['distinct_iso']}  "
          f"distinct_H={sig['H_multiset']}")
    print(f"   transitive_cells={sig['n_transitive_cells']} "
          f"(times the speed-order is realized as a clean ranking)")
    # H histogram
    from collections import Counter
    hist = Counter(sig["H_per_cell"])
    print(f"   H histogram over cells: {dict(sorted(hist.items()))}")
    return sig


def main():
    print("Tournament clock — phase comparator (oracle-2026-06-01-S24)\n")

    print("=" * 68)
    print("A. n=5 (basketball-sized) speed families")
    print("=" * 68)
    fams5 = {
        "arithmetic 0..4": (0, 1, 2, 3, 4),
        "LRC initial 0,1,2,3,4": (0, 1, 2, 3, 4),   # observer + 1..4
        "geometric 0,1,2,4,8": (0, 1, 2, 4, 8),
        "primes 0,2,3,5,7": (0, 2, 3, 5, 7),
        "random-ish 0,1,4,9,11": (0, 1, 4, 9, 11),
        "resonant 0,1,2,3,6": (0, 1, 2, 3, 6),
    }
    sigs5 = {}
    for label, s in fams5.items():
        sigs5[label] = report(label, s)
    print()

    print("=" * 68)
    print("B. n=6 speed families")
    print("=" * 68)
    fams6 = {
        "arithmetic 0..5": (0, 1, 2, 3, 4, 5),
        "geometric 0,1,2,4,8,16": (0, 1, 2, 4, 8, 16),
        "primes 0,2,3,5,7,11": (0, 2, 3, 5, 7, 11),
        "lonely extremal 0..5": (0, 1, 2, 3, 4, 5),
        "coprime spread 0,1,5,7,11,13": (0, 1, 5, 7, 11, 13),
    }
    for label, s in fams6.items():
        report(label, s)
    print()

    print("=" * 68)
    print("C. n=7 (a couple)")
    print("=" * 68)
    for label, s in {
        "arithmetic 0..6": (0, 1, 2, 3, 4, 5, 6),
        "primes 0,2,3,5,7,11,13": (0, 2, 3, 5, 7, 11, 13),
    }.items():
        report(label, s)
    print()

    print("SUMMARY")
    print(" The tournament clock turns a runner system into a closed walk in G_n.")
    print(" Pattern probes: #cells vs sum|s_i-s_j|; distinct iso-classes vs speed")
    print(" arithmetic structure; whether the transitive tournament is visited")
    print(" (a clean global ranking moment) and how often.")


if __name__ == "__main__":
    main()
