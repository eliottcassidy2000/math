#!/usr/bin/env python3
"""
lrc_observer_source_tournament_s511.py

oracle-2026-06-01-S511

LRC as a MARKED tournament-iso-class problem -- fixing the "safety is not a class
function" wall (codex S509/HYP-1977).

codex lifted runner movies to half-turn tournaments and found LRC safety is NOT a
function of the (even observer-marked) half-turn class.  The reason is the WRONG
WALLS: the half-turn comparator meters the 1/2 gap, not the LRC 1/n gap, so it
cannot see loneliness.

THIS construction uses the LRC walls and an observer-anchored comparator so that
loneliness becomes a PURE marked-class property:

  vertices: observer (speed 0) + runners (speeds v_1..v_{n-1}); threshold 1/n.
  observer--runner i:  observer -> i   iff   ||v_i t|| >= 1/n   (runner i is FAR,
                       i.e. in the safe arc [1/n, 1-1/n]); else i -> observer.
  runner i--runner j:  half-turn  (i -> j iff frac((v_i - v_j) t) in (0,1/2)).

Then, exactly:   observer is LONELY at t  <=>  observer beats every runner
                                          <=>  observer is a SOURCE (out-deg n-1).

So LRC( observer )  <=>  the observer-marked tournament clock reaches a cell whose
class has the observer as a source.  That is a statement purely about marked
tournament isomorphism classes (A000568 with one marked vertex).  This script:
 (1) verifies observer-source <=> LRC-safe (0 mismatches),
 (2) records, per cell, the marked iso-class and whether observer is a source,
 (3) asks the new question: is observer-source a function of the marked class?
     (yes, by construction) and of the UNMARKED runner-runner class? (the real
     content -- does the half-turn runner sub-tournament constrain loneliness?).
"""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations, permutations
from math import gcd
from collections import defaultdict


ONE = Fraction(1)


def frac(x: Fraction) -> Fraction:
    return x - (x.numerator // x.denominator)


def dist0(x: Fraction) -> Fraction:
    f = frac(x)
    return min(f, ONE - f)


def marked_tournament(speeds, n, t):
    """speeds = (0, v_1, ..., v_{n-1}); vertex 0 = observer. Return adj (n x n)."""
    thr = Fraction(1, n)
    adj = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            if i == 0 or j == 0:
                # observer edge: observer -> runner iff runner is FAR (||v t||>=1/n)
                runner = j if i == 0 else i
                far = dist0(Fraction(speeds[runner]) * t) >= thr
                if i == 0:
                    adj[0][j] = 1 if far else 0
                else:
                    adj[i][0] = 0 if far else 1
            else:
                f = frac(Fraction(speeds[i] - speeds[j]) * t)
                adj[i][j] = 1 if 0 < f < Fraction(1, 2) else 0
    return tuple(tuple(r) for r in adj)


def walls(speeds, n):
    thr = Fraction(1, n)
    W = set()
    # observer-runner LRC endpoint walls: ||v t|| = 1/n
    for v in speeds[1:]:
        if v == 0:
            continue
        # v t = m +/- 1/n
        for m in range(0, v):
            W.add(Fraction(m, v) + thr / v)
            W.add(Fraction(m, v) - thr / v)
    # runner-runner half-turn walls: frac((vi-vj)t) in {0,1/2}
    for i, j in combinations(range(1, n), 2):
        d = abs(speeds[i] - speeds[j])
        if d == 0:
            continue
        for k in range(0, 2 * d):
            W.add(Fraction(k, 2 * d))
    W = sorted(frac(w) for w in W)
    # dedupe within [0,1)
    out = []
    for w in W:
        if 0 <= w < 1 and (not out or out[-1] != w):
            out.append(w)
    return out


def observer_is_source(adj):
    return all(adj[0][j] == 1 for j in range(1, len(adj)))


def lonely_direct(speeds, n, t):
    thr = Fraction(1, n)
    return all(dist0(Fraction(v) * t) >= thr for v in speeds[1:])


def canon_marked(adj):
    """canonical form fixing vertex 0 (observer); minimize over perms of 1..n-1."""
    n = len(adj)
    best = None
    for perm in permutations(range(1, n)):
        p = (0,) + perm
        flat = tuple(adj[p[i]][p[j]] for i in range(n) for j in range(n) if i != j)
        if best is None or flat < best:
            best = flat
    return best


def canon_unmarked_runners(adj):
    """canonical form of the runner-runner sub-tournament (drop observer)."""
    n = len(adj)
    idx = list(range(1, n))
    best = None
    for perm in permutations(idx):
        flat = tuple(adj[perm[a]][perm[b]] for a in range(len(idx)) for b in range(len(idx)) if a != b)
        if best is None or flat < best:
            best = flat
    return best


def Hcount(adj):
    n = len(adj)
    full = (1 << n) - 1
    from functools import lru_cache

    @lru_cache(None)
    def dp(mask, last):
        if mask == full:
            return 1
        s = 0
        for nx in range(n):
            if not (mask >> nx) & 1 and adj[last][nx]:
                s += dp(mask | (1 << nx), nx)
        return s
    return sum(dp(1 << s, s) for s in range(n))


def analyze(speeds):
    speeds = tuple(speeds)
    n = len(speeds)
    W = walls(speeds, n)
    W2 = W + [ONE]
    cells = []
    for a, b in zip(W2, W2[1:]):
        mid = (a + b) / 2
        cells.append(mid)
    mism = 0
    source_classes = set()
    nonsource_classes = set()
    # map: unmarked runner class -> set of source bits seen
    runnerclass_source = defaultdict(set)
    n_source_cells = 0
    for t in cells:
        adj = marked_tournament(speeds, n, t)
        src = observer_is_source(adj)
        lon = lonely_direct(speeds, n, t)
        if src != lon:
            mism += 1
        if src:
            n_source_cells += 1
        mc = canon_marked(adj)
        (source_classes if src else nonsource_classes).add(mc)
        rc = canon_unmarked_runners(adj)
        runnerclass_source[rc].add(src)
    # is observer-source determined by the UNMARKED runner class?
    runner_mixed = sum(1 for v in runnerclass_source.values() if len(v) > 1)
    return {
        "n": n, "cells": len(cells), "mism": mism,
        "source_cells": n_source_cells,
        "marked_source_classes": len(source_classes),
        "marked_nonsource_classes": len(nonsource_classes),
        "marked_overlap": len(source_classes & nonsource_classes),
        "runner_classes": len(runnerclass_source),
        "runner_mixed": runner_mixed,
    }


def main():
    print("LRC observer-source marked tournament (oracle-2026-06-01-S511)\n")
    print("Claim: observer is a SOURCE  <=>  observer is LONELY (LRC-safe).")
    print("=> LRC is a reachability question on observer-MARKED tournament classes.\n")
    fams = {
        "n4 initial (0,1,2,3)": (0, 1, 2, 3),
        "n4 (0,1,2,5)": (0, 1, 2, 5),
        "n5 initial (0,1,2,3,4)": (0, 1, 2, 3, 4),
        "n5 (0,2,3,5,7)": (0, 2, 3, 5, 7),
        "n6 initial (0,1,2,3,4,5)": (0, 1, 2, 3, 4, 5),
        "n6 (0,1,3,4,5,9)": (0, 1, 3, 4, 5, 9),
    }
    print(f"{'family':<26} cells mism src_cells | mk_src mk_nonsrc mk_overlap | runnerC mixed")
    for label, s in fams.items():
        r = analyze(s)
        flag = "OK" if r["mism"] == 0 else f"MISMATCH={r['mism']}"
        print(f"{label:<26} {r['cells']:<5} {flag:<5} {r['source_cells']:<9} | "
              f"{r['marked_source_classes']:<6} {r['marked_nonsource_classes']:<9} "
              f"{r['marked_overlap']:<10} | {r['runner_classes']:<7} {r['runner_mixed']}")
    print()
    print("Reading:")
    print(" mism=0 everywhere  => observer-source = LRC-safe exactly (by design):")
    print("   LRC = 'the observer-marked clock reaches a SOURCE class'. A pure")
    print("   marked tournament-iso-class (A000568-with-mark) reachability problem.")
    print(" mk_overlap=0       => observer-source IS a function of the marked class")
    print("   (sanity: source-ness is a marked-class invariant).")
    print(" runner_mixed>0     => loneliness is NOT a function of the UNMARKED")
    print("   runner-runner half-turn class: the half-turn sub-tournament does not")
    print("   determine loneliness (matches codex's negative result). The MARK +")
    print("   the LRC walls are exactly what carries the loneliness information.")


if __name__ == "__main__":
    main()
