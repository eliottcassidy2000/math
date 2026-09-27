#!/usr/bin/env python3
"""Whole-orbit landing multiplicities as chain lengths of the D-coarsened excursion
order (session note section 6.6).  For an orbit (x_j) and depth D, the dipper i lands
at the first j > i with x_j < 2^-D x_i; m(j) = number of dippers landing at j.
THM-4506 (S): all dippers of j lie in (2^D x_j, 2^(D+1) x_j].  We report, per depth,
the number of landing points, the mean and maximum multiplicity, and (for the 5x+1
control) how the mean grows with the log-size of the values.

Run: python 04-computation/experiments/collatz_posets_dags_20260927_chains.py
"""
from __future__ import annotations

import math
from collections import Counter


def orbit(n: int, q: int, steps: int | None) -> list[int]:
    out = [n]
    x = n
    while True:
        if steps is not None and len(out) > steps:
            break
        if steps is None and x == 1:
            break
        x = (q * x + 1) // 2 if x % 2 else x // 2
        out.append(x)
    return out


def multiplicities(vals: list[int], D: int) -> tuple[Counter, list[int]]:
    """Landing counts m(j) using a monotone stack: dippers waiting for their first
    D-drop form a stack ordered by threshold; when x_j < 2^-D x_i for the stack top
    it also holds for every deeper entry with a larger threshold?  Not monotone in
    general, so do it directly but prune with the shell lemma: a dipper i can only
    land while x_j > 2^-(D+1) x_i ... we simply scan forward with a list of pending
    dippers, which is fast enough for 2*10^4 steps."""
    n = len(vals)
    land = [None] * n
    pending: list[int] = []
    for j in range(n):
        # resolve pending dippers whose first D-drop is j
        still = []
        for i in pending:
            if vals[j] * (1 << D) < vals[i]:
                land[i] = j
            else:
                still.append(i)
        pending = still
        pending.append(j)
    counts = Counter(l for l in land if l is not None)
    # shell lemma check
    for i, l in enumerate(land):
        if l is not None:
            assert vals[i] > vals[l] * (1 << D) and vals[i] <= vals[l] * (1 << (D + 1))
    return counts, land


def report(name: str, vals: list[int], depths: list[int]) -> None:
    print(f"-- {name}: {len(vals) - 1} steps, log2(max) = {math.log2(max(vals)):.1f}")
    for D in depths:
        counts, land = multiplicities(vals, D)
        if not counts:
            print(f"   D={D}: no landings")
            continue
        ms = list(counts.values())
        dist = Counter(ms)
        print(f"   D={D}: landings {len(ms):>6}  dippers {sum(ms):>6}  mean mult {sum(ms)/len(ms):.3f}  max {max(ms):>3}  "
              f"dist {dict(sorted(dist.items()))}" if len(dist) <= 12 else
              f"   D={D}: landings {len(ms):>6}  dippers {sum(ms):>6}  mean mult {sum(ms)/len(ms):.3f}  max {max(ms):>3}")


if __name__ == "__main__":
    print("== whole-orbit landing multiplicities (chain lengths of the D-coarsened excursion order) ==")
    for n in [27, 703, 6171, 77031, 837799, 8400511, 63728127, 670617279]:
        report(f"3x+1 orbit of {n}", orbit(n, 3, None), [1, 2, 3, 4, 6])
    v = orbit(7, 5, 20000)
    report("5x+1 orbit of 7 (20000 steps)", v, [1, 2, 3, 4, 6, 8, 12])
    # growth of the mean with scale: split the 5x+1 orbit into stretches by log2 of the value
    print("-- 5x+1 orbit of 7: mean multiplicity per stretch of 2500 steps (D = 2, 4, 8)")
    for D in (2, 4, 8):
        counts, land = multiplicities(v, D)
        row = []
        for s in range(0, 20000, 2500):
            js = [j for j in counts if s <= j < s + 2500]
            if js:
                mean = sum(counts[j] for j in js) / len(js)
                L = math.log2(v[s + 1250])
                row.append(f"[{s},{s+2500}): L={L:.0f} mean={mean:.2f} sqrtL={math.sqrt(L):.1f}")
        print(f"   D={D}: " + "; ".join(row))
    print("DONE")
