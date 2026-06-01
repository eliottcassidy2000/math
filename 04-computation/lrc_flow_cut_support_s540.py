#!/usr/bin/env python3
"""Support-flow and cut-flow probes for LRC.

codex-2026-06-01 S540

This continues the nowhere-zero-flow attack beyond HYP-2025's full-support
speed-dipole resonance.

Two directed-graph views are audited side by side:

1. Fourier/support-flow side.
   The exact lonely-time measure is

       sum_{m : sum m_i v_i = 0} prod_i ghat_n(m_i),

   where the support S={i:m_i!=0} is a nowhere-zero integer flow on the
   speed-weighted sub-dipole induced by S.  HYP-2025 studied the full support;
   here we keep every support layer and measure cancellation by layer.

2. Cover/cut-flow side.
   Each danger interval B_v={t: ||v t|| < 1/n} is a directed arc on the time
   circle.  On the endpoint arrangement, the coverage count is a nonnegative
   flow on the directed cell cycle.  A positive-length lonely interval is a
   zero-flow cut.  A wall-only extremal such as the AP has nowhere-zero
   coverage on open cells but an empty strict endpoint-protection core.

Tournament Analysis declaration:
   vertices:
       support layers 0..n-1 of the speed-dipole flow enumerator;
   pairwise observable:
       absolute Fourier contribution, modular NZ-flow support count, and
       bridge-free subset count;
   switch/gauge:
       lexicographic dominance of support-layer pressure;
   tie Hamiltonian path:
       increasing support size;
   fingerprints:
       score histogram, directed triangles, SCCs, and Hamiltonian paths for
       the support-layer tournament.

Stored output:
   05-knowledge/results/lrc_flow_cut_support_s540.out
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from math import gcd, pi, sin


Adj = tuple[tuple[int, ...], ...]
ZERO = Fraction(0)
ONE = Fraction(1)


@dataclass(frozen=True)
class Case:
    name: str
    n: int
    speeds: tuple[int, ...]
    fourier_m: int


CASES = (
    Case("n4_AP_wall", 4, (1, 2, 3), 18),
    Case("n4_odd_sum_open", 4, (1, 2, 4), 18),
    Case("n5_AP_wall", 5, (1, 2, 3, 4), 14),
    Case("n5_open_sample", 5, (1, 2, 3, 5), 14),
    Case("n6_AP_wall", 6, (1, 2, 3, 4, 5), 10),
    Case("n6_Z3_bridge", 6, (1, 3, 6, 9, 12), 10),
    Case("n6_two_active", 6, (1, 2, 3, 6, 9), 10),
    Case("n7_AP_wall", 7, (1, 2, 3, 4, 5, 6), 7),
    Case("n7_open_sample", 7, (1, 2, 3, 4, 6, 8), 7),
    Case("n14_AP_wall", 14, tuple(range(1, 14)), 3),
)


def frac(x: Fraction) -> Fraction:
    return x - Fraction(x.numerator // x.denominator)


def dist0(x: Fraction) -> Fraction:
    f = frac(x)
    return min(f, ONE - f)


def nstar(n: int) -> int:
    return n // 2 if n % 2 == 0 else n


def ghat(m: int, n: int) -> float:
    """Fourier coefficient of 1[||x|| >= 1/n]."""
    if m == 0:
        return 1.0 - 2.0 / n
    return -sin(2 * pi * m / n) / (pi * m)


def exact_safe_measure_and_best(speeds: tuple[int, ...], n: int) -> tuple[Fraction, Fraction]:
    cuts = {ZERO, ONE}
    threshold = Fraction(1, n)
    for s in speeds:
        for k in range(s + 1):
            for sign in (-1, 1):
                t = Fraction(n * k + sign, n * s)
                if ZERO <= t <= ONE:
                    cuts.add(t)
    pts = sorted(cuts)
    total = ZERO
    best = ZERO
    probes = set(pts)
    probes.update((a + b) / 2 for a, b in zip(pts, pts[1:]) if a < b)
    for t in probes:
        margin = min(dist0(Fraction(s) * t) for s in speeds)
        if margin > best:
            best = margin
    for a, b in zip(pts, pts[1:]):
        if a >= b:
            continue
        mid = (a + b) / 2
        if all(dist0(Fraction(s) * mid) >= threshold for s in speeds):
            total += b - a
    return total, best


def add_interval_mod(out: list[tuple[Fraction, Fraction, int]], left: Fraction, right: Fraction, owner: int) -> None:
    """Add an interval modulo 1, split into non-wrapping pieces."""
    while left < ZERO:
        out.append((left + ONE, ONE, owner))
        left += ONE
        right += ONE
    while right > ONE:
        out.append((ZERO, right - ONE, owner))
        right -= ONE
        left -= ONE
    if left < right:
        out.append((left, right, owner))


def danger_intervals(speeds: tuple[int, ...], n: int) -> list[tuple[Fraction, Fraction, int]]:
    threshold = Fraction(1, n)
    intervals: list[tuple[Fraction, Fraction, int]] = []
    for idx, s in enumerate(speeds):
        for k in range(s):
            left = (Fraction(k) - threshold) / s
            right = (Fraction(k) + threshold) / s
            add_interval_mod(intervals, left, right, idx)
    return intervals


def inside_open_interval(p: Fraction, interval: tuple[Fraction, Fraction, int]) -> bool:
    left, right, _ = interval
    return left < p < right


def cover_cut_profile(speeds: tuple[int, ...], n: int) -> dict[str, object]:
    intervals = danger_intervals(speeds, n)
    endpoints = sorted({ZERO, ONE, *(x for a, b, _ in intervals for x in (a, b))})
    coverage_lengths: Counter[int] = Counter()
    zero_cells = 0
    min_cov: int | None = None
    for a, b in zip(endpoints, endpoints[1:]):
        if a >= b:
            continue
        mid = (a + b) / 2
        cov = sum(1 for interval in intervals if inside_open_interval(mid, interval))
        coverage_lengths[cov] += int((b - a) * 10**9)  # scaled only for a compact histogram
        zero_cells += int(cov == 0)
        min_cov = cov if min_cov is None else min(min_cov, cov)

    active = set(range(len(intervals)))
    changed = True
    while changed:
        changed = False
        remove = set()
        for idx in active:
            left, right, _ = intervals[idx]
            left_protected = any(
                j != idx and j in active and inside_open_interval(left, intervals[j])
                for j in active
            )
            right_protected = any(
                j != idx and j in active and inside_open_interval(right, intervals[j])
                for j in active
            )
            if not (left_protected and right_protected):
                remove.add(idx)
        if remove:
            active -= remove
            changed = True

    return {
        "intervals": len(intervals),
        "cells": max(0, len(endpoints) - 1),
        "zero_cells": zero_cells,
        "min_coverage": min_cov if min_cov is not None else 0,
        "coverage_hist": tuple(sorted((k, v) for k, v in coverage_lengths.items())),
        "strict_core_intervals": len(active),
    }


def support_flow_ledger(speeds: tuple[int, ...], n: int, bound: int) -> tuple[dict[int, float], dict[int, int]]:
    """Truncated Fourier sum by support size and bounded integer NZ-flow count."""
    vals: dict[tuple[int, int], float] = {(0, 0): 1.0}
    counts: dict[tuple[int, int], int] = {(0, 0): 1}
    for v in speeds:
        new_vals: defaultdict[tuple[int, int], float] = defaultdict(float)
        new_counts: defaultdict[tuple[int, int], int] = defaultdict(int)
        for (total, supp), value in vals.items():
            new_vals[(total, supp)] += value * ghat(0, n)
            new_counts[(total, supp)] += counts[(total, supp)]
            for m in range(-bound, bound + 1):
                if m == 0:
                    continue
                key = (total + m * v, supp + 1)
                new_vals[key] += value * ghat(m, n)
                new_counts[key] += counts[(total, supp)]
        vals = dict(new_vals)
        counts = dict(new_counts)
    ledger = {s: value for (total, s), value in vals.items() if total == 0}
    count_ledger = {s: value for (total, s), value in counts.items() if total == 0}
    return ledger, count_ledger


def nz_flow_count_mod(weights: tuple[int, ...], k: int) -> int:
    if not weights:
        return 1
    dist = {0: 1}
    for w in weights:
        nd: defaultdict[int, int] = defaultdict(int)
        for residue, count in dist.items():
            for m in range(1, k):
                nd[(residue + m * w) % k] += count
        dist = dict(nd)
    return dist.get(0, 0)


def modular_subset_profile(speeds: tuple[int, ...], n: int) -> dict[int, tuple[int, int, int]]:
    """support size -> (subsets, subsets with NZ mod-n* flow, subsets with >=2 active)."""
    k = nstar(n)
    out: dict[int, list[int]] = {s: [0, 0, 0] for s in range(len(speeds) + 1)}
    for r in range(len(speeds) + 1):
        for idxs in combinations(range(len(speeds)), r):
            weights = tuple(speeds[i] for i in idxs)
            active = sum(1 for w in weights if w % k != 0)
            out[r][0] += 1
            out[r][1] += int(nz_flow_count_mod(weights, k) > 0)
            out[r][2] += int(active >= 2)
    return {key: tuple(value) for key, value in out.items()}


def tournament_from_layer_scores(scores: list[tuple[int, int, int]]) -> Adj:
    q = len(scores)
    adj = [[0] * q for _ in range(q)]
    for i, j in combinations(range(q), 2):
        if scores[i] > scores[j] or (scores[i] == scores[j] and i < j):
            adj[i][j] = 1
        else:
            adj[j][i] = 1
    return tuple(tuple(row) for row in adj)


def score_hist(adj: Adj) -> tuple[tuple[int, int], ...]:
    return tuple(sorted(Counter(sum(row) for row in adj).items()))


def directed_triangles(adj: Adj) -> int:
    total = 0
    for a, b, c in combinations(range(len(adj)), 3):
        if (adj[a][b] and adj[b][c] and adj[c][a]) or (
            adj[a][c] and adj[c][b] and adj[b][a]
        ):
            total += 1
    return total


def scc_sizes(adj: Adj) -> tuple[int, ...]:
    n = len(adj)
    graph = [[j for j in range(n) if adj[i][j]] for i in range(n)]
    rgraph = [[] for _ in range(n)]
    for i, row in enumerate(graph):
        for j in row:
            rgraph[j].append(i)
    seen: set[int] = set()
    order: list[int] = []

    def dfs(v: int) -> None:
        seen.add(v)
        for w in graph[v]:
            if w not in seen:
                dfs(w)
        order.append(v)

    for v in range(n):
        if v not in seen:
            dfs(v)
    seen.clear()
    sizes = []

    def rdfs(v: int) -> int:
        seen.add(v)
        total = 1
        for w in rgraph[v]:
            if w not in seen:
                total += rdfs(w)
        return total

    for v in reversed(order):
        if v not in seen:
            sizes.append(rdfs(v))
    return tuple(sorted(sizes, reverse=True))


def hamiltonian_paths(adj: Adj) -> int:
    n = len(adj)
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for mask in range(1 << n):
        for last in range(n):
            cur = dp[mask][last]
            if not cur:
                continue
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                if adj[last][nxt]:
                    dp[mask | (1 << nxt)][nxt] += cur
    return sum(dp[full])


def fmt_fraction(x: Fraction) -> str:
    return f"{x} ({float(x):.6f})"


def main() -> None:
    print("LRC support-flow and cut-flow probe -- codex-2026-06-01 S540")
    print("=" * 78)
    print("Flow side: |SAFE| is a support-decomposed NZ-flow enumerator over sub-dipoles.")
    print("Cut side: danger arcs induce a nonnegative cover flow on the time-cell cycle.")
    print()

    for case in CASES:
        safe, best = exact_safe_measure_and_best(case.speeds, case.n)
        cut = cover_cut_profile(case.speeds, case.n)
        ledger, counts = support_flow_ledger(case.speeds, case.n, case.fourier_m)
        mod_profile = modular_subset_profile(case.speeds, case.n)
        total = sum(ledger.values())
        main = ledger.get(0, 0.0)
        nonzero = total - main
        layer_scores = []
        print(f"CASE {case.name}: n={case.n}, speeds={case.speeds}, M={case.fourier_m}")
        print("-" * 78)
        print(
            f"exact open |SAFE|={fmt_fraction(safe)}; best_margin={fmt_fraction(best)}; "
            f"closed_LRC={best >= Fraction(1, case.n)}"
        )
        print(
            "cover-cut: intervals={intervals} cells={cells} zero_cells={zero_cells} "
            "min_coverage={min_coverage} strict_core_intervals={strict_core_intervals}".format(**cut)
        )
        print(f"support-flow truncated total={total:+.8f}; main={main:+.8f}; nonzero={nonzero:+.8f}")
        print("support layers: k -> contribution, bounded-flow-count, modNZ/subsets, active>=2")
        for support in range(len(case.speeds) + 1):
            contribution = ledger.get(support, 0.0)
            bounded_count = counts.get(support, 0)
            subsets, mod_nz, active_ge2 = mod_profile.get(support, (0, 0, 0))
            layer_scores.append(
                (
                    int(round(abs(contribution) * 10**12)),
                    mod_nz,
                    active_ge2,
                )
            )
            if contribution or bounded_count or subsets:
                print(
                    f"  {support:2d}: {contribution:+.8f} "
                    f"count={bounded_count:7d} modNZ={mod_nz:4d}/{subsets:<4d} "
                    f"active>=2={active_ge2:4d}"
                )
        adj = tournament_from_layer_scores(layer_scores)
        print(
            f"layer-tournament: score_hist={score_hist(adj)} c3={directed_triangles(adj)} "
            f"SCC={scc_sizes(adj)} H={hamiltonian_paths(adj)}"
        )
        if safe == 0 and best >= Fraction(1, case.n):
            print("reading: wall-only cover flow; open cells have no zero cut, but closed LRC survives on a boundary.")
        elif safe > 0:
            print("reading: positive zero-flow cut exists on the cover cycle; nonzero support flows fail to cancel main.")
        else:
            print("reading: this would be counterexample-shaped only if best_margin also fell below threshold.")
        print()

    print("SYNTHESIS")
    print("=" * 78)
    print("1. HYP-2025's full-support NZ flow is one layer of a larger support-flow")
    print("   enumerator: every nonzero Fourier support is a NZ flow on a sub-dipole.")
    print("2. On the dual cover graph, a positive lonely interval is a zero-flow cut")
    print("   in the danger-arc coverage flow.  AP walls are nowhere-zero on open")
    print("   cells but have empty strict endpoint-protection cores.")
    print("3. A true counterexample would need both sides at once: full cancellation")
    print("   of the support-flow enumerator and a nonpeeling nowhere-zero cover-flow")
    print("   core with best margin below 1/n.  That is much more rigid than merely")
    print("   having many full-support resonances.")


if __name__ == "__main__":
    main()
