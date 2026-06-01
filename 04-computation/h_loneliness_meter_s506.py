#!/usr/bin/env python3
"""
h_loneliness_meter_s506.py

codex-2026-06-01 S506

Investigate H(T) as a loneliness meter for runner-clock tournaments.

There are two related but different notions:

1. Unanchored half-turn gap:
   the half-turn circular tournament is transitive (H=1) exactly when all
   points lie in an open semicircle, i.e. max circular gap > 1/2.

2. LRC marked loneliness:
   for n total runners including the stationary runner 0, runner 0 is lonely
   at time t when every moving runner is at distance >= 1/n from 0.  More
   generally, a runner is locally lonely iff the two adjacent circular gaps
   around it are both >= 1/n (THM-370).

The script compares H to both measurements over exact clock cells and selected
LRC witness times.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
from math import sqrt
from pathlib import Path
from importlib.machinery import SourceFileLoader


ROOT = Path(__file__).resolve().parents[1]
clock = SourceFileLoader(
    "tournament_clock_s24",
    str(ROOT / "04-computation" / "tournament_clock_s24.py"),
).load_module()


def mod1(x: Fraction) -> Fraction:
    return x - (x.numerator // x.denominator)


def norm_dist(x: Fraction) -> Fraction:
    y = mod1(x)
    return min(y, 1 - y)


def fmt_frac(x: Fraction) -> str:
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def tournament_at_with_tie_path(speeds: tuple[int, ...], t: Fraction) -> tuple[int, ...]:
    """Return adjacency bitmasks. Ties at 0 or 1/2 use label path i -> j for i<j."""
    n = len(speeds)
    out = [0] * n
    for i, j in combinations(range(n), 2):
        f = mod1(Fraction(speeds[i] - speeds[j]) * t)
        if Fraction(0) < f < Fraction(1, 2):
            out[i] |= 1 << j
        elif f in (Fraction(0), Fraction(1, 2)):
            out[i] |= 1 << j
        else:
            out[j] |= 1 << i
    return tuple(out)


@lru_cache(maxsize=None)
def h_count(adj_bits: tuple[int, ...]) -> int:
    n = len(adj_bits)
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        unvisited_mask = full ^ mask
        for last in range(n):
            val = dp[mask][last]
            if not val:
                continue
            nxts = adj_bits[last] & unvisited_mask
            while nxts:
                bit = nxts & -nxts
                nxt = bit.bit_length() - 1
                dp[mask | bit][nxt] += val
                nxts ^= bit
    return sum(dp[full])


def score_sequence(adj_bits: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sorted(x.bit_count() for x in adj_bits))


def positions(speeds: tuple[int, ...], t: Fraction) -> list[tuple[int, Fraction]]:
    return [(i, mod1(Fraction(v) * t)) for i, v in enumerate(speeds)]


def circular_gaps(speeds: tuple[int, ...], t: Fraction) -> tuple[list[tuple[int, Fraction]], list[Fraction]]:
    ordered = sorted(positions(speeds, t), key=lambda item: (item[1], item[0]))
    gaps = []
    for idx, (_, here) in enumerate(ordered):
        nxt = ordered[(idx + 1) % len(ordered)][1]
        gaps.append(mod1(nxt - here))
    return ordered, gaps


def geometry_record(speeds: tuple[int, ...], t: Fraction) -> dict[str, object]:
    n = len(speeds)
    threshold = Fraction(1, n)
    ordered, gaps = circular_gaps(speeds, t)
    label_to_idx = {label: i for i, (label, _) in enumerate(ordered)}
    lonely = []
    for label, _ in ordered:
        idx = label_to_idx[label]
        if gaps[(idx - 1) % n] >= threshold and gaps[idx] >= threshold:
            lonely.append(label)
    stationary_min = min(norm_dist(Fraction(v) * t) for v in speeds[1:])
    return {
        "max_gap": max(gaps),
        "min_gap": min(gaps),
        "lonely_vertices": tuple(sorted(lonely)),
        "stationary_lonely": stationary_min >= threshold,
        "stationary_min": stationary_min,
        "threshold": threshold,
        "safe_gap_count": sum(1 for g in gaps if g >= threshold),
    }


def pearson(xs: list[float], ys: list[float]) -> float:
    if len(xs) < 2:
        return 0.0
    mx = sum(xs) / len(xs)
    my = sum(ys) / len(ys)
    vx = sum((x - mx) ** 2 for x in xs)
    vy = sum((y - my) ** 2 for y in ys)
    if vx == 0 or vy == 0:
        return 0.0
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sqrt(vx * vy)


def ranks(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=lambda i: values[i])
    out = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i + 1
        while j < len(order) and values[order[j]] == values[order[i]]:
            j += 1
        rank = (i + j - 1) / 2
        for k in range(i, j):
            out[order[k]] = rank
        i = j
    return out


def spearman(xs: list[float], ys: list[float]) -> float:
    return pearson(ranks(xs), ranks(ys))


def exact_cells(speeds: tuple[int, ...]) -> list[tuple[Fraction, tuple[int, ...], dict[str, object]]]:
    cells = []
    for t, _ in clock.clock_cells(speeds):
        adj = tournament_at_with_tie_path(speeds, t)
        cells.append((t, adj, geometry_record(speeds, t)))
    return cells


def bucket_report(label: str, speeds: tuple[int, ...]) -> None:
    cells = exact_cells(speeds)
    by_h: dict[int, list[dict[str, object]]] = defaultdict(list)
    hs = []
    max_gaps = []
    stationary_mins = []
    lonely_counts = []
    for _, adj, geom in cells:
        h = h_count(adj)
        by_h[h].append(geom)
        hs.append(float(h))
        max_gaps.append(float(geom["max_gap"]))
        stationary_mins.append(float(geom["stationary_min"]))
        lonely_counts.append(float(len(geom["lonely_vertices"])))

    print(f"[{label}] speeds={speeds} n={len(speeds)} cells={len(cells)}")
    print(
        "  correlations: "
        f"Pearson(H,max_gap)={pearson(hs, max_gaps):+.3f} "
        f"Spearman={spearman(hs, max_gaps):+.3f}; "
        f"Pearson(H,stationary_min)={pearson(hs, stationary_mins):+.3f}; "
        f"Pearson(H,#lonely_vertices)={pearson(hs, lonely_counts):+.3f}"
    )
    print("  H buckets: H cells max_gap_mean[min,max] stat_lonely avg_lonely avg_safe_gaps")
    prev_mean: float | None = None
    inversions = 0
    for h in sorted(by_h):
        rows = by_h[h]
        mg = [float(r["max_gap"]) for r in rows]
        stat = sum(1 for r in rows if r["stationary_lonely"])
        avg_lonely = sum(len(r["lonely_vertices"]) for r in rows) / len(rows)
        avg_safe = sum(int(r["safe_gap_count"]) for r in rows) / len(rows)
        mean_mg = sum(mg) / len(mg)
        if prev_mean is not None and mean_mg > prev_mean + 1e-12:
            inversions += 1
        prev_mean = mean_mg
        print(
            f"    {h:<8} {len(rows):<5} "
            f"{mean_mg:.4f}[{min(mg):.4f},{max(mg):.4f}] "
            f"{stat:<4} {avg_lonely:>6.2f} {avg_safe:>6.2f}"
        )
    print(f"  bucket max-gap mean inversions as H increases: {inversions}")


def stationary_candidate_times(speeds: tuple[int, ...]) -> list[Fraction]:
    n = len(speeds)
    threshold = Fraction(1, n)
    times = {Fraction(0), Fraction(1)}
    for v in speeds[1:]:
        for m in range(v):
            times.add(mod1((Fraction(m) + threshold) / v))
            times.add(mod1((Fraction(m) - threshold) / v))
    ordered = sorted(times)
    mids = set()
    wrap = ordered + [ordered[0] + 1]
    for a, b in zip(wrap, wrap[1:]):
        mids.add(mod1((a + b) / 2))
    return sorted(times | mids)


def best_stationary_time(speeds: tuple[int, ...]) -> tuple[Fraction, Fraction, tuple[int, ...], dict[str, object]]:
    best_t = Fraction(0)
    best_dist = Fraction(-1)
    best_adj: tuple[int, ...] | None = None
    best_geom: dict[str, object] | None = None
    for t in stationary_candidate_times(speeds):
        geom = geometry_record(speeds, t)
        dist = geom["stationary_min"]
        if dist > best_dist:
            best_dist = dist
            best_t = t
            best_adj = tournament_at_with_tie_path(speeds, t)
            best_geom = geom
    assert best_adj is not None and best_geom is not None
    return best_t, best_dist, best_adj, best_geom


def lrc_witness_report(label: str, speeds: tuple[int, ...]) -> None:
    t, dist, adj, geom = best_stationary_time(speeds)
    n = len(speeds)
    threshold = Fraction(1, n)
    h = h_count(adj)
    print(f"[{label}] speeds={speeds} n={n}")
    print(
        f"  best stationary time t={fmt_frac(t)} "
        f"min_dist0={fmt_frac(dist)} threshold={fmt_frac(threshold)} "
        f"margin={fmt_frac(dist - threshold)} stationary_lonely={geom['stationary_lonely']}"
    )
    print(
        f"  half-turn H={h} score={score_sequence(adj)} "
        f"max_gap={fmt_frac(geom['max_gap'])} "
        f"lonely_vertices={geom['lonely_vertices']} "
        f"safe_gap_count={geom['safe_gap_count']}"
    )


def selected_lrc_time_report(label: str, speeds: tuple[int, ...], times: list[Fraction]) -> None:
    print(f"[{label}] selected times, speeds={speeds} n={len(speeds)}")
    for t in times:
        adj = tournament_at_with_tie_path(speeds, t)
        geom = geometry_record(speeds, t)
        print(
            f"  t={fmt_frac(t):>10} H={h_count(adj):>18} "
            f"max_gap={fmt_frac(geom['max_gap']):>8} "
            f"stat_min={fmt_frac(geom['stationary_min']):>8} "
            f"stat_lonely={str(geom['stationary_lonely']):<5} "
            f"lonely={geom['lonely_vertices']}"
        )


def main() -> None:
    print("H as a loneliness meter for runner-clock tournaments")
    print("=" * 88)
    print("Part A: exact clock-cell audits")
    print()
    families = [
        ("initial n=5", (0, 1, 2, 3, 4)),
        ("primes n=5", (0, 2, 3, 5, 7)),
        ("spread n=5", (0, 1, 4, 9, 11)),
        ("initial n=6", (0, 1, 2, 3, 4, 5)),
        ("primes n=6", (0, 2, 3, 5, 7, 11)),
        ("initial n=7", (0, 1, 2, 3, 4, 5, 6)),
        ("spread n=7", (0, 1, 3, 7, 12, 18, 27)),
    ]
    for label, speeds in families:
        bucket_report(label, speeds)
        print()

    print("=" * 88)
    print("Part B: LRC marked stationary-witness times")
    print()
    for n in range(5, 10):
        lrc_witness_report(f"initial segment total n={n}", tuple(range(n)))
    print()
    lrc_witness_report("prime speeds total n=7", (0, 2, 3, 5, 7, 11, 13))
    lrc_witness_report("spread speeds total n=7", (0, 1, 3, 7, 12, 18, 27))
    print()

    print("=" * 88)
    print("Part C: selected larger LRC rows")
    print()
    selected_lrc_time_report(
        "n=14 initial/ladders",
        tuple(range(14)),
        [Fraction(1, 14), Fraction(1, 28), Fraction(3053, 25872), Fraction(4339, 51744)],
    )
    print()
    selected_lrc_time_report(
        "n=18 initial/ladders",
        tuple(range(18)),
        [Fraction(1, 18), Fraction(1, 36), Fraction(3991, 57024), Fraction(8681, 114048)],
    )
    print()

    print("SYNTHESIS")
    print("=" * 88)
    print("H has two loneliness readings.")
    print("  Low H, especially H=1, detects unanchored bunching: a huge empty semicircle.")
    print("  High H detects LRC-style even separation: many vertices have two safe gaps.")
    print("Plain H is therefore not the LRC invariant by itself; LRC needs a marked")
    print("tournament structure: circular order, stationary vertex, safe-gap mask, and")
    print("pressure/deletion data.  H is the scalar shadow of that richer tournament movie.")


if __name__ == "__main__":
    main()
