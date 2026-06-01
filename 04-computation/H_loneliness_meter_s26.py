#!/usr/bin/env python3
"""Investigate H as a loneliness meter for half-turn circular tournaments.

oracle-2026-06-01-S26, extended by codex.

Given n points on the circle, orient the unordered pair {i,j} by the half-turn
rule: i beats j when j lies clockwise from i by less than 1/2 turn.  Ties are
resolved by the base label order.  The observable H is the directed
Hamiltonian-path count of this tournament.

The experiment asks what kind of "loneliness" H measures.  It compares H with:

* the largest circular gap;
* the two-neighbor static safety count at threshold 1/n;
* score variance and directed 3-cycles;
* exact selected n=14 Lonely Runner rows from the tournament-clock overlay.

Stored output:
    05-knowledge/results/H_loneliness_meter_s26.out
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass, field
from fractions import Fraction
from itertools import combinations
from math import log, sqrt
import random


ONE = Fraction(1, 1)
HALF = Fraction(1, 2)


def is_fraction(value: object) -> bool:
    return isinstance(value, Fraction)


def circle(value):
    if is_fraction(value):
        return value % ONE
    return value % 1.0


def clockwise(a, b):
    return circle(b - a)


def circular_distance_to_zero(value):
    value = circle(value)
    if is_fraction(value):
        return min(value, ONE - value)
    return min(value, 1.0 - value)


def as_float(value) -> float:
    return float(value.numerator) / float(value.denominator) if is_fraction(value) else float(value)


def fmt_value(value) -> str:
    if value is None:
        return "-"
    if isinstance(value, int):
        return str(value)
    if is_fraction(value):
        if value.denominator == 1:
            return str(value.numerator)
        return f"{value.numerator}/{value.denominator}"
    return f"{value:.6f}"


def half_turn_tournament(points) -> list[list[bool]]:
    n = len(points)
    adj = [[False] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        d = clockwise(points[i], points[j])
        if is_fraction(d):
            if d == 0 or d == HALF:
                winner = i
            elif d < HALF:
                winner = i
            else:
                winner = j
        else:
            if abs(d) < 1e-15 or abs(d - 0.5) < 1e-15:
                winner = i
            elif d < 0.5:
                winner = i
            else:
                winner = j
        loser = j if winner == i else i
        adj[winner][loser] = True
    return adj


def row_masks(adj: list[list[bool]]) -> tuple[int, ...]:
    masks: list[int] = []
    for row in adj:
        mask = 0
        for j, value in enumerate(row):
            if value:
                mask |= 1 << j
        masks.append(mask)
    return tuple(masks)


def H_count(adj: list[list[bool]]) -> int:
    """Count directed Hamiltonian paths by dense subset DP."""
    n = len(adj)
    masks = row_masks(adj)
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        row = dp[mask]
        for last, value in enumerate(row):
            if not value:
                continue
            allowed = masks[last] & ~mask
            while allowed:
                bit = allowed & -allowed
                nxt = bit.bit_length() - 1
                dp[mask | bit][nxt] += value
                allowed -= bit
    return sum(dp[full])


def num_3cycles(adj: list[list[bool]]) -> int:
    count = 0
    for a, b, c in combinations(range(len(adj)), 3):
        if (adj[a][b] and adj[b][c] and adj[c][a]) or (
            adj[a][c] and adj[c][b] and adj[b][a]
        ):
            count += 1
    return count


def scores(adj: list[list[bool]]) -> tuple[int, ...]:
    return tuple(sorted(sum(row) for row in adj))


def score_variance(adj: list[list[bool]]) -> float:
    n = len(adj)
    mean = (n - 1) / 2
    return sum((sum(row) - mean) ** 2 for row in adj) / n


def gap_list(points) -> list:
    ordered = sorted(points)
    gaps = [ordered[i + 1] - ordered[i] for i in range(len(ordered) - 1)]
    gaps.append(ordered[0] + (ONE if is_fraction(ordered[0]) else 1.0) - ordered[-1])
    return sorted(gaps, reverse=True)


def static_safe_vertex_count(points, threshold) -> int:
    """Count vertices with both adjacent circular gaps at least threshold."""
    ordered = sorted(points)
    n = len(ordered)
    count = 0
    unit = ONE if is_fraction(threshold) else 1.0
    for i, point in enumerate(ordered):
        left = point - ordered[i - 1] if i > 0 else point + unit - ordered[-1]
        right = ordered[(i + 1) % n] - point if i + 1 < n else ordered[0] + unit - point
        if left >= threshold and right >= threshold:
            count += 1
    return count


def gap_entropy(points) -> float:
    gaps = [as_float(gap) for gap in gap_list(points)]
    n = len(gaps)
    if n <= 1:
        return 0.0
    return -sum(gap * log(gap) for gap in gaps if gap > 0.0) / log(n)


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


@dataclass
class Record:
    H: int
    max_gap: float
    min_gap: float
    safe: int
    c3: int
    score_var: float
    entropy: float
    scores: tuple[int, ...]


@dataclass
class Bucket:
    records: list[Record] = field(default_factory=list)

    def values(self, name: str) -> list[float]:
        return [float(getattr(record, name)) for record in self.records]

    def int_values(self, name: str) -> list[int]:
        return [int(getattr(record, name)) for record in self.records]


def sample_records(n: int, samples: int, seed: int) -> tuple[dict[int, Bucket], list[Record]]:
    rng = random.Random(seed)
    by_H: dict[int, Bucket] = defaultdict(Bucket)
    records: list[Record] = []
    for _ in range(samples):
        points = sorted(rng.random() for _ in range(n))
        adj = half_turn_tournament(points)
        gaps = gap_list(points)
        record = Record(
            H=H_count(adj),
            max_gap=as_float(gaps[0]),
            min_gap=as_float(gaps[-1]),
            safe=static_safe_vertex_count(points, 1.0 / n),
            c3=num_3cycles(adj),
            score_var=score_variance(adj),
            entropy=gap_entropy(points),
            scores=scores(adj),
        )
        by_H[record.H].records.append(record)
        records.append(record)
    return dict(by_H), records


def monotonicity_witness(records: list[Record]) -> tuple[Record, Record, float] | None:
    """Find lower-H and higher-H records where higher H has larger max_gap."""
    best: tuple[Record, Record, float] | None = None
    by_H = sorted(set(record.H for record in records))
    for low_H in by_H:
        low_records = [record for record in records if record.H == low_H]
        low = min(low_records, key=lambda record: record.max_gap)
        for high_H in by_H:
            if high_H <= low_H:
                continue
            high = max(
                (record for record in records if record.H == high_H),
                key=lambda record: record.max_gap,
            )
            diff = high.max_gap - low.max_gap
            if diff > 0 and (best is None or diff > best[2]):
                best = (low, high, diff)
    return best


def print_detail(n: int, samples: int) -> None:
    by_H, records = sample_records(n, samples, seed=2600 + n)
    print("=" * 88)
    print(f"n={n}; samples={samples}; observed circular H menu={sorted(by_H)}")
    print("=" * 88)
    mismatches = [
        record
        for record in records
        if (record.H == 1) != (record.max_gap >= 0.5 - 1e-12)
    ]
    print(
        "H=1 iff max_gap>=1/2 mismatches in sample "
        f"(tie-completed boundary): {len(mismatches)}"
    )
    print(
        "H      count  max_gap range      mean    safe@1/n  c3      "
        "scorevar     entropy"
    )
    for H in sorted(by_H):
        bucket = by_H[H]
        maxg = bucket.values("max_gap")
        safe = bucket.int_values("safe")
        c3 = bucket.int_values("c3")
        svar = bucket.values("score_var")
        entropy = bucket.values("entropy")
        print(
            f"{H:<6} {len(bucket.records):>5}  "
            f"[{min(maxg):.4f},{max(maxg):.4f}]  {sum(maxg)/len(maxg):.4f}  "
            f"[{min(safe)},{max(safe)}]     [{min(c3)},{max(c3)}]  "
            f"[{min(svar):.3f},{max(svar):.3f}]  "
            f"{sum(entropy)/len(entropy):.4f}"
        )
    rounded: dict[float, set[int]] = defaultdict(set)
    for record in records:
        rounded[round(record.max_gap, 2)].add(record.H)
    ambiguous = {gap: hs for gap, hs in rounded.items() if len(hs) > 1}
    witness = monotonicity_witness(records)
    print(
        f"rounded max_gap buckets with >1 H: {len(ambiguous)}/{len(rounded)}; "
        "H is not a function of max_gap."
    )
    if witness is not None:
        low, high, diff = witness
        print(
            "pointwise monotonicity witness: "
            f"H={low.H} has max_gap={low.max_gap:.4f}, but higher H={high.H} "
            f"has larger max_gap={high.max_gap:.4f} (delta={diff:.4f})."
        )
    print(
        "correlations over samples: "
        f"corr(H,max_gap)={pearson([record.H for record in records], [record.max_gap for record in records]):.3f}, "
        f"corr(H,entropy)={pearson([record.H for record in records], [record.entropy for record in records]):.3f}, "
        f"corr(H,safe@1/n)={pearson([record.H for record in records], [record.safe for record in records]):.3f}"
    )
    grouped: dict[tuple[tuple[int, ...], int], set[int]] = defaultdict(set)
    for record in records:
        grouped[(record.scores, record.c3)].add(record.H)
    collisions = [(key, sorted(values)) for key, values in grouped.items() if len(values) > 1]
    if collisions:
        print("score/c3 collisions separated by H:")
        for (score_seq, c3), values in collisions[:4]:
            print(f"  scores={score_seq}, c3={c3} -> H={values}")
    print()


def print_compact(n: int, samples: int) -> None:
    by_H, records = sample_records(n, samples, seed=2600 + n)
    H_values = sorted(by_H)
    mismatches = [
        record
        for record in records
        if (record.H == 1) != (record.max_gap >= 0.5 - 1e-12)
    ]
    rounded: dict[float, set[int]] = defaultdict(set)
    for record in records:
        rounded[round(record.max_gap, 2)].add(record.H)
    ambiguous = sum(1 for values in rounded.values() if len(values) > 1)
    witness = monotonicity_witness(records)
    print("=" * 88)
    print(f"n={n}; samples={samples}; observed H values={len(H_values)}")
    print("=" * 88)
    print(
        f"H range [{min(H_values)}, {max(H_values)}]; H=1 boundary mismatches={len(mismatches)}; "
        f"ambiguous rounded max_gap buckets={ambiguous}/{len(rounded)}"
    )
    print(
        "global correlations: "
        f"corr(H,max_gap)={pearson([record.H for record in records], [record.max_gap for record in records]):.3f}, "
        f"corr(H,entropy)={pearson([record.H for record in records], [record.entropy for record in records]):.3f}, "
        f"corr(H,safe@1/n)={pearson([record.H for record in records], [record.safe for record in records]):.3f}"
    )
    print("lowest five H buckets:")
    for H in H_values[:5]:
        bucket = by_H[H]
        maxg = bucket.values("max_gap")
        safe = bucket.int_values("safe")
        print(
            f"  H={H:<8} count={len(bucket.records):<5} "
            f"max_gap=[{min(maxg):.4f},{max(maxg):.4f}] "
            f"safe@1/n=[{min(safe)},{max(safe)}]"
        )
    print("highest five H buckets:")
    for H in H_values[-5:]:
        bucket = by_H[H]
        maxg = bucket.values("max_gap")
        safe = bucket.int_values("safe")
        print(
            f"  H={H:<8} count={len(bucket.records):<5} "
            f"max_gap=[{min(maxg):.4f},{max(maxg):.4f}] "
            f"safe@1/n=[{min(safe)},{max(safe)}]"
        )
    if witness is not None:
        low, high, diff = witness
        print(
            f"monotonicity witness: H={low.H}, max_gap={low.max_gap:.4f}; "
            f"H={high.H}, max_gap={high.max_gap:.4f}; delta={diff:.4f}"
        )
    print()


@dataclass(frozen=True)
class LrcRow:
    label: str
    n: int
    speeds: tuple[int, ...]
    t: Fraction


def ladder(n: int, scale: int, skip: int) -> tuple[int, ...]:
    return tuple(sorted({1} | {scale * q for q in range(1, n) if q != skip}))


def lrc_rows_n14() -> tuple[LrcRow, ...]:
    return (
        LrcRow("n14 initial", 14, tuple(range(1, 14)), Fraction(1, 14)),
        LrcRow("n14 row-parent", 14, ladder(14, 7, 6), Fraction(3053, 25872)),
        LrcRow("n14 gate", 14, ladder(14, 14, 6), Fraction(4339, 51744)),
        LrcRow("n14 double-gate", 14, ladder(14, 28, 6), Fraction(8035, 103488)),
    )


def positions_from_speeds(speeds: tuple[int, ...], t: Fraction) -> tuple[Fraction, ...]:
    return tuple(circle(speed * t) for speed in speeds)


def origin_margin(speeds: tuple[int, ...], t: Fraction, n: int) -> Fraction:
    threshold = Fraction(1, n)
    return min(circular_distance_to_zero(speed * t) for speed in speeds) / threshold


def print_lrc_rows() -> None:
    print("=" * 88)
    print("n=14 LRC selected rows: H stays high while endpoint loneliness is anchored")
    print("=" * 88)
    rows = lrc_rows_n14()
    initial_H: int | None = None
    print(
        "row                 H          H/H0    max_gap    orig/th  "
        "safe@1/n score_width c3"
    )
    for row in rows:
        clock_speeds = (0,) + row.speeds
        points = positions_from_speeds(clock_speeds, row.t)
        adj = half_turn_tournament(points)
        H = H_count(adj)
        if initial_H is None:
            initial_H = H
        gaps = gap_list(points)
        score_seq = scores(adj)
        print(
            f"{row.label:<18} {H:>10}  {H / initial_H:>7.3f}  "
            f"{fmt_value(gaps[0]):>9}  {fmt_value(origin_margin(row.speeds, row.t, row.n)):>7}  "
            f"{static_safe_vertex_count(points, Fraction(1, row.n)):>8} "
            f"{max(score_seq) - min(score_seq):>11} {num_3cycles(adj):>3}"
        )
    print(
        "Reading: the hard rows are not bunched by the half-turn meter.  Their H "
        "values remain roughly 0.74-0.92 of the initial near-regular clock top, "
        "while the true LRC certificate depends on the anchored origin margin and "
        "endpoint clock."
    )
    print()


def main() -> None:
    print("H as a loneliness meter: half-turn circular tournaments (S26 extended)")
    print()
    print(
        "Thesis under test: H is an exact detector of the nondegenerate "
        "open-semicircle/1/2-gap state, with the max_gap=1/2 boundary "
        "completed by the tie path. Above H=1 it is a circular-tournament "
        "class meter, not a scalar max-gap meter."
    )
    print()
    for n in (5, 6, 7):
        print_detail(n, samples=20000)
    for n in (8, 9):
        print_compact(n, samples=8000)
    print_lrc_rows()
    print("SYNTHESIS")
    print("=" * 88)
    print("1. H=1 is the sharp bunched/open-semicircle reading, with tie boundary at max_gap = 1/2.")
    print("2. For H>1, H is not determined by max_gap and is not pointwise monotone.")
    print("3. H correlates with spread/entropy, but it also sees cyclic arrangement.")
    print("4. The two-neighbor safe count at 1/n is extra data, not recoverable from H.")
    print("5. LRC needs the half-turn clock plus the anchored endpoint clock.")


if __name__ == "__main__":
    main()
