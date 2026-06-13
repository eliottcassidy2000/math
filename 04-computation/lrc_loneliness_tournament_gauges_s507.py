#!/usr/bin/env python3
"""Candidate LRC tournament gauges and loneliness-meter fingerprints.

codex-2026-06-01 S507

The S26 half-turn tournament makes Hamiltonian-path count a sharp meter for the
1/2 loneliness threshold, but LRC needs the much smaller 1/(k+1) threshold.
This script tries several arc-assignment criteria that keep the tournament
shape while injecting endpoint/slack information from the original LRC.

Tournament Analysis declaration:

* vertices: moving runners in a fixed speed set;
* pairwise observable: phase separation, origin distance, endpoint slack,
  endpoint-wall nearness, side of the origin, or deletion relief;
* switch/gauge: each named gauge below turns the pair observable into one
  directed edge;
* tie Hamiltonian path: increasing speed/index order.

For each gauge, the script compares tournament fingerprints
H, score sequence, score width, score variance, 3-cycles, and SCC shape against
the actual normalized loneliness depth

    (k+1) * min_i ||s_i t||.

Values >= 1 are LRC-safe at that time.
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from functools import lru_cache
from itertools import combinations
from math import sqrt
from statistics import mean


Q = 840
EPS = 1e-12


@dataclass(frozen=True)
class Case:
    name: str
    speeds: tuple[int, ...]


@dataclass(frozen=True)
class Context:
    case: Case
    tick: int
    positions: tuple[int, ...]
    distances: tuple[int, ...]
    signed: tuple[int, ...]
    unsafe: tuple[bool, ...]
    bands: tuple[int, ...]
    wall_nearness: tuple[int, ...]
    relief: tuple[int, ...]
    min_distance: int
    unsafe_count: int

    @property
    def k(self) -> int:
        return len(self.case.speeds)

    @property
    def threshold_den(self) -> int:
        return self.k + 1

    @property
    def lonely_norm(self) -> float:
        return self.min_distance * self.threshold_den / Q

    @property
    def safe(self) -> bool:
        return self.min_distance * self.threshold_den >= Q


@dataclass(frozen=True)
class Gauge:
    name: str
    description: str
    winner: object


@dataclass(frozen=True)
class Record:
    gauge: str
    case: str
    safe: bool
    lonely_norm: float
    unsafe_count: int
    hamiltonian_paths: int
    cyclic_triples: int
    score_sequence: tuple[int, ...]
    score_width: int
    score_variance: float
    scc_count: int
    largest_scc: int


CASES = (
    Case("consecutive-4", (1, 2, 3, 4)),
    Case("consecutive-5", (1, 2, 3, 4, 5)),
    Case("consecutive-6", (1, 2, 3, 4, 5, 6)),
    Case("consecutive-7", (1, 2, 3, 4, 5, 6, 7)),
    Case("prime-6", (2, 3, 5, 7, 11, 13)),
    Case("mixed-7", (1, 4, 6, 9, 10, 15, 21)),
)


def base_winner(_ctx: Context, i: int, j: int) -> int:
    return i if i < j else j


def reverse_base(ctx: Context, i: int, j: int) -> int:
    return j if base_winner(ctx, i, j) == i else i


def clockwise_delta(a: int, b: int) -> int:
    return (b - a) % Q


def circular_distance(a: int, b: int) -> int:
    d = clockwise_delta(a, b)
    return min(d, Q - d)


def half_turn_winner(ctx: Context, i: int, j: int) -> int:
    delta = (ctx.positions[i] - ctx.positions[j]) % Q
    if delta == 0 or 2 * delta == Q:
        return base_winner(ctx, i, j)
    return i if 2 * delta < Q else j


def distance_to_origin(position: int) -> int:
    r = position % Q
    return min(r, Q - r)


def signed_distance(position: int) -> int:
    r = position % Q
    if 2 * r <= Q:
        return r
    return r - Q


def band(distance: int, den: int) -> int:
    """Unsafe, near-safe, or comfortable relative to theta=1/den."""
    if distance * den < Q:
        return 0
    if distance * den < 2 * Q:
        return 1
    return 2


def context_for(case: Case, tick: int) -> Context:
    positions = tuple((speed * tick) % Q for speed in case.speeds)
    distances = tuple(distance_to_origin(position) for position in positions)
    signed = tuple(signed_distance(position) for position in positions)
    den = len(case.speeds) + 1
    unsafe = tuple(distance * den < Q for distance in distances)
    bands = tuple(band(distance, den) for distance in distances)
    wall_nearness = tuple(abs(distance * den - Q) for distance in distances)
    min_distance = min(distances)
    relief = []
    for idx in range(len(case.speeds)):
        without = [d for j, d in enumerate(distances) if j != idx]
        relief.append((min(without) if without else Q // 2) - min_distance)
    return Context(
        case=case,
        tick=tick,
        positions=positions,
        distances=distances,
        signed=signed,
        unsafe=unsafe,
        bands=bands,
        wall_nearness=wall_nearness,
        relief=tuple(relief),
        min_distance=min_distance,
        unsafe_count=sum(unsafe),
    )


def theta_close_switch(ctx: Context, i: int, j: int) -> int:
    dist = circular_distance(ctx.positions[i], ctx.positions[j])
    if dist * ctx.threshold_den <= 2 * Q:
        return reverse_base(ctx, i, j)
    return base_winner(ctx, i, j)


def antipodal_open_switch(ctx: Context, i: int, j: int) -> int:
    dist = circular_distance(ctx.positions[i], ctx.positions[j])
    if abs(2 * dist - Q) * ctx.threshold_den <= 2 * Q:
        return reverse_base(ctx, i, j)
    return base_winner(ctx, i, j)


def unsafe_dominance(ctx: Context, i: int, j: int) -> int:
    if ctx.unsafe[i] != ctx.unsafe[j]:
        return i if ctx.unsafe[i] else j
    return half_turn_winner(ctx, i, j)


def safe_dominance(ctx: Context, i: int, j: int) -> int:
    if ctx.unsafe[i] != ctx.unsafe[j]:
        return i if not ctx.unsafe[i] else j
    return half_turn_winner(ctx, i, j)


def lrc_band_dominance(ctx: Context, i: int, j: int) -> int:
    if ctx.bands[i] != ctx.bands[j]:
        return i if ctx.bands[i] > ctx.bands[j] else j
    return half_turn_winner(ctx, i, j)


def endpoint_wall_pressure(ctx: Context, i: int, j: int) -> int:
    if ctx.wall_nearness[i] != ctx.wall_nearness[j]:
        return i if ctx.wall_nearness[i] < ctx.wall_nearness[j] else j
    return half_turn_winner(ctx, i, j)


def origin_bracket_flip(ctx: Context, i: int, j: int) -> int:
    si = ctx.signed[i]
    sj = ctx.signed[j]
    brackets_origin = si * sj < 0 and (abs(si) + abs(sj)) * ctx.threshold_den <= 2 * Q
    if brackets_origin:
        return j if half_turn_winner(ctx, i, j) == i else i
    return half_turn_winner(ctx, i, j)


def same_side_escape(ctx: Context, i: int, j: int) -> int:
    si = ctx.signed[i]
    sj = ctx.signed[j]
    if si and sj and si * sj > 0 and ctx.distances[i] != ctx.distances[j]:
        return i if ctx.distances[i] > ctx.distances[j] else j
    return half_turn_winner(ctx, i, j)


def relief_bottleneck(ctx: Context, i: int, j: int) -> int:
    if ctx.relief[i] != ctx.relief[j]:
        return i if ctx.relief[i] > ctx.relief[j] else j
    return half_turn_winner(ctx, i, j)


GAUGES = (
    Gauge(
        "half_turn",
        "i beats j when j lies in i's trailing/leading semicircle; the S26 1/2 meter.",
        half_turn_winner,
    ),
    Gauge(
        "theta_close_switch",
        "reverse the base path on pairs closer than 2/(k+1).",
        theta_close_switch,
    ),
    Gauge(
        "antipodal_open_switch",
        "reverse the base path on pairs within 1/(k+1) of antipodal separation.",
        antipodal_open_switch,
    ),
    Gauge(
        "unsafe_dominance",
        "a runner inside the forbidden origin interval beats a safe runner; ties use half-turn.",
        unsafe_dominance,
    ),
    Gauge(
        "safe_dominance",
        "a safe runner beats an unsafe runner; ties use half-turn.",
        safe_dominance,
    ),
    Gauge(
        "lrc_band_dominance",
        "comfortable-safe beats near-safe beats unsafe; equal bands use half-turn.",
        lrc_band_dominance,
    ),
    Gauge(
        "endpoint_wall_pressure",
        "runner nearer the LRC endpoint wall ||st||=1/(k+1) wins; ties use half-turn.",
        endpoint_wall_pressure,
    ),
    Gauge(
        "origin_bracket_flip",
        "flip half-turn on pairs bracketing the origin inside the 2/(k+1) corridor.",
        origin_bracket_flip,
    ),
    Gauge(
        "same_side_escape",
        "on one side of the origin, farther runner wins; opposite-side pairs use half-turn.",
        same_side_escape,
    ),
    Gauge(
        "relief_bottleneck",
        "runner whose deletion most improves the current lonely distance wins.",
        relief_bottleneck,
    ),
)


def tournament_bits(ctx: Context, gauge: Gauge) -> tuple[int, ...]:
    n = ctx.k
    rows = [0] * n
    for i, j in combinations(range(n), 2):
        winner = gauge.winner(ctx, i, j)
        loser = j if winner == i else i
        rows[winner] |= 1 << loser
    return tuple(rows)


def signature(rows: tuple[int, ...]) -> str:
    bits = []
    for i, j in combinations(range(len(rows)), 2):
        bits.append("1" if rows[i] & (1 << j) else "0")
    return "".join(bits)


def score_sequence(rows: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sorted(row.bit_count() for row in rows))


def cyclic_triples(rows: tuple[int, ...]) -> int:
    total = 0
    for a, b, c in combinations(range(len(rows)), 3):
        if (
            (rows[a] & (1 << b) and rows[b] & (1 << c) and rows[c] & (1 << a))
            or (rows[a] & (1 << c) and rows[c] & (1 << b) and rows[b] & (1 << a))
        ):
            total += 1
    return total


def scc_sizes(rows: tuple[int, ...]) -> tuple[int, ...]:
    n = len(rows)
    adj = [[] for _ in range(n)]
    radj = [[] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if rows[i] & (1 << j):
                adj[i].append(j)
                radj[j].append(i)

    seen: set[int] = set()
    order: list[int] = []

    def dfs(v: int) -> None:
        seen.add(v)
        for nxt in adj[v]:
            if nxt not in seen:
                dfs(nxt)
        order.append(v)

    for v in range(n):
        if v not in seen:
            dfs(v)

    seen.clear()
    sizes: list[int] = []
    for start in reversed(order):
        if start in seen:
            continue
        todo = deque([start])
        seen.add(start)
        size = 0
        while todo:
            v = todo.pop()
            size += 1
            for nxt in radj[v]:
                if nxt not in seen:
                    seen.add(nxt)
                    todo.append(nxt)
        sizes.append(size)
    return tuple(sorted(sizes, reverse=True))


@lru_cache(maxsize=None)
def hamiltonian_paths(rows: tuple[int, ...]) -> int:
    n = len(rows)
    full = (1 << n) - 1

    @lru_cache(maxsize=None)
    def dp(mask: int, last: int) -> int:
        if mask == full:
            return 1
        total = 0
        available = (~mask) & full & rows[last]
        while available:
            bit = available & -available
            nxt = bit.bit_length() - 1
            total += dp(mask | bit, nxt)
            available -= bit
        return total

    return sum(dp(1 << start, start) for start in range(n))


def score_variance(seq: tuple[int, ...]) -> float:
    center = (len(seq) - 1) / 2
    return sum((score - center) ** 2 for score in seq) / len(seq)


def make_record(ctx: Context, gauge: Gauge) -> Record:
    rows = tournament_bits(ctx, gauge)
    seq = score_sequence(rows)
    sizes = scc_sizes(rows)
    return Record(
        gauge=gauge.name,
        case=ctx.case.name,
        safe=ctx.safe,
        lonely_norm=ctx.lonely_norm,
        unsafe_count=ctx.unsafe_count,
        hamiltonian_paths=hamiltonian_paths(rows),
        cyclic_triples=cyclic_triples(rows),
        score_sequence=seq,
        score_width=max(seq) - min(seq),
        score_variance=score_variance(seq),
        scc_count=len(sizes),
        largest_scc=sizes[0],
    )


def rank_values(values: list[float]) -> list[float]:
    indexed = sorted(enumerate(values), key=lambda item: item[1])
    ranks = [0.0] * len(values)
    i = 0
    while i < len(indexed):
        j = i + 1
        while j < len(indexed) and indexed[j][1] == indexed[i][1]:
            j += 1
        rank = (i + j - 1) / 2
        for pos in range(i, j):
            ranks[indexed[pos][0]] = rank
        i = j
    return ranks


def pearson(xs: list[float], ys: list[float]) -> float:
    mx = mean(xs)
    my = mean(ys)
    vx = sum((x - mx) ** 2 for x in xs)
    vy = sum((y - my) ** 2 for y in ys)
    if vx <= EPS or vy <= EPS:
        return 0.0
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sqrt(vx * vy)


def spearman(xs: list[float], ys: list[float]) -> float:
    return pearson(rank_values(xs), rank_values(ys))


def bucket_stats(
    records: list[Record], key
) -> tuple[float, float, float, int, int, float, float]:
    buckets: dict[object, list[Record]] = defaultdict(list)
    for record in records:
        buckets[key(record)].append(record)

    pure_weight = 0
    pure_safe_weight = 0
    pure_unsafe_weight = 0
    weighted_span = 0.0
    max_span = 0.0
    ambiguous = 0
    total_safe = sum(record.safe for record in records)
    total_unsafe = len(records) - total_safe
    for rows in buckets.values():
        safe_count = sum(1 for record in rows if record.safe)
        unsafe_count = len(rows) - safe_count
        pure_weight += max(safe_count, unsafe_count)
        if safe_count and unsafe_count:
            ambiguous += 1
        elif safe_count:
            pure_safe_weight += safe_count
        else:
            pure_unsafe_weight += unsafe_count
        values = [record.lonely_norm for record in rows]
        span = max(values) - min(values)
        weighted_span += span * len(rows)
        max_span = max(max_span, span)
    return (
        pure_weight / len(records),
        pure_safe_weight / total_safe if total_safe else 1.0,
        pure_unsafe_weight / total_unsafe if total_unsafe else 1.0,
        ambiguous,
        len(buckets),
        weighted_span / len(records),
        max_span,
    )


def feature_values(records: list[Record], feature: str) -> list[float]:
    return [float(getattr(record, feature)) for record in records]


def summarize_gauge(records: list[Record]) -> dict[str, object]:
    numeric_features = (
        "hamiltonian_paths",
        "cyclic_triples",
        "score_width",
        "score_variance",
        "scc_count",
        "largest_scc",
    )
    y = [record.lonely_norm for record in records]
    correlations = {
        feature: spearman(feature_values(records, feature), y)
        for feature in numeric_features
    }
    best_feature = max(correlations, key=lambda item: abs(correlations[item]))
    h_purity = bucket_stats(records, lambda record: record.hamiltonian_paths)
    score_purity = bucket_stats(records, lambda record: record.score_sequence)
    hs_purity = bucket_stats(
        records, lambda record: (record.hamiltonian_paths, record.score_sequence)
    )
    usefulness = (
        0.35 * ((hs_purity[1] + hs_purity[2]) / 2)
        + 0.25 * hs_purity[1]
        + 0.20 * abs(correlations[best_feature])
        + 0.20 * (1.0 - min(1.0, hs_purity[5]))
    )
    return {
        "best_feature": best_feature,
        "best_rho": correlations[best_feature],
        "correlations": correlations,
        "h_purity": h_purity,
        "score_purity": score_purity,
        "hs_purity": hs_purity,
        "usefulness": usefulness,
    }


def fmt_float(value: float) -> str:
    return f"{value: .3f}"


def top_distributions(records: list[Record], attr: str, safe: bool) -> str:
    counter = Counter(getattr(record, attr) for record in records if record.safe == safe)
    total = sum(counter.values())
    if not total:
        return "-"
    parts = []
    for value, count in counter.most_common(4):
        parts.append(f"{value}:{count / total:.2f}")
    return ", ".join(parts)


def case_summary(case: Case) -> str:
    contexts = [context_for(case, tick) for tick in range(1, Q)]
    safe_count = sum(ctx.safe for ctx in contexts)
    best = max(ctx.lonely_norm for ctx in contexts)
    mean_lonely = mean(ctx.lonely_norm for ctx in contexts)
    unsafe_hist = Counter(ctx.unsafe_count for ctx in contexts)
    hist = " ".join(f"{count}:{unsafe_hist[count]}" for count in sorted(unsafe_hist))
    return (
        f"{case.name:<15} k={len(case.speeds)} theta=1/{len(case.speeds)+1:<2} "
        f"safe={safe_count:>3}/{Q-1:<3} max_norm={best:5.3f} "
        f"mean_norm={mean_lonely:5.3f} unsafe_count_hist={hist}"
    )


def main() -> None:
    print("LRC loneliness tournament gauges (codex-2026-06-01-S507)")
    print(f"sample grid: t=a/{Q}, a=1..{Q - 1}")
    print()
    print("Tournament Analysis declaration")
    print("  vertices: moving runners in the speed set")
    print("  target: normalized loneliness depth (k+1)*min_i ||s_i t||")
    print("  tie Hamiltonian path: increasing speed/index order")
    print("  gauges:")
    for gauge in GAUGES:
        print(f"    {gauge.name:<24} {gauge.description}")
    print()

    print("Case target summaries")
    for case in CASES:
        print("  " + case_summary(case))
    print()

    records_by_gauge: dict[str, list[Record]] = {gauge.name: [] for gauge in GAUGES}
    by_case_gauge: dict[tuple[str, str], list[Record]] = defaultdict(list)
    for case in CASES:
        for tick in range(1, Q):
            ctx = context_for(case, tick)
            for gauge in GAUGES:
                record = make_record(ctx, gauge)
                records_by_gauge[gauge.name].append(record)
                by_case_gauge[(case.name, gauge.name)].append(record)

    summaries = {
        gauge_name: summarize_gauge(records)
        for gauge_name, records in records_by_gauge.items()
    }
    ranked = sorted(
        summaries.items(), key=lambda item: item[1]["usefulness"], reverse=True
    )

    print("Gauge leaderboard across all cases")
    header = (
        "  gauge                    use   best rho(feature)       "
        "pure safeIso unsafeIso  span(H+score)"
    )
    print(header)
    for gauge_name, summary in ranked:
        h_purity = summary["h_purity"]
        score_purity = summary["score_purity"]
        hs_purity = summary["hs_purity"]
        print(
            f"  {gauge_name:<24}"
            f"{summary['usefulness']:5.3f}  "
            f"{summary['best_rho']:>6.3f}({summary['best_feature']:<17})  "
            f"{hs_purity[0]:5.3f}  {hs_purity[1]:7.3f}  {hs_purity[2]:9.3f}  "
            f"{hs_purity[5]:6.3f}"
        )
    print()

    print("Per-case best gauges by H+score safe/unsafe purity")
    for case in CASES:
        rows = []
        for gauge in GAUGES:
            stats = bucket_stats(
                by_case_gauge[(case.name, gauge.name)],
                lambda record: (record.hamiltonian_paths, record.score_sequence),
            )
            balanced = (stats[1] + stats[2]) / 2
            rows.append((balanced, stats[1], -stats[5], gauge.name, stats))
        rows.sort(reverse=True)
        best = rows[:4]
        parts = []
        for _balanced, _safe_iso, _neg_span, name, stats in best:
            parts.append(
                f"{name} safeIso={stats[1]:.3f} amb={stats[3]}/{stats[4]} span={stats[5]:.3f}"
            )
        print(f"  {case.name:<15} " + " | ".join(parts))
    print()

    print("Safe vs unsafe fingerprint sketches")
    for gauge_name in [name for name, _summary in ranked[:5]]:
        records = records_by_gauge[gauge_name]
        print(f"  {gauge_name}")
        print(
            f"    H unsafe: {top_distributions(records, 'hamiltonian_paths', False)}"
        )
        print(f"    H safe:   {top_distributions(records, 'hamiltonian_paths', True)}")
        print(
            f"    score unsafe: {top_distributions(records, 'score_sequence', False)}"
        )
        print(
            f"    score safe:   {top_distributions(records, 'score_sequence', True)}"
        )
    print()

    print("Interpretation")
    print("  1. Pure half-turn remains a spread/1/2-gap meter, not an LRC-safe meter.")
    print("  2. Unsafe/safe dominance turns the score sequence into a coarse LRC")
    print("     witness detector: zero unsafe runners becomes a visible block event.")
    print("  3. Theta-close and antipodal-open switches are the best shape-only")
    print("     candidates: weaker than status dominance, but less tautological and")
    print("     especially good on the mixed-7 case.")
    print("  4. LRC band dominance records unsafe/near/deep endpoint slack, but by")
    print("     itself is a pressure diagnostic rather than a safe-witness classifier.")
    print("  5. Endpoint-wall pressure is a boundary meter rather than a witness meter:")
    print("     it spikes near ||s_i t|| = 1/(k+1), useful for corridor/peel scripts.")
    print("  6. Origin-bracket and same-side escape are shape gauges: weaker as direct")
    print("     classifiers, but they preserve pair geometry that scalar rank gauges")
    print("     erase.")


if __name__ == "__main__":
    main()
