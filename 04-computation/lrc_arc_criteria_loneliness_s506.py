#!/usr/bin/env python3
"""
lrc_arc_criteria_loneliness_s506.py

codex-2026-06-01 S506

Search for LRC tournament-analysis arc rules whose tournament shape gives a
useful loneliness-related metric.

Each criterion below is a Tournament Analysis switchboard:

* objects: runners, either all runners including the observer 0 or only moving
  runners;
* pairwise observable: phase, local safe-gap scores, deletion relief, or
  threshold deficit;
* switch/gauge: which side of a pair wins the directed arc;
* tie path: increasing runner label, i -> j for i < j.

The script tests two uses.

1. Small exact clock cells: rank every arc criterion against LRC targets
   (origin clearance, origin two-sided bracket, safe-gap count, local lonely
   vertices, max circular gap).

2. n=14/n=18 hard rows: record score histograms, directed triangles, SCCs,
   edge-local tie rates, and Hamiltonian path counts when feasible.
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from importlib.machinery import SourceFileLoader
from itertools import combinations
from math import comb, sqrt
from pathlib import Path
from statistics import mean


ROOT = Path(__file__).resolve().parents[1]
S356 = SourceFileLoader(
    "lonely_runner_residue_probe_s356",
    str(ROOT / "04-computation" / "lonely_runner_residue_probe_s356.py"),
).load_module()
S360 = SourceFileLoader(
    "lonely_runner_endpoint_protection_s360",
    str(ROOT / "04-computation" / "lonely_runner_endpoint_protection_s360.py"),
).load_module()
S492 = SourceFileLoader(
    "lrc_n14_n18_tournament_pingpong_s492",
    str(ROOT / "04-computation" / "lrc_n14_n18_tournament_pingpong_s492.py"),
).load_module()
clock = SourceFileLoader(
    "tournament_clock_s24",
    str(ROOT / "04-computation" / "tournament_clock_s24.py"),
).load_module()


ONE = Fraction(1, 1)
HALF = Fraction(1, 2)
H_LIMIT = 14


@dataclass(frozen=True)
class Geometry:
    threshold: Fraction
    stationary_min: Fraction
    origin_left: Fraction
    origin_right: Fraction
    max_gap: Fraction
    min_gap: Fraction
    safe_gap_count: int
    lonely_vertices: tuple[int, ...]


@dataclass(frozen=True)
class Criterion:
    name: str
    scope: str
    observable: str
    switch: str
    tie_path: str


@dataclass(frozen=True)
class Snapshot:
    criterion: Criterion
    labels: tuple[int, ...]
    out: tuple[int, ...]
    strict_out: tuple[int, ...]
    ties: int
    strict_arcs: int


@dataclass(frozen=True)
class Fingerprint:
    vertex_count: int
    hamiltonian_paths: int | None
    score_hist: tuple[tuple[int, int], ...]
    score_width: int
    origin_score: int | None
    strict_triangles: int
    complete_triangles: int
    largest_strict_scc: int
    ties: int
    strict_arcs: int


CRITERIA: tuple[Criterion, ...] = (
    Criterion(
        "phase_half",
        "all",
        "pair phase delta frac(x_j-x_i)",
        "i -> j iff j is clockwise from i by less than 1/2",
        "collisions/antipodes use increasing labels",
    ),
    Criterion(
        "lrc_close_sector",
        "all",
        "pair phase delta frac(x_j-x_i)",
        "strict arc only inside the LRC danger sector of width 1/n",
        "non-close pairs use increasing labels",
    ),
    Criterion(
        "local_moat_min",
        "all",
        "vertex score min(left_gap,right_gap)",
        "larger two-sided local moat points to smaller local moat",
        "equal moats use increasing labels",
    ),
    Criterion(
        "local_moat_sum",
        "all",
        "vertex score left_gap+right_gap",
        "larger adjacent gap sum points to smaller sum",
        "equal sums use increasing labels",
    ),
    Criterion(
        "safe_deficit",
        "all",
        "vertex score local deficit below threshold on two adjacent gaps",
        "smaller deficit points to larger deficit",
        "equal deficits use increasing labels",
    ),
    Criterion(
        "origin_distance",
        "moving",
        "moving runner distance from stationary runner 0",
        "farther from 0 points to closer to 0",
        "equal distances use increasing labels",
    ),
    Criterion(
        "origin_blocker_relief",
        "moving",
        "origin nearest distance after deleting one moving runner",
        "runner whose deletion most improves origin moat points to the other",
        "equal relief uses increasing labels",
    ),
    Criterion(
        "origin_deficit_relief",
        "moving",
        "origin threshold-deficit relief after deleting one moving runner",
        "runner whose deletion removes more origin deficit points to the other",
        "equal relief uses increasing labels",
    ),
    Criterion(
        "k1_relief_pressure",
        "all",
        "pairwise nearest-neighbor deletion relief",
        "blocker points to runner whose k=1 nearest moat it improves more",
        "equal relief uses increasing labels",
    ),
    Criterion(
        "k2_relief_pressure",
        "all",
        "pairwise two-nearest-neighbor deletion relief",
        "blocker points to runner whose k=2 nearest moat it improves more",
        "equal relief uses increasing labels",
    ),
    Criterion(
        "threshold_deficit_pressure",
        "all",
        "pairwise relief of two-neighbor deficit below 1/n",
        "blocker points to runner whose threshold deficit it removes more",
        "equal relief uses increasing labels",
    ),
    Criterion(
        "bracket_relief_pressure",
        "all",
        "pairwise relief of circular two-sided local bracket",
        "blocker points to runner whose adjacent safe bracket it improves more",
        "equal relief uses increasing labels",
    ),
)


def fmt(value: Fraction | int | None) -> str:
    if value is None:
        return "-"
    if isinstance(value, int):
        return str(value)
    return S356.fmt_frac(value)


def mod1(value: Fraction) -> Fraction:
    return value % ONE


def circular_distance(a: Fraction, b: Fraction) -> Fraction:
    delta = mod1(b - a)
    return min(delta, ONE - delta)


def clockwise_delta(a: Fraction, b: Fraction) -> Fraction:
    return mod1(b - a)


def positions(speeds_with_observer: tuple[int, ...], t: Fraction) -> tuple[Fraction, ...]:
    return tuple(mod1(Fraction(speed) * t) for speed in speeds_with_observer)


def local_gaps(
    pos: tuple[Fraction, ...], labels: tuple[int, ...] | None = None
) -> tuple[list[int], list[Fraction]]:
    if labels is None:
        labels = tuple(range(len(pos)))
    ordered = [label for label, _ in sorted(((label, pos[label]) for label in labels), key=lambda x: (x[1], x[0]))]
    gaps: list[Fraction] = []
    for idx, label in enumerate(ordered):
        nxt = ordered[(idx + 1) % len(ordered)]
        gaps.append(mod1(pos[nxt] - pos[label]))
    return ordered, gaps


def bracket_values(
    pos: tuple[Fraction, ...], labels: tuple[int, ...] | None = None
) -> tuple[dict[int, Fraction], dict[int, Fraction], dict[int, tuple[Fraction, Fraction]]]:
    ordered, gaps = local_gaps(pos, labels)
    min_bracket: dict[int, Fraction] = {}
    sum_bracket: dict[int, Fraction] = {}
    pair_bracket: dict[int, tuple[Fraction, Fraction]] = {}
    for idx, label in enumerate(ordered):
        left = gaps[idx - 1]
        right = gaps[idx]
        min_bracket[label] = min(left, right)
        sum_bracket[label] = left + right
        pair_bracket[label] = (left, right)
    return min_bracket, sum_bracket, pair_bracket


def geometry(speeds_with_observer: tuple[int, ...], t: Fraction) -> Geometry:
    pos = positions(speeds_with_observer, t)
    n = len(pos)
    threshold = Fraction(1, n)
    ordered, gaps = local_gaps(pos)
    order_index = {label: idx for idx, label in enumerate(ordered)}
    safe = tuple(gap >= threshold for gap in gaps)
    lonely = tuple(
        sorted(label for label in ordered if safe[order_index[label] - 1] and safe[order_index[label]])
    )
    origin_idx = order_index[0]
    stationary_min = min(circular_distance(pos[0], pos[label]) for label in range(1, n))
    return Geometry(
        threshold=threshold,
        stationary_min=stationary_min,
        origin_left=gaps[origin_idx - 1],
        origin_right=gaps[origin_idx],
        max_gap=max(gaps),
        min_gap=min(gaps),
        safe_gap_count=sum(1 for value in safe if value),
        lonely_vertices=lonely,
    )


def distance_matrix(pos: tuple[Fraction, ...]) -> tuple[tuple[Fraction, ...], ...]:
    return tuple(
        tuple(Fraction(0) if i == j else circular_distance(pos[i], pos[j]) for j in range(len(pos)))
        for i in range(len(pos))
    )


def nearest_sum(
    dist: tuple[tuple[Fraction, ...], ...], i: int, k: int, removed: int | None = None
) -> Fraction:
    values = sorted(dist[i][j] for j in range(len(dist)) if j != i and j != removed)
    return sum(values[:k], Fraction(0))


def first_distances(
    dist: tuple[tuple[Fraction, ...], ...], i: int, k: int, removed: int | None = None
) -> tuple[Fraction, ...]:
    return tuple(sorted(dist[i][j] for j in range(len(dist)) if j != i and j != removed)[:k])


def deficit_score(values: tuple[Fraction, ...], threshold: Fraction) -> Fraction:
    return sum((threshold - value for value in values if value < threshold), Fraction(0))


def bracket_after_delete(
    pos: tuple[Fraction, ...], label: int, removed: int
) -> Fraction:
    labels = tuple(v for v in range(len(pos)) if v != removed)
    min_bracket, _, _ = bracket_values(pos, labels)
    return min_bracket[label]


def compare_fraction(left: Fraction, right: Fraction) -> int:
    if left > right:
        return 1
    if right > left:
        return -1
    return 0


def criterion_labels(criterion: Criterion, speeds_with_observer: tuple[int, ...]) -> tuple[int, ...]:
    if criterion.scope == "all":
        return tuple(range(len(speeds_with_observer)))
    if criterion.scope == "moving":
        return tuple(range(1, len(speeds_with_observer)))
    raise ValueError(f"unknown criterion scope {criterion.scope}")


def compare_for_criterion(
    criterion: Criterion,
    i: int,
    j: int,
    speeds_with_observer: tuple[int, ...],
    t: Fraction,
) -> int:
    pos = positions(speeds_with_observer, t)
    n = len(pos)
    threshold = Fraction(1, n)
    dist = distance_matrix(pos)
    min_bracket, sum_bracket, pair_bracket = bracket_values(pos)

    if criterion.name == "phase_half":
        delta = clockwise_delta(pos[i], pos[j])
        if delta == 0 or delta == HALF:
            return 0
        return 1 if delta < HALF else -1

    if criterion.name == "lrc_close_sector":
        delta = clockwise_delta(pos[i], pos[j])
        if delta == 0:
            return 0
        if delta < threshold:
            return 1
        if ONE - delta < threshold:
            return -1
        return 0

    if criterion.name == "local_moat_min":
        return compare_fraction(min_bracket[i], min_bracket[j])

    if criterion.name == "local_moat_sum":
        return compare_fraction(sum_bracket[i], sum_bracket[j])

    if criterion.name == "safe_deficit":
        left_i, right_i = pair_bracket[i]
        left_j, right_j = pair_bracket[j]
        deficit_i = deficit_score((left_i, right_i), threshold)
        deficit_j = deficit_score((left_j, right_j), threshold)
        return compare_fraction(deficit_j, deficit_i)

    if criterion.name == "origin_distance":
        return compare_fraction(dist[0][i], dist[0][j])

    if criterion.name == "origin_blocker_relief":
        moving = tuple(label for label in range(1, n))
        base = min(dist[0][label] for label in moving)
        after_i = min(dist[0][label] for label in moving if label != i)
        after_j = min(dist[0][label] for label in moving if label != j)
        return compare_fraction(after_i - base, after_j - base)

    if criterion.name == "origin_deficit_relief":
        moving = tuple(label for label in range(1, n))
        base = sum((threshold - dist[0][label] for label in moving if dist[0][label] < threshold), Fraction(0))
        after_i = sum(
            (threshold - dist[0][label] for label in moving if label != i and dist[0][label] < threshold),
            Fraction(0),
        )
        after_j = sum(
            (threshold - dist[0][label] for label in moving if label != j and dist[0][label] < threshold),
            Fraction(0),
        )
        return compare_fraction(base - after_i, base - after_j)

    if criterion.name == "k1_relief_pressure":
        base_i = nearest_sum(dist, i, 1)
        base_j = nearest_sum(dist, j, 1)
        relief_i_from_j = nearest_sum(dist, i, 1, removed=j) - base_i
        relief_j_from_i = nearest_sum(dist, j, 1, removed=i) - base_j
        return compare_fraction(relief_j_from_i, relief_i_from_j)

    if criterion.name == "k2_relief_pressure":
        base_i = nearest_sum(dist, i, 2)
        base_j = nearest_sum(dist, j, 2)
        relief_i_from_j = nearest_sum(dist, i, 2, removed=j) - base_i
        relief_j_from_i = nearest_sum(dist, j, 2, removed=i) - base_j
        return compare_fraction(relief_j_from_i, relief_i_from_j)

    if criterion.name == "threshold_deficit_pressure":
        base_i = deficit_score(first_distances(dist, i, 2), threshold)
        base_j = deficit_score(first_distances(dist, j, 2), threshold)
        after_i = deficit_score(first_distances(dist, i, 2, removed=j), threshold)
        after_j = deficit_score(first_distances(dist, j, 2, removed=i), threshold)
        relief_i_from_j = base_i - after_i
        relief_j_from_i = base_j - after_j
        return compare_fraction(relief_j_from_i, relief_i_from_j)

    if criterion.name == "bracket_relief_pressure":
        relief_i_from_j = bracket_after_delete(pos, i, j) - min_bracket[i]
        relief_j_from_i = bracket_after_delete(pos, j, i) - min_bracket[j]
        return compare_fraction(relief_j_from_i, relief_i_from_j)

    raise ValueError(f"missing comparator for {criterion.name}")


def snapshot(
    criterion: Criterion, speeds_with_observer: tuple[int, ...], t: Fraction
) -> Snapshot:
    labels = criterion_labels(criterion, speeds_with_observer)
    label_to_local = {label: idx for idx, label in enumerate(labels)}
    out = [0] * len(labels)
    strict_out = [0] * len(labels)
    ties = 0
    strict_arcs = 0
    for a, b in combinations(labels, 2):
        cmp = compare_for_criterion(criterion, a, b, speeds_with_observer, t)
        if cmp > 0:
            winner, loser = a, b
            strict_arcs += 1
        elif cmp < 0:
            winner, loser = b, a
            strict_arcs += 1
        else:
            winner, loser = a, b
            ties += 1
        wi = label_to_local[winner]
        li = label_to_local[loser]
        out[wi] |= 1 << li
        if cmp != 0:
            strict_out[wi] |= 1 << li
    return Snapshot(
        criterion=criterion,
        labels=labels,
        out=tuple(out),
        strict_out=tuple(strict_out),
        ties=ties,
        strict_arcs=strict_arcs,
    )


@lru_cache(maxsize=None)
def h_count(out: tuple[int, ...]) -> int:
    n = len(out)
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        unvisited = full ^ mask
        for last in range(n):
            val = dp[mask][last]
            if not val:
                continue
            nxts = out[last] & unvisited
            while nxts:
                bit = nxts & -nxts
                nxt = bit.bit_length() - 1
                dp[mask | bit][nxt] += val
                nxts ^= bit
    return sum(dp[full])


def score_list(out: tuple[int, ...]) -> list[int]:
    return [row.bit_count() for row in out]


def triangle_count(out: tuple[int, ...]) -> int:
    total = 0
    for a, b, c in combinations(range(len(out)), 3):
        if ((out[a] >> b) & 1 and (out[b] >> c) & 1 and (out[c] >> a) & 1) or (
            (out[a] >> c) & 1 and (out[c] >> b) & 1 and (out[b] >> a) & 1
        ):
            total += 1
    return total


def scc_sizes(out: tuple[int, ...]) -> tuple[int, ...]:
    n = len(out)
    graph = [[] for _ in range(n)]
    reverse = [[] for _ in range(n)]
    for i in range(n):
        bits = out[i]
        while bits:
            bit = bits & -bits
            j = bit.bit_length() - 1
            graph[i].append(j)
            reverse[j].append(i)
            bits ^= bit

    seen: set[int] = set()
    order: list[int] = []

    def dfs(v: int) -> None:
        seen.add(v)
        for nxt in graph[v]:
            if nxt not in seen:
                dfs(nxt)
        order.append(v)

    for v in range(n):
        if v not in seen:
            dfs(v)

    sizes: list[int] = []
    seen.clear()
    for start in reversed(order):
        if start in seen:
            continue
        todo = deque([start])
        seen.add(start)
        size = 0
        while todo:
            v = todo.pop()
            size += 1
            for nxt in reverse[v]:
                if nxt not in seen:
                    seen.add(nxt)
                    todo.append(nxt)
        sizes.append(size)
    return tuple(sorted(sizes, reverse=True))


def edge_key(out: tuple[int, ...]) -> tuple[int, ...]:
    bits: list[int] = []
    for i, j in combinations(range(len(out)), 2):
        bits.append(1 if (out[i] >> j) & 1 else 0)
    return tuple(bits)


def edge_flips(left: tuple[int, ...], right: tuple[int, ...]) -> int:
    return sum(1 for a, b in zip(left, right) if a != b)


def fingerprint(snap: Snapshot, compute_h: bool) -> Fingerprint:
    scores = score_list(snap.out)
    origin_score = None
    if 0 in snap.labels:
        origin_score = scores[snap.labels.index(0)]
    strict_scc = scc_sizes(snap.strict_out)
    return Fingerprint(
        vertex_count=len(snap.labels),
        hamiltonian_paths=h_count(snap.out) if compute_h and len(snap.labels) <= H_LIMIT else None,
        score_hist=tuple(sorted(Counter(scores).items())),
        score_width=max(scores) - min(scores) if scores else 0,
        origin_score=origin_score,
        strict_triangles=triangle_count(snap.strict_out),
        complete_triangles=triangle_count(snap.out),
        largest_strict_scc=strict_scc[0] if strict_scc else 0,
        ties=snap.ties,
        strict_arcs=snap.strict_arcs,
    )


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
    order = sorted(range(len(values)), key=lambda idx: values[idx])
    out = [0.0] * len(values)
    cursor = 0
    while cursor < len(order):
        end = cursor + 1
        while end < len(order) and values[order[end]] == values[order[cursor]]:
            end += 1
        rank = (cursor + end - 1) / 2
        for item in range(cursor, end):
            out[order[item]] = rank
        cursor = end
    return out


def spearman(xs: list[float], ys: list[float]) -> float:
    return pearson(ranks(xs), ranks(ys))


def exact_cell_times(speeds_with_observer: tuple[int, ...]) -> tuple[Fraction, ...]:
    return tuple(t for t, _ in clock.clock_cells(speeds_with_observer))


def feature_series(snapshots: list[Snapshot], geoms: list[Geometry]) -> dict[str, list[float]]:
    out: dict[str, list[float]] = defaultdict(list)
    for snap, geom in zip(snapshots, geoms):
        fp = fingerprint(snap, compute_h=True)
        denom = max(1, fp.vertex_count - 1)
        total_pairs = max(1, comb(fp.vertex_count, 2))
        total_triples = max(1, comb(fp.vertex_count, 3)) if fp.vertex_count >= 3 else 1
        out["H"].append(float(fp.hamiltonian_paths or 0))
        out["score_width_norm"].append(float(fp.score_width / denom))
        out["tie_rate"].append(float(fp.ties / total_pairs))
        out["strict_arc_rate"].append(float(fp.strict_arcs / total_pairs))
        out["strict_triangle_density"].append(float(fp.strict_triangles / total_triples))
        out["complete_triangle_density"].append(float(fp.complete_triangles / total_triples))
        out["strict_scc_norm"].append(float(fp.largest_strict_scc / fp.vertex_count))
        if fp.origin_score is not None:
            out["origin_score_norm"].append(float(fp.origin_score / denom))
        else:
            out["origin_score_norm"].append(0.0)
        out["top_score_norm"].append(float(max(score_list(snap.out)) / denom))
    return dict(out)


def target_series(geoms: list[Geometry]) -> dict[str, list[float]]:
    return {
        "origin_clearance": [float(g.stationary_min / g.threshold) for g in geoms],
        "origin_bracket": [float(min(g.origin_left, g.origin_right) / g.threshold) for g in geoms],
        "safe_gap_count": [float(g.safe_gap_count) for g in geoms],
        "lonely_vertices": [float(len(g.lonely_vertices)) for g in geoms],
        "max_gap": [float(g.max_gap) for g in geoms],
    }


def small_families() -> tuple[tuple[str, tuple[int, ...]], ...]:
    return (
        ("initial n=5", (0, 1, 2, 3, 4)),
        ("primes n=5", (0, 2, 3, 5, 7)),
        ("spread n=5", (0, 1, 4, 9, 11)),
        ("initial n=6", (0, 1, 2, 3, 4, 5)),
        ("primes n=6", (0, 2, 3, 5, 7, 11)),
        ("initial n=7", (0, 1, 2, 3, 4, 5, 6)),
        ("spread n=7", (0, 1, 3, 7, 12, 18, 27)),
    )


def criterion_scorecard() -> dict[str, dict[str, object]]:
    summary: dict[str, dict[str, object]] = {}
    for criterion in CRITERIA:
        best_scores: list[tuple[float, str, str]] = []
        h_maxgap: list[float] = []
        h_lonely: list[float] = []
        origin_score_origin_bracket: list[float] = []
        flip_rates: list[float] = []
        tie_rates: list[float] = []

        for _, speeds in small_families():
            times = exact_cell_times(speeds)
            geoms = [geometry(speeds, t) for t in times]
            snaps = [snapshot(criterion, speeds, t) for t in times]
            features = feature_series(snaps, geoms)
            targets = target_series(geoms)
            local_best = (0.0, "-", "-")
            for feature_name, xs in features.items():
                for target_name, ys in targets.items():
                    rho = spearman(xs, ys)
                    if abs(rho) > abs(local_best[0]):
                        local_best = (rho, feature_name, target_name)
            best_scores.append(local_best)
            h_maxgap.append(spearman(features["H"], targets["max_gap"]))
            h_lonely.append(spearman(features["H"], targets["lonely_vertices"]))
            origin_score_origin_bracket.append(
                spearman(features["origin_score_norm"], targets["origin_bracket"])
            )
            keys = [edge_key(snap.out) for snap in snaps]
            flips = [edge_flips(a, b) for a, b in zip(keys, keys[1:])]
            pair_count = max(1, comb(len(snaps[0].labels), 2))
            if flips:
                flip_rates.append(mean(flips) / pair_count)
            total_pairs = pair_count * len(snaps)
            tie_rates.append(sum(snap.ties for snap in snaps) / total_pairs)

        best_abs = mean(abs(row[0]) for row in best_scores)
        feature_votes = Counter((row[1], row[2]) for row in best_scores)
        best_pair, _ = feature_votes.most_common(1)[0]
        summary[criterion.name] = {
            "criterion": criterion,
            "best_abs": best_abs,
            "best_pair": best_pair,
            "h_maxgap": mean(h_maxgap),
            "h_lonely": mean(h_lonely),
            "origin_score_origin_bracket": mean(origin_score_origin_bracket),
            "flip_rate": mean(flip_rates),
            "tie_rate": mean(tie_rates),
        }
    return summary


def print_methodology() -> None:
    print("LRC arc criteria for tournament-analysis loneliness metrics")
    print("=" * 100)
    print("Tie Hamiltonian path for every criterion: increasing runner label i -> j for i < j.")
    print("Strict arcs are reported separately from tie-path completion.")
    print()
    print("Criteria switchboard:")
    for item in CRITERIA:
        print(f"  {item.name:<28} vertices={item.scope:<6} observable={item.observable}")
        print(f"    switch={item.switch}")
    print()


def print_scorecard(summary: dict[str, dict[str, object]]) -> None:
    print("SMALL EXACT CLOCK SCORECARD")
    print("=" * 100)
    print(
        "criterion                    best_abs  common_best_feature->target"
        "                  H~maxgap H~lonely originScore~bracket flips ties"
    )
    for name, row in sorted(summary.items(), key=lambda item: -float(item[1]["best_abs"])):
        best_feature, best_target = row["best_pair"]  # type: ignore[misc]
        print(
            f"{name:<28} "
            f"{float(row['best_abs']):>8.3f}  "
            f"{best_feature}->{best_target:<30} "
            f"{float(row['h_maxgap']):>+8.3f} "
            f"{float(row['h_lonely']):>+8.3f} "
            f"{float(row['origin_score_origin_bracket']):>+10.3f} "
            f"{float(row['flip_rate']):>5.3f} "
            f"{float(row['tie_rate']):>5.3f}"
        )
    print()
    print("Reading the scorecard:")
    print("  best_abs is the mean absolute within-family Spearman score for the best")
    print("  fingerprint/target pairing.  H~maxgap and H~lonely are shown separately")
    print("  because H is only useful when the switch itself preserves cyclic structure.")
    print()


def initial(n: int) -> tuple[int, ...]:
    return tuple(range(1, n))


def selected_examples() -> tuple[tuple[str, int, tuple[int, ...]], ...]:
    n14_lpd_skip, n14_lpd = S492.best_ladder(14, S492.largest_proper_divisor(14))
    n18_lpd_skip, n18_lpd = S492.best_ladder(18, S492.largest_proper_divisor(18))
    n14_gate_skip, n14_gate = S492.best_ladder(14, 14)
    n18_gate_skip, n18_gate = S492.best_ladder(18, 18)
    return (
        ("initial n=14", 14, initial(14)),
        ("initial n=18", 18, initial(18)),
        (f"n=14 lpd ladder skip {n14_lpd_skip}", 14, n14_lpd),
        (f"n=18 lpd ladder skip {n18_lpd_skip}", 18, n18_lpd),
        (f"n=14 gate ladder skip {n14_gate_skip}", 14, n14_gate),
        (f"n=18 gate ladder skip {n18_gate_skip}", 18, n18_gate),
    )


def selected_times(moving_speeds: tuple[int, ...], n: int) -> tuple[tuple[str, Fraction], ...]:
    report = S356.report("selected", list(moving_speeds))
    times: list[tuple[str, Fraction]] = [
        ("unit", Fraction(1, n)),
        ("half-unit", Fraction(1, 2 * n)),
    ]
    if report.witness is not None:
        times.append(("gap-mid", report.witness))
    if report.boundary_witness is not None:
        times.append(("boundary", report.boundary_witness))
    seen: set[Fraction] = set()
    out: list[tuple[str, Fraction]] = []
    for tag, t in times:
        if t not in seen:
            seen.add(t)
            out.append((tag, t))
    return tuple(out)


def short_score_hist(hist: tuple[tuple[int, int], ...]) -> str:
    return " ".join(f"{score}:{count}" for score, count in hist)


def print_large_rows() -> None:
    shortlist = {
        "phase_half",
        "local_moat_min",
        "origin_blocker_relief",
        "k2_relief_pressure",
        "threshold_deficit_pressure",
        "bracket_relief_pressure",
    }
    criteria = tuple(item for item in CRITERIA if item.name in shortlist)

    print("SELECTED n=14/n=18 LRC ROWS")
    print("=" * 100)
    print("H is omitted above n=14; score histograms, strict cycles, SCCs, and ties remain exact.")
    print()
    for label, n, moving in selected_examples():
        speeds = (0,) + moving
        report = S356.report(label, list(moving))
        summary = S360.summarize(list(moving))
        print(f"[{label}]")
        print(
            f"  class={summary.classification} gap/th={fmt(report.max_gap / report.threshold)} "
            f"unprotected={summary.unprotected_count} speeds={moving}"
        )
        for tag, t in selected_times(moving, n):
            geom = geometry(speeds, t)
            print(
                f"  time={tag:<9} t={fmt(t):>14} "
                f"origin_clear={fmt(geom.stationary_min / geom.threshold):>7} "
                f"bracket={fmt(min(geom.origin_left, geom.origin_right) / geom.threshold):>7} "
                f"safe={geom.safe_gap_count:>2} lonely={len(geom.lonely_vertices):>2} "
                f"max_gap/th={fmt(geom.max_gap / geom.threshold):>7}"
            )
            for criterion in criteria:
                snap = snapshot(criterion, speeds, t)
                fp = fingerprint(snap, compute_h=True)
                h_text = fmt(fp.hamiltonian_paths)
                origin_text = fmt(fp.origin_score)
                print(
                    f"    {criterion.name:<28} "
                    f"v={fp.vertex_count:>2} H={h_text:>12} "
                    f"score=[{short_score_hist(fp.score_hist)}] "
                    f"w={fp.score_width:>2} o={origin_text:>3} "
                    f"strict={fp.strict_arcs:>3} ties={fp.ties:>3} "
                    f"tri={fp.strict_triangles:>3} scc={fp.largest_strict_scc:>2}"
                )
        print()


def print_synthesis(summary: dict[str, dict[str, object]]) -> None:
    print("SYNTHESIS")
    print("=" * 100)
    ordered = sorted(summary.items(), key=lambda item: -float(item[1]["best_abs"]))
    print("Most useful criteria in this audit:")
    for name, row in ordered[:5]:
        feature, target = row["best_pair"]  # type: ignore[misc]
        print(
            f"  {name}: {feature}->{target}, "
            f"mean_abs_spearman={float(row['best_abs']):.3f}, tie_rate={float(row['tie_rate']):.3f}"
        )
    print()
    print("Operational conclusion:")
    print("  1. Half-turn phase H remains the global spread meter: low H is the")
    print("     unanchored open-semicircle signal, while high H often means many")
    print("     locally lonely vertices.")
    print("  2. Scalar rank gauges such as local_moat_min and safe_deficit are better")
    print("     marked LRC meters, but their H is mostly tie-path or transitive noise.")
    print("     Use the marked origin score and score histogram instead.")
    print("  3. Origin blocker relief is not a good clock-wide scalar by itself, but")
    print("     it is a useful selected-snapshot diagnostic for 'why is runner 0 not")
    print("     lonely?': read its strict/tie pattern before trusting its completed H.")
    print("  4. k2, threshold-deficit, and bracket-relief pressure gauges should be read")
    print("     by strict SCCs and directed triangles.  In the selected n=14/n=18 rows")
    print("     they still peel rather than forming a cyclic pressure core.")
    print()
    print("Proposed LRC tournament metric vector:")
    print("  (phase_H, phase_score_width, origin_marked_score, safe_deficit_score_hist,")
    print("   origin_blocker_score_width, pressure_largest_SCC, pressure_triangles,")
    print("   strict_tie_rate, edge_flip_rate)")
    print()
    print("This is the useful replacement for a single loneliness scalar: H supplies")
    print("the cyclic spread coordinate, while marked scores and pressure-core shape")
    print("supply the LRC threshold coordinate.")


def main() -> None:
    print_methodology()
    summary = criterion_scorecard()
    print_scorecard(summary)
    print_large_rows()
    print_synthesis(summary)


if __name__ == "__main__":
    main()
