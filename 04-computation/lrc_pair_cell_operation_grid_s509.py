#!/usr/bin/env python3
"""
lrc_pair_cell_operation_grid_s509.py

codex-2026-06-01 S509

Extend the S506b LRC arc-criteria search by lifting one level: make the
unordered runner-pairs themselves into tournament vertices.

Methodology:

* objects: pair-cells {i,j} of runners, including the stationary observer 0;
* pairwise observable: a scalar attached to each pair-cell, then compared
  between pair-cells;
* switch/gauge: larger scalar points to smaller scalar;
* tie Hamiltonian path: lexicographic order on pair-cells.

The tested pair-cell scalars deliberately mix geometry and arithmetic:

* current circular distance;
* current LRC danger deficit below 1/n;
* half-turn wall frequency |v_i-v_j|;
* dyadic row v2(|v_i-v_j|);
* odd core of |v_i-v_j| along the x+2 row;
* same-odd-chain status in the x*2 branch grid;
* product-sum proximity to (a-1)(b-1)=n-1.

This is the place where A000568-like tournament structure enters the LRC
search: A000568 counts unlabeled complete binary relations, and Burnside's
formula only keeps odd cycle partitions.  LRC pair-cells have the same
complete-relation surface, but their useful gauges live on the two operation
directions: odd-core x+2 chains and dyadic x*2 branches.
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from importlib.machinery import SourceFileLoader
from itertools import combinations
from math import comb, factorial, gcd, log2, prod, sqrt
from pathlib import Path
from statistics import mean
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
S506B = SourceFileLoader(
    "lrc_arc_criteria_loneliness_s506",
    str(ROOT / "04-computation" / "lrc_arc_criteria_loneliness_s506.py"),
).load_module()
S365 = SourceFileLoader(
    "natural_operation_digraphs_s365",
    str(ROOT / "04-computation" / "natural_operation_digraphs_s365.py"),
).load_module()


H_LIMIT = 15


@dataclass(frozen=True)
class EdgeCriterion:
    name: str
    observable: str
    switch: str
    dynamic: bool


@dataclass(frozen=True)
class EdgeSnapshot:
    criterion: EdgeCriterion
    labels: tuple[tuple[int, int], ...]
    out: tuple[int, ...]
    strict_out: tuple[int, ...]
    ties: int
    strict_arcs: int
    values: tuple[Any, ...]


@dataclass(frozen=True)
class EdgeFingerprint:
    vertex_count: int
    hamiltonian_paths: int | None
    score_hist: tuple[tuple[int, int], ...]
    score_width: int
    score_variance: float
    top_score: int
    origin_pair_max_score: int
    origin_pair_mean_score: float
    strict_triangles: int
    complete_triangles: int
    largest_strict_scc: int
    ties: int
    strict_arcs: int


EDGE_CRITERIA: tuple[EdgeCriterion, ...] = (
    EdgeCriterion(
        "edge_danger_deficit",
        "pair-cell deficit max(0,1/n-dist(i,j,t))",
        "larger current LRC danger deficit points to smaller deficit",
        True,
    ),
    EdgeCriterion(
        "edge_phase_distance",
        "pair-cell circular distance dist(i,j,t)",
        "larger distance points to smaller distance",
        True,
    ),
    EdgeCriterion(
        "edge_wall_frequency",
        "pair-cell speed gap |v_i-v_j|",
        "larger half-turn wall frequency points to smaller frequency",
        False,
    ),
    EdgeCriterion(
        "edge_dyadic_row",
        "dyadic row v2(|v_i-v_j|)",
        "higher x*2 branch row points to lower row",
        False,
    ),
    EdgeCriterion(
        "edge_odd_core",
        "odd core of |v_i-v_j| on the x+2 row",
        "smaller odd core points to larger odd core",
        False,
    ),
    EdgeCriterion(
        "edge_same_odd_chain",
        "whether both speeds lie on the same odd-core doubling chain",
        "same-chain and deeper-row pairs point to other pairs",
        False,
    ),
    EdgeCriterion(
        "edge_product_sum_defect",
        "residue-pair distance from (a-1)(b-1)=n-1",
        "closer product-sum critical pair points to farther pair",
        False,
    ),
)


def fmt(value: Fraction | int | float | None) -> str:
    if value is None:
        return "-"
    if isinstance(value, Fraction):
        return S506B.fmt(value)
    if isinstance(value, float):
        return f"{value:.3f}"
    return str(value)


def v2(value: int) -> int:
    value = abs(value)
    if value == 0:
        return 0
    out = 0
    while value % 2 == 0:
        out += 1
        value //= 2
    return out


def odd_core(value: int) -> int:
    value = abs(value)
    if value == 0:
        return 0
    return value >> v2(value)


def compare_values(left: Any, right: Any) -> int:
    return (left > right) - (left < right)


def pair_labels(n: int) -> tuple[tuple[int, int], ...]:
    return tuple(combinations(range(n), 2))


def circular_distances(pos: tuple[Fraction, ...]) -> dict[tuple[int, int], Fraction]:
    out: dict[tuple[int, int], Fraction] = {}
    for i, j in pair_labels(len(pos)):
        out[(i, j)] = S506B.circular_distance(pos[i], pos[j])
    return out


def edge_value(
    criterion: EdgeCriterion,
    speeds_with_observer: tuple[int, ...],
    t: Fraction,
    pair: tuple[int, int],
) -> Any:
    n = len(speeds_with_observer)
    i, j = pair
    speed_i = speeds_with_observer[i]
    speed_j = speeds_with_observer[j]
    gap = abs(speed_j - speed_i)

    if criterion.name == "edge_danger_deficit":
        pos = S506B.positions(speeds_with_observer, t)
        dist = S506B.circular_distance(pos[i], pos[j])
        threshold = Fraction(1, n)
        return max(Fraction(0), threshold - dist)

    if criterion.name == "edge_phase_distance":
        pos = S506B.positions(speeds_with_observer, t)
        return S506B.circular_distance(pos[i], pos[j])

    if criterion.name == "edge_wall_frequency":
        return gap

    if criterion.name == "edge_dyadic_row":
        return (v2(gap), gap)

    if criterion.name == "edge_odd_core":
        return (-odd_core(gap), v2(gap))

    if criterion.name == "edge_same_odd_chain":
        core_i = odd_core(speed_i)
        core_j = odd_core(speed_j)
        same = int(core_i != 0 and core_i == core_j)
        return (same, min(v2(speed_i), v2(speed_j)), -abs(v2(speed_i) - v2(speed_j)))

    if criterion.name == "edge_product_sum_defect":
        # Two-factor product-sum witnesses at total arity n satisfy
        # (A-1)(B-1)=n-1.  Use residues as the shifted A-1,B-1 coordinates.
        # Smaller defect is more critical, so return the negative defect.
        a = abs(speed_i) % n
        b = abs(speed_j) % n
        return (-abs(a * b - (n - 1)), -(a + b))

    raise ValueError(f"unknown edge criterion {criterion.name}")


def edge_snapshot(
    criterion: EdgeCriterion,
    speeds_with_observer: tuple[int, ...],
    t: Fraction,
) -> EdgeSnapshot:
    labels = pair_labels(len(speeds_with_observer))
    values = tuple(edge_value(criterion, speeds_with_observer, t, pair) for pair in labels)
    out = [0] * len(labels)
    strict_out = [0] * len(labels)
    ties = 0
    strict_arcs = 0
    for a_idx, b_idx in combinations(range(len(labels)), 2):
        cmp = compare_values(values[a_idx], values[b_idx])
        if cmp > 0:
            winner, loser = a_idx, b_idx
            strict_arcs += 1
        elif cmp < 0:
            winner, loser = b_idx, a_idx
            strict_arcs += 1
        else:
            winner, loser = a_idx, b_idx
            ties += 1
        out[winner] |= 1 << loser
        if cmp != 0:
            strict_out[winner] |= 1 << loser
    return EdgeSnapshot(
        criterion=criterion,
        labels=labels,
        out=tuple(out),
        strict_out=tuple(strict_out),
        ties=ties,
        strict_arcs=strict_arcs,
        values=values,
    )


@lru_cache(maxsize=None)
def h_count(out: tuple[int, ...]) -> int:
    n = len(out)
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for vertex in range(n):
        dp[1 << vertex][vertex] = 1
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
    for i, bits in enumerate(out):
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

    for vertex in range(n):
        if vertex not in seen:
            dfs(vertex)

    sizes: list[int] = []
    seen.clear()
    for start in reversed(order):
        if start in seen:
            continue
        todo = deque([start])
        seen.add(start)
        size = 0
        while todo:
            vertex = todo.pop()
            size += 1
            for nxt in reverse[vertex]:
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


def fingerprint(snap: EdgeSnapshot, compute_h: bool) -> EdgeFingerprint:
    scores = score_list(snap.out)
    origin_scores = [
        scores[idx]
        for idx, pair in enumerate(snap.labels)
        if 0 in pair
    ]
    avg = mean(scores) if scores else 0.0
    variance = mean([(score - avg) ** 2 for score in scores]) if scores else 0.0
    strict_scc = scc_sizes(snap.strict_out)
    return EdgeFingerprint(
        vertex_count=len(snap.labels),
        hamiltonian_paths=h_count(snap.out) if compute_h and len(snap.labels) <= H_LIMIT else None,
        score_hist=tuple(sorted(Counter(scores).items())),
        score_width=max(scores) - min(scores) if scores else 0,
        score_variance=variance,
        top_score=max(scores) if scores else 0,
        origin_pair_max_score=max(origin_scores) if origin_scores else 0,
        origin_pair_mean_score=mean(origin_scores) if origin_scores else 0.0,
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


def feature_series(snaps: list[EdgeSnapshot]) -> dict[str, list[float]]:
    out: dict[str, list[float]] = defaultdict(list)
    for snap in snaps:
        fp = fingerprint(snap, compute_h=True)
        denom = max(1, fp.vertex_count - 1)
        total_pairs = max(1, comb(fp.vertex_count, 2))
        total_triples = max(1, comb(fp.vertex_count, 3)) if fp.vertex_count >= 3 else 1
        out["H"].append(float(fp.hamiltonian_paths or 0))
        out["score_width_norm"].append(float(fp.score_width / denom))
        out["score_variance_norm"].append(float(fp.score_variance / (denom * denom)))
        out["top_score_norm"].append(float(fp.top_score / denom))
        out["origin_pair_max_norm"].append(float(fp.origin_pair_max_score / denom))
        out["origin_pair_mean_norm"].append(float(fp.origin_pair_mean_score / denom))
        out["tie_rate"].append(float(fp.ties / total_pairs))
        out["strict_arc_rate"].append(float(fp.strict_arcs / total_pairs))
        out["strict_triangle_density"].append(float(fp.strict_triangles / total_triples))
        out["complete_triangle_density"].append(float(fp.complete_triangles / total_triples))
        out["strict_scc_norm"].append(float(fp.largest_strict_scc / fp.vertex_count))
    return dict(out)


def target_series(
    speeds_with_observer: tuple[int, ...],
    times: tuple[Fraction, ...],
    geoms: list[Any],
) -> dict[str, list[float]]:
    close_counts: list[float] = []
    threshold = Fraction(1, len(speeds_with_observer))
    for t in times:
        pos = S506B.positions(speeds_with_observer, t)
        dists = circular_distances(pos)
        close_counts.append(float(sum(1 for dist in dists.values() if dist < threshold)))
    return {
        "origin_clearance": [float(g.stationary_min / g.threshold) for g in geoms],
        "origin_bracket": [float(min(g.origin_left, g.origin_right) / g.threshold) for g in geoms],
        "safe_gap_count": [float(g.safe_gap_count) for g in geoms],
        "lonely_vertices": [float(len(g.lonely_vertices)) for g in geoms],
        "max_gap": [float(g.max_gap) for g in geoms],
        "close_pair_count": close_counts,
    }


def criterion_scorecard() -> dict[str, dict[str, object]]:
    summary: dict[str, dict[str, object]] = {}
    for criterion in EDGE_CRITERIA:
        best_scores: list[tuple[float, str, str]] = []
        h_safe: list[float] = []
        h_close: list[float] = []
        origin_pair_clearance: list[float] = []
        flip_rates: list[float] = []
        tie_rates: list[float] = []
        scc_rates: list[float] = []

        for _, speeds in S506B.small_families():
            times = S506B.exact_cell_times(speeds)
            geoms = [S506B.geometry(speeds, t) for t in times]
            snaps = [edge_snapshot(criterion, speeds, t) for t in times]
            features = feature_series(snaps)
            targets = target_series(speeds, times, geoms)
            local_best = (0.0, "-", "-")
            for feature_name, xs in features.items():
                for target_name, ys in targets.items():
                    rho = spearman(xs, ys)
                    if abs(rho) > abs(local_best[0]):
                        local_best = (rho, feature_name, target_name)
            best_scores.append(local_best)
            h_safe.append(spearman(features["H"], targets["safe_gap_count"]))
            h_close.append(spearman(features["H"], targets["close_pair_count"]))
            origin_pair_clearance.append(
                spearman(features["origin_pair_max_norm"], targets["origin_clearance"])
            )
            keys = [edge_key(snap.out) for snap in snaps]
            flips = [edge_flips(a, b) for a, b in zip(keys, keys[1:])]
            pair_count = max(1, comb(len(snaps[0].labels), 2))
            if flips:
                flip_rates.append(mean(flips) / pair_count)
            total_pairs = pair_count * len(snaps)
            tie_rates.append(sum(snap.ties for snap in snaps) / total_pairs)
            scc_rates.append(mean(feature_series(snaps)["strict_scc_norm"]))

        best_abs = mean(abs(row[0]) for row in best_scores)
        feature_votes = Counter((row[1], row[2]) for row in best_scores)
        best_pair, _ = feature_votes.most_common(1)[0]
        summary[criterion.name] = {
            "criterion": criterion,
            "best_abs": best_abs,
            "best_pair": best_pair,
            "h_safe": mean(h_safe),
            "h_close": mean(h_close),
            "origin_pair_clearance": mean(origin_pair_clearance),
            "flip_rate": mean(flip_rates),
            "tie_rate": mean(tie_rates),
            "scc_rate": mean(scc_rates),
        }
    return summary


def short_score_hist(hist: tuple[tuple[int, int], ...], limit: int = 12) -> str:
    parts = [f"{score}:{count}" for score, count in hist[:limit]]
    if len(hist) > limit:
        parts.append("...")
    return " ".join(parts)


def print_methodology() -> None:
    print("LRC pair-cell operation-grid Tournament Analysis")
    print("=" * 100)
    print("Tie Hamiltonian path: lexicographic pair-cell order (0,1)->(0,2)->... .")
    print("Strict arcs compare pair-cell scalar values; equal values use the tie path.")
    print()
    print("Pair-cell criteria:")
    for item in EDGE_CRITERIA:
        print(f"  {item.name:<25} dynamic={str(item.dynamic):<5} observable={item.observable}")
        print(f"    switch={item.switch}")
    print()


def print_scorecard(summary: dict[str, dict[str, object]]) -> None:
    print("SMALL EXACT CLOCK PAIR-CELL SCORECARD")
    print("=" * 100)
    print(
        "criterion                  best_abs  common_best_feature->target"
        "                  H~safe H~close originPair~clear flips ties scc"
    )
    for name, row in sorted(summary.items(), key=lambda item: -float(item[1]["best_abs"])):
        best_feature, best_target = row["best_pair"]  # type: ignore[misc]
        print(
            f"{name:<25} "
            f"{float(row['best_abs']):>8.3f}  "
            f"{best_feature}->{best_target:<30} "
            f"{float(row['h_safe']):>+7.3f} "
            f"{float(row['h_close']):>+7.3f} "
            f"{float(row['origin_pair_clearance']):>+10.3f} "
            f"{float(row['flip_rate']):>5.3f} "
            f"{float(row['tie_rate']):>5.3f} "
            f"{float(row['scc_rate']):>5.3f}"
        )
    print()
    print("Reading: this is a tournament on runner-pairs, so H is computed only")
    print(f"when the pair-cell count is at most {H_LIMIT}.  For larger rows, score")
    print("shape, SCCs, triangles, and tie rates are the intended fingerprints.")
    print()


def print_selected_rows() -> None:
    shortlist = {
        "edge_danger_deficit",
        "edge_phase_distance",
        "edge_dyadic_row",
        "edge_odd_core",
        "edge_same_odd_chain",
        "edge_product_sum_defect",
    }
    criteria = tuple(item for item in EDGE_CRITERIA if item.name in shortlist)
    print("SELECTED n=14/n=18 PAIR-CELL ROWS")
    print("=" * 100)
    print("H is omitted for these pair-cell tournaments; pair-cell counts are 91 or 153.")
    print()
    for label, n, moving in S506B.selected_examples():
        speeds = (0,) + moving
        report = S506B.S356.report(label, list(moving))
        print(f"[{label}]")
        print(
            f"  gap/th={fmt(report.max_gap / report.threshold)} "
            f"pair_cells={comb(len(speeds), 2)} speeds={moving}"
        )
        for tag, t in S506B.selected_times(moving, n):
            geom = S506B.geometry(speeds, t)
            print(
                f"  time={tag:<9} t={fmt(t):>14} "
                f"origin_clear={fmt(geom.stationary_min / geom.threshold):>7} "
                f"safe={geom.safe_gap_count:>2} lonely={len(geom.lonely_vertices):>2} "
                f"max_gap/th={fmt(geom.max_gap / geom.threshold):>7}"
            )
            for criterion in criteria:
                snap = edge_snapshot(criterion, speeds, t)
                fp = fingerprint(snap, compute_h=False)
                total_pairs = comb(fp.vertex_count, 2)
                print(
                    f"    {criterion.name:<25} "
                    f"score=[{short_score_hist(fp.score_hist)}] "
                    f"w={fp.score_width:>3} top={fp.top_score:>3} "
                    f"opair={fp.origin_pair_max_score:>3}/{fmt(fp.origin_pair_mean_score):>5} "
                    f"strict={fp.strict_arcs:>5} ties={fp.ties:>5} "
                    f"tri={fp.strict_triangles:>6} scc={fp.largest_strict_scc:>3} "
                    f"tie%={fp.ties / max(1, total_pairs):.3f}"
                )
        print()


def odd_partitions(n: int) -> list[tuple[int, ...]]:
    out: list[tuple[int, ...]] = []

    def rec(remaining: int, max_part: int, parts: list[int]) -> None:
        if remaining == 0:
            out.append(tuple(parts))
            return
        start = min(max_part, remaining)
        if start % 2 == 0:
            start -= 1
        for part in range(start, 0, -2):
            rec(remaining - part, part, parts + [part])

    rec(n, n, [])
    return out


def partition_counts(parts: tuple[int, ...]) -> Counter[int]:
    return Counter(parts)


def a000568_term(parts: tuple[int, ...]) -> Fraction:
    counts = partition_counts(parts)
    exponent = (
        sum(counts[r] * counts[s] * gcd(r, s) for r in counts for s in counts)
        - sum(counts.values())
    ) // 2
    denominator = prod((part ** mult) * factorial(mult) for part, mult in counts.items())
    return Fraction(2**exponent, denominator)


def a000568_burnside(n: int) -> int:
    return int(sum(a000568_term(parts) for parts in odd_partitions(n)))


def print_a000568_operation_grid() -> None:
    print("A000568 AND THE x+2 / x*2 GRID")
    print("=" * 100)
    print("A000568 Burnside terms use only odd cycle partitions.")
    print("That is the tournament-counting version of an odd-core rule: even")
    print("cycle orbits cannot carry a coherent complete orientation.")
    print()
    print("n  A000568(n)  #odd_partitions  largest odd partitions contributing")
    print("-" * 100)
    for n in range(1, 11):
        parts = odd_partitions(n)
        examples = ", ".join("+".join(map(str, item)) for item in parts[:4])
        print(f"{n:>2} {a000568_burnside(n):>11} {len(parts):>16}  {examples}")
    print()
    print("2-adic operation grid up to 40:")
    print("  top row moves by x+2; each top-row odd core branches by x*2.")
    for odd in range(1, 20, 2):
        chain: list[int] = []
        value = odd
        while value <= 40:
            chain.append(value)
            value *= 2
        print(f"  odd core {odd:>2}: {' -> '.join(map(str, chain))}")
    print()
    print("Product-sum two-factor critical pairs for LRC denominators:")
    print("  total arity k has r=k-2 ones and (a-1)(b-1)=k-1.")
    for k in (14, 18, 21, 24):
        sols = S365.two_factor_solutions_for_arity(k)
        text = ", ".join(f"({a},{b})->{p}" for a, b, p in sols)
        seeds = S365.enumerate_product_sum_seeds(max_k=k, max_product=4 * k).get(k, [])
        best = seeds[0] if seeds else None
        best_text = S365.compact_witness(best) if best else "-"
        print(f"  k={k:>2}: two-factor {text or '-'}; best seed {best_text}")
    print()


def print_synthesis(summary: dict[str, dict[str, object]]) -> None:
    print("SYNTHESIS")
    print("=" * 100)
    ordered = sorted(summary.items(), key=lambda item: -float(item[1]["best_abs"]))
    print("Most useful pair-cell criteria in this audit:")
    for name, row in ordered[:5]:
        feature, target = row["best_pair"]  # type: ignore[misc]
        print(
            f"  {name}: {feature}->{target}, "
            f"mean_abs_spearman={float(row['best_abs']):.3f}, "
            f"ties={float(row['tie_rate']):.3f}, scc={float(row['scc_rate']):.3f}"
        )
    print()
    print("Interpretation:")
    print("  1. Runner-level Tournament Analysis asks which runner wins a relation.")
    print("     Pair-cell Tournament Analysis asks which arc/distance cell is the")
    print("     dominant obstruction.  LRC needs both levels.")
    print("  2. A000568's odd-partition Burnside formula is the unlabeled tournament")
    print("     version of odd-core survival.  The LRC operation grid splits this")
    print("     into horizontal x+2 odd-core movement and vertical x*2 doubling.")
    print("  3. Product-sum defects are not a replacement for geometry.  They are a")
    print("     static arithmetic pressure coordinate saying which pair-cells sit near")
    print("     addition/multiplication critical pairs.")
    print("  4. The useful loneliness metric is therefore two-layered:")
    print("       runner vector from S506b")
    print("       + pair-cell operation vector from this script.")


def main() -> None:
    print_methodology()
    summary = criterion_scorecard()
    print_scorecard(summary)
    print_selected_rows()
    print_a000568_operation_grid()
    print_synthesis(summary)


if __name__ == "__main__":
    main()
