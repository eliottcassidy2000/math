#!/usr/bin/env python3
"""
lrc_operation_grid_arc_criteria_s511.py

codex-2026-06-01 S511

Explore LRC Tournament Analysis arc criteria that mix three layers:

1. runner-level phase/threshold gauges from S506;
2. pair-cell operation-grid labels from S509;
3. runner aggregates induced by pair-cells, so static arithmetic branch data
   can become a moving loneliness pressure signal.

Tournament Analysis declaration:

* objects: runners including stationary observer 0;
* pairwise observable: either an existing pairwise phase/pressure switch, or a
  runner scalar induced by local gaps, incident pair deficits, and arithmetic
  operation-grid labels;
* switch/gauge: larger runner scalar points to smaller runner scalar, while
  imported S506 criteria keep their native switches;
* tie Hamiltonian path: increasing runner label 0 -> 1 -> ... -> N-1.

The goal is not to collapse the LRC to one magic H value.  The goal is to find
which tournament fingerprints become useful loneliness-related coordinates,
and which arithmetic fingerprints are better read as A000568-style fiber
labels over the moving LRC quotient path.

Stored output:
    05-knowledge/results/lrc_operation_grid_arc_criteria_s511.out
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from fractions import Fraction
from importlib.machinery import SourceFileLoader
from itertools import combinations
from math import comb, log2, sqrt
from pathlib import Path
from statistics import mean
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
S506 = SourceFileLoader(
    "lrc_arc_criteria_loneliness_s506",
    str(ROOT / "04-computation" / "lrc_arc_criteria_loneliness_s506.py"),
).load_module()
S509 = SourceFileLoader(
    "lrc_pair_cell_operation_grid_s509",
    str(ROOT / "04-computation" / "lrc_pair_cell_operation_grid_s509.py"),
).load_module()


H_LIMIT = 14


@dataclass(frozen=True)
class ArcCriterion:
    name: str
    layer: str
    observable: str
    switch: str
    dynamic: bool
    source: str = "local"


@dataclass(frozen=True)
class RunnerSnapshot:
    criterion: ArcCriterion
    labels: tuple[int, ...]
    out: tuple[int, ...]
    strict_out: tuple[int, ...]
    ties: int
    strict_arcs: int
    values: tuple[Any, ...] = ()


@dataclass(frozen=True)
class ScorecardRow:
    criterion: ArcCriterion
    best_abs: float
    best_pair: tuple[str, str]
    h_origin: float
    h_close: float
    origin_score_origin: float
    top_score_danger: float
    flip_rate: float
    tie_rate: float
    scc_rate: float


@dataclass(frozen=True)
class CriteriaProfile:
    name: str
    signal: int
    threshold: int
    operation_grid: int
    dynamic_motion: int
    novelty: int
    projection_risk: int
    score: int


LOCAL_CRITERIA: tuple[ArcCriterion, ...] = (
    ArcCriterion(
        "origin_clearance_rank",
        "runner geometry",
        "runner distance from observer 0",
        "farther from 0 points to closer to 0",
        True,
    ),
    ArcCriterion(
        "local_two_guard",
        "runner geometry",
        "min(left adjacent gap, right adjacent gap)",
        "larger two-sided guard points to smaller guard",
        True,
    ),
    ArcCriterion(
        "incident_danger_sum",
        "pair-cell aggregate",
        "sum_j max(0,1/N-dist(i,j,t))",
        "more incident threshold danger points to less danger",
        True,
    ),
    ArcCriterion(
        "incident_close_count",
        "pair-cell aggregate",
        "number of incident pair-cells below the LRC threshold",
        "more close incident pair-cells points to fewer",
        True,
    ),
    ArcCriterion(
        "dyadic_branch_pressure",
        "operation grid",
        "sum_j v2(|v_i-v_j|)",
        "higher x*2 branch pressure points to lower",
        False,
    ),
    ArcCriterion(
        "odd_core_branch_count",
        "operation grid",
        "number of distinct odd cores among incident speed gaps",
        "more x+2 branch diversity points to less",
        False,
    ),
    ArcCriterion(
        "same_odd_chain_degree",
        "operation grid",
        "count of runners sharing the same odd-core doubling chain",
        "more same-chain incidence points to less",
        False,
    ),
    ArcCriterion(
        "additive_shadow_degree",
        "operation grid",
        "participation in speed equations a+b=c or |a-b|=c",
        "more additive-shadow incidence points to less",
        False,
    ),
    ArcCriterion(
        "multiplicative_shadow_degree",
        "operation grid",
        "divisibility incidence among nonzero speeds",
        "more multiplicative-shadow incidence points to less",
        False,
    ),
    ArcCriterion(
        "product_sum_interface",
        "operation grid",
        "sum_j 1/(1+|(v_i mod N)(v_j mod N)-(N-1)|)",
        "closer to product-sum interface points to farther",
        False,
    ),
    ArcCriterion(
        "dyadic_danger_curvature",
        "hybrid",
        "incident danger deficit weighted by 1+v2(|v_i-v_j|)",
        "more deep-branch danger points to less",
        True,
    ),
    ArcCriterion(
        "odd_chain_danger",
        "hybrid",
        "incident danger deficit restricted to same odd-core chains",
        "more same-chain danger points to less",
        True,
    ),
    ArcCriterion(
        "additive_danger_interface",
        "hybrid",
        "incident danger deficit weighted by additive-shadow pair incidence",
        "more additive danger points to less",
        True,
    ),
    ArcCriterion(
        "multiplicative_danger_interface",
        "hybrid",
        "incident danger deficit weighted by divisibility pair incidence",
        "more multiplicative danger points to less",
        True,
    ),
    ArcCriterion(
        "product_sum_danger",
        "hybrid",
        "incident danger deficit weighted by product-sum closeness",
        "more product-sum danger points to less",
        True,
    ),
)


IMPORTED_NAMES = (
    "phase_half",
    "lrc_close_sector",
    "local_moat_min",
    "safe_deficit",
    "threshold_deficit_pressure",
)


def imported_criteria() -> tuple[ArcCriterion, ...]:
    by_name = {item.name: item for item in S506.CRITERIA}
    out: list[ArcCriterion] = []
    for name in IMPORTED_NAMES:
        item = by_name[name]
        out.append(
            ArcCriterion(
                name=f"s506_{name}",
                layer="imported runner gauge",
                observable=item.observable,
                switch=item.switch,
                dynamic=True,
                source=name,
            )
        )
    return tuple(out)


CRITERIA: tuple[ArcCriterion, ...] = imported_criteria() + LOCAL_CRITERIA


def fmt(value: Fraction | int | float | None) -> str:
    if value is None:
        return "-"
    if isinstance(value, Fraction):
        return S506.fmt(value)
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


def nonzero_speed_set(speeds: tuple[int, ...]) -> set[int]:
    return {abs(speed) for speed in speeds if speed != 0}


def pair_product_closeness(
    speeds_with_observer: tuple[int, ...], i: int, j: int
) -> Fraction:
    n = len(speeds_with_observer)
    a = abs(speeds_with_observer[i]) % n
    b = abs(speeds_with_observer[j]) % n
    defect = abs(a * b - (n - 1))
    return Fraction(1, 1 + defect)


def additive_pair_weight(speeds_with_observer: tuple[int, ...], i: int, j: int) -> int:
    si = abs(speeds_with_observer[i])
    sj = abs(speeds_with_observer[j])
    if si == 0 or sj == 0:
        return 0
    speeds = nonzero_speed_set(speeds_with_observer)
    weight = 0
    if si + sj in speeds:
        weight += 1
    if abs(si - sj) in speeds:
        weight += 1
    return weight


def multiplicative_pair_weight(
    speeds_with_observer: tuple[int, ...], i: int, j: int
) -> int:
    si = abs(speeds_with_observer[i])
    sj = abs(speeds_with_observer[j])
    if si == 0 or sj == 0 or si == sj:
        return 0
    big = max(si, sj)
    small = min(si, sj)
    return int(small != 0 and big % small == 0)


def additive_degrees(speeds_with_observer: tuple[int, ...]) -> tuple[int, ...]:
    speeds = nonzero_speed_set(speeds_with_observer)
    out = [0] * len(speeds_with_observer)
    for i, j in combinations(range(len(speeds_with_observer)), 2):
        weight = additive_pair_weight(speeds_with_observer, i, j)
        if weight:
            out[i] += weight
            out[j] += weight
    for i, speed in enumerate(speeds_with_observer):
        if speed == 0:
            continue
        # Count direct decompositions speed = a+b with a,b present.
        target = abs(speed)
        for a in speeds:
            b = target - a
            if b in speeds and a <= b:
                out[i] += 1
    return tuple(out)


def multiplicative_degrees(speeds_with_observer: tuple[int, ...]) -> tuple[int, ...]:
    out = [0] * len(speeds_with_observer)
    for i, j in combinations(range(len(speeds_with_observer)), 2):
        weight = multiplicative_pair_weight(speeds_with_observer, i, j)
        if weight:
            out[i] += weight
            out[j] += weight
    return tuple(out)


def runner_values(
    criterion: ArcCriterion,
    speeds_with_observer: tuple[int, ...],
    t: Fraction,
) -> tuple[Any, ...]:
    pos = S506.positions(speeds_with_observer, t)
    n = len(pos)
    threshold = Fraction(1, n)
    dist = S506.distance_matrix(pos)
    min_bracket, _sum_bracket, _pair_bracket = S506.bracket_values(pos)
    add_deg = additive_degrees(speeds_with_observer)
    mult_deg = multiplicative_degrees(speeds_with_observer)

    incident_deficit: list[Fraction] = []
    incident_close: list[int] = []
    incident_open: list[Fraction] = []
    dyadic_pressure: list[int] = []
    odd_core_count: list[int] = []
    same_chain: list[int] = []
    product_interface: list[Fraction] = []
    dyadic_danger: list[Fraction] = []
    odd_danger: list[Fraction] = []
    additive_danger: list[Fraction] = []
    multiplicative_danger: list[Fraction] = []
    product_danger: list[Fraction] = []

    for i in range(n):
        deficits: list[Fraction] = []
        close = 0
        open_mass = Fraction(0)
        dyadic_total = 0
        odd_cores: set[int] = set()
        same_total = 0
        product_total = Fraction(0)
        dyadic_weighted = Fraction(0)
        odd_weighted = Fraction(0)
        additive_weighted = Fraction(0)
        multiplicative_weighted = Fraction(0)
        product_weighted = Fraction(0)

        for j in range(n):
            if i == j:
                continue
            gap = abs(speeds_with_observer[i] - speeds_with_observer[j])
            deficit = max(Fraction(0), threshold - dist[i][j])
            deficits.append(deficit)
            if dist[i][j] < threshold:
                close += 1
            if dist[i][j] > threshold:
                open_mass += dist[i][j] - threshold
            if gap:
                dyadic_total += v2(gap)
                odd_cores.add(odd_core(gap))

            core_i = odd_core(speeds_with_observer[i])
            core_j = odd_core(speeds_with_observer[j])
            same = int(core_i != 0 and core_i == core_j)
            same_total += same

            product_close = pair_product_closeness(speeds_with_observer, i, j)
            add_weight = additive_pair_weight(speeds_with_observer, i, j)
            mult_weight = multiplicative_pair_weight(speeds_with_observer, i, j)
            product_total += product_close
            dyadic_weighted += deficit * (1 + v2(gap))
            odd_weighted += deficit * same
            additive_weighted += deficit * (1 + add_weight)
            multiplicative_weighted += deficit * (1 + mult_weight)
            product_weighted += deficit * product_close

        incident_deficit.append(sum(deficits, Fraction(0)))
        incident_close.append(close)
        incident_open.append(open_mass)
        dyadic_pressure.append(dyadic_total)
        odd_core_count.append(len(odd_cores))
        same_chain.append(same_total)
        product_interface.append(product_total)
        dyadic_danger.append(dyadic_weighted)
        odd_danger.append(odd_weighted)
        additive_danger.append(additive_weighted)
        multiplicative_danger.append(multiplicative_weighted)
        product_danger.append(product_weighted)

    if criterion.name == "origin_clearance_rank":
        return tuple(dist[0][i] for i in range(n))
    if criterion.name == "local_two_guard":
        return tuple(min_bracket[i] for i in range(n))
    if criterion.name == "incident_danger_sum":
        return tuple(incident_deficit)
    if criterion.name == "incident_close_count":
        return tuple(incident_close)
    if criterion.name == "dyadic_branch_pressure":
        return tuple(dyadic_pressure)
    if criterion.name == "odd_core_branch_count":
        return tuple(odd_core_count)
    if criterion.name == "same_odd_chain_degree":
        return tuple(same_chain)
    if criterion.name == "additive_shadow_degree":
        return tuple(add_deg)
    if criterion.name == "multiplicative_shadow_degree":
        return tuple(mult_deg)
    if criterion.name == "product_sum_interface":
        return tuple(product_interface)
    if criterion.name == "dyadic_danger_curvature":
        return tuple(dyadic_danger)
    if criterion.name == "odd_chain_danger":
        return tuple(odd_danger)
    if criterion.name == "additive_danger_interface":
        return tuple(additive_danger)
    if criterion.name == "multiplicative_danger_interface":
        return tuple(multiplicative_danger)
    if criterion.name == "product_sum_danger":
        return tuple(product_danger)
    raise ValueError(f"unknown local criterion {criterion.name}")


def snapshot(
    criterion: ArcCriterion,
    speeds_with_observer: tuple[int, ...],
    t: Fraction,
) -> RunnerSnapshot:
    if criterion.source != "local":
        by_name = {item.name: item for item in S506.CRITERIA}
        base = by_name[criterion.source]
        snap = S506.snapshot(base, speeds_with_observer, t)
        return RunnerSnapshot(
            criterion=criterion,
            labels=snap.labels,
            out=snap.out,
            strict_out=snap.strict_out,
            ties=snap.ties,
            strict_arcs=snap.strict_arcs,
        )

    labels = tuple(range(len(speeds_with_observer)))
    values = runner_values(criterion, speeds_with_observer, t)
    out = [0] * len(labels)
    strict_out = [0] * len(labels)
    ties = 0
    strict_arcs = 0
    for i, j in combinations(range(len(labels)), 2):
        cmp = compare_values(values[i], values[j])
        if cmp > 0:
            winner, loser = i, j
            strict_arcs += 1
        elif cmp < 0:
            winner, loser = j, i
            strict_arcs += 1
        else:
            winner, loser = i, j
            ties += 1
        out[winner] |= 1 << loser
        if cmp != 0:
            strict_out[winner] |= 1 << loser
    return RunnerSnapshot(
        criterion=criterion,
        labels=labels,
        out=tuple(out),
        strict_out=tuple(strict_out),
        ties=ties,
        strict_arcs=strict_arcs,
        values=values,
    )


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


def h_count(out: tuple[int, ...]) -> int:
    n = len(out)
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for vertex in range(n):
        dp[1 << vertex][vertex] = 1
    for mask in range(1 << n):
        unvisited = full ^ mask
        for last in range(n):
            value = dp[mask][last]
            if not value:
                continue
            nxts = out[last] & unvisited
            while nxts:
                bit = nxts & -nxts
                nxt = bit.bit_length() - 1
                dp[mask | bit][nxt] += value
                nxts ^= bit
    return sum(dp[full])


def fingerprint(snap: RunnerSnapshot, compute_h: bool) -> dict[str, Any]:
    scores = score_list(snap.out)
    strict_scc = scc_sizes(snap.strict_out)
    origin_score = scores[snap.labels.index(0)] if 0 in snap.labels else None
    vertex_count = len(snap.labels)
    return {
        "vertex_count": vertex_count,
        "H": h_count(snap.out) if compute_h and vertex_count <= H_LIMIT else None,
        "score_hist": tuple(sorted(Counter(scores).items())),
        "score_width": max(scores) - min(scores) if scores else 0,
        "score_variance": variance(scores),
        "origin_score": origin_score,
        "top_score": max(scores) if scores else 0,
        "strict_triangles": triangle_count(snap.strict_out),
        "complete_triangles": triangle_count(snap.out),
        "largest_strict_scc": strict_scc[0] if strict_scc else 0,
        "ties": snap.ties,
        "strict_arcs": snap.strict_arcs,
    }


def variance(values: list[int]) -> float:
    if not values:
        return 0.0
    avg = mean(values)
    return mean((value - avg) ** 2 for value in values)


def edge_key(out: tuple[int, ...]) -> tuple[int, ...]:
    bits: list[int] = []
    for i, j in combinations(range(len(out)), 2):
        bits.append(1 if (out[i] >> j) & 1 else 0)
    return tuple(bits)


def edge_flips(left: tuple[int, ...], right: tuple[int, ...]) -> int:
    return sum(1 for a, b in zip(left, right) if a != b)


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


def pair_danger_targets(
    speeds_with_observer: tuple[int, ...], times: tuple[Fraction, ...]
) -> dict[str, list[float]]:
    close_pair_count: list[float] = []
    danger_mass: list[float] = []
    origin_incident_deficit: list[float] = []
    threshold = Fraction(1, len(speeds_with_observer))
    for t in times:
        pos = S506.positions(speeds_with_observer, t)
        dist = S506.distance_matrix(pos)
        close = 0
        mass = Fraction(0)
        origin_mass = Fraction(0)
        for i, j in combinations(range(len(speeds_with_observer)), 2):
            deficit = max(Fraction(0), threshold - dist[i][j])
            if dist[i][j] < threshold:
                close += 1
            mass += deficit
            if i == 0 or j == 0:
                origin_mass += deficit
        close_pair_count.append(float(close))
        danger_mass.append(float(mass / threshold))
        origin_incident_deficit.append(float(origin_mass / threshold))
    return {
        "close_pair_count": close_pair_count,
        "danger_mass": danger_mass,
        "origin_incident_deficit": origin_incident_deficit,
    }


def feature_series(snaps: list[RunnerSnapshot]) -> dict[str, list[float]]:
    out: dict[str, list[float]] = defaultdict(list)
    for snap in snaps:
        fp = fingerprint(snap, compute_h=True)
        n = fp["vertex_count"]
        denom = max(1, n - 1)
        total_pairs = max(1, comb(n, 2))
        total_triples = max(1, comb(n, 3)) if n >= 3 else 1
        h_value = fp["H"] or 0
        out["H"].append(float(h_value))
        out["logH"].append(log2(h_value + 1))
        out["score_width_norm"].append(float(fp["score_width"] / denom))
        out["score_variance_norm"].append(float(fp["score_variance"] / (denom * denom)))
        out["origin_score_norm"].append(
            float((fp["origin_score"] or 0) / denom)
            if fp["origin_score"] is not None
            else 0.0
        )
        out["top_score_norm"].append(float(fp["top_score"] / denom))
        out["tie_rate"].append(float(fp["ties"] / total_pairs))
        out["strict_arc_rate"].append(float(fp["strict_arcs"] / total_pairs))
        out["strict_triangle_density"].append(float(fp["strict_triangles"] / total_triples))
        out["complete_triangle_density"].append(float(fp["complete_triangles"] / total_triples))
        out["strict_scc_norm"].append(float(fp["largest_strict_scc"] / n))
    return dict(out)


def target_series(
    speeds_with_observer: tuple[int, ...],
    times: tuple[Fraction, ...],
    geoms: list[Any],
) -> dict[str, list[float]]:
    out = {
        "origin_clearance": [float(g.stationary_min / g.threshold) for g in geoms],
        "origin_bracket": [float(min(g.origin_left, g.origin_right) / g.threshold) for g in geoms],
        "safe_gap_count": [float(g.safe_gap_count) for g in geoms],
        "lonely_vertices": [float(len(g.lonely_vertices)) for g in geoms],
        "max_gap": [float(g.max_gap / g.threshold) for g in geoms],
    }
    out.update(pair_danger_targets(speeds_with_observer, times))
    return out


def criterion_scorecard() -> dict[str, ScorecardRow]:
    summary: dict[str, ScorecardRow] = {}
    for criterion in CRITERIA:
        best_scores: list[tuple[float, str, str]] = []
        h_origin: list[float] = []
        h_close: list[float] = []
        origin_score_origin: list[float] = []
        top_score_danger: list[float] = []
        flip_rates: list[float] = []
        tie_rates: list[float] = []
        scc_rates: list[float] = []

        for _family_name, speeds in S506.small_families():
            times = S506.exact_cell_times(speeds)
            geoms = [S506.geometry(speeds, t) for t in times]
            snaps = [snapshot(criterion, speeds, t) for t in times]
            features = feature_series(snaps)
            targets = target_series(speeds, times, geoms)
            local_best = (0.0, "-", "-")
            for feature_name, xs in features.items():
                for target_name, ys in targets.items():
                    rho = spearman(xs, ys)
                    if abs(rho) > abs(local_best[0]):
                        local_best = (rho, feature_name, target_name)
            best_scores.append(local_best)
            h_origin.append(spearman(features["logH"], targets["origin_clearance"]))
            h_close.append(spearman(features["logH"], targets["close_pair_count"]))
            origin_score_origin.append(
                spearman(features["origin_score_norm"], targets["origin_clearance"])
            )
            top_score_danger.append(spearman(features["top_score_norm"], targets["danger_mass"]))

            keys = [edge_key(snap.out) for snap in snaps]
            flips = [edge_flips(a, b) for a, b in zip(keys, keys[1:])]
            pair_count = max(1, comb(len(snaps[0].labels), 2))
            if flips:
                flip_rates.append(mean(flips) / pair_count)
            total_pairs = pair_count * len(snaps)
            tie_rates.append(sum(snap.ties for snap in snaps) / total_pairs)
            scc_rates.append(mean(features["strict_scc_norm"]))

        feature_votes = Counter((row[1], row[2]) for row in best_scores)
        best_pair, _ = feature_votes.most_common(1)[0]
        summary[criterion.name] = ScorecardRow(
            criterion=criterion,
            best_abs=mean(abs(row[0]) for row in best_scores),
            best_pair=best_pair,
            h_origin=mean(h_origin),
            h_close=mean(h_close),
            origin_score_origin=mean(origin_score_origin),
            top_score_danger=mean(top_score_danger),
            flip_rate=mean(flip_rates),
            tie_rate=mean(tie_rates),
            scc_rate=mean(scc_rates),
        )
    return summary


def short_score_hist(hist: tuple[tuple[int, int], ...], limit: int = 12) -> str:
    parts = [f"{score}:{count}" for score, count in hist[:limit]]
    if len(hist) > limit:
        parts.append("...")
    return " ".join(parts)


def print_methodology() -> None:
    print("LRC operation-grid arc-criteria Tournament Analysis (S511)")
    print("=" * 104)
    print("Tie Hamiltonian path: increasing runner label 0->1->... .")
    print("Imported S506 criteria keep their native pairwise switch; local criteria")
    print("compare runner scalar values induced by pair-cells and operation labels.")
    print()
    print("Criteria switchboard:")
    for item in CRITERIA:
        dynamic_text = str(item.dynamic)
        print(
            f"  {item.name:<34} layer={item.layer:<23} "
            f"dynamic={dynamic_text}"
        )
        print(f"    observable={item.observable}")
        print(f"    switch={item.switch}")
    print()


def print_scorecard(summary: dict[str, ScorecardRow]) -> None:
    print("SMALL EXACT CLOCK SCORECARD")
    print("=" * 104)
    print(
        "criterion                          best_abs  common_best_feature->target"
        "                  H~origin H~close originScore~origin topScore~danger flips ties scc"
    )
    for name, row in sorted(summary.items(), key=lambda item: -item[1].best_abs):
        best_feature, best_target = row.best_pair
        print(
            f"{name:<34} "
            f"{row.best_abs:>8.3f}  "
            f"{best_feature}->{best_target:<30} "
            f"{row.h_origin:>+8.3f} "
            f"{row.h_close:>+7.3f} "
            f"{row.origin_score_origin:>+10.3f} "
            f"{row.top_score_danger:>+10.3f} "
            f"{row.flip_rate:>5.3f} "
            f"{row.tie_rate:>5.3f} "
            f"{row.scc_rate:>5.3f}"
        )
    print()
    print("Reading: dynamic threshold aggregates should score in the time direction;")
    print("static operation-grid criteria can score only through fixed branch shape,")
    print("so zero Spearman here means coordinate label, not mathematical uselessness.")
    print()


def selected_examples() -> tuple[tuple[str, int, tuple[int, ...]], ...]:
    base = list(S506.selected_examples())
    n14_double_skip, n14_double = S506.S492.best_ladder(14, 28)
    n18_double_skip, n18_double = S506.S492.best_ladder(18, 36)
    base.extend(
        [
            (f"n=14 double-gate ladder skip {n14_double_skip}", 14, n14_double),
            (f"n=18 double-gate ladder skip {n18_double_skip}", 18, n18_double),
        ]
    )
    return tuple(base)


def print_selected_rows() -> None:
    shortlist = {
        "s506_phase_half",
        "s506_lrc_close_sector",
        "incident_danger_sum",
        "dyadic_danger_curvature",
        "additive_danger_interface",
        "multiplicative_danger_interface",
        "product_sum_danger",
        "dyadic_branch_pressure",
        "product_sum_interface",
    }
    criteria = tuple(item for item in CRITERIA if item.name in shortlist)
    print("SELECTED n=14/n=18 RUNNER ROWS")
    print("=" * 104)
    print("H is computed through N=14 and omitted for N=18.  Large-row fingerprints")
    print("should be read as score shape, observer score, ties, triangles, and SCCs.")
    print()
    for label, n, moving in selected_examples():
        speeds = (0,) + moving
        report = S506.S356.report(label, list(moving))
        summary = S506.S360.summarize(list(moving))
        print(f"[{label}]")
        print(
            f"  class={summary.classification} gap/th={fmt(report.max_gap / report.threshold)} "
            f"unprotected={summary.unprotected_count} speeds={moving}"
        )
        for tag, t in S506.selected_times(moving, n):
            geom = S506.geometry(speeds, t)
            pair_targets = pair_danger_targets(speeds, (t,))
            print(
                f"  time={tag:<9} t={fmt(t):>14} "
                f"origin_clear={fmt(geom.stationary_min / geom.threshold):>7} "
                f"bracket={fmt(min(geom.origin_left, geom.origin_right) / geom.threshold):>7} "
                f"safe={geom.safe_gap_count:>2} lonely={len(geom.lonely_vertices):>2} "
                f"max_gap/th={fmt(geom.max_gap / geom.threshold):>7} "
                f"close_pairs={int(pair_targets['close_pair_count'][0]):>3} "
                f"danger/th={pair_targets['danger_mass'][0]:>6.2f}"
            )
            for criterion in criteria:
                snap = snapshot(criterion, speeds, t)
                fp = fingerprint(snap, compute_h=True)
                total_pairs = max(1, comb(fp["vertex_count"], 2))
                print(
                    f"    {criterion.name:<34} "
                    f"v={fp['vertex_count']:>2} H={fmt(fp['H']):>14} "
                    f"score=[{short_score_hist(fp['score_hist'])}] "
                    f"w={fp['score_width']:>2} o={fmt(fp['origin_score']):>3} "
                    f"strict={fp['strict_arcs']:>3} ties={fp['ties']:>3} "
                    f"tri={fp['strict_triangles']:>3} scc={fp['largest_strict_scc']:>2} "
                    f"tie%={fp['ties'] / total_pairs:.3f}"
                )
        print()


def operation_grid_lines(limit: int = 48) -> list[str]:
    lines: list[str] = []
    for odd in range(1, min(limit, 24), 2):
        chain: list[int] = []
        value = odd
        while value <= limit:
            chain.append(value)
            value *= 2
        lines.append(f"  odd core {odd:>2}: {' -> '.join(map(str, chain))}")
    return lines


def print_operation_dictionary() -> None:
    print("A000568, x+2/x*2, AND OPERATION INTERFACES")
    print("=" * 104)
    print("A000568 terms from the odd-partition Burnside formula:")
    print("  " + ", ".join(str(S509.a000568_burnside(n)) for n in range(1, 11)))
    print("Only odd cycle partitions contribute; this is the quotient-level analogue")
    print("of an odd-core survival rule for complete orientations.")
    print()
    print("2-adic operation grid up to 48:")
    print("  horizontal move: odd core x -> x+2")
    print("  vertical move:   branch height x -> 2x")
    for line in operation_grid_lines():
        print(line)
    print()
    print("How the script turns operations into arc criteria:")
    print("  addition: speed equations a+b=c or |a-b|=c give additive-shadow degree.")
    print("  multiplication: divisibility and doubling branches give x*2 pressure.")
    print("  product-sum: (a-1)(b-1)=N-1 gives a critical interface score.")
    print("  hybrid: multiply each static operation label by current 1/N danger.")
    print()
    print("Two-factor product-sum landmarks:")
    for k in (14, 18, 21, 24, 30):
        sols = S509.S365.two_factor_solutions_for_arity(k)
        text = ", ".join(f"({a},{b})->{p}" for a, b, p in sols) or "-"
        seeds = S509.S365.enumerate_product_sum_seeds(max_k=k, max_product=4 * k).get(k, [])
        best = S509.S365.compact_witness(seeds[0]) if seeds else "-"
        print(f"  N={k:>2}: two-factor {text}; seed {best}")
    print()


def profile_for(row: ScorecardRow) -> CriteriaProfile:
    name = row.criterion.name
    signal = min(10, max(0, round(10 * row.best_abs)))
    threshold_words = ("danger", "threshold", "origin", "close", "safe", "guard", "lrc")
    operation_words = ("dyadic", "odd", "additive", "multiplicative", "product", "branch")
    threshold = 5 if any(word in name for word in threshold_words) else 2
    operation_grid = 5 if any(word in name for word in operation_words) else 1
    dynamic_motion = 5 if row.criterion.dynamic else 1
    novelty = 5 if row.criterion.layer == "hybrid" else 4 if operation_grid >= 5 else 2
    projection_risk = 5
    if row.criterion.layer == "hybrid":
        projection_risk = 1
    elif row.criterion.dynamic and threshold >= 5:
        projection_risk = 2
    elif row.criterion.dynamic:
        projection_risk = 3
    elif operation_grid >= 5:
        projection_risk = 4
    score = signal + 2 * threshold + 2 * operation_grid + dynamic_motion + novelty - 2 * projection_risk
    return CriteriaProfile(
        name=name,
        signal=signal,
        threshold=threshold,
        operation_grid=operation_grid,
        dynamic_motion=dynamic_motion,
        novelty=novelty,
        projection_risk=projection_risk,
        score=score,
    )


def criteria_tournament(profiles: list[CriteriaProfile]) -> tuple[tuple[int, ...], list[int]]:
    n = len(profiles)
    out = [0] * n
    for i, j in combinations(range(n), 2):
        votes_i = 0
        votes_j = 0
        for attr in ("signal", "threshold", "operation_grid", "dynamic_motion", "novelty"):
            ai = getattr(profiles[i], attr)
            aj = getattr(profiles[j], attr)
            if ai > aj:
                votes_i += 1
            elif aj > ai:
                votes_j += 1
        if profiles[i].projection_risk < profiles[j].projection_risk:
            votes_i += 1
        elif profiles[j].projection_risk < profiles[i].projection_risk:
            votes_j += 1
        winner = i if votes_i >= votes_j else j
        loser = j if winner == i else i
        out[winner] |= 1 << loser
    return tuple(out), score_list(tuple(out))


def print_criteria_tournament(summary: dict[str, ScorecardRow]) -> None:
    profiles = [profile_for(row) for row in summary.values()]
    out, scores = criteria_tournament(profiles)
    print("CRITERIA-CHOICE TOURNAMENT")
    print("=" * 104)
    print("Vertices are arc criteria. Observable is a six-attribute profile:")
    print("signal, threshold fidelity, operation-grid content, dynamic motion, novelty,")
    print("and lower projection risk. Switch is majority vote; tie path is list order.")
    print(
        f"  H={h_count(out)} score_hist={dict(sorted(Counter(scores).items()))} "
        f"c3={triangle_count(out)} largest_scc={scc_sizes(out)[0]}"
    )
    ranked = sorted(zip(profiles, scores), key=lambda pair: (-pair[0].score, -pair[1], pair[0].name))
    print("  top criteria by profile:")
    for profile, tournament_score in ranked[:10]:
        print(
            f"    {profile.name:<34} profile={profile.score:>2} "
            f"t_score={tournament_score:>2} sig={profile.signal} "
            f"thr={profile.threshold} op={profile.operation_grid} "
            f"dyn={profile.dynamic_motion} risk={profile.projection_risk}"
        )
    print()


def print_synthesis(summary: dict[str, ScorecardRow]) -> None:
    print("SYNTHESIS")
    print("=" * 104)
    ordered = sorted(summary.values(), key=lambda row: -row.best_abs)
    print("Most useful arc criteria in this audit:")
    for row in ordered[:7]:
        feature, target = row.best_pair
        print(
            f"  {row.criterion.name}: {feature}->{target}, "
            f"mean_abs_spearman={row.best_abs:.3f}, "
            f"ties={row.tie_rate:.3f}, scc={row.scc_rate:.3f}"
        )
    print()
    print("Interpretation:")
    print("  1. Static x+2/x*2 labels are not time-loneliness meters by themselves.")
    print("     Their job is to name the operation-grid branch of each obstruction.")
    print("  2. Once those labels weight current pair danger, they become moving")
    print("     loneliness coordinates: dyadic, additive, multiplicative, and")
    print("     product-sum danger are all different projections of the same fiber.")
    print("  3. A000568 supplies the quotient/base analogy: complete orientations are")
    print("     classified after forgetting labels, while LRC needs the labelled")
    print("     operation-grid and endpoint-pressure fiber above that base.")
    print("  4. Addition behaves like horizontal transport among odd-core cells;")
    print("     multiplication behaves like vertical branch refinement; product-sum")
    print("     equations are the interfaces where arity/addition is traded for")
    print("     multiplicative factor structure.")
    print()
    print("Candidate loneliness-related metric vector:")
    print("  phase_H, lrc_close_tie_rate, incident_danger_origin_score,")
    print("  dyadic_danger_score_hist, additive_danger_score_hist,")
    print("  multiplicative_danger_score_hist, product_sum_danger_score_hist,")
    print("  strict pressure SCCs, endpoint-labelled origin bracket.")
    print()
    print("Conjectural next move: project hard rows to the A000568 base, but carry this")
    print("operation-weighted danger vector as fiber data; mixed fibers are repair")
    print("opportunities and all-forbidden fibers are the counterexample-shaped object.")


def main() -> None:
    print_methodology()
    summary = criterion_scorecard()
    print_scorecard(summary)
    print_selected_rows()
    print_operation_dictionary()
    print_criteria_tournament(summary)
    print_synthesis(summary)


if __name__ == "__main__":
    main()
