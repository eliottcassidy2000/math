#!/usr/bin/env python3
"""Audit a three-layer LRC add/multiply stack on hard rows.

codex-2026-06-01 S514

This script recovers and extends the S511-S513 line.  S513 predicted that a
counterexample-shaped LRC row would have to align several bad coordinates at
once:

* dynamic runner/endpoint danger,
* pair-cell operation-weighted danger,
* denominator add/multiply/A000568 labels.

Tournament Analysis declaration:

Row tournaments
    vertices: selected LRC rows, represented by exact sampled corridor times.
    pairwise observables:
        runner layer: LRC gap, close-pair mass, threshold danger, safe-gap
            shortage, and pressure SCC size;
        operation layer: dyadic/additive/multiplicative/product-sum danger
            mass induced by current pair-cell deficits;
        denominator layer: phi(N), Goldbach/Lemoine scarcity, dyadic height,
            product-sum two-factor collisions, odd A000568 survival, and the
            S513 scalar grid-stack score.
    switches/gauges: a row points to another row when it is more
        obstruction-shaped by majority vote in the chosen layer.
    tie Hamiltonian path: the listed row order below.

State tournaments
    vertices: selected row-time states.
    observable/switch: the full-stack majority vote on state-level features.
    tie Hamiltonian path: row order, then selected-time order.

Stored output:
    05-knowledge/results/lrc_three_layer_stack_audit_s514.out
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from fractions import Fraction
from importlib.machinery import SourceFileLoader
from itertools import combinations
from math import comb, log2
from pathlib import Path
from typing import Callable


ROOT = Path(__file__).resolve().parents[1]
S506 = SourceFileLoader(
    "lrc_arc_criteria_loneliness_s506",
    str(ROOT / "04-computation" / "lrc_arc_criteria_loneliness_s506.py"),
).load_module()
S511 = SourceFileLoader(
    "lrc_operation_grid_arc_criteria_s511",
    str(ROOT / "04-computation" / "lrc_operation_grid_arc_criteria_s511.py"),
).load_module()
S513 = SourceFileLoader(
    "lrc_add_mult_gauge_stack_s513",
    str(ROOT / "04-computation" / "lrc_add_mult_gauge_stack_s513.py"),
).load_module()


H_LIMIT = 16


def fmt_frac(value: Fraction | int | None) -> str:
    if value is None:
        return "-"
    if isinstance(value, int):
        return str(value)
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def fmt_float(value: float | Fraction | None, places: int = 3) -> str:
    if value is None:
        return "-"
    return f"{float(value):.{places}f}"


def initial(n: int) -> tuple[int, ...]:
    return tuple(range(1, n))


def selected_rows() -> tuple[tuple[str, int, tuple[int, ...]], ...]:
    n14_lpd_skip, n14_lpd = S506.S492.best_ladder(14, S506.S492.largest_proper_divisor(14))
    n18_lpd_skip, n18_lpd = S506.S492.best_ladder(18, S506.S492.largest_proper_divisor(18))
    n14_gate_skip, n14_gate = S506.S492.best_ladder(14, 14)
    n18_gate_skip, n18_gate = S506.S492.best_ladder(18, 18)
    n14_double_skip, n14_double = S506.S492.best_ladder(14, 28)
    n18_double_skip, n18_double = S506.S492.best_ladder(18, 36)
    n16_lpd_skip, n16_lpd = S506.S492.best_ladder(16, S506.S492.largest_proper_divisor(16))
    n16_gate_skip, n16_gate = S506.S492.best_ladder(16, 16)
    return (
        ("n14 initial", 14, initial(14)),
        (f"n14 row-parent {S506.S492.largest_proper_divisor(14)}* skip {n14_lpd_skip}", 14, n14_lpd),
        (f"n14 gate 14* skip {n14_gate_skip}", 14, n14_gate),
        (f"n14 double-gate 28* skip {n14_double_skip}", 14, n14_double),
        ("n16 initial", 16, initial(16)),
        (f"n16 row-parent {S506.S492.largest_proper_divisor(16)}* skip {n16_lpd_skip}", 16, n16_lpd),
        (f"n16 gate 16* skip {n16_gate_skip}", 16, n16_gate),
        ("n18 initial", 18, initial(18)),
        (f"n18 row-parent {S506.S492.largest_proper_divisor(18)}* skip {n18_lpd_skip}", 18, n18_lpd),
        (f"n18 gate 18* skip {n18_gate_skip}", 18, n18_gate),
        (f"n18 double-gate 36* skip {n18_double_skip}", 18, n18_double),
    )


LOCAL_BY_NAME = {criterion.name: criterion for criterion in S511.LOCAL_CRITERIA}
CRITERIA_BY_NAME = {criterion.name: criterion for criterion in S511.CRITERIA}


@dataclass(frozen=True)
class StateFeature:
    row_index: int
    label: str
    n: int
    speeds: tuple[int, ...]
    time_tag: str
    time: Fraction
    origin_ratio: Fraction
    bracket_ratio: Fraction
    max_gap_ratio: Fraction
    safe_gaps: int
    lonely_vertices: int
    close_pair_count: int
    danger_density: float
    origin_danger_density: float
    pressure_scc_norm: float
    pressure_tri_density: float
    phase_log_h: float | None
    add_danger_density: float
    mult_danger_density: float
    product_danger_density: float
    dyadic_danger_density: float


@dataclass(frozen=True)
class RowFeature:
    index: int
    label: str
    n: int
    speeds: tuple[int, ...]
    gap_ratio: Fraction
    endpoint_debt: int
    add_count: int
    dyadic_height: int
    odd_core: int
    twofactor_count: int
    odd_survival: Fraction
    scalar_stack: int
    max_close_fraction: float
    max_danger_density: float
    max_origin_danger_density: float
    min_safe_fraction: float
    max_pressure_scc: float
    max_pressure_tri_density: float
    max_add_danger: float
    max_mult_danger: float
    max_product_danger: float
    max_dyadic_danger: float
    best_origin_ratio: Fraction
    states: tuple[StateFeature, ...]


@dataclass(frozen=True)
class Dimension:
    name: str
    getter: Callable[[RowFeature], float]
    high_wins: bool
    note: str


@dataclass(frozen=True)
class StateDimension:
    name: str
    getter: Callable[[StateFeature], float]
    high_wins: bool
    note: str


def operation_density(
    speeds: tuple[int, ...],
    time: Fraction,
    criterion_name: str,
) -> float:
    criterion = LOCAL_BY_NAME[criterion_name]
    values = S511.runner_values(criterion, speeds, time)
    threshold = Fraction(1, len(speeds))
    pair_count = max(1, comb(len(speeds), 2))
    return float(sum(values, Fraction(0)) / 2 / threshold / pair_count)


def state_feature(row_index: int, label: str, n: int, moving: tuple[int, ...], tag: str, time: Fraction) -> StateFeature:
    speeds = (0,) + moving
    geom = S506.geometry(speeds, time)
    pair_targets = S511.pair_danger_targets(speeds, (time,))
    pressure_snap = S511.snapshot(CRITERIA_BY_NAME["s506_threshold_deficit_pressure"], speeds, time)
    pressure_fp = S511.fingerprint(pressure_snap, compute_h=False)
    phase_snap = S511.snapshot(CRITERIA_BY_NAME["s506_phase_half"], speeds, time)
    phase_fp = S511.fingerprint(phase_snap, compute_h=True)
    triples = max(1, comb(len(speeds), 3))
    pair_count = max(1, comb(len(speeds), 2))
    phase_h = phase_fp["H"]
    return StateFeature(
        row_index=row_index,
        label=label,
        n=n,
        speeds=speeds,
        time_tag=tag,
        time=time,
        origin_ratio=geom.stationary_min / geom.threshold,
        bracket_ratio=min(geom.origin_left, geom.origin_right) / geom.threshold,
        max_gap_ratio=geom.max_gap / geom.threshold,
        safe_gaps=geom.safe_gap_count,
        lonely_vertices=len(geom.lonely_vertices),
        close_pair_count=int(pair_targets["close_pair_count"][0]),
        danger_density=pair_targets["danger_mass"][0] / pair_count,
        origin_danger_density=pair_targets["origin_incident_deficit"][0] / max(1, n - 1),
        pressure_scc_norm=pressure_fp["largest_strict_scc"] / len(speeds),
        pressure_tri_density=pressure_fp["strict_triangles"] / triples,
        phase_log_h=None if phase_h is None else log2(phase_h + 1),
        add_danger_density=operation_density(speeds, time, "additive_danger_interface"),
        mult_danger_density=operation_density(speeds, time, "multiplicative_danger_interface"),
        product_danger_density=operation_density(speeds, time, "product_sum_danger"),
        dyadic_danger_density=operation_density(speeds, time, "dyadic_danger_curvature"),
    )


def row_feature(index: int, label: str, n: int, moving: tuple[int, ...]) -> RowFeature:
    report = S506.S356.report(label, list(moving))
    record = S513.denominator_record(n)
    states = tuple(
        state_feature(index, label, n, moving, tag, time)
        for tag, time in S506.selected_times(moving, n)
    )
    pair_count = max(1, comb(n, 2))
    return RowFeature(
        index=index,
        label=label,
        n=n,
        speeds=moving,
        gap_ratio=report.max_gap / report.threshold,
        endpoint_debt=record.phi,
        add_count=record.add_count,
        dyadic_height=record.dyadic_height,
        odd_core=record.odd_core,
        twofactor_count=record.twofactor_count,
        odd_survival=record.odd_partition_fraction,
        scalar_stack=S513.scalar_grid_stack_score(record),
        max_close_fraction=max(state.close_pair_count / pair_count for state in states),
        max_danger_density=max(state.danger_density for state in states),
        max_origin_danger_density=max(state.origin_danger_density for state in states),
        min_safe_fraction=min(state.safe_gaps / n for state in states),
        max_pressure_scc=max(state.pressure_scc_norm for state in states),
        max_pressure_tri_density=max(state.pressure_tri_density for state in states),
        max_add_danger=max(state.add_danger_density for state in states),
        max_mult_danger=max(state.mult_danger_density for state in states),
        max_product_danger=max(state.product_danger_density for state in states),
        max_dyadic_danger=max(state.dyadic_danger_density for state in states),
        best_origin_ratio=max(state.origin_ratio for state in states),
        states=states,
    )


RUNNER_DIMS = (
    Dimension("small_LRC_gap", lambda r: float(r.gap_ratio), False, "smaller certified gap wins"),
    Dimension("close_pair_mass", lambda r: r.max_close_fraction, True, "more close pair-cells wins"),
    Dimension("threshold_danger", lambda r: r.max_danger_density, True, "more threshold deficit wins"),
    Dimension("safe_gap_shortage", lambda r: r.min_safe_fraction, False, "fewer safe gaps wins"),
    Dimension("pressure_SCC", lambda r: r.max_pressure_scc, True, "larger pressure SCC wins"),
    Dimension("pressure_triangles", lambda r: r.max_pressure_tri_density, True, "more pressure cycles wins"),
)

OPERATION_DIMS = (
    Dimension("dyadic_danger", lambda r: r.max_dyadic_danger, True, "more dyadic-weighted danger wins"),
    Dimension("additive_danger", lambda r: r.max_add_danger, True, "more additive-weighted danger wins"),
    Dimension("multiplicative_danger", lambda r: r.max_mult_danger, True, "more multiplicative-weighted danger wins"),
    Dimension("product_sum_danger", lambda r: r.max_product_danger, True, "more product-sum danger wins"),
    Dimension("scalar_stack", lambda r: r.scalar_stack, True, "larger S513 stack score wins"),
)

DENOM_DIMS = (
    Dimension("endpoint_debt_phi", lambda r: r.endpoint_debt, True, "larger phi(N) wins"),
    Dimension("additive_scarcity", lambda r: r.add_count, False, "fewer prime representations wins"),
    Dimension("dyadic_height", lambda r: r.dyadic_height, True, "higher x*2 branch wins"),
    Dimension("product_sum_collisions", lambda r: r.twofactor_count, True, "more two-factor product-sum gates wins"),
    Dimension("odd_A000568_survival", lambda r: float(r.odd_survival), True, "larger odd partition survival wins"),
    Dimension("small_best_origin", lambda r: float(r.best_origin_ratio), False, "smaller best witness margin wins"),
)

FULL_DIMS = RUNNER_DIMS + OPERATION_DIMS + DENOM_DIMS

STATE_FULL_DIMS = (
    StateDimension("small_origin_ratio", lambda s: float(s.origin_ratio), False, "closer observer wins"),
    StateDimension("small_bracket", lambda s: float(s.bracket_ratio), False, "smaller two-sided bracket wins"),
    StateDimension("close_pairs", lambda s: s.close_pair_count, True, "more close pairs wins"),
    StateDimension("danger_density", lambda s: s.danger_density, True, "more pair danger wins"),
    StateDimension("origin_danger", lambda s: s.origin_danger_density, True, "more observer incident danger wins"),
    StateDimension("safe_gaps", lambda s: s.safe_gaps, False, "fewer safe gaps wins"),
    StateDimension("pressure_SCC", lambda s: s.pressure_scc_norm, True, "larger pressure SCC wins"),
    StateDimension("additive_danger", lambda s: s.add_danger_density, True, "more additive danger wins"),
    StateDimension("multiplicative_danger", lambda s: s.mult_danger_density, True, "more multiplicative danger wins"),
    StateDimension("product_sum_danger", lambda s: s.product_danger_density, True, "more product-sum danger wins"),
)


def compare(left: float, right: float) -> int:
    return (left > right) - (left < right)


def row_tournament(rows: tuple[RowFeature, ...], dims: tuple[Dimension, ...]) -> list[list[int]]:
    n = len(rows)
    adj = [[0] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        votes_i = 0
        votes_j = 0
        for dim in dims:
            cmp = compare(dim.getter(rows[i]), dim.getter(rows[j]))
            if cmp == 0:
                continue
            if dim.high_wins:
                votes_i += int(cmp > 0)
                votes_j += int(cmp < 0)
            else:
                votes_i += int(cmp < 0)
                votes_j += int(cmp > 0)
        winner = i if votes_i >= votes_j else j
        loser = j if winner == i else i
        adj[winner][loser] = 1
    return adj


def state_tournament(states: tuple[StateFeature, ...]) -> list[list[int]]:
    n = len(states)
    adj = [[0] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        votes_i = 0
        votes_j = 0
        for dim in STATE_FULL_DIMS:
            cmp = compare(dim.getter(states[i]), dim.getter(states[j]))
            if cmp == 0:
                continue
            if dim.high_wins:
                votes_i += int(cmp > 0)
                votes_j += int(cmp < 0)
            else:
                votes_i += int(cmp < 0)
                votes_j += int(cmp > 0)
        winner = i if votes_i >= votes_j else j
        loser = j if winner == i else i
        adj[winner][loser] = 1
    return adj


def scores(adj: list[list[int]]) -> tuple[int, ...]:
    return tuple(sum(row) for row in adj)


def score_hist(adj: list[list[int]]) -> str:
    hist = Counter(scores(adj))
    return ",".join(f"{score}:{hist[score]}" for score in sorted(hist))


def directed_triangles(adj: list[list[int]]) -> int:
    total = 0
    for a, b, c in combinations(range(len(adj)), 3):
        if (adj[a][b] and adj[b][c] and adj[c][a]) or (
            adj[a][c] and adj[c][b] and adj[b][a]
        ):
            total += 1
    return total


def scc_sizes(adj: list[list[int]]) -> tuple[int, ...]:
    n = len(adj)
    graph = [[j for j in range(n) if adj[i][j]] for i in range(n)]
    reverse = [[i for i in range(n) if adj[i][j]] for j in range(n)]
    seen = [False] * n
    order: list[int] = []

    def dfs(vertex: int) -> None:
        seen[vertex] = True
        for nxt in graph[vertex]:
            if not seen[nxt]:
                dfs(nxt)
        order.append(vertex)

    for vertex in range(n):
        if not seen[vertex]:
            dfs(vertex)

    seen = [False] * n
    sizes: list[int] = []
    for start in reversed(order):
        if seen[start]:
            continue
        todo = deque([start])
        seen[start] = True
        size = 0
        while todo:
            vertex = todo.pop()
            size += 1
            for nxt in reverse[vertex]:
                if not seen[nxt]:
                    seen[nxt] = True
                    todo.append(nxt)
        sizes.append(size)
    return tuple(sorted(sizes, reverse=True))


def hamiltonian_paths(adj: list[list[int]]) -> int | None:
    n = len(adj)
    if n > H_LIMIT:
        return None
    dp: dict[tuple[int, int], int] = defaultdict(int)
    for vertex in range(n):
        dp[(1 << vertex, vertex)] = 1
    for mask in range(1, 1 << n):
        for vertex in range(n):
            ways = dp.get((mask, vertex), 0)
            if not ways:
                continue
            for nxt in range(n):
                if not (mask & (1 << nxt)) and adj[vertex][nxt]:
                    dp[(mask | (1 << nxt), nxt)] += ways
    full = (1 << n) - 1
    return sum(dp.get((full, vertex), 0) for vertex in range(n))


def edge_flip_rate(left: list[list[int]], right: list[list[int]]) -> Fraction:
    flips = 0
    total = 0
    for i, j in combinations(range(len(left)), 2):
        total += 1
        flips += int(left[i][j] != right[i][j])
    return Fraction(flips, total)


def top_rows(rows: tuple[RowFeature, ...], adj: list[list[int]], limit: int = 5) -> str:
    ranked = sorted(((scores(adj)[idx], rows[idx].label) for idx in range(len(rows))), reverse=True)
    return " | ".join(f"{label}:{score}" for score, label in ranked[:limit])


def top_states(states: tuple[StateFeature, ...], adj: list[list[int]], limit: int = 8) -> str:
    ranked = sorted(
        ((scores(adj)[idx], -idx, states[idx]) for idx in range(len(states))),
        reverse=True,
    )
    return " | ".join(f"{state.label}/{state.time_tag}:{score}" for score, _idx, state in ranked[:limit])


def print_methodology() -> None:
    print("LRC THREE-LAYER STACK AUDIT - S514")
    print("=" * 92)
    print("Recovered line: S511 operation-weighted danger, S512 decorated fibers,")
    print("S513 add/multiply gauge stack.  This pass asks whether selected hard")
    print("rows align the whole obstruction stack at once.")
    print()
    print("Tie Hamiltonian paths:")
    print("  row tournaments: listed row order")
    print("  state tournament: row order, then selected-time order")
    print()
    print("Layer switches:")
    for title, dims in (
        ("runner dynamic", RUNNER_DIMS),
        ("operation-weighted", OPERATION_DIMS),
        ("denominator static", DENOM_DIMS),
    ):
        print(f"  {title}:")
        for dim in dims:
            direction = "high" if dim.high_wins else "low"
            print(f"    {dim.name:<24} {direction:<4} wins  {dim.note}")
    print()


def print_rows(rows: tuple[RowFeature, ...]) -> None:
    print("ROW STACK FEATURES")
    print("=" * 92)
    print(
        "row                                  N gap/th bestO phi add h core ps oddsurv "
        "close danger safeSCC addD multD prodD"
    )
    for row in rows:
        print(
            f"{row.label:<36} {row.n:>2} "
            f"{fmt_frac(row.gap_ratio):>7} "
            f"{fmt_frac(row.best_origin_ratio):>5} "
            f"{row.endpoint_debt:>3} "
            f"{row.add_count:>3} "
            f"{row.dyadic_height:>1} "
            f"{row.odd_core:>4} "
            f"{row.twofactor_count:>2} "
            f"{fmt_float(row.odd_survival):>7} "
            f"{row.max_close_fraction:>5.3f} "
            f"{row.max_danger_density:>6.3f} "
            f"{row.max_pressure_scc:>6.3f} "
            f"{row.max_add_danger:>5.3f} "
            f"{row.max_mult_danger:>5.3f} "
            f"{row.max_product_danger:>5.3f}"
        )
    print()


def print_row_tournaments(rows: tuple[RowFeature, ...]) -> dict[str, list[list[int]]]:
    tournaments = {
        "runner_dynamic": row_tournament(rows, RUNNER_DIMS),
        "operation_weighted": row_tournament(rows, OPERATION_DIMS),
        "denominator_static": row_tournament(rows, DENOM_DIMS),
        "full_stack": row_tournament(rows, FULL_DIMS),
    }
    print("ROW TOURNAMENT FINGERPRINTS")
    print("=" * 92)
    print("gauge                  H       c3  SCCs                 score_hist              top obstruction-shaped rows")
    for name, adj in tournaments.items():
        h_value = hamiltonian_paths(adj)
        print(
            f"{name:<22} {fmt_frac(h_value):>8} "
            f"{directed_triangles(adj):>4} "
            f"{str(scc_sizes(adj)):<20} "
            f"{score_hist(adj):<22} "
            f"{top_rows(rows, adj)}"
        )
    print()
    print("ROW LAYER EDGE-FLIP RATES")
    print("=" * 92)
    names = tuple(tournaments)
    print(" " * 24 + " ".join(f"{name[:8]:>8}" for name in names))
    for left in names:
        cells = []
        for right in names:
            if left == right:
                cells.append("      --")
            else:
                cells.append(f"{float(edge_flip_rate(tournaments[left], tournaments[right])):>8.2f}")
        print(f"{left:<24}" + " ".join(cells))
    print()
    return tournaments


def print_state_tournament(states: tuple[StateFeature, ...]) -> None:
    selected = tuple(
        state
        for state in states
        if state.time_tag in {"unit", "half-unit", "gap-mid", "boundary"}
    )
    adj = state_tournament(selected)
    print("STATE TOURNAMENT FINGERPRINT")
    print("=" * 92)
    print(
        f"states={len(selected)} H={fmt_frac(hamiltonian_paths(adj))} "
        f"c3={directed_triangles(adj)} SCCs={scc_sizes(adj)} score_hist={score_hist(adj)}"
    )
    print(f"top states: {top_states(selected, adj)}")
    print()
    print("Selected top-state diagnostics:")
    ranked = sorted(
        ((scores(adj)[idx], -idx, selected[idx]) for idx in range(len(selected))),
        reverse=True,
    )
    for score, _idx, state in ranked[:10]:
        print(
            f"  {state.label:<36} {state.time_tag:<9} score={score:>2} "
            f"t={fmt_frac(state.time):>14} "
            f"origin/th={fmt_frac(state.origin_ratio):>7} "
            f"close={state.close_pair_count:>3} "
            f"danger={state.danger_density:>5.3f} "
            f"safe={state.safe_gaps:>2} "
            f"opD=({state.add_danger_density:.3f},"
            f"{state.mult_danger_density:.3f},"
            f"{state.product_danger_density:.3f})"
        )
    print()


def median(values: list[float]) -> float:
    ordered = sorted(values)
    mid = len(ordered) // 2
    if len(ordered) % 2:
        return ordered[mid]
    return (ordered[mid - 1] + ordered[mid]) / 2


def print_conjunction_audit(rows: tuple[RowFeature, ...]) -> None:
    add_med = median([row.add_count for row in rows])
    phi_med = median([row.endpoint_debt for row in rows])
    ps_med = median([row.twofactor_count for row in rows])
    odd_med = median([float(row.odd_survival) for row in rows])
    pressure_cut = 1.0 / min(row.n for row in rows)
    print("COUNTEREXAMPLE-SHAPED CONJUNCTION AUDIT")
    print("=" * 92)
    print(
        "Flags: scarce(add<=median), debt(phi>=median), ps(twofactor>=median), "
        "A568(odd survival>=median), pressure(nontrivial SCC)."
    )
    print(
        f"medians: add={add_med:.1f}, phi={phi_med:.1f}, "
        f"twofactor={ps_med:.1f}, oddsurv={odd_med:.3f}, pressure>{pressure_cut:.3f}"
    )
    print("row                                  scarce debt ps A568 pressure  flags")
    best_flags: list[tuple[int, str]] = []
    for row in rows:
        flags = {
            "scarce": row.add_count <= add_med,
            "debt": row.endpoint_debt >= phi_med,
            "ps": row.twofactor_count >= ps_med,
            "A568": float(row.odd_survival) >= odd_med,
            "pressure": row.max_pressure_scc > pressure_cut,
        }
        count = sum(flags.values())
        best_flags.append((count, row.label))
        print(
            f"{row.label:<36} "
            f"{'Y' if flags['scarce'] else '-':>7} "
            f"{'Y' if flags['debt'] else '-':>4} "
            f"{'Y' if flags['ps'] else '-':>2} "
            f"{'Y' if flags['A568'] else '-':>4} "
            f"{'Y' if flags['pressure'] else '-':>8} "
            f"{count}/5"
        )
    top = sorted(best_flags, reverse=True)[:4]
    print()
    print("Closest sampled rows: " + " | ".join(f"{label}:{count}/5" for count, label in top))
    print()


def print_synthesis(rows: tuple[RowFeature, ...], tournaments: dict[str, list[list[int]]]) -> None:
    full = tournaments["full_stack"]
    runner = tournaments["runner_dynamic"]
    denom = tournaments["denominator_static"]
    op = tournaments["operation_weighted"]
    print("SYNTHESIS")
    print("=" * 92)
    print(
        "1. The full-stack row tournament is not a scalar ledger: it has "
        f"c3={directed_triangles(full)}, SCCs={scc_sizes(full)}, and H={fmt_frac(hamiltonian_paths(full))}."
    )
    print(
        "2. Runner dynamics and denominator labels disagree on "
        f"{fmt_float(edge_flip_rate(runner, denom), 2)} of row edges; operation-weighted danger sits "
        f"{fmt_float(edge_flip_rate(op, runner), 2)} from runner dynamics and "
        f"{fmt_float(edge_flip_rate(op, denom), 2)} from denominator labels."
    )
    pressure_rows = [row.label for row in rows if row.max_pressure_scc > 1 / row.n]
    if pressure_rows:
        print(
            "3. Nontrivial pressure SCCs appear in sampled rows: "
            + ", ".join(pressure_rows)
            + "."
        )
    else:
        print(
            "3. No sampled row has a nontrivial pressure SCC at the selected times. "
            "The S513 counterexample stack currently fails at the pressure-peeling coordinate."
        )
    print(
        "4. Practical next test: keep this row-level stack but replace the pressure coordinate "
        "by owner-compatible endpoint cores from THM-380, because the coarse row-time SCC "
        "coordinate is still too easy to peel."
    )


def main() -> None:
    rows = tuple(row_feature(index, label, n, moving) for index, (label, n, moving) in enumerate(selected_rows()))
    states = tuple(state for row in rows for state in row.states)
    print_methodology()
    print_rows(rows)
    tournaments = print_row_tournaments(rows)
    print_state_tournament(states)
    print_conjunction_audit(rows)
    print_synthesis(rows, tournaments)


if __name__ == "__main__":
    main()
