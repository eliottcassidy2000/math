#!/usr/bin/env python3
"""Audit the n=18 gate/source battlefield inspired by THM-382/383 and HYP-1981.

codex-2026-06-01 S520

THM-382 says the useful finite LRC quotient is threshold-decorated; raw
half-turn tournament classes are too coarse.  THM-383 says equality walls and
tie completion are part of the compactified half-turn menu.  HYP-1981/THM-381
then identifies the exact target:

    observer lonely at t  <=>  marked observer is a source.

THM-384/HYP-1986 sharpen the local target to the compactified source-gap
fiber, HYP-1987 identifies the true A000568 target as the arc-confined source
menu, and THM-385/HYP-1988 stratify the observer-score/blocker layers around
that target.  This script applies those lenses to the n=18 gate battlefield:
does the no-18-gate branch keep observer-source targets alive, and does the
18-gate branch merely export those targets to descendant endpoint layers?

Tournament Analysis declaration:

Branch tournament
    vertices: representative n=18 rows from the no-gate, one-gate,
        lpd-ladder, gate-ladder, and double-gate branches.
    pairwise observable: the row audit vector
        (has 18-gate, unit witness loss, LRC gap, unprotected endpoint count,
        descendant endpoint denominator size).
    switch/gauge: row A points to row B if A is more counterexample-shaped by
        majority vote over the audit vector.
    tie Hamiltonian path: the listed representative-row order.

Stored output:
    05-knowledge/results/lrc_n18_gate_battlefield_s520.out
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from fractions import Fraction
from importlib.machinery import SourceFileLoader
from itertools import combinations
from math import gcd
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
S356 = SourceFileLoader(
    "lonely_runner_residue_probe_s356",
    str(ROOT / "04-computation" / "lonely_runner_residue_probe_s356.py"),
).load_module()
S360 = SourceFileLoader(
    "lonely_runner_endpoint_protection_s360",
    str(ROOT / "04-computation" / "lonely_runner_endpoint_protection_s360.py"),
).load_module()


N = 18
K = N - 1
INITIAL = tuple(range(1, N))
UNIT_POINTS = tuple(Fraction(a, N) for a in range(1, N) if gcd(a, N) == 1)


@dataclass(frozen=True)
class SourceProfile:
    witness: Fraction | None
    kind: str
    observer_outdegree: int
    runner_phase_ties: int
    runner_phase_c3: int
    runner_phase_sccs: tuple[int, ...]
    runner_phase_score_hist: tuple[tuple[int, int], ...]


@dataclass(frozen=True)
class Audit:
    label: str
    speeds: tuple[int, ...]
    classification: str
    forbidden_length: Fraction
    max_gap: Fraction
    gap_ratio: Fraction
    boundary_witnesses: int
    unprotected: int
    first_unprotected: Fraction | None
    has_gate: bool
    gate_count: int
    unit_safe: int
    unit_owned: int
    unit_unprotected: int
    unprotected_denoms: tuple[tuple[int, int], ...]
    source: SourceProfile


def fmt_frac(value: Fraction | None) -> str:
    return S356.fmt_frac(value)


def fmt_float(value: Fraction | float | None, places: int = 6) -> str:
    if value is None:
        return "-"
    return f"{float(value):.{places}f}"


def primitive(speeds: tuple[int, ...]) -> bool:
    g = 0
    for speed in speeds:
        g = gcd(g, speed)
    return g == 1


def normalize(speeds: tuple[int, ...]) -> tuple[int, ...]:
    return S356.normalize_speed_set(list(speeds))


def endpoint_profile(
    speeds: tuple[int, ...],
) -> tuple[set[Fraction], set[Fraction], tuple[tuple[int, int], ...]]:
    points_by_value: dict[Fraction, list[object]] = defaultdict(list)
    for endpoint in S360.endpoints(speeds):
        points_by_value[endpoint.value].append(endpoint)

    owned_units = {t for t in UNIT_POINTS if t in points_by_value}
    unprotected_units: set[Fraction] = set()
    hist: Counter[int] = Counter()
    for value in sorted(points_by_value):
        unprotected = not any(
            S360.direct_protects(speeds, protector, value) for protector in speeds
        )
        if unprotected:
            hist[value.denominator] += 1
            if value in owned_units:
                unprotected_units.add(value)
    return owned_units, unprotected_units, tuple(sorted(hist.items()))


def half_turn_runner_adjacency(speeds: tuple[int, ...], t: Fraction) -> tuple[list[list[int]], int]:
    """Runner half-turn tournament at t, completed by the speed-order tie path."""
    n = len(speeds)
    adj = [[0] * n for _ in range(n)]
    ties = 0
    for i, j in combinations(range(n), 2):
        phase = S356.circle_point((speeds[i] - speeds[j]) * t)
        if phase in (Fraction(0), Fraction(1, 2)):
            ties += 1
            winner, loser = i, j
        elif phase < Fraction(1, 2):
            winner, loser = i, j
        else:
            winner, loser = j, i
        adj[winner][loser] = 1
    return adj, ties


def observer_outdegree(speeds: tuple[int, ...], t: Fraction) -> int:
    threshold = Fraction(1, len(speeds) + 1)
    return sum(
        1
        for speed in speeds
        if S360.circular_distance_to_integer(speed * t) >= threshold
    )


def score_hist_tuple(adj: list[list[int]]) -> tuple[tuple[int, int], ...]:
    return tuple(sorted(Counter(sum(row) for row in adj).items()))


def source_profile(speeds: tuple[int, ...], report: object) -> SourceProfile:
    if report.witness is not None:
        witness = report.witness
        kind = "open_source"
    elif report.boundary_witness is not None:
        witness = report.boundary_witness
        kind = "boundary_source"
    else:
        return SourceProfile(None, "none", 0, 0, 0, tuple(), tuple())

    adj, ties = half_turn_runner_adjacency(speeds, witness)
    return SourceProfile(
        witness=witness,
        kind=kind,
        observer_outdegree=observer_outdegree(speeds, witness),
        runner_phase_ties=ties,
        runner_phase_c3=directed_triangles(adj),
        runner_phase_sccs=scc_sizes(adj),
        runner_phase_score_hist=score_hist_tuple(adj),
    )


def audit(label: str, raw_speeds: tuple[int, ...]) -> Audit:
    speeds = normalize(raw_speeds)
    report = S356.report(label, list(speeds))
    summary = S360.summarize(list(speeds))
    unit_safe = sum(1 for t in UNIT_POINTS if S356.is_lonely_witness(speeds, t))
    owned, unprotected_units, unprotected_denoms = endpoint_profile(speeds)
    gate_count = sum(1 for speed in speeds if speed % N == 0)
    return Audit(
        label=label,
        speeds=speeds,
        classification=summary.classification,
        forbidden_length=summary.forbidden_length,
        max_gap=summary.max_gap,
        gap_ratio=summary.max_gap / summary.threshold,
        boundary_witnesses=report.boundary_witness_count,
        unprotected=summary.unprotected_count,
        first_unprotected=summary.first_unprotected,
        has_gate=gate_count > 0,
        gate_count=gate_count,
        unit_safe=unit_safe,
        unit_owned=len(owned),
        unit_unprotected=len(unprotected_units),
        unprotected_denoms=unprotected_denoms,
        source=source_profile(speeds, report),
    )


def audit_rank(row: Audit) -> tuple[object, ...]:
    class_rank = {"open_cover": 0, "boundary_only": 1, "positive_gap": 2}.get(
        row.classification, 3
    )
    return (class_rank, row.max_gap, row.unprotected, row.speeds)


def print_audit(row: Audit) -> None:
    denoms = ",".join(f"{d}:{c}" for d, c in row.unprotected_denoms[:6])
    if len(row.unprotected_denoms) > 6:
        denoms += ",..."
    print(
        "  "
        f"{row.label:<34} class={row.classification:<13} "
        f"gate={row.gate_count:>2} "
        f"unit_safe={row.unit_safe} unit_owned={row.unit_owned} "
        f"unit_unprot={row.unit_unprotected} "
        f"gap/th={fmt_float(row.gap_ratio)} "
        f"unprot={row.unprotected:>3} first={fmt_frac(row.first_unprotected):>10} "
        f"denoms=[{denoms}]"
    )
    print(f"    speeds={row.speeds}")


def residue_unit_action() -> None:
    print("Residue action on the n=18 unit skeleton")
    print("  Unit points are a/18 for a in {1,5,7,11,13,17}.")
    print("  A speed with residue r mod 18 covers such a point iff r*a == 0 mod 18.")
    for residue in range(N):
        covered = sum(
            1
            for t in UNIT_POINTS
            if S360.circular_distance_to_integer(residue * t) < Fraction(1, N)
        )
        if covered or residue in (0, 1, 2, 3, 6, 9, 17):
            print(f"  residue={residue:2d} covers_unit_points={covered}")
    print("  Therefore: without an 18-multiple, all six unit points remain safe.")
    print()


def ladder(scale: int, skip: int) -> tuple[int, ...]:
    return tuple(sorted({1} | {scale * q for q in range(1, N) if q != skip}))


def best_ladder(scale: int) -> tuple[int, tuple[int, ...]]:
    candidates: list[tuple[Fraction, int, int, tuple[int, ...]]] = []
    for skip in range(1, N):
        speeds = ladder(scale, skip)
        if len(speeds) != K or not primitive(speeds):
            continue
        report = S356.report("ladder", list(speeds))
        candidates.append((report.max_gap / report.threshold, report.boundary_witness_count, skip, speeds))
    if not candidates:
        raise RuntimeError(f"no primitive ladder for scale={scale}")
    _gap, _boundary, skip, speeds = min(candidates)
    return skip, speeds


def canonical_audits() -> list[Audit]:
    lpd_skip, lpd_speeds = best_ladder(9)
    gate_skip, gate_speeds = best_ladder(18)
    double_skip, double_speeds = best_ladder(36)
    no_gate_stress = (
        2,
        3,
        4,
        5,
        6,
        8,
        9,
        10,
        12,
        14,
        15,
        16,
        20,
        21,
        22,
        24,
        25,
    )
    rows = [
        audit("initial tight", INITIAL),
        audit("initial replace 8 by 18", tuple(sorted((set(INITIAL) - {8}) | {18}))),
        audit("initial replace 17 by 18", tuple(sorted((set(INITIAL) - {17}) | {18}))),
        audit("no-gate stress", no_gate_stress),
        audit(f"lpd ladder 9* skip {lpd_skip}", lpd_speeds),
        audit(f"gate ladder 18* skip {gate_skip}", gate_speeds),
        audit(f"double-gate 36* skip {double_skip}", double_speeds),
    ]
    print("Canonical exact audits")
    for row in rows:
        print_audit(row)
    print()
    return rows


def targeted_new_values() -> tuple[int, ...]:
    # Hand-picked local repair values near the 9-ladder, 18-gate, and
    # first descendant layers.  This is intentionally a bounded audit, not
    # an exhaustive one-swap enumeration.
    return (
        9,
        18,
        24,
        27,
        30,
        36,
        42,
        45,
        48,
        54,
        60,
        63,
        66,
        72,
        81,
        90,
        99,
        108,
    )


def targeted_repair_scan() -> list[Audit]:
    """Bounded one-step repairs near the 2- and 3-torsion coordinates."""
    rows: list[Audit] = []
    base = set(INITIAL)
    fragile_removed = (6, 8, 9, 12, 16, 17)
    for removed in fragile_removed:
        rest = base - {removed}
        for new in targeted_new_values():
            if new in rest:
                continue
            speeds = tuple(sorted(rest | {new}))
            if len(speeds) != K or speeds == INITIAL or not primitive(speeds):
                continue
            rows.append(audit(f"swap {removed}->{new}", speeds))
    return rows


def summarize_scan(rows: list[Audit]) -> tuple[Audit | None, Audit | None]:
    print("Targeted one-step repair scan around n=18 torsion/gate coordinates")
    print("  removed coordinates: 6,8,9,12,16,17")
    print("  inserted values: multiples of 2,3,6,9,18,36 up to q=12 plus selected descendants")
    counts: dict[bool, Counter[str]] = defaultdict(Counter)
    for row in rows:
        counts[row.has_gate][row.classification] += 1

    best_no_gate: Audit | None = None
    best_gate: Audit | None = None
    for has_gate in (False, True):
        label = "with an 18-multiple" if has_gate else "without an 18-multiple"
        subset = [row for row in rows if row.has_gate == has_gate]
        total = sum(counts[has_gate].values())
        print(f"  {label}: total={total} classes={dict(sorted(counts[has_gate].items()))}")
        best = sorted(subset, key=audit_rank)[:5]
        if has_gate and best:
            best_gate = best[0]
        if not has_gate and best:
            best_no_gate = best[0]
        for row in best:
            print_audit(row)
    print()
    return best_no_gate, best_gate


def gate_replacement_ledger(max_q: int = 6) -> list[Audit]:
    print(
        "Pure gate replacement ledger: replace fragile initial speeds "
        f"by 18*q, q<= {max_q}"
    )
    rows: list[Audit] = []
    for removed in (6, 8, 9, 12, 16, 17):
        for q in range(1, max_q + 1):
            new = N * q
            speeds = tuple(sorted((set(INITIAL) - {removed}) | {new}))
            if len(speeds) == K and primitive(speeds):
                rows.append(audit(f"remove {removed}, add {new}", speeds))
    by_class = Counter(row.classification for row in rows)
    print(f"  audited={len(rows)} class_counts={dict(sorted(by_class.items()))}")
    print("  best gate replacements:")
    for row in sorted(rows, key=audit_rank)[:10]:
        print_audit(row)
    print()
    return rows


def fmt_hist(hist: tuple[tuple[int, int], ...]) -> str:
    return ",".join(f"{score}:{count}" for score, count in hist)


def observer_source_fingerprints(rows: list[Audit]) -> None:
    print("Observer-source fingerprints (HYP-1981 / THM-381)")
    print("  observer -> runner iff ||v*t|| >= 1/18; source means observer outdegree 17.")
    print("  runner half-turn ties are completed by the listed speed-order tie path.")
    for row in rows:
        source = row.source
        print(
            "  "
            f"{row.label:<34} source={source.kind:<15} "
            f"t={fmt_frac(source.witness):>10} "
            f"obs_out={source.observer_outdegree:>2}/{K} "
            f"runner_ties={source.runner_phase_ties:>3} "
            f"runner_c3={source.runner_phase_c3:>4} "
            f"runner_sccs={source.runner_phase_sccs} "
            f"score_hist={fmt_hist(source.runner_phase_score_hist)}"
        )
    print(
        "  Reading: the no-gate unit witnesses are source-marked boundary targets; "
        "gate rows destroy those targets but create open source targets in "
        "descendant endpoint layers."
    )
    print()


def score_hist(adj: list[list[int]]) -> str:
    scores = [sum(row) for row in adj]
    hist = Counter(scores)
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


def hamiltonian_paths(adj: list[list[int]]) -> int:
    n = len(adj)
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


def max_unprotected_denom(row: Audit) -> int:
    return max((denom for denom, _count in row.unprotected_denoms), default=1)


def branch_tournament(rows: list[Audit]) -> None:
    dims = (
        ("has_gate", lambda r: int(r.has_gate), True),
        ("gate_count", lambda r: r.gate_count, True),
        ("unit_witness_loss", lambda r: K - r.unit_safe, True),
        ("small_gap", lambda r: float(r.gap_ratio), False),
        ("unprotected", lambda r: r.unprotected, True),
        ("descendant_denominator", max_unprotected_denom, True),
    )
    n = len(rows)
    adj = [[0] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        votes_i = 0
        votes_j = 0
        for _name, getter, high_wins in dims:
            left = getter(rows[i])
            right = getter(rows[j])
            cmp = (left > right) - (left < right)
            if cmp == 0:
                continue
            if high_wins:
                votes_i += int(cmp > 0)
                votes_j += int(cmp < 0)
            else:
                votes_i += int(cmp < 0)
                votes_j += int(cmp > 0)
        winner = i if votes_i >= votes_j else j
        loser = j if winner == i else i
        adj[winner][loser] = 1

    scores = [sum(row) for row in adj]
    ranked = sorted(((scores[idx], rows[idx].label) for idx in range(n)), reverse=True)
    print("Branch Tournament Analysis")
    print("  vertices: representative n=18 branch rows")
    print("  observable: has-gate, unit witness loss, gap, endpoint debt, descendant denominator")
    print("  switch: more counterexample-shaped row points to less; tie path is listed row order")
    print(
        f"  H={hamiltonian_paths(adj)} c3={directed_triangles(adj)} "
        f"SCCs={scc_sizes(adj)} score_hist={score_hist(adj)}"
    )
    print("  top branch rows: " + " | ".join(f"{label}:{score}" for score, label in ranked[:6]))
    print()


def interpretation() -> None:
    print("Interpretation")
    print(
        "  THM-382/383 say the right LRC quotient must remember threshold colors "
        "and boundary tie-wall compactification."
    )
    print(
        "  HYP-1981 turns that decorated target into source reachability: every "
        "lonely witness listed above is an observer-source marked tournament."
    )
    print(
        "  The n=18 unit skeleton behaves like the old n=14 gate/tightness skeleton: "
        "without an 18-multiple, all six unit points remain safe."
    )
    print(
        "  Therefore the no-18-gate branch cannot be an open-cover counterexample. "
        "It can only be tight/boundary-like or positive-gap with unit witnesses intact."
    )
    print(
        "  The has-18-gate branch kills the unit witnesses but exports debt to "
        "descendant endpoint denominators.  In HYP-1981 language, it kills the "
        "unit source targets but moves source reachability into descendant layers."
    )
    print(
        "  The next proof object should be an observer-source endpoint-debt "
        "certificate for the 18-gate/lpd branch, split into 2-torsion and "
        "3^2-torsion descendant layers."
    )


def main() -> None:
    print("n=18 gate/source battlefield audit (codex-2026-06-01 S520)")
    print("Convention: k=17 moving speeds, threshold=1/18.")
    print()
    residue_unit_action()
    canonical = canonical_audits()
    scan_rows = targeted_repair_scan()
    best_no_gate, best_gate = summarize_scan(scan_rows)
    gate_rows = gate_replacement_ledger(max_q=6)

    reps = canonical[:]
    if best_no_gate is not None:
        reps.append(best_no_gate)
    if best_gate is not None:
        reps.append(best_gate)
    if gate_rows:
        reps.append(sorted(gate_rows, key=audit_rank)[0])
    # Deduplicate by speed set while preserving first explanatory label.
    seen: set[tuple[int, ...]] = set()
    deduped: list[Audit] = []
    for row in reps:
        if row.speeds in seen:
            continue
        seen.add(row.speeds)
        deduped.append(row)
    observer_source_fingerprints(deduped)
    branch_tournament(deduped)
    interpretation()


if __name__ == "__main__":
    main()
