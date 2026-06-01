#!/usr/bin/env python3
"""
tournament_analysis_switchboard_s454.py

codex-2026-06-01 S454

Second Tournament Analysis pass.

S23 established the ranker/analyzer split.  This script makes the user's
"switch each of the N choose 2 things by a metric" idea into a first-class
object: a switchboard path through the tournament cube.

Each pairwise measurement supplies either:

* an antisymmetric flux, already a binary relation after tie-break;
* a scalar rank shadow, which usually collapses to a transitive tournament;
* a symmetric metric switch, which toggles a fixed Hamiltonian-path edge;
* a lens, where pair i,j is compared through a label, anchor, volume, or
  information-geometric asymmetry.

The output tracks live edges, flip density, Hamming diameter, cycle ranges, and
Hamiltonian-path ranges for basketball, runners, cuboids, spheres, and simplex
probability points.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass
from importlib.machinery import SourceFileLoader
from itertools import combinations
from math import acos, cos, exp, gcd, log, pi, sin, sqrt
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
S23 = SourceFileLoader(
    "tournament_analysis_metric_lifts_s23",
    str(ROOT / "04-computation" / "tournament_analysis_metric_lifts_s23.py"),
).load_module()

EPS = 1.0e-12


@dataclass(frozen=True)
class PathSummary:
    family: str
    lift: str
    kind: str
    n: int
    samples: int
    states: int
    live_edges: int
    total_flips: int
    max_step_flips: int
    hamming_diameter: int
    transitive: int
    cycle_range: tuple[int, int]
    h_range: tuple[int | None, int | None]
    scc_range: tuple[int, int]


def fmt_range(values: tuple[int | None, int | None]) -> str:
    lo, hi = values
    if lo is None or hi is None:
        return "-"
    return str(lo) if lo == hi else f"{lo}..{hi}"


def edge_count(n: int) -> int:
    return n * (n - 1) // 2


def bit_at(bits: int, idx: int) -> int:
    return 1 if bits & (1 << idx) else 0


def summarize_path(family: str, lift: str, kind: str, adjs: list[list[list[int]]]) -> PathSummary:
    rows = [S23.stats(adj) for adj in adjs]
    bits = [row.bits for row in rows]
    flips = [(bits[i] ^ bits[i - 1]).bit_count() for i in range(1, len(bits))]
    m = edge_count(rows[0].n)
    live_edges = sum(
        len({bit_at(value, edge) for value in bits}) > 1
        for edge in range(m)
    )
    diameter = 0
    for i, left in enumerate(bits):
        for right in bits[i + 1:]:
            diameter = max(diameter, (left ^ right).bit_count())
    h_values = [row.hamiltonian_paths for row in rows if row.hamiltonian_paths is not None]
    return PathSummary(
        family=family,
        lift=lift,
        kind=kind,
        n=rows[0].n,
        samples=len(rows),
        states=len(set(bits)),
        live_edges=live_edges,
        total_flips=sum(flips),
        max_step_flips=max(flips) if flips else 0,
        hamming_diameter=diameter,
        transitive=sum(1 for row in rows if row.cyclic_triples == 0),
        cycle_range=(min(row.cyclic_triples for row in rows), max(row.cyclic_triples for row in rows)),
        h_range=(min(h_values), max(h_values)) if h_values else (None, None),
        scc_range=(min(row.scc_count for row in rows), max(row.scc_count for row in rows)),
    )


def print_summary(title: str, rows: list[PathSummary]) -> None:
    print(title)
    print("=" * 126)
    print(
        "  family      lift                      kind        n samp states live flips "
        "maxF diam trans cycles   H        scc"
    )
    for row in rows:
        print(
            f"  {row.family:<11} {row.lift:<25} {row.kind:<10} "
            f"{row.n:>2} {row.samples:>4} {row.states:>6} {row.live_edges:>4} "
            f"{row.total_flips:>5} {row.max_step_flips:>4} {row.hamming_diameter:>4} "
            f"{row.transitive:>5} {fmt_range(row.cycle_range):<8} "
            f"{fmt_range(row.h_range):<8} {fmt_range(row.scc_range)}"
        )
    print()


def print_clusters(family: str, sequences: dict[str, tuple[int, ...]]) -> None:
    buckets: defaultdict[tuple[int, ...], list[str]] = defaultdict(list)
    for name, seq in sequences.items():
        buckets[seq].append(name)
    print(f"{family} exact switchboard clusters")
    print("=" * 126)
    found = False
    for names in buckets.values():
        if len(names) > 1:
            found = True
            print(f"  same path: {', '.join(names)}")
    if not found:
        print("  no exact coincidences")
    print()


def orient_by_antisym(n: int, value) -> list[list[int]]:
    def beats(i: int, j: int) -> bool:
        out = value(i, j)
        if abs(out) <= EPS:
            return S23.tie_beats(i, j)
        return out > 0

    return S23.orient_by_edge_rule(n, beats)


def orient_switch(n: int, switch_value) -> list[list[int]]:
    """Switch a fixed Hamiltonian path edge by the sign of a pair metric."""

    def beats(i: int, j: int) -> bool:
        value = switch_value(i, j)
        if abs(value) <= EPS:
            return S23.tie_beats(i, j)
        keep_base = value > 0
        return S23.tie_beats(i, j) if keep_base else S23.tie_beats(j, i)

    return S23.orient_by_edge_rule(n, beats)


def bits_sequence(adjs: list[list[list[int]]]) -> tuple[int, ...]:
    return tuple(S23.bits_from_adj(adj) for adj in adjs)


def basketball_window(k: int) -> list[list[int]]:
    """Synthetic five-starter directed pass counts for one game window."""

    n = 5
    matrix = [[0] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        base = 6 + ((i + 2) * (j + 3) + 2 * k) % 5
        pulse = int(3 * (1 + sin((k + 1) * (i + 1) + (j + 2))))
        counter = int(3 * (1 + cos((k + 2) * (j + 1) - i)))
        matrix[i][j] = base + pulse + ((i + k) % 3)
        matrix[j][i] = base + counter + ((j + 2 * k) % 3)
        if (i + j + k) % 7 == 0:
            matrix[j][i] = matrix[i][j]
    return matrix


def basketball_series() -> tuple[list[PathSummary], dict[str, tuple[int, ...]]]:
    matrices = [basketball_window(k) for k in range(8)]
    n = 5
    adjs_by_lift: dict[str, tuple[str, list[list[list[int]]]]] = {
        "pass-flux": ("flux", []),
        "assist-rank": ("rank", []),
        "reciprocity-switch": ("switch", []),
        "two-hop-lens": ("lens", []),
        "pressure-switch": ("switch", []),
    }
    for matrix in matrices:
        totals = [sum(matrix[i]) - sum(matrix[j][i] for j in range(n)) for i in range(n)]
        pair_totals = [matrix[i][j] + matrix[j][i] for i, j in combinations(range(n), 2)]
        median_pair = sorted(pair_totals)[len(pair_totals) // 2]

        adjs_by_lift["pass-flux"][1].append(
            orient_by_antisym(n, lambda i, j, m=matrix: m[i][j] - m[j][i])
        )
        adjs_by_lift["assist-rank"][1].append(S23.orient_by_score(totals))
        adjs_by_lift["reciprocity-switch"][1].append(
            orient_switch(n, lambda i, j, m=matrix, med=median_pair: (m[i][j] + m[j][i]) - med)
        )
        adjs_by_lift["two-hop-lens"][1].append(
            orient_by_antisym(
                n,
                lambda i, j, m=matrix: sum(
                    min(m[i][k], m[k][j]) - min(m[j][k], m[k][i])
                    for k in range(n)
                    if k not in (i, j)
                ),
            )
        )
        adjs_by_lift["pressure-switch"][1].append(
            orient_switch(
                n,
                lambda i, j, m=matrix: abs(m[i][j] - m[j][i])
                - ((m[i][j] + m[j][i]) / 4.0),
            )
        )
    summaries = [
        summarize_path("basket5", lift, kind, adjs)
        for lift, (kind, adjs) in adjs_by_lift.items()
    ]
    sequences = {lift: bits_sequence(adjs) for lift, (_kind, adjs) in adjs_by_lift.items()}
    return summaries, sequences


def circle_positions(speeds: list[int], t: float) -> list[float]:
    return [(speed * t) % 1.0 for speed in speeds]


def circ_delta(a: float, b: float) -> float:
    return (b - a) % 1.0


def circ_dist(a: float, b: float) -> float:
    delta = circ_delta(a, b)
    return min(delta, 1.0 - delta)


def chord_dist(a: float, b: float) -> float:
    return sqrt(max(0.0, 2.0 - 2.0 * cos(2.0 * pi * circ_delta(a, b))))


def circle_derivative_value(speeds: list[int], positions: list[float], i: int, j: int) -> float:
    dt = 1.0e-5
    now = chord_dist(positions[i], positions[j])
    later_positions = [((positions[k] + speeds[k] * dt) % 1.0) for k in range(len(speeds))]
    later = chord_dist(later_positions[i], later_positions[j])
    return later - now


def circle_series() -> tuple[list[PathSummary], dict[str, tuple[int, ...]]]:
    speeds = [0, 1, 2, 3, 5, 8, 13, 21]
    times = [k / 37.0 for k in range(1, 37)]
    lifts: dict[str, tuple[str, list[list[list[int]]]]] = {
        "phase-halfturn": ("phase", []),
        "chord-median-switch": ("switch", []),
        "chord-annulus-switch": ("switch", []),
        "chord-resonance-7": ("switch", []),
        "approach-switch": ("switch", []),
        "relative-speed-lens": ("lens", []),
        "local-gap-rank": ("rank", []),
    }
    for t in times:
        pos = circle_positions(speeds, t)
        pair_dists = [chord_dist(pos[i], pos[j]) for i, j in combinations(range(len(pos)), 2)]
        sorted_dists = sorted(pair_dists)
        q1 = sorted_dists[len(sorted_dists) // 4]
        q2 = sorted_dists[len(sorted_dists) // 2]
        q3 = sorted_dists[3 * len(sorted_dists) // 4]
        lifts["phase-halfturn"][1].append(S23.phase_halfturn_tournament(pos))
        lifts["chord-median-switch"][1].append(
            orient_switch(len(pos), lambda i, j, p=pos, med=q2: chord_dist(p[i], p[j]) - med)
        )
        lifts["chord-annulus-switch"][1].append(
            orient_switch(
                len(pos),
                lambda i, j, p=pos, lo=q1, hi=q3: min(chord_dist(p[i], p[j]) - lo, hi - chord_dist(p[i], p[j])),
            )
        )
        lifts["chord-resonance-7"][1].append(
            orient_switch(len(pos), lambda i, j, p=pos: sin(7.0 * pi * chord_dist(p[i], p[j])))
        )
        lifts["approach-switch"][1].append(
            orient_switch(len(pos), lambda i, j, p=pos, s=speeds: -circle_derivative_value(s, p, i, j))
        )
        lifts["relative-speed-lens"][1].append(
            orient_by_antisym(
                len(pos),
                lambda i, j, p=pos, s=speeds: (
                    circ_dist(p[i], p[(j + 1) % len(p)]) / (1 + abs(s[i] - s[(j + 1) % len(p)]))
                    - circ_dist(p[j], p[(i + 1) % len(p)]) / (1 + abs(s[j] - s[(i + 1) % len(p)]))
                ),
            )
        )
        gaps = []
        order = sorted(range(len(pos)), key=lambda idx: (pos[idx], idx))
        score = [0.0] * len(pos)
        for idx, vertex in enumerate(order):
            prev_vertex = order[(idx - 1) % len(pos)]
            next_vertex = order[(idx + 1) % len(pos)]
            score[vertex] = min(circ_delta(pos[prev_vertex], pos[vertex]), circ_delta(pos[vertex], pos[next_vertex]))
            gaps.append(score[vertex])
        lifts["local-gap-rank"][1].append(S23.orient_by_score(score))
    summaries = [
        summarize_path("circle8", lift, kind, adjs)
        for lift, (kind, adjs) in lifts.items()
    ]
    sequences = {lift: bits_sequence(adjs) for lift, (_kind, adjs) in lifts.items()}
    return summaries, sequences


def cube_points(speed_vectors: list[tuple[int, int, int]], t: float) -> list[tuple[float, float, float]]:
    return [tuple((component * t) % 1.0 for component in vector) for vector in speed_vectors]


def lp_dist(p: tuple[float, ...], q: tuple[float, ...], power: float) -> float:
    diffs = [abs(a - b) for a, b in zip(p, q)]
    if power == float("inf"):
        return max(diffs)
    return sum(value**power for value in diffs) ** (1.0 / power)


def cube_series() -> tuple[list[PathSummary], dict[str, tuple[int, ...]]]:
    speed_vectors = [
        (0, 0, 0),
        (1, 2, 3),
        (2, 3, 5),
        (3, 5, 8),
        (5, 8, 13),
        (8, 13, 21),
        (13, 21, 34),
        (21, 34, 55),
    ]
    times = [k / 41.0 for k in range(1, 41)]
    lifts: dict[str, tuple[str, list[list[list[int]]]]] = {
        "L2-row-rank": ("rank", []),
        "L2-median-switch": ("switch", []),
        "Linf-resonance-5": ("switch", []),
        "dominant-axis-lens": ("lens", []),
        "volume-drift-lens": ("lens", []),
    }
    drift = (1.0, sqrt(2.0), (1.0 + sqrt(5.0)) / 2.0)

    def det3(a, b, c) -> float:
        return (
            a[0] * (b[1] * c[2] - b[2] * c[1])
            - a[1] * (b[0] * c[2] - b[2] * c[0])
            + a[2] * (b[0] * c[1] - b[1] * c[0])
        )

    for t in times:
        pts = cube_points(speed_vectors, t)
        l2_pairs = [lp_dist(pts[i], pts[j], 2.0) for i, j in combinations(range(len(pts)), 2)]
        med = sorted(l2_pairs)[len(l2_pairs) // 2]
        row_scores = [
            sum(lp_dist(pts[i], pts[j], 2.0) for j in range(len(pts)) if i != j)
            for i in range(len(pts))
        ]
        lifts["L2-row-rank"][1].append(S23.orient_by_score(row_scores))
        lifts["L2-median-switch"][1].append(
            orient_switch(len(pts), lambda i, j, p=pts, m=med: lp_dist(p[i], p[j], 2.0) - m)
        )
        lifts["Linf-resonance-5"][1].append(
            orient_switch(len(pts), lambda i, j, p=pts: sin(5.0 * pi * lp_dist(p[i], p[j], float("inf"))))
        )
        lifts["dominant-axis-lens"][1].append(
            orient_by_antisym(
                len(pts),
                lambda i, j, p=pts: (
                    p[i][max(range(3), key=lambda a: abs(p[i][a] - p[j][a]))]
                    - p[j][max(range(3), key=lambda a: abs(p[i][a] - p[j][a]))]
                ),
            )
        )
        lifts["volume-drift-lens"][1].append(
            orient_by_antisym(
                len(pts),
                lambda i, j, p=pts: det3(
                    (p[i][0] - 0.5, p[i][1] - 0.5, p[i][2] - 0.5),
                    (p[j][0] - 0.5, p[j][1] - 0.5, p[j][2] - 0.5),
                    drift,
                ),
            )
        )
    summaries = [
        summarize_path("cube8", lift, kind, adjs)
        for lift, (kind, adjs) in lifts.items()
    ]
    sequences = {lift: bits_sequence(adjs) for lift, (_kind, adjs) in lifts.items()}
    return summaries, sequences


def sphere_points(speed_vectors: list[tuple[int, int]], t: float) -> list[tuple[float, float, float]]:
    points = []
    for a, b in speed_vectors:
        theta = 2.0 * pi * ((a * t) % 1.0)
        z = 2.0 * ((b * t) % 1.0) - 1.0
        radius = sqrt(max(0.0, 1.0 - z * z))
        points.append((radius * cos(theta), radius * sin(theta), z))
    return points


def dot(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def cross(a: tuple[float, float, float], b: tuple[float, float, float]) -> tuple[float, float, float]:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def great_circle(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return acos(max(-1.0, min(1.0, dot(a, b))))


def sphere_series() -> tuple[list[PathSummary], dict[str, tuple[int, ...]]]:
    speeds = [(0, 0), (1, 2), (2, 3), (3, 5), (5, 8), (8, 13), (13, 21), (21, 34)]
    times = [k / 43.0 for k in range(1, 43)]
    north = (0.0, 0.0, 1.0)
    lifts: dict[str, tuple[str, list[list[list[int]]]]] = {
        "greatcircle-rank": ("rank", []),
        "greatcircle-switch": ("switch", []),
        "oriented-normal-lens": ("lens", []),
        "cap-membership-switch": ("switch", []),
    }
    for t in times:
        pts = sphere_points(speeds, t)
        pair_d = [great_circle(pts[i], pts[j]) for i, j in combinations(range(len(pts)), 2)]
        med = sorted(pair_d)[len(pair_d) // 2]
        row_scores = [
            sum(great_circle(pts[i], pts[j]) for j in range(len(pts)) if i != j)
            for i in range(len(pts))
        ]
        lifts["greatcircle-rank"][1].append(S23.orient_by_score(row_scores))
        lifts["greatcircle-switch"][1].append(
            orient_switch(len(pts), lambda i, j, p=pts, m=med: great_circle(p[i], p[j]) - m)
        )
        lifts["oriented-normal-lens"][1].append(
            orient_by_antisym(len(pts), lambda i, j, p=pts: dot(cross(p[i], p[j]), north))
        )
        lifts["cap-membership-switch"][1].append(
            orient_switch(
                len(pts),
                lambda i, j, p=pts: (1 if p[i][2] * p[j][2] >= 0 else -1) * abs(p[i][2] - p[j][2]),
            )
        )
    summaries = [
        summarize_path("sphere8", lift, kind, adjs)
        for lift, (kind, adjs) in lifts.items()
    ]
    sequences = {lift: bits_sequence(adjs) for lift, (_kind, adjs) in lifts.items()}
    return summaries, sequences


def simplex_points(speed_vectors: list[tuple[int, int, int, int]], t: float) -> list[tuple[float, ...]]:
    points = []
    for vec in speed_vectors:
        raw = [exp(sin(2.0 * pi * ((component * t) % 1.0))) for component in vec]
        total = sum(raw)
        points.append(tuple(value / total for value in raw))
    return points


def kl(p: tuple[float, ...], q: tuple[float, ...]) -> float:
    return sum(pi_value * log(pi_value / qi_value) for pi_value, qi_value in zip(p, q))


def l1(p: tuple[float, ...], q: tuple[float, ...]) -> float:
    return sum(abs(a - b) for a, b in zip(p, q))


def simplex_series() -> tuple[list[PathSummary], dict[str, tuple[int, ...]]]:
    speeds = [
        (0, 1, 2, 3),
        (1, 2, 3, 5),
        (2, 3, 5, 8),
        (3, 5, 8, 13),
        (5, 8, 13, 21),
        (8, 13, 21, 34),
        (13, 21, 34, 55),
        (21, 34, 55, 89),
    ]
    times = [k / 47.0 for k in range(1, 47)]
    lifts: dict[str, tuple[str, list[list[list[int]]]]] = {
        "entropy-rank": ("rank", []),
        "L1-switch": ("switch", []),
        "KL-skew-lens": ("lens", []),
        "dominant-coordinate-lens": ("lens", []),
    }
    for t in times:
        pts = simplex_points(speeds, t)
        pair_l1 = [l1(pts[i], pts[j]) for i, j in combinations(range(len(pts)), 2)]
        med = sorted(pair_l1)[len(pair_l1) // 2]
        entropy_scores = [-sum(value * log(value) for value in p) for p in pts]
        lifts["entropy-rank"][1].append(S23.orient_by_score(entropy_scores))
        lifts["L1-switch"][1].append(
            orient_switch(len(pts), lambda i, j, p=pts, m=med: l1(p[i], p[j]) - m)
        )
        lifts["KL-skew-lens"][1].append(
            orient_by_antisym(len(pts), lambda i, j, p=pts: kl(p[j], p[i]) - kl(p[i], p[j]))
        )
        lifts["dominant-coordinate-lens"][1].append(
            orient_by_antisym(
                len(pts),
                lambda i, j, p=pts: (
                    p[i][max(range(len(p[i])), key=lambda a: abs(p[i][a] - p[j][a]))]
                    - p[j][max(range(len(p[i])), key=lambda a: abs(p[i][a] - p[j][a]))]
                ),
            )
        )
    summaries = [
        summarize_path("simplex8", lift, kind, adjs)
        for lift, (kind, adjs) in lifts.items()
    ]
    sequences = {lift: bits_sequence(adjs) for lift, (_kind, adjs) in lifts.items()}
    return summaries, sequences


def primitive_gcd(values: tuple[int, ...]) -> int:
    out = 0
    for value in values:
        out = gcd(out, value)
    return out


def lrc_speed_sets() -> list[tuple[str, tuple[int, ...]]]:
    return [
        ("n14-initial", tuple(range(14))),
        ("n14-seven", (0,) + tuple(sorted({1} | {7 * q for q in range(1, 14) if q != 6}))),
        ("n16-initial", tuple(range(16))),
        ("n16-dyadic", (0,) + tuple(sorted({1} | {8 * q for q in range(1, 16) if q != 14}))),
    ]


def lrc_snapshot_table() -> None:
    print("E. LRC SNAPSHOTS AS TOURNAMENT ANALYSIS SWITCHBOARDS")
    print("=" * 126)
    print("  sample       n  t       phase cycles  safe-switch cycles  relspeed cycles  live warning")
    for label, speeds in lrc_speed_sets():
        n = len(speeds)
        t = 1.0 / n
        pos = circle_positions(list(speeds), t)
        phase = S23.stats(S23.phase_halfturn_tournament(pos), compute_h=False)
        safe = S23.stats(
            orient_switch(n, lambda i, j, p=pos, nn=n: circ_dist(p[i], p[j]) - 1.0 / nn),
            compute_h=False,
        )
        relspeed = S23.stats(
            orient_switch(
                n,
                lambda i, j, p=pos, s=speeds: circ_dist(p[i], p[j]) / (1 + abs(s[i] - s[j])) - (1.0 / (2 * n)),
            ),
            compute_h=False,
        )
        warning = (
            "antipodal ties" if n % 2 == 0 else "odd no antipodes"
        )
        print(
            f"  {label:<12} {n:>2} {t:<7.4f} {phase.cyclic_triples:>12} "
            f"{safe.cyclic_triples:>19} {relspeed.cyclic_triples:>16}  {warning}"
        )
    print()


def print_switchboard_synthesis(all_rows: list[PathSummary]) -> None:
    rank_rows = [row for row in all_rows if row.kind == "rank"]
    analyzer_rows = [row for row in all_rows if row.kind != "rank"]
    rank_trans = sum(row.transitive for row in rank_rows)
    rank_samples = sum(row.samples for row in rank_rows)
    analyzer_trans = sum(row.transitive for row in analyzer_rows)
    analyzer_samples = sum(row.samples for row in analyzer_rows)
    live_edge_mean = sum(row.live_edges for row in analyzer_rows) / len(analyzer_rows)
    print("F. SWITCHBOARD SYNTHESIS")
    print("=" * 126)
    print("  Tournament Analysis is not just extracting a tournament; it is choosing a switchboard.")
    print("  A switchboard has one comparator channel per unordered pair and a fixed tie-break path.")
    print()
    print(f"  rank shadows:     transitive states {rank_trans}/{rank_samples}")
    print(f"  analyzer shadows: transitive states {analyzer_trans}/{analyzer_samples}")
    print(f"  analyzer mean live edges per path: {live_edge_mean:.2f}")
    print()
    print("  human dictionary:")
    print("    basketball pass-flux       = directed social energy")
    print("    reciprocity switch         = symmetric pair intensity toggles a base playbook path")
    print("    circle/sphere phase        = chirality, a tournament-native geometric relation")
    print("    annulus/resonance switch   = arbitrary metric banding on each pair")
    print("    dominant-axis lens         = cuboid pair chooses its own coordinate judge")
    print("    KL-skew lens               = simplex information geometry gives an asymmetric edge")
    print("    LRC safe-distance switch   = pairwise moat relation beside the marked lonely bracket")
    print()
    print("  working principle:")
    print("    rankers summarize objects; switchboards summarize relations.")
    print("    If the problem lives in pairwise tension, the switchboard is closer to the object.")


def sphere_simplex_series() -> tuple[list[PathSummary], dict[str, tuple[int, ...]]]:
    sphere_rows, sphere_sequences = sphere_series()
    simplex_rows, simplex_sequences = simplex_series()
    return sphere_rows + simplex_rows, {**sphere_sequences, **simplex_sequences}


def main() -> None:
    print("TOURNAMENT ANALYSIS SWITCHBOARD (codex-2026-06-01-S454)")
    print("=" * 126)
    print("A metric becomes tournament-relevant only after we choose how each pair gets one bit.")
    print("This pass studies those pairwise bit channels as a switchboard path through the tournament cube.")
    print()

    all_rows: list[PathSummary] = []
    for title, maker in [
        ("A. BASKETBALL: FLUX, RANK, AND PLAYBOOK SWITCHES", basketball_series),
        ("B. CIRCLE RUNNERS: PHASE, CHORD, ANNULUS, AND APPROACH SWITCHES", circle_series),
        ("C. CUBOID SWITCHBOARDS: DISTANCE, AXIS, AND VOLUME LENSES", cube_series),
        ("D. SPHERE AND SIMPLEX SWITCHBOARDS", sphere_simplex_series),
    ]:
        rows, sequences = maker()
        all_rows.extend(rows)
        print_summary(title, rows)
        print_clusters(title.split(":")[0].split(". ")[-1].lower(), sequences)

    lrc_snapshot_table()
    print_switchboard_synthesis(all_rows)


if __name__ == "__main__":
    main()
