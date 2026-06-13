#!/usr/bin/env python3
"""
lrc_neighbor_tournament_lift_s431.py

codex-2026-05-31 S431

Audit what extra tournament-like structure appears when a Lonely Runner
configuration is not collapsed to the single number min_i ||v_i t||.

At a fixed time t, the stationary runner plus moving runners give:

  * a distance-to-zero rank tournament (usually transitive);
  * a circular half-distance tournament, oriented i -> j when j lies within
    the clockwise open semicircle from i;
  * a two-nearest-neighbor directed graph;
  * a pairwise-distance ledger that sees collisions and quotient compression
    invisible to the lonely scalar.

The point of this script is exploratory.  It gives concrete vocabulary for
using tournament tools on LRC: SCC/condensation, good-cut protection,
incomplete tournaments from antipodal ties, and nearest-neighbor cores.
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from fractions import Fraction
from importlib.machinery import SourceFileLoader
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
S356 = SourceFileLoader(
    "lonely_runner_residue_probe_s356",
    str(ROOT / "04-computation" / "lonely_runner_residue_probe_s356.py"),
).load_module()

ONE = Fraction(1, 1)
HALF = Fraction(1, 2)


@dataclass(frozen=True)
class Vertex:
    label: str
    speed: int
    pos: Fraction


@dataclass(frozen=True)
class HalfTournament:
    oriented_arcs: tuple[tuple[str, str], ...]
    antipodal_pairs: tuple[tuple[str, str], ...]
    collision_pairs: tuple[tuple[str, str], ...]
    score_sequence: tuple[int, ...]
    score_hist: tuple[tuple[int, int], ...]
    scc_count: int
    orientation_density: Fraction


@dataclass(frozen=True)
class TwoNeighborGraph:
    arcs: tuple[tuple[str, str], ...]
    mutual_pairs: tuple[tuple[str, str], ...]
    weak_component_sizes: tuple[int, ...]
    outdegree_hist: tuple[tuple[int, int], ...]
    excess_tie_arcs: int


@dataclass(frozen=True)
class LiftSnapshot:
    label: str
    speeds: tuple[int, ...]
    threshold: Fraction
    time: Fraction
    source: str
    lonely_gap: Fraction
    second_zero_gap: Fraction
    pairwise_min_gap: Fraction
    closest_zero: tuple[tuple[str, Fraction, str], ...]
    side_bracket: tuple[str, str]
    residue_multiplicities: tuple[tuple[int, int], ...]
    circular_order: tuple[str, ...]
    half_tournament: HalfTournament
    neighbor_graph: TwoNeighborGraph


def fmt_frac(x: Fraction | None) -> str:
    if x is None:
        return "-"
    if x.denominator == 1:
        return str(x.numerator)
    return f"{x.numerator}/{x.denominator}"


def fmt_ratio(x: Fraction, base: Fraction) -> str:
    if base == 0:
        return "-"
    return f"{float(x / base):.6f}"


def mod1(x: Fraction) -> Fraction:
    return x % ONE


def clockwise(a: Fraction, b: Fraction) -> Fraction:
    return (b - a) % ONE


def circle_distance(a: Fraction, b: Fraction) -> Fraction:
    d = clockwise(a, b)
    return min(d, ONE - d)


def signed_from_zero(x: Fraction) -> tuple[Fraction, str]:
    x = mod1(x)
    if x == 0:
        return Fraction(0), "0"
    if x <= HALF:
        return x, "+"
    return ONE - x, "-"


def signed_text(pos: Fraction) -> str:
    dist, side = signed_from_zero(pos)
    if side == "0":
        return "0"
    return f"{side}{fmt_frac(dist)}"


def vertices_at_time(speeds: tuple[int, ...], t: Fraction) -> tuple[Vertex, ...]:
    return (Vertex("0", 0, Fraction(0)),) + tuple(
        Vertex(f"v{v}", v, mod1(v * t)) for v in speeds
    )


def scc_count(labels: tuple[str, ...], arcs: tuple[tuple[str, str], ...]) -> int:
    adj = {label: [] for label in labels}
    radj = {label: [] for label in labels}
    for a, b in arcs:
        adj[a].append(b)
        radj[b].append(a)

    seen: set[str] = set()
    order: list[str] = []

    def dfs(v: str) -> None:
        seen.add(v)
        for w in adj[v]:
            if w not in seen:
                dfs(w)
        order.append(v)

    for label in labels:
        if label not in seen:
            dfs(label)

    seen.clear()
    components = 0

    def rdfs(v: str) -> None:
        seen.add(v)
        for w in radj[v]:
            if w not in seen:
                rdfs(w)

    for label in reversed(order):
        if label not in seen:
            components += 1
            rdfs(label)
    return components


def weak_component_sizes(
    labels: tuple[str, ...], arcs: tuple[tuple[str, str], ...]
) -> tuple[int, ...]:
    graph = {label: set() for label in labels}
    for a, b in arcs:
        graph[a].add(b)
        graph[b].add(a)
    seen: set[str] = set()
    sizes: list[int] = []
    for label in labels:
        if label in seen:
            continue
        queue = deque([label])
        seen.add(label)
        size = 0
        while queue:
            node = queue.popleft()
            size += 1
            for nxt in graph[node]:
                if nxt not in seen:
                    seen.add(nxt)
                    queue.append(nxt)
        sizes.append(size)
    return tuple(sorted(sizes, reverse=True))


def half_distance_tournament(vertices: tuple[Vertex, ...]) -> HalfTournament:
    scores = Counter({v.label: 0 for v in vertices})
    arcs: list[tuple[str, str]] = []
    antipodal: list[tuple[str, str]] = []
    collisions: list[tuple[str, str]] = []

    for a, b in combinations(vertices, 2):
        cw = clockwise(a.pos, b.pos)
        if cw == 0:
            collisions.append(tuple(sorted((a.label, b.label))))
        elif cw == HALF:
            antipodal.append(tuple(sorted((a.label, b.label))))
        elif cw < HALF:
            arcs.append((a.label, b.label))
            scores[a.label] += 1
        else:
            arcs.append((b.label, a.label))
            scores[b.label] += 1

    total_pairs = len(vertices) * (len(vertices) - 1) // 2
    score_values = tuple(sorted(scores.values()))
    return HalfTournament(
        oriented_arcs=tuple(sorted(arcs)),
        antipodal_pairs=tuple(sorted(antipodal)),
        collision_pairs=tuple(sorted(collisions)),
        score_sequence=score_values,
        score_hist=tuple(sorted(Counter(score_values).items())),
        scc_count=scc_count(tuple(v.label for v in vertices), tuple(arcs)),
        orientation_density=Fraction(len(arcs), total_pairs) if total_pairs else ONE,
    )


def two_nearest_neighbor_graph(vertices: tuple[Vertex, ...]) -> TwoNeighborGraph:
    arcs: list[tuple[str, str]] = []
    outdegrees = Counter()
    excess = 0

    for v in vertices:
        distances = sorted(
            (circle_distance(v.pos, w.pos), w.label)
            for w in vertices
            if w.label != v.label
        )
        if not distances:
            cutoff = Fraction(0)
        else:
            # "Two nearest" means the second order statistic, with all ties
            # at that distance retained.  On a regular polygon this selects
            # the two adjacent neighbors, not the whole next shell.
            cutoff = distances[min(1, len(distances) - 1)][0]
        neighbors = [label for d, label in distances if d <= cutoff]
        outdegrees[v.label] = len(neighbors)
        excess += max(0, len(neighbors) - 2)
        for w in neighbors:
            arcs.append((v.label, w))

    arc_set = set(arcs)
    mutual = []
    for a, b in combinations(sorted(v.label for v in vertices), 2):
        if (a, b) in arc_set and (b, a) in arc_set:
            mutual.append((a, b))

    return TwoNeighborGraph(
        arcs=tuple(sorted(arc_set)),
        mutual_pairs=tuple(mutual),
        weak_component_sizes=weak_component_sizes(
            tuple(v.label for v in vertices), tuple(sorted(arc_set))
        ),
        outdegree_hist=tuple(sorted(Counter(outdegrees.values()).items())),
        excess_tie_arcs=excess,
    )


def closest_zero_records(
    vertices: tuple[Vertex, ...], limit: int = 4
) -> tuple[tuple[str, Fraction, str], ...]:
    moving = [v for v in vertices if v.speed != 0]
    records = []
    for v in moving:
        dist, side = signed_from_zero(v.pos)
        records.append((dist, v.label, side))
    records.sort(key=lambda item: (item[0], item[1]))
    return tuple((label, dist, side) for dist, label, side in records[:limit])


def side_bracket(vertices: tuple[Vertex, ...]) -> tuple[str, str]:
    left: list[tuple[Fraction, str]] = []
    right: list[tuple[Fraction, str]] = []
    for v in vertices:
        if v.speed == 0:
            continue
        dist, side = signed_from_zero(v.pos)
        if side == "+":
            right.append((dist, v.label))
        elif side == "-":
            left.append((dist, v.label))
        elif side == "0":
            left.append((dist, v.label))
            right.append((dist, v.label))
    left_label = min(left)[1] if left else "-"
    right_label = min(right)[1] if right else "-"
    return left_label, right_label


def residue_multiplicities(speeds: tuple[int, ...], t: Fraction) -> tuple[tuple[int, int], ...]:
    q = t.denominator
    p = t.numerator
    residues = [0] + [(v * p) % q for v in speeds]
    return tuple(sorted(Counter(residues).items()))


def circular_order(vertices: tuple[Vertex, ...]) -> tuple[str, ...]:
    groups: defaultdict[Fraction, list[str]] = defaultdict(list)
    for v in vertices:
        groups[v.pos].append(v.label)
    out = []
    for pos in sorted(groups):
        labels = "/".join(sorted(groups[pos]))
        out.append(f"{labels}@{fmt_frac(pos)}")
    return tuple(out)


def pairwise_min_gap(vertices: tuple[Vertex, ...]) -> Fraction:
    return min(
        circle_distance(a.pos, b.pos) for a, b in combinations(vertices, 2)
    )


def choose_time(label: str, speeds: tuple[int, ...]) -> tuple[Fraction, str]:
    row = S356.report(label, list(speeds))
    if row.witness is not None:
        return row.witness, "positive-gap midpoint"
    if row.boundary_witness is not None:
        return row.boundary_witness, "boundary witness"
    raise ValueError(f"no lonely witness found by S356 for {label}")


def snapshot(label: str, raw_speeds: list[int]) -> LiftSnapshot:
    speeds = S356.normalize_speed_set(raw_speeds)
    t, source = choose_time(label, speeds)
    vertices = vertices_at_time(speeds, t)
    zero_distances = sorted(
        circle_distance(Fraction(0), v.pos) for v in vertices if v.speed != 0
    )
    threshold = Fraction(1, len(speeds) + 1)
    return LiftSnapshot(
        label=label,
        speeds=speeds,
        threshold=threshold,
        time=t,
        source=source,
        lonely_gap=zero_distances[0],
        second_zero_gap=zero_distances[1] if len(zero_distances) > 1 else zero_distances[0],
        pairwise_min_gap=pairwise_min_gap(vertices),
        closest_zero=closest_zero_records(vertices),
        side_bracket=side_bracket(vertices),
        residue_multiplicities=residue_multiplicities(speeds, t),
        circular_order=circular_order(vertices),
        half_tournament=half_distance_tournament(vertices),
        neighbor_graph=two_nearest_neighbor_graph(vertices),
    )


def print_methodology() -> None:
    print("LRC TOURNAMENT-LIFT METHODOLOGY")
    print("=" * 96)
    rows = [
        (
            "scalar lonely gap",
            "min_i ||v_i t||",
            "height of a lower envelope",
            "too lossy: forgets all pairwise runner relations",
        ),
        (
            "danger rank tournament",
            "orient i->j if i is closer to 0 than j",
            "momentary threat order",
            "transitive, useful mostly as a braid/order-change trace",
        ),
        (
            "circular half-distance tournament",
            "i->j if j is clockwise within 1/2 from i",
            "round/incomplete tournament from positions",
            "sees antipodal ties, SCCs, and circular order strata",
        ),
        (
            "two-nearest graph",
            "each runner points to its two nearest neighbors",
            "local contact graph",
            "possible LRC analogue of tournament good-cut cores",
        ),
        (
            "pairwise distance ledger",
            "all ||(v_i-v_j)t||",
            "difference-speed LRC shadow",
            "detects quotient collisions invisible to zero-distance scalar",
        ),
        (
            "residue quotient tournament",
            "positions at t=p/q are v_i p mod q",
            "labelled cyclic tournament with multiplicities",
            "tight examples often collapse to small residue polygons",
        ),
    ]
    print(f"{'lens':<34} {'data':<38} {'tournament use'}")
    print("-" * 96)
    for name, data, use, warning in rows:
        print(f"{name:<34} {data:<38} {use}")
        print(f"{'':<34} {'':<38} {warning}")
    print()


def print_snapshot_table(snapshots: tuple[LiftSnapshot, ...]) -> None:
    print("SAFE/TIGHT TIME SNAPSHOTS")
    print("=" * 96)
    header = (
        "label",
        "n",
        "t",
        "gap/th",
        "2nd/th",
        "pair/th",
        "half",
        "2NN",
    )
    print(
        f"{header[0]:<26} {header[1]:>2} {header[2]:>12} "
        f"{header[3]:>9} {header[4]:>9} {header[5]:>9} "
        f"{header[6]:<20} {header[7]}"
    )
    print("-" * 96)
    for s in snapshots:
        h = s.half_tournament
        g = s.neighbor_graph
        half = (
            f"dens={float(h.orientation_density):.3f} "
            f"scc={h.scc_count} tie={len(h.antipodal_pairs)} col={len(h.collision_pairs)}"
        )
        nn = (
            f"mut={len(g.mutual_pairs)} comp={g.weak_component_sizes} "
            f"out={dict(g.outdegree_hist)}"
        )
        print(
            f"{s.label:<26} {len(s.speeds)+1:>2} {fmt_frac(s.time):>12} "
            f"{fmt_ratio(s.lonely_gap, s.threshold):>9} "
            f"{fmt_ratio(s.second_zero_gap, s.threshold):>9} "
            f"{fmt_ratio(s.pairwise_min_gap, s.threshold):>9} "
            f"{half:<20} {nn}"
        )
    print()


def print_snapshot_details(snapshots: tuple[LiftSnapshot, ...]) -> None:
    print("DETAILS: WHAT THE SCALAR GAP FORGETS")
    print("=" * 96)
    for s in snapshots:
        print(f"[{s.label}]")
        print(f"  speeds={s.speeds}")
        print(
            "  "
            f"time={fmt_frac(s.time)} ({s.source}); "
            f"threshold={fmt_frac(s.threshold)}; "
            f"lonely_gap={fmt_frac(s.lonely_gap)}; "
            f"second_zero_gap={fmt_frac(s.second_zero_gap)}"
        )
        closest = ", ".join(
            f"{label}:{side}{fmt_frac(dist)}"
            for label, dist, side in s.closest_zero
        )
        print(f"  closest_to_zero={closest}")
        print(
            "  "
            f"side_bracket=(left {s.side_bracket[0]}, right {s.side_bracket[1]}); "
            f"pairwise_min={fmt_frac(s.pairwise_min_gap)}"
        )
        print(f"  residue_multiplicities={s.residue_multiplicities}")
        print(f"  circular_order={s.circular_order}")
        h = s.half_tournament
        print(
            "  "
            f"half_tournament score_seq={h.score_sequence} "
            f"score_hist={dict(h.score_hist)} "
            f"antipodal={h.antipodal_pairs[:6]} "
            f"collisions={h.collision_pairs[:6]}"
        )
        g = s.neighbor_graph
        print(
            "  "
            f"two_nearest mutual_sample={g.mutual_pairs[:8]} "
            f"weak_components={g.weak_component_sizes} "
            f"excess_tie_arcs={g.excess_tie_arcs}"
        )
        print()


def print_generated_routes() -> None:
    print("GENERATED PROOF ROUTES")
    print("=" * 96)
    routes = [
        (
            "Two-neighbor bracket lemma",
            "At any interior maximum of min_i ||v_i t||, the nearest left and nearest right runners bracket 0. If both nearest runners are on one side, a small motion improves the gap. Turn this into a local cut-protection rule.",
        ),
        (
            "Incomplete tournament defect",
            "At tight boundary witnesses, antipodal ties and runner collisions are not noise; they are quotient-compression data. Treat them like missing/identified arcs in an incomplete tournament and ask for a condensation theorem.",
        ),
        (
            "Difference-speed shadow",
            "The full pairwise ledger is an LRC instance for differences v_i-v_j. A counterexample to the original problem should force a highly compressed difference tournament at every would-be witness time.",
        ),
        (
            "Residue polygon lifting",
            "For t=p/q, the configuration is a labelled cyclic tournament on residues mod q with multiplicities. Tight examples collapse to q=n polygons; near-counterexamples should be classified by how far their residue tournament is from a regular polygon.",
        ),
        (
            "Two-nearest core peeling",
            "Build a peel process on the two-nearest graph: delete vertices whose two-neighbor contact is not reciprocated across a zero-bracketing cut. A true open-cover obstruction should leave a nonempty nearest-neighbor core.",
        ),
        (
            "Omega of danger packets",
            "Make vertices be compatible pairs of near-zero runners or endpoint protectors. Independent sets are simultaneously realizable danger packets. This imports the tournament Omega/OCF language without pretending raw LRC is a tournament.",
        ),
    ]
    for name, body in routes:
        print(f"- {name}: {body}")
    print()


def main() -> None:
    examples = (
        ("initial k=4", [1, 2, 3, 4]),
        ("sporadic n=5", [1, 3, 4, 7]),
        ("initial k=5", [1, 2, 3, 4, 5]),
        ("sporadic n=6", [1, 3, 4, 5, 9]),
        ("near-tight k=6", [1, 4, 5, 6, 7, 11]),
        ("n14 seven-ladder", [1, 7, 14, 21, 28, 35, 49, 56, 63, 70, 77, 84, 91]),
        (
            "n14 exported-debt",
            [1] + [14 * q for q in range(1, 14) if q != 6],
        ),
    )
    snapshots = tuple(snapshot(label, speeds) for label, speeds in examples)

    print("LRC neighbor/tournament lift (codex-2026-05-31 S431)")
    print("All times and distances are exact Fractions inherited from S356.\n")
    print_methodology()
    print_snapshot_table(snapshots)
    print_snapshot_details(snapshots)
    print_generated_routes()
    print("Takeaway")
    print(
        "  The lonely scalar is only the height of one vertex, 0, against the "
        "configuration.  The circular half-distance tournament and the "
        "two-nearest graph remember quotient collisions, antipodal ties, "
        "bracketing, and local cores.  These are exactly the kinds of residue "
        "and protection structures where the tournament toolkit has been most "
        "productive."
    )


if __name__ == "__main__":
    main()
