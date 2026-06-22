#!/usr/bin/env python3
"""
Octahedral current realizability for the LRC(14) repeated-packet graph.

S101 found an exact finite local balance on the HYP-2632 repeated-residue
packets.  S102 showed that the balance does not lift for free to actual
integer reciprocal tails: a nonzero local divergence remains.

This scout tests a more rigid realizability object.  The nonzero finite packet
edges are not an arbitrary graph and not naturally a tournament: they are
K_6 minus the affine-zero perfect matching, i.e. the octahedron graph, i.e.
the line graph L(K_4).  Integer lifts are therefore height 0-cochains on this
octahedron.  The possible lifted currents are much smaller than all abstract
signed graph currents.

The computation below:

1. verifies the L(K_4) / octahedron identification;
2. enumerates all layer cochains V -> {0,1,2};
3. measures the resulting reciprocal-lift divergence at H=10;
4. compares it with wall counts and octahedral stretch invariants.

This is evidence and proof-target sharpening, not a proof of LRC(14).
"""
from __future__ import annotations

import itertools
import math
import sys
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "04-computation"))

import lrc14_support6_residue_cusp_codex_s12 as s12  # noqa: E402

AMBIENT_D = 9
HMAX = 10
VERTICES = (0, 2, 3, 4, 5, 6)
CORE = (1, 8, 15, 22)
LOW_HEIGHT = 2

# The affine-zero pairs are the three opposite edge-pairs in K_4.
K4_EDGE_OF_RESIDUE = {
    0: (1, 2),
    2: (3, 4),
    3: (1, 3),
    6: (2, 4),
    4: (1, 4),
    5: (2, 3),
}


def chi7(x: int) -> int:
    x %= 7
    if x == 0:
        return 0
    return 1 if x in {1, 2, 4} else -1


def q_selector(a: int, b: int) -> int:
    return (a * b * (1 + 3 * ((a + b) % 7)) - 1) % 7


def loop_weight(residue: int) -> int:
    if residue == 0:
        return -4
    return (-43 - 7 * chi7(residue)) // 2


def edge_weight(a: int, b: int) -> int:
    if (a + b) % 7 == 2:
        return 0
    return 8 if chi7(q_selector(a, b)) == 1 else 1


def lift(residue: int, layer: int) -> int:
    if residue == 0:
        return 7 * (layer + 1)
    return residue + 7 * layer


def layer_of(layers: tuple[int, ...], residue: int) -> int:
    return layers[VERTICES.index(residue)]


def loop_support(layers: tuple[int, ...], residue: int) -> tuple[int, ...]:
    layer = layer_of(layers, residue)
    a = lift(residue, layer)
    b = lift(residue, layer + 1)
    return tuple(sorted(CORE + (a, b)))


def edge_support(layers: tuple[int, ...], a: int, b: int) -> tuple[int, ...]:
    x = lift(a, layer_of(layers, a))
    y = lift(b, layer_of(layers, b))
    return tuple(sorted(CORE + (x, y)))


def nonzero_edges() -> list[tuple[int, int]]:
    out = []
    for a, b in itertools.combinations(VERTICES, 2):
        if edge_weight(a, b):
            out.append((a, b))
    return out


def zero_edges() -> list[tuple[int, int]]:
    out = []
    for a, b in itertools.combinations(VERTICES, 2):
        if edge_weight(a, b) == 0:
            out.append((a, b))
    return out


def section(title: str) -> None:
    print("\n" + "=" * 92)
    print(title)
    print("=" * 92)


@lru_cache(None)
def packet_value(support: tuple[int, ...]) -> float:
    by_shell, _ = s12.six_support_shells(support, AMBIENT_D, HMAX)
    rows = s12.cumulative_rows(by_shell, (HMAX,))
    if not rows:
        return 0.0
    return float(rows[-1][2].real)


@lru_cache(None)
def low_height_relation_count(support: tuple[int, ...], height: int = LOW_HEIGHT) -> int:
    count = 0
    coeffs = range(-height, height + 1)

    def rec(i: int, dot: int, nonzero: bool) -> None:
        nonlocal count
        if i == len(support):
            if nonzero and dot == 0:
                count += 1
            return
        for c in coeffs:
            rec(i + 1, dot + c * support[i], nonzero or c != 0)

    rec(0, 0, False)
    return count // 2


def k4_adjacent(a: int, b: int) -> bool:
    ea = set(K4_EDGE_OF_RESIDUE[a])
    eb = set(K4_EDGE_OF_RESIDUE[b])
    return bool(ea & eb)


def k4_triangle_type(tri: tuple[int, int, int]) -> str:
    edges = [set(K4_EDGE_OF_RESIDUE[v]) for v in tri]
    common = set.intersection(*edges)
    if common:
        return f"star@{next(iter(common))}"
    used_vertices = set().union(*edges)
    if len(used_vertices) == 3:
        return "face"
    return "other"


def verify_octahedron() -> None:
    section("OCTAHEDRAL PACKET GRAPH")
    print("Residue vertex -> K4 edge:")
    for v in VERTICES:
        print(f"  {v}: {K4_EDGE_OF_RESIDUE[v]}")
    print(f"\naffine-zero matching: {zero_edges()}")
    print(f"nonzero packet edges: {nonzero_edges()}")

    mismatches = []
    for a, b in itertools.combinations(VERTICES, 2):
        packet_nonzero = edge_weight(a, b) != 0
        line_graph_nonzero = k4_adjacent(a, b)
        if packet_nonzero != line_graph_nonzero:
            mismatches.append((a, b, packet_nonzero, line_graph_nonzero))
    print(f"\nL(K4) incidence mismatches: {mismatches}")
    print("So the finite packet graph is the octahedron graph L(K4), not a free tournament.")

    print("\nOctahedron triangles (4 K4 vertex-stars + 4 K4 faces):")
    print(f"{'triangle':>14} {'type':>8} {'edge weights':>18} {'high count':>10}")
    for tri in itertools.combinations(VERTICES, 3):
        if not all(edge_weight(*sorted(e)) for e in itertools.combinations(tri, 2)):
            continue
        weights = [edge_weight(*sorted(e)) for e in itertools.combinations(tri, 2)]
        print(
            f"{str(tri):>14} {k4_triangle_type(tri):>8} "
            f"{str(weights):>18} {sum(1 for w in weights if w == 8):>10}"
        )


@dataclass(frozen=True)
class GaugeScore:
    layers: tuple[int, ...]
    l1_div: float
    max_div: float
    net_div: float
    positive_vertices: int
    negative_vertices: int
    layer_sum: int
    diameter: int
    zero_stretch: int
    edge_stretch: int
    weighted_edge_stretch: int
    lap_l1: int
    wall_incidence: int
    wall_max: int

    @property
    def word(self) -> str:
        return "".join(str(x) for x in self.layers)


def divergence(layers: tuple[int, ...]) -> dict[int, float]:
    divs: dict[int, float] = {}
    for v in VERTICES:
        total = packet_value(loop_support(layers, v))
        for u in VERTICES:
            if u == v:
                continue
            a, b = sorted((u, v))
            if edge_weight(a, b):
                total += packet_value(edge_support(layers, a, b))
        divs[v] = total
    return divs


def finite_weighted_laplacian_l1(layers: tuple[int, ...]) -> int:
    total = 0
    for v in VERTICES:
        lv = layer_of(layers, v)
        loc = 0
        for u in VERTICES:
            if u == v:
                continue
            a, b = sorted((u, v))
            w = edge_weight(a, b)
            if w:
                loc += w * (lv - layer_of(layers, u))
        total += abs(loc)
    return total


def wall_incidence(layers: tuple[int, ...]) -> tuple[int, int]:
    total = 0
    wall_max = 0
    for v in VERTICES:
        walls = low_height_relation_count(loop_support(layers, v))
        total += walls
        wall_max = max(wall_max, walls)
    for a, b in nonzero_edges():
        walls = low_height_relation_count(edge_support(layers, a, b))
        total += 2 * walls
        wall_max = max(wall_max, walls)
    return total, wall_max


def score_gauge(layers: tuple[int, ...]) -> GaugeScore:
    divs = divergence(layers)
    vals = list(divs.values())
    zs = zero_edges()
    es = nonzero_edges()
    wall_total, wall_max = wall_incidence(layers)
    return GaugeScore(
        layers=layers,
        l1_div=sum(abs(x) for x in vals),
        max_div=max(abs(x) for x in vals),
        net_div=sum(vals),
        positive_vertices=sum(1 for x in vals if x > 0),
        negative_vertices=sum(1 for x in vals if x < 0),
        layer_sum=sum(layers),
        diameter=max(layers) - min(layers),
        zero_stretch=sum(abs(layer_of(layers, a) - layer_of(layers, b)) for a, b in zs),
        edge_stretch=sum(abs(layer_of(layers, a) - layer_of(layers, b)) for a, b in es),
        weighted_edge_stretch=sum(
            edge_weight(a, b) * abs(layer_of(layers, a) - layer_of(layers, b))
            for a, b in es
        ),
        lap_l1=finite_weighted_laplacian_l1(layers),
        wall_incidence=wall_total,
        wall_max=wall_max,
    )


def corr(xs: list[float], ys: list[float]) -> float:
    mx = sum(xs) / len(xs)
    my = sum(ys) / len(ys)
    vx = sum((x - mx) ** 2 for x in xs)
    vy = sum((y - my) ** 2 for y in ys)
    if vx == 0 or vy == 0:
        return float("nan")
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / math.sqrt(vx * vy)


def all_scores() -> list[GaugeScore]:
    return [score_gauge(layers) for layers in itertools.product((0, 1, 2), repeat=len(VERTICES))]


def print_score_table(title: str, rows: list[GaugeScore], limit: int = 12) -> None:
    print(f"\n{title}")
    print(
        f"{'word':>8} {'L1div':>11} {'maxdiv':>11} {'net':>11} "
        f"{'+/-':>5} {'sum':>3} {'diam':>4} {'zstr':>4} {'estr':>4} "
        f"{'wstr':>5} {'lap':>5} {'walls':>6} {'wmax':>5}"
    )
    for s in rows[:limit]:
        signs = f"{s.positive_vertices}/{s.negative_vertices}"
        print(
            f"{s.word:>8} {s.l1_div:>11.6g} {s.max_div:>11.6g} {s.net_div:>11.6g} "
            f"{signs:>5} {s.layer_sum:>3} {s.diameter:>4} {s.zero_stretch:>4} "
            f"{s.edge_stretch:>4} {s.weighted_edge_stretch:>5} {s.lap_l1:>5} "
            f"{s.wall_incidence:>6} {s.wall_max:>5}"
        )


def gauge_scan() -> None:
    section("REALIZABLE LAYER-COCHAIN SCAN")
    scores = all_scores()
    print(f"Enumerated {len(scores)} layer cochains V(L(K4)) -> {{0,1,2}} at H={HMAX}.")
    print(f"unique reciprocal packet supports evaluated: {packet_value.cache_info().currsize}")
    print(f"unique low-height wall supports evaluated: {low_height_relation_count.cache_info().currsize}")

    print_score_table("Best gauges by L1 divergence", sorted(scores, key=lambda s: s.l1_div))
    print_score_table("Worst gauges by L1 divergence", sorted(scores, key=lambda s: -s.l1_div), 8)

    constants = [s for s in scores if len(set(s.layers)) == 1]
    print_score_table("Constant-layer gauges (S102 start/raised plus one more)", constants, 3)

    print("\nBest gauge at each layer-sum:")
    print(f"{'sum':>3} {'word':>8} {'L1div':>11} {'maxdiv':>11} {'zero':>5} {'wall':>6}")
    for layer_sum in sorted({s.layer_sum for s in scores}):
        row = min((s for s in scores if s.layer_sum == layer_sum), key=lambda s: s.l1_div)
        print(
            f"{layer_sum:>3} {row.word:>8} {row.l1_div:>11.6g} "
            f"{row.max_div:>11.6g} {row.zero_stretch:>5} {row.wall_incidence:>6}"
        )

    section("CORRELATIONS WITH L1 DIVERGENCE")
    metrics = [
        ("layer_sum", [s.layer_sum for s in scores]),
        ("diameter", [s.diameter for s in scores]),
        ("zero_stretch", [s.zero_stretch for s in scores]),
        ("edge_stretch", [s.edge_stretch for s in scores]),
        ("weighted_edge_stretch", [s.weighted_edge_stretch for s in scores]),
        ("weighted_laplacian_L1", [s.lap_l1 for s in scores]),
        ("wall_incidence", [s.wall_incidence for s in scores]),
        ("wall_max", [s.wall_max for s in scores]),
        ("negative_vertices", [s.negative_vertices for s in scores]),
    ]
    ys = [s.l1_div for s in scores]
    print(f"{'metric':>24} {'corr(all)':>12} {'mean within-sum corr':>22}")
    for name, vals in metrics:
        within = []
        for layer_sum in sorted({s.layer_sum for s in scores}):
            idx = [i for i, s in enumerate(scores) if s.layer_sum == layer_sum]
            if len(idx) < 3:
                continue
            c = corr([vals[i] for i in idx], [ys[i] for i in idx])
            if not math.isnan(c):
                within.append(c)
        mean_within = sum(within) / len(within) if within else float("nan")
        print(f"{name:>24} {corr([float(v) for v in vals], ys):>12.6g} {mean_within:>22.6g}")

    section("S103 REALIZABILITY READING")
    print(
        "The useful object is not an arbitrary tournament orientation.  The six "
        "finite packet vertices are the six edges of K4, the affine-zero lane is "
        "the opposite-edge matching, and nonzero packet coupling is octahedral "
        "incidence.  Integer lifting chooses a height 0-cochain on this "
        "octahedron; the reciprocal defect is the divergence of the lifted "
        "current after this cochain is applied."
    )
    print(
        "This suggests the missing LRC structure is an octahedral Hodge problem: "
        "finite Kirchhoff balance on L(K4), coherent curl on the eight triangular "
        "faces routed to low-height walls, and the remaining spread current "
        "routed to HYP-2636 additive-frequency Abel summation."
    )


def main() -> None:
    section("LRC14 OCTAHEDRAL CURRENT REALIZABILITY - CODEX S103")
    verify_octahedron()
    gauge_scan()


if __name__ == "__main__":
    main()
