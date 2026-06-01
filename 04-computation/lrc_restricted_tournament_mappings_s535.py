#!/usr/bin/env python3
"""Restricted LRC-to-tournament quotient maps.

codex-2026-06-01 S535

The user asked for alternate mappings of LRC into tournament space where the
question is still "which isomorphism classes can be exhibited?", but the class
set is meaningfully smaller than raw A000568.

This script audits several such quotient maps on exact small primitive clocks.
Each map sends sampled LRC states to a colored or uncolored tournament class.
The main metrics are:

* image size: how many classes the clock actually visits;
* target image: how many classes are visited at lonely/source states;
* mixed fibers: classes containing both good and bad states;
* certification: how many primitive speed sets have a sampled state in a
  good-only class;
* compression entropy: class distribution entropy and states-per-class bits.

Tournament Analysis declaration:
    vertices: quotient maps in MAPS
    pairwise observable: aggregate profile
        (certification rate, purity rate, compression bits, target bits,
         negative image size)
    switch/gauge: lexicographic dominance of the profile vector
    tie Hamiltonian path: map declaration order below

Stored output:
    05-knowledge/results/lrc_restricted_tournament_mappings_s535.out
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, permutations
from math import gcd, log2


ONE = Fraction(1, 1)
HALF = Fraction(1, 2)


Adj = tuple[tuple[int, ...], ...]
ClassKey = tuple[tuple[int, ...], tuple[int, ...]]


A000568_TABLE = {
    1: 1,
    2: 1,
    3: 2,
    4: 4,
    5: 12,
    6: 56,
    7: 456,
    8: 6880,
    9: 191536,
}


BOUNDS = {
    4: 12,
    5: 10,
    6: 9,
    7: 7,
}


@dataclass(frozen=True)
class State:
    speeds: tuple[int, ...]
    t: Fraction
    good: bool
    pos: tuple[Fraction, ...]


@dataclass(frozen=True)
class MapSpec:
    name: str
    kind: str
    vertices: str
    threshold_decorated: bool
    conditional_good_only: bool


@dataclass
class ClassStat:
    good: int = 0
    bad: int = 0
    exemplar_adj: Adj | None = None
    speed_sets_good: set[tuple[int, ...]] | None = None

    def __post_init__(self) -> None:
        if self.speed_sets_good is None:
            self.speed_sets_good = set()


MAPS: tuple[MapSpec, ...] = (
    MapSpec(
        "phase_runner",
        "raw circular body",
        "moving runners",
        False,
        False,
    ),
    MapSpec(
        "source_deleted_phase",
        "target menu",
        "moving runners, good states only",
        False,
        True,
    ),
    MapSpec(
        "observer_source_marked",
        "exact source lift",
        "observer plus moving runners",
        True,
        False,
    ),
    MapSpec(
        "gap_threshold_necklace",
        "outside clasp necklace",
        "circular gaps",
        True,
        False,
    ),
    MapSpec(
        "gap_kinetic_flow",
        "THM-387 flow necklace",
        "circular gaps",
        True,
        False,
    ),
    MapSpec(
        "blocker_deficit_shadow",
        "endpoint blocker fiber",
        "moving runners",
        True,
        False,
    ),
    MapSpec(
        "apex_boundary_runner",
        "source-sink/apex marked body",
        "moving runners",
        True,
        False,
    ),
)


def a000568(n: int) -> int | None:
    return A000568_TABLE.get(n)


def mod1(x: Fraction) -> Fraction:
    return x % ONE


def circular_distance(a: Fraction, b: Fraction) -> Fraction:
    delta = mod1(b - a)
    return min(delta, ONE - delta)


def positions(speeds: tuple[int, ...], t: Fraction) -> tuple[Fraction, ...]:
    return tuple(mod1(Fraction(v) * t) for v in speeds)


def primitive_speed_sets(total_n: int, max_speed: int) -> list[tuple[int, ...]]:
    out: list[tuple[int, ...]] = []
    for moving in combinations(range(1, max_speed + 1), total_n - 1):
        g = 0
        for speed in moving:
            g = gcd(g, speed)
        if g == 1:
            out.append((0,) + moving)
    return out


def lrc_good(pos: tuple[Fraction, ...]) -> bool:
    threshold = Fraction(1, len(pos))
    return all(circular_distance(pos[0], pos[i]) >= threshold for i in range(1, len(pos)))


def event_times(speeds: tuple[int, ...]) -> list[Fraction]:
    total_n = len(speeds)
    threshold = Fraction(1, total_n)
    out: set[Fraction] = {Fraction(0), Fraction(1)}

    for speed in speeds[1:]:
        for k in range(speed):
            out.add(mod1((Fraction(k) + threshold) / speed))
            out.add(mod1((Fraction(k) - threshold) / speed))

    for i, j in combinations(range(total_n), 2):
        d = abs(speeds[i] - speeds[j])
        if d == 0:
            continue
        for k in range(2 * d + 1):
            out.add(Fraction(k, 2 * d))

    return sorted(out)


def sampled_states(speeds: tuple[int, ...]) -> list[State]:
    walls = event_times(speeds)
    times: set[Fraction] = set(walls)
    for left, right in zip(walls, walls[1:]):
        if left < right:
            times.add((left + right) / 2)
    out = []
    for t in sorted(times):
        pos = positions(speeds, t)
        out.append(State(speeds=speeds, t=t, good=lrc_good(pos), pos=pos))
    return out


@lru_cache(maxsize=None)
def _perms(n: int) -> tuple[tuple[int, ...], ...]:
    return tuple(permutations(range(n)))


@lru_cache(maxsize=None)
def canonical_colored(adj: Adj, colors: tuple[int, ...]) -> ClassKey:
    n = len(adj)
    best: ClassKey | None = None
    for p in _perms(n):
        pc = tuple(colors[p[i]] for i in range(n))
        bits = tuple(adj[p[i]][p[j]] for i in range(n) for j in range(i + 1, n))
        key = (pc, bits)
        if best is None or key < best:
            best = key
    assert best is not None
    return best


def tournament_from_winner(n: int, winner) -> Adj:
    adj = [[0] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        w = winner(i, j)
        if w == i:
            adj[i][j] = 1
        else:
            adj[j][i] = 1
    return tuple(tuple(row) for row in adj)


def phase_adj_for_vertices(pos: tuple[Fraction, ...], vertices: tuple[int, ...]) -> Adj:
    def winner(i: int, j: int) -> int:
        a = vertices[i]
        b = vertices[j]
        delta = mod1(pos[b] - pos[a])
        if delta == 0 or delta == HALF:
            return i
        return i if delta < HALF else j

    return tournament_from_winner(len(vertices), winner)


def observer_source_adj(state: State) -> tuple[Adj, tuple[int, ...]]:
    n = len(state.speeds)
    threshold = Fraction(1, n)
    adj = [[0] * n for _ in range(n)]
    for i in range(1, n):
        if circular_distance(state.pos[0], state.pos[i]) >= threshold:
            adj[0][i] = 1
        else:
            adj[i][0] = 1
    runner_adj = phase_adj_for_vertices(state.pos, tuple(range(1, n)))
    for a in range(1, n):
        for b in range(1, n):
            if a != b:
                adj[a][b] = runner_adj[a - 1][b - 1]
    colors = (1,) + (0,) * (n - 1)
    return tuple(tuple(row) for row in adj), colors


def gap_data(state: State) -> tuple[list[int], list[Fraction], list[int], list[int]]:
    order = sorted(range(len(state.pos)), key=lambda idx: (state.pos[idx], idx))
    gaps: list[Fraction] = []
    obs_adj: list[int] = []
    derivs: list[int] = []
    for k, idx in enumerate(order):
        nxt = order[(k + 1) % len(order)]
        gaps.append(mod1(state.pos[nxt] - state.pos[idx]))
        obs_adj.append(1 if idx == 0 or nxt == 0 else 0)
        derivs.append(state.speeds[nxt] - state.speeds[idx])
    return order, gaps, obs_adj, derivs


def gap_threshold_map(state: State) -> tuple[Adj, tuple[int, ...]]:
    n = len(state.speeds)
    threshold = Fraction(1, n)
    _, gaps, obs_adj, _ = gap_data(state)

    def winner(i: int, j: int) -> int:
        if gaps[i] != gaps[j]:
            return i if gaps[i] > gaps[j] else j
        return i

    colors = tuple((2 if obs_adj[i] else 0) + (1 if gaps[i] >= threshold else 0) for i in range(n))
    return tournament_from_winner(n, winner), colors


def sign(x: int) -> int:
    if x > 0:
        return 1
    if x < 0:
        return -1
    return 0


def gap_kinetic_map(state: State) -> tuple[Adj, tuple[int, ...]]:
    n = len(state.speeds)
    threshold = Fraction(1, n)
    _, gaps, obs_adj, derivs = gap_data(state)
    potentials = tuple((1 if gaps[i] >= threshold else 0, sign(derivs[i]), gaps[i]) for i in range(n))

    def winner(i: int, j: int) -> int:
        if potentials[i] != potentials[j]:
            return i if potentials[i] > potentials[j] else j
        return i

    colors = tuple((2 if obs_adj[i] else 0) + (1 if gaps[i] >= threshold else 0) for i in range(n))
    return tournament_from_winner(n, winner), colors


def blocker_deficit_map(state: State) -> tuple[Adj, tuple[int, ...]]:
    n = len(state.speeds)
    threshold = Fraction(1, n)
    data = []
    colors = []
    for i in range(1, n):
        p = mod1(state.pos[i] - state.pos[0])
        dist = min(p, ONE - p)
        deficit = max(Fraction(0), threshold - dist)
        if deficit == 0:
            side = 0
        elif p < threshold:
            side = 1
        else:
            side = 2
        colors.append(side)
        data.append((deficit, side, -dist, state.speeds[i]))

    def winner(i: int, j: int) -> int:
        if data[i] != data[j]:
            return i if data[i] > data[j] else j
        return i

    return tournament_from_winner(n - 1, winner), tuple(colors)


def apex_boundary_runner_map(state: State) -> tuple[Adj, tuple[int, ...]]:
    n = len(state.speeds)
    threshold = Fraction(1, n)
    moving = tuple(range(1, n))
    adj = phase_adj_for_vertices(state.pos, moving)
    colors = [0] * (n - 1)
    order = sorted(range(n), key=lambda idx: (state.pos[idx], idx))
    obs_pos = order.index(0)
    left_runner = order[(obs_pos - 1) % n]
    right_runner = order[(obs_pos + 1) % n]
    left_gap = mod1(state.pos[0] - state.pos[left_runner])
    right_gap = mod1(state.pos[right_runner] - state.pos[0])
    if left_runner != 0:
        colors[left_runner - 1] = 2 if left_gap >= threshold else 1
    if right_runner != 0:
        colors[right_runner - 1] = 4 if right_gap >= threshold else 3
    return adj, tuple(colors)


def emit_map(spec: MapSpec, state: State) -> tuple[ClassKey, Adj] | None:
    n = len(state.speeds)
    if spec.name == "phase_runner":
        adj = phase_adj_for_vertices(state.pos, tuple(range(1, n)))
        key = canonical_colored(adj, (0,) * (n - 1))
        return key, adj
    if spec.name == "source_deleted_phase":
        if not state.good:
            return None
        adj = phase_adj_for_vertices(state.pos, tuple(range(1, n)))
        key = canonical_colored(adj, (0,) * (n - 1))
        return key, adj
    if spec.name == "observer_source_marked":
        adj, colors = observer_source_adj(state)
        return canonical_colored(adj, colors), adj
    if spec.name == "gap_threshold_necklace":
        adj, colors = gap_threshold_map(state)
        return canonical_colored(adj, colors), adj
    if spec.name == "gap_kinetic_flow":
        adj, colors = gap_kinetic_map(state)
        return canonical_colored(adj, colors), adj
    if spec.name == "blocker_deficit_shadow":
        adj, colors = blocker_deficit_map(state)
        return canonical_colored(adj, colors), adj
    if spec.name == "apex_boundary_runner":
        adj, colors = apex_boundary_runner_map(state)
        return canonical_colored(adj, colors), adj
    raise ValueError(spec.name)


def hamiltonian_paths(adj: Adj) -> int:
    n = len(adj)
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for mask in range(1 << n):
        for last in range(n):
            cur = dp[mask][last]
            if not cur:
                continue
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                if adj[last][nxt]:
                    dp[mask | (1 << nxt)][nxt] += cur
    return sum(dp[full])


def directed_triangles(adj: Adj) -> int:
    out = 0
    for a, b, c in combinations(range(len(adj)), 3):
        if (adj[a][b] and adj[b][c] and adj[c][a]) or (
            adj[a][c] and adj[c][b] and adj[b][a]
        ):
            out += 1
    return out


def scc_sizes(adj: Adj) -> tuple[int, ...]:
    n = len(adj)
    graph = [[j for j in range(n) if adj[i][j]] for i in range(n)]
    rgraph = [[] for _ in range(n)]
    for i, row in enumerate(graph):
        for j in row:
            rgraph[j].append(i)

    seen: set[int] = set()
    order: list[int] = []

    def dfs(v: int) -> None:
        seen.add(v)
        for w in graph[v]:
            if w not in seen:
                dfs(w)
        order.append(v)

    for v in range(n):
        if v not in seen:
            dfs(v)

    seen.clear()
    sizes = []

    def rdfs(v: int) -> int:
        seen.add(v)
        size = 1
        for w in rgraph[v]:
            if w not in seen:
                size += rdfs(w)
        return size

    for v in reversed(order):
        if v not in seen:
            sizes.append(rdfs(v))
    return tuple(sorted(sizes, reverse=True))


def score_hist(adj: Adj) -> tuple[tuple[int, int], ...]:
    return tuple(sorted(Counter(sum(row) for row in adj).items()))


def entropy(counts: list[int]) -> float:
    total = sum(counts)
    if total == 0:
        return 0.0
    out = 0.0
    for count in counts:
        p = count / total
        out -= p * log2(p)
    return out


def summarize_stats(
    total_n: int,
    spec: MapSpec,
    class_stats: dict[ClassKey, ClassStat],
    speed_sets: list[tuple[int, ...]],
    state_count: int,
) -> dict[str, object]:
    classes = len(class_stats)
    good_classes = sum(1 for stat in class_stats.values() if stat.good)
    pure_good = sum(1 for stat in class_stats.values() if stat.good and not stat.bad)
    mixed = sum(1 for stat in class_stats.values() if stat.good and stat.bad)
    bad_only = sum(1 for stat in class_stats.values() if stat.bad and not stat.good)
    pure_good_keys = {key for key, stat in class_stats.items() if stat.good and not stat.bad}
    certified = 0
    if not spec.conditional_good_only:
        for speeds in speed_sets:
            if any(speeds in class_stats[key].speed_sets_good for key in pure_good_keys):
                certified += 1
    else:
        certified = len({sp for stat in class_stats.values() for sp in (stat.speed_sets_good or set())})
    counts = [stat.good + stat.bad for stat in class_stats.values()]
    ent = entropy(counts)
    compression_bits = log2(state_count / classes) if classes else 0.0
    good_h_values = sorted(
        {
            hamiltonian_paths(stat.exemplar_adj)
            for stat in class_stats.values()
            if stat.good and stat.exemplar_adj is not None
        }
    )
    ambient = a000568(len(next(iter(class_stats.values())).exemplar_adj)) if classes else None
    ambient_bits = None
    if ambient and classes:
        ambient_bits = log2(ambient / classes)
    return {
        "total_n": total_n,
        "map": spec.name,
        "classes": classes,
        "good_classes": good_classes,
        "pure_good": pure_good,
        "mixed": mixed,
        "bad_only": bad_only,
        "certified": certified,
        "speed_sets": len(speed_sets),
        "entropy": ent,
        "compression_bits": compression_bits,
        "ambient": ambient,
        "ambient_bits": ambient_bits,
        "good_h_values": good_h_values,
    }


def bool_tournament_from_profile(profile: dict[str, tuple[int, ...]]) -> Adj:
    names = list(profile)
    n = len(names)
    adj = [[0] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        a = profile[names[i]]
        b = profile[names[j]]
        if a > b or (a == b and i < j):
            adj[i][j] = 1
        else:
            adj[j][i] = 1
    return tuple(tuple(row) for row in adj)


def profile_tuple(rows: list[dict[str, object]]) -> tuple[int, ...]:
    speed_total = sum(int(row["speed_sets"]) for row in rows)
    certified = sum(int(row["certified"]) for row in rows)
    classes = sum(int(row["classes"]) for row in rows)
    pure_good = sum(int(row["pure_good"]) for row in rows)
    mixed = sum(int(row["mixed"]) for row in rows)
    good_classes = sum(int(row["good_classes"]) for row in rows)
    cert_rate = round(1000 * certified / speed_total) if speed_total else 0
    purity = round(1000 * pure_good / good_classes) if good_classes else 0
    compression = round(1000 * sum(float(row["compression_bits"]) for row in rows))
    target_bits = round(1000 * sum(log2(max(1, int(row["good_classes"]))) for row in rows))
    # Higher is better: certify, stay pure, compress, keep target image small, keep total image small.
    return (cert_rate, purity, compression, -target_bits, -classes, -mixed)


def main() -> None:
    print("LRC restricted tournament mappings -- codex-2026-06-01 S535")
    print()
    print("Audit bounds:")
    for total_n, max_speed in BOUNDS.items():
        print(f"  total_n={total_n} max_speed={max_speed}")
    print()

    all_rows: list[dict[str, object]] = []

    for total_n, max_speed in BOUNDS.items():
        speed_sets = primitive_speed_sets(total_n, max_speed)
        map_stats: dict[str, dict[ClassKey, ClassStat]] = {spec.name: {} for spec in MAPS}
        state_count = 0
        good_state_count = 0
        for speeds in speed_sets:
            states = sampled_states(speeds)
            state_count += len(states)
            good_state_count += sum(1 for state in states if state.good)
            for state in states:
                for spec in MAPS:
                    emitted = emit_map(spec, state)
                    if emitted is None:
                        continue
                    key, adj = emitted
                    stat = map_stats[spec.name].setdefault(key, ClassStat(exemplar_adj=adj))
                    if state.good:
                        stat.good += 1
                        assert stat.speed_sets_good is not None
                        stat.speed_sets_good.add(speeds)
                    else:
                        stat.bad += 1

        print(f"PART total_n={total_n}")
        print("=" * 78)
        print(
            f"primitive_speed_sets={len(speed_sets)} sampled_states={state_count} "
            f"good_states={good_state_count}"
        )
        for spec in MAPS:
            row = summarize_stats(total_n, spec, map_stats[spec.name], speed_sets, state_count)
            all_rows.append(row)
            ambient_bits = row["ambient_bits"]
            ambient_text = "-" if ambient_bits is None else f"{float(ambient_bits):+.3f}"
            h_values = row["good_h_values"]
            h_text = ",".join(str(x) for x in h_values[:8])
            if len(h_values) > 8:
                h_text += ",..."
            print(
                "{:<25s} kind={:<24s} classes={:>5d} good={:>4d} pure={:>4d} "
                "mixed={:>4d} cert={:>4d}/{:<4d} Hgood=[{}] "
                "entropy={:>7.3f} comp_bits={:>6.3f} ambient_bits={:>8s}".format(
                    spec.name,
                    spec.kind,
                    int(row["classes"]),
                    int(row["good_classes"]),
                    int(row["pure_good"]),
                    int(row["mixed"]),
                    int(row["certified"]),
                    int(row["speed_sets"]),
                    h_text,
                    float(row["entropy"]),
                    float(row["compression_bits"]),
                    ambient_text,
                )
            )
        print()

    print("PART meta-tournament on quotient maps")
    print("=" * 78)
    rows_by_map: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in all_rows:
        rows_by_map[str(row["map"])].append(row)
    profiles = {name: profile_tuple(rows) for name, rows in rows_by_map.items()}
    for name in [spec.name for spec in MAPS]:
        print(f"  {name:<25s} profile={profiles[name]}")
    meta_adj = bool_tournament_from_profile(profiles)
    print(f"score_hist={score_hist(meta_adj)}")
    print(f"directed_triangles={directed_triangles(meta_adj)}")
    print(f"SCCs={scc_sizes(meta_adj)}")
    print(f"Hamiltonian_paths={hamiltonian_paths(meta_adj)}")
    print("dominance_order:")
    for score, name in sorted(((sum(meta_adj[i]), name) for i, name in enumerate(profiles)), reverse=True):
        print(f"  score={score} {name}")
    print()

    print("PART metric glossary")
    print("=" * 78)
    print("image_size: number of colored/uncolored tournament isomorphism classes visited.")
    print("target_image: classes visited at lonely/source states.")
    print("mixed_fiber: one class contains both lonely and non-lonely states.")
    print("certification_rate: speed sets with a sampled state in a good-only class.")
    print("compression_bits: log2(sampled_states / image_size), a coarse quotient strength.")
    print("ambient_bits: log2(A000568(vertex_count) / image_size), only a rough uncolored base gauge.")
    print("label_tax: colored maps can have negative ambient_bits; that is not failure, it is the")
    print("  price paid to make the source/threshold predicate class-local.")
    print()

    print("SYNTHESIS")
    print("=" * 78)
    print("1. Raw phase classes are still too coarse: they mix good and bad states.")
    print("2. The source-deleted map is the clean target-menu language: it asks which")
    print("   runner-subtournament classes can be exhibited when the observer is a source.")
    print("3. Threshold-decorated maps trade a label tax for purity: observer-source,")
    print("   gap-necklace, blocker-shadow, and apex-boundary classes make goodness")
    print("   local to the class in these bounded audits.")
    print("4. The most meaningful restriction is not one quotient but a stack:")
    print("   circular body -> source target -> gap/apex threshold fiber -> kinetic")
    print("   flow fiber -> endpoint blocker fiber.")
    print("5. A future proof should show that every primitive arithmetic clock reaches")
    print("   the source target or a compactified boundary target inside this stack,")
    print("   while a counterexample would have to be a closed walk avoiding all pure")
    print("   target fibers and carrying nonempty labelled endpoint pressure.")


if __name__ == "__main__":
    main()
