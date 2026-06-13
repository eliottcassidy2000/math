#!/usr/bin/env python3
"""Circle-section and boundary-event tournament functors for LRC.

codex-2026-06-01 S539

This script challenges the implicit assumption that tournament vertices should
be runners.  It tests functors whose vertices are fixed sections of the time
circle, or fixed section boundaries, and whose edges change when runners enter
or leave those sections.

For total LRC denominator n, the q=n section grid makes the open unsafe set
exactly the two observer-adjacent sections:

    [0,1/n) and (1-1/n,1).

Thus LRC can be read as an empty-danger-section exhibition problem.  The
q=2n grid refines this into four danger half-sections and tests whether extra
metric memory gives a better restricted class space.

Tournament Analysis declaration:
    vertices:
        fixed circle sections or fixed section boundaries;
    pairwise observable:
        occupancy count, runner-speed flux, empty/vacuum status, or boundary
        crossing flux;
    switch/gauge:
        lexicographic pressure comparison, with cyclic-section colors marking
        danger sections/boundaries and optional occupancy bits;
    tie Hamiltonian path:
        increasing cyclic section index;
    fingerprints:
        image size, good-only/bad-only/mixed fibers, certification rate,
        score histograms, directed triangles, SCCs, Hamiltonian paths, and a
        meta-tournament over functors.

Stored output:
    05-knowledge/results/lrc_section_boundary_functors_s539.out
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from math import gcd, log2


ONE = Fraction(1, 1)
ZERO = Fraction(0, 1)


Adj = tuple[tuple[int, ...], ...]
Key = tuple[tuple[int, ...], tuple[int, ...]]


BOUNDS = {
    4: 10,
    5: 8,
    6: 7,
}


@dataclass(frozen=True)
class State:
    speeds: tuple[int, ...]
    t: Fraction
    good: bool
    sections: tuple[int, ...]


@dataclass(frozen=True)
class Functor:
    name: str
    q_factor: int
    kind: str


@dataclass
class ClassStat:
    good: int = 0
    bad: int = 0
    exemplar: Adj | None = None
    good_speed_sets: set[tuple[int, ...]] | None = None

    def __post_init__(self) -> None:
        if self.good_speed_sets is None:
            self.good_speed_sets = set()


FUNCTORS = (
    Functor("section_pressure", 1, "sections by occupancy+speed flux"),
    Functor("section_pressure", 2, "half-sections by occupancy+speed flux"),
    Functor("section_empty_colored", 1, "sections colored by danger+occupied"),
    Functor("section_empty_colored", 2, "half-sections colored by danger+occupied"),
    Functor("boundary_flux", 1, "boundaries by incoming section flux"),
    Functor("boundary_flux", 2, "half-boundaries by incoming section flux"),
    Functor("void_pressure", 1, "empty sections beat occupied sections"),
    Functor("void_pressure", 2, "empty half-sections beat occupied half-sections"),
)


def frac(x: Fraction) -> Fraction:
    return x - (x.numerator // x.denominator)


def dist0(x: Fraction) -> Fraction:
    f = frac(x)
    return min(f, ONE - f)


def primitive_speed_sets(total_n: int, max_speed: int) -> list[tuple[int, ...]]:
    out = []
    for moving in combinations(range(1, max_speed + 1), total_n - 1):
        g = 0
        for v in moving:
            g = gcd(g, v)
        if g == 1:
            out.append(moving)
    return out


def good_at(speeds: tuple[int, ...], total_n: int, t: Fraction) -> bool:
    threshold = Fraction(1, total_n)
    return all(dist0(Fraction(v) * t) >= threshold for v in speeds)


def section_of(x: Fraction, q: int) -> int:
    f = frac(x)
    idx = (f * q).numerator // (f * q).denominator
    if idx >= q:
        return q - 1
    return idx


def danger_sections(total_n: int, q: int) -> set[int]:
    threshold = Fraction(1, total_n)
    danger = set()
    for s in range(q):
        left = Fraction(s, q)
        right = Fraction(s + 1, q)
        # Open section intersects [0,threshold) or (1-threshold,1).
        if left < threshold and right > 0:
            danger.add(s)
        if right > ONE - threshold and left < ONE:
            danger.add(s)
    return danger


def danger_boundaries(total_n: int, q: int) -> set[int]:
    points = {ZERO, Fraction(1, total_n), ONE - Fraction(1, total_n)}
    out = set()
    for b in range(q):
        x = Fraction(b, q)
        if any(x == p for p in points):
            out.add(b)
    return out


def event_times(speeds: tuple[int, ...], total_n: int, q_values: tuple[int, ...]) -> list[Fraction]:
    out = {ZERO, ONE}
    threshold = Fraction(1, total_n)
    for v in speeds:
        for q in q_values:
            for k in range(q * v):
                out.add(Fraction(k, q * v))
        for k in range(v):
            out.add(frac((Fraction(k) + threshold) / v))
            out.add(frac((Fraction(k) - threshold) / v))
    return sorted(out)


def sampled_states(speeds: tuple[int, ...], total_n: int, q_values: tuple[int, ...]) -> list[tuple[Fraction, bool]]:
    walls = event_times(speeds, total_n, q_values)
    out = []
    for left, right in zip(walls, walls[1:]):
        if left < right:
            t = (left + right) / 2
            out.append((t, good_at(speeds, total_n, t)))
    return out


def dihedral_perms(q: int) -> list[tuple[int, ...]]:
    perms = []
    for shift in range(q):
        perms.append(tuple((i + shift) % q for i in range(q)))
    for shift in range(q):
        perms.append(tuple((shift - i) % q for i in range(q)))
    return perms


def canonical_cyclic(adj: Adj, colors: tuple[int, ...]) -> Key:
    q = len(colors)
    best: Key | None = None
    for p in dihedral_perms(q):
        pc = tuple(colors[p[i]] for i in range(q))
        bits = tuple(adj[p[i]][p[j]] for i in range(q) for j in range(i + 1, q))
        key = (pc, bits)
        if best is None or key < best:
            best = key
    assert best is not None
    return best


def tournament_from_scores(scores: tuple[tuple, ...]) -> Adj:
    q = len(scores)
    adj = [[0] * q for _ in range(q)]
    for i, j in combinations(range(q), 2):
        if scores[i] > scores[j] or (scores[i] == scores[j] and i < j):
            adj[i][j] = 1
        else:
            adj[j][i] = 1
    return tuple(tuple(row) for row in adj)


def section_vectors(speeds: tuple[int, ...], t: Fraction, q: int) -> tuple[tuple[int, ...], tuple[int, ...]]:
    counts = [0] * q
    flux = [0] * q
    for v in speeds:
        s = section_of(Fraction(v) * t, q)
        counts[s] += 1
        flux[s] += v
    return tuple(counts), tuple(flux)


def section_distance_rank(s: int, q: int) -> int:
    # Distance from the observer boundary at 0 in section steps.
    return min(s, q - s)


def emit(functor: Functor, speeds: tuple[int, ...], total_n: int, t: Fraction) -> tuple[Key, Adj]:
    q = total_n * functor.q_factor
    danger = danger_sections(total_n, q)
    counts, flux = section_vectors(speeds, t, q)

    if functor.name == "section_pressure":
        scores = tuple((counts[s], flux[s], -section_distance_rank(s, q)) for s in range(q))
        colors = tuple(1 if s in danger else 0 for s in range(q))
    elif functor.name == "section_empty_colored":
        scores = tuple((counts[s], flux[s], -section_distance_rank(s, q)) for s in range(q))
        colors = tuple((2 if s in danger else 0) + (1 if counts[s] else 0) for s in range(q))
    elif functor.name == "void_pressure":
        scores = tuple((1 if counts[s] == 0 else 0, -counts[s], -section_distance_rank(s, q)) for s in range(q))
        colors = tuple(1 if s in danger else 0 for s in range(q))
    elif functor.name == "boundary_flux":
        bdanger = danger_boundaries(total_n, q)
        scores = []
        for b in range(q):
            incoming_section = (b - 1) % q
            outgoing_section = b
            scores.append(
                (
                    flux[incoming_section],
                    counts[incoming_section],
                    -flux[outgoing_section],
                    -section_distance_rank(b, q),
                )
            )
        colors = tuple(1 if b in bdanger else 0 for b in range(q))
    else:
        raise ValueError(functor.name)

    adj = tournament_from_scores(tuple(scores))
    return canonical_cyclic(adj, colors), adj


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
    total = 0
    for a, b, c in combinations(range(len(adj)), 3):
        if (adj[a][b] and adj[b][c] and adj[c][a]) or (
            adj[a][c] and adj[c][b] and adj[b][a]
        ):
            total += 1
    return total


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
        total = 1
        for w in rgraph[v]:
            if w not in seen:
                total += rdfs(w)
        return total

    for v in reversed(order):
        if v not in seen:
            sizes.append(rdfs(v))
    return tuple(sorted(sizes, reverse=True))


def score_hist(adj: Adj) -> tuple[tuple[int, int], ...]:
    return tuple(sorted(Counter(sum(row) for row in adj).items()))


def summarize(classes: dict[Key, ClassStat], speed_sets: list[tuple[int, ...]]) -> dict[str, object]:
    pure_good_keys = {key for key, stat in classes.items() if stat.good and not stat.bad}
    certified = 0
    for speeds in speed_sets:
        if any(speeds in (classes[key].good_speed_sets or set()) for key in pure_good_keys):
            certified += 1
    h_values = sorted(
        {
            hamiltonian_paths(stat.exemplar)
            for stat in classes.values()
            if stat.good and stat.exemplar is not None and len(stat.exemplar) <= 8
        }
    )
    c3_values = sorted(
        {
            directed_triangles(stat.exemplar)
            for stat in classes.values()
            if stat.good and stat.exemplar is not None
        }
    )
    scc_values = sorted(
        {
            scc_sizes(stat.exemplar)
            for stat in classes.values()
            if stat.good and stat.exemplar is not None
        }
    )
    return {
        "classes": len(classes),
        "good_classes": sum(1 for stat in classes.values() if stat.good),
        "pure_good": len(pure_good_keys),
        "mixed": sum(1 for stat in classes.values() if stat.good and stat.bad),
        "bad_only": sum(1 for stat in classes.values() if stat.bad and not stat.good),
        "certified": certified,
        "h_values": h_values,
        "c3_values": c3_values,
        "scc_values": scc_values,
    }


def meta_tournament(profile: dict[str, tuple[int, ...]]) -> Adj:
    names = list(profile)
    n = len(names)
    adj = [[0] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        if profile[names[i]] > profile[names[j]] or (
            profile[names[i]] == profile[names[j]] and i < j
        ):
            adj[i][j] = 1
        else:
            adj[j][i] = 1
    return tuple(tuple(row) for row in adj)


def profile_tuple(rows: list[dict[str, object]]) -> tuple[int, ...]:
    speed_total = sum(int(row["speed_sets"]) for row in rows)
    certified = sum(int(row["certified"]) for row in rows)
    good_classes = sum(int(row["good_classes"]) for row in rows)
    pure_good = sum(int(row["pure_good"]) for row in rows)
    classes = sum(int(row["classes"]) for row in rows)
    mixed = sum(int(row["mixed"]) for row in rows)
    cert_rate = round(1000 * certified / speed_total) if speed_total else 0
    purity = round(1000 * pure_good / good_classes) if good_classes else 0
    compression = -classes
    anti_mixed = -mixed
    return (cert_rate, purity, anti_mixed, compression)


def short_list(values: list | tuple, limit: int = 6) -> str:
    if len(values) <= limit:
        return str(values)
    return str(values[:limit])[:-1] + ", ...]"


def main() -> None:
    print("LRC section/boundary tournament functors -- codex-2026-06-01 S539")
    print()
    print("Assumption being challenged: tournament vertices need not be runners.")
    print("Here vertices are fixed circle sections or fixed section boundaries.")
    print()

    all_rows: list[dict[str, object]] = []
    for total_n, max_speed in BOUNDS.items():
        speed_sets = primitive_speed_sets(total_n, max_speed)
        q_values = tuple(sorted({total_n * f.q_factor for f in FUNCTORS}))
        class_maps: dict[str, dict[Key, ClassStat]] = {
            f"{f.name}_q{total_n * f.q_factor}": {} for f in FUNCTORS
        }
        state_count = 0
        good_count = 0

        for speeds in speed_sets:
            states = sampled_states(speeds, total_n, q_values)
            state_count += len(states)
            good_count += sum(1 for _, good in states if good)
            for t, good in states:
                for functor in FUNCTORS:
                    name = f"{functor.name}_q{total_n * functor.q_factor}"
                    key, adj = emit(functor, speeds, total_n, t)
                    stat = class_maps[name].setdefault(key, ClassStat(exemplar=adj))
                    if good:
                        stat.good += 1
                        assert stat.good_speed_sets is not None
                        stat.good_speed_sets.add(speeds)
                    else:
                        stat.bad += 1

        print(f"PART total_n={total_n} max_speed={max_speed}")
        print("=" * 78)
        print(
            f"primitive_speed_sets={len(speed_sets)} open_section_states={state_count} "
            f"good_states={good_count}"
        )
        for functor in FUNCTORS:
            name = f"{functor.name}_q{total_n * functor.q_factor}"
            summary = summarize(class_maps[name], speed_sets)
            row = dict(summary)
            row["name"] = name
            row["speed_sets"] = len(speed_sets)
            all_rows.append(row)
            print(
                "{:<28s} classes={:>5d} good={:>4d} pure={:>4d} mixed={:>4d} "
                "cert={:>4d}/{:<4d} Hgood={} c3good={} SCCgood={}".format(
                    name,
                    int(summary["classes"]),
                    int(summary["good_classes"]),
                    int(summary["pure_good"]),
                    int(summary["mixed"]),
                    int(summary["certified"]),
                    len(speed_sets),
                    short_list(summary["h_values"]),
                    short_list(summary["c3_values"]),
                    short_list(summary["scc_values"], 4),
                )
            )
        print()

    print("PART meta-tournament over section/boundary functors")
    print("=" * 78)
    by_name: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in all_rows:
        by_name[str(row["name"])].append(row)
    profiles = {name: profile_tuple(rows) for name, rows in by_name.items()}
    for name in sorted(profiles):
        print(f"  {name:<28s} profile={profiles[name]}")
    meta = meta_tournament(profiles)
    print(f"score_hist={score_hist(meta)}")
    print(f"directed_triangles={directed_triangles(meta)}")
    print(f"SCCs={scc_sizes(meta)}")
    print(f"Hamiltonian_paths={hamiltonian_paths(meta)}")
    print("dominance_order:")
    names = list(profiles)
    for score, name in sorted(((sum(meta[i]), names[i]) for i in range(len(names))), reverse=True):
        print(f"  score={score} {name}")
    print()

    print("SYNTHESIS")
    print("=" * 78)
    print("1. The section-node idea is exact at q=n in open cells: LRC is")
    print("   emptiness of the two danger sections adjacent to the observer.")
    print("2. Pure class quotients appear when occupancy is remembered as a color")
    print("   on danger sections.  Without that color, pressure and void maps can")
    print("   compress aggressively but may mix safe and unsafe states.")
    print("3. Boundary-flux tournaments are dynamic: their edges change exactly when")
    print("   runners enter or leave fixed section boundaries.  They are less pure")
    print("   alone, but they look like a useful clock derivative for BLEX/near-graph.")
    print("4. q=2n half-sections pay a label tax but separate left/right boundary")
    print("   events.  This is a concrete bridge from section nodes to sentinel maps.")
    print("5. Future agents should challenge the default vertex choice before coding:")
    print("   runners, sections, boundaries, gaps, events, covers, residues, and proof")
    print("   obligations are all possible tournament vertex sets.")


if __name__ == "__main__":
    main()
