#!/usr/bin/env python3
"""
LRC section check-off and switch atlas.

Prompt: use regions of the loop as the primary object instead of runners.

For an n-runner speed set in the canonical slowest-runner gauge

    V = (0, v_1, ..., v_{n-1}),

divide the loop into n fixed sections [s/n,(s+1)/n).  For each runner r,
collect the sections in which r has a lonely witness:

    ||(v_r-v_j)t|| >= 1/n for every j != r.

This gives a bipartite "check-off" graph

    runners  <-->  sections.

The dream certificate is a perfect matching: every runner can be assigned its
own section.  Exotic cases are not discarded; they are reduced to Hall-deficit
packets and section-switch graphs.

Gauge warning:
  Absolute sections are not Galilean invariant under adding the same speed to
  every runner.  The slowest-runner gauge is a canonical section frame for this
  experiment.  The invariant observer-relative version is S539's section
  functor.  This script is about the user's "runner uses a quadrant/section"
  picture, so the gauge issue is recorded explicitly rather than hidden.

Tournament Analysis contract:
  * Vertices: fixed loop sections, not runners.
  * Pairwise observable: section support vector
      (number of witnessed runners, total open witness measure, private runners).
  * Switch/gauge: orient section a -> b when a has the larger support vector;
    ties use cyclic section order as the Hamiltonian path.
  * Tie Hamiltonian path: increasing section index.
  * Fingerprints: score histogram, directed 3-cycles, SCCs, Hamiltonian path
    count, plus edge flips under an alternate "private-first" switch.

Assumption challenge:
  Considered vertices = runners, fixed sections, section boundaries, lonely
  events, runner-section edges, Hall-deficit packets, cover arcs, residues, and
  proof obligations.  Chosen primary vertices are fixed sections.  This quotient
  preserves the predicate "which regions can host a runner's lonely moment" and
  the Hall obstruction to checking off all sections, but it destroys exact
  timing, pairwise runner identity away from shared section support, and the
  observer-relative invariance unless the gauge is fixed.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import combinations
from math import gcd


@dataclass(frozen=True)
class Checkoff:
    speeds: tuple[int, ...]
    open_measure: tuple[tuple[F, ...], ...]
    wall_edges: tuple[tuple[int, ...], ...]


def frac(x: F) -> F:
    return x - (x.numerator // x.denominator)


def dist0(x: F) -> F:
    y = frac(x)
    return min(y, 1 - y)


def gcd_all(values):
    g = 0
    for value in values:
        g = gcd(g, value)
    return g


def speed_sets(n: int, max_speed: int) -> list[tuple[int, ...]]:
    out = []
    for tail in combinations(range(1, max_speed + 1), n - 1):
        if gcd_all(tail) == 1:
            out.append((0,) + tail)
    return out


def event_times(speeds: tuple[int, ...], n: int) -> list[F]:
    out = {F(0), F(1)}
    threshold = F(1, n)

    for v in speeds:
        if v:
            for k in range(n * v + 1):
                out.add(F(k, n * v))

    for a, b in combinations(speeds, 2):
        d = abs(a - b)
        if d == 0:
            continue
        for k in range(d):
            out.add(frac((F(k) + threshold) / d))
            out.add(frac((F(k) - threshold) / d))
    return sorted(out)


def is_lonely_runner(speeds: tuple[int, ...], n: int, r: int, t: F) -> bool:
    threshold = F(1, n)
    vr = speeds[r]
    for j, vj in enumerate(speeds):
        if j == r:
            continue
        if dist0(F(vr - vj) * t) < threshold:
            return False
    return True


def section_set(x: F, n: int) -> tuple[int, ...]:
    y = frac(x)
    z = y * n
    if z.denominator == 1:
        boundary = z.numerator % n
        return tuple(sorted({boundary, (boundary - 1) % n}))
    idx = z.numerator // z.denominator
    return (idx,)


def section_of_open(x: F, n: int) -> int:
    y = frac(x) * n
    return min(n - 1, y.numerator // y.denominator)


def checkoff_for(speeds: tuple[int, ...], n: int) -> Checkoff:
    walls = event_times(speeds, n)
    open_measure = [[F(0) for _ in range(n)] for _ in range(n)]
    wall_edges: list[set[int]] = [set() for _ in range(n)]

    for left, right in zip(walls, walls[1:]):
        if left == right:
            continue
        mid = (left + right) / 2
        length = right - left
        for r, v in enumerate(speeds):
            if is_lonely_runner(speeds, n, r, mid):
                s = section_of_open(F(v) * mid, n)
                open_measure[r][s] += length

    for t in walls:
        for r, v in enumerate(speeds):
            if is_lonely_runner(speeds, n, r, t):
                for s in section_set(F(v) * t, n):
                    wall_edges[r].add(s)

    return Checkoff(
        speeds=speeds,
        open_measure=tuple(tuple(row) for row in open_measure),
        wall_edges=tuple(tuple(sorted(row)) for row in wall_edges),
    )


def edge_sets(chk: Checkoff, include_walls: bool = True) -> tuple[set[int], ...]:
    edges = []
    for r, row in enumerate(chk.open_measure):
        sset = {s for s, value in enumerate(row) if value > 0}
        if include_walls:
            sset.update(chk.wall_edges[r])
        edges.append(sset)
    return tuple(edges)


def has_perfect_matching(edges: tuple[set[int], ...], n: int) -> bool:
    match_to_runner = [-1] * n

    def augment(r: int, seen: set[int]) -> bool:
        for s in sorted(edges[r]):
            if s in seen:
                continue
            seen.add(s)
            if match_to_runner[s] == -1 or augment(match_to_runner[s], seen):
                match_to_runner[s] = r
                return True
        return False

    for r in range(n):
        if not augment(r, set()):
            return False
    return True


def hall_deficits(edges: tuple[set[int], ...], n: int):
    deficits = []
    runners = range(n)
    for size in range(1, n + 1):
        for subset in combinations(runners, size):
            union = set()
            for r in subset:
                union.update(edges[r])
            deficit = size - len(union)
            if deficit > 0:
                deficits.append((deficit, subset, tuple(sorted(union))))
    return sorted(deficits, key=lambda item: (-item[0], len(item[1]), item[1]))


def section_support(chk: Checkoff, include_walls: bool = True):
    n = len(chk.speeds)
    edges = edge_sets(chk, include_walls)
    runners_by_section = [set() for _ in range(n)]
    private_by_section = [set() for _ in range(n)]
    measure_by_section = [F(0) for _ in range(n)]

    section_count_by_runner = [len(edges[r]) for r in range(n)]
    for r in range(n):
        for s in edges[r]:
            runners_by_section[s].add(r)
        for s, value in enumerate(chk.open_measure[r]):
            measure_by_section[s] += value
    for r, sections in enumerate(edges):
        if len(sections) == 1:
            private_by_section[next(iter(sections))].add(r)
        else:
            for s in sections:
                if sum(1 for rr in runners_by_section[s] if s in edges[rr]) == 1:
                    private_by_section[s].add(r)
    return {
        "edges": edges,
        "runners_by_section": tuple(tuple(sorted(x)) for x in runners_by_section),
        "private_by_section": tuple(tuple(sorted(x)) for x in private_by_section),
        "measure_by_section": tuple(measure_by_section),
        "section_count_by_runner": tuple(section_count_by_runner),
    }


def tournament_from_scores(scores):
    n = len(scores)
    adj = [[False] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        if scores[i] > scores[j] or (scores[i] == scores[j] and i < j):
            adj[i][j] = True
        else:
            adj[j][i] = True
    return adj


def tournament_scores(adj):
    return [sum(row) for row in adj]


def directed_triangles(adj):
    n = len(adj)
    cycles = 0
    transitive = 0
    for i, j, k in combinations(range(n), 3):
        out = []
        for a in (i, j, k):
            out.append(sum(1 for b in (i, j, k) if a != b and adj[a][b]))
        if sorted(out) == [1, 1, 1]:
            cycles += 1
        else:
            transitive += 1
    return cycles, transitive


def sccs(adj):
    n = len(adj)
    graph = [[j for j in range(n) if adj[i][j]] for i in range(n)]
    rev = [[i for i in range(n) if adj[i][j]] for j in range(n)]
    seen = set()
    order = []

    def dfs(v):
        seen.add(v)
        for u in graph[v]:
            if u not in seen:
                dfs(u)
        order.append(v)

    for v in range(n):
        if v not in seen:
            dfs(v)
    seen.clear()
    comps = []

    def rdfs(v, comp):
        seen.add(v)
        comp.append(v)
        for u in rev[v]:
            if u not in seen:
                rdfs(u, comp)

    for v in reversed(order):
        if v not in seen:
            comp = []
            rdfs(v, comp)
            comps.append(comp)
    return comps


def hamiltonian_path_count(adj):
    n = len(adj)
    dp = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for mask in range(1 << n):
        for last in range(n):
            val = dp[mask][last]
            if not val:
                continue
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                if adj[last][nxt]:
                    dp[mask | (1 << nxt)][nxt] += val
    return sum(dp[-1])


def edge_flips(a, b):
    n = len(a)
    total = 0
    for i, j in combinations(range(n), 2):
        if a[i][j] != b[i][j]:
            total += 1
    return total


def section_tournament(chk: Checkoff):
    support = section_support(chk)
    n = len(chk.speeds)
    scores = []
    private_first = []
    for s in range(n):
        runners = support["runners_by_section"][s]
        private = support["private_by_section"][s]
        measure = support["measure_by_section"][s]
        scores.append((len(runners), measure, len(private), -s))
        private_first.append((len(private), len(runners), measure, -s))
    adj = tournament_from_scores(tuple(scores))
    alt = tournament_from_scores(tuple(private_first))
    return adj, alt, support


def fingerprint(adj):
    cyc, trans = directed_triangles(adj)
    return {
        "scores": tournament_scores(adj),
        "score_hist": dict(sorted(Counter(tournament_scores(adj)).items())),
        "directed_cycles": cyc,
        "transitive_triples": trans,
        "scc_sizes": sorted((len(c) for c in sccs(adj)), reverse=True),
        "hamiltonian_paths": hamiltonian_path_count(adj),
    }


def fmt(fr):
    return f"{fr} ({float(fr):.6f})"


def summarize_scan(n: int, max_speed: int):
    sets = speed_sets(n, max_speed)
    rows = []
    hall_counter = Counter()
    open_hall_counter = Counter()
    for speeds in sets:
        chk = checkoff_for(speeds, n)
        edges = edge_sets(chk, include_walls=True)
        open_edges = edge_sets(chk, include_walls=False)
        matching = has_perfect_matching(edges, n)
        open_matching = has_perfect_matching(open_edges, n)
        section_cover = len(set().union(*edges)) if edges else 0
        open_section_cover = len(set().union(*open_edges)) if open_edges else 0
        all_runners = all(edges[r] for r in range(n))
        deficits = hall_deficits(edges, n)
        open_deficits = hall_deficits(open_edges, n)
        if deficits:
            d, subset, union = deficits[0]
            hall_counter[(d, len(subset), len(union))] += 1
        if open_deficits:
            d, subset, union = open_deficits[0]
            open_hall_counter[(d, len(subset), len(union))] += 1
        wall_debt_sections = tuple(sorted(set().union(*edges) - set().union(*open_edges)))
        wall_only_runners = tuple(r for r in range(n) if edges[r] and not open_edges[r])
        min_runner_measure = min(sum(chk.open_measure[r]) for r in range(n))
        total_measure = sum(sum(row) for row in chk.open_measure)
        adj, alt, support = section_tournament(chk)
        rows.append(
            {
                "speeds": speeds,
                "checkoff": chk,
                "matching": matching,
                "open_matching": open_matching,
                "section_cover": section_cover,
                "open_section_cover": open_section_cover,
                "all_runners": all_runners,
                "deficits": deficits,
                "open_deficits": open_deficits,
                "wall_debt_sections": wall_debt_sections,
                "wall_only_runners": wall_only_runners,
                "min_runner_measure": min_runner_measure,
                "total_measure": total_measure,
                "support": support,
                "fingerprint": fingerprint(adj),
                "edge_flips_private_first": edge_flips(adj, alt),
            }
        )
    return sets, rows, hall_counter, open_hall_counter


def print_example(row, label):
    chk = row["checkoff"]
    n = len(chk.speeds)
    print(f"  {label}: speeds={chk.speeds}")
    print(
        f"    perfect_matching={row['matching']} open_matching={row['open_matching']} "
        f"section_cover={row['section_cover']}/{n} "
        f"open_cover={row['open_section_cover']}/{n} "
        f"min_runner_open_measure={fmt(row['min_runner_measure'])}"
    )
    print(f"    runner->sections={tuple(tuple(sorted(x)) for x in row['support']['edges'])}")
    print(f"    section->runners={row['support']['runners_by_section']}")
    print(f"    section_open_measure={tuple(str(x) for x in row['support']['measure_by_section'])}")
    if row["deficits"]:
        d, subset, union = row["deficits"][0]
        print(f"    Hall packet: deficit={d}, runners={subset}, sections={union}")
    if row["open_deficits"]:
        d, subset, union = row["open_deficits"][0]
        print(f"    open Hall packet: deficit={d}, runners={subset}, sections={union}")
    print(
        f"    wall_debt_sections={row['wall_debt_sections']} "
        f"wall_only_runners={row['wall_only_runners']}"
    )
    print(
        f"    section tournament={row['fingerprint']}; "
        f"private-first edge flips={row['edge_flips_private_first']}"
    )


def main():
    print("=" * 78)
    print("LRC section check-off and switch atlas")
    print("=" * 78)
    print("Canonical gauge: include speed 0 and divide the loop into n fixed sections.")
    print("Dream certificate: perfect matching in runner<->section lonely-witness graph.")
    print("Exotic cases reduce to Hall packets plus section-support tournaments.")
    print()

    bounds = {4: 14, 5: 10, 6: 8}
    all_profiles = {}

    for n, bound in bounds.items():
        sets, rows, hall_counter, open_hall_counter = summarize_scan(n, bound)
        total = len(rows)
        matching = sum(1 for row in rows if row["matching"])
        open_matching = sum(1 for row in rows if row["open_matching"])
        all_runners = sum(1 for row in rows if row["all_runners"])
        full_cover = sum(1 for row in rows if row["section_cover"] == n)
        open_full_cover = sum(1 for row in rows if row["open_section_cover"] == n)
        no_hall = sum(1 for row in rows if not row["deficits"])
        open_no_hall = sum(1 for row in rows if not row["open_deficits"])
        all_profiles[n] = (matching, open_matching, full_cover, no_hall, total)

        print(f"(A) Exhaustive canonical speed sets: n={n}, max_speed={bound}")
        print(f"  primitive sets={total}")
        print(f"  every runner has some witnessed section={all_runners}/{total}")
        print(f"  full section cover={full_cover}/{total}")
        print(f"  open full section cover={open_full_cover}/{total}")
        print(f"  perfect section check-off matching={matching}/{total}")
        print(f"  open-cell-only matching={open_matching}/{total}")
        print(f"  no Hall deficit={no_hall}/{total}")
        print(f"  open no Hall deficit={open_no_hall}/{total}")
        print(f"  Hall packet histogram={(dict(sorted(hall_counter.items())))}")
        print(f"  open Hall packet histogram={(dict(sorted(open_hall_counter.items())))}")

        dream = next((row for row in rows if row["matching"]), None)
        exotic = next((row for row in rows if row["matching"] and row["open_deficits"]), None)
        tight = min(rows, key=lambda row: (row["min_runner_measure"], row["section_cover"]))
        if dream:
            print_example(dream, "dream check-off example")
        if exotic:
            print_example(exotic, "exotic wall-switch example")
        print_example(tight, "smallest min-runner-measure example")
        print()

    print("(B) Named rows")
    named = [
        ("LRC4 AP", (0, 1, 2, 3)),
        ("LRC4 uneven", (0, 1, 2, 4)),
        ("LRC5 AP", (0, 1, 2, 3, 4)),
        ("LRC6 cluster", (0, 1, 2, 3, 5, 8)),
        ("LRC14 AP", tuple(range(14))),
        ("LRC14 GW", (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 24)),
    ]
    for label, speeds in named:
        chk = checkoff_for(speeds, len(speeds))
        edges = edge_sets(chk)
        row = {
            "checkoff": chk,
            "matching": has_perfect_matching(edges, len(speeds)),
            "open_matching": has_perfect_matching(edge_sets(chk, False), len(speeds)),
            "section_cover": len(set().union(*edges)),
            "open_section_cover": len(set().union(*edge_sets(chk, False))),
            "all_runners": all(edges),
            "deficits": hall_deficits(edges, len(speeds)),
            "open_deficits": hall_deficits(edge_sets(chk, False), len(speeds)),
            "wall_debt_sections": tuple(
                sorted(set().union(*edges) - set().union(*edge_sets(chk, False)))
            ),
            "wall_only_runners": tuple(
                r for r in range(len(speeds)) if edges[r] and not edge_sets(chk, False)[r]
            ),
            "min_runner_measure": min(sum(chk.open_measure[r]) for r in range(len(speeds))),
            "support": section_support(chk),
        }
        adj, alt, _ = section_tournament(chk)
        row["fingerprint"] = fingerprint(adj)
        row["edge_flips_private_first"] = edge_flips(adj, alt)
        print_example(row, label)
    print()

    print("(C) Meta Tournament Analysis over n-level profiles")
    names = list(all_profiles)
    adj = tournament_from_scores(tuple(all_profiles[n] for n in names))
    print(f"  profiles={all_profiles}")
    print(f"  n-order={names}")
    print(f"  fingerprint={fingerprint(adj)}")
    print()

    print("SYNTHESIS")
    print("  1. A section proof can be exact after replacing 'runner i' by the")
    print("     runner-section witness graph.  The compactified dream case is Hall")
    print("     matching; strict-open matching exposes the wall debt.")
    print("  2. The exotic cases are not mysterious: their compressed certificate is")
    print("     an open Hall packet (runner subset, section union), wall-debt sections,")
    print("     and a switch tournament on the contested sections.")
    print("  3. The tournament analogy improves when vertices are sections: section")
    print("     dominance records which regions have private support, abundant open")
    print("     witness measure, or only wall witnesses.")
    print("  4. The proof route suggested by this data is local: show every Hall packet")
    print("     has a boundary-crossing switch that strictly enlarges its section union,")
    print("     or descends to an already-proved observer-source/endpoint-mouth case.")


if __name__ == "__main__":
    main()
