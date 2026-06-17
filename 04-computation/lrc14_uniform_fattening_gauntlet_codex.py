#!/usr/bin/env python3
"""
LRC(14) uniform-fattening gauntlet.

Target: OPEN-Q-108.  For a 12-speed core C, let

    G_C = {t in [0,1): ||v t|| > 1/14 for every v in C}.

The remaining singular-series proof of LRC(14) wants a uniform lower bound
meas(G_C) >= c > 0 over every primitive 12-subset C.  This script stress-tests
that statement around the known near-tight AP cores and records a tournament
fingerprint for the load-bearing structure.

Tournament Analysis contract:
  * Vertices: speeds in the 12-core C.
  * Pairwise observable: with a,b removed, compare how much of G_{C minus {a,b}} is
    covered by D_a versus D_b, where D_v={t: ||v t|| <= 1/14}.
  * Switch/gauge: orient a -> b when a covers more; ties use the Hamiltonian
    path sorted by speed.  We compare edge flips under scale by 2 and under
    threshold switch 1/14 -> 1/13.
  * Fingerprints: score histogram, directed 3-cycles, SCCs, edge flips, and
    Hamiltonian-path count.

Assumption challenge:
  Speeds are only one possible vertex set.  Other natural choices are safe-set
  components, band endpoints, q-witness grids, deleted AP positions, residues
  mod 14, pair-sum wall crossings, and proof obligations.  The speed tournament
  preserves the predicate "which runner most reduces meas(G_C) in a pairwise
  context", but it destroys endpoint alignment, component adjacency, and the
  q-grid witness logic.  If this speed tournament stays transitive, the next
  quotient should move to safe components or endpoint events.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from functools import lru_cache, reduce
from fractions import Fraction as F
from itertools import combinations
from math import gcd
import random


N = 14
THRESH = F(1, N)


def gcd_all(values):
    return reduce(gcd, values, 0)


def primitive(values):
    g = gcd_all(values)
    return tuple(sorted(v // g for v in values))


@lru_cache(maxsize=None)
def danger_arcs(v, q=N):
    """Intervals where ||v t|| <= 1/q, split inside [0,1)."""
    h = F(1, q * v)
    arcs = []
    for k in range(v):
        c = F(k, v)
        lo = c - h
        hi = c + h
        if lo < 0:
            arcs.append((lo + 1, F(1)))
            arcs.append((F(0), hi))
        elif hi > 1:
            arcs.append((lo, F(1)))
            arcs.append((F(0), hi - 1))
        else:
            arcs.append((lo, hi))
    return tuple(arcs)


def merge_intervals(arcs):
    arcs = sorted((a, b) for a, b in arcs if b > a)
    if not arcs:
        return []
    merged = []
    lo, hi = arcs[0]
    for a, b in arcs[1:]:
        if a <= hi:
            if b > hi:
                hi = b
        else:
            merged.append((lo, hi))
            lo, hi = a, b
    merged.append((lo, hi))
    return merged


def interval_measure(arcs):
    return sum(b - a for a, b in arcs)


def complement_intervals(merged):
    if not merged:
        return [(F(0), F(1))]
    comps = []
    cursor = F(0)
    for lo, hi in merged:
        if lo > cursor:
            comps.append((cursor, lo))
        if hi > cursor:
            cursor = hi
    if cursor < 1:
        comps.append((cursor, F(1)))
    return comps


@lru_cache(maxsize=None)
def safe_components(core, q=N):
    arcs = []
    for v in core:
        arcs.extend(danger_arcs(v, q))
    danger = merge_intervals(arcs)
    safe = complement_intervals(danger)
    return interval_measure(safe), tuple(safe)


def intersect_measure(a_intervals, b_intervals):
    i = j = 0
    total = F(0)
    a_intervals = list(a_intervals)
    b_intervals = list(b_intervals)
    while i < len(a_intervals) and j < len(b_intervals):
        a0, a1 = a_intervals[i]
        b0, b1 = b_intervals[j]
        lo = max(a0, b0)
        hi = min(a1, b1)
        if hi > lo:
            total += hi - lo
        if a1 < b1:
            i += 1
        else:
            j += 1
    return total


def cover_measure(safe, v, q=N):
    return intersect_measure(safe, merge_intervals(danger_arcs(v, q)))


def fmt(fr):
    return f"{fr} = {float(fr):.9f}"


def ap_drop_one_table():
    rows = []
    base = tuple(range(1, N))
    for missing in base:
        core = tuple(v for v in base if v != missing)
        measure, comps = safe_components(core, N)
        max_comp = max((b - a for a, b in comps), default=F(0))
        rows.append((measure, missing, core, len(comps), max_comp))
    return sorted(rows)


def two_drop_one_replacement(limit=180):
    """Remove two AP speeds and add one outside speed; keep the best positive core."""
    base = tuple(range(1, N))
    best = (F(99), None)
    tested = 0
    zero = 0
    for missing in combinations(base, 2):
        stem = [v for v in base if v not in missing]
        for w in range(N, limit + 1):
            if w in stem:
                continue
            core = primitive(stem + [w])
            if len(core) != 12:
                continue
            measure, comps = safe_components(core, N)
            tested += 1
            if measure == 0:
                zero += 1
            elif measure < best[0]:
                max_comp = max((b - a for a, b in comps), default=F(0))
                best = (measure, (missing, w, core, len(comps), max_comp))
    return tested, zero, best


def random_core_search(trials=3000, bound=90, seed=20260617):
    rng = random.Random(seed)
    best = (F(99), None)
    zero = 0
    tested = 0
    for _ in range(trials):
        core = primitive(rng.sample(range(1, bound + 1), 12))
        if len(core) != 12:
            continue
        measure, comps = safe_components(core, N)
        tested += 1
        if measure == 0:
            zero += 1
        elif measure < best[0]:
            max_comp = max((b - a for a, b in comps), default=F(0))
            best = (measure, (core, len(comps), max_comp))
    return tested, zero, best


def single_swap_greedy(seed_core, limit=180):
    current = tuple(seed_core)
    history = []
    for _ in range(2):
        start_measure, _ = safe_components(current, N)
        step_best = (start_measure, current, None)
        for old in current:
            stem = [v for v in current if v != old]
            for w in range(1, limit + 1):
                if w in stem:
                    continue
                candidate = primitive(stem + [w])
                if len(candidate) != 12:
                    continue
                measure, _ = safe_components(candidate, N)
                if measure > 0 and measure < step_best[0]:
                    step_best = (measure, candidate, (old, w))
        history.append(step_best)
        if step_best[1] == current:
            break
        current = step_best[1]
    return history


def load_tournament(core, q=N):
    core = tuple(sorted(core))
    n = len(core)
    adj = [[False] * n for _ in range(n)]
    weights = {}
    for i, a in enumerate(core):
        for j, b in enumerate(core):
            if i >= j:
                continue
            base = tuple(v for k, v in enumerate(core) if k not in (i, j))
            _, safe = safe_components(base, q)
            ca = cover_measure(safe, a, q)
            cb = cover_measure(safe, b, q)
            weights[(i, j)] = (ca, cb)
            if ca > cb or (ca == cb and a < b):
                adj[i][j] = True
            else:
                adj[j][i] = True
    return adj, weights


def tournament_scores(adj):
    return [sum(row) for row in adj]


def directed_triangles(adj):
    n = len(adj)
    cycles = 0
    transitive = 0
    for i, j, k in combinations(range(n), 3):
        verts = (i, j, k)
        outdeg = []
        for a in verts:
            outdeg.append(sum(1 for b in verts if a != b and adj[a][b]))
        if sorted(outdeg) == [1, 1, 1]:
            cycles += 1
        else:
            transitive += 1
    return cycles, transitive


def sccs(adj):
    n = len(adj)
    graph = [[j for j in range(n) if adj[i][j]] for i in range(n)]
    rev = [[i for i in range(n) if adj[i][j]] for j in range(n)]
    seen = [False] * n
    order = []

    def dfs(v):
        seen[v] = True
        for u in graph[v]:
            if not seen[u]:
                dfs(u)
        order.append(v)

    for v in range(n):
        if not seen[v]:
            dfs(v)
    comps = []
    seen = [False] * n

    def rdfs(v, comp):
        seen[v] = True
        comp.append(v)
        for u in rev[v]:
            if not seen[u]:
                rdfs(u, comp)

    for v in reversed(order):
        if not seen[v]:
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
            row = adj[last]
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                if row[nxt]:
                    dp[mask | (1 << nxt)][nxt] += val
    return sum(dp[-1])


def edge_flips(adj_a, adj_b):
    n = len(adj_a)
    flips = 0
    for i in range(n):
        for j in range(i + 1, n):
            if adj_a[i][j] != adj_b[i][j]:
                flips += 1
    return flips


def tournament_report(core):
    core = tuple(sorted(core))
    adj14, _ = load_tournament(core, 14)
    adj_scale, _ = load_tournament(tuple(2 * v for v in core), 14)
    adj13, _ = load_tournament(core, 13)
    scores = tournament_scores(adj14)
    cyc, trans = directed_triangles(adj14)
    comps = sccs(adj14)
    return {
        "core": core,
        "scores": scores,
        "score_hist": dict(sorted(Counter(scores).items())),
        "cycles": cyc,
        "transitive_triples": trans,
        "scc_sizes": sorted([len(c) for c in comps], reverse=True),
        "hp_count": hamiltonian_path_count(adj14),
        "scale2_flips": edge_flips(adj14, adj_scale),
        "threshold13_flips": edge_flips(adj14, adj13),
    }


def print_tournament_report(rep):
    print(f"  core={rep['core']}")
    print(f"    scores={rep['scores']}")
    print(f"    score_hist={rep['score_hist']}")
    print(
        "    triples: directed_cycles="
        f"{rep['cycles']}, transitive={rep['transitive_triples']}"
    )
    print(f"    SCC sizes={rep['scc_sizes']}")
    print(f"    Hamiltonian paths={rep['hp_count']}")
    print(f"    edge flips under scale-by-2 gauge={rep['scale2_flips']}")
    print(f"    edge flips under threshold switch 1/14 -> 1/13={rep['threshold13_flips']}")


def main():
    print("=" * 78)
    print("LRC(14) uniform-fattening gauntlet for OPEN-Q-108")
    print("=" * 78)
    print("Threshold: 1/14.  Core size: 12 speeds.")
    print()

    print("(A) AP drop-one 12-cores: exact meas(G_C)")
    drop_rows = ap_drop_one_table()
    for measure, missing, core, comp_count, max_comp in drop_rows:
        flag = "  <-- smallest in AP drop-one" if missing == 6 else ""
        print(
            f"  missing {missing:2d}: meas(G_C)={fmt(measure)}; "
            f"safe_components={comp_count}; max_safe_component={max_comp}{flag}"
        )

    print()
    print("(B) Two AP deletions plus one replacement, w <= 180")
    tested, zero, best = two_drop_one_replacement(limit=180)
    measure, detail = best
    print(f"  tested={tested}; zero-measure cores={zero}")
    print(f"  best positive meas(G_C)={fmt(measure)}")
    print(
        "  best detail: missing_pair={}, replacement={}, core={}, "
        "safe_components={}, max_safe_component={}".format(*detail)
    )
    print("  result: no two-drop/one-replacement core beat the AP drop-6 core.")

    print()
    print("(C) Random primitive 12-cores, values <= 90")
    tested, zero, best = random_core_search(trials=3000, bound=90)
    measure, detail = best
    print(f"  tested={tested}; zero-measure cores={zero}")
    print(f"  best positive meas(G_C)={fmt(measure)}")
    print(
        f"  best random detail: core={detail[0]}, safe_components={detail[1]}, "
        f"max_safe_component={detail[2]}"
    )

    print()
    print("(D) Single-swap greedy stress from three near-tight seeds")
    seeds = [
        tuple(v for v in range(1, N) if v != 6),
        (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13),
        detail[0],
    ]
    for seed in seeds:
        start_measure, _ = safe_components(tuple(seed), N)
        print(f"  seed={tuple(seed)}; start meas(G_C)={fmt(start_measure)}")
        for measure, core, swap in single_swap_greedy(tuple(seed), limit=180):
            print(f"    best step meas={fmt(measure)}; swap={swap}; core={core}")

    print()
    print("(E) Tournament Analysis fingerprints")
    drop6_core = tuple(v for v in range(1, N) if v != 6)
    best_replacement_core = detail[0]
    print("  Speed-load tournament on AP drop-6 core:")
    print_tournament_report(tournament_report(drop6_core))
    print("  Speed-load tournament on best random core:")
    print_tournament_report(tournament_report(best_replacement_core))

    print()
    print("Conclusion")
    print("  The AP drop-6 12-core has meas(G_C)=7/858.")
    print("  The searched coordinated replacements did not beat 7/858.")
    print("  This supports a sharper OPEN-Q-108 subtarget: prove the drop-6")
    print("  core is extremal, or prove every non-AP coordinated core has")
    print("  meas(G_C) >= 7/858.  The speed tournament is only a quotient;")
    print("  component/end-point tournaments should be tested next.")


if __name__ == "__main__":
    main()
