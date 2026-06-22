#!/usr/bin/env python3
"""S127: exact tournament-realizability atlas for the LRC14 summit.

The user prompt is methodological: if an LRC proof step has n comparable
objects, make the pairwise relation exact, quotient to tournament isomorphism
classes, and ask which classes are actually achievable.

This script keeps that discipline narrow and testable.  It builds two carriers:

1. The apex clock carrier: vertices are selected points on Z/14Z.  For points
   x,y, orient x -> y iff 0 < y-x < 7 on the 14-clock.  Diameter ties x-y=7
   are either excluded or broken by the listed Hamiltonian path.
2. The terminal runner-phase carrier: vertices are the thirteen speeds in an
   AP/Goddyn-Wong terminal row, placed at s mod 14.  The same half-clock cutoff
   orients non-tied pairs; collisions and diameters are tie-broken by increasing
   speed.  Exact M is still reported, because the tournament quotient discards
   off-apex metric escape data.

The result is not a proof of LRC14.  It is a realizability interface: a future
state-lift theorem must say why a bad LRC row lands in one of these achieved
tournament classes, and why the forbidden tournament/OCF class cannot occur.
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from functools import lru_cache
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, permutations, product
from pathlib import Path

import networkx as nx


REPO = Path(__file__).resolve().parents[1]
S124_PATH = REPO / "04-computation" / "lrc14_ap_gw_common_conditions_codex_s124.py"
spec = spec_from_file_location("s124_ap_gw", S124_PATH)
assert spec and spec.loader
s124 = module_from_spec(spec)
spec.loader.exec_module(s124)

D = 14
HALF = D // 2
AP = tuple(range(1, 14))
GW = tuple(list(range(1, 12)) + [13, 24])
PETAL_8_16 = tuple(sorted((set(AP) - {8}) | {16}))
PETAL_10_20 = tuple(sorted((set(AP) - {10}) | {20}))
RESIDUE_LIAR_26 = tuple(list(range(1, 12)) + [13, 26])
NEAR_MISS_36 = tuple(list(range(1, 12)) + [13, 36])


def orient_by_clock(a: int, b: int, modulus: int = D) -> int | None:
    """Return 1 for a->b, -1 for b->a, None for collision/diameter tie."""
    diff = (b - a) % modulus
    if diff == 0 or 2 * diff == modulus:
        return None
    return 1 if 2 * diff < modulus else -1


def mask_from_clock_points(points: tuple[int, ...], tie_policy: str) -> int | None:
    """Tournament mask on sorted clock points.

    Bits are indexed by i<j in vertex order.  Bit 1 means i->j; bit 0 means
    j->i.  The tie policy is:
      * exclude: any collision/diameter tie makes the row non-achievable;
      * path: orient every tie forward along the listed Hamiltonian path.
    """
    n = len(points)
    mask = 0
    bit = 0
    for i in range(n):
        for j in range(i + 1, n):
            ori = orient_by_clock(points[i], points[j])
            if ori is None:
                if tie_policy == "exclude":
                    return None
                ori = 1
            if ori == 1:
                mask |= 1 << bit
            bit += 1
    return mask


def mask_from_runner_row(row: tuple[int, ...]) -> int:
    """Tie-broken half-clock tournament for speeds ordered by increasing speed."""
    speeds = tuple(sorted(row))
    residues = tuple(s % D for s in speeds)
    n = len(speeds)
    mask = 0
    bit = 0
    for i in range(n):
        for j in range(i + 1, n):
            ori = orient_by_clock(residues[i], residues[j])
            if ori is None:
                ori = 1  # speed-order Hamiltonian path for collision/diameter ties
            if ori == 1:
                mask |= 1 << bit
            bit += 1
    return mask


def edge(mask: int, n: int, i: int, j: int) -> bool:
    """Return True iff i -> j in the tournament mask."""
    if i == j:
        raise ValueError("no loop")
    if i < j:
        bit = i * n - i * (i + 1) // 2 + (j - i - 1)
        return bool(mask & (1 << bit))
    return not edge(mask, n, j, i)


def encode_order(mask: int, order: tuple[int, ...]) -> int:
    """Encode the relabelled tournament in the supplied old-vertex order."""
    n = len(order)
    out = 0
    bit = 0
    for a in range(n):
        for b in range(a + 1, n):
            if edge(mask, n, order[a], order[b]):
                out |= 1 << bit
            bit += 1
    return out


def scores(mask: int, n: int) -> tuple[int, ...]:
    return tuple(sum(1 for j in range(n) if j != i and edge(mask, n, i, j)) for i in range(n))


@lru_cache(maxsize=None)
def canonical_mask(mask: int, n: int) -> int:
    """Exact isomorphism canonical form for n<=8, using score-block pruning."""
    sc = scores(mask, n)
    groups: list[list[int]] = []
    for value in sorted(set(sc), reverse=True):
        groups.append([i for i, s in enumerate(sc) if s == value])

    best: int | None = None
    for blocks in product(*(permutations(group) for group in groups)):
        order = tuple(v for block in blocks for v in block)
        enc = encode_order(mask, order)
        if best is None or enc < best:
            best = enc
    assert best is not None
    return best


def tournament_fingerprint(mask: int, n: int) -> dict[str, object]:
    sc = scores(mask, n)
    c3 = 0
    for tri in combinations(range(n), 3):
        local = [sum(1 for j in tri if i != j and edge(mask, n, i, j)) for i in tri]
        if sorted(local) == [1, 1, 1]:
            c3 += 1

    adj = {i: {j for j in range(n) if i != j and edge(mask, n, i, j)} for i in range(n)}
    radj = {i: set() for i in range(n)}
    for i, outs in adj.items():
        for j in outs:
            radj[j].add(i)

    seen: set[int] = set()
    order: list[int] = []

    def dfs(v: int) -> None:
        seen.add(v)
        for w in adj[v]:
            if w not in seen:
                dfs(w)
        order.append(v)

    for v in range(n):
        if v not in seen:
            dfs(v)

    seen.clear()
    scc: list[int] = []

    def rdfs(v: int, comp: list[int]) -> None:
        seen.add(v)
        comp.append(v)
        for w in radj[v]:
            if w not in seen:
                rdfs(w, comp)

    for v in reversed(order):
        if v not in seen:
            comp: list[int] = []
            rdfs(v, comp)
            scc.append(len(comp))

    hp = hamiltonian_path_count(mask, n)
    return {
        "score_hist": tuple(sorted(Counter(sc).items())),
        "score_sorted": tuple(sorted(sc, reverse=True)),
        "c3": c3,
        "scc": tuple(sorted(scc, reverse=True)),
        "hp": hp,
    }


def hamiltonian_path_count(mask: int, n: int) -> int:
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for state in range(1 << n):
        for last in range(n):
            val = dp[state][last]
            if not val:
                continue
            for nxt in range(n):
                if state & (1 << nxt):
                    continue
                if edge(mask, n, last, nxt):
                    dp[state | (1 << nxt)][nxt] += val
    return sum(dp[full])


def digraph_from_mask(mask: int, n: int) -> nx.DiGraph:
    graph = nx.DiGraph()
    graph.add_nodes_from(range(n))
    for i in range(n):
        for j in range(i + 1, n):
            if edge(mask, n, i, j):
                graph.add_edge(i, j)
            else:
                graph.add_edge(j, i)
    return graph


def classify_large_masks(masks: list[int], n: int) -> list[int]:
    """Class IDs by exact DiGraph isomorphism, for small lists with n>8."""
    reps: list[nx.DiGraph] = []
    ids: list[int] = []
    for mask in masks:
        graph = digraph_from_mask(mask, n)
        found = None
        for idx, rep in enumerate(reps):
            if nx.is_isomorphic(graph, rep):
                found = idx
                break
        if found is None:
            found = len(reps)
            reps.append(graph)
        ids.append(found)
    return ids


@dataclass(frozen=True)
class ClassSummary:
    class_id: int
    count: int
    representative: tuple[int, ...]
    fp: dict[str, object]


def clock_class_summaries(tie_policy: str) -> tuple[int, int, list[ClassSummary], Counter[int]]:
    classes: dict[int, list[tuple[int, ...]]] = defaultdict(list)
    diameter_hist: Counter[int] = Counter()
    labelled = 0
    for points in combinations(range(D), 7):
        diameter_count = sum(1 for x in points if (x + HALF) % D in points) // 2
        mask = mask_from_clock_points(points, tie_policy)
        if mask is None:
            continue
        labelled += 1
        diameter_hist[diameter_count] += 1
        can = canonical_mask(mask, 7)
        classes[can].append(points)

    summaries = []
    for idx, (can, reps) in enumerate(sorted(classes.items(), key=lambda kv: (-len(kv[1]), kv[0])), 1):
        summaries.append(
            ClassSummary(
                class_id=idx,
                count=len(reps),
                representative=reps[0],
                fp=tournament_fingerprint(can, 7),
            )
        )
    return labelled, len(classes), summaries, diameter_hist


def print_clock_section() -> None:
    print("[Apex-clock tournament realizability]")
    print("  vertices: seven selected points on Z/14Z")
    print("  observable: clockwise displacement")
    print("  switch/gauge: x->y iff 0 < y-x < 7; diameter/collision ties either excluded or")
    print("    broken by the Hamiltonian path induced by increasing residue")
    print("  quotient preserves: half-clock order class and diameter-tie count")
    print("  quotient destroys: actual speed divisibility, runner labels, off-apex M")
    print()

    for tie_policy in ("exclude", "path"):
        labelled, class_count, summaries, diameter_hist = clock_class_summaries(tie_policy)
        print(f"  tie_policy={tie_policy}: labelled rows={labelled}, isomorphism classes={class_count}")
        print(f"    diameter_pair_hist={dict(sorted(diameter_hist.items()))}")
        print("    largest achieved classes:")
        for item in summaries[:8]:
            print(
                f"      class {item.class_id:02d}: count={item.count:4d}, "
                f"rep={item.representative}, score_hist={item.fp['score_hist']}, "
                f"c3={item.fp['c3']}, scc={item.fp['scc']}, hp={item.fp['hp']}"
            )
        print()


def row_signature(row: tuple[int, ...]) -> str:
    M, pts = s124.M_exact(row)
    denoms = sorted({t.denominator for t in pts})
    return f"q={s124.q_threshold(row)}, M={M}, denoms={denoms}"


def print_terminal_runner_section() -> None:
    rows = [
        ("AP", AP),
        ("GW 12->24", GW),
        ("loose 8->16", PETAL_8_16),
        ("loose 10->20", PETAL_10_20),
        ("residue-liar 12->26", RESIDUE_LIAR_26),
        ("Farey near-miss 12->36", NEAR_MISS_36),
    ]
    masks = [mask_from_runner_row(row) for _, row in rows]
    class_ids = classify_large_masks(masks, 13)
    print("[Terminal runner-phase tournament classes]")
    print("  vertices: the thirteen speeds in increasing-speed order")
    print("  observable: residues s mod 14 at the denominator-14 apex")
    print("  switch/gauge: half-clock cutoff; collisions and diameters broken by speed order")
    print("  tie Hamiltonian path: increasing speed")
    print("  quotient preserves: apex winding/tie class")
    print("  quotient destroys: q-threshold divisibility and off-apex escape size")
    print()
    by_class: dict[int, list[str]] = defaultdict(list)
    for (name, row), mask, cid in zip(rows, masks, class_ids, strict=True):
        fp = tournament_fingerprint(mask, 13)
        by_class[cid].append(name)
        print(
            f"  {name:24s} class=T{cid} {row_signature(row)} "
            f"score_hist={fp['score_hist']} c3={fp['c3']} scc={fp['scc']} hp={fp['hp']}"
        )
    print("  class groups:")
    for cid, names in sorted(by_class.items()):
        print(f"    T{cid}: {', '.join(names)}")
    print()


def print_assumption_challenge() -> None:
    print("[Assumption challenge]")
    print("  considered vertex sets:")
    print("    runners/speeds, residues, gaps, fixed circle sections, section boundaries,")
    print("    wall-crossing events, cover arcs, Fourier modes, matroid/cycle circuits,")
    print("    and proof obligations/filters.")
    print("  chosen here:")
    print("    Z/14 apex points and denominator-14 runner phases, because these preserve")
    print("    the exact half-clock cutoff where AP/GW bind.")
    print("  challenged assumption:")
    print("    raw runner tournaments are not enough.  The residue-liar and Farey")
    print("    near-miss can share coarse winding data with tight rows while q(S) and")
    print("    off-apex M separate them.  A proof must state both the tournament")
    print("    realizability theorem and the arithmetic data the quotient intentionally")
    print("    forgets.")
    print()


def main() -> None:
    print("S127 LRC14 TOURNAMENT-REALIZABILITY SUMMIT ATLAS")
    print("=" * 72)
    print_assumption_challenge()
    print_clock_section()
    print_terminal_runner_section()
    print("[Proof readout]")
    print("  Tournament Analysis is useful only after three exact choices are fixed:")
    print("    (1) vertex set, (2) pairwise observable/cutoff, (3) tie Hamiltonian path.")
    print("  For the apex clock, the achieved tournament classes are a small finite")
    print("  atlas inside the 2x7 clock.  For terminal AP/GW rows, the apex-winding")
    print("  quotient separates some geometry but not the full LRC statement.")
    print("  The next summit lemma should be a state lift:")
    print("    bad LRC14 row -> achieved tournament/OCF class -> forbidden class absent,")
    print("  with q-threshold and off-apex witness data carried outside the quotient.")


if __name__ == "__main__":
    main()
