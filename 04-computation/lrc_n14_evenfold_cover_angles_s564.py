#!/usr/bin/env python3
"""
lrc_n14_evenfold_cover_angles_s564.py

codex-2026-06-02 S564

New n=14 proof angles from the recursive multi-sieve ledger.

Core reduction:
  For n=14, split speeds into even and odd.  The even-fold route says the set

      G = {t in [0,1): ||v t|| >= 1/14 for every even speed v}

  is the free region supplied by lower-n LRC.  The remaining n=14 problem is
  whether the odd danger arcs cover G.

This script computes the interval cover exactly for deterministic frontier
rows, then records:

  * exact measure of G and of the fully safe set G minus odd danger;
  * first closed rational witness t=a/q up to Q;
  * greedy odd-cover chains on G components and their first gaps;
  * a Tournament Analysis fingerprint on odd runners using pairwise exclusive
    pressure on G.

Proof-angle reading:
  A counterexample must make every positive G component covered by odd arcs and
  must also kill all boundary endpoints.  The AP/V* rows are wall covers with
  zero safe measure but closed rational witnesses.  The residual rows audited
  here have positive exact safe measure, often with low-denominator witnesses
  in the first fine window.  This suggests an n=14 route:

      coarse sieve witness
      or exact even-good cover gap
      or wall endpoint witness
      or owner-labelled residual frontier export.

Tournament Analysis declaration:
  Vertices: odd runners inside a fixed speed-set obligation, not all runners.
  Pairwise observable: exclusive odd-danger pressure on the even-good set G.
  Switch/gauge: x -> y when x has more pressure on G not already covered by y;
    ties follow increasing speed as the Hamiltonian path.
  Fingerprints: score histograms, directed 3-cycles, SCC sizes, edge flips
    against scalar pressure order, and Hamiltonian-path counts when feasible.

Assumption challenge:
  Alternative vertices considered: runners, G components, danger intervals,
  greedy handoff events, endpoint walls, residues, CRT channels, cover arcs,
  fine denominators, and proof obligations.  This pass uses odd runners because
  the even-fold quotient preserves exactly the predicate "do odd arcs cover G?"
  It destroys endpoint ownership within a runner's many arcs, so the handoff is
  to refine odd runners into owner-labelled danger intervals next.
"""

from __future__ import annotations

from collections import Counter, deque
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import combinations, permutations
from math import gcd


N = 14
TH = F(1, N)
QMAX = 220
Interval = tuple[F, F]
LabelledInterval = tuple[F, F, int]


@dataclass(frozen=True)
class Sample:
    label: str
    speeds: tuple[int, ...]
    note: str


@dataclass(frozen=True)
class Witness:
    a: int
    q: int
    margin: F


def fmt(x: F) -> str:
    if x.denominator == 1:
        return str(x.numerator)
    return f"{x.numerator}/{x.denominator}"


def fmt_float(x: F) -> str:
    return f"{float(x):.6f}"


def normalize(speeds: tuple[int, ...]) -> tuple[int, ...]:
    g = 0
    for s in speeds:
        g = gcd(g, abs(s))
    return tuple(sorted({abs(s) // g for s in speeds if s}))


def dist_frac(x: F) -> F:
    x %= 1
    return min(x, 1 - x)


def collar_at(speeds: tuple[int, ...], t: F) -> F:
    return min(dist_frac(v * t) for v in speeds)


def measure(intervals: list[Interval]) -> F:
    return sum((hi - lo for lo, hi in intervals), F(0))


def merge(intervals: list[Interval]) -> list[Interval]:
    intervals = sorted((lo, hi) for lo, hi in intervals if lo <= hi)
    out: list[Interval] = []
    for lo, hi in intervals:
        if not out or lo > out[-1][1]:
            out.append((lo, hi))
        elif hi > out[-1][1]:
            out[-1] = (out[-1][0], hi)
    return out


def intersect(a: list[Interval], b: list[Interval]) -> list[Interval]:
    out: list[Interval] = []
    i = j = 0
    while i < len(a) and j < len(b):
        lo = max(a[i][0], b[j][0])
        hi = min(a[i][1], b[j][1])
        if lo <= hi:
            out.append((lo, hi))
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return merge(out)


def subtract(a: list[Interval], b: list[Interval]) -> list[Interval]:
    b = merge(b)
    out: list[Interval] = []
    for lo, hi in a:
        cursor = lo
        for blo, bhi in b:
            if bhi <= cursor:
                continue
            if blo >= hi:
                break
            if blo > cursor:
                out.append((cursor, min(blo, hi)))
            cursor = max(cursor, bhi)
            if cursor >= hi:
                break
        if cursor < hi:
            out.append((cursor, hi))
    return merge(out)


def safe_intervals_for_speed(v: int) -> list[Interval]:
    return merge([(F(r, v) + TH / v, F(r + 1, v) - TH / v) for r in range(v)])


def danger_intervals_for_speed(v: int) -> list[Interval]:
    intervals: list[Interval] = []
    for r in range(v):
        intervals.append((F(r, v), F(r, v) + TH / v))
        intervals.append((F(r + 1, v) - TH / v, F(r + 1, v)))
    return merge(intervals)


def even_good_set(speeds: tuple[int, ...]) -> list[Interval]:
    g = [(F(0), F(1))]
    for v in speeds:
        if v % 2 == 0:
            g = intersect(g, safe_intervals_for_speed(v))
    return g


def odd_danger_on_g(speeds: tuple[int, ...], g_set: list[Interval]) -> list[Interval]:
    danger: list[Interval] = []
    for v in speeds:
        if v % 2 == 1:
            danger.extend(intersect(g_set, danger_intervals_for_speed(v)))
    return merge(danger)


def labelled_odd_danger_on_component(speeds: tuple[int, ...], comp: Interval) -> list[LabelledInterval]:
    base = [comp]
    out: list[LabelledInterval] = []
    for v in speeds:
        if v % 2 == 0:
            continue
        for lo, hi in intersect(base, danger_intervals_for_speed(v)):
            if lo < hi:
                out.append((lo, hi, v))
    return sorted(out)


def greedy_cover_component(speeds: tuple[int, ...], comp: Interval) -> tuple[bool, list[LabelledInterval], Interval | None]:
    lo, hi = comp
    intervals = labelled_odd_danger_on_component(speeds, comp)
    cursor = lo
    chain: list[LabelledInterval] = []
    while cursor < hi:
        candidates = [seg for seg in intervals if seg[0] <= cursor and seg[1] > cursor]
        if not candidates:
            next_lo = min((seg[0] for seg in intervals if seg[0] > cursor), default=hi)
            return False, chain, (cursor, min(next_lo, hi))
        best = max(candidates, key=lambda seg: (seg[1], -seg[2]))
        chain.append(best)
        cursor = best[1]
    return True, chain, None


def first_closed_witness(speeds: tuple[int, ...], qmax: int = QMAX) -> Witness | None:
    for q in range(2, qmax + 1):
        for a in range(1, q):
            margin = collar_at(speeds, F(a, q))
            if margin >= TH:
                return Witness(a, q, margin)
    return None


def odd_pressure_tournament(speeds: tuple[int, ...], g_set: list[Interval]) -> dict[str, object]:
    odds = [v for v in speeds if v % 2 == 1]
    n = len(odds)
    per = {
        v: intersect(g_set, danger_intervals_for_speed(v))
        for v in odds
    }
    pressure = {v: measure(per[v]) for v in odds}
    adj = [[False] * n for _ in range(n)]
    flips_vs_pressure = 0
    for i, x in enumerate(odds):
        for j, y in enumerate(odds):
            if i == j:
                continue
            x_excl = measure(subtract(per[x], per[y]))
            y_excl = measure(subtract(per[y], per[x]))
            if x_excl != y_excl:
                adj[i][j] = x_excl > y_excl
            else:
                adj[i][j] = i < j
            if i < j and pressure[x] != pressure[y]:
                pressure_order = pressure[x] > pressure[y]
                if adj[i][j] != pressure_order:
                    flips_vs_pressure += 1

    scores = [sum(row) for row in adj]
    c3 = 0
    for i, j, k in combinations(range(n), 3):
        if (adj[i][j] and adj[j][k] and adj[k][i]) or (adj[i][k] and adj[k][j] and adj[j][i]):
            c3 += 1

    def reach(start: int) -> set[int]:
        seen = {start}
        todo = deque([start])
        while todo:
            u = todo.popleft()
            for v, edge in enumerate(adj[u]):
                if edge and v not in seen:
                    seen.add(v)
                    todo.append(v)
        return seen

    remaining = set(range(n))
    sccs: list[int] = []
    while remaining:
        u = next(iter(remaining))
        ru = reach(u)
        comp = {v for v in remaining if v in ru and u in reach(v)}
        sccs.append(len(comp))
        remaining -= comp

    hp: int | str = 0
    if n <= 8:
        for perm in permutations(range(n)):
            if all(adj[perm[i]][perm[i + 1]] for i in range(n - 1)):
                hp += 1
    else:
        hp = "skipped"

    top_pressure = sorted(((pressure[v], v) for v in odds), reverse=True)[:3]
    return {
        "vertices": odds,
        "top_pressure": [(v, fmt(p)) for p, v in top_pressure],
        "score_hist": dict(sorted(Counter(scores).items())),
        "directed_3_cycles": c3,
        "sccs": sorted(sccs, reverse=True),
        "hamiltonian_paths": hp,
        "edge_flips_vs_pressure_order": flips_vs_pressure,
    }


def packet(n: int, scale: int, skip: int) -> tuple[int, ...]:
    return (1,) + tuple(scale * q for q in range(1, n) if q != skip)


def samples() -> list[Sample]:
    return [
        Sample("AP_wall", tuple(range(1, 14)), "wall row; zero safe measure but closed q=14 witness"),
        Sample("V_star_wall", tuple(list(range(1, 12)) + [13, 24]), "sporadic tight wall row"),
        Sample("sieve_minimal_lonely", tuple(range(2, 15)), "contains apex 14 and omits runner 1"),
        Sample("near_AP_apex", (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 14), "AP with 12 replaced by apex 14"),
        Sample("S562_packet_n14", packet(14, 7, 6), "HYP-2073 n=14 residual packet"),
        Sample("S562_packet_n14_lift", packet(14, 14, 6), "dyadic lift of HYP-2073 n=14 packet"),
        Sample("no_small_pinch_proxy", (1, 2, 9, 26, 110, 153, 166, 170, 178, 190, 192, 196, 201), "THM-396 proxy: lonely but no small-pinch witness"),
    ]


def analyze(sample: Sample) -> dict[str, object]:
    speeds = normalize(sample.speeds)
    g_set = even_good_set(speeds)
    odd_cover = odd_danger_on_g(speeds, g_set)
    safe = subtract(g_set, odd_cover)
    positive_g = [(lo, hi) for lo, hi in g_set if lo < hi]
    positive_safe = [(lo, hi) for lo, hi in safe if lo < hi]
    witness = first_closed_witness(speeds)
    chains = [greedy_cover_component(speeds, comp) for comp in positive_g]
    complete = sum(1 for ok, _, _ in chains if ok)
    gaps = [(gap, chain) for ok, chain, gap in chains if not ok and gap is not None]
    max_chain = max((len(chain) for _, chain, _ in chains), default=0)
    first_gap = gaps[0][0] if gaps else None
    return {
        "speeds": speeds,
        "g_measure": measure(positive_g),
        "safe_measure": measure(positive_safe),
        "g_components": len(positive_g),
        "safe_components": len(positive_safe),
        "first_witness": witness,
        "complete_components": complete,
        "gap_components": len(gaps),
        "first_gap": first_gap,
        "max_chain": max_chain,
        "tournament": odd_pressure_tournament(speeds, positive_g),
    }


def main() -> None:
    print("S564 n=14 even-fold cover angles")
    print("=" * 78)
    print("Exact reduction: even-good set G, odd danger cover, closed witnesses.")
    print(f"threshold=1/{N}, witness cutoff Q={QMAX}\n")

    rows = []
    for sample in samples():
        row = analyze(sample)
        rows.append((sample, row))
        w = row["first_witness"]
        w_s = "-" if w is None else f"{w.a}/{w.q} margin={fmt(w.margin)}"
        gap = row["first_gap"]
        gap_s = "-" if gap is None else f"{fmt(gap[0])}..{fmt(gap[1])} len={fmt(gap[1] - gap[0])}"
        print(sample.label)
        print(f"  note: {sample.note}")
        print(f"  speeds={row['speeds']}")
        print(
            f"  |G|={fmt(row['g_measure'])} ({fmt_float(row['g_measure'])})  "
            f"|G\\D_odd|={fmt(row['safe_measure'])} ({fmt_float(row['safe_measure'])})"
        )
        print(
            f"  components: G={row['g_components']} safe={row['safe_components']} "
            f"complete_odd_cover={row['complete_components']} gap_components={row['gap_components']} "
            f"max_greedy_chain={row['max_chain']}"
        )
        print(f"  first closed witness <=Q: {w_s}")
        print(f"  first greedy cover gap: {gap_s}")
        print(f"  odd-pressure tournament: {row['tournament']}")

    print("\nProof-angle synthesis")
    print("  Angle A: exact even-fold cover.  Prove that, outside wall rows, odd danger")
    print("    cannot cover every positive component of G.  S564 gives exact cover gaps")
    print("    for residual rows where coarse q<=14 sieves are blind.")
    print("  Angle B: wall endpoint fallback.  AP/V* have zero safe measure but retain")
    print("    closed rational witnesses; a counterexample must also kill all such")
    print("    endpoints, not merely cover G in measure.")
    print("  Angle C: recursive fine denominators.  The first witnesses for HYP-2073")
    print("    residual packets sit in the first fine window (23) or nearby fine tiers,")
    print("    matching HYP-2076's local-tier branch.")
    print("  Angle D: no-return cover chains.  A full odd cover of a G component is a")
    print("    left-to-right greedy handoff chain.  To use the user's hidden")
    print("    transitivity fact, refine odd runners into interval owner events and")
    print("    forbid return triangles in those handoffs.")


if __name__ == "__main__":
    main()
