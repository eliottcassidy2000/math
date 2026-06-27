#!/usr/bin/env python3
"""S155: taut bridge graph curvature for the LRC14 endpoint arrangement.

This is a local refinement of HYP-2970's endpoint-credit winding graph.  The
new object is not another scalar mass.  It is the exact boundary layer left by
the open danger arcs:

  * positive safe intervals are directed bridges between endpoint owners;
  * isolated safe points are zero-length taut vertices;
  * the owner-current records which mod-14/gcd labels are paid at the boundary.

The theorem-facing question is whether every non-AP/GW row creates a positive
bridge, violates a dual certificate such as HYP-2974, or carries a named
K33/state-lift transfer.  A strict counterexample would have neither positive
bridges nor taut vertices.

Tournament Analysis is over boundary/proof carriers rather than runners.  The
pairwise observable is how much of the LRC predicate survives before
scalarization: open witness, boundary-only equality, endpoint credit, Fourier
duality, C27/K33 state labels, and exact Farey scale.
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from math import gcd
from pathlib import Path
import argparse
import sys


REPO = Path(__file__).resolve().parents[1]
AP = tuple(range(1, 14))
THRESHOLD = Fraction(1, 14)


def load_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    assert spec and spec.loader
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


s138 = load_module(
    "s155_s138_frontier",
    REPO / "04-computation" / "lrc14_c27_two_swap_frontier_codex_s138.py",
)
s147 = load_module(
    "s155_s147_haar",
    REPO / "04-computation" / "lrc14_baire_haar_anyangle_codex_s147.py",
)
s148 = load_module(
    "s155_s148_gauntlet",
    REPO / "04-computation" / "lrc14_adversarial_counterexample_gauntlet_codex_s148.py",
)


@dataclass(frozen=True)
class Arc:
    speed: int
    center: int
    left: Fraction
    right: Fraction
    length: Fraction

    @property
    def owner(self) -> str:
        return owner_label(self.speed)


@dataclass(frozen=True)
class PositiveBridge:
    lo: Fraction
    hi: Fraction
    length: Fraction
    left_end: tuple[Arc, ...]
    right_start: tuple[Arc, ...]

    @property
    def transition(self) -> str:
        return f"{owners_text(self.left_end)} -> {owners_text(self.right_start)}"


@dataclass(frozen=True)
class TautVertex:
    t: Fraction
    left_depth: int
    right_depth: int
    point_depth: int
    left_end: tuple[Arc, ...]
    right_start: tuple[Arc, ...]
    zero_sum_pairs: tuple[tuple[int, int], ...]
    pair_sum_mod14: tuple[int, ...]

    @property
    def transition(self) -> str:
        return f"{owners_text(self.left_end)} -> {owners_text(self.right_start)}"


@dataclass(frozen=True)
class RowTautAudit:
    name: str
    family: str
    speeds: tuple[int, ...]
    qdiv: int
    M: Fraction
    denoms: tuple[int, ...]
    safe_mu: Fraction
    positive_bridges: tuple[PositiveBridge, ...]
    taut_vertices: tuple[TautVertex, ...]
    owner_current_l1: int
    owner_support: tuple[str, ...]
    curvature: str
    route: str


def fmt(x: Fraction) -> str:
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def mod1(x: Fraction) -> Fraction:
    return x - (x.numerator // x.denominator)


def primitive(row: tuple[int, ...]) -> bool:
    g = 0
    for v in row:
        g = gcd(g, v)
    return g == 1


def qdiv(row: tuple[int, ...]) -> int:
    for d in range(2, 200):
        if all(v % d != 0 for v in row):
            return d
    return 200


def owner_label(speed: int) -> str:
    return f"{speed % 14}:g{gcd(speed, 14)}"


def owners_text(arcs: tuple[Arc, ...]) -> str:
    if not arcs:
        return "none"
    counts = Counter(arc.owner for arc in arcs)
    return ",".join(
        f"{label}" if count == 1 else f"{label}x{count}"
        for label, count in sorted(counts.items())
    )


def arcs_for(row: tuple[int, ...]) -> tuple[Arc, ...]:
    arcs: list[Arc] = []
    for speed in row:
        length = Fraction(1, 7 * speed)
        radius = THRESHOLD / speed
        for center in range(speed):
            c = Fraction(center, speed)
            left = mod1(c - radius)
            right = mod1(c + radius)
            arcs.append(Arc(speed, center, left, right, length))
    return tuple(arcs)


def arc_contains_open(arc: Arc, t: Fraction) -> bool:
    delta = mod1(t - arc.left)
    return Fraction(0) < delta < arc.length


def endpoint_map(arcs: tuple[Arc, ...]) -> dict[Fraction, dict[str, list[Arc]]]:
    events: dict[Fraction, dict[str, list[Arc]]] = defaultdict(lambda: {"start": [], "end": []})
    for arc in arcs:
        events[arc.left]["start"].append(arc)
        events[arc.right]["end"].append(arc)
    for event in events.values():
        event["start"].sort(key=lambda a: (a.speed, a.center))
        event["end"].sort(key=lambda a: (a.speed, a.center))
    return events


def event_depths(arcs: tuple[Arc, ...], t: Fraction, event: dict[str, list[Arc]]) -> tuple[int, int, int]:
    point_depth = sum(1 for arc in arcs if arc_contains_open(arc, t))
    left_depth = point_depth + len(event["end"])
    right_depth = point_depth + len(event["start"])
    return left_depth, point_depth, right_depth


def safe_component_owners(
    components: list[tuple[Fraction, Fraction]], events: dict[Fraction, dict[str, list[Arc]]]
) -> tuple[PositiveBridge, ...]:
    bridges: list[PositiveBridge] = []
    for lo, hi in components:
        left_event = events.get(mod1(lo), {"start": [], "end": []})
        right_event = events.get(mod1(hi), {"start": [], "end": []})
        bridges.append(
            PositiveBridge(
                lo=lo,
                hi=hi,
                length=hi - lo,
                left_end=tuple(left_event["end"]),
                right_start=tuple(right_event["start"]),
            )
        )
    return tuple(bridges)


def taut_vertices(arcs: tuple[Arc, ...], events: dict[Fraction, dict[str, list[Arc]]]) -> tuple[TautVertex, ...]:
    out: list[TautVertex] = []
    for t in sorted(events):
        event = events[t]
        left_depth, point_depth, right_depth = event_depths(arcs, t, event)
        if point_depth != 0 or left_depth == 0 or right_depth == 0:
            continue
        pair_mods = []
        zero_pairs = []
        for left in event["end"]:
            for right in event["start"]:
                pair_mod = (left.speed + right.speed) % 14
                pair_mods.append(pair_mod)
                if pair_mod == 0:
                    zero_pairs.append((left.speed, right.speed))
        out.append(
            TautVertex(
                t=t,
                left_depth=left_depth,
                right_depth=right_depth,
                point_depth=point_depth,
                left_end=tuple(event["end"]),
                right_start=tuple(event["start"]),
                zero_sum_pairs=tuple(sorted(zero_pairs)),
                pair_sum_mod14=tuple(sorted(pair_mods)),
            )
        )
    return tuple(out)


def owner_current(
    bridges: tuple[PositiveBridge, ...], tauts: tuple[TautVertex, ...]
) -> tuple[int, tuple[str, ...]]:
    current: Counter[str] = Counter()
    support: set[str] = set()
    for bridge in bridges:
        for arc in bridge.left_end:
            current[arc.owner] += 1
            support.add(arc.owner)
        for arc in bridge.right_start:
            current[arc.owner] -= 1
            support.add(arc.owner)
    for taut in tauts:
        for arc in taut.left_end:
            current[arc.owner] += 1
            support.add(arc.owner)
        for arc in taut.right_start:
            current[arc.owner] -= 1
            support.add(arc.owner)
    return sum(abs(v) for v in current.values()), tuple(sorted(support))


def exact_m(row: tuple[int, ...]) -> tuple[Fraction, tuple[int, ...]]:
    _, M, denoms = s138.exact_m_worker(row)
    return M, denoms


def classify_curvature(safe_mu: Fraction, bridges: tuple[PositiveBridge, ...], tauts: tuple[TautVertex, ...]) -> str:
    if safe_mu > 0:
        return "positive-open-bridge"
    if tauts:
        zero_count = sum(len(t.zero_sum_pairs) for t in tauts)
        if zero_count:
            return "zero-curvature-taut"
        return "nonzero-taut"
    return "fully-covered-no-taut"


def route_for(name: str, family: str, row: tuple[int, ...], curvature: str, M: Fraction) -> str:
    if M < THRESHOLD:
        return "COUNTEREXAMPLE"
    gw = tuple(list(range(1, 12)) + [13, 24])
    if curvature == "zero-curvature-taut" and row in {AP, gw}:
        return "AP/GW-TAUT-EQUALITY"
    if "K33" in name or "12->36" in name:
        return "K33-STATE-LIFT"
    if "petal" in family or "P10" in name or "P13" in name:
        return "C27-PETAL"
    if qdiv(row) > 14:
        return "COVERING-POSITIVE"
    if curvature == "positive-open-bridge":
        return "OPEN-BRIDGE"
    return "UNKNOWN-TAUT"


def audit_row(name: str, family: str, row: tuple[int, ...], compute_m: bool = True) -> RowTautAudit:
    row = tuple(sorted(row))
    data = s147.exact_row_measure(row)
    safe_mu = data["safe_measure"]
    arcs = arcs_for(row)
    events = endpoint_map(arcs)
    bridges = safe_component_owners(data["safe_components"], events)
    tauts = taut_vertices(arcs, events) if safe_mu == 0 else tuple()
    current_l1, support = owner_current(bridges, tauts)
    M, denoms = exact_m(row) if compute_m else (Fraction(-1), tuple())
    curvature = classify_curvature(safe_mu, bridges, tauts)
    route = route_for(name, family, row, curvature, M)
    return RowTautAudit(
        name=name,
        family=family,
        speeds=row,
        qdiv=qdiv(row),
        M=M,
        denoms=denoms,
        safe_mu=safe_mu,
        positive_bridges=bridges,
        taut_vertices=tauts,
        owner_current_l1=current_l1,
        owner_support=support,
        curvature=curvature,
        route=route,
    )


def named_rows() -> list[tuple[str, str, tuple[int, ...]]]:
    ap_set = set(AP)
    return [
        ("AP", "boundary", AP),
        ("GW 12->24", "boundary", tuple(list(range(1, 12)) + [13, 24])),
        ("near/K33 12->36", "K33", tuple(list(range(1, 12)) + [13, 36])),
        ("petal 10->20", "petal", tuple(sorted((ap_set - {10}) | {20}))),
        ("petal 13->26", "petal", tuple(sorted((ap_set - {13}) | {26}))),
        ("P10+GW", "petal+boundary", tuple(sorted((ap_set - {10, 12}) | {20, 24}))),
        ("P10+K33", "petal+K33", tuple(sorted((ap_set - {10, 12}) | {20, 36}))),
        ("residue liar 12->26", "q-witness", tuple(list(range(1, 12)) + [13, 26])),
        ("magnitude liar 12->96", "q14-loose", tuple(list(range(1, 12)) + [13, 96])),
        ("covering 12->84", "covering", tuple(list(range(1, 12)) + [13, 84])),
        ("covering drop13 add182", "covering", tuple(list(range(1, 13)) + [182])),
    ]


def unique_rows(rows: list[tuple[str, str, tuple[int, ...]]]) -> list[tuple[str, str, tuple[int, ...]]]:
    seen: set[tuple[int, ...]] = set()
    out = []
    for name, family, row in rows:
        row = tuple(sorted(set(row)))
        if len(row) != 13 or not primitive(row) or row in seen:
            continue
        seen.add(row)
        out.append((name, family, row))
    return out


def bank_rows(one_limit: int, two_limit: int) -> list[tuple[str, str, tuple[int, ...]]]:
    rows: list[tuple[str, str, tuple[int, ...]]] = []
    ap = set(AP)
    for remove in AP:
        base = ap - {remove}
        for add in range(14, one_limit + 1):
            rows.append((f"single {remove}->{add}", "single-swap", tuple(sorted(base | {add}))))
    add_values = range(14, two_limit + 1)
    for rem in combinations(AP, 2):
        base = ap - set(rem)
        for add in combinations(add_values, 2):
            rows.append((f"two {rem}->{add}", "two-swap", tuple(sorted(base | set(add)))))
    return unique_rows(rows)


def scan_bank(one_limit: int, two_limit: int, max_zero_details: int) -> dict[str, object]:
    rows = bank_rows(one_limit, two_limit)
    curvature_counts: Counter[str] = Counter()
    zero_rows: list[RowTautAudit] = []
    min_positive: tuple[Fraction, str, tuple[int, ...]] | None = None
    for name, family, row in rows:
        data = s147.exact_row_measure(row)
        safe_mu = data["safe_measure"]
        if safe_mu > 0:
            curvature_counts["positive-open-bridge"] += 1
            if min_positive is None or safe_mu < min_positive[0]:
                min_positive = (safe_mu, name, row)
            continue
        audit = audit_row(name, family, row, compute_m=True)
        curvature_counts[audit.curvature] += 1
        zero_rows.append(audit)
    return {
        "rows": rows,
        "curvature_counts": curvature_counts,
        "zero_rows": zero_rows[:max_zero_details],
        "zero_total": len(zero_rows),
        "min_positive": min_positive,
    }


PROOF_CARRIERS = {
    "positive_open_bridge": (8, 8, 5, 3, 3, 7),
    "taut_vertex_current": (7, 7, 8, 5, 5, 7),
    "endpoint_credit_winding": (7, 6, 8, 7, 5, 8),
    "Toeplitz_PSD_dual": (6, 5, 4, 8, 8, 8),
    "multiplicity_moment_dual": (5, 5, 4, 6, 7, 7),
    "C27_K33_state_labels": (5, 5, 7, 8, 6, 8),
    "exact_M_Farey_scale": (8, 8, 3, 4, 4, 7),
    "raw_runner_vertices": (2, 2, 2, 1, 1, 1),
}

TIE_PATH = (
    "exact_M_Farey_scale",
    "positive_open_bridge",
    "taut_vertex_current",
    "endpoint_credit_winding",
    "Toeplitz_PSD_dual",
    "multiplicity_moment_dual",
    "C27_K33_state_labels",
    "raw_runner_vertices",
)


def carrier_winner(a: str, b: str) -> str:
    va = PROOF_CARRIERS[a]
    vb = PROOF_CARRIERS[b]
    score = sum(x > y for x, y in zip(va, vb)) - sum(x < y for x, y in zip(va, vb))
    if score > 0:
        return a
    if score < 0:
        return b
    return a if TIE_PATH.index(a) < TIE_PATH.index(b) else b


def tournament_fingerprint() -> dict[str, object]:
    vertices = tuple(PROOF_CARRIERS)
    wins = Counter()
    edges: dict[str, set[str]] = {v: set() for v in vertices}
    for a, b in combinations(vertices, 2):
        winner = carrier_winner(a, b)
        loser = b if winner == a else a
        wins[winner] += 1
        edges[winner].add(loser)
    hist = Counter(wins[v] for v in vertices)
    c3 = 0
    for a, b, c in combinations(vertices, 3):
        if b in edges[a] and c in edges[b] and a in edges[c]:
            c3 += 1
        if c in edges[a] and b in edges[c] and a in edges[b]:
            c3 += 1
    sccs = strongly_connected_components(vertices, edges)
    return {
        "score_hist": dict(sorted(hist.items())),
        "directed_3cycles": c3,
        "scc_sizes": sorted((len(s) for s in sccs), reverse=True),
        "order": sorted(vertices, key=lambda v: (-wins[v], TIE_PATH.index(v))),
    }


def strongly_connected_components(vertices: tuple[str, ...], edges: dict[str, set[str]]) -> list[set[str]]:
    def reach(start: str, graph: dict[str, set[str]]) -> set[str]:
        seen = {start}
        q = deque([start])
        while q:
            v = q.popleft()
            for w in graph[v]:
                if w not in seen:
                    seen.add(w)
                    q.append(w)
        return seen

    rev = {v: set() for v in vertices}
    for v, outs in edges.items():
        for w in outs:
            rev[w].add(v)
    unseen = set(vertices)
    sccs = []
    while unseen:
        v = next(iter(unseen))
        comp = reach(v, edges) & reach(v, rev)
        sccs.append(comp)
        unseen -= comp
    return sccs


def print_named_audits(audits: list[RowTautAudit]) -> None:
    print("[1] Named-row taut bridge audit")
    print(
        f"  {'row':24s} {'M':>8s} {'qdiv':>4s} {'safe_mu':>10s} "
        f"{'bridges':>7s} {'taut':>5s} {'curvature':25s} route"
    )
    for audit in audits:
        print(
            f"  {audit.name:24s} {fmt(audit.M):>8s} {audit.qdiv:4d} "
            f"{fmt(audit.safe_mu):>10s} {len(audit.positive_bridges):7d} "
            f"{len(audit.taut_vertices):5d} {audit.curvature:25s} {audit.route}"
        )
    print()
    print("  Taut equality atoms:")
    for audit in audits:
        if not audit.taut_vertices:
            continue
        zero_pairs = sum(len(t.zero_sum_pairs) for t in audit.taut_vertices)
        print(
            f"    {audit.name}: taut_vertices={len(audit.taut_vertices)} "
            f"zero_sum_pairs={zero_pairs} owner_current_l1={audit.owner_current_l1}"
        )
        for taut in audit.taut_vertices[:6]:
            print(
                f"      t={fmt(taut.t):>5s} left={owners_text(taut.left_end):20s} "
                f"right={owners_text(taut.right_start):20s} "
                f"pair_mod14={list(taut.pair_sum_mod14)[:8]}"
            )
    print()
    print("  Smallest positive bridges in named rows:")
    for audit in audits:
        if not audit.positive_bridges:
            continue
        best = min(audit.positive_bridges, key=lambda b: b.length)
        print(
            f"    {audit.name}: min_bridge={fmt(best.length)} "
            f"{fmt(best.lo)}->{fmt(best.hi)} {best.transition}"
        )
    print()


def print_bank_scan(scan: dict[str, object]) -> None:
    rows = scan["rows"]
    counts = scan["curvature_counts"]
    print("[2] AP-neighborhood exact open-mass bank")
    print(f"  rows_scanned={len(rows)}")
    print(f"  curvature_counts={dict(sorted(counts.items()))}")
    print(f"  zero_open_total={scan['zero_total']}")
    min_positive = scan["min_positive"]
    if min_positive:
        mu, name, row = min_positive
        print(f"  smallest_positive_safe_mu={fmt(mu)} at {name} row={row}")
    zero_rows = scan["zero_rows"]
    if zero_rows:
        print("  zero-open rows detailed:")
        for audit in zero_rows:
            print(
                f"    {audit.name}: M={fmt(audit.M)} qdiv={audit.qdiv} "
                f"taut={len(audit.taut_vertices)} route={audit.route} row={audit.speeds}"
            )
    print()


def print_tournament_analysis() -> None:
    fp = tournament_fingerprint()
    print("[3] Tournament Analysis")
    print("  vertices considered:")
    print("    runners, raw endpoints, endpoint-owner labels, positive safe intervals,")
    print("    isolated taut points, boundary-current states, missed-depth sectors,")
    print("    HYP-2970 endpoint-credit cycles, HYP-2974 Toeplitz dual sections,")
    print("    C27/K33 proof obligations.")
    print("  chosen vertices:")
    print("    boundary/proof carriers, especially positive_open_bridge and taut_vertex_current.")
    print("  pairwise observable:")
    print("    retention of open-witness status, boundary equality, endpoint ownership,")
    print("    dual-certifiability, state-lift address, and scalar-decoy resistance.")
    print("  switch/gauge:")
    print("    componentwise majority of the six retention scores; ties follow the")
    print("    Hamiltonian path exact_M -> open bridge -> taut current -> endpoint credit.")
    print(
        f"  fingerprint: score_hist={fp['score_hist']} "
        f"directed_3cycles={fp['directed_3cycles']} scc_sizes={fp['scc_sizes']}"
    )
    print("  Hamiltonian order:")
    for idx, name in enumerate(fp["order"], start=1):
        print(f"    {idx}. {name}")
    print()


def print_readout() -> None:
    print("[4] Readout")
    print("  HYP-2975 is not a standalone proof.  It contributes a local boundary")
    print("  curvature audit inside HYP-2970's endpoint-credit graph.")
    print("  Positive safe intervals are visible as bridges; AP/GW are zero-open but")
    print("  leave isolated taut vertices with mod-14 zero-sum endpoint transfers.")
    print("  The falsifier to search for is a primitive qdiv>=14 row with no positive")
    print("  bridge, no AP/GW taut current, and no K33/state-lift or Toeplitz-negative")
    print("  exit.  The bounded bank scanned here found none.")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--one-limit", type=int, default=160)
    parser.add_argument("--two-limit", type=int, default=36)
    parser.add_argument("--max-zero-details", type=int, default=12)
    args = parser.parse_args()

    print("S155 LRC14 TAUT BRIDGE GRAPH CURVATURE")
    print("=" * 78)
    print("[0] Scope")
    print("  This is a HYP-2975 local refinement of HYP-2970 endpoint-credit winding.")
    print("  It keeps exact open intervals and isolated endpoint witnesses before")
    print("  scalarizing to Haar mass or Fourier moments.")
    print()

    audits = [audit_row(name, family, row) for name, family, row in named_rows()]
    print_named_audits(audits)
    scan = scan_bank(args.one_limit, args.two_limit, args.max_zero_details)
    print_bank_scan(scan)
    print_tournament_analysis()
    print_readout()


if __name__ == "__main__":
    main()
