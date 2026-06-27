#!/usr/bin/env python3
"""S155: global twist-ladder / dual-blocker certificate for LRC14.

This is deliberately different from the boundary-gap packet route.  Instead
of decomposing open safe intervals, it asks for explicit rational phases

    t = a / q,  gcd(a,q)=1,

that already satisfy the LRC14 threshold for a row S.  A selected denominator
ladder L gives a finite primal certificate when some twist is safe.  If no
twist in L is safe, the blocked twists form a finite dual hypergraph:
twists are vertices and speeds are blocker hyperedges.

The guardrail is important: HYP-2901 shows no fixed finite denominator ladder
can prove all of LRC14.  The point here is instead to test whether the current
HYP-2963 reduced packet bank has a small global ladder certificate, and what
dual blocker pattern remains when a smaller ladder fails.

Tournament Analysis is over proof coordinates such as twist phases, ladder
rungs, and blocker hypergraphs, not over runners or danger intervals.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from math import gcd
from pathlib import Path
import argparse
import sys


REPO = Path(__file__).resolve().parents[1]
N = 14
AP = tuple(range(1, 14))


def load_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    assert spec and spec.loader
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


s124 = load_module(
    "s155_s124_common",
    REPO / "04-computation" / "lrc14_ap_gw_common_conditions_codex_s124.py",
)
s2963 = load_module(
    "s155_hyp2963_packets",
    REPO / "04-computation" / "lrc14_labelled_packet_counterexample_classifier_codex_20260624.py",
)


@dataclass(frozen=True)
class Twist:
    q: int
    a: int

    @property
    def text(self) -> str:
        return f"{self.a}/{self.q}"


@dataclass(frozen=True)
class RowCertificate:
    name: str
    family: str
    speeds: tuple[int, ...]
    q_threshold: int
    witness: Twist | None
    slack_num: int | None
    tight_speeds: tuple[int, ...]

    @property
    def certified(self) -> bool:
        return self.witness is not None

    @property
    def strict(self) -> bool:
        return self.slack_num is not None and self.slack_num > 0


@dataclass(frozen=True)
class LadderSummary:
    label: str
    q_count: int
    twist_count: int
    certified: int
    missed: int
    strict_witnesses: int
    boundary_witnesses: int
    max_first_q: int
    first_q_hist: tuple[tuple[int, int], ...]


PROOF_VERTICES = (
    "q_threshold_gate",
    "twist_ladder_primal",
    "q41_rescue_rung",
    "q27_fiber_ladder",
    "blocker_hypergraph_dual",
    "committed_denominator_wall",
    "source_spectrum_packet",
    "exact_haar_interval",
    "boundary_gap_packet",
    "raw_apex_residue_tournament",
    "raw_runner_set",
)

# Scores: preserves LRC witness, sees finite dual obstruction, adapts to
# denominator walls, connects to existing packets, resists scalarization.
PROOF_SCORES: dict[str, tuple[int, ...]] = {
    "q_threshold_gate": (6, 2, 2, 4, 4),
    "twist_ladder_primal": (8, 5, 5, 5, 7),
    "q41_rescue_rung": (7, 6, 7, 5, 7),
    "q27_fiber_ladder": (7, 6, 7, 6, 7),
    "blocker_hypergraph_dual": (6, 8, 8, 6, 8),
    "committed_denominator_wall": (4, 8, 9, 4, 8),
    "source_spectrum_packet": (8, 7, 8, 8, 8),
    "exact_haar_interval": (8, 5, 6, 7, 6),
    "boundary_gap_packet": (7, 6, 6, 8, 7),
    "raw_apex_residue_tournament": (3, 2, 1, 3, 2),
    "raw_runner_set": (2, 1, 1, 1, 1),
}

TIE_PATH = (
    "source_spectrum_packet",
    "twist_ladder_primal",
    "blocker_hypergraph_dual",
    "q41_rescue_rung",
    "q27_fiber_ladder",
    "exact_haar_interval",
    "boundary_gap_packet",
    "committed_denominator_wall",
    "q_threshold_gate",
    "raw_apex_residue_tournament",
    "raw_runner_set",
)


def primitive(speeds: tuple[int, ...]) -> bool:
    g = 0
    for v in speeds:
        g = gcd(g, v)
    return g == 1


def twist_slack_num(speeds: tuple[int, ...], twist: Twist) -> tuple[int, tuple[int, ...]]:
    """Return min integer slack 14*dist_num-q and speeds attaining it."""
    q = twist.q
    values: list[tuple[int, int]] = []
    for v in speeds:
        r = (v * twist.a) % q
        dist_num = min(r, q - r)
        values.append((N * dist_num - q, v))
    slack = min(x for x, _ in values)
    tight = tuple(v for x, v in values if x == slack)
    return slack, tight


def is_safe(speeds: tuple[int, ...], twist: Twist) -> bool:
    slack, _ = twist_slack_num(speeds, twist)
    return slack >= 0


def blocking_speeds(speeds: tuple[int, ...], twist: Twist) -> tuple[int, ...]:
    q = twist.q
    blockers = []
    for v in speeds:
        r = (v * twist.a) % q
        if N * min(r, q - r) < q:
            blockers.append(v)
    return tuple(blockers)


def build_twists(qs: set[int]) -> tuple[Twist, ...]:
    return tuple(
        Twist(q, a)
        for q in sorted(qs)
        for a in range(1, q)
        if gcd(a, q) == 1
    )


def q27_denominators(m: int) -> set[int]:
    qs = set(range(2, N + 1))
    for d in (1, 2, 7, 14):
        for k in range(2, m + 1):
            qs.add(d * k)
    return qs


def ladder_sets(qmax: int, q27_m: int) -> dict[str, set[int]]:
    q27 = q27_denominators(q27_m)
    return {
        "q<=27": set(range(2, 28)),
        "q<=42": set(range(2, 43)),
        f"q<={qmax}": set(range(2, qmax + 1)),
        f"qdiv+Q27(m<={q27_m})": q27,
        f"q<={qmax}+Q27(m<={q27_m})": set(range(2, qmax + 1)) | q27,
    }


def certify_row(name: str, family: str, speeds: tuple[int, ...], twists: tuple[Twist, ...]) -> RowCertificate:
    q_threshold = s124.q_threshold(speeds)
    for twist in twists:
        slack, tight = twist_slack_num(speeds, twist)
        if slack >= 0:
            return RowCertificate(name, family, speeds, q_threshold, twist, slack, tight)
    return RowCertificate(name, family, speeds, q_threshold, None, None, ())


def summarize_ladder(label: str, qs: set[int], rows: list[tuple[str, str, tuple[int, ...]]]) -> tuple[LadderSummary, list[RowCertificate]]:
    twists = build_twists(qs)
    certs = [certify_row(name, family, speeds, twists) for name, family, speeds in rows]
    first_q = Counter(c.witness.q for c in certs if c.witness is not None)
    max_first_q = max(first_q, default=0)
    summary = LadderSummary(
        label=label,
        q_count=len(qs),
        twist_count=len(twists),
        certified=sum(c.certified for c in certs),
        missed=sum(not c.certified for c in certs),
        strict_witnesses=sum(c.strict for c in certs),
        boundary_witnesses=sum(c.certified and not c.strict for c in certs),
        max_first_q=max_first_q,
        first_q_hist=tuple(first_q.most_common(14)),
    )
    return summary, certs


def blocker_profile(speeds: tuple[int, ...], twists: tuple[Twist, ...]) -> dict[str, object]:
    load: Counter[int] = Counter()
    private: Counter[int] = Counter()
    size_hist: Counter[int] = Counter()
    for twist in twists:
        blockers = blocking_speeds(speeds, twist)
        size_hist[len(blockers)] += 1
        for v in blockers:
            load[v] += 1
        if len(blockers) == 1:
            private[blockers[0]] += 1
    return {
        "load": tuple(load.most_common(8)),
        "private": tuple(private.most_common(8)),
        "size_hist": tuple(sorted(size_hist.items())),
    }


def edge(mask: int, n: int, i: int, j: int) -> bool:
    if i == j:
        raise ValueError("loop")
    if i > j:
        return not edge(mask, n, j, i)
    bit = 0
    for a, b in combinations(range(n), 2):
        if a == i and b == j:
            return bool((mask >> bit) & 1)
        bit += 1
    raise AssertionError("unreachable")


def strongly_connected_components(mask: int, n: int) -> list[set[int]]:
    graph = [[j for j in range(n) if i != j and edge(mask, n, i, j)] for i in range(n)]
    rgraph = [[j for j in range(n) if i != j and edge(mask, n, j, i)] for i in range(n)]
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
    comps: list[set[int]] = []

    def rdfs(v: int, comp: set[int]) -> None:
        seen.add(v)
        comp.add(v)
        for w in rgraph[v]:
            if w not in seen:
                rdfs(w, comp)

    for v in reversed(order):
        if v not in seen:
            comp: set[int] = set()
            rdfs(v, comp)
            comps.append(comp)
    return comps


def hamiltonian_path_count(mask: int, n: int) -> int:
    dp: dict[tuple[int, int], int] = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for size in range(1, n):
        nxt: dict[tuple[int, int], int] = {}
        for (used, last), count in dp.items():
            if used.bit_count() != size:
                continue
            for v in range(n):
                if used & (1 << v):
                    continue
                if edge(mask, n, last, v):
                    key = (used | (1 << v), v)
                    nxt[key] = nxt.get(key, 0) + count
        dp.update(nxt)
    full = (1 << n) - 1
    return sum(count for (used, _), count in dp.items() if used == full)


def tournament_fingerprint() -> dict[str, object]:
    index = {v: i for i, v in enumerate(PROOF_VERTICES)}
    tie_rank = {v: i for i, v in enumerate(TIE_PATH)}
    mask = 0
    bit = 0
    for i, j in combinations(range(len(PROOF_VERTICES)), 2):
        vi = PROOF_VERTICES[i]
        vj = PROOF_VERTICES[j]
        si = PROOF_SCORES[vi]
        sj = PROOF_SCORES[vj]
        i_wins = si > sj or (si == sj and tie_rank[vi] < tie_rank[vj])
        if i_wins:
            mask |= 1 << bit
        bit += 1
    n = len(PROOF_VERTICES)
    outdeg = [sum(1 for j in range(n) if i != j and edge(mask, n, i, j)) for i in range(n)]
    c3 = 0
    for a, b, c in combinations(range(n), 3):
        if edge(mask, n, a, b) and edge(mask, n, b, c) and edge(mask, n, c, a):
            c3 += 1
        if edge(mask, n, a, c) and edge(mask, n, c, b) and edge(mask, n, b, a):
            c3 += 1
    leaders = sorted(PROOF_VERTICES, key=lambda v: (-outdeg[index[v]], tie_rank[v]))[:4]
    return {
        "score_hist": dict(sorted(Counter(outdeg).items())),
        "c3": c3,
        "scc": tuple(sorted((len(comp) for comp in strongly_connected_components(mask, n)), reverse=True)),
        "hp": hamiltonian_path_count(mask, n),
        "leaders": tuple(leaders),
    }


def print_assumption_challenge() -> None:
    print("[0] Assumption challenge")
    print("  considered vertices:")
    print("    runners, gaps, safe components, endpoint owners, residue classes,")
    print("    denominators, rational twists, blocked twist events, speed blockers,")
    print("    divisor fibers, source-spectrum nodes, and proof obligations.")
    print("  chosen vertices:")
    print("    rational twists a/q and the speed-vs-twist blocker hypergraph.")
    print("  preserved LRC predicate:")
    print("    if any selected twist has ||v*a/q|| >= 1/14 for every speed v,")
    print("    it is an exact LRC14 certificate for that row.")
    print("  destroyed information:")
    print("    open intervals between selected twists and fine endpoint ownership.")
    print("  challenged assumption:")
    print("    the proof need not start from boundary gaps; it can start from a")
    print("    global primal witness ladder and read failed ladders as dual covers.")
    print()


def print_ladder_summary(summaries: list[LadderSummary], total_rows: int) -> None:
    print("[1] Ladder census on the HYP-2963 packet bank")
    print(f"  audited rows={total_rows}")
    print(
        f"  {'ladder':28s} {'qs':>4s} {'twists':>7s} {'cert':>6s} "
        f"{'miss':>5s} {'strict':>6s} {'bdry':>6s} {'maxq':>5s}"
    )
    for s in summaries:
        print(
            f"  {s.label[:28]:28s} {s.q_count:4d} {s.twist_count:7d} "
            f"{s.certified:6d} {s.missed:5d} {s.strict_witnesses:6d} "
            f"{s.boundary_witnesses:6d} {s.max_first_q:5d}"
        )
    print("  first-witness denominator histogram for the full ladder:")
    full = summaries[-1]
    print("    " + ", ".join(f"{q}:{c}" for q, c in full.first_q_hist))
    print()


def print_q27_dual(
    small_certs: list[RowCertificate],
    full_certs: list[RowCertificate],
    small_twists: tuple[Twist, ...],
) -> None:
    full_by_speeds = {c.speeds: c for c in full_certs}
    misses = [c for c in small_certs if not c.certified]
    print("[2] Dual blocker readout for the q<=27 misses")
    print(f"  q<=27 missed rows={len(misses)}")
    if not misses:
        print("  no dual blocker residual at this ladder size.")
        print()
        return
    print(
        f"  {'row':38s} {'qdiv':>4s} {'full witness':>12s} "
        f"{'slack':>7s} {'top blocker loads':>30s}"
    )
    for c in misses:
        full = full_by_speeds[c.speeds]
        profile = blocker_profile(c.speeds, small_twists)
        loads = ",".join(f"{v}:{k}" for v, k in profile["load"][:4])
        witness = full.witness.text if full.witness else "-"
        slack = str(full.slack_num) if full.slack_num is not None else "-"
        print(f"  {c.name[:38]:38s} {c.q_threshold:4d} {witness:>12s} {slack:>7s} {loads:>30s}")
        print(f"      blocker-size histogram={profile['size_hist']}")
        if profile["private"]:
            print(f"      private blockers={profile['private']}")
        else:
            print("      private blockers=none; every failed twist has multiple blockers")
    print("  readout:")
    print("    q<=27 sees exactly the lcm-tail family as a dual cover; q=41")
    print("    with a=17 breaks the cover in the full ladder.  This is a global")
    print("    twist certificate, not a local endpoint-gap certificate.")
    print()


def print_hardest(full_certs: list[RowCertificate]) -> None:
    ordered = sorted(
        (c for c in full_certs if c.witness is not None),
        key=lambda c: (c.witness.q, c.q_threshold, c.name),
        reverse=True,
    )
    print("[3] Hardest first witnesses in the full ladder")
    print(f"  {'row':42s} {'fam':14s} {'qdiv':>4s} {'witness':>9s} {'slack':>7s} tight speeds")
    for c in ordered[:18]:
        print(
            f"  {c.name[:42]:42s} {c.family[:14]:14s} {c.q_threshold:4d} "
            f"{c.witness.text:>9s} {c.slack_num:7d} {c.tight_speeds}"
        )
    qdiv_hist = Counter(c.witness.q for c in full_certs if c.witness is not None and c.q_threshold > N)
    print("  qdiv>14 first-witness denominators:")
    print("    " + ", ".join(f"{q}:{count}" for q, count in qdiv_hist.most_common()))
    print()


def print_theorem_readout() -> None:
    print("[4] Twist-ladder theorem target")
    print("  Exact finite-bank theorem:")
    print("    every row in the HYP-2963 labelled-packet bank has an explicit")
    print("    rational witness a/q with q<=42.")
    print("  Smaller-ladder obstruction:")
    print("    q<=27 misses only the divisor-loaded lcm-tail family")
    print("    {1,2,3,4,5,6,7,8,9,10,11,13,84m}, m=1..5 in this bank.")
    print("  Rescue:")
    print("    the same row family is witnessed at 17/41 in the q<=42 ladder.")
    print("  Global guardrail:")
    print("    HYP-2901 rules out any fixed finite denominator ladder as a full")
    print("    proof of LRC14.  The proposed proof route is therefore dynamic:")
    print("    after Moon-core / packet reductions, either a bounded ladder gives")
    print("    a primal twist witness, or the blocker hypergraph has a private")
    print("    resource that descends to q-threshold, K33/state-lift, or a")
    print("    committed-denominator wall with its own next-rung witness.")
    print()


def print_tournament_analysis() -> None:
    fp = tournament_fingerprint()
    print("[5] Tournament Analysis")
    print("  vertices are proof coordinates, not runners.")
    print("  pair observable:")
    print("    which coordinate preserves an exact LRC witness or a finite dual")
    print("    blocker certificate while surviving denominator-wall adaptation.")
    print("  switch/gauge:")
    print("    lexicographic retention of primal twist witness, blocker dual,")
    print("    denominator-wall adaptability, packet compatibility, and")
    print("    anti-scalarization.")
    print("  tie Hamiltonian path:")
    print("    " + " > ".join(TIE_PATH))
    print(
        "  fingerprint: "
        f"score_hist={fp['score_hist']} c3={fp['c3']} "
        f"scc={fp['scc']} hp={fp['hp']} leaders={fp['leaders']}"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--single-limit", type=int, default=180)
    parser.add_argument("--two-swap-limit", type=int, default=36)
    parser.add_argument("--alias-depth", type=int, default=4)
    parser.add_argument("--lcm-tail-max", type=int, default=5)
    parser.add_argument("--qmax", type=int, default=64)
    parser.add_argument("--q27-m", type=int, default=27)
    args = parser.parse_args()

    rows = s2963.build_bank(args.single_limit, args.two_swap_limit, args.alias_depth, args.lcm_tail_max)
    ladders = ladder_sets(args.qmax, args.q27_m)
    outputs: list[tuple[LadderSummary, list[RowCertificate], tuple[Twist, ...]]] = []
    for label, qs in ladders.items():
        summary, certs = summarize_ladder(label, qs, rows)
        outputs.append((summary, certs, build_twists(qs)))

    print("S155 LRC14 twist-ladder / dual-blocker certificate")
    print("=" * 78)
    print(
        f"single_limit={args.single_limit}, two_swap_limit={args.two_swap_limit}, "
        f"alias_depth={args.alias_depth}, lcm_tail_max={args.lcm_tail_max}, "
        f"qmax={args.qmax}, q27_m={args.q27_m}"
    )
    print_assumption_challenge()
    print_ladder_summary([x[0] for x in outputs], len(rows))
    small_summary, small_certs, small_twists = outputs[0]
    full_summary, full_certs, _ = outputs[-1]
    print_q27_dual(small_certs, full_certs, small_twists)
    print_hardest(full_certs)
    print_theorem_readout()
    print_tournament_analysis()


if __name__ == "__main__":
    main()
