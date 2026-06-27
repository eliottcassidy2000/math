#!/usr/bin/env python3
"""HYP-3133 scout: A000568 as a middle edge-witness quotient.

Executable synthesis, not an LRC14 proof.

Prompt cue: 12 sits between 10 and 16, and 56 sits between 20 and 80.
Here 10,20 are four-sector edge words from HYP-3124, 16,80 are the
sector-plus-paired-tail/tip-child signatures, and 12,56 are A000568
unlabeled tournament counts one vertex higher.

Tournament Analysis declaration:
  vertices: proof carriers and quotient shadows, not runners or scalar counts;
  pairwise observable: majority comparison over retained proof-payload axes;
  switch/gauge: orient A->B when A wins more axes, ties follow the declared
                Hamiltonian path;
  tie path: SPEC certificate, edge packet, A000568 middle shadow, paired child
            deck, sector word, even-graph parity, Asano legality, constants,
            raw counts.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from functools import lru_cache
from itertools import combinations
import sys
from typing import Dict, List, Sequence, Tuple


sys.path.append("04-computation")
import lrc14_tournament_edge_witness_recursion_codex_s268 as edge_scout  # noqa: E402


# OEIS A000568: number of unlabeled tournaments on n vertices, offset n=0.
# Values used here are the small exact values needed for the shifted sandwich.
A000568 = {
    0: 1,
    1: 1,
    2: 1,
    3: 2,
    4: 4,
    5: 12,
    6: 56,
    7: 456,
    8: 6880,
}

AXES = (
    "lrc_predicate_retention",
    "edge_sandwich_explains",
    "tail_tip_recursion",
    "equidecomposition_power",
    "equinumerosity_guard",
    "spec_resonance_bridge",
    "finite_constant_chase",
    "formalization_readiness",
    "scalar_guardrail",
)


@dataclass(frozen=True)
class Carrier:
    name: str
    axes: Dict[str, int]
    preserves: str
    destroys: str
    next_hook: str


def carrier(
    name: str,
    scores: Sequence[int],
    preserves: str,
    destroys: str,
    next_hook: str,
) -> Carrier:
    if len(scores) != len(AXES):
        raise ValueError(name)
    return Carrier(name, dict(zip(AXES, scores)), preserves, destroys, next_hook)


@lru_cache(maxsize=None)
def one_sided_child(mask: int, n: int, deleted: int) -> Tuple[object, ...]:
    child, child_n = edge_scout.delete_vertex(mask, n, deleted)
    return edge_scout.tournament_sig(child, child_n)


def edge_sandwich_census(n: int) -> Dict[str, object]:
    """Count sector words and paired child signatures on n-vertex edge witnesses."""
    total_tournaments = 1 << (n * (n - 1) // 2)
    sectors: Counter[Tuple[int, int, int, int]] = Counter()
    child_pairs: Counter[Tuple[object, ...]] = Counter()
    sector_to_child_pairs: Dict[Tuple[int, int, int, int], set] = {}

    for mask in range(total_tournaments):
        for tail, tip in edge_scout.directed_edges(mask, n):
            sector = edge_scout.sector_word(mask, n, tail, tip)
            tail_child = one_sided_child(mask, n, tail)
            tip_child = one_sided_child(mask, n, tip)
            child_pair = (sector, tail_child, tip_child)
            sectors[sector] += 1
            child_pairs[child_pair] += 1
            sector_to_child_pairs.setdefault(sector, set()).add(child_pair)

    top_splits = sorted(
        (len(children), sectors[sector], sector)
        for sector, children in sector_to_child_pairs.items()
    )[-5:]
    top_splits.reverse()
    shifted = A000568[n + 1]
    lower = len(sectors)
    upper = len(child_pairs)
    total_gap = upper - lower
    lower_gap = shifted - lower
    upper_gap = upper - shifted
    return {
        "edge_n": n,
        "a000568_n": n + 1,
        "labelled_tournaments": total_tournaments,
        "directed_edge_instances": total_tournaments * n * (n - 1) // 2,
        "sector_words": lower,
        "a000568_shifted": shifted,
        "paired_child_signatures": upper,
        "lower_gap": lower_gap,
        "upper_gap": upper_gap,
        "gap_position": (lower_gap / total_gap) if total_gap else None,
        "sector_groups_split_by_child_pair": sum(
            1 for children in sector_to_child_pairs.values() if len(children) > 1
        ),
        "max_child_pairs_inside_one_sector": max(len(v) for v in sector_to_child_pairs.values()),
        "top_sector_splits": top_splits,
    }


CARRIERS: List[Carrier] = [
    carrier(
        "SPEC_resonance_lattice_certificate",
        (10, 7, 6, 6, 5, 10, 10, 8, 9),
        "HYP-3129 exact-low plus Parseval-tail floor for the multi-far SPEC term",
        "tail/tip edge ownership if used as a scalar Fourier bound alone",
        "attach sector/extension shadows to the finite low-frequency constant chase",
    ),
    carrier(
        "edge_floor_two_ended_packet",
        (9, 8, 10, 9, 6, 8, 7, 8, 9),
        "HYP-3125 R-safe -> Q-safe edge packet with both endpoint deletion children",
        "free middle quotient unless the A000568 shadow is stored",
        "add `a000568_extension_shadow` to the edge_floor_packet schema",
    ),
    carrier(
        "A000568_unlabeled_extension_shadow",
        (7, 10, 6, 7, 10, 6, 6, 8, 9),
        "one-extra-vertex isomorphism quotient sitting between sector words and paired children",
        "which endpoint deletion card owns the obstruction",
        "use as a middle quotient before escalating to the full paired child deck",
    ),
    carrier(
        "paired_tail_tip_child_deck",
        (8, 8, 10, 10, 6, 5, 6, 7, 9),
        "equidecomposition data from both endpoint-deletion cards",
        "unrooted global extension economy; can over-specify the quotient",
        "keep as the repair layer when A000568 shadow is not enough",
    ),
    carrier(
        "four_sector_tetrahedral_word",
        (5, 7, 3, 4, 8, 4, 5, 8, 7),
        "weak-composition sector sizes around a directed edge",
        "outside-vertex internal tournament and both deletion cards",
        "use as the equinumerosity floor, never as the full witness",
    ),
    carrier(
        "even_graph_cycle_space_parity",
        (6, 7, 4, 7, 9, 5, 5, 6, 8),
        "OEIS-linked equinumerosity between tournaments and even graphs",
        "orientation/chirality unless an edge packet reattaches it",
        "test whether HYP-3129 resonance rows have a cycle-space parity sidecar",
    ),
    carrier(
        "Asano_tip_tail_legality_guard",
        (7, 5, 8, 7, 5, 5, 6, 6, 8),
        "HYP-3128/HYP-3127 tip legality and the warning that the crowded tail is not zero-free",
        "the elementary SPEC floor if used as the main engine",
        "use only as a legality sidecar after the SPEC certificate is named",
    ),
    carrier(
        "closed_form_constant_chase",
        (8, 5, 4, 4, 5, 8, 10, 7, 8),
        "the remaining HYP-3129 task: uniform lower SPEC_low and upper tail constants",
        "which quotient caused a bad finite row",
        "stratify constant-chase rows by sector, A000568 shadow, and child pair",
    ),
    carrier(
        "raw_sequence_numerology",
        (1, 2, 0, 1, 2, 0, 1, 1, 0),
        "the visible numbers 10,12,16 and 20,56,80",
        "the LRC predicate, quotient map, recursion, and proof obligations",
        "keep only as the observation that triggered the middle quotient",
    ),
]

TIE_PATH = [c.name for c in CARRIERS]


def compare(a: Carrier, b: Carrier) -> str:
    a_votes = b_votes = 0
    for axis in AXES:
        if a.axes[axis] > b.axes[axis]:
            a_votes += 1
        elif b.axes[axis] > a.axes[axis]:
            b_votes += 1
    if a_votes > b_votes:
        return a.name
    if b_votes > a_votes:
        return b.name
    return a.name if TIE_PATH.index(a.name) < TIE_PATH.index(b.name) else b.name


def tournament_edges(vertices: Sequence[Carrier]) -> Dict[str, set]:
    adj: Dict[str, set] = {v.name: set() for v in vertices}
    for a, b in combinations(vertices, 2):
        winner = compare(a, b)
        loser = b.name if winner == a.name else a.name
        adj[winner].add(loser)
    return adj


def directed_3cycles(adj: Dict[str, set]) -> int:
    count = 0
    for a, b, c in combinations(sorted(adj), 3):
        if (b in adj[a] and c in adj[b] and a in adj[c]) or (
            a in adj[b] and b in adj[c] and c in adj[a]
        ):
            count += 1
    return count


def scc_sizes(adj: Dict[str, set]) -> List[int]:
    index = 0
    stack: List[str] = []
    on_stack = set()
    indices: Dict[str, int] = {}
    low: Dict[str, int] = {}
    sizes: List[int] = []

    def visit(v: str) -> None:
        nonlocal index
        indices[v] = low[v] = index
        index += 1
        stack.append(v)
        on_stack.add(v)
        for w in adj[v]:
            if w not in indices:
                visit(w)
                low[v] = min(low[v], low[w])
            elif w in on_stack:
                low[v] = min(low[v], indices[w])
        if low[v] == indices[v]:
            size = 0
            while True:
                w = stack.pop()
                on_stack.remove(w)
                size += 1
                if w == v:
                    break
            sizes.append(size)

    for v in adj:
        if v not in indices:
            visit(v)
    return sorted(sizes, reverse=True)


def hamiltonian_path_count(adj: Dict[str, set]) -> int:
    names = list(adj)
    idx = {name: i for i, name in enumerate(names)}
    n = len(names)
    dp: Dict[Tuple[int, int], int] = {}
    for i in range(n):
        dp[(1 << i, i)] = 1
    for mask in range(1 << n):
        for last in range(n):
            paths = dp.get((mask, last), 0)
            if not paths:
                continue
            for nxt_name in adj[names[last]]:
                nxt = idx[nxt_name]
                if mask & (1 << nxt):
                    continue
                dp[(mask | (1 << nxt), nxt)] = dp.get((mask | (1 << nxt), nxt), 0) + paths
    full = (1 << n) - 1
    return sum(dp.get((full, last), 0) for last in range(n))


def one_hamiltonian_path(adj: Dict[str, set]) -> List[str]:
    names = list(adj)

    def search(path: List[str], unused: set) -> List[str] | None:
        if not unused:
            return path
        last = path[-1]
        for nxt in sorted(unused, key=lambda x: -len(adj[x])):
            if nxt in adj[last]:
                found = search(path + [nxt], unused - {nxt})
                if found:
                    return found
        return None

    for start in sorted(names, key=lambda x: -len(adj[x])):
        found = search([start], set(names) - {start})
        if found:
            return found
    return []


def sandwich_report() -> List[str]:
    lines = [
        "2. SHIFTED A000568 EDGE-WITNESS SANDWICH",
        "A000568(n)=number of unlabeled tournaments on n vertices.",
        "Compare edge-local signatures on m vertices with A000568(m+1).",
        (
            "m  sector_words  A000568(m+1)  paired_child_signatures  "
            "lower_gap  upper_gap  gap_position  status"
        ),
    ]
    summaries = [edge_sandwich_census(n) for n in range(2, 7)]
    for s in summaries:
        lower = s["sector_words"]
        middle = s["a000568_shifted"]
        upper = s["paired_child_signatures"]
        if lower < middle < upper:
            status = "proper_sandwich"
        elif lower == middle == upper:
            status = "trivial_equality"
        else:
            status = "outside_small_boundary"
        position = s["gap_position"]
        position_text = "NA" if position is None else f"{position:.3f}"
        lines.append(
            f"{s['edge_n']}  {lower:12d}  {middle:13d}  {upper:23d}  "
            f"{s['lower_gap']:9d}  {s['upper_gap']:9d}  {position_text:>12}  {status}"
        )
    lines.append("")
    lines.append(
        "Readout: the user's 10<12<16 and 20<56<80 pattern is the first two "
        "proper cases of a shifted sandwich that continues at m=6 as "
        "35<456<632.  The m=3 row is equality, and m=2 is a small boundary "
        "where the free one-vertex quotient is larger than the available "
        "two-ended edge deck."
    )
    lines.append("")
    lines.append("m=6 largest sector refinements:")
    for count, instances, sector in summaries[-1]["top_sector_splits"]:
        lines.append(f"- sector={sector}: child_pairs={count}, instances={instances}")
    return lines


def tournament_report() -> List[str]:
    adj = tournament_edges(CARRIERS)
    path = one_hamiltonian_path(adj)
    return [
        "3. TOURNAMENT ANALYSIS OVER QUOTIENT CARRIERS",
        "vertices=proof carriers and quotient shadows, not runners or scalar sequences",
        "pairwise_observable=majority over axes " + ",".join(AXES),
        "switch=orient A->B when A wins more axes; ties use the declared Hamiltonian path",
        "tie_hamiltonian_path=" + " -> ".join(TIE_PATH),
        f"score_hist={dict(sorted(Counter(len(v) for v in adj.values()).items()))}",
        f"directed_3cycles={directed_3cycles(adj)}",
        f"scc_sizes={scc_sizes(adj)}",
        f"hamiltonian_path_count={hamiltonian_path_count(adj)}",
        "selected_hamiltonian_path=" + " -> ".join(path),
        "",
        "Top carrier hooks:",
        *[
            f"- {next(c for c in CARRIERS if c.name == name).name}: "
            f"preserves={next(c for c in CARRIERS if c.name == name).preserves}; "
            f"next={next(c for c in CARRIERS if c.name == name).next_hook}"
            for name in path[:5]
        ],
    ]


def synthesis_report() -> List[str]:
    return [
        "1. CLAIMED MIDDLE QUOTIENT",
        "namespace=HYP-3133/T1200/LTI-261/LTT-159",
        "candidate invariant:",
        "  edge_extension_sandwich_certificate = (",
        "    four_sector_tetrahedral_word,",
        "    A000568_unlabeled_extension_shadow,",
        "    paired_tail_tip_child_deck,",
        "    SPEC_resonance_lattice_sidecar_or_named_debt",
        "  )",
        "",
        "Interpretation:",
        "- equinumerosity layer: sector words are weak compositions of m-2 outside vertices into four edge sectors.",
        "- equidistribution layer: A000568(m+1) is the unrooted one-extra-vertex tournament quotient.",
        "- equidecomposability layer: paired child signatures remember both tail-deletion and tip-deletion cards.",
        "",
        "Preserved LRC predicate: in the HYP-3125/HYP-3129 multi-far floor, the quotient is legal only if the",
        "R-safe -> Q-safe witness edge still routes to positive SPEC floor, finite constant chase, endpoint recursion,",
        "or a named H7/Asano/Lee-Yang/phi4/Cech debt.",
        "Destroyed by scalarization: the free middle extension class, endpoint owner, which deletion child improves",
        "the floor, and which finite low-frequency row in HYP-3129 caused the constant-chase burden.",
        "",
        "Challenged assumption: A000568 is not just numerology.  It is the natural unrooted middle quotient between",
        "sector equinumerosity and paired endpoint-deletion equidecomposition.",
    ]


def proof_route_report() -> List[str]:
    return [
        "4. LRC14 PROOF ROUTE CONSEQUENCES",
        "- HYP-3129 already moves the multi-far SPEC bound away from EH/BV and into exact-low plus Parseval-tail harmonic analysis.",
        "- HYP-3133 adds a finite quotient stratifier for that constant chase: sector word, A000568 middle shadow, then paired child deck.",
        "- HYP-3128 warns that naive Asano over the crowded R-tail fails; the A000568 shadow should therefore be a diagnostic, not a zero-free proof engine.",
        "- HYP-3125 edge-floor packets should add `a000568_extension_shadow` between `edge_tail_tip_sector_word` and deletion-child Rprime ratios.",
        "- Next executable test: for HYP-3129's worst finite SPEC rows, group rows by this sandwich certificate and ask whether every bad middle shadow either improves under one endpoint deletion or lands in a named resonance-lattice debt class.",
    ]


def main() -> None:
    print("HYP-3133 / codex-2026-06-27-A000568")
    print("A000568 middle quotient for edge-witness recursion; executable synthesis, not a proof.")
    print()
    for block in (synthesis_report(), sandwich_report(), tournament_report(), proof_route_report()):
        for line in block:
            print(line)
        print()


if __name__ == "__main__":
    main()
