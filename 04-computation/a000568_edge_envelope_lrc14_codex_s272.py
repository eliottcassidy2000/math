#!/usr/bin/env python3
"""HYP-3134 scout: A000568 inside the edge-witness envelope.

This is an executable synthesis scout, not a proof.

The prompt observation is:

  12 sits between 10 and 16,
  56 sits between 20 and 80.

HYP-3124 produced the two envelope sequences from directed-edge witnesses:

  lower sector-count envelope          1, 4, 10, 20, ...
  upper sector+paired-child envelope   1, 4, 16, 80, ...

The A000568 tournament class count appears one vertex later:

  U(n+1) sits between lower(n) and upper(n).

This scout extends the wedge to n=6 and treats the pattern as a controlled
forgetting principle for LRC14: raw sectors are too coarse, paired tail/tip
children are too fine, and the useful proof object is the global-consistency
quotient between them.

Tournament Analysis declaration:
  vertices: proof carriers / quotient operators, not runners or scalar counts;
  pairwise observable: majority comparison over retained proof-payload axes;
  switch/gauge: A beats B when it wins more axes; ties use the declared path.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from itertools import combinations, permutations
from math import comb
from typing import Dict, Iterable, List, Sequence, Tuple


AXES = (
    "a000568_wedge_explains",
    "lrc_predicate_retention",
    "controlled_forgetting",
    "tail_tip_payload",
    "global_gluing",
    "spec_floor_transfer",
    "finite_checkability",
    "failure_guard",
)


@dataclass(frozen=True)
class Carrier:
    name: str
    anchors: Tuple[str, ...]
    axes: Dict[str, int]
    preserves: str
    destroys: str
    next_hook: str


def carrier(
    name: str,
    anchors: Sequence[str],
    scores: Sequence[int],
    preserves: str,
    destroys: str,
    next_hook: str,
) -> Carrier:
    assert len(scores) == len(AXES)
    return Carrier(
        name=name,
        anchors=tuple(anchors),
        axes=dict(zip(AXES, scores)),
        preserves=preserves,
        destroys=destroys,
        next_hook=next_hook,
    )


CARRIERS: List[Carrier] = [
    carrier(
        "global_consistency_quotient",
        ("HYP-3134", "HYP-3133", "HYP-3132", "HYP-3124", "HYP-3054", "HYP-3106"),
        (10, 8, 10, 8, 10, 7, 8, 9),
        "the middle quotient between raw edge sectors and paired child signatures",
        "which local witnesses are equivalent unless the gluing rule is stated",
        "define the quotient map from paired edge children to global A000568-like packets",
    ),
    carrier(
        "paired_tail_tip_child_envelope",
        ("HYP-3124", "LTI-259", "LTT-157"),
        (9, 8, 9, 10, 7, 6, 8, 8),
        "four-sector deck plus both endpoint-deletion children",
        "global equivalences among locally distinct children",
        "use as the safe upper envelope before quotienting by observer/gluing legality",
    ),
    carrier(
        "A000568_middle_orbit_count",
        ("HYP-3047", "HYP-3054", "A000568"),
        (10, 5, 7, 5, 9, 3, 9, 7),
        "global isomorphism classes that sit inside the edge-witness envelope",
        "endpoint role, edge sector, and proof-route sidecars",
        "treat A000568 as a consistency check, not as a standalone proof invariant",
    ),
    carrier(
        "observer_extension_cut_payload",
        ("HYP-3054", "HYP-3059", "HYP-3065"),
        (8, 7, 9, 8, 9, 5, 8, 8),
        "the 12/48/56 cut-payload arithmetic and cross-sector orientation repair",
        "the paired child recursion unless joined to HYP-3124",
        "reuse the 56=P(5)+U(5)-U(4) splice as a gluing normal form",
    ),
    carrier(
        "resonance_lattice_SPEC_certificate",
        ("HYP-3129", "HYP-2861", "HYP-2606"),
        (5, 10, 8, 4, 7, 10, 8, 9),
        "the elementary SPEC floor through exact low modes plus Parseval tail",
        "edge-local witness geometry if it is only a Fourier scalar",
        "map edge-envelope gluing to resonance-lattice sidecar fields",
    ),
    carrier(
        "edge_floor_packet_schema",
        ("HYP-3125", "HYP-3127", "HYP-3129"),
        (7, 10, 9, 9, 8, 9, 8, 9),
        "R-safe -> Q-safe edge packet with deletion-child Rprime ratios and SPEC status",
        "global quotient minimality unless the A000568 wedge is retained",
        "add envelope_position and global_consistency_class to edge_floor_packet rows",
    ),
    carrier(
        "asano_lee_yang_dichotomy",
        ("HYP-3127", "HYP-3128", "HYP-3122"),
        (6, 8, 7, 7, 6, 8, 7, 10),
        "zero-free tip factors and the warning that naive joint Asano can fail",
        "the elementary SPEC floor if zero-free status is overpromoted",
        "separate Asano packaging from load-bearing SPEC certificate columns",
    ),
    carrier(
        "raw_four_sector_composition",
        ("HYP-3124", "HYP-3049"),
        (7, 3, 4, 5, 2, 2, 10, 4),
        "the lower envelope C(n+1,3) of possible edge-sector sizes",
        "paired endpoint children and global gluing",
        "keep only as the first compression alarm",
    ),
    carrier(
        "raw_count_numerology",
        ("negative-control",),
        (2, 1, 1, 1, 1, 1, 10, 1),
        "the visible numbers 10,12,16,20,56,80",
        "the proof predicate and all sidecars",
        "never accept without a carrier, quotient map, and terminal exit",
    ),
]

TIE_PATH = [c.name for c in CARRIERS]

A000568_SMALL = {
    1: 1,
    2: 1,
    3: 2,
    4: 4,
    5: 12,
    6: 56,
    7: 456,
}


def edge_index(i: int, j: int, n: int) -> int:
    if i > j:
        i, j = j, i
    idx = 0
    for a in range(n):
        for b in range(a + 1, n):
            if a == i and b == j:
                return idx
            idx += 1
    raise ValueError("bad edge")


def orient(mask: int, n: int, i: int, j: int) -> bool:
    """Return True iff i -> j in the labelled tournament mask."""
    if i < j:
        return bool((mask >> edge_index(i, j, n)) & 1)
    return not bool((mask >> edge_index(j, i, n)) & 1)


CANON_CACHE: Dict[Tuple[int, int], str] = {}


def canonical(mask: int, n: int) -> str:
    key = (mask, n)
    cached = CANON_CACHE.get(key)
    if cached is not None:
        return cached

    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    best = None
    for perm in permutations(range(n)):
        bits = []
        for a, b in pairs:
            old_a, old_b = perm[a], perm[b]
            bits.append("1" if orient(mask, n, old_a, old_b) else "0")
        word = "".join(bits)
        if best is None or word < best:
            best = word
    assert best is not None
    CANON_CACHE[key] = best
    return best


def delete_vertex(mask: int, n: int, deleted: int) -> int:
    verts = [v for v in range(n) if v != deleted]
    new_mask = 0
    bit_index = 0
    for ai in range(n - 1):
        for bi in range(ai + 1, n - 1):
            if orient(mask, n, verts[ai], verts[bi]):
                new_mask |= 1 << bit_index
            bit_index += 1
    return new_mask


def edge_signature_counts(n: int) -> Tuple[int, int, int, int]:
    sectors = set()
    tail_child = set()
    tip_child = set()
    child_pair = set()

    edge_count = n * (n - 1) // 2
    for mask in range(1 << edge_count):
        child_words = [canonical(delete_vertex(mask, n, v), n - 1) for v in range(n)]
        for tail in range(n):
            for tip in range(n):
                if tail == tip or not orient(mask, n, tail, tip):
                    continue
                deck = [0, 0, 0, 0]
                for w in range(n):
                    if w in (tail, tip):
                        continue
                    tail_to_w = orient(mask, n, tail, w)
                    tip_to_w = orient(mask, n, tip, w)
                    if tail_to_w and tip_to_w:
                        deck[0] += 1
                    elif tail_to_w and not tip_to_w:
                        deck[1] += 1
                    elif (not tail_to_w) and tip_to_w:
                        deck[2] += 1
                    else:
                        deck[3] += 1
                d = tuple(deck)
                sectors.add(d)
                tail_child.add((d, child_words[tail]))
                tip_child.add((d, child_words[tip]))
                child_pair.add((d, child_words[tail], child_words[tip]))
    return len(sectors), len(tail_child), len(tip_child), len(child_pair)


def unrooted_tournament_count(n: int) -> int:
    if n in A000568_SMALL:
        return A000568_SMALL[n]
    return len({canonical(mask, n) for mask in range(1 << (n * (n - 1) // 2))})


def compare(a: Carrier, b: Carrier) -> int:
    wins_a = wins_b = 0
    for axis in AXES:
        if a.axes[axis] > b.axes[axis]:
            wins_a += 1
        elif b.axes[axis] > a.axes[axis]:
            wins_b += 1
    if wins_a > wins_b:
        return 1
    if wins_b > wins_a:
        return -1
    return 1 if TIE_PATH.index(a.name) < TIE_PATH.index(b.name) else -1


def adjacency(carriers: Sequence[Carrier]) -> Dict[str, List[str]]:
    out = {carrier.name: [] for carrier in carriers}
    for a, b in combinations(carriers, 2):
        if compare(a, b) > 0:
            out[a.name].append(b.name)
        else:
            out[b.name].append(a.name)
    return out


def directed_three_cycles(out: Dict[str, List[str]]) -> List[Tuple[str, str, str]]:
    names = list(out)
    edge = {(a, b) for a, bs in out.items() for b in bs}
    cycles = []
    for a, b, c in combinations(names, 3):
        if (a, b) in edge and (b, c) in edge and (c, a) in edge:
            cycles.append((a, b, c))
        elif (a, c) in edge and (c, b) in edge and (b, a) in edge:
            cycles.append((a, c, b))
    return cycles


def scc_sizes(out: Dict[str, List[str]]) -> List[int]:
    reverse = {name: [] for name in out}
    for a, bs in out.items():
        for b in bs:
            reverse[b].append(a)

    seen = set()
    order: List[str] = []

    def dfs(v: str) -> None:
        seen.add(v)
        for w in out[v]:
            if w not in seen:
                dfs(w)
        order.append(v)

    for name in out:
        if name not in seen:
            dfs(name)

    seen.clear()
    sizes = []

    def rdfs(v: str, acc: List[str]) -> None:
        seen.add(v)
        acc.append(v)
        for w in reverse[v]:
            if w not in seen:
                rdfs(w, acc)

    for name in reversed(order):
        if name not in seen:
            acc: List[str] = []
            rdfs(name, acc)
            sizes.append(len(acc))
    return sorted(sizes, reverse=True)


def hamiltonian_path_count(out: Dict[str, List[str]]) -> int:
    names = list(out)
    n = len(names)
    index = {name: i for i, name in enumerate(names)}
    edge = [[False] * n for _ in range(n)]
    for a, bs in out.items():
        for b in bs:
            edge[index[a]][index[b]] = True

    dp = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for mask in range(1 << n):
        for last in range(n):
            count = dp[mask][last]
            if not count:
                continue
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                if edge[last][nxt]:
                    dp[mask | (1 << nxt)][nxt] += count
    return sum(dp[(1 << n) - 1])


def selected_path(out: Dict[str, List[str]]) -> List[str]:
    remaining = set(out)
    path = []
    while remaining:
        best = max(
            remaining,
            key=lambda name: (
                sum(1 for b in out[name] if b in remaining),
                -TIE_PATH.index(name),
            ),
        )
        path.append(best)
        remaining.remove(best)
    return path


def main() -> None:
    rows = []
    for n in range(2, 7):
        lower, tail, tip, upper = edge_signature_counts(n)
        unext = unrooted_tournament_count(n + 1)
        span = upper - lower
        position = None if span == 0 else (unext - lower) / span
        rows.append((n, lower, tail, tip, upper, unext, span, position))

    out = adjacency(CARRIERS)
    scores = {name: len(bs) for name, bs in out.items()}
    cycles = directed_three_cycles(out)

    print("HYP-3134 / codex-2026-06-27-S272")
    print("A000568 edge-envelope scout; executable synthesis, not a proof.")
    print()
    print("1. EDGE-SIGNATURE ENVELOPE")
    print("n lower=C(n+1,3) sector+tail sector+tip upper=sector+child_pair U(n+1) wedge_span wedge_position")
    for n, lower, tail, tip, upper, unext, span, position in rows:
        pos = "NA" if position is None else f"{position:.6f}"
        print(f"{n} {lower} {tail} {tip} {upper} {unext} {span} {pos}")
    print()
    print("Key rows:")
    print("  n=4: 10 < A000568(5)=12 < 16")
    print("  n=5: 20 < A000568(6)=56 < 80")
    print("  n=6: 35 < A000568(7)=456 < 632")
    print()

    print("2. INTERPRETATION")
    print("lower envelope=four-sector composition around a directed edge; it is local equinumerosity/equidistribution only")
    print("upper envelope=four-sector deck plus both endpoint-deletion children; it is the safe local two-ended witness")
    print("A000568(n+1)=global tournament classes; it sits between as a global-consistency quotient of local edge witnesses")
    print("LRC translation: raw scalar/SPEC sectors are too coarse, full edge packets are deliberately overcomplete, and the proof certificate should be the quotient where local edge witnesses glue globally without losing the LRC predicate.")
    print()

    print("3. TOURNAMENT ANALYSIS OVER QUOTIENT CARRIERS")
    print("vertices=proof carriers / quotient operators, not runners or scalar counts")
    print("pairwise_observable=majority over axes " + ",".join(AXES))
    print("switch=A->B when A wins more axes; ties use " + " -> ".join(TIE_PATH))
    print(f"score_hist={dict(sorted(Counter(scores.values()).items()))}")
    print(f"directed_3cycles={len(cycles)}")
    print(f"scc_sizes={scc_sizes(out)}")
    print(f"hamiltonian_path_count={hamiltonian_path_count(out)}")
    print("selected_hamiltonian_path=" + " -> ".join(selected_path(out)))
    print()
    print("Top carrier hooks:")
    for name in selected_path(out)[:6]:
        c = next(x for x in CARRIERS if x.name == name)
        print(f"- {name}: preserves={c.preserves}; next={c.next_hook}")
    print()

    print("4. PROOF-FRONTIER CONSEQUENCE")
    print("The count wedge suggests a proof tactic: start with the HYP-3124 paired child edge packet, then quotient only by a named global-consistency rule.  For LRC14 this means adding envelope_position, global_consistency_class, resonance_lattice_class, SPEC_bound_status, and edge_child_gluing_status to the HYP-3125 edge-floor packet.  HYP-3129 supplies the load-bearing SPEC certificate; the A000568 wedge supplies the controlled-forgetting discipline that prevents both under-refined scalar sectors and over-refined local child packets from becoming false proof endpoints.")


if __name__ == "__main__":
    main()
