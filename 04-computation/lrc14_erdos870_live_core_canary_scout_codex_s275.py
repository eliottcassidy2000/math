#!/usr/bin/env python3
"""HYP-3148 scout: Erdos-870 style filler/canary lens for n=4 tournaments.

Executable synthesis, not an LRC(14) proof.

The external trigger is davidturturean/erdos-870: a negative answer to
Erdos Problem #870 formalized in Lean.  The transferable proof pattern is:
large representation counts do not force a minimal subbasis; a proof should
separate an order-two core, deterministic fillers, canaries for shifted
exceptions, and a no-minimality/deletability audit.

Here we test the user's two n=4 tournament models exactly.  The objects are
colored Cayley tables of flip masks under xor; the colors are the four
unlabeled n=4 tournament classes:

  T: score sequence (0,1,2,3)
  S: score sequence (1,1,2,2)
  +: score sequence (0,2,2,2)
  -: score sequence (1,1,1,3)

Tournament Analysis declaration:
  vertices: proof-carrier lenses, not runners or raw arcs;
  pairwise observable: majority over retained-payload axes;
  switch/gauge: A -> B when A wins more axes; ties use the declared path.
"""

from __future__ import annotations

import itertools as it
import math
from collections import Counter
from dataclasses import dataclass
from typing import Dict, Iterable, List, Sequence, Tuple


V = tuple(range(4))
EDGES = tuple(it.combinations(V, 2))
CLASS_FROM_SCORE = {
    (0, 1, 2, 3): "T",
    (1, 1, 2, 2): "S",
    (0, 2, 2, 2): "+",
    (1, 1, 1, 3): "-",
}


def out_scores(bits: Sequence[int]) -> Tuple[int, ...]:
    out = [0] * 4
    for idx, (i, j) in enumerate(EDGES):
        if bits[idx] == 0:
            out[i] += 1
        else:
            out[j] += 1
    return tuple(sorted(out))


def cls(bits: Sequence[int]) -> str:
    return CLASS_FROM_SCORE[out_scores(bits)]


def partial_scores(assign: Dict[int, int]) -> Tuple[int, ...]:
    out = [0] * 4
    for idx, bit in assign.items():
        i, j = EDGES[idx]
        if bit == 0:
            out[i] += 1
        else:
            out[j] += 1
    return tuple(sorted(out))


def flip(bits: Sequence[int], indices: Iterable[int]) -> Tuple[int, ...]:
    out = list(bits)
    for i in indices:
        out[i] ^= 1
    return tuple(out)


def entropy(counts: Counter) -> float:
    total = sum(counts.values())
    return -sum((c / total) * math.log(c / total, 2) for c in counts.values())


def format_table(labels: Sequence[str], class_by_mask: Dict[int, str]) -> List[str]:
    rows = ["      " + " ".join(f"{x:>2}" for x in labels)]
    for i, a in enumerate(labels):
        vals = []
        for j, _b in enumerate(labels):
            vals.append(class_by_mask[i ^ j])
        rows.append(f"{a:>4} | " + " ".join(f"{x:>2}" for x in vals))
    return rows


def format_table_masks(labels: Sequence[str], masks: Sequence[int], class_by_mask: Dict[int, str]) -> List[str]:
    rows = ["      " + " ".join(f"{x:>2}" for x in labels)]
    for a, ma in zip(labels, masks):
        vals = [class_by_mask[ma ^ mb] for mb in masks]
        rows.append(f"{a:>4} | " + " ".join(f"{x:>2}" for x in vals))
    return rows


def scheme_a() -> Tuple[Dict[int, str], Counter, List[Tuple[str, Tuple[int, ...], Tuple[int, ...], str]]]:
    """Fixed Hamiltonian path 0->1->2->3, live skips a=(0,2), b=(1,3), c=(0,3)."""
    a = EDGES.index((0, 2))
    b = EDGES.index((1, 3))
    c = EDGES.index((0, 3))
    live = (a, b, c)
    base = (0, 0, 0, 0, 0, 0)
    names = ("E", "a", "b", "ab", "c", "ac", "bc", "abc")
    rows = []
    class_by_mask = {}
    for mask in range(8):
        fs = [live[k] for k in range(3) if (mask >> k) & 1]
        bits = flip(base, fs)
        class_by_mask[mask] = cls(bits)
        rows.append((names[mask], bits, out_scores(bits), cls(bits)))
    return class_by_mask, Counter(class_by_mask.values()), rows


def scheme_b_anchor() -> Tuple[Dict[int, str], Counter, Tuple[int, ...], Dict[int, int]]:
    """Canonical two-bit anchor: fix path plus c=(0,3); live x=a=(0,2), y=b=(1,3)."""
    x = EDGES.index((0, 2))
    y = EDGES.index((1, 3))
    fixed = {
        EDGES.index((0, 1)): 0,
        EDGES.index((0, 3)): 0,
        EDGES.index((1, 2)): 0,
        EDGES.index((2, 3)): 0,
    }
    base = [0] * 6
    live = (x, y)
    class_by_mask = {}
    for mask in range(4):
        fs = [live[k] for k in range(2) if (mask >> k) & 1]
        bits = flip(base, fs)
        class_by_mask[mask] = cls(bits)
    return class_by_mask, Counter(class_by_mask.values()), tuple(base), fixed


def two_bit_anchors() -> List[Tuple[Tuple[int, int], int, int, Tuple[int, ...], Dict[int, int]]]:
    out = []
    for free in it.combinations(range(6), 2):
        fixed_edges = [i for i in range(6) if i not in free]
        for fixed_bits in it.product((0, 1), repeat=4):
            fixed = dict(zip(fixed_edges, fixed_bits))
            if partial_scores(fixed) != (0, 1, 1, 2):
                continue
            for base_free in it.product((0, 1), repeat=2):
                base = [None] * 6
                for idx, bit in fixed.items():
                    base[idx] = bit
                for idx, bit in zip(free, base_free):
                    base[idx] = bit
                for x, y in (free, tuple(reversed(free))):
                    bx = flip(base, (x,))
                    by = flip(base, (y,))
                    bxy = flip(base, (x, y))
                    if (cls(base), cls(bx), cls(by), cls(bxy)) == ("T", "+", "-", "S"):
                        out.append((tuple(free), x, y, tuple(base), fixed))
    return out


def minimal_class_covers(live_count: int, class_by_mask: Dict[int, str]) -> List[Tuple[str, ...]]:
    names = ("a", "b", "c") if live_count == 3 else ("x", "y")
    covers = []
    for r in range(1, live_count + 1):
        for subset in it.combinations(range(live_count), r):
            reachable = {0}
            for mask in range(1, 1 << live_count):
                if all(((mask >> i) & 1) == 0 for i in range(live_count) if i not in subset):
                    reachable.add(mask)
            if {class_by_mask[m] for m in reachable} == {"T", "+", "-", "S"}:
                covers.append(tuple(names[i] for i in subset))
        if covers:
            break
    return covers


AXES = (
    "keeps_live_core",
    "names_filler",
    "detects_deletable_coordinate",
    "separates_class_from_rep_count",
    "preserves_edge_witness_payload",
    "lrc_packet_handoff",
    "lean_formalization_shape",
    "failure_guardrail",
)


@dataclass(frozen=True)
class Carrier:
    name: str
    scores: Tuple[int, ...]
    preserves: str
    destroys: str
    next_hook: str


CARRIERS = [
    Carrier(
        "two_bit_canary_anchor",
        (10, 10, 9, 10, 9, 9, 8, 9),
        "a minimal x/y core with one representation of each n=4 class",
        "the redundant tiling multiplicity that created the question",
        "promote every LRC quotient to a live-core plus filler/canary packet",
    ),
    Carrier(
        "erdos870_deletable_subbasis_audit",
        (8, 9, 10, 10, 8, 9, 10, 10),
        "the warning that many representations need not imply minimal proof support",
        "the concrete n=4 class table unless attached to a finite model",
        "attach deletable-coordinate tests to edge/fiber/U4 ledgers",
    ),
    Carrier(
        "fixed_path_tiling_cube",
        (8, 7, 8, 7, 8, 8, 7, 8),
        "the original Hamiltonian-path tiling model and its three free skips",
        "which skip is filler versus live core if only class counts are read",
        "freeze the long diagonal c when the proof needs a minimal core",
    ),
    Carrier(
        "edge_witness_tip_tail_packet",
        (7, 8, 8, 8, 10, 10, 7, 9),
        "HYP-3141's Tail/Tip/Orbit/Comm/Exit fields for directed edges",
        "the small F2 table if no anchor/filler coordinates are named",
        "add core/filler/canary words to edge_bounded_core_floor_exit",
    ),
    Carrier(
        "score_sequence_color_table",
        (7, 6, 5, 8, 6, 6, 8, 7),
        "the four n=4 classes as exact score-sequence colors",
        "representation multiplicity and live-coordinate minimality",
        "use only as a typed color after the flip mask is retained",
    ),
    Carrier(
        "raw_A000568_class_count",
        (4, 4, 3, 3, 3, 4, 5, 5),
        "that there are four unlabeled n=4 tournaments",
        "every coordinate needed for proof transfer",
        "negative control only",
    ),
]

TIE_PATH = [c.name for c in CARRIERS]


def tournament() -> Tuple[Counter, int, List[str]]:
    order = {name: i for i, name in enumerate(TIE_PATH)}
    score = Counter({c.name: 0 for c in CARRIERS})
    edges = set()
    for a, b in it.combinations(CARRIERS, 2):
        aw = sum(x > y for x, y in zip(a.scores, b.scores))
        bw = sum(y > x for x, y in zip(a.scores, b.scores))
        if aw > bw or (aw == bw and order[a.name] < order[b.name]):
            score[a.name] += 1
            edges.add((a.name, b.name))
        else:
            score[b.name] += 1
            edges.add((b.name, a.name))
    cycles = 0
    names = [c.name for c in CARRIERS]
    for a, b, c in it.combinations(names, 3):
        if (a, b) in edges and (b, c) in edges and (c, a) in edges:
            cycles += 1
        if (a, c) in edges and (c, b) in edges and (b, a) in edges:
            cycles += 1
    path = sorted(names, key=lambda n: -score[n])
    return score, cycles, path


def main() -> None:
    class_a, counts_a, rows_a = scheme_a()
    class_b, counts_b, base_b, fixed_b = scheme_b_anchor()
    anchors = two_bit_anchors()
    anchor_pair_hist = Counter(
        "disjoint" if len(set(EDGES[i][j] for i in free for j in (0, 1))) == 4 else "incident"
        for free, _x, _y, _base, _fixed in anchors
    )
    anchor_free_hist = Counter(tuple(EDGES[i] for i in free) for free, *_ in anchors)

    print("HYP-3148 / codex-2026-06-27-S275")
    print("Erdos-870 filler/canary lens for two n=4 tournament tables")
    print()
    print("1. EXTERNAL PATTERN IMPORT")
    print("Erdos-870 negative answer pattern:")
    print("  large representation lower bound does not force a minimal subbasis")
    print("  proof shape = order-two core + deterministic filler + clustered canary + deletability audit")
    print("LRC/tournament translation:")
    print("  class multiplicity is not proof support; keep live bits, filler bits, and deletion tests")
    print()
    print("2. N=4 CLASSES")
    for score, name in CLASS_FROM_SCORE.items():
        print(f"  {name}: score_sequence={score}")
    print()
    print("3. SCHEME A: FIXED HAMILTONIAN PATH, LIVE SKIPS a,b,c")
    for name, bits, score, color in rows_a:
        print(f"  {name:>3}: bits={bits} score={score} class={color}")
    print(f"  class_counts_full_cube={dict(counts_a)} entropy_bits={entropy(counts_a):.6f}")
    print(f"  minimal_all_class_covers={minimal_class_covers(3, class_a)}")
    print("  user_table_on_Eabc:")
    for line in format_table_masks(("E", "a", "b", "c"), (0, 1, 2, 4), class_a):
        print("   " + line)
    print()
    print("4. SCHEME B: FIX FOUR ARCS, LIVE x,y")
    print(f"  canonical_live: x=(0,2), y=(1,3); fixed={{{', '.join(str(EDGES[i])+':'+str(b) for i,b in fixed_b.items())}}}")
    print(f"  fixed_partial_score={partial_scores(fixed_b)} base_bits={base_b} base_class={cls(base_b)}")
    print(f"  class_counts_full_square={dict(counts_b)} entropy_bits={entropy(counts_b):.6f}")
    print(f"  minimal_all_class_covers={minimal_class_covers(2, class_b)}")
    print("  user_table_on_Exy:")
    for line in format_table(("E", "x", "y"), class_b):
        print("   " + line)
    print(f"  all_labelled_two_bit_anchors={len(anchors)}")
    print(f"  anchor_pair_type_hist={dict(anchor_pair_hist)}")
    print(f"  anchor_free_pair_hist={dict(anchor_free_hist)}")
    print()
    print("5. SYNTHESIS")
    print("  Scheme B is Scheme A with the long diagonal c frozen as deterministic filler.")
    print("  This turns the skewed class distribution T:+:-:S = 1:1:1:5 into 1:1:1:1.")
    print("  In Scheme A, c is class-cover-deletable because {a,b} already reaches T,+,-,S.")
    print("  In Scheme B, both live bits are load-bearing; deleting x or y loses class coverage.")
    print("  Erdos-870 warns that many representations can hide non-minimality; the LRC ledger should")
    print("  therefore record live_core_bits, filler_bits, canary_bits, deletable_coordinates,")
    print("  class_distribution, and terminal_exit_or_named_debt before using a tournament quotient.")
    print()
    print("6. TOURNAMENT ANALYSIS OVER PROOF CARRIERS")
    scores, cycles, path = tournament()
    print(f"  vertices={len(CARRIERS)} axes={','.join(AXES)}")
    print(f"  score_hist={dict(sorted(Counter(scores.values()).items()))}")
    print(f"  directed_3cycles={cycles}")
    print("  selected_hamiltonian_path=" + " -> ".join(path))
    print("  top_hooks:")
    for name in path[:4]:
        c = next(x for x in CARRIERS if x.name == name)
        print(f"   - {c.name}: preserves={c.preserves}; next={c.next_hook}")


if __name__ == "__main__":
    main()
