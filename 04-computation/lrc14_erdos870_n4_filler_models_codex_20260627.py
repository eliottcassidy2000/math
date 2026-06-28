#!/usr/bin/env python3
"""HYP-3145 scout: Erdos-870 filler interfaces and two n=4 tournament models.

This is a finite modeling scout, not a proof.  It compares two ways to present
the four unlabeled tournaments on four vertices:

1. fixed Hamiltonian path / tiling cube: three off-path arcs a,b,c remain;
2. partial-score filler interface: four arcs are fixed with partial score
   profile (0,1,1,2), leaving two core arcs x,y.

The first model is a useful atlas but the quotient to isomorphism classes is
not a congruence: the S class has many bit representatives.  The second model
is an exact four-state interface: the two remaining arcs form a Klein square
over the classes T,+,-,S.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from itertools import combinations, product


PAIRS = [(i, j) for i in range(4) for j in range(i + 1, 4)]
CLASS_BY_SCORE = {
    (0, 1, 2, 3): "T",
    (1, 1, 2, 2): "S",
    (0, 2, 2, 2): "+",
    (1, 1, 1, 3): "-",
}


def scores(bits: int) -> tuple[int, ...]:
    out = [0] * 4
    for k, (i, j) in enumerate(PAIRS):
        if (bits >> k) & 1:
            out[i] += 1
        else:
            out[j] += 1
    return tuple(sorted(out))


def class_name(bits: int) -> str:
    return CLASS_BY_SCORE[scores(bits)]


def word_from_mask(mask: int, arc_labels: dict[int, str]) -> str:
    if mask == 0:
        return "E"
    labels = [label for bit, label in arc_labels.items() if mask & bit]
    return "".join(labels)


def fixed_path_model() -> dict[str, object]:
    """Use base path 0->1->2->3 and transitive orientation as E."""
    transitive = (1 << len(PAIRS)) - 1
    # Label the three non-path arcs so the user's one-flip names are matched.
    labels = {
        1 << 1: "a",  # arc (0,2), one flip gives +
        1 << 4: "b",  # arc (1,3), one flip gives -
        1 << 2: "c",  # arc (0,3), one flip gives S
    }
    generators = [0, 1 << 1, 1 << 4, 1 << 2]
    table = []
    for row in generators:
        table_row = []
        for col in generators:
            table_row.append(class_name(transitive ^ row ^ col))
        table.append(table_row)

    fiber: dict[str, list[str]] = {name: [] for name in ["T", "+", "-", "S"]}
    nonpath_bits = list(labels)
    for choices in product([0, 1], repeat=3):
        mask = 0
        for bit, flag in zip(nonpath_bits, choices):
            if flag:
                mask |= bit
        fiber[class_name(transitive ^ mask)].append(word_from_mask(mask, labels))

    ambiguity = {}
    for cls, words in fiber.items():
        images = {}
        for gen in nonpath_bits:
            gen_label = labels[gen]
            images[gen_label] = sorted(
                {class_name(transitive ^ mask ^ gen) for mask_word in words for mask in [mask_from_word(mask_word, labels)]}
            )
        ambiguity[cls] = images

    return {
        "base_bits": transitive,
        "labels": labels,
        "generators": generators,
        "table": table,
        "fiber": fiber,
        "ambiguity": ambiguity,
    }


def mask_from_word(word: str, labels: dict[int, str]) -> int:
    if word == "E":
        return 0
    out = 0
    reverse = {v: k for k, v in labels.items()}
    for ch in word:
        out |= reverse[ch]
    return out


def find_partial_score_model() -> dict[str, object]:
    target = ("T", "+", "-", "S")
    for free in combinations(range(len(PAIRS)), 2):
        fixed = [k for k in range(len(PAIRS)) if k not in free]
        for fixed_vals in product([0, 1], repeat=len(fixed)):
            partial = [0] * 4
            fixed_bits = 0
            for k, value in zip(fixed, fixed_vals):
                if value:
                    fixed_bits |= 1 << k
                i, j = PAIRS[k]
                partial[i if value else j] += 1
            if tuple(sorted(partial)) != (0, 1, 1, 2):
                continue
            for base_free in product([0, 1], repeat=2):
                base = fixed_bits
                for k, value in zip(free, base_free):
                    if value:
                        base |= 1 << k
                    else:
                        base &= ~(1 << k)
                x_bit, y_bit = (1 << free[0], 1 << free[1])
                classes = (
                    class_name(base),
                    class_name(base ^ x_bit),
                    class_name(base ^ y_bit),
                    class_name(base ^ x_bit ^ y_bit),
                )
                if classes == target:
                    return {
                        "base_bits": base,
                        "free": free,
                        "x_bit": x_bit,
                        "y_bit": y_bit,
                        "fixed": fixed,
                        "fixed_bits": fixed_bits,
                        "partial": tuple(sorted(partial)),
                    }
    raise RuntimeError("no partial-score model found")


def table_for_states(base: int, masks: list[int]) -> list[list[str]]:
    return [[class_name(base ^ row ^ col) for col in masks] for row in masks]


def print_table(title: str, headers: list[str], table: list[list[str]]) -> None:
    print(title)
    print("      " + " ".join(f"{h:>3}" for h in headers))
    for h, row in zip(headers, table):
        print(f"{h:>3} | " + " ".join(f"{x:>3}" for x in row))


@dataclass(frozen=True)
class Carrier:
    name: str
    exactness: int
    quotient_congruence: int
    filler_core_split: int
    deletion_alarm: int
    lrc_transfer: int
    formal_interface: int
    overclaim_safety: int


CARRIERS = [
    Carrier("partial_score_two_arc_core", 10, 10, 9, 6, 8, 9, 8),
    Carrier("erdos870_filler_interface", 9, 9, 10, 9, 8, 10, 9),
    Carrier("quotient_congruence_audit", 8, 10, 8, 9, 9, 9, 10),
    Carrier("fixed_path_tiling_cube", 9, 4, 5, 8, 7, 6, 7),
    Carrier("nonminimal_deletable_fiber_alarm", 7, 5, 7, 10, 8, 7, 9),
    Carrier("edge_witness_SPEC_packet", 7, 8, 8, 8, 10, 8, 8),
    Carrier("raw_n4_class_table", 5, 3, 3, 2, 3, 4, 4),
]


def score(carrier: Carrier) -> int:
    return (
        2 * carrier.exactness
        + 3 * carrier.quotient_congruence
        + 2 * carrier.filler_core_split
        + 2 * carrier.deletion_alarm
        + 3 * carrier.lrc_transfer
        + 2 * carrier.formal_interface
        + 2 * carrier.overclaim_safety
    )


def tournament_fingerprint() -> None:
    names = [c.name for c in CARRIERS]
    wins = Counter()
    edges = []
    tie_path = {name: i for i, name in enumerate(names)}
    for i, a in enumerate(CARRIERS):
        for j, b in enumerate(CARRIERS):
            if i >= j:
                continue
            sa, sb = score(a), score(b)
            if sa > sb or (sa == sb and tie_path[a.name] < tie_path[b.name]):
                winner, loser = a.name, b.name
            else:
                winner, loser = b.name, a.name
            wins[winner] += 1
            edges.append((winner, loser))
    score_hist = Counter(wins.get(name, 0) for name in names)
    # Directed 3-cycles by brute force.
    edge_set = set(edges)
    cycles = 0
    for a, b, c in combinations(names, 3):
        if (
            ((a, b) in edge_set and (b, c) in edge_set and (c, a) in edge_set)
            or ((a, c) in edge_set and (c, b) in edge_set and (b, a) in edge_set)
        ):
            cycles += 1
    order = sorted(names, key=lambda n: (-wins.get(n, 0), tie_path[n]))
    print("TOURNAMENT ANALYSIS")
    print(f"vertices={len(names)}")
    print(f"score_hist={dict(sorted(score_hist.items()))}")
    print(f"directed_3cycles={cycles}")
    print("selected_path=" + " -> ".join(order))


def main() -> None:
    print("HYP-3145 / codex-2026-06-27")
    print("ERDOS-870 FILLER INTERFACES AND TWO N=4 TOURNAMENT MODELS")
    print()
    print("Class names:")
    for score_tuple, name in sorted(CLASS_BY_SCORE.items(), key=lambda item: item[1]):
        print(f"  {name}: score={score_tuple}")

    model1 = fixed_path_model()
    headers1 = ["E", "a", "b", "c"]
    print()
    print("MODEL 1: FIXED HAMILTONIAN PATH / THREE-FREE-ARC TILING")
    print("Base path 0->1->2->3, transitive E, with a=(0,2), b=(1,3), c=(0,3).")
    print_table("Representative xor table:", headers1, model1["table"])
    print("Full fixed-path fiber over all 2^3 off-path states:")
    for cls in ["T", "+", "-", "S"]:
        print(f"  {cls}: {model1['fiber'][cls]}")
    print("Noncongruence alarm: multiplying from an S representative is ambiguous.")
    for gen, images in model1["ambiguity"]["S"].items():
        print(f"  S * {gen} can land in {images}")

    model2 = find_partial_score_model()
    x_bit = model2["x_bit"]
    y_bit = model2["y_bit"]
    print()
    print("MODEL 2: PARTIAL-SCORE FILLER / TWO-FREE-ARC CORE")
    print(f"Fixed partial score profile: {model2['partial']}")
    print("One exact realization:")
    for k in model2["fixed"]:
        i, j = PAIRS[k]
        direction = f"{i}->{j}" if (model2["fixed_bits"] >> k) & 1 else f"{j}->{i}"
        print(f"  fixed arc {PAIRS[k]} as {direction}")
    print(f"  x toggles arc {PAIRS[model2['free'][0]]}")
    print(f"  y toggles arc {PAIRS[model2['free'][1]]}")
    print(f"  E full score={scores(model2['base_bits'])}")
    headers2 = ["E", "x", "y"]
    print_table("User-facing partial table:", headers2, table_for_states(model2["base_bits"], [0, x_bit, y_bit]))
    headers2_full = ["E", "x", "y", "xy"]
    print_table(
        "Closed Klein-square table, with xy representing S:",
        headers2_full,
        table_for_states(model2["base_bits"], [0, x_bit, y_bit, x_bit ^ y_bit]),
    )

    print()
    print("ERDOS-870 TRANSFER")
    print("- Treat the four fixed arcs as finite fillers that force the residue/score scaffold.")
    print("- Treat x,y as the order-two core; all four class states are reached without quotient ambiguity.")
    print("- The fixed-path tiling cube is an atlas, but its S fiber is nonminimal/deletable data.")
    print("- LRC14 analogue: choose filler sidecars first, then leave a small signed core for SPEC/resolvent checks.")
    print("- Formalization lesson: expose the interface boundary before proving the terminal theorem.")
    print()
    tournament_fingerprint()


if __name__ == "__main__":
    main()
