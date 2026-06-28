#!/usr/bin/env python3
"""Creative LRC14 reframing scout: Bring/SC/BDH/Menger/charal recursion.

This is an exact, dependency-free scaffold over the mixed fibers exposed by
HYP-3406.  It does not use famous theorems as black boxes.  Instead it turns
their proof shapes into packet fields:

* Bring radical: a five-exit first-failure branch resolver.
* Schwarz-Christoffel: contact polygon turning word plus accessory debt.
* Barban-Davenport-Halberstam: finite channel-discrepancy variance.
* Menger: minimum owner-channel cuts separating theorem exits.
* charal signature: characteristic/chiral/arc-lift recursive sidecar.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations, permutations
import re


EXIT_ALPHABET = (
    "q-witness-discharge",
    "AP/GW-zero-open-equality",
    "unit-petal-named",
    "named-K33-state-lift",
    "positive-Haar-open",
)


@dataclass(frozen=True)
class Row:
    name: str
    exit: str
    route: str
    unit_slot: tuple[int, ...]
    owner_support: tuple[str, ...]
    v2: tuple[int, ...]


FIBERS: dict[str, tuple[Row, ...]] = {
    "height_leak_12_family": (
        Row(
            "GW-shell alias 12->132",
            "positive-Haar-open",
            "COVERING-MOMENT",
            (2, 2, 2),
            ("11:g1", "13:g1", "5:g1", "6:g2", "7:g7"),
            (0, 1, 1, 1, 2, 2, 3),
        ),
        Row(
            "single swap 12->48",
            "positive-Haar-open",
            "COVERING-MOMENT",
            (2, 2, 2),
            ("5:g1", "6:g2", "7:g7"),
            (0, 1, 1, 1, 2, 3, 4),
        ),
        Row(
            "P10+GW",
            "unit-petal-named",
            "BOUNDARY-PETAL-SPORADIC",
            (2, 2, 2),
            ("6:g2", "7:g7"),
            (0, 1, 1, 2, 2, 3, 3),
        ),
    ),
    "persistent_owner_leak_26_40_54_family": (
        Row(
            "single swap 1->26",
            "positive-Haar-open",
            "COVERING-MOMENT",
            (1, 2, 2),
            ("12:g2", "13:g1"),
            (0, 1, 1, 1, 1, 2, 2, 3),
        ),
        Row(
            "single swap 1->40",
            "positive-Haar-open",
            "COVERING-MOMENT",
            (1, 2, 2),
            ("12:g2", "13:g1", "2:g2"),
            (0, 1, 1, 1, 2, 2, 3, 3),
        ),
        Row(
            "single swap 1->54",
            "positive-Haar-open",
            "COVERING-MOMENT",
            (1, 2, 2),
            ("12:g2", "13:g1"),
            (0, 1, 1, 1, 1, 2, 2, 3),
        ),
        Row(
            "single swap 3->26",
            "positive-Haar-open",
            "COVERING-MOMENT",
            (2, 1, 2),
            ("11:g1", "12:g2", "13:g1", "6:g2"),
            (0, 1, 1, 1, 1, 2, 2, 3),
        ),
        Row(
            "single swap 3->40",
            "positive-Haar-open",
            "COVERING-MOMENT",
            (2, 1, 2),
            ("11:g1", "12:g2", "13:g1", "6:g2"),
            (0, 1, 1, 1, 2, 2, 3, 3),
        ),
        Row(
            "single swap 3->54",
            "positive-Haar-open",
            "COVERING-MOMENT",
            (2, 1, 2),
            ("11:g1", "12:g2", "13:g1", "6:g2"),
            (0, 1, 1, 1, 1, 2, 2, 3),
        ),
        Row(
            "single swap 5->26",
            "positive-Haar-open",
            "COVERING-MOMENT",
            (2, 2, 1),
            ("10:g2", "11:g1", "12:g2", "13:g1", "9:g1"),
            (0, 1, 1, 1, 1, 2, 2, 3),
        ),
        Row(
            "single swap 5->40",
            "positive-Haar-open",
            "COVERING-MOMENT",
            (2, 2, 1),
            ("10:g2", "11:g1", "12:g2", "13:g1", "9:g1"),
            (0, 1, 1, 1, 2, 2, 3, 3),
        ),
        Row(
            "single swap 5->54",
            "positive-Haar-open",
            "COVERING-MOMENT",
            (2, 2, 1),
            ("10:g2", "11:g1", "12:g2", "9:g1"),
            (0, 1, 1, 1, 1, 2, 2, 3),
        ),
        Row(
            "single swap 9->54",
            "positive-Haar-open",
            "COVERING-MOMENT",
            (2, 2, 1),
            ("10:g2", "11:g1", "12:g2", "13:g1", "5:g1", "7:g7", "8:g2"),
            (0, 1, 1, 1, 1, 2, 2, 3),
        ),
        Row(
            "petal 13->26",
            "unit-petal-named",
            "BOUNDARY-PETAL-SPORADIC",
            (1, 2, 2),
            ("12:g2", "1:g1"),
            (0, 1, 1, 1, 1, 2, 2, 3),
        ),
    ),
    "height_persistent_owner_leak_10_20_drop_add_family": (
        Row(
            "two drop(1, 10)->add(15, 20)",
            "positive-Haar-open",
            "COVERING-MOMENT",
            (2, 2, 2),
            ("1:g1", "2:g2", "6:g2", "7:g7"),
            (0, 1, 1, 2, 2, 2, 3),
        ),
        Row(
            "two drop(1, 10)->add(17, 20)",
            "positive-Haar-open",
            "COVERING-MOMENT",
            (1, 3, 2),
            ("13:g1", "2:g2", "3:g1", "6:g2"),
            (0, 1, 1, 2, 2, 2, 3),
        ),
        Row(
            "two drop(1, 10)->add(19, 20)",
            "positive-Haar-open",
            "COVERING-MOMENT",
            (1, 2, 3),
            ("13:g1", "2:g2", "5:g1", "6:g2", "7:g7"),
            (0, 1, 1, 2, 2, 2, 3),
        ),
        Row(
            "two drop(3, 10)->add(15, 20)",
            "positive-Haar-open",
            "COVERING-MOMENT",
            (3, 1, 2),
            ("11:g1", "13:g1", "6:g2", "7:g7"),
            (0, 1, 1, 2, 2, 2, 3),
        ),
        Row(
            "two drop(3, 10)->add(17, 20)",
            "positive-Haar-open",
            "COVERING-MOMENT",
            (2, 2, 2),
            ("13:g1", "6:g2"),
            (0, 1, 1, 2, 2, 2, 3),
        ),
        Row(
            "two drop(3, 10)->add(19, 20)",
            "positive-Haar-open",
            "COVERING-MOMENT",
            (2, 1, 3),
            ("11:g1", "5:g1", "6:g2", "7:g7"),
            (0, 1, 1, 2, 2, 2, 3),
        ),
        Row(
            "two drop(5, 10)->add(15, 20)",
            "positive-Haar-open",
            "COVERING-MOMENT",
            (3, 2, 1),
            ("11:g1", "12:g2", "13:g1", "1:g1", "6:g2", "7:g7", "9:g1"),
            (0, 1, 1, 2, 2, 2, 3),
        ),
        Row(
            "two drop(5, 10)->add(17, 20)",
            "positive-Haar-open",
            "COVERING-MOMENT",
            (2, 3, 1),
            ("11:g1", "13:g1", "3:g1", "6:g2", "9:g1"),
            (0, 1, 1, 2, 2, 2, 3),
        ),
        Row(
            "two drop(5, 10)->add(19, 20)",
            "positive-Haar-open",
            "COVERING-MOMENT",
            (2, 2, 2),
            ("11:g1", "12:g2", "13:g1", "5:g1", "6:g2", "7:g7"),
            (0, 1, 1, 2, 2, 2, 3),
        ),
        Row(
            "petal 10->20",
            "unit-petal-named",
            "BOUNDARY-PETAL-SPORADIC",
            (2, 2, 2),
            ("6:g2", "7:g7"),
            (0, 1, 1, 2, 2, 2, 3),
        ),
    ),
}


def owner_group(label: str) -> str:
    return label.split(":")[1]


def owner_index(label: str) -> int:
    return int(label.split(":")[0])


def sc_turn_word(row: Row) -> tuple[str, ...]:
    """Schwarz-Christoffel-style turn word around owner channels."""
    labels = sorted(row.owner_support, key=lambda x: (owner_index(x), owner_group(x)))
    out = []
    for label in labels:
        group = owner_group(label)
        if group == "g7":
            out.append("A")  # apex turn
        elif group == "g2":
            out.append("E")  # even-cover turn
        else:
            out.append("U")  # unit-owner turn
    return tuple(out)


def bdh_channel_vector(row: Row) -> tuple[tuple[str, int], ...]:
    counts = Counter(owner_group(x) for x in row.owner_support)
    return tuple(sorted(counts.items()))


def parse_replacement(name: str) -> tuple[str, int | None, int | None, int | None]:
    m = re.search(r"(\d+)->(\d+)", name)
    if not m:
        return ("named", None, None, None)
    old = int(m.group(1))
    new = int(m.group(2))
    return ("swap", old, new, new % 14)


def charal_signature(row: Row) -> tuple[object, ...]:
    """Characteristic/chiral/arc-lift signature for recursive comparisons."""
    kind, old, new, new_mod = parse_replacement(row.name)
    labels = tuple(sorted(row.owner_support, key=lambda x: (owner_group(x), owner_index(x))))
    return (
        ("replace", kind, old, new_mod),
        ("unit_slot", row.unit_slot),
        ("turn", sc_turn_word(row)),
        ("owner_groups", bdh_channel_vector(row)),
        ("apex", "7:g7" in row.owner_support),
        ("v2_span", min(row.v2), max(row.v2), sum(row.v2)),
        ("owner_labels", labels),
    )


def is_exit_pure(rows: tuple[Row, ...], key_fn) -> bool:
    fibers: dict[object, set[str]] = defaultdict(set)
    for row in rows:
        fibers[key_fn(row)].add(row.exit)
    return all(len(exits) == 1 for exits in fibers.values())


def min_owner_cut(rows: tuple[Row, ...]) -> tuple[int, list[tuple[str, ...]]]:
    labels = sorted({label for row in rows for label in row.owner_support})
    for r in range(1, len(labels) + 1):
        winners = []
        for subset in combinations(labels, r):
            if is_exit_pure(rows, lambda row, s=subset: tuple(x for x in s if x in row.owner_support)):
                winners.append(subset)
        if winners:
            return r, winners
    return 0, []


def channel_variance(rows: tuple[Row, ...]) -> list[tuple[Fraction, str]]:
    """Finite BDH analogue: channel indicator variance against theorem exits."""
    labels = sorted({label for row in rows for label in row.owner_support})
    exits = sorted({row.exit for row in rows})
    n = Fraction(len(rows), 1)
    scores = []
    for label in labels:
        total = Fraction(sum(label in row.owner_support for row in rows), 1)
        score = Fraction(0, 1)
        for exit_name in exits:
            exit_rows = [row for row in rows if row.exit == exit_name]
            observed = Fraction(sum(label in row.owner_support for row in exit_rows), 1)
            expected = total * Fraction(len(exit_rows), 1) / n
            score += (observed - expected) ** 2
        scores.append((score, label))
    return sorted(scores, reverse=True)


def recursive_ladders(rows: tuple[Row, ...]) -> dict[tuple[int, int], list[Row]]:
    ladders: dict[tuple[int, int], list[Row]] = defaultdict(list)
    for row in rows:
        kind, old, new, new_mod = parse_replacement(row.name)
        if kind == "swap" and old is not None and new_mod is not None:
            ladders[(old, new_mod)].append(row)
    for key in ladders:
        ladders[key].sort(key=lambda row: parse_replacement(row.name)[2] or 0)
    return dict(sorted(ladders.items()))


FEATURE_WEIGHTS = {
    "first_failure_branch": 11,
    "quintic_exit_alphabet": 7,
    "contact_polygon_turns": 9,
    "accessory_parameter_debt": 10,
    "channel_variance": 9,
    "owner_cut_separator": 12,
    "recursive_signature": 11,
    "height_wall": 8,
    "endpoint_owner": 10,
    "finite_bank_exactness": 10,
    "globalization_path": 8,
    "raw_analogy": 1,
}


@dataclass(frozen=True)
class Carrier:
    name: str
    preserves: frozenset[str]
    destroys: frozenset[str]
    priority: int

    def score(self) -> int:
        kept = sum(FEATURE_WEIGHTS[x] for x in self.preserves)
        lost = sum(max(1, FEATURE_WEIGHTS[x] // 3) for x in self.destroys)
        return kept - lost


def carriers() -> list[Carrier]:
    return [
        Carrier(
            "menger_owner_cut_recursion",
            frozenset(
                {
                    "owner_cut_separator",
                    "endpoint_owner",
                    "recursive_signature",
                    "finite_bank_exactness",
                    "globalization_path",
                    "first_failure_branch",
                }
            ),
            frozenset({"raw_analogy"}),
            90,
        ),
        Carrier(
            "charal_recursive_signature",
            frozenset(
                {
                    "recursive_signature",
                    "contact_polygon_turns",
                    "endpoint_owner",
                    "height_wall",
                    "finite_bank_exactness",
                    "first_failure_branch",
                }
            ),
            frozenset({"accessory_parameter_debt"}),
            85,
        ),
        Carrier(
            "bdh_channel_discrepancy_packet",
            frozenset(
                {
                    "channel_variance",
                    "endpoint_owner",
                    "finite_bank_exactness",
                    "globalization_path",
                    "first_failure_branch",
                }
            ),
            frozenset({"contact_polygon_turns", "height_wall"}),
            80,
        ),
        Carrier(
            "schwarz_christoffel_contact_polygon",
            frozenset(
                {
                    "contact_polygon_turns",
                    "accessory_parameter_debt",
                    "endpoint_owner",
                    "first_failure_branch",
                }
            ),
            frozenset({"channel_variance", "height_wall"}),
            70,
        ),
        Carrier(
            "bring_branch_resolvent",
            frozenset(
                {
                    "first_failure_branch",
                    "quintic_exit_alphabet",
                    "globalization_path",
                }
            ),
            frozenset({"endpoint_owner", "height_wall", "contact_polygon_turns"}),
            65,
        ),
        Carrier(
            "tropical_height_wall_backend",
            frozenset({"height_wall", "recursive_signature", "finite_bank_exactness"}),
            frozenset({"endpoint_owner", "channel_variance", "accessory_parameter_debt"}),
            55,
        ),
        Carrier(
            "raw_cross_discipline_analogy",
            frozenset({"raw_analogy"}),
            frozenset(
                {
                    "first_failure_branch",
                    "contact_polygon_turns",
                    "accessory_parameter_debt",
                    "channel_variance",
                    "owner_cut_separator",
                    "recursive_signature",
                    "height_wall",
                    "endpoint_owner",
                    "finite_bank_exactness",
                    "globalization_path",
                }
            ),
            5,
        ),
    ]


def orient(a: Carrier, b: Carrier) -> tuple[str, str]:
    ka = (a.score(), -len(a.destroys), a.priority, a.name)
    kb = (b.score(), -len(b.destroys), b.priority, b.name)
    return (a.name, b.name) if ka >= kb else (b.name, a.name)


def adjacency(cs: list[Carrier]) -> dict[str, set[str]]:
    adj = {c.name: set() for c in cs}
    for a, b in combinations(cs, 2):
        u, v = orient(a, b)
        adj[u].add(v)
    return adj


def directed_3cycles(adj: dict[str, set[str]]) -> list[tuple[str, str, str]]:
    out = []
    for a, b, c in combinations(sorted(adj), 3):
        if b in adj[a] and c in adj[b] and a in adj[c]:
            out.append((a, b, c))
        if c in adj[a] and b in adj[c] and a in adj[b]:
            out.append((a, c, b))
    return out


def hamiltonian_paths(adj: dict[str, set[str]]) -> list[tuple[str, ...]]:
    names = tuple(adj)
    out = []
    for perm in permutations(names):
        if all(perm[i + 1] in adj[perm[i]] for i in range(len(perm) - 1)):
            out.append(perm)
    return out


def print_row_signatures() -> None:
    print("\n## Row charal signatures")
    for fiber_name, rows in FIBERS.items():
        print(f"\n### {fiber_name}")
        for row in rows:
            print(f"{row.exit:19s} | {row.name:28s} | turn={sc_turn_word(row)} "
                  f"| groups={bdh_channel_vector(row)} | charal={charal_signature(row)}")


def print_cut_and_variance() -> None:
    print("\n## Menger cuts and finite BDH variance")
    for fiber_name, rows in FIBERS.items():
        cut_size, cuts = min_owner_cut(rows)
        print(f"\n### {fiber_name}")
        print(f"rows={len(rows)} exit_hist={dict(Counter(row.exit for row in rows))}")
        print(f"minimum_owner_label_cut_size={cut_size}")
        print(f"first_min_cuts={cuts[:8]}")
        print("top_channel_variance_labels=")
        for score, label in channel_variance(rows)[:8]:
            print(f"  {label}: {score}")


def print_recursive_patterns() -> None:
    print("\n## Recursive ladder patterns")
    rows = tuple(row for fiber in FIBERS.values() for row in fiber)
    for key, ladder in recursive_ladders(rows).items():
        names = [row.name for row in ladder]
        turns = [sc_turn_word(row) for row in ladder]
        owners = [row.owner_support for row in ladder]
        exits = [row.exit for row in ladder]
        print(f"ladder old_mod_key={key}: names={names}")
        print(f"  exits={exits}")
        print(f"  turns={turns}")
        print(f"  owners={owners}")


def print_theorem_routes() -> None:
    print("\n## Proof-route reframes")
    routes = [
        (
            "Bring branch resolver",
            "The five theorem exits form a quintic-like branch alphabet; "
            "normalization is useful only if the accessory owner/height fields "
            "make exit_status a function on packet fibers.",
        ),
        (
            "Schwarz-Christoffel polygon",
            "Unit contacts are polygon vertices; owner_support is the accessory "
            "parameter word.  A polygon turn word without accessory debt can hide "
            "the exact off-unit endpoint owner.",
        ),
        (
            "Barban-Davenport-Halberstam channel variance",
            "Use mean-square owner-channel imbalance on finite packet fibers as a "
            "route-pure discrepancy certificate; do not import prime-distribution "
            "theorems as black boxes.",
        ),
        (
            "Menger owner cut",
            "Treat endpoint-owner labels as a cut system separating theorem exits; "
            "a proof must show every mixed fiber has a bounded owner cut or a named "
            "residual debt.",
        ),
        (
            "Charal recursion",
            "Track characteristic/chiral/arc-lift signatures under +14 cusp "
            "ladders and endpoint deletion; persistent equality of this signature "
            "is the next candidate for a recursive no-new-kernel theorem.",
        ),
    ]
    for name, text in routes:
        print(f"- {name}: {text}")


def print_tournament() -> None:
    cs = carriers()
    adj = adjacency(cs)
    paths = hamiltonian_paths(adj)
    best_path = max(paths, key=lambda p: tuple(next(c for c in cs if c.name == x).score() for x in p))
    print("\n## Tournament Analysis")
    print("vertices=proof carriers / sidecar transformations, not runners")
    print("pairwise_observable=retained first-failure payload minus destroyed sidecars")
    print("switch_gauge=higher weighted retained payload; ties by fewer destroyed sidecars")
    print(f"vertex_count={len(cs)}")
    print(f"score_hist={dict(sorted(Counter(c.score() for c in cs).items()))}")
    print(f"directed_3cycles={len(directed_3cycles(adj))}")
    print(f"hamiltonian_path_count={len(paths)}")
    print("priority_hamiltonian_path=")
    for name in best_path:
        c = next(x for x in cs if x.name == name)
        print(f"  {name}: score={c.score()} keeps={sorted(c.preserves)} destroys={sorted(c.destroys)}")


def main() -> None:
    print("HYP-3410 Bring/Schwarz-Christoffel/BDH/Menger charal recursion scout")
    print("=" * 78)
    print("substrate=HYP-3406 expanded-bank mixed fibers")
    print(f"exit_alphabet={EXIT_ALPHABET}")
    print("charal=characteristic/chiral/arc-lift signature, not raw row identity")
    print_row_signatures()
    print_cut_and_variance()
    print_recursive_patterns()
    print_theorem_routes()
    print_tournament()
    print("\n## Assumption challenge")
    print("Considered vertices: runners, residues, owner labels, polygon turns, "
          "BDH channels, Menger cuts, Bring branches, deletion events, and proof obligations.")
    print("Chosen vertices: proof carriers and sidecar transformations.  The quotient "
          "preserves LRC route status only when exit_status is constant on fibers or "
          "owner/height/accessory debt is restored.")


if __name__ == "__main__":
    main()
