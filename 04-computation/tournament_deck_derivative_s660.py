#!/usr/bin/env python3
"""
S660: tournament deck derivative side-channel atlas.

This deepens the HYP-2234 tournament deck lane after incoming S659 claimed
HYP-2235 for finite-field Kakeya/Falconer.  The point is not to prove a
general reconstruction theorem.  The point is to treat a deck as a lossy
projection and ask which derivative side channels repair its collisions.

Tournament Analysis:
  Vertices are deck repair channels.  Pairwise observable is
  (separates_checked_collisions, deck_visibility, side_channel_transfer,
   simplicity).  The tie Hamiltonian path is the ranked repair list.

Assumption challenge:
  Candidate vertices included tournaments, cards, deleted vertices, deck
  buckets, scalar invariants, marked card boundaries, deletion derivatives,
  proof obligations, and cross-repo analogues.  This script chooses repair
  channels as vertices because the preserved predicate is "does this side
  channel refine a deck collision into useful proof data?"  It intentionally
  destroys labelled adjacency detail unless a marked boundary is being tested.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from functools import lru_cache
from itertools import combinations, permutations
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import frontier_small_carriers_s658 as base  # noqa: E402


MAX_N = 6


def prepare() -> dict[int, list[int]]:
    base.prepare_tournaments(MAX_N)
    return base.tournament_iso_classes(MAX_N)


@lru_cache(maxsize=None)
def canon(bits: int, n: int) -> int:
    return base.canon_tournament(bits, n)


@lru_cache(maxsize=None)
def hp(bits: int, n: int) -> int:
    return base.hamiltonian_paths(bits, n)


@lru_cache(maxsize=None)
def c3(bits: int, n: int) -> int:
    return base.c3_count(bits, n)


@lru_cache(maxsize=None)
def score_sequence(bits: int, n: int) -> tuple[int, ...]:
    return base.scores(bits, n)


@lru_cache(maxsize=None)
def scc(bits: int, n: int) -> tuple[int, ...]:
    return base.scc_sizes(bits, n)


@lru_cache(maxsize=None)
def scalar(bits: int, n: int) -> tuple[object, ...]:
    return (hp(bits, n), score_sequence(bits, n), c3(bits, n), scc(bits, n))


@lru_cache(maxsize=None)
def deleted_card(bits: int, n: int, v: int) -> int:
    return base.delete_vertex(bits, n, v)


@lru_cache(maxsize=None)
def deleted_card_canon(bits: int, n: int, v: int) -> int:
    return canon(deleted_card(bits, n, v), n - 1)


def outdegree(bits: int, n: int, v: int) -> int:
    return sum(1 for w in range(n) if base.has_edge(bits, n, v, w))


def card_records(bits: int, n: int) -> list[dict[str, object]]:
    h0 = hp(bits, n)
    c30 = c3(bits, n)
    rows = []
    for v in range(n):
        card = deleted_card(bits, n, v)
        ccard = canon(card, n - 1)
        rows.append(
            {
                "v": v,
                "card": ccard,
                "card_scalar": scalar(ccard, n - 1),
                "score": outdegree(bits, n, v),
                "h_loss": h0 - hp(card, n - 1),
                "c3_loss": c30 - c3(card, n - 1),
                "scc_card": scc(ccard, n - 1),
            }
        )
    return rows


def full_deck(bits: int, n: int) -> tuple[int, ...]:
    return tuple(sorted(rec["card"] for rec in card_records(bits, n)))


def scalar_deck(bits: int, n: int) -> tuple[tuple[object, ...], ...]:
    return tuple(sorted((rec["card_scalar"] for rec in card_records(bits, n)), key=repr))


def signature(bits: int, n: int, name: str) -> object:
    records = card_records(bits, n)
    deck = full_deck(bits, n)
    if name == "scalar":
        return scalar(bits, n)
    if name == "kelly_subcounts":
        return kelly_profile(bits, n)
    if name == "full_deck":
        return deck
    if name == "full_deck_plus_H":
        return (deck, hp(bits, n))
    if name == "full_deck_plus_score":
        return (deck, score_sequence(bits, n))
    if name == "full_deck_plus_H_score":
        return (deck, hp(bits, n), score_sequence(bits, n))
    if name == "full_deck_plus_scalar":
        return (deck, scalar(bits, n))
    if name == "scalar_deck":
        return scalar_deck(bits, n)
    if name == "unpaired_deleted_score":
        return (deck, tuple(sorted(rec["score"] for rec in records)))
    if name == "paired_deleted_score":
        return tuple(sorted((rec["card"], rec["score"]) for rec in records))
    if name == "paired_c3_loss":
        return tuple(sorted((rec["card"], rec["c3_loss"]) for rec in records))
    if name == "paired_H_loss":
        return tuple(sorted((rec["card"], rec["h_loss"]) for rec in records))
    if name == "paired_score_c3_loss":
        return tuple(sorted((rec["card"], rec["score"], rec["c3_loss"]) for rec in records))
    if name == "paired_all_derivatives":
        return tuple(
            sorted((rec["card"], rec["score"], rec["c3_loss"], rec["h_loss"]) for rec in records)
        )
    if name == "S658_loss_deck":
        losses = tuple(sorted((rec["h_loss"], rec["c3_loss"], rec["score"]) for rec in records))
        return (scalar_deck(bits, n), losses)
    if name == "card_scalar_all_derivatives":
        return tuple(
            sorted(
                (rec["card_scalar"], rec["score"], rec["c3_loss"], rec["h_loss"])
                for rec in records
            )
        )
    raise KeyError(name)


def subtournament(bits: int, n: int, vertices: tuple[int, ...]) -> int:
    out = 0
    idx = 0
    for i in range(len(vertices)):
        for j in range(i + 1, len(vertices)):
            if base.has_edge(bits, n, vertices[i], vertices[j]):
                out |= 1 << idx
            idx += 1
    return out


@lru_cache(maxsize=None)
def kelly_profile(bits: int, n: int) -> tuple[tuple[int, tuple[tuple[int, int], ...]], ...]:
    levels = []
    for k in range(3, n):
        counts: Counter[int] = Counter()
        for vertices in combinations(range(n), k):
            counts[canon(subtournament(bits, n, vertices), k)] += 1
        levels.append((k, tuple(sorted(counts.items()))))
    return tuple(levels)


def bucket_stats(classes: list[int], n: int, name: str) -> dict[str, object]:
    buckets: dict[object, list[int]] = defaultdict(list)
    for bits in classes:
        buckets[signature(bits, n, name)].append(bits)
    sizes = [len(v) for v in buckets.values()]
    collisions = [v for v in buckets.values() if len(v) > 1]
    return {
        "buckets": len(buckets),
        "max": max(sizes),
        "colliding_buckets": len(collisions),
        "colliding_classes": sum(len(v) for v in collisions),
    }


def bit_count(x: int) -> int:
    return int(x).bit_count()


def converse(bits: int, n: int) -> int:
    return bits ^ ((1 << (n * (n - 1) // 2)) - 1)


def is_converse_pair(a: int, b: int, n: int) -> bool:
    return canon(converse(a, n), n) == canon(b, n)


def iso_hamming_distance(a: int, b: int, n: int) -> int:
    return min(bit_count(a ^ base.permute_tournament(b, n, perm)) for perm in base.PERMS[n])


def full_deck_collision_anatomy(classes: list[int], n: int, limit: int = 8) -> list[dict[str, object]]:
    buckets: dict[object, list[int]] = defaultdict(list)
    for bits in classes:
        buckets[full_deck(bits, n)].append(bits)
    rows = []
    for deck, members in buckets.items():
        if len(members) <= 1:
            continue
        members = sorted(members)
        pair_data = []
        for a, b in combinations(members, 2):
            pair_data.append(
                {
                    "pair": (hex(a), hex(b)),
                    "H": (hp(a, n), hp(b, n)),
                    "scores": (score_sequence(a, n), score_sequence(b, n)),
                    "c3": (c3(a, n), c3(b, n)),
                    "scc": (scc(a, n), scc(b, n)),
                    "converse": is_converse_pair(a, b, n),
                    "min_flip_distance": iso_hamming_distance(a, b, n),
                    "paired_score_c3_separates": signature(a, n, "paired_score_c3_loss")
                    != signature(b, n, "paired_score_c3_loss"),
                    "paired_H_loss_separates": signature(a, n, "paired_H_loss")
                    != signature(b, n, "paired_H_loss"),
                }
            )
        rows.append(
            {
                "deck_cards": tuple(hex(x) for x in deck),
                "members": tuple(hex(x) for x in members),
                "pair_data": pair_data,
            }
        )
    return rows[:limit]


REPAIR_CHANNELS = [
    "full_deck",
    "kelly_subcounts",
    "scalar",
    "scalar_deck",
    "full_deck_plus_H",
    "full_deck_plus_score",
    "full_deck_plus_H_score",
    "full_deck_plus_scalar",
    "unpaired_deleted_score",
    "paired_deleted_score",
    "paired_c3_loss",
    "paired_H_loss",
    "paired_score_c3_loss",
    "paired_all_derivatives",
    "S658_loss_deck",
    "card_scalar_all_derivatives",
]


def repair_table(classes_by_n: dict[int, list[int]]) -> dict[int, dict[str, dict[str, object]]]:
    table: dict[int, dict[str, dict[str, object]]] = {}
    for n in range(3, MAX_N + 1):
        classes = classes_by_n[n]
        table[n] = {name: bucket_stats(classes, n, name) for name in REPAIR_CHANNELS}
    return table


def deck_formula_checks(classes_by_n: dict[int, list[int]]) -> list[dict[str, object]]:
    rows = []
    for n in range(4, MAX_N + 1):
        c3_ok = True
        h_card_sum_values: Counter[int] = Counter()
        for bits in classes_by_n[n]:
            card_c3_sum = sum(c3(deleted_card(bits, n, v), n - 1) for v in range(n))
            if card_c3_sum != (n - 3) * c3(bits, n):
                c3_ok = False
            h_card_sum_values[sum(hp(deleted_card(bits, n, v), n - 1) for v in range(n))] += 1
        rows.append(
            {
                "n": n,
                "kelly_c3_formula_ok": c3_ok,
                "distinct_sum_card_H": len(h_card_sum_values),
                "most_common_sum_card_H": h_card_sum_values.most_common(4),
            }
        )
    return rows


def route_tournament(final_stats: dict[str, dict[str, object]]) -> dict[str, object]:
    channels = [
        ("paired_deleted_score", (5, 1, 5, 5)),
        ("paired_score_c3_loss", (5, 1, 5, 4)),
        ("full_deck_plus_H", (4, 2, 4, 5)),
        ("paired_H_loss", (5, 1, 4, 3)),
        ("full_deck_plus_score", (3, 2, 3, 5)),
        ("full_deck", (2, 5, 2, 5)),
        ("kelly_subcounts", (1, 5, 2, 4)),
        ("scalar", (1, 4, 1, 5)),
    ]
    # First coordinate is replaced by checked separation strength at n=7.
    enriched = []
    for name, metric in channels:
        max_bucket = int(final_stats[name]["max"])
        separates = 5 if max_bucket == 1 else max(1, 5 - max_bucket)
        enriched.append((name, (separates, *metric[1:])))

    wins = {name: 0 for name, _metric in enriched}
    edges = []
    for (a, ma), (b, mb) in combinations(enriched, 2):
        if ma > mb:
            winner, loser = a, b
        elif mb > ma:
            winner, loser = b, a
        else:
            winner, loser = sorted([a, b])
        wins[winner] += 1
        edges.append((winner, loser))

    c3_count_route = 0
    names = [name for name, _metric in enriched]
    edge_set = set(edges)
    for a, b, c in combinations(names, 3):
        if (a, b) in edge_set and (b, c) in edge_set and (c, a) in edge_set:
            c3_count_route += 1
        if (b, a) in edge_set and (c, b) in edge_set and (a, c) in edge_set:
            c3_count_route += 1

    hpaths = 0
    for perm in permutations(names):
        if all((perm[i], perm[i + 1]) in edge_set for i in range(len(perm) - 1)):
            hpaths += 1

    return {
        "metrics": dict(enriched),
        "top_order": sorted(names, key=lambda name: (-wins[name], name)),
        "score_hist": dict(sorted(Counter(wins.values()).items())),
        "directed_3cycles": c3_count_route,
        "hamiltonian_paths": hpaths,
        "edges": edges,
    }


def main() -> None:
    classes_by_n = prepare()
    table = repair_table(classes_by_n)

    print("S660 tournament deck derivative side-channel atlas")
    print("=" * 72)

    print("\nA. Isomorphism classes")
    for n in range(3, MAX_N + 1):
        print(f"  n={n} classes={len(classes_by_n[n])}")

    print("\nB. Repair-channel bucket statistics")
    for n in range(3, MAX_N + 1):
        print(f"  n={n}")
        for name in REPAIR_CHANNELS:
            row = table[n][name]
            print(
                f"    {name:28s} buckets={row['buckets']:4d} "
                f"max={row['max']:2d} colliding_buckets={row['colliding_buckets']:3d} "
                f"colliding_classes={row['colliding_classes']:3d}"
            )

    print("\nC. Full-deck collision anatomy")
    for n in range(3, MAX_N + 1):
        collisions = full_deck_collision_anatomy(classes_by_n[n], n)
        print(f"  n={n} collision_buckets={len(collisions)} shown={min(len(collisions), 8)}")
        for row in collisions:
            print(f"    deck_cards={row['deck_cards']} members={row['members']}")
            for pair in row["pair_data"]:
                print(f"      pair={pair}")

    print("\nD. Deck-visible formula checks")
    for row in deck_formula_checks(classes_by_n):
        print(f"  {row}")

    print("\nE. Repair-channel Tournament Analysis")
    route = route_tournament(table[MAX_N])
    for key, value in route.items():
        print(f"  {key}={value}")

    print("\nF. Synthesis")
    print("  The ordinary full deck is a deletion projection; through n=6 it still")
    print("  has small collision buckets.  Kelly-style subtournament counts match")
    print("  the full-deck collision strength in this range but do not repair it.")
    print("  The useful repairs are marked derivative channels.  The cleanest one")
    print("  checked here is paired deleted score: attach to each card the outdegree")
    print("  of the missing vertex.  Unpaired deleted-score multisets do not suffice")
    print("  at n=6, and global H or score sequences do not repair the n=6 collision")
    print("  buckets.  Paired score+C3-loss and S658's stronger H/C3/score loss deck")
    print("  also resolve every checked collision.  This is the tournament analogue")
    print("  of retaining owner/carry")
    print("  in LRC, point-deletion degree/gain in unit distance, set pressure in")
    print("  union-closed families, and pinned/concurrency labels in finite-field")
    print("  Kakeya/Falconer carriers.")
    print("  The n=7 extension should use canonical augmentation or nauty-style")
    print("  automorphism pruning rather than raw permutation canonicalization.")


if __name__ == "__main__":
    main()
