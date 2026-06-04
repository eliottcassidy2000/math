#!/usr/bin/env python3
"""
cauldron_complexity_variants_s619.py

S619: difficulty and variant atlas for the cauldron game.

This script does three things:

1. Measures exact small last-boil state counts and a bounded k=4 pressure test.
2. Measures first-boil frontier growth, showing where Schur-frontier explosion
   starts before any A000568-style graph isomorphism enters.
3. Runs a tiny adversarial "attack-only" cauldron game, where players take
   turns placing the next natural number into each other's cauldrons.

Tournament Analysis is over problem variants/proof routes rather than raw
cauldrons.  This keeps the quotient predicate explicit.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from itertools import combinations
from math import log10
import sys

from cauldron_game_s618 import RelationSpec, canonical, max_all_boiled, relation_witness


State = tuple[tuple[int, ...], ...]


WEAK = RelationSpec("two-term cauldron", 2, 2, False, "literal distinct values")
CLASSICAL = RelationSpec("two-term classical", 2, 2, True, "allows A=A")
FINITE = RelationSpec("finite-sums cauldron", 2, None, False)
TWO_THREE = RelationSpec("two-or-three-term cauldron", 2, 3, False)


def first_boil_frontier_counts(
    k: int,
    spec: RelationSpec,
    n_limit: int,
    cap: int = 250_000,
) -> list[tuple[int, int, bool]]:
    """Count canonical safe color-frontiers after placing each n."""
    states: set[State] = {canonical(tuple(() for _ in range(k)))}
    rows: list[tuple[int, int, bool]] = []
    for n in range(1, n_limit + 1):
        next_states: set[State] = set()
        for state in states:
            for idx, color in enumerate(state):
                candidate = tuple(sorted(color + (n,)))
                if relation_witness(candidate, spec) is not None:
                    continue
                new_state = list(state)
                new_state[idx] = candidate
                next_states.add(canonical(new_state))
                if len(next_states) > cap:
                    rows.append((n, len(next_states), True))
                    return rows
        states = next_states
        rows.append((n, len(states), False))
        if not states:
            return rows
    return rows


def frontier_milestones(rows: list[tuple[int, int, bool]], wanted: list[int]) -> str:
    by_n = {n: (count, capped) for n, count, capped in rows}
    parts: list[str] = []
    for n in wanted:
        if n in by_n:
            count, capped = by_n[n]
            marker = "+" if capped else ""
            parts.append(f"n={n}:{count}{marker}")
    if rows and rows[-1][2]:
        n, count, _ = rows[-1]
        if n not in wanted:
            parts.append(f"n={n}:{count}+")
    return ", ".join(parts)


def place_in(state: State, idx: int, n: int, spec: RelationSpec) -> State:
    colors = [tuple(c) for c in state]
    candidate = tuple(sorted(colors[idx] + (n,)))
    if relation_witness(candidate, spec) is not None:
        del colors[idx]
        return canonical(colors)
    colors[idx] = candidate
    return canonical(colors)


@dataclass(frozen=True)
class AdversarialResult:
    spec: RelationSpec
    a_cauldrons: int
    b_cauldrons: int
    winner: str
    at_n: int
    states: int


def attack_only_game(
    a_cauldrons: int,
    b_cauldrons: int,
    spec: RelationSpec,
    cap: int = 40,
) -> AdversarialResult:
    """Player A attacks B on odd n; player B attacks A on even n."""

    @lru_cache(maxsize=None)
    def rec(n: int, turn: int, a_state: State, b_state: State) -> tuple[int, int]:
        if not a_state:
            return -1, n - 1
        if not b_state:
            return 1, n - 1
        if n > cap:
            return 0, cap

        target = b_state if turn == 0 else a_state
        choices: list[tuple[int, int]] = []
        for idx in range(len(target)):
            if turn == 0:
                choices.append(rec(n + 1, 1, a_state, place_in(b_state, idx, n, spec)))
            else:
                choices.append(rec(n + 1, 0, place_in(a_state, idx, n, spec), b_state))

        if turn == 0:
            wins = [choice for choice in choices if choice[0] > 0]
            if wins:
                return min(wins, key=lambda item: item[1])
            draws = [choice for choice in choices if choice[0] == 0]
            if draws:
                return 0, cap
            return max(choices, key=lambda item: item[1])

        wins = [choice for choice in choices if choice[0] < 0]
        if wins:
            return min(wins, key=lambda item: item[1])
        draws = [choice for choice in choices if choice[0] == 0]
        if draws:
            return 0, cap
        return max(choices, key=lambda item: item[1])

    start_a = canonical(tuple(() for _ in range(a_cauldrons)))
    start_b = canonical(tuple(() for _ in range(b_cauldrons)))
    value, at_n = rec(1, 0, start_a, start_b)
    winner = "A" if value > 0 else "B" if value < 0 else "draw"
    return AdversarialResult(spec, a_cauldrons, b_cauldrons, winner, at_n, rec.cache_info().currsize)


@dataclass(frozen=True)
class Variant:
    name: str
    difficulty: int
    repo_connection: int
    novelty: int
    computability: int
    proof_payload: int
    note: str


VARIANTS: tuple[Variant, ...] = (
    Variant(
        "attack-only adversarial parity cauldrons",
        5,
        5,
        5,
        4,
        5,
        "turn order splits odd/even attack streams; a correlated residue can beat density",
    ),
    Variant(
        "gift-or-poison adversarial cauldrons",
        5,
        4,
        5,
        3,
        5,
        "players may place n into own or opponent cauldrons; PSPACE-style finite game",
    ),
    Variant(
        "modular CRT cauldrons",
        4,
        5,
        4,
        5,
        4,
        "cauldrons live in Z/mZ and boil on zero-sum/subset-sum residues",
    ),
    Variant(
        "LRC depth cauldrons",
        5,
        5,
        4,
        2,
        5,
        "a cauldron boils when its danger arcs cover the observer circle, p0=0",
    ),
    Variant(
        "unit-distance geometric cauldrons",
        5,
        5,
        4,
        2,
        5,
        "points are placed into cauldrons; boil on unit edges/triangles or unfaithful embeddings",
    ),
    Variant(
        "OCF odd-cycle cauldrons",
        4,
        5,
        4,
        4,
        4,
        "odd cycles are ingredients; boil on forbidden independence-polynomial values",
    ),
    Variant(
        "hidden-quotient A000568 cauldrons",
        5,
        4,
        4,
        3,
        4,
        "state is observed only up to relabelling/orbit type, forcing Burnside-style side data",
    ),
    Variant(
        "base last-boil Schur removal",
        3,
        3,
        3,
        5,
        3,
        "finite active-state DP over unlabeled sum-free color classes",
    ),
    Variant(
        "finite-sums frontier cauldrons",
        4,
        4,
        3,
        4,
        4,
        "pair sums become subset-sum bitmasks, an all-orders overlap toy model",
    ),
    Variant(
        "upper-cap certificate cauldrons",
        4,
        4,
        4,
        4,
        4,
        "a boil creates a reusable certificate bounding all later safe placements",
    ),
    Variant(
        "Collatz residue cauldrons",
        4,
        5,
        4,
        3,
        5,
        "ingredients are 2/3 residue blocks; boil on correlated blind-density residues",
    ),
    Variant(
        "misere cauldrons",
        3,
        3,
        4,
        5,
        3,
        "players try to force their own cauldrons to boil last/first under reversed payoff",
    ),
)


def variant_score(variant: Variant) -> int:
    return (
        2 * variant.repo_connection
        + 2 * variant.proof_payload
        + variant.difficulty
        + variant.novelty
        + variant.computability
    )


def variant_tournament() -> dict[str, object]:
    edges: dict[tuple[str, str], int] = {}
    for a, b in combinations(VARIANTS, 2):
        score_a = variant_score(a)
        score_b = variant_score(b)
        if score_a == score_b:
            score_a += -len(a.name)
            score_b += -len(b.name)
        winner, loser = (a, b) if score_a > score_b else (b, a)
        edges[(winner.name, loser.name)] = 1

    names = [variant.name for variant in VARIANTS]
    adjacency = {name: set() for name in names}
    out_scores = {name: 0 for name in names}
    for winner, loser in edges:
        adjacency[winner].add(loser)
        out_scores[winner] += 1

    cycles3 = 0
    for a, b, c in combinations(names, 3):
        if b in adjacency[a] and c in adjacency[b] and a in adjacency[c]:
            cycles3 += 1
        if c in adjacency[a] and b in adjacency[c] and a in adjacency[b]:
            cycles3 += 1

    index = {name: idx for idx, name in enumerate(names)}
    out_masks = [0] * len(names)
    for name in names:
        mask = 0
        for target in adjacency[name]:
            mask |= 1 << index[target]
        out_masks[index[name]] = mask

    @lru_cache(maxsize=None)
    def hp_dp(mask: int, last: int) -> int:
        if mask == (1 << last):
            return 1
        prior_mask = mask ^ (1 << last)
        total = 0
        for prior in range(len(names)):
            if (prior_mask >> prior) & 1 and ((out_masks[prior] >> last) & 1):
                total += hp_dp(prior_mask, prior)
        return total

    full_mask = (1 << len(names)) - 1
    hp_count = sum(hp_dp(full_mask, last) for last in range(len(names)))

    score_hist: dict[int, int] = {}
    for score in out_scores.values():
        score_hist[score] = score_hist.get(score, 0) + 1

    return {
        "ranking": sorted(VARIANTS, key=lambda item: (-variant_score(item), item.name)),
        "out_scores": out_scores,
        "score_hist": dict(sorted(score_hist.items())),
        "cycles3": cycles3,
        "hamiltonian_path_count": hp_count,
    }


def print_difficulty_verdict() -> None:
    print("DIFFICULTY VERDICT")
    print("------------------")
    print("Base fixed-k last-boil is Schur-frontier hard, not A000568/LRC/unit-distance hard.")
    print("  Reason: the state is a canonical tuple of live sum-free color classes; no S_n graph")
    print("  isomorphism, no continuous p0 coimage, and no geometric embedding side-channel is present.")
    print("  But k as a growing parameter inherits Schur/Rado-style frontier explosion.")
    print("Adversarial, hidden-quotient, LRC, and geometric cauldron variants can be as hard")
    print("  as the repo's larger problems because they add alternating quantifiers, orbit")
    print("  side-data, all-orders coverage, or embedding obstructions.")
    print()


def print_one_player_evidence() -> None:
    print("ONE-PLAYER LAST-BOIL EVIDENCE")
    print("-----------------------------")
    for spec in [WEAK, CLASSICAL, FINITE]:
        row = []
        for k in [1, 2, 3]:
            result = max_all_boiled(k, spec)
            row.append(f"k={k}:last={result.last_boil_at},states={result.states}")
        print(f"{spec.label()}: " + "; ".join(row))
    bounded = max_all_boiled(4, WEAK, n_limit=90, state_limit=250_000)
    print(
        f"bounded k=4 weak-removal pressure test: last={bounded.last_boil_at}, "
        f"states={bounded.states}, capped={bounded.capped}"
    )
    print()

    print("FIRST-BOIL FRONTIER GROWTH")
    print("--------------------------")
    rows3 = first_boil_frontier_counts(3, WEAK, 24)
    rows4 = first_boil_frontier_counts(4, WEAK, 35)
    print("weak two-term k=3 frontiers:", frontier_milestones(rows3, [5, 10, 15, 20, 23, 24]))
    print("weak two-term k=4 frontiers:", frontier_milestones(rows4, [5, 10, 12, 13]))
    if rows4[-1][2]:
        n, count, _ = rows4[-1]
        print(f"k=4 crossed the {count - 1} frontier cap already at n={n}; no last-boil value claimed.")
    print()


def print_adversarial_evidence() -> None:
    print("ADVERSARIAL ATTACK-ONLY GAME")
    print("----------------------------")
    print("Rule: on odd n, A places n into one of B's cauldrons; on even n, B places n into one of A's.")
    print("Terminal: if your cauldrons are all boiled, you lose.")
    for spec in [WEAK, TWO_THREE, FINITE]:
        print(f"{spec.label()}:")
        for a, b in [(1, 1), (1, 2), (2, 1), (2, 2)]:
            result = attack_only_game(a, b, spec)
            print(
                f"  A={a},B={b}: winner={result.winner}, at_n={result.at_n}, states={result.states}"
            )
    print("Reading: with two-term sums, A attacks only with odd numbers and cannot make")
    print("  an odd target from odd+odd, while B's even stream can make 2+4=6 and beyond.")
    print("  Allowing three/finite terms lets A use 1+3+5=9, flipping the asymmetric 2-vs-1 case.")
    print()


def print_variant_atlas() -> None:
    print("VARIANT ATLAS / TOURNAMENT ANALYSIS")
    print("-----------------------------------")
    print("Vertices: cauldron variants and proof routes.")
    print("Pairwise observable: weighted score = 2*repo_connection + 2*proof_payload +")
    print("  difficulty + novelty + computability; ties use shorter certificate names.")
    print("Switch/gauge: prefer variants that retain side-channel payload while still being computable.")
    analysis = variant_tournament()
    print(f"Score histogram: {analysis['score_hist']}")
    print(f"Directed 3-cycles: {analysis['cycles3']}")
    print(f"Hamiltonian-path count: {analysis['hamiltonian_path_count']}")
    print("Ranking:")
    for idx, variant in enumerate(analysis["ranking"], start=1):
        print(
            f"  {idx:2d}. {variant.name} "
            f"(score={variant_score(variant)}, log10 score={log10(variant_score(variant)):.2f})"
        )
        print(f"      {variant.note}")
    print()


def main() -> None:
    sys.stdout.reconfigure(line_buffering=True)
    print("S619 CAULDRON COMPLEXITY AND VARIANT ATLAS")
    print("==========================================")
    print_difficulty_verdict()
    print_one_player_evidence()
    print_adversarial_evidence()
    print_variant_atlas()


if __name__ == "__main__":
    main()
