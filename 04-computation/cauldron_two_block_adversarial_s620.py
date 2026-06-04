#!/usr/bin/env python3
"""
cauldron_two_block_adversarial_s620.py

S620: two-block adversarial cauldrons.

The user-proposed schedule is:

    A, B, B, A, A, B, B, A, A, ...

The first turn is a one-move handicap for A, then the players alternate
two-move blocks.  This script compares that schedule with ordinary parity
alternation in the attack-only cauldron game from S619.

Tournament Analysis is over schedules and proof routes, not over raw cauldron
labels.  The pairwise observable is retained proof payload: residue balance,
small-game depth, repo connection, and computability.  The quotient preserves
who can force a boil under a schedule and destroys most ingredient histories
except the explicit additive witnesses found by the cauldron predicate.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from itertools import combinations_with_replacement
from typing import Callable

from cauldron_game_s618 import RelationSpec, canonical, relation_witness


State = tuple[tuple[int, ...], ...]
TurnFn = Callable[[int], int]


WEAK = RelationSpec("two-term cauldron", 2, 2, False, "literal distinct values")
TWO_THREE = RelationSpec("two-or-three-term cauldron", 2, 3, False)
FINITE = RelationSpec("finite-sums cauldron", 2, None, False)


@dataclass(frozen=True)
class Schedule:
    name: str
    word: str
    turn: TurnFn
    note: str


def parity_turn(n: int) -> int:
    """A attacks on odd n, B attacks on even n."""
    return 0 if n % 2 == 1 else 1


def two_block_handicap_turn(n: int) -> int:
    """A, B, B, A, A, ...; A residues 0,1 mod 4 and B residues 2,3."""
    return 0 if n % 4 in (0, 1) else 1


def two_block_no_handicap_turn(n: int) -> int:
    """A, A, B, B, A, A, ...; included as the strict block baseline."""
    return 0 if n % 4 in (1, 2) else 1


SCHEDULES: tuple[Schedule, ...] = (
    Schedule(
        "parity alternation",
        "A,B,A,B,...",
        parity_turn,
        "pure odd/even attack streams",
    ),
    Schedule(
        "one-handicap two-block",
        "A,B,B,A,A,B,B,...",
        two_block_handicap_turn,
        "period-four square wave with a one-move A start",
    ),
    Schedule(
        "strict two-block",
        "A,A,B,B,A,A,B,B,...",
        two_block_no_handicap_turn,
        "same period-four square wave without the first-turn handicap",
    ),
)


def player_name(turn: int) -> str:
    return "A" if turn == 0 else "B"


def residues_for(schedule: Schedule, player: int, modulus: int = 4) -> tuple[int, ...]:
    period = 4 if modulus == 4 else modulus
    return tuple(sorted({n % modulus for n in range(1, period + 1) if schedule.turn(n) == player}))


def closure_residues(residues: tuple[int, ...], modulus: int, min_terms: int, max_terms: int) -> tuple[int, ...]:
    hits: set[int] = set()
    for arity in range(min_terms, max_terms + 1):
        for terms in combinations_with_replacement(residues, arity):
            hits.add(sum(terms) % modulus)
    return tuple(sorted(hits))


def schedule_residue_report(schedule: Schedule) -> list[str]:
    lines: list[str] = []
    for player in (0, 1):
        stream = residues_for(schedule, player)
        pair = closure_residues(stream, 4, 2, 2)
        two_three = closure_residues(stream, 4, 2, 3)
        pair_self_hits = tuple(sorted(set(stream) & set(pair)))
        missing_pair_targets = tuple(sorted(set(stream) - set(pair)))
        lines.append(
            f"{player_name(player)} stream {stream}: pair sums {pair}, "
            f"self-hit {pair_self_hits}, pair-missing {missing_pair_targets}, "
            f"2/3-term closure {two_three}"
        )
    return lines


def place_in(state: State, idx: int, n: int, spec: RelationSpec) -> State:
    colors = [tuple(c) for c in state]
    candidate = tuple(sorted(colors[idx] + (n,)))
    if relation_witness(candidate, spec) is not None:
        del colors[idx]
        return canonical(colors)
    colors[idx] = candidate
    return canonical(colors)


@dataclass(frozen=True)
class GameResult:
    schedule: str
    spec: str
    a_cauldrons: int
    b_cauldrons: int
    winner: str
    at_n: int
    states: int
    cap: int


def choose_for_a(choices: list[tuple[int, int]], cap: int) -> tuple[int, int]:
    wins = [choice for choice in choices if choice[0] > 0]
    if wins:
        return min(wins, key=lambda item: item[1])
    draws = [choice for choice in choices if choice[0] == 0]
    if draws:
        return 0, cap
    return max(choices, key=lambda item: item[1])


def choose_for_b(choices: list[tuple[int, int]], cap: int) -> tuple[int, int]:
    wins = [choice for choice in choices if choice[0] < 0]
    if wins:
        return min(wins, key=lambda item: item[1])
    draws = [choice for choice in choices if choice[0] == 0]
    if draws:
        return 0, cap
    return max(choices, key=lambda item: item[1])


def attack_only_game(
    schedule: Schedule,
    a_cauldrons: int,
    b_cauldrons: int,
    spec: RelationSpec,
    cap: int = 48,
) -> GameResult:
    """Players place n into an opponent cauldron; boiled cauldrons disappear."""

    @lru_cache(maxsize=None)
    def rec(n: int, a_state: State, b_state: State) -> tuple[int, int]:
        if not a_state:
            return -1, n - 1
        if not b_state:
            return 1, n - 1
        if n > cap:
            return 0, cap

        mover = schedule.turn(n)
        choices: list[tuple[int, int]] = []
        if mover == 0:
            for idx in range(len(b_state)):
                choices.append(rec(n + 1, a_state, place_in(b_state, idx, n, spec)))
            return choose_for_a(choices, cap)

        for idx in range(len(a_state)):
            choices.append(rec(n + 1, place_in(a_state, idx, n, spec), b_state))
        return choose_for_b(choices, cap)

    start_a = canonical(tuple(() for _ in range(a_cauldrons)))
    start_b = canonical(tuple(() for _ in range(b_cauldrons)))
    value, at_n = rec(1, start_a, start_b)
    winner = "A" if value > 0 else "B" if value < 0 else f"draw-to-{cap}"
    return GameResult(
        schedule.name,
        spec.name,
        a_cauldrons,
        b_cauldrons,
        winner,
        at_n,
        rec.cache_info().currsize,
        cap,
    )


@dataclass(frozen=True)
class Route:
    name: str
    residue_balance: int
    depth_signal: int
    repo_connection: int
    computability: int
    novelty: int
    note: str

    def score(self) -> int:
        return self.residue_balance + self.depth_signal + self.repo_connection + self.computability + self.novelty


ROUTES: tuple[Route, ...] = (
    Route(
        "residue-channel sumset calculus",
        5,
        4,
        5,
        5,
        4,
        "turn schedules are subsets S of Z/4Z and attack power is S+S hitting S",
    ),
    Route(
        "automatic/Walsh schedule analysis",
        5,
        3,
        5,
        5,
        4,
        "A,BB,AA is a period-four square wave / Walsh bit",
    ),
    Route(
        "Maker-Breaker Schur hypergraph game",
        4,
        5,
        4,
        3,
        5,
        "cauldrons are bins in a positional game on triples (a,b,a+b)",
    ),
    Route(
        "Collatz two-block residue analogy",
        5,
        4,
        5,
        3,
        4,
        "density is blind when moves arrive in correlated residue blocks",
    ),
    Route(
        "LRC depth/Helly cauldron analogy",
        3,
        4,
        5,
        2,
        4,
        "boil events mimic retained overlap order rather than scalar density",
    ),
    Route(
        "adaptive adversarial cauldrons",
        3,
        5,
        4,
        2,
        5,
        "players choose own/opponent and schedule control becomes PSPACE-like",
    ),
    Route(
        "raw cauldron-label minimax",
        2,
        3,
        2,
        4,
        2,
        "exact for small k but forgets the reusable residue mechanism",
    ),
)


def route_tournament(routes: tuple[Route, ...]) -> tuple[dict[str, int], int, int, list[tuple[str, str]]]:
    wins = {route.name: 0 for route in routes}
    flips: list[tuple[str, str]] = []
    cycles = 0
    for i, left in enumerate(routes):
        for right in routes[i + 1 :]:
            if left.score() >= right.score():
                wins[left.name] += 1
                flips.append((left.name, right.name))
            else:
                wins[right.name] += 1
                flips.append((right.name, left.name))

    for i in range(len(routes)):
        for j in range(i + 1, len(routes)):
            for k in range(j + 1, len(routes)):
                trio = (routes[i], routes[j], routes[k])
                names = [route.name for route in trio]
                edge_set = set(flips)
                cyclic = (
                    (names[0], names[1]) in edge_set
                    and (names[1], names[2]) in edge_set
                    and (names[2], names[0]) in edge_set
                ) or (
                    (names[1], names[0]) in edge_set
                    and (names[2], names[1]) in edge_set
                    and (names[0], names[2]) in edge_set
                )
                if cyclic:
                    cycles += 1

    route_index = {route.name: idx for idx, route in enumerate(routes)}
    ordered = sorted(routes, key=lambda route: (-route.score(), route_index[route.name]))
    hamiltonian_paths = 1 if cycles == 0 else 0
    return wins, cycles, hamiltonian_paths, [(ordered[i].name, ordered[i + 1].name) for i in range(len(ordered) - 1)]


def print_schedule_reports() -> None:
    print("TWO-BLOCK ADVERSARIAL CAULDRONS S620")
    print()
    print("Schedule residue filters mod 4")
    for schedule in SCHEDULES:
        print(f"- {schedule.name}: {schedule.word} ({schedule.note})")
        for line in schedule_residue_report(schedule):
            print(f"  {line}")
    print()


def print_game_table() -> list[GameResult]:
    specs = (WEAK, TWO_THREE, FINITE)
    pairs = ((1, 1), (1, 2), (2, 1), (2, 2))
    results: list[GameResult] = []
    print("Exact attack-only minimax table")
    print("schedule | relation | A,B | winner | last/forced n | states")
    print("--- | --- | --- | --- | --- | ---")
    for schedule in SCHEDULES:
        for spec in specs:
            for a_cauldrons, b_cauldrons in pairs:
                result = attack_only_game(schedule, a_cauldrons, b_cauldrons, spec)
                results.append(result)
                print(
                    f"{result.schedule} | {result.spec} | "
                    f"{a_cauldrons},{b_cauldrons} | {result.winner} | "
                    f"{result.at_n} | {result.states}"
                )
    print()
    return results


def print_interpretation(results: list[GameResult]) -> None:
    print("Interpretation")
    print("- Parity alternation gives A the odd stream and B the even stream.")
    print("  In the weak two-term rule, odd+odd is even and cannot hit A's odd target stream; even+even remains even, so B has a first-order residue advantage.")
    print("- The one-handicap two-block schedule turns the streams into {0,1} and {2,3} mod 4.")
    print("  Now each side has at least one pair-sum self-hit residue, so the game is no longer decided by density alone.")
    print("- The schedule is a period-four Walsh/Rademacher square wave. Generalizing the word turns the problem into additive combinatorics of a schedule set S: can r in S be represented by allowed sums from S?")
    print("- As a positional game, the board is the Schur hypergraph with edges (a,b,a+b); cauldrons are surviving bins and a boil is a claimed hyperedge.")
    print("- This is the cauldron analogue of the Collatz two-block lesson: correlated residue blocks can carry information that scalar density erases.")
    print()

    weak = [result for result in results if result.spec == WEAK.name and result.a_cauldrons == 2 and result.b_cauldrons == 1]
    print("Focused 2v1 weak-sum comparison")
    for result in weak:
        print(f"- {result.schedule}: {result.winner} at {result.at_n} using {result.states} cached states")
    print()


def print_tournament_analysis() -> None:
    wins, cycles, hamiltonian_paths, path_edges = route_tournament(ROUTES)
    print("Tournament Analysis over proof routes")
    print("Route | score | outdegree | note")
    print("--- | --- | --- | ---")
    route_index = {route.name: idx for idx, route in enumerate(ROUTES)}
    for route in sorted(ROUTES, key=lambda item: (-item.score(), route_index[item.name])):
        print(f"{route.name} | {route.score()} | {wins[route.name]} | {route.note}")
    print()
    print(f"directed 3-cycles: {cycles}")
    print(f"Hamiltonian path count in this score quotient: {hamiltonian_paths}")
    print("tie Hamiltonian path:")
    for left, right in path_edges:
        print(f"- {left} -> {right}")
    print()
    print("Challenged assumption")
    print("- Do not use cauldrons, players, or turns as the only vertices.")
    print("  The useful quotient vertices here are residue channels, schedule languages, hypergraph proof obligations, and retained witness types.")


def main() -> None:
    print_schedule_reports()
    results = print_game_table()
    print_interpretation(results)
    print_tournament_analysis()


if __name__ == "__main__":
    main()
