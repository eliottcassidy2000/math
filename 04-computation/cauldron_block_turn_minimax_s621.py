#!/usr/bin/env python3
"""
cauldron_block_turn_minimax_s621.py

S621: adversarial block-turn cauldrons.

This is not the S619 attack-only parity game or the S620 two-block
schedule-channel game.  Here both players act on one shared active cauldron
pool under the all-boiled/removal rule.  The Delayer chooses placements to
maximize the final boil time; the Spoiler chooses placements to minimize it.

Schedules checked:

  - single-start: D, S, S, D, D, S, S, ...
  - double-start: D, D, S, S, D, D, ...

The single-start schedule matches the user's "the first turn only be one move"
normalization.  No resonance-energy scalar is used: the value is the literal
minimax final-boil time.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from itertools import combinations

from cauldron_game_s618 import RelationSpec, canonical, max_all_boiled, relation_witness


State = tuple[tuple[int, ...], ...]
Witness = tuple[tuple[int, ...], int]
Trace = tuple[tuple, ...]

DELAYER = "D"
SPOILER = "S"

WEAK = RelationSpec("two-term cauldron", 2, 2, False, "literal distinct values")
CLASSICAL = RelationSpec("two-term classical", 2, 2, True, "allows A=A")
TWO_THREE = RelationSpec("two-or-three-term cauldron", 2, 3, False)
FINITE = RelationSpec("finite-sums cauldron", 2, None, False)


@dataclass(frozen=True)
class BlockTurnResult:
    k: int
    spec: RelationSpec
    schedule: str
    last_boil_at: int | None
    trace: Trace
    states: int
    capped: bool = False


def controller_at(n: int, first_single: bool) -> str:
    """Return the controller for placement n."""
    if first_single:
        if n == 1:
            return DELAYER
        block = (n - 2) // 2
        return SPOILER if block % 2 == 0 else DELAYER
    block = (n - 1) // 2
    return DELAYER if block % 2 == 0 else SPOILER


def schedule_label(first_single: bool) -> str:
    if first_single:
        return "single-start D,S,S,D,D,..."
    return "double-start D,D,S,S,..."


def maxmin_all_boiled(
    k: int,
    spec: RelationSpec,
    first_single: bool = True,
    n_limit: int = 120,
    state_limit: int = 2_000_000,
) -> BlockTurnResult:
    """Solve the shared-pool removal game by exact minimax for small k."""
    state_count = 0
    capped = False

    @lru_cache(maxsize=None)
    def rec(n: int, active: State) -> tuple[int | None, Trace]:
        nonlocal state_count, capped
        state_count += 1
        if state_count > state_limit or n > n_limit:
            capped = True
            return None, (("CAP", n, state_count),)
        if not active:
            return n - 1, ()

        turn = controller_at(n, first_single)
        best_value: int | None = None
        best_trace: Trace = ()

        for idx, color in enumerate(active):
            candidate = tuple(sorted(color + (n,)))
            witness = relation_witness(candidate, spec)
            if witness is not None:
                new_active = list(active)
                new_active.pop(idx)
                next_state = canonical(new_active)
                event = ("boil", turn, n, candidate, witness)
            else:
                new_active = list(active)
                new_active[idx] = candidate
                next_state = canonical(new_active)
                event = ("keep", turn, n, candidate, None)

            value, suffix = rec(n + 1, next_state)
            if value is None:
                return None, suffix
            trace = (event,) + suffix

            if best_value is None:
                best_value = value
                best_trace = trace
                continue

            improves = value > best_value if turn == DELAYER else value < best_value
            if improves or (value == best_value and repr(trace) < repr(best_trace)):
                best_value = value
                best_trace = trace

        return best_value, best_trace

    value, trace = rec(1, tuple(() for _ in range(k)))
    return BlockTurnResult(k, spec, schedule_label(first_single), value, trace, state_count, capped)


def boil_events(trace: Trace) -> list[tuple]:
    return [event for event in trace if event[0] == "boil"]


def format_witness(witness: Witness) -> str:
    terms, target = witness
    return "+".join(map(str, terms)) + f"={target}"


def print_block_turn_table() -> None:
    print("ADVERSARIAL BLOCK-TURN MINIMAX")
    print("------------------------------")
    print("Shared active cauldron pool; a placement that creates the relation removes that cauldron.")
    print("D chooses placements maximizing the final boil time; S chooses placements minimizing it.")
    print("The single-start schedule is D,S,S,D,D,...; no resonance-energy scalar is used.")
    print()

    specs = [WEAK, CLASSICAL, TWO_THREE, FINITE]
    for first_single in [True, False]:
        print(f"SCHEDULE: {schedule_label(first_single)}")
        for spec in specs:
            print(f"{spec.label()}:")
            for k in [1, 2, 3]:
                result = maxmin_all_boiled(k, spec, first_single=first_single)
                cap = " [CAPPED]" if result.capped else ""
                print(
                    f"  k={k}: last_boil_at={result.last_boil_at}, "
                    f"states={result.states}{cap}"
                )
                events = boil_events(result.trace)
                if events:
                    compact = []
                    for _, turn, n, cauldron, witness in events:
                        compact.append(f"{turn}@{n}:{cauldron} via {format_witness(witness)}")
                    print("    boils: " + "; ".join(compact))
        print()


def print_comparison_reading() -> None:
    print("COMPARISON TO ONE-PLAYER REMOVAL")
    print("--------------------------------")
    print("The block-turn game is an alternating-quantifier policy on the S618 removal DP.")
    print("For k=3, S's first double move cuts the weak two-term final boil from 27 to 13.")
    print()

    rows = []
    for spec in [WEAK, CLASSICAL, FINITE]:
        one_player = max_all_boiled(3, spec)
        single = maxmin_all_boiled(3, spec, first_single=True)
        double = maxmin_all_boiled(3, spec, first_single=False)
        rows.append(
            (
                spec.name,
                one_player.last_boil_at,
                single.last_boil_at,
                double.last_boil_at,
            )
        )

    print("variant                         one-player  single-start  double-start")
    for name, one, single, double in rows:
        print(f"{name:30s} {one:10d} {single:13d} {double:13d}")
    print()
    print("Reading: the 'go twice' rule adds real strategic content, not just a")
    print("renaming of the one-player sacrifice search.  The single-start normalization")
    print("is especially hostile to the first player because S immediately receives the")
    print("first full two-move block.")
    print()


@dataclass(frozen=True)
class Route:
    name: str
    prompt_match: int
    exactness: int
    strategic_payload: int
    repo_connection: int
    computability: int
    note: str


ROUTES: tuple[Route, ...] = (
    Route(
        "single-start block minimax",
        5,
        5,
        5,
        5,
        4,
        "directly matches D,S,S,D,D,... over the shared removal pool",
    ),
    Route(
        "double-start block minimax",
        4,
        5,
        5,
        5,
        4,
        "checks the unnormalized D,D,S,S,... comparison",
    ),
    Route(
        "attack-only parity cauldrons",
        2,
        4,
        5,
        5,
        5,
        "S619 separate-pool odd/even attack game; useful but not this rule",
    ),
    Route(
        "gift-or-poison cauldrons",
        3,
        2,
        5,
        4,
        2,
        "future richer game where players may choose own or opponent pools",
    ),
    Route(
        "one-player removal DP",
        1,
        5,
        0,
        4,
        5,
        "baseline sacrifice dynamic without adversarial quantifiers",
    ),
    Route(
        "raw first-boil Schur",
        1,
        5,
        0,
        3,
        5,
        "first-boil threshold only; destroys removal timing",
    ),
    Route(
        "raw cauldron labels",
        0,
        1,
        0,
        1,
        5,
        "retains labels but not the quotient predicate that matters",
    ),
)


def route_score(route: Route) -> int:
    return (
        2 * route.prompt_match
        + 2 * route.strategic_payload
        + route.exactness
        + route.repo_connection
        + route.computability
    )


def build_route_edges(use_computability_only: bool = False) -> dict[tuple[str, str], int]:
    edges: dict[tuple[str, str], int] = {}
    tie_path = {route.name: idx for idx, route in enumerate(ROUTES)}
    for a, b in combinations(ROUTES, 2):
        score_a = a.computability if use_computability_only else route_score(a)
        score_b = b.computability if use_computability_only else route_score(b)
        if score_a == score_b:
            score_a += tie_path[b.name] - tie_path[a.name]
        winner, loser = (a.name, b.name) if score_a > score_b else (b.name, a.name)
        edges[(winner, loser)] = 1
    return edges


def strongly_connected_components(names: list[str], adjacency: dict[str, set[str]]) -> list[list[str]]:
    remaining = set(names)
    components: list[list[str]] = []
    while remaining:
        start = next(iter(remaining))
        forward = reachable(start, adjacency)
        reverse_adj = {name: set() for name in names}
        for src, targets in adjacency.items():
            for target in targets:
                reverse_adj[target].add(src)
        backward = reachable(start, reverse_adj)
        component = sorted(forward & backward)
        components.append(component)
        remaining -= set(component)
    return sorted(components, key=lambda comp: (-len(comp), comp))


def reachable(start: str, adjacency: dict[str, set[str]]) -> set[str]:
    seen = {start}
    stack = [start]
    while stack:
        node = stack.pop()
        for nxt in adjacency[node]:
            if nxt not in seen:
                seen.add(nxt)
                stack.append(nxt)
    return seen


def route_tournament() -> dict[str, object]:
    names = [route.name for route in ROUTES]
    edges = build_route_edges()
    alt_edges = build_route_edges(use_computability_only=True)
    adjacency = {name: set() for name in names}
    out_scores = {name: 0 for name in names}
    for winner, loser in edges:
        adjacency[winner].add(loser)
        out_scores[winner] += 1

    cycles3 = []
    for a, b, c in combinations(names, 3):
        if b in adjacency[a] and c in adjacency[b] and a in adjacency[c]:
            cycles3.append((a, b, c))
        if c in adjacency[a] and b in adjacency[c] and a in adjacency[b]:
            cycles3.append((a, c, b))

    index = {name: idx for idx, name in enumerate(names)}
    out_masks = [0] * len(names)
    for name in names:
        for target in adjacency[name]:
            out_masks[index[name]] |= 1 << index[target]

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

    edge_flips = 0
    for a, b in combinations(names, 2):
        base = (a, b) in edges
        alt = (a, b) in alt_edges
        if base != alt:
            edge_flips += 1

    ranking = sorted(ROUTES, key=lambda item: (-route_score(item), item.name))
    return {
        "ranking": ranking,
        "out_scores": out_scores,
        "score_hist": dict(sorted(score_hist.items())),
        "directed_3_cycles": cycles3,
        "sccs": strongly_connected_components(names, adjacency),
        "hamiltonian_path_count": hp_count,
        "edge_flips_vs_computability": edge_flips,
    }


def print_tournament_analysis() -> None:
    print("TOURNAMENT ANALYSIS")
    print("-------------------")
    print("Vertices: schedules/proof routes, not runners, arcs, or raw cauldrons.")
    print("Pairwise observable: which route preserves the requested adversarial block predicate.")
    print("Switch/gauge: 2*prompt_match + 2*strategic_payload + exactness + repo_connection + computability.")
    print("Tie Hamiltonian path: listed route order in the script.")
    analysis = route_tournament()
    print(f"Score histogram: {analysis['score_hist']}")
    print(f"Directed 3-cycles: {len(analysis['directed_3_cycles'])}")
    print(f"SCCs: {analysis['sccs']}")
    print(f"Hamiltonian-path count: {analysis['hamiltonian_path_count']}")
    print(f"Edge flips vs computability-only gauge: {analysis['edge_flips_vs_computability']}")
    print("Ranking:")
    for idx, route in enumerate(analysis["ranking"], start=1):
        print(f"  {idx}. {route.name} (score={route_score(route)}): {route.note}")
    print()


def print_assumption_challenge() -> None:
    print("ASSUMPTION CHALLENGE")
    print("--------------------")
    print("Alternate vertex sets considered: active cauldron states, controllers, two-move")
    print("blocks, schedules, boil witnesses, parity residue streams, proof routes, and")
    print("raw cauldron labels.")
    print("Preserved predicate: exact minimax final-boil time under the selected block")
    print("schedule and additive boil relation.")
    print("Destroyed information: canonical active-state quotients forget interchangeable")
    print("cauldron labels and all off-equilibrium non-optimal continuations; the trace")
    print("retains only one deterministic optimal line.")
    print("Challenged assumption: the user's block-turn adversarial game is not the same")
    print("as S619's attack-only parity game.  It is an alternating-quantifier policy")
    print("on the shared S618 removal DP.")
    print()


def main() -> None:
    print("S621 CAULDRON BLOCK-TURN MINIMAX")
    print("================================")
    print_block_turn_table()
    print_comparison_reading()
    print_tournament_analysis()
    print_assumption_challenge()


if __name__ == "__main__":
    main()
