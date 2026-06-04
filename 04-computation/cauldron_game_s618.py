#!/usr/bin/env python3
"""
cauldron_game_s618.py

Finite scouts for the cauldron game:

  - First-boil variant: place 1,2,3,... into k active cauldrons. A placement is
    legal only if the chosen cauldron remains free of the additive relation.
  - All-boiled/removal variant: a placement that creates the additive relation
    removes that cauldron; play continues until every cauldron has boiled.

The literal rule "a cauldron contains three values A,B,C with C=A+B" uses
distinct summand values, since each natural number is placed once. The script
also reports the repeated-summand comparison, which is the classical Schur
number convention.

Tournament Analysis is included, but the vertices are proof/computation lenses
and active-state quotients, not runners/arcs/cauldrons by default.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from itertools import combinations, combinations_with_replacement, permutations
from typing import Iterable


@dataclass(frozen=True)
class RelationSpec:
    name: str
    min_terms: int = 2
    max_terms: int | None = 2
    repeated_terms: bool = False
    note: str = ""

    def label(self) -> str:
        max_label = "unbounded" if self.max_terms is None else str(self.max_terms)
        repeat_label = "repeated" if self.repeated_terms else "distinct"
        return f"{self.name} ({repeat_label}, terms {self.min_terms}..{max_label})"


@dataclass(frozen=True)
class FirstBoilResult:
    k: int
    spec: RelationSpec
    safe_prefix: int
    forced_at: int | None
    coloring: tuple[tuple[int, ...], ...] | None
    states: int
    capped: bool = False


@dataclass(frozen=True)
class RemovalResult:
    k: int
    spec: RelationSpec
    last_boil_at: int | None
    trace: tuple[tuple, ...]
    states: int
    capped: bool = False


def canonical(colors: Iterable[Iterable[int]]) -> tuple[tuple[int, ...], ...]:
    return tuple(sorted(tuple(sorted(c)) for c in colors))


def relation_witness(values: tuple[int, ...], spec: RelationSpec) -> tuple[tuple[int, ...], int] | None:
    """Return one additive witness `(summands, target)` if the cauldron boils."""
    vals = tuple(sorted(values))
    for target in vals:
        small = [v for v in vals if v < target]
        if not small:
            continue
        if spec.repeated_terms:
            if spec.max_terms is None:
                hi = target // min(small)
            else:
                hi = spec.max_terms
            iterator = combinations_with_replacement
        else:
            hi = len(small) if spec.max_terms is None else min(spec.max_terms, len(small))
            iterator = combinations
        for term_count in range(spec.min_terms, hi + 1):
            for terms in iterator(small, term_count):
                total = sum(terms)
                if total == target:
                    return tuple(terms), target
    return None


def find_safe_coloring(
    k: int,
    n_max: int,
    spec: RelationSpec,
    state_limit: int = 2_000_000,
) -> tuple[tuple[tuple[int, ...], ...] | None, int, bool]:
    """Find one legal first-boil coloring of [1,n_max], modulo cauldron labels."""
    state_count = 0
    capped = False

    @lru_cache(maxsize=None)
    def rec(n: int, colors: tuple[tuple[int, ...], ...]) -> tuple[tuple[int, ...], ...] | None:
        nonlocal state_count, capped
        state_count += 1
        if state_count > state_limit:
            capped = True
            return None
        if n > n_max:
            return colors
        for idx, color in enumerate(colors):
            candidate = tuple(sorted(color + (n,)))
            if relation_witness(candidate, spec) is not None:
                continue
            new_colors = list(colors)
            new_colors[idx] = candidate
            found = rec(n + 1, canonical(new_colors))
            if found is not None:
                return found
        return None

    answer = rec(1, tuple(() for _ in range(k)))
    return answer, state_count, capped


def max_first_boil_prefix(
    k: int,
    spec: RelationSpec,
    search_limit: int,
    state_limit: int = 2_000_000,
) -> FirstBoilResult:
    last_coloring: tuple[tuple[int, ...], ...] | None = canonical(tuple(() for _ in range(k)))
    total_states = 0
    for n in range(0, search_limit + 1):
        coloring, states, capped = find_safe_coloring(k, n, spec, state_limit=state_limit)
        total_states += states
        if capped:
            return FirstBoilResult(k, spec, n - 1, None, last_coloring, total_states, capped=True)
        if coloring is None:
            return FirstBoilResult(k, spec, n - 1, n, last_coloring, total_states)
        last_coloring = coloring
    return FirstBoilResult(k, spec, search_limit, None, last_coloring, total_states, capped=True)


def max_all_boiled(
    k: int,
    spec: RelationSpec,
    n_limit: int = 200,
    state_limit: int = 2_000_000,
) -> RemovalResult:
    """Maximize the time of the final boil in the removal variant."""
    state_count = 0
    capped = False

    @lru_cache(maxsize=None)
    def rec(n: int, active: tuple[tuple[int, ...], ...]) -> tuple[int | None, tuple[tuple, ...]]:
        nonlocal state_count, capped
        state_count += 1
        if state_count > state_limit or n > n_limit:
            capped = True
            return None, (("CAP", n_limit, state_count),)
        if not active:
            return n - 1, ()

        best_value = -1
        best_trace: tuple[tuple, ...] = ()
        for idx, color in enumerate(active):
            candidate = tuple(sorted(color + (n,)))
            witness = relation_witness(candidate, spec)
            if witness is not None:
                new_active = list(active)
                new_active.pop(idx)
                new_active = list(canonical(new_active))
                event = ("boil", n, candidate, witness)
                if not new_active:
                    value = n
                    trace = (event,)
                else:
                    value, suffix = rec(n + 1, tuple(new_active))
                    if value is None:
                        return None, suffix
                    trace = (event,) + suffix
            else:
                new_active = list(active)
                new_active[idx] = candidate
                value, suffix = rec(n + 1, canonical(new_active))
                if value is None:
                    return None, suffix
                trace = (("keep", n, candidate, None),) + suffix
            if value > best_value:
                best_value = value
                best_trace = trace
        return best_value, best_trace

    value, trace = rec(1, tuple(() for _ in range(k)))
    return RemovalResult(k, spec, value, trace, state_count, capped)


def score_lenses() -> tuple[list[str], dict[tuple[str, str], int], list[str]]:
    lenses = [
        "weak Schur first-boil",
        "classical Schur comparison",
        "removal active-state DP",
        "finite-sums subset obligations",
        "fixed-arity generalized Schur",
        "greedy residue coloring",
        "raw cauldron labels",
    ]
    criteria = {
        "exact base predicate": {
            "weak Schur first-boil": 5,
            "classical Schur comparison": 3,
            "removal active-state DP": 5,
            "finite-sums subset obligations": 4,
            "fixed-arity generalized Schur": 4,
            "greedy residue coloring": 2,
            "raw cauldron labels": 1,
        },
        "variant separation": {
            "weak Schur first-boil": 3,
            "classical Schur comparison": 2,
            "removal active-state DP": 5,
            "finite-sums subset obligations": 4,
            "fixed-arity generalized Schur": 4,
            "greedy residue coloring": 2,
            "raw cauldron labels": 1,
        },
        "computed witness retention": {
            "weak Schur first-boil": 4,
            "classical Schur comparison": 3,
            "removal active-state DP": 5,
            "finite-sums subset obligations": 4,
            "fixed-arity generalized Schur": 2,
            "greedy residue coloring": 2,
            "raw cauldron labels": 1,
        },
        "connection to repo coimage work": {
            "weak Schur first-boil": 4,
            "classical Schur comparison": 3,
            "removal active-state DP": 5,
            "finite-sums subset obligations": 5,
            "fixed-arity generalized Schur": 3,
            "greedy residue coloring": 2,
            "raw cauldron labels": 1,
        },
    }
    totals = {lens: sum(c[lens] for c in criteria.values()) for lens in lenses}
    tie_path = {
        lens: rank
        for rank, lens in enumerate(
            [
                "removal active-state DP",
                "finite-sums subset obligations",
                "weak Schur first-boil",
                "fixed-arity generalized Schur",
                "classical Schur comparison",
                "greedy residue coloring",
                "raw cauldron labels",
            ]
        )
    }
    edges: dict[tuple[str, str], int] = {}
    for a, b in combinations(lenses, 2):
        delta = totals[a] - totals[b]
        if delta == 0:
            delta = tie_path[b] - tie_path[a]
        winner, loser = (a, b) if delta > 0 else (b, a)
        edges[(winner, loser)] = 1
    return lenses, edges, sorted(lenses, key=lambda x: (-totals[x], tie_path[x]))


def tournament_fingerprints(lenses: list[str], edges: dict[tuple[str, str], int]) -> dict:
    out_scores = {v: 0 for v in lenses}
    adjacency = {v: set() for v in lenses}
    for winner, loser in edges:
        out_scores[winner] += 1
        adjacency[winner].add(loser)

    cycles3 = []
    for a, b, c in combinations(lenses, 3):
        if b in adjacency[a] and c in adjacency[b] and a in adjacency[c]:
            cycles3.append((a, b, c))
        if c in adjacency[a] and b in adjacency[c] and a in adjacency[b]:
            cycles3.append((a, c, b))

    hp_count = 0
    for order in permutations(lenses):
        if all(order[i + 1] in adjacency[order[i]] for i in range(len(order) - 1)):
            hp_count += 1

    score_hist: dict[int, int] = {}
    for score in out_scores.values():
        score_hist[score] = score_hist.get(score, 0) + 1
    return {
        "out_scores": out_scores,
        "score_hist": dict(sorted(score_hist.items())),
        "directed_3_cycles": cycles3,
        "hamiltonian_path_count": hp_count,
    }


def fmt_coloring(coloring: tuple[tuple[int, ...], ...] | None) -> str:
    if coloring is None:
        return "(none)"
    return " | ".join("{" + ",".join(map(str, c)) + "}" for c in coloring)


def print_first_boil_section() -> None:
    weak2 = RelationSpec("two-term cauldron", 2, 2, False, "literal distinct values")
    schur2 = RelationSpec("two-term classical", 2, 2, True, "allows A=A")
    up_to_3 = RelationSpec("two-or-three-term cauldron", 2, 3, False)
    finite_sums = RelationSpec("finite-sums cauldron", 2, None, False)

    jobs = [
        (weak2, 3, 30),
        (schur2, 3, 20),
        (up_to_3, 3, 30),
        (finite_sums, 3, 30),
    ]
    print("FIRST-BOIL EXACT SCOUTS")
    print("-----------------------")
    for spec, k, limit in jobs:
        result = max_first_boil_prefix(k, spec, limit)
        forced = "not reached" if result.forced_at is None else str(result.forced_at)
        cap = " [CAPPED]" if result.capped else ""
        print(f"{spec.label()}, k={k}: safe_prefix={result.safe_prefix}, forced_at={forced}, states={result.states}{cap}")
        print(f"  witness coloring at safe prefix: {fmt_coloring(result.coloring)}")
    print()

    print("SMALL k TABLES FOR TWO-TERM RULE")
    print("--------------------------------")
    for spec in [weak2, schur2]:
        row = []
        for k in [1, 2, 3]:
            limit = 30 if not spec.repeated_terms else 20
            result = max_first_boil_prefix(k, spec, limit)
            row.append(f"k={k}:safe={result.safe_prefix},force={result.forced_at}")
        print(f"{spec.label()}: " + "; ".join(row))
    print()


def print_removal_section() -> None:
    print("ALL-BOILED REMOVAL EXACT SCOUTS")
    print("-------------------------------")
    for spec in [
        RelationSpec("two-term cauldron", 2, 2, False, "literal distinct values"),
        RelationSpec("two-term classical", 2, 2, True, "allows A=A"),
        RelationSpec("finite-sums cauldron", 2, None, False),
    ]:
        for k in [1, 2, 3]:
            result = max_all_boiled(k, spec)
            cap = " [CAPPED]" if result.capped else ""
            print(f"{spec.label()}, k={k}: last_boil_at={result.last_boil_at}, states={result.states}{cap}")
            boil_events = [event for event in result.trace if event[0] == "boil"]
            for event in boil_events:
                _, n, cauldron, witness = event
                terms, target = witness
                print(f"  boil at {n}: cauldron={cauldron}, witness={terms}->{target}")
        print()


def print_variant_notes() -> None:
    print("VARIANT NOTES")
    print("-------------")
    print("- Literal cauldrons use distinct summands: weak Schur, not classical Schur.")
    print("- Repeated-summand output is included because many Schur-number tables use that convention.")
    print("- Exactly three summands for three colors is a higher generalized Schur-number problem; this scout does not brute-force its threshold.")
    print("- Finite-sums means any distinct finite subset of size at least 2 summing to another value; for k=3 it forces earlier than the weak two-term rule in the computed prefix.")
    print("- Finite-sums removal is exact here for k<=3 and gives last boils 3,10,25.")
    print()


def print_tournament_analysis() -> None:
    lenses, edges, ranking = score_lenses()
    fingerprints = tournament_fingerprints(lenses, edges)
    print("TOURNAMENT ANALYSIS")
    print("-------------------")
    print("Vertices: proof/computation lenses, not cauldrons.")
    print("Pairwise observable: which lens better preserves the cauldron survival predicate while retaining useful witnesses.")
    print("Switch/gauge: sum of exactness, variant separation, witness retention, and repo coimage connection scores.")
    print("Tie Hamiltonian path: removal DP -> finite-sums -> weak Schur -> fixed-arity -> classical comparison -> residues -> raw labels.")
    print("Ranking: " + " > ".join(ranking))
    print(f"Score histogram: {fingerprints['score_hist']}")
    print(f"Directed 3-cycles: {len(fingerprints['directed_3_cycles'])}")
    print(f"Hamiltonian-path count: {fingerprints['hamiltonian_path_count']}")
    print("Out-scores:")
    for lens, score in sorted(fingerprints["out_scores"].items(), key=lambda item: (-item[1], item[0])):
        print(f"  {lens}: {score}")
    print()


def print_assumption_challenge() -> None:
    print("ASSUMPTION CHALLENGE")
    print("--------------------")
    print("Alternate vertex sets considered: cauldrons/colors, active-state signatures, first-boil witnesses, removal order types, fixed arities, subset-sum obligations, residue classes, and proof routes.")
    print("Preserved predicate: survival of the selected additive boil relation over the current prefix.")
    print("Destroyed information: canonical active-state quotients forget cauldron labels; first-boil Schur quotients forget removal timing; finite-sums quotients forget which subset size caused failure unless witnesses are retained.")
    print("Challenged assumption: the literal game is not automatically classical Schur, because each natural is placed once and repeated summands require an extra rule.")
    print()


def main() -> None:
    print("CAULDRON GAME S618")
    print("==================")
    print_first_boil_section()
    print_removal_section()
    print_variant_notes()
    print_tournament_analysis()
    print_assumption_challenge()


if __name__ == "__main__":
    main()

