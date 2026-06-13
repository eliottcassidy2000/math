#!/usr/bin/env python3
"""Observer-marked predicate zoo inspired by HYP-1981.

codex-2026-06-01 S517

HYP-1981 and THM-381 identify the exact LRC witness:

    observer lonely  <=>  marked observer is a source.

This script asks what the nearby marked-tournament predicates represent.  The
main point is that the observer score stratifies the whole movie by exact
LRC blocker count:

    outdeg(observer) = number of safe runners
    indeg(observer)  = number of near/blocking runners

Thus "almost source" means exactly one runner blocks the observer, and the
observer's directed distance to other vertices gives a possible repair/pressure
language: a 2-king observer means every blocker is beaten by some currently
safe runner in the half-turn phase tournament.

Tournament Analysis declaration:

Predicate tournament
    vertices: source, almost-source, observer score layer, side-defect layer,
        observer 2-king, runner-subclass menu, and endpoint-pressure core.
    pairwise observable: exactness, locality, boundary fidelity, quotient
        compression, non-tautological content, pressure compatibility, and
        proof leverage.
    switch/gauge: a predicate points to another when it wins a majority of
        those dimensions.
    tie Hamiltonian path: the listed predicate order.

Stored output:
    05-knowledge/results/lrc_observer_predicate_zoo_s517.out
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, permutations
from math import gcd
from typing import Callable


ONE = Fraction(1)
INF = 10**9


def frac(value: Fraction) -> Fraction:
    return value - (value.numerator // value.denominator)


def dist0(value: Fraction) -> Fraction:
    x = frac(value)
    return min(x, ONE - x)


def fmt_frac(value: Fraction | int | None) -> str:
    if value is None:
        return "-"
    if isinstance(value, int):
        return str(value)
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def primitive_speed_sets(total_n: int, max_speed: int) -> list[tuple[int, ...]]:
    out: list[tuple[int, ...]] = []
    for moving in combinations(range(1, max_speed + 1), total_n - 1):
        g = 0
        for speed in moving:
            g = gcd(g, speed)
        if g == 1:
            out.append((0,) + moving)
    return out


def positions(speeds: tuple[int, ...], time: Fraction) -> tuple[Fraction, ...]:
    return tuple(frac(Fraction(speed) * time) for speed in speeds)


def near_count(speeds: tuple[int, ...], time: Fraction) -> int:
    threshold = Fraction(1, len(speeds))
    return sum(1 for speed in speeds[1:] if dist0(Fraction(speed) * time) < threshold)


def side_defect_code(speeds: tuple[int, ...], time: Fraction) -> tuple[int, int, int]:
    """Return (left blockers, right blockers, observer ties)."""
    threshold = Fraction(1, len(speeds))
    left = right = tied = 0
    for pos in positions(speeds, time)[1:]:
        if pos == 0:
            tied += 1
        elif 0 < pos < threshold:
            right += 1
        elif ONE - threshold < pos < ONE:
            left += 1
    return left, right, tied


def observer_walls(speeds: tuple[int, ...]) -> set[Fraction]:
    total_n = len(speeds)
    out: set[Fraction] = {Fraction(0)}
    for speed in speeds[1:]:
        for m in range(speed):
            out.add(frac(Fraction(m * total_n + 1, total_n * speed)))
            out.add(frac(Fraction(m * total_n + total_n - 1, total_n * speed)))
    return out


def half_turn_walls(speeds: tuple[int, ...]) -> set[Fraction]:
    out: set[Fraction] = set()
    for i, j in combinations(range(1, len(speeds)), 2):
        d = abs(speeds[i] - speeds[j])
        if d == 0:
            continue
        for k in range(2 * d):
            out.add(frac(Fraction(k, 2 * d)))
    return out


def event_walls(speeds: tuple[int, ...]) -> list[Fraction]:
    return sorted(observer_walls(speeds) | half_turn_walls(speeds))


def sampled_states(speeds: tuple[int, ...]) -> list[tuple[Fraction, str]]:
    walls = event_walls(speeds)
    out: list[tuple[Fraction, str]] = [(time, "wall") for time in walls]
    cyclic = walls + [ONE]
    for left, right in zip(cyclic, cyclic[1:]):
        if left < right:
            out.append(((left + right) / 2, "open"))
    return sorted(out, key=lambda item: (item[0], item[1]))


def marked_tournament(speeds: tuple[int, ...], time: Fraction) -> tuple[tuple[int, ...], ...]:
    total_n = len(speeds)
    threshold = Fraction(1, total_n)
    adj = [[0] * total_n for _ in range(total_n)]
    for i, j in combinations(range(total_n), 2):
        if i == 0 or j == 0:
            runner = j if i == 0 else i
            observer_wins = dist0(Fraction(speeds[runner]) * time) >= threshold
            winner = 0 if observer_wins else runner
        else:
            gap = frac(Fraction(speeds[i] - speeds[j]) * time)
            if 0 < gap < Fraction(1, 2):
                winner = i
            elif Fraction(1, 2) < gap < ONE:
                winner = j
            else:
                # Fixed Hamiltonian tie path on runner labels.
                winner = min(i, j)
        loser = j if winner == i else i
        adj[winner][loser] = 1
    return tuple(tuple(row) for row in adj)


@lru_cache(maxsize=None)
def canonical_marked(adj: tuple[tuple[int, ...], ...]) -> tuple[int, ...]:
    total_n = len(adj)
    best: tuple[int, ...] | None = None
    for perm_tail in permutations(range(1, total_n)):
        perm = (0,) + perm_tail
        bits = tuple(adj[perm[i]][perm[j]] for i in range(total_n) for j in range(i + 1, total_n))
        if best is None or bits < best:
            best = bits
    assert best is not None
    return best


@lru_cache(maxsize=None)
def canonical_runners(adj: tuple[tuple[int, ...], ...]) -> tuple[int, ...]:
    runner_count = len(adj) - 1
    best: tuple[int, ...] | None = None
    for perm in permutations(range(1, len(adj))):
        bits = tuple(adj[perm[i]][perm[j]] for i in range(runner_count) for j in range(i + 1, runner_count))
        if best is None or bits < best:
            best = bits
    assert best is not None
    return best


def runner_adj(adj: tuple[tuple[int, ...], ...]) -> tuple[tuple[int, ...], ...]:
    return tuple(tuple(adj[i][j] for j in range(1, len(adj))) for i in range(1, len(adj)))


def score_sequence(adj: tuple[tuple[int, ...], ...]) -> tuple[int, ...]:
    return tuple(sorted(sum(row) for row in adj))


def hamiltonian_paths(adj: tuple[tuple[int, ...], ...]) -> int:
    n = len(adj)
    full = (1 << n) - 1

    @lru_cache(maxsize=None)
    def dp(mask: int, last: int) -> int:
        if mask == full:
            return 1
        return sum(
            dp(mask | (1 << nxt), nxt)
            for nxt in range(n)
            if not ((mask >> nxt) & 1) and adj[last][nxt]
        )

    return sum(dp(1 << start, start) for start in range(n))


def observer_radius(adj: tuple[tuple[int, ...], ...]) -> int | None:
    n = len(adj)
    dist = [INF] * n
    dist[0] = 0
    queue = deque([0])
    while queue:
        vertex = queue.popleft()
        for nxt, bit in enumerate(adj[vertex]):
            if bit and dist[nxt] == INF:
                dist[nxt] = dist[vertex] + 1
                queue.append(nxt)
    if any(value == INF for value in dist):
        return None
    return max(dist)


def observer_outdegree(adj: tuple[tuple[int, ...], ...]) -> int:
    return sum(adj[0])


@dataclass
class Audit:
    total_n: int
    max_speed: int
    systems: int = 0
    states: int = 0
    score_mismatches: int = 0
    marked_mixed: int = 0
    runner_mixed: int = 0
    system_profile: Counter[str] | None = None
    min_open_hist: Counter[int] | None = None
    min_all_hist: Counter[int] | None = None
    layer_states: Counter[int] | None = None
    layer_marked_classes: dict[int, set[tuple[int, ...]]] | None = None
    layer_runner_classes: dict[int, set[tuple[int, ...]]] | None = None
    layer_radius: Counter[tuple[int, str]] | None = None
    side_codes: Counter[tuple[int, tuple[int, int, int]]] | None = None
    source_runner_classes: dict[tuple[int, ...], tuple[int, tuple[int, ...]]] | None = None


def audit_window(total_n: int, max_speed: int) -> Audit:
    audit = Audit(
        total_n=total_n,
        max_speed=max_speed,
        system_profile=Counter(),
        min_open_hist=Counter(),
        min_all_hist=Counter(),
        layer_states=Counter(),
        layer_marked_classes=defaultdict(set),
        layer_runner_classes=defaultdict(set),
        layer_radius=Counter(),
        side_codes=Counter(),
        source_runner_classes={},
    )
    class_to_near: dict[tuple[int, ...], set[int]] = defaultdict(set)
    runner_to_near: dict[tuple[int, ...], set[int]] = defaultdict(set)
    for speeds in primitive_speed_sets(total_n, max_speed):
        audit.systems += 1
        min_open: int | None = None
        min_all: int | None = None
        open_source = False
        wall_source = False
        for time, kind in sampled_states(speeds):
            audit.states += 1
            adj = marked_tournament(speeds, time)
            blockers = near_count(speeds, time)
            outdeg = observer_outdegree(adj)
            if outdeg != total_n - 1 - blockers:
                audit.score_mismatches += 1
            if kind == "open":
                min_open = blockers if min_open is None else min(min_open, blockers)
                open_source = open_source or blockers == 0
            else:
                wall_source = wall_source or blockers == 0
            min_all = blockers if min_all is None else min(min_all, blockers)

            marked = canonical_marked(adj)
            runners = canonical_runners(adj)
            class_to_near[marked].add(blockers)
            runner_to_near[runners].add(blockers)
            audit.layer_states[blockers] += 1
            audit.layer_marked_classes[blockers].add(marked)
            audit.layer_runner_classes[blockers].add(runners)
            radius = observer_radius(adj)
            audit.layer_radius[(blockers, "unreachable" if radius is None else str(radius))] += 1
            audit.side_codes[(blockers, side_defect_code(speeds, time))] += 1
            if blockers == 0:
                r_adj = runner_adj(adj)
                audit.source_runner_classes.setdefault(
                    runners,
                    (hamiltonian_paths(r_adj), score_sequence(r_adj)),
                )
        assert min_all is not None
        audit.min_all_hist[min_all] += 1
        if min_open is not None:
            audit.min_open_hist[min_open] += 1
        if open_source:
            audit.system_profile["open-source"] += 1
        elif wall_source:
            audit.system_profile["wall-only-source"] += 1
        else:
            audit.system_profile["source-avoiding"] += 1
    audit.marked_mixed = sum(1 for values in class_to_near.values() if len(values) > 1)
    audit.runner_mixed = sum(1 for values in runner_to_near.values() if len(values) > 1)
    return audit


def compact_counter(counter: Counter[object], limit: int = 6) -> str:
    items = counter.most_common(limit)
    return ", ".join(f"{key}:{value}" for key, value in items)


def print_audit(audit: Audit) -> None:
    assert audit.system_profile is not None
    assert audit.min_open_hist is not None
    assert audit.min_all_hist is not None
    assert audit.layer_states is not None
    assert audit.layer_marked_classes is not None
    assert audit.layer_runner_classes is not None
    assert audit.layer_radius is not None
    assert audit.side_codes is not None
    assert audit.source_runner_classes is not None
    print(
        f"N={audit.total_n} max_speed={audit.max_speed} systems={audit.systems} "
        f"states={audit.states} score_mismatches={audit.score_mismatches}"
    )
    print(
        f"  marked mixed near-count classes={audit.marked_mixed}; "
        f"runner-subclass mixed near-count classes={audit.runner_mixed}"
    )
    print(f"  system profiles: {dict(sorted(audit.system_profile.items()))}")
    print(f"  min blockers over all states: {dict(sorted(audit.min_all_hist.items()))}")
    print(f"  min blockers over open cells: {dict(sorted(audit.min_open_hist.items()))}")
    print("  observer-score layers (blockers = indeg(observer)):")
    for blockers in sorted(audit.layer_states):
        radius_hist = Counter(
            radius for (layer, radius), count in audit.layer_radius.items() if layer == blockers for _ in range(count)
        )
        side_hist = Counter(
            code for (layer, code), count in audit.side_codes.items() if layer == blockers for _ in range(count)
        )
        print(
            f"    blockers={blockers:<2} states={audit.layer_states[blockers]:<6} "
            f"marked_classes={len(audit.layer_marked_classes[blockers]):<5} "
            f"runner_classes={len(audit.layer_runner_classes[blockers]):<5} "
            f"radius={compact_counter(radius_hist, 4)} "
            f"sides={compact_counter(side_hist, 4)}"
        )
    source_details = sorted(set(audit.source_runner_classes.values()))
    print(
        f"  source runner menu: {len(audit.source_runner_classes)} classes; "
        f"H/score list={source_details[:10]}"
    )
    print()


@dataclass(frozen=True)
class Predicate:
    name: str
    exactness: int
    locality: int
    boundary: int
    compression: int
    non_taut: int
    pressure: int
    leverage: int


PREDICATES = (
    Predicate("source", 5, 5, 5, 5, 2, 2, 5),
    Predicate("almost_source", 5, 5, 5, 4, 3, 3, 4),
    Predicate("observer_score_layer", 5, 5, 5, 4, 2, 3, 4),
    Predicate("side_defect_layer", 5, 5, 5, 3, 3, 4, 4),
    Predicate("observer_2_king", 4, 4, 3, 4, 5, 4, 3),
    Predicate("runner_subclass_menu", 3, 2, 3, 5, 5, 2, 3),
    Predicate("endpoint_pressure_core", 5, 3, 4, 3, 5, 5, 5),
)


DIMENSIONS: tuple[tuple[str, Callable[[Predicate], int]], ...] = (
    ("exactness", lambda pred: pred.exactness),
    ("locality", lambda pred: pred.locality),
    ("boundary", lambda pred: pred.boundary),
    ("compression", lambda pred: pred.compression),
    ("non_taut", lambda pred: pred.non_taut),
    ("pressure", lambda pred: pred.pressure),
    ("leverage", lambda pred: pred.leverage),
)


def predicate_tournament() -> list[list[int]]:
    adj = [[0] * len(PREDICATES) for _ in PREDICATES]
    for i, j in combinations(range(len(PREDICATES)), 2):
        wins_i = 0
        wins_j = 0
        for _, getter in DIMENSIONS:
            vi = getter(PREDICATES[i])
            vj = getter(PREDICATES[j])
            if vi > vj:
                wins_i += 1
            elif vj > vi:
                wins_j += 1
        winner = i if wins_i >= wins_j else j
        loser = j if winner == i else i
        adj[winner][loser] = 1
    return adj


def h_paths_matrix(adj: list[list[int]]) -> int:
    n = len(adj)
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for start in range(n):
        dp[1 << start][start] = 1
    for mask in range(1 << n):
        for last in range(n):
            val = dp[mask][last]
            if not val:
                continue
            for nxt in range(n):
                if not ((mask >> nxt) & 1) and adj[last][nxt]:
                    dp[mask | (1 << nxt)][nxt] += val
    return sum(dp[full])


def triangle_count(adj: list[list[int]]) -> int:
    total = 0
    for a, b, c in combinations(range(len(adj)), 3):
        if (adj[a][b] and adj[b][c] and adj[c][a]) or (
            adj[a][c] and adj[c][b] and adj[b][a]
        ):
            total += 1
    return total


def scc_sizes(adj: list[list[int]]) -> tuple[int, ...]:
    n = len(adj)
    seen = [False] * n
    order: list[int] = []

    def dfs(vertex: int) -> None:
        seen[vertex] = True
        for nxt, bit in enumerate(adj[vertex]):
            if bit and not seen[nxt]:
                dfs(nxt)
        order.append(vertex)

    for vertex in range(n):
        if not seen[vertex]:
            dfs(vertex)
    rev = [[adj[j][i] for j in range(n)] for i in range(n)]
    seen = [False] * n
    sizes: list[int] = []
    for start in reversed(order):
        if seen[start]:
            continue
        queue = deque([start])
        seen[start] = True
        size = 0
        while queue:
            vertex = queue.popleft()
            size += 1
            for nxt, bit in enumerate(rev[vertex]):
                if bit and not seen[nxt]:
                    seen[nxt] = True
                    queue.append(nxt)
        sizes.append(size)
    return tuple(sorted(sizes, reverse=True))


def print_predicate_tournament() -> None:
    print("Predicate Tournament Analysis")
    print("-" * 88)
    adj = predicate_tournament()
    scores = [sum(row) for row in adj]
    print(
        f"predicates={len(PREDICATES)} H={h_paths_matrix(adj)} "
        f"c3={triangle_count(adj)} SCCs={scc_sizes(adj)} "
        f"score_hist={dict(sorted(Counter(scores).items()))}"
    )
    ranked = sorted(((scores[idx], -idx, PREDICATES[idx]) for idx in range(len(PREDICATES))), reverse=True)
    for score, _, pred in ranked:
        print(
            f"  {score:>2} {pred.name:<24} exact={pred.exactness} local={pred.locality} "
            f"boundary={pred.boundary} compress={pred.compression} non_taut={pred.non_taut} "
            f"pressure={pred.pressure} leverage={pred.leverage}"
        )
    print()


def main() -> None:
    print("LRC observer-marked predicate zoo (S517)")
    print("=" * 88)
    print("Part A. What source-like predicates represent")
    print("-" * 88)
    print("source:        blockers=0; exact LRC witness (THM-381)")
    print("almost-source: blockers=1; one observer incident edge flip from source")
    print("score layer:   blockers=k; distance-to-source in observer incident edges")
    print("side layer:    blockers split by left/right forbidden caps and observer ties")
    print("2-king:        every blocker is reachable via a safe runner in two steps")
    print()
    print("Part B. Exact bounded audits with tie-completed runner clocks")
    print("-" * 88)
    for total_n, max_speed in ((4, 12), (5, 10), (6, 9), (7, 8)):
        print_audit(audit_window(total_n, max_speed))
    print_predicate_tournament()
    print("SYNTHESIS")
    print("=" * 88)
    print(
        "1. HYP-1981's source predicate has a whole exact stratification around it: "
        "observer indegree is precisely the number of near/blocking runners."
    )
    print(
        "2. Marked classes never mix blocker counts, because marked isomorphism "
        "preserves observer score.  Runner-subtournament classes do mix them, "
        "so forgetting the observer incident threshold edges loses the LRC layer."
    )
    print(
        "3. Almost-source and side-defect layers are the natural next targets for "
        "tight and near-tight systems: they record one-blocker or two-cap debt "
        "instead of only the zero-debt source event."
    )
    print(
        "4. Observer 2-king status is not equivalent to loneliness.  It is a "
        "repair/pressure predicate: every current blocker is beaten by some "
        "safe runner, suggesting a bridge to labelled endpoint-pressure cores."
    )


if __name__ == "__main__":
    main()
