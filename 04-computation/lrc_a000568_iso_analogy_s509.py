#!/usr/bin/env python3
"""Probe the LRC -> tournament-isoclass / A000568 analogy.

codex-2026-06-01 S509

The user's hypothesis is that the Lonely Runner Conjecture is ultimately
analogous to a problem on tournament isomorphism classes and A000568.  This
script makes that precise as a projection experiment.

For a speed set with the stationary observer included, time traces a movie of
phase configurations on the circle.  Between endpoint and half-turn walls, the
following data are constant enough to sample at the cell midpoint:

    LRC safe bit          are all moving runners at distance >= 1/N from 0?
    phase tournament      orient i -> j if j is clockwise from i by < 1/2
    unmarked iso class    quotient by all vertex relabelings (A000568 vertex)
    pointed iso class     quotient by relabelings that keep observer 0 fixed

Tournament Analysis declaration:

    pairwise observable: half-turn circular phase difference
    switch/gauge: delta in (0, 1/2), ties by increasing vertex label
    tie Hamiltonian path: 0 -> 1 -> ... -> N-1

The key defect measured here is whether the safe/unsafe LRC bit is a function
of the unmarked or pointed tournament isomorphism class.  If not, plain
A000568 is not the final state space; it is the base of a richer pointed,
threshold-labelled sheaf.

Stored output:
    05-knowledge/results/lrc_a000568_iso_analogy_s509.out
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations, permutations
from math import factorial, gcd


ONE = Fraction(1, 1)
HALF = Fraction(1, 2)

A000568 = {
    1: 1,
    2: 1,
    3: 2,
    4: 4,
    5: 12,
    6: 56,
    7: 456,
    8: 6880,
    9: 191536,
    10: 9733056,
}


@dataclass(frozen=True)
class Family:
    name: str
    speeds: tuple[int, ...]  # moving speeds; observer 0 is added by the script


@dataclass(frozen=True)
class CellRecord:
    family: str
    t: Fraction
    safe: bool
    unmarked: tuple[int, ...]
    pointed: tuple[int, ...]
    hamiltonian_paths: int
    score_sequence: tuple[int, ...]
    c3: int
    scc_count: int
    stationary_score: int


def circle(x: Fraction) -> Fraction:
    return x % ONE


def dist0(x: Fraction) -> Fraction:
    x = circle(x)
    return min(x, ONE - x)


def positions(speeds0: tuple[int, ...], t: Fraction) -> tuple[Fraction, ...]:
    return tuple(circle(s * t) for s in speeds0)


def phase_tournament(pos: tuple[Fraction, ...]) -> list[list[int]]:
    n = len(pos)
    adj = [[0] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        delta = circle(pos[j] - pos[i])
        if delta == 0 or delta == HALF:
            winner = i
        elif delta < HALF:
            winner = i
        else:
            winner = j
        loser = j if winner == i else i
        adj[winner][loser] = 1
    return adj


def canonical(adj: list[list[int]], fixed_zero: bool = False) -> tuple[int, ...]:
    n = len(adj)
    if fixed_zero:
        perms = ((0,) + tail for tail in permutations(range(1, n)))
    else:
        perms = permutations(range(n))
    best: tuple[int, ...] | None = None
    for perm in perms:
        form = tuple(adj[perm[i]][perm[j]] for i in range(n) for j in range(n))
        if best is None or form < best:
            best = form
    assert best is not None
    return best


def count_hp(adj: list[list[int]]) -> int:
    n = len(adj)
    dp: dict[tuple[int, int], int] = defaultdict(int)
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            ways = dp.get((mask, v), 0)
            if not ways:
                continue
            for w in range(n):
                if not (mask & (1 << w)) and adj[v][w]:
                    dp[(mask | (1 << w), w)] += ways
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))


def directed_triangles(adj: list[list[int]]) -> int:
    out = 0
    for a, b, c in combinations(range(len(adj)), 3):
        if (adj[a][b] and adj[b][c] and adj[c][a]) or (
            adj[a][c] and adj[c][b] and adj[b][a]
        ):
            out += 1
    return out


def scc_count(adj: list[list[int]]) -> int:
    n = len(adj)
    graph = [[j for j in range(n) if adj[i][j]] for i in range(n)]
    rev = [[i for i in range(n) if adj[i][j]] for j in range(n)]
    seen = [False] * n
    order: list[int] = []

    def dfs(v: int) -> None:
        seen[v] = True
        for w in graph[v]:
            if not seen[w]:
                dfs(w)
        order.append(v)

    for v in range(n):
        if not seen[v]:
            dfs(v)

    seen = [False] * n
    count = 0
    for start in reversed(order):
        if seen[start]:
            continue
        count += 1
        stack = [start]
        seen[start] = True
        while stack:
            v = stack.pop()
            for w in rev[v]:
                if not seen[w]:
                    seen[w] = True
                    stack.append(w)
    return count


def wall_events(speeds: tuple[int, ...]) -> tuple[Fraction, ...]:
    """Endpoint and half-turn walls in [0,1]."""
    n_total = len(speeds) + 1
    speeds0 = (0,) + speeds
    events: set[Fraction] = {Fraction(0), Fraction(1)}

    # Half-turn walls for all pairs, including stationary observer.
    for a, b in combinations(speeds0, 2):
        d = abs(a - b)
        if d == 0:
            continue
        for m in range(d):
            events.add(Fraction(2 * m + 1, 2 * d))

    # LRC endpoint walls ||s_i t|| = 1/N.
    for speed in speeds:
        for m in range(speed + 1):
            left = Fraction(m, speed) - Fraction(1, n_total * speed)
            right = Fraction(m, speed) + Fraction(1, n_total * speed)
            if 0 <= left <= 1:
                events.add(left)
            if 0 <= right <= 1:
                events.add(right)

    return tuple(sorted(events))


def primitive(speeds: tuple[int, ...]) -> bool:
    g = 0
    for speed in speeds:
        g = gcd(g, speed)
    return g == 1


def cells_for_family(family: Family) -> list[CellRecord]:
    speeds = family.speeds
    assert primitive(speeds)
    n_total = len(speeds) + 1
    threshold = Fraction(1, n_total)
    speeds0 = (0,) + speeds
    records: list[CellRecord] = []
    events = wall_events(speeds)
    for left, right in zip(events, events[1:]):
        if right <= left:
            continue
        t = (left + right) / 2
        pos = positions(speeds0, t)
        safe = all(dist0(p) >= threshold for p in pos[1:])
        adj = phase_tournament(pos)
        scores = tuple(sorted(sum(row) for row in adj))
        records.append(
            CellRecord(
                family=family.name,
                t=t,
                safe=safe,
                unmarked=canonical(adj, fixed_zero=False),
                pointed=canonical(adj, fixed_zero=True),
                hamiltonian_paths=count_hp(adj),
                score_sequence=scores,
                c3=directed_triangles(adj),
                scc_count=scc_count(adj),
                stationary_score=sum(adj[0]),
            )
        )
    return records


def records_at_events(family: Family) -> list[CellRecord]:
    """Sample exact event times, where tight LRC witnesses can live."""
    speeds = family.speeds
    n_total = len(speeds) + 1
    threshold = Fraction(1, n_total)
    speeds0 = (0,) + speeds
    records: list[CellRecord] = []
    for t in wall_events(speeds):
        if t == 0 or t == 1:
            continue
        pos = positions(speeds0, t)
        safe = all(dist0(p) >= threshold for p in pos[1:])
        adj = phase_tournament(pos)
        scores = tuple(sorted(sum(row) for row in adj))
        records.append(
            CellRecord(
                family=family.name,
                t=t,
                safe=safe,
                unmarked=canonical(adj, fixed_zero=False),
                pointed=canonical(adj, fixed_zero=True),
                hamiltonian_paths=count_hp(adj),
                score_sequence=scores,
                c3=directed_triangles(adj),
                scc_count=scc_count(adj),
                stationary_score=sum(adj[0]),
            )
        )
    return records


def fmt_frac(x: Fraction) -> str:
    if x.denominator == 1:
        return str(x.numerator)
    return f"{x.numerator}/{x.denominator}"


def mixed_classes(records: list[CellRecord], attr: str) -> tuple[int, int, int]:
    buckets: dict[tuple[int, ...], Counter] = defaultdict(Counter)
    for rec in records:
        buckets[getattr(rec, attr)][rec.safe] += 1
    safe_only = sum(1 for c in buckets.values() if c[True] and not c[False])
    unsafe_only = sum(1 for c in buckets.values() if c[False] and not c[True])
    mixed = sum(1 for c in buckets.values() if c[True] and c[False])
    return safe_only, unsafe_only, mixed


def transition_stats(records: list[CellRecord], attr: str) -> tuple[int, int, int]:
    transitions = []
    for a, b in zip(records, records[1:]):
        if getattr(a, attr) != getattr(b, attr):
            transitions.append((getattr(a, attr), getattr(b, attr)))
    unique_undirected = {frozenset(pair) for pair in transitions}
    loops = sum(1 for a, b in zip(records, records[1:]) if getattr(a, attr) == getattr(b, attr))
    return len(transitions), len(unique_undirected), loops


def summarize_family(
    family: Family,
    records: list[CellRecord],
    boundary_records: list[CellRecord],
) -> None:
    n_total = len(family.speeds) + 1
    safe_cells = sum(1 for r in records if r.safe)
    unmarked = {r.unmarked for r in records}
    pointed = {r.pointed for r in records}
    safe_unmarked = {r.unmarked for r in records if r.safe}
    safe_pointed = {r.pointed for r in records if r.safe}
    su, uu, mu = mixed_classes(records, "unmarked")
    sp, up, mp = mixed_classes(records, "pointed")
    tr_u, edge_u, loop_u = transition_stats(records, "unmarked")
    tr_p, edge_p, loop_p = transition_stats(records, "pointed")
    hp_safe = Counter(r.hamiltonian_paths for r in records if r.safe)
    hp_unsafe = Counter(r.hamiltonian_paths for r in records if not r.safe)
    stat_safe = Counter(r.stationary_score for r in records if r.safe)
    stat_unsafe = Counter(r.stationary_score for r in records if not r.safe)
    safe_boundaries = [r for r in boundary_records if r.safe]
    boundary_unmarked = {r.unmarked for r in safe_boundaries}
    boundary_pointed = {r.pointed for r in safe_boundaries}

    print(f"\nFAMILY {family.name}")
    print(f"  speeds={family.speeds} total_N={n_total} A000568={A000568.get(n_total)}")
    print(f"  cells={len(records)} safe_cells={safe_cells} unsafe_cells={len(records)-safe_cells}")
    print(
        "  safe boundary events=%d boundary unmarked/pointed classes=%d/%d"
        % (len(safe_boundaries), len(boundary_unmarked), len(boundary_pointed))
    )
    print(
        "  visited unmarked iso classes=%d (%.3f of A000568); pointed=%d"
        % (
            len(unmarked),
            len(unmarked) / A000568[n_total],
            len(pointed),
        )
    )
    print(
        "  safe-hit classes: unmarked=%d pointed=%d"
        % (len(safe_unmarked), len(safe_pointed))
    )
    print(
        "  projection defect unmarked safe-only/unsafe-only/mixed = %d/%d/%d"
        % (su, uu, mu)
    )
    print(
        "  projection defect pointed   safe-only/unsafe-only/mixed = %d/%d/%d"
        % (sp, up, mp)
    )
    print(
        "  quotient walk unmarked: transitions=%d unique_edges=%d self_steps=%d"
        % (tr_u, edge_u, loop_u)
    )
    print(
        "  quotient walk pointed:   transitions=%d unique_edges=%d self_steps=%d"
        % (tr_p, edge_p, loop_p)
    )
    print(f"  H among safe cells:   {sorted(hp_safe.items())[:8]}")
    print(f"  H among unsafe cells: {sorted(hp_unsafe.items())[:8]}")
    print(f"  observer score safe:   {sorted(stat_safe.items())}")
    print(f"  observer score unsafe: {sorted(stat_unsafe.items())}")

    # Give one concrete mixed witness if present.
    by_class: dict[tuple[int, ...], list[CellRecord]] = defaultdict(list)
    for rec in records:
        by_class[rec.unmarked].append(rec)
    for cls, items in by_class.items():
        if any(r.safe for r in items) and any(not r.safe for r in items):
            safe = next(r for r in items if r.safe)
            unsafe = next(r for r in items if not r.safe)
            print(
                "  example unmarked-mixed class:"
                f" safe_t={fmt_frac(safe.t)} unsafe_t={fmt_frac(unsafe.t)}"
                f" H={safe.hamiltonian_paths}/{unsafe.hamiltonian_paths}"
                f" score={safe.score_sequence}/{unsafe.score_sequence}"
                f" stationary_score={safe.stationary_score}/{unsafe.stationary_score}"
            )
            break


@dataclass(frozen=True)
class Route:
    name: str
    closeness_to_a000568: int
    keeps_lrc_threshold: int
    computational_handle: int
    proof_potential: int
    novelty: int
    projection_risk: int

    @property
    def score(self) -> int:
        return (
            2 * self.closeness_to_a000568
            + 3 * self.keeps_lrc_threshold
            + 2 * self.computational_handle
            + 3 * self.proof_potential
            + self.novelty
            - 2 * self.projection_risk
        )


def route_tournament(routes: list[Route]) -> list[list[int]]:
    n = len(routes)
    adj = [[0] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        votes_i = 0
        votes_j = 0
        for attr in (
            "closeness_to_a000568",
            "keeps_lrc_threshold",
            "computational_handle",
            "proof_potential",
            "novelty",
        ):
            ai = getattr(routes[i], attr)
            aj = getattr(routes[j], attr)
            if ai > aj:
                votes_i += 1
            elif aj > ai:
                votes_j += 1
        if routes[i].projection_risk < routes[j].projection_risk:
            votes_i += 1
        elif routes[j].projection_risk < routes[i].projection_risk:
            votes_j += 1

        if votes_i == votes_j:
            winner = i
        else:
            winner = i if votes_i > votes_j else j
        loser = j if winner == i else i
        adj[winner][loser] = 1
    return adj


def route_analysis() -> None:
    routes = [
        Route("clock-word in G_N", 5, 2, 5, 3, 4, 4),
        Route("pointed A000568 lift", 4, 3, 4, 4, 4, 3),
        Route("alpha=1/N gauge bundle", 3, 5, 4, 5, 5, 2),
        Route("endpoint-pressure sheaf", 2, 5, 3, 5, 5, 2),
        Route("Burnside fiber mass", 5, 1, 5, 3, 4, 5),
        Route("G_N wall-crossing corridors", 5, 3, 4, 4, 5, 3),
        Route("A000568 entropy barrier", 5, 1, 5, 2, 5, 5),
        Route("even-graph cycle-space dual", 3, 3, 3, 4, 5, 3),
        Route("Hamiltonian-path H meter", 4, 2, 5, 3, 3, 4),
        Route("projection-defect obstruction", 4, 5, 4, 5, 5, 1),
        Route("residue bucket transport", 3, 5, 3, 5, 4, 2),
        Route("staircase extremal geodesic", 5, 3, 4, 4, 5, 3),
    ]
    adj = route_tournament(routes)
    scores = [sum(row) for row in adj]
    print("\nROUTE TOURNAMENT")
    print("  observable = six route attributes")
    print("  switch = majority of attributes wins; lower projection risk wins risk vote")
    print("  tie path = listed order")
    print(f"  H={count_hp(adj)}")
    print(f"  score_hist={dict(sorted(Counter(scores).items()))}")
    print(f"  c3={directed_triangles(adj)} scc_count={scc_count(adj)}")
    ranked = sorted(zip(routes, scores), key=lambda rs: (-rs[0].score, -rs[1], rs[0].name))
    print("  ranking:")
    for route, outdeg in ranked:
        print(
            "    %-31s route_score=%2d tournament_score=%2d"
            % (route.name, route.score, outdeg)
        )


def main() -> None:
    print("LRC -> A000568 / tournament-isoclass analogy probe (S509)")
    print("=" * 72)
    print("Interpretation: if safe/unsafe is mixed inside an iso-class bucket,")
    print("plain A000568 is the base quotient, not the whole LRC state.")

    families = [
        Family("N4 consecutive", (1, 2, 3)),
        Family("N4 sparse", (1, 3, 4)),
        Family("N5 consecutive", (1, 2, 3, 4)),
        Family("N5 powers-ish", (1, 2, 4, 8)),
        Family("N5 prime-ish", (2, 3, 5, 7)),
        Family("N6 consecutive", (1, 2, 3, 4, 5)),
        Family("N6 mixed", (1, 3, 4, 7, 9)),
        Family("N7 consecutive", (1, 2, 3, 4, 5, 6)),
        Family("N7 mixed", (1, 4, 6, 9, 10, 15)),
    ]

    all_records: list[CellRecord] = []
    for family in families:
        records = cells_for_family(family)
        boundary_records = records_at_events(family)
        summarize_family(family, records, boundary_records)
        all_records.extend(records)

    print("\nGLOBAL PROJECTION SUMMARY")
    for n_total in sorted({len(f.speeds) + 1 for f in families}):
        recs = [r for r in all_records if len(next(f for f in families if f.name == r.family).speeds) + 1 == n_total]
        su, uu, mu = mixed_classes(recs, "unmarked")
        sp, up, mp = mixed_classes(recs, "pointed")
        print(
            "  N=%d A000568=%d cells=%d unmarked_mixed=%d pointed_mixed=%d"
            " safe_only_unmarked=%d unsafe_only_unmarked=%d"
            % (n_total, A000568[n_total], len(recs), mu, mp, su, uu)
        )

    route_analysis()

    print("\nSYNTHESIS")
    print("  1. The LRC movie naturally walks on tournament isomorphism classes.")
    print("  2. The safe bit is often not a function of the unmarked A000568 class.")
    print("  3. Even pointed classes can mix safe and unsafe cells; gap lengths and")
    print("     endpoint labels are extra sheaf data over the A000568 base.")
    print("  4. Therefore the most plausible analogy is convoluted but sharp:")
    print("     LRC is a forbidden-cover / section-existence problem over the")
    print("     A000568 isoclass quotient, with projection defect as the obstruction.")


if __name__ == "__main__":
    main()
