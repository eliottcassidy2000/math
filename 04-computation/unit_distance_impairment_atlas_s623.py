#!/usr/bin/env python3
"""
unit_distance_impairment_atlas_s623.py

S623: small-size impairment atlas for unit-distance construction/proof methods.

The point is to deliberately damage a method and measure the damage.  If a
small dense construction collapses when a side-channel is removed, that
side-channel is a candidate proof invariant rather than implementation detail.

Inspirations imported from other repo threads:

  - LRC: quotient compression works only with retained side channels.
  - Tournaments: scalar counts hide forbidden packet obstructions.
  - Cauldrons: adversarial schedules expose the resource a one-player search
    forgot.
  - Unit-distance S617: frontier gains are the state-local observable.

This script keeps the runs small and reproducible.  It is not a proof of any
unit-distance optimum.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from itertools import combinations, permutations, product
import sys


Point2 = tuple[int, int]
Point4 = tuple[int, int, int, int]


TRI_UNITS: tuple[Point2, ...] = (
    (1, 0),
    (0, 1),
    (1, -1),
    (-1, 0),
    (0, -1),
    (-1, 1),
)


def add2(a: Point2, b: Point2) -> Point2:
    return (a[0] + b[0], a[1] + b[1])


def sub2(a: Point2, b: Point2) -> Point2:
    return (a[0] - b[0], a[1] - b[1])


def canon2(points: tuple[Point2, ...]) -> tuple[Point2, ...]:
    pts = sorted(points)
    base = pts[0]
    return tuple(sorted(sub2(p, base) for p in pts))


def span2(points: tuple[Point2, ...]) -> int:
    qs = [p[0] for p in points]
    rs = [p[1] for p in points]
    ss = [-p[0] - p[1] for p in points]
    return (max(qs) - min(qs)) + (max(rs) - min(rs)) + (max(ss) - min(ss))


def triple(p: Point2) -> tuple[int, int, int]:
    return (p[0], p[1], -p[0] - p[1])


@lru_cache(maxsize=None)
def d6_maps() -> tuple[tuple[int, int, int, int], ...]:
    maps: list[tuple[int, int, int, int]] = []
    # encoded as sign, i, j, k where transformed triple is sign*(old[i],old[j],old[k])
    for perm in permutations(range(3)):
        for sign in (1, -1):
            maps.append((sign, perm[0], perm[1], perm[2]))
    return tuple(maps)


def canon2_d6(points: tuple[Point2, ...]) -> tuple[Point2, ...]:
    reps = []
    triples = [triple(p) for p in points]
    for sign, i, j, _ in d6_maps():
        transformed = tuple((sign * t[i], sign * t[j]) for t in triples)
        reps.append(canon2(transformed))
    return min(reps)


def moser_unit_vectors() -> tuple[Point4, ...]:
    units: list[Point4] = []
    for a, b, c, d in product(range(-4, 5), repeat=4):
        if a == b == c == d == 0:
            continue
        if a * d != b * c:
            continue
        value = (
            6 * a * a
            + 6 * a * b
            + 10 * a * c
            + 5 * a * d
            + 6 * b * b
            + 5 * b * c
            + 10 * b * d
            + 6 * c * c
            + 6 * c * d
            + 6 * d * d
        )
        if value == 6:
            units.append((a, b, c, d))
    return tuple(sorted(units))


MOSER_UNITS = moser_unit_vectors()


def add4(a: Point4, b: Point4) -> Point4:
    return tuple(a[i] + b[i] for i in range(4))  # type: ignore[return-value]


def sub4(a: Point4, b: Point4) -> Point4:
    return tuple(a[i] - b[i] for i in range(4))  # type: ignore[return-value]


def canon4(points: tuple[Point4, ...]) -> tuple[Point4, ...]:
    pts = sorted(points)
    base = pts[0]
    return tuple(sorted(sub4(p, base) for p in pts))


def span4(points: tuple[Point4, ...]) -> int:
    return sum(max(p[i] for p in points) - min(p[i] for p in points) for i in range(4))


def antipodal_reps(units: tuple) -> tuple:
    reps = []
    seen = set()
    for u in units:
        neg = tuple(-x for x in u)
        if u in seen or neg in seen:
            continue
        reps.append(min(u, neg))
        seen.add(u)
        seen.add(neg)
    return tuple(sorted(reps))


def drop_direction(units: tuple, rep: tuple) -> tuple:
    neg = tuple(-x for x in rep)
    return tuple(u for u in units if u != rep and u != neg)


@dataclass(frozen=True)
class State:
    points: tuple
    score_edges: int
    true_edges: int
    span: int
    last_gain: int
    gain_hist: tuple[tuple[int, int], ...]


@dataclass(frozen=True)
class BeamResult:
    carrier: str
    policy: str
    target: int
    width: int
    score_edges: int
    true_edges: int
    span: int
    kept: int
    unique_children: int
    frontier_evals: int
    gain_hist: tuple[tuple[int, int], ...]
    cluster: tuple


def edge_count(points: tuple, units: tuple, add) -> int:
    s = set(points)
    return sum(1 for p in points for u in units if add(p, u) in s) // 2


def frontier_gains(points: tuple, units: tuple, add) -> dict[tuple, int]:
    s = set(points)
    gains: dict[tuple, int] = {}
    for p in points:
        for u in units:
            q = add(p, u)
            if q not in s:
                gains[q] = gains.get(q, 0) + 1
    return gains


def gain_at(q: tuple, points: tuple, units: tuple, add) -> int:
    s = set(points)
    return sum(1 for u in units if add(q, u) in s)


def hist_tuple(values: list[int]) -> tuple[tuple[int, int], ...]:
    hist: dict[int, int] = {}
    for value in values:
        hist[value] = hist.get(value, 0) + 1
    return tuple(sorted(hist.items()))


def shell_direction_profile(points: tuple, reps: tuple, add) -> tuple[int, ...]:
    s = set(points)
    profile = []
    for rep in reps:
        neg = tuple(-x for x in rep)
        count = 0
        for p in points:
            if add(p, rep) in s:
                count += 1
            if add(p, neg) in s:
                count += 1
        profile.append(count // 2)
    return tuple(profile)


def profile_imbalance(profile: tuple[int, ...]) -> int:
    if not profile:
        return 0
    return max(profile) - min(profile)


def beam_search(
    carrier: str,
    target: int,
    width: int,
    score_units: tuple,
    true_units: tuple,
    add,
    canon,
    span,
    policy: str,
    gain_ceiling: int | None = None,
) -> BeamResult:
    zero = (0, 0) if carrier == "triangular" else (0, 0, 0, 0)
    start = (zero,)
    start_state = State(start, 0, 0, 0, 0, ())
    states: dict[tuple, State] = {start: start_state}
    unique_children = 0
    frontier_evals = 0

    reps = antipodal_reps(true_units)

    for size in range(2, target + 1):
        children: dict[tuple, State] = {}
        frontier_evals = 0
        for state in states.values():
            gains = frontier_gains(state.points, score_units, add)
            frontier_evals += len(gains)
            for q, gain in gains.items():
                if gain_ceiling is not None and gain > gain_ceiling and size > 3:
                    continue
                child = canon(state.points + (q,))
                score_edges = state.score_edges + gain
                true_edges = state.true_edges + gain_at(q, state.points, true_units, add)
                child_state = State(
                    child,
                    score_edges,
                    true_edges,
                    span(child),
                    gain,
                    (),
                )
                old = children.get(child)
                if old is None or rank_key(child_state, reps, add, policy) > rank_key(old, reps, add, policy):
                    children[child] = child_state
        if not children:
            break
        unique_children = len(children)
        ranked = sorted(children.values(), key=lambda st: rank_key(st, reps, add, policy), reverse=True)
        states = {st.points: st for st in ranked[:width]}

    best = max(states.values(), key=lambda st: (st.true_edges, st.score_edges, -st.span, st.points))
    best_gain_hist = hist_tuple(list(frontier_gains(best.points, score_units, add).values()))
    return BeamResult(
        carrier,
        policy if gain_ceiling is None else f"{policy}+gain_ceiling_{gain_ceiling}",
        target,
        width,
        best.score_edges,
        best.true_edges,
        best.span,
        len(states),
        unique_children,
        frontier_evals,
        best_gain_hist,
        best.points,
    )


def rank_key(state: State, reps: tuple, add, policy: str) -> tuple:
    profile = shell_direction_profile(state.points, reps, add)
    imbalance = profile_imbalance(profile)
    max_frontier_gain = state.last_gain
    frontier_high = 1 if state.last_gain >= 3 else 0
    if policy == "healthy":
        return (state.score_edges, max_frontier_gain, frontier_high, -state.span, -imbalance, state.points)
    if policy == "edge_only":
        return (state.score_edges, state.points)
    if policy == "sprawl_bias":
        return (state.score_edges, state.span, max_frontier_gain, state.points)
    if policy == "balance_bias":
        return (state.score_edges, -imbalance, max_frontier_gain, -state.span, state.points)
    if policy == "future_gain_bias":
        return (state.score_edges + max_frontier_gain, frontier_high, -state.span, state.points)
    raise ValueError(policy)


def print_inspiration() -> None:
    print("INSPIRATION IMPORTED FROM OTHER THREADS")
    print("---------------------------------------")
    print("- LRC: a quotient is useful only when side channels are retained; impairments test which labels are load-bearing.")
    print("- Tournaments: scalar counts miss forbidden packet structure; impairments should target packets, not only totals.")
    print("- Cauldrons: hostile turn schedules expose which resource a one-player process forgot to protect.")
    print("- Equidecomposability: equal edge counts can be non-equivalent; keep construction class and deletion deck.")
    print("- S617 frontier gains: update the local observable state-locally; impaired methods recompute or forget gain shape.")
    print()


def print_width_sweep() -> None:
    print("SMALL WIDTH IMPAIRMENT: MOSER CARRIER")
    print("-------------------------------------")
    print("target=14, policy=healthy; width is the impaired resource.")
    for width in [1, 3, 10, 30, 60]:
        result = beam_search(
            "moser",
            14,
            width,
            MOSER_UNITS,
            MOSER_UNITS,
            add4,
            canon4,
            span4,
            "healthy",
        )
        print(
            f"  width={width:3d}: true_edges={result.true_edges:2d}, "
            f"span={result.span:2d}, frontier_gain_hist={dict(result.gain_hist)}"
        )
    print("Reading: small widths are a deliberate adversary.  If the best edge count")
    print("stabilizes early but the gain histogram does not, the proof object should")
    print("retain future-gain shape rather than only current edge count.")
    print()


def print_policy_impairments() -> None:
    print("POLICY IMPAIRMENTS AT SMALL SIZE")
    print("--------------------------------")
    jobs = [
        ("triangular", 12, 80, TRI_UNITS, TRI_UNITS, add2, canon2, span2),
        ("moser", 12, 60, MOSER_UNITS, MOSER_UNITS, add4, canon4, span4),
    ]
    policies = ["healthy", "edge_only", "sprawl_bias", "balance_bias", "future_gain_bias"]
    for carrier, target, width, units, true_units, add, canon, span in jobs:
        print(f"{carrier} target={target} width={width}")
        baseline = None
        for policy in policies:
            result = beam_search(carrier, target, width, units, true_units, add, canon, span, policy)
            if baseline is None:
                baseline = result.true_edges
            loss = baseline - result.true_edges
            print(
                f"  {policy:16s} true_edges={result.true_edges:2d} "
                f"loss_vs_healthy={loss:2d} span={result.span:2d} "
                f"gain_hist={dict(result.gain_hist)}"
            )
        print()


def print_direction_drop() -> None:
    print("DIRECTION-DROP JACKKNIFE")
    print("------------------------")
    print("Build with one antipodal direction removed; evaluate edges in the full shell.")
    for carrier, target, width, units, add, canon, span in [
        ("triangular", 12, 50, TRI_UNITS, add2, canon2, span2),
        ("moser", 10, 40, MOSER_UNITS, add4, canon4, span4),
    ]:
        full = beam_search(carrier, target, width, units, units, add, canon, span, "healthy")
        print(f"{carrier}: full-shell true_edges={full.true_edges}")
        for rep in antipodal_reps(units):
            kept = drop_direction(units, rep)
            result = beam_search(carrier, target, width, kept, units, add, canon, span, "healthy")
            print(
                f"  drop {rep}: score_edges={result.score_edges:2d}, "
                f"true_edges={result.true_edges:2d}, loss={full.true_edges - result.true_edges:2d}, "
                f"span={result.span:2d}"
            )
        print()


def print_gain_ceiling_impairment() -> None:
    print("GAIN-CEILING IMPAIRMENT")
    print("-----------------------")
    print("Forbid high-gain extensions after size 3; this sabotages dense-core growth.")
    for ceiling in [1, 2, 3, 4]:
        result = beam_search(
            "moser",
            12,
            50,
            MOSER_UNITS,
            MOSER_UNITS,
            add4,
            canon4,
            span4,
            "healthy",
            gain_ceiling=ceiling,
        )
        print(
            f"  ceiling={ceiling}: true_edges={result.true_edges:2d}, "
            f"span={result.span:2d}, gain_hist={dict(result.gain_hist)}"
        )
    print("Repair idea: a proof search should enumerate where gain-4/5 extensions are")
    print("possible, then attach embeddability and unfaithful-obstruction labels there.")
    print()


def print_canonicalization_waste() -> None:
    print("CANONICALIZATION IMPAIRMENT")
    print("---------------------------")
    print("Triangular carrier target=10: translation-only canonicalization versus D6.")
    trans = beam_search("triangular", 10, 200, TRI_UNITS, TRI_UNITS, add2, canon2, span2, "healthy")
    d6 = beam_search("triangular", 10, 200, TRI_UNITS, TRI_UNITS, add2, canon2_d6, span2, "healthy")
    ratio = trans.unique_children / max(1, d6.unique_children)
    print(
        f"  translation-only: unique_children={trans.unique_children}, "
        f"kept={trans.kept}, edges={trans.true_edges}, span={trans.span}"
    )
    print(
        f"  D6 canonical:     unique_children={d6.unique_children}, "
        f"kept={d6.kept}, edges={d6.true_edges}, span={d6.span}"
    )
    print(f"  duplicate-work ratio at last layer: {ratio:.2f}x")
    print("Repair idea: before scaling, separate mathematical impairment from symmetry")
    print("waste.  Otherwise a beam-width failure may just be duplicate orbit traffic.")
    print()


@dataclass(frozen=True)
class Lens:
    name: str
    small_test: int
    proof_payload: int
    construction_payload: int
    side_channel: int
    scaling_risk: int
    note: str


LENSES = (
    Lens("extension gain ledger", 5, 5, 5, 5, 2, "local degree/gain packets before embedding"),
    Lens("direction-drop jackknife", 5, 4, 4, 4, 2, "unit-shell indispensability by leave-one-direction-out tests"),
    Lens("deletion-core impairment", 4, 5, 3, 5, 2, "which dense cores survive hostile vertex deletion"),
    Lens("canonical-orbit repair", 5, 3, 4, 4, 1, "distinguish true search scarcity from duplicate state traffic"),
    Lens("gain-ceiling adversary", 5, 4, 4, 4, 3, "force search to reveal where high-gain steps are necessary"),
    Lens("count-only quotient", 3, 1, 2, 1, 5, "same edge count but missing geometry and obstruction labels"),
    Lens("triangular-only baseline", 5, 1, 2, 1, 4, "good sanity carrier but too thin for n=22"),
    Lens("raw wider beam", 4, 2, 4, 2, 4, "more compute without a new certificate invariant"),
)


def lens_edges() -> list[list[int]]:
    n = len(LENSES)
    adj = [[0] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        a = LENSES[i]
        b = LENSES[j]
        av = sum(
            [
                a.small_test > b.small_test,
                a.proof_payload > b.proof_payload,
                a.construction_payload > b.construction_payload,
                a.side_channel > b.side_channel,
                a.scaling_risk < b.scaling_risk,
            ]
        )
        bv = sum(
            [
                b.small_test > a.small_test,
                b.proof_payload > a.proof_payload,
                b.construction_payload > a.construction_payload,
                b.side_channel > a.side_channel,
                b.scaling_risk < a.scaling_risk,
            ]
        )
        if av > bv or (av == bv and i < j):
            adj[i][j] = 1
        else:
            adj[j][i] = 1
    return adj


def tournament_fingerprints(adj: list[list[int]]) -> tuple[dict[int, int], int, list[int], int]:
    hist: dict[int, int] = {}
    for row in adj:
        score = sum(row)
        hist[score] = hist.get(score, 0) + 1
    c3 = 0
    for i, j, k in combinations(range(len(adj)), 3):
        if adj[i][j] and adj[j][k] and adj[k][i]:
            c3 += 1
        if adj[i][k] and adj[k][j] and adj[j][i]:
            c3 += 1

    n = len(adj)
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        for last in range(n):
            if not dp[mask][last]:
                continue
            for nxt in range(n):
                if not (mask >> nxt) & 1 and adj[last][nxt]:
                    dp[mask | (1 << nxt)][nxt] += dp[mask][last]
    hp = sum(dp[-1])
    scores = [sum(row) for row in adj]
    return dict(sorted(hist.items())), c3, scores, hp


def print_tournament_analysis() -> None:
    print("TOURNAMENT ANALYSIS: IMPAIRMENT LENSES")
    print("--------------------------------------")
    print("Vertices: impairment/repair lenses, not points or raw graphs.")
    print("Pairwise observable: which lens best preserves small-test signal, proof payload,")
    print("construction payload, side-channel retention, and low scaling risk.")
    adj = lens_edges()
    hist, c3, scores, hp = tournament_fingerprints(adj)
    print(f"Score histogram: {hist}")
    print(f"Directed 3-cycles: {c3}")
    print(f"Hamiltonian-path count: {hp}")
    for idx, lens in sorted(enumerate(LENSES), key=lambda item: (-scores[item[0]], item[1].name)):
        print(f"  score={scores[idx]} {lens.name}: {lens.note}")
    print()


def print_assumption_challenge() -> None:
    print("ASSUMPTION CHALLENGE")
    print("--------------------")
    print("Alternate vertex sets considered: points, unit directions, frontier candidates,")
    print("gain packets, deletion cores, obstruction filters, construction classes, and")
    print("proof obligations.  The script uses impairment lenses as Tournament Analysis")
    print("vertices because the preserved predicate is method robustness, not a graph.")
    print("Preserved predicate: small-size response of dense unit-distance construction")
    print("or proof routes to controlled loss of side information.")
    print("Destroyed information: exact planar embedding outside the carrier, full graph")
    print("isomorphism, and proof-grade unfaithful-subgraph certificates.")
    print("Challenged assumption: the next improvement is just a wider beam.  The small")
    print("tests instead point to gain ledgers, direction jackknives, deletion cores, and")
    print("canonical-orbit repair as separate load-bearing invariants.")
    print()


def print_technique_program() -> None:
    print("NOVEL TECHNIQUE PROGRAM")
    print("-----------------------")
    print("1. Damage-response ledger: run the same small carrier under controlled")
    print("   impairments and record which missing label causes the largest loss.")
    print("2. Direction jackknife certificates: leave out one antipodal unit direction,")
    print("   build anyway, then evaluate in the full shell to identify essential")
    print("   direction packets.")
    print("3. Gain-threshold extension solver: for a dense n-core, enumerate candidate")
    print("   gain-4/5 extensions first; only then pay for embedding and obstruction")
    print("   checks.")
    print("4. Deletion-core resilience score: rank clusters by the multiset of edge")
    print("   losses under vertex deletion, not only by total edges.")
    print("5. Orbit-budget accounting: measure duplicate work under weak canonicalization")
    print("   before interpreting a beam-width failure as mathematical scarcity.")
    print()


def main() -> None:
    sys.stdout.reconfigure(line_buffering=True)
    print("S623 UNIT-DISTANCE IMPAIRMENT ATLAS")
    print("===================================")
    print_inspiration()
    print_width_sweep()
    print_policy_impairments()
    print_direction_drop()
    print_gain_ceiling_impairment()
    print_canonicalization_waste()
    print_tournament_analysis()
    print_assumption_challenge()
    print_technique_program()


if __name__ == "__main__":
    main()
