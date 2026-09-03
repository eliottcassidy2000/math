#!/usr/bin/env python3
"""Exact third-tooth competition audit for THM-4335 renewal chains.

The probe has three purposes.

1.  It rewrites farthest-reach greedy as a prefix-record envelope and checks
    an integer centered-residue criterion for the successor of an addressed
    tooth.
2.  It factors the large transition capacities of the two frozen h=420
    controls into obligation, local-successor, and reached-edge stages.
3.  It tests a sharp obstruction: odd integer multiples of one participant
    are nested third teeth and preserve every two-tooth cover.  It then counts
    exactly how many signed walls one residual can kill, audits the sharp
    quotient-three two-residual equality pattern, and constructs a primitive
    twelve-tail family with 420|h and unbounded reuse of one coprime pair.

All endpoint, activity, and tie decisions use exact Fraction arithmetic.
This is not an LRC(14) census and the counterfamily consists of safe rows.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from math import ceil, floor, gcd


@dataclass(frozen=True)
class Tooth:
    w: int
    n: int
    left: Fraction
    right: Fraction


@dataclass(frozen=True)
class Event:
    u: int
    n: int
    v: int
    m: int
    determinant: int
    q: int
    component: int | None


def tooth(w: int, n: int) -> Tooth:
    return Tooth(
        w,
        n,
        Fraction(14 * n - 1, 14 * w),
        Fraction(14 * n + 1, 14 * w),
    )


def anchor_component(h: int, k: int) -> tuple[Fraction, Fraction]:
    assert 0 <= k < 2 * h
    return Fraction(14 * k + 1, 28 * h), Fraction(14 * k + 13, 28 * h)


def selection_key(t: Tooth) -> tuple[Fraction, Fraction, int, int]:
    """The deterministic farthest-reach order, including exact ties."""
    return (t.right, -t.left, -t.w, -t.n)


def meeting_teeth(w: int, left: Fraction, right: Fraction) -> tuple[Tooth, ...]:
    lo = floor(w * left - Fraction(1, 14))
    hi = ceil(w * right + Fraction(1, 14))
    return tuple(
        t
        for n in range(lo - 1, hi + 2)
        if (t := tooth(w, n)).left < right and left < t.right
    )


def component_teeth(
    h: int, k: int, speeds: tuple[int, ...]
) -> tuple[Tooth, ...]:
    left, right = anchor_component(h, k)
    return tuple(t for w in speeds for t in meeting_teeth(w, left, right))


def active_trace(
    h: int, k: int, speeds: tuple[int, ...]
) -> tuple[str, tuple[Tooth, ...], tuple[Fraction, ...]]:
    """Literal open-tooth greedy trace; retain a partial trace on failure."""
    left, right = anchor_component(h, k)
    teeth = component_teeth(h, k, speeds)
    cursor = left
    chain: list[Tooth] = []
    frontiers: list[Fraction] = []
    while True:
        eligible = [t for t in teeth if t.left < cursor < t.right]
        if not eligible:
            return "missing", tuple(chain), tuple(frontiers)
        winner = max(eligible, key=selection_key)
        frontiers.append(cursor)
        chain.append(winner)
        cursor = winner.right
        if cursor > right:
            status = "span" if len(chain) == 1 else "renew"
            return status, tuple(chain), tuple(frontiers)


def prefix_trace(
    h: int, k: int, speeds: tuple[int, ...]
) -> tuple[str, tuple[Tooth, ...], tuple[Fraction, ...]]:
    """Independent prefix-record implementation of the same greedy walk."""
    left, right = anchor_component(h, k)
    teeth = component_teeth(h, k, speeds)
    cursor = left
    chain: list[Tooth] = []
    frontiers: list[Fraction] = []
    while True:
        started = [t for t in teeth if t.left < cursor]
        if not started:
            return "missing", tuple(chain), tuple(frontiers)
        winner = max(started, key=selection_key)
        if winner.right <= cursor:
            return "missing", tuple(chain), tuple(frontiers)
        frontiers.append(cursor)
        chain.append(winner)
        cursor = winner.right
        if cursor > right:
            status = "span" if len(chain) == 1 else "renew"
            return status, tuple(chain), tuple(frontiers)


def visibility_cap(t: Tooth, teeth: tuple[Tooth, ...]) -> Fraction:
    """Right edge of the frontier interval on which t is the prefix winner."""
    outranking_starts = [
        z.left
        for z in teeth
        if selection_key(z) > selection_key(t) and z.left < t.right
    ]
    return min((t.right, *outranking_starts))


def centered_residue(value: int, modulus: int) -> int:
    """Representative in (-modulus/2, modulus/2], adequate at strict walls."""
    r = value % modulus
    if 2 * r > modulus:
        r -= modulus
    return r


def circle_distance(x: Fraction) -> Fraction:
    residue = x - floor(x)
    return min(residue, 1 - residue)


def wall_danger_count(u: int, s: int) -> int:
    """Number of the 2u signed u-walls strictly killed by speed s."""
    count = 0
    for n in range(u):
        for sigma in (-1, 1):
            e = centered_residue(s * (14 * n + sigma), 14 * u)
            count += int(abs(e) < u)
    return count


def wall_capacity_bound(u: int, s: int) -> int:
    g = gcd(u, s)
    U = u // g
    return 2 * g * ceil(Fraction(U, 7))


def congruence_interval_count(q: int, residue: int) -> int:
    return sum((y - residue) % 14 == 0 for y in range(-q + 1, q))


def wall_danger_count_formula(u: int, s: int) -> int:
    d = gcd(u, s)
    q, reduced = u // d, s // d
    return d * (
        congruence_interval_count(q, reduced)
        + congruence_interval_count(q, -reduced)
    )


def wall_mask(u: int, s: int) -> int:
    mask = 0
    bit = 0
    for n in range(u):
        for sigma in (-1, 1):
            e = centered_residue(s * (14 * n + sigma), 14 * u)
            if abs(e) < u:
                mask |= 1 << bit
            bit += 1
    return mask


def audit_wall_capacity() -> tuple[object, ...]:
    checks = 0
    proper_quotient_checks = 0
    for u in range(1, 64, 2):
        for s in range(1, 14 * u + 1):
            direct = wall_danger_count(u, s)
            assert direct == wall_danger_count_formula(u, s)
            assert direct <= wall_capacity_bound(u, s)
            if u // gcd(u, s) > 1:
                assert 3 * direct <= 2 * u
                proper_quotient_checks += 1
            checks += 1
    sharp = (31, 433, wall_danger_count(31, 433), wall_capacity_bound(31, 433))
    assert sharp[2:] == (10, 10)

    # Complete quotient-three equality audit.  The factor d merely repeats
    # each of the six quotient walls, so testing odd d through 15 also checks
    # the claimed scaling rather than only the primitive u=3 picture.
    equality_checks = 0
    cover_rows = 0
    odd_units_mod_3 = tuple(r for r in range(1, 42, 2) if r % 3)
    for d in range(1, 16, 2):
        u = 3 * d
        full = (1 << (2 * u)) - 1
        for H in range(1, 43):
            if H % 3 == 0:
                continue
            anchor = wall_mask(u, 2 * d * H)
            assert anchor != full
            for R1, R2 in combinations(odd_units_mod_3, 2):
                r1, r2 = d * R1, d * R2
                covered = anchor | wall_mask(u, r1) | wall_mask(u, r2) == full
                class_1 = "A" if R1 in (1, 41) else "B" if R1 in (13, 29) else "-"
                class_2 = "A" if R2 in (1, 41) else "B" if R2 in (13, 29) else "-"
                predicted = H % 21 in (7, 14) and {class_1, class_2} == {"A", "B"}
                assert covered == predicted, (d, H, R1, R2)
                cover_rows += int(covered)
                equality_checks += 1

    # The least positive representatives are 1 and 13.  Replacing them by
    # 41=-1 and 29=-13 modulo 42 leaves the masks unchanged and gives the
    # least representatives which both exceed the anchor 2h=14.
    assert wall_mask(3, 1) == wall_mask(3, 41)
    assert wall_mask(3, 13) == wall_mask(3, 29)
    hostile = (
        3,
        7,
        (1, 13),
        (29, 41),
        wall_mask(3, 14),
        wall_mask(3, 29),
        wall_mask(3, 41),
    )
    assert hostile[-3] | hostile[-2] | hostile[-1] == (1 << 6) - 1
    return checks, proper_quotient_checks, sharp, equality_checks, cover_rows, hostile


def residue_winner(u: int, n: int, speeds: tuple[int, ...]) -> tuple[int, int]:
    """Winner at b(u,n), using only centered integer residues.

    If e_z is z(14n+1) modulo 14u in the centered interval, z is active at
    b(u,n) iff |e_z|<u.  Its right extension, after a common factor 1/(14u),
    is (u-e_z)/z.  Equal right walls are won by the smaller speed.
    """
    candidates: list[tuple[Fraction, int, int]] = []
    value_factor = 14 * n + 1
    for z in speeds:
        value = z * value_factor
        e = centered_residue(value, 14 * u)
        if abs(e) < u:
            ell = (value - e) // (14 * u)
            candidates.append((Fraction(u - e, z), -z, ell))
    assert candidates
    _, minus_z, ell = max(candidates)
    return -minus_z, ell


def geometric_winner(u: int, n: int, speeds: tuple[int, ...]) -> tuple[int, int]:
    cursor = tooth(u, n).right
    candidates: list[Tooth] = []
    for z in speeds:
        centre = z * cursor
        base = floor(centre)
        for ell in range(base - 1, base + 3):
            t = tooth(z, ell)
            if t.left < cursor < t.right:
                candidates.append(t)
    assert candidates
    winner = max(candidates, key=selection_key)
    return winner.w, winner.n


def determinant_indices(u: int, v: int) -> tuple[int, ...]:
    g = gcd(u, v)
    U, V = u // g, v // g
    lo = floor(abs(U - V) / 14) + 1
    hi = floor(Fraction(U + V - 1, 14))
    return tuple(range(lo, hi + 1))


def locate_event(h: int, u: int, n: int, v: int, m: int) -> Event:
    D = u * m - v * n
    assert abs(u - v) < 14 * D < u + v
    q = u + v - 14 * D
    assert q > 0 and q % (2 * gcd(u, v)) == 0

    # [a(v,m), b(u,n)] is allowed to touch either anchor endpoint.
    incoming = 2 * h * (14 * m - 1)
    outgoing = 2 * h * (14 * n + 1)
    k_in, r_in = divmod(incoming, 14 * v)
    k_out, r_out = divmod(outgoing, 14 * u)
    if (
        k_in == k_out
        and v <= r_in < 13 * v
        and u < r_out <= 13 * u
    ):
        assert 0 <= k_out < 2 * h
        component: int | None = k_out
    else:
        component = None
    return Event(u, n, v, m, D, q, component)


def addressed_events(h: int, u: int, v: int) -> tuple[Event, ...]:
    """All period-one proper crossings u->v, with the anchor location filter."""
    g = gcd(u, v)
    U, V = u // g, v // g
    inverse = pow(V, -1, U) if U > 1 else 0
    events: list[Event] = []
    for d in determinant_indices(u, v):
        n0 = (-d * inverse) % U if U > 1 else 0
        m0 = (d + V * n0) // U
        for s in range(g):
            event = locate_event(h, u, n0 + U * s, v, m0 + V * s)
            if event.component is not None:
                events.append(event)
    return tuple(events)


def edge_tuple(k: int, a: Tooth, b: Tooth) -> tuple[int, int, int, int, int, int]:
    return (k, a.w, a.n, b.w, b.n, a.w * b.n - b.w * a.n)


def actual_data(
    h: int, speeds: tuple[int, ...]
) -> tuple[Counter[str], set[tuple[int, int, int, int, int, int]]]:
    statuses: Counter[str] = Counter()
    edges: set[tuple[int, int, int, int, int, int]] = set()
    for k in range(2 * h):
        status, chain, _ = active_trace(h, k, speeds)
        statuses[status] += 1
        if status == "renew":
            for a, b in zip(chain, chain[1:]):
                edges.add(edge_tuple(k, a, b))
    return statuses, edges


def audit_small_prefix_and_residue() -> tuple[int, int, int, int]:
    trace_checks = 0
    selected_edges = 0
    candidate_checks = 0
    visibility_checks = 0
    odds = (1, 3, 5, 7, 9, 11, 13)
    for h in (2, 3, 5, 8, 13):
        for speeds in combinations(odds, 3):
            for k in range(2 * h):
                literal = active_trace(h, k, speeds)
                prefix = prefix_trace(h, k, speeds)
                assert literal == prefix
                trace_checks += 1
                status, chain, frontiers = literal
                teeth = component_teeth(h, k, speeds)
                for t, x in zip(chain, frontiers):
                    cap = visibility_cap(t, teeth)
                    assert t.left < x <= cap and x < t.right
                    visibility_checks += 1
                if status == "renew":
                    for a, b in zip(chain, chain[1:]):
                        event = locate_event(h, a.w, a.n, b.w, b.n)
                        assert event.component == k
                        assert residue_winner(a.w, a.n, speeds) == (b.w, b.n)
                        assert geometric_winner(a.w, a.n, speeds) == (b.w, b.n)
                        selected_edges += 1

            for u in speeds:
                for v in speeds:
                    if u == v:
                        continue
                    for event in addressed_events(h, u, v):
                        assert residue_winner(u, event.n, speeds) == geometric_winner(
                            u, event.n, speeds
                        )
                        candidate_checks += 1
    return trace_checks, selected_edges, candidate_checks, visibility_checks


def audit_nested_multiple_preservation() -> tuple[int, tuple[object, ...]]:
    checks = 0
    multipliers = (3, 5, 7, 9, 11, 13, 15)
    odds = tuple(range(1, 26, 2))
    for h in range(2, 16):
        for u, v in combinations(odds, 2):
            for k in range(2 * h):
                status, chain, _ = active_trace(h, k, (u, v))
                if status != "renew" or len(chain) != 2:
                    continue
                for root in (u, v):
                    for c in multipliers:
                        z = c * root
                        if z in (u, v):
                            continue
                        extended = active_trace(h, k, (u, v, z))
                        assert extended[0] == "renew" and extended[1] == chain
                        checks += 1

    # Divisibility by 14 is the exact boundary-crossing obstruction.  This
    # even-speed example lies outside the odd-tail branch, but shows why the
    # hypothesis cannot be dropped in the general interval lemma.
    base = active_trace(26, 35, (3, 29))
    boundary = active_trace(26, 35, (3, 29, 42))
    assert tuple((t.w, t.n) for t in base[1]) == ((3, 2), (29, 20))
    assert tuple((t.w, t.n) for t in boundary[1]) == ((3, 2), (42, 29))
    witness = (
        26,
        35,
        tuple((t.w, t.n) for t in base[1]),
        tuple((t.w, t.n) for t in boundary[1]),
    )
    return checks, witness


def control_sieve(extra: int) -> None:
    h = 420
    common = tuple(11 + 1680 * k for k in range(7)) + (525, 945, 1365, 1575)
    speeds = common + (extra,)

    statuses, actual = actual_data(h, speeds)
    candidates: list[Event] = []
    for u in speeds:
        for v in speeds:
            if u != v:
                candidates.extend(addressed_events(h, u, v))

    survivor: list[Event] = []
    killer_counts: Counter[int] = Counter()
    groups: defaultdict[tuple[int, int, int], list[Event]] = defaultdict(list)
    status_by_component: dict[int, str] = {}
    chains_by_component: dict[int, tuple[Tooth, ...]] = {}
    for k in range(2 * h):
        status, chain, _ = active_trace(h, k, speeds)
        status_by_component[k] = status
        chains_by_component[k] = chain

    for event in candidates:
        assert event.component is not None
        groups[(event.component, event.u, event.n)].append(event)
        arithmetic = residue_winner(event.u, event.n, speeds)
        geometric = geometric_winner(event.u, event.n, speeds)
        assert arithmetic == geometric
        if arithmetic == (event.v, event.m):
            survivor.append(event)
        else:
            killer_counts[arithmetic[0]] += 1

    # At most one successor survives at any addressed outgoing frontier.
    survivor_by_group: Counter[tuple[int, int, int]] = Counter(
        (e.component, e.u, e.n) for e in survivor
    )
    assert max(survivor_by_group.values(), default=0) <= 1

    candidate_keys = {
        (e.component, e.u, e.n, e.v, e.m, e.determinant) for e in candidates
    }
    survivor_keys = {
        (e.component, e.u, e.n, e.v, e.m, e.determinant) for e in survivor
    }
    assert actual <= candidate_keys
    assert actual <= survivor_keys

    candidate_status = Counter(status_by_component[e.component] for e in candidates)
    survivor_status = Counter(status_by_component[e.component] for e in survivor)
    renewal_candidates = [
        e for e in candidates if status_by_component[e.component] == "renew"
    ]
    renewal_survivors = [
        e for e in survivor if status_by_component[e.component] == "renew"
    ]

    # The exact capacity loss splits into collisions among candidate
    # successors with one outgoing address, and groups whose true winner is
    # not itself a located candidate from that address.
    shared_frontier_excess = len(candidates) - len(groups)
    empty_winner_groups = len(groups) - len(survivor)
    assert len(candidates) - len(survivor) == (
        shared_frontier_excess + empty_winner_groups
    )

    print(f"h420_P{extra}")
    print(" component_status", tuple(sorted(statuses.items())))
    print(
        " sieve_all_obligation_successor_reached",
        (len(candidates), len(renewal_candidates), len(renewal_survivors), len(actual)),
    )
    print(
        " local_partition_candidate_groups_survivors_shared_excess_empty_groups",
        (
            len(candidates),
            len(groups),
            len(survivor),
            shared_frontier_excess,
            empty_winner_groups,
        ),
    )
    print(" candidate_status", tuple(sorted(candidate_status.items())))
    print(" successor_status", tuple(sorted(survivor_status.items())))
    print(" top_chosen_killers", tuple(killer_counts.most_common()))

    for u, v in ((525, 1365), (1365, 525)):
        pair_candidates = [e for e in candidates if (e.u, e.v) == (u, v)]
        pair_survivors = [e for e in survivor if (e.u, e.v) == (u, v)]
        pair_actual = [row for row in actual if (row[1], row[3]) == (u, v)]
        pair_killers: Counter[int] = Counter()
        for event in pair_candidates:
            winner = residue_winner(event.u, event.n, speeds)
            if winner != (event.v, event.m):
                pair_killers[winner[0]] += 1
        print(
            f" pair_{u}_{v}_candidate_obligation_successor_reached",
            (
                len(pair_candidates),
                sum(status_by_component[e.component] == "renew" for e in pair_candidates),
                sum(status_by_component[e.component] == "renew" for e in pair_survivors),
                len(pair_actual),
            ),
        )
        print(
            f" pair_{u}_{v}_candidate_status",
            tuple(
                sorted(Counter(status_by_component[e.component] for e in pair_candidates).items())
            ),
        )
        print(f" pair_{u}_{v}_chosen_killers", tuple(sorted(pair_killers.items())))

    print(" actual_edges", tuple(sorted(actual)))


def full_twelve_tail_counterfamily() -> None:
    # 9u,11u,13u pay the remaining THM-4330 divisibility necessities;
    # 39u,...,47u lie in the forbidden half-turn residue band.
    multipliers = (9, 11, 13, 39, 41, 43, 45, 47, 49, 51)
    print("full_twelve_tail_nested_counterfamily")
    for M in (1, 2, 5):
        N = 10 * M
        h = 420 * M
        u, v = 140 * M + 1, 140 * M + 3
        speeds = (u, v) + tuple(c * u for c in multipliers)
        assert len(speeds) == len(set(speeds)) == 12
        assert all(w % 2 == 1 for w in speeds)
        assert gcd(u, v) == 1
        assert gcd(2 * h, *speeds) == 1
        assert h % 420 == 0
        assert all(any(w % modulus == 0 for w in speeds) for modulus in (9, 11, 13))

        base_status, base_edges = actual_data(h, (u, v))
        full_status, full_edges = actual_data(h, speeds)
        base_pair = {row for row in base_edges if {row[1], row[3]} == {u, v}}
        full_pair = {row for row in full_edges if {row[1], row[3]} == {u, v}}
        assert base_pair <= full_pair

        # The inherited explicit determinant range d=N+1,...,2N occurs in
        # each orientation and gives the rigorous 2N lower bound.
        for d in range(N + 1, 2 * N + 1):
            assert any(row[1] == u and row[3] == v and row[5] == d for row in base_pair)
            assert any(row[1] == v and row[3] == u and row[5] == d for row in base_pair)

        bad_multipliers = tuple(
            c
            for c in multipliers
            if 12 * h < (c * u) % (28 * h) < 16 * h
        )
        assert bad_multipliers
        # Every non-14-divisible multiple of u is equality-safe or better at
        # an u-wall.  Here the first positive u-wall also leaves 2h and v
        # safe, so it is a uniform exact witness for the whole family.
        wall_witness = Fraction(1, 14 * u)
        clearance = min(circle_distance(w * wall_witness) for w in (2 * h, *speeds))
        assert clearance == Fraction(1, 14)
        residual_wall_capacity = wall_capacity_bound(u, 2 * h) + wall_capacity_bound(u, v)
        assert residual_wall_capacity < 2 * u
        print(
            " M_h_u_v",
            (M, h, u, v),
            "base_pair_full_pair_lower",
            (len(base_pair), len(full_pair), 2 * N),
            "full_status",
            tuple(sorted(full_status.items())),
            "halfturn_bad_multipliers",
            bad_multipliers,
            "wall_witness_clearance",
            (wall_witness, clearance),
            "anchor_plus_residual_wall_capacity",
            (residual_wall_capacity, 2 * u),
        )


def main() -> None:
    print("small_prefix_edge_candidate_visibility_checks", audit_small_prefix_and_residue())
    print("nested_preservation_checks_and_c14_boundary", audit_nested_multiple_preservation())
    print("wall_capacity_checks_and_sharp_control", audit_wall_capacity())
    control_sieve(1287)
    control_sieve(9009)
    full_twelve_tail_counterfamily()
    print("PASS")


if __name__ == "__main__":
    main()
