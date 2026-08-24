#!/usr/bin/env python3
"""Exact independent THM-3910 audit of the compact-to-open response filter.

The computation is intentionally standalone.  It rebuilds AP11's closed
safe-set wall graph and every scale-one pair-danger component from rational
endpoints.  The analytic implication being audited is:

    lambda_+(u) = longest positive-length closed component of
                  {y: min_i ||u_i y|| >= 1/14}.

If t >= U=max(u), a THM-3878 scale-one counterexample must obey

    t*lambda_+(u) < beta_1(p,q),

and the scale-two odd-pair branch obeys the analogous strict inequality with
beta_2.  Indeed a closed component maps continuously to a compact connected
subset of the relevant *open* danger set.  It must fit strictly inside one
open component.  Equality of lengths already rules out containment.  Isolated
closed-safe equality walls have length zero and supply no such response.

This is a necessary filter for a counterexample, not a proof of LRC(14).
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
import json
import math
import sys


sys.dont_write_bytecode = True
sys.stdout.reconfigure(newline="\n")

H = Q(1, 14)
AP11 = tuple(range(1, 12))

SCALE1 = (
    (1, 3), (1, 4), (1, 9), (1, 10),
    (2, 3), (2, 9), (2, 15), (2, 21), (2, 23),
    (3, 7), (3, 8), (3, 14), (3, 17), (3, 19), (3, 20),
    (3, 22), (3, 26), (3, 31), (3, 38),
    (4, 7), (4, 13), (4, 19), (4, 21), (4, 25), (4, 37),
    (4, 43), (4, 49), (4, 51),
    (5, 6), (5, 12), (5, 17), (5, 18), (5, 24), (5, 29),
    (5, 36), (5, 39), (5, 41), (5, 42), (5, 48), (5, 53),
    (5, 54), (5, 63),
    (6, 11), (6, 17), (6, 19), (6, 23), (6, 41), (6, 47),
    (6, 53), (6, 65),
    (7, 10), (7, 13), (7, 15), (7, 22),
    (8, 9), (8, 21), (9, 11),
)

# (beta, multiplicity, cumulative number closed at U*lambda >= beta).
EXPECTED_BETA_STAIRCASE = (
    (Q(19, 784), 1, 1), (Q(8, 329), 1, 2),
    (Q(13, 532), 1, 3), (Q(6, 245), 1, 4),
    (Q(17, 693), 1, 5), (Q(47, 1820), 1, 6),
    (Q(115, 4452), 1, 7), (Q(37, 1428), 1, 8),
    (Q(57, 2156), 1, 9), (Q(17, 637), 1, 10),
    (Q(8, 287), 1, 11), (Q(1, 35), 4, 15),
    (Q(37, 1260), 1, 16), (Q(19, 644), 1, 17),
    (Q(18, 595), 1, 18), (Q(89, 2940), 1, 19),
    (Q(23, 756), 1, 20), (Q(31, 1015), 1, 21),
    (Q(44, 1435), 1, 22), (Q(57, 1855), 1, 23),
    (Q(31, 1008), 1, 24), (Q(43, 1365), 1, 25),
    (Q(31, 980), 1, 26), (Q(31, 924), 1, 27),
    (Q(19, 560), 1, 28), (Q(1, 28), 6, 34),
    (Q(31, 840), 1, 35), (Q(89, 2408), 1, 36),
    (Q(31, 728), 1, 37), (Q(8, 175), 1, 38),
    (Q(1, 21), 8, 46), (Q(19, 364), 1, 47),
    (Q(31, 588), 1, 48), (Q(1, 14), 4, 52),
    (Q(8, 105), 1, 53), (Q(1, 7), 4, 57),
)

COHERENT_Q25 = tuple(26 * k for k in range(1, 13)) + (339,)
INCOHERENT_Q25 = (81, 91, 131, 151, 157, 196, 258, 274, 313, 328, 330, 339, 348)

CHECKS = 0


def require(condition: bool, label: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(label)


def fmt(x: Q) -> str:
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def circle_distance(x: Q) -> Q:
    r = x % 1
    return min(r, 1 - r)


def safe_closed(base: tuple[int, ...], x: Q) -> bool:
    return all(circle_distance(v * x) >= H for v in base)


def safe_strict(base: tuple[int, ...], x: Q) -> bool:
    return all(circle_distance(v * x) > H for v in base)


def equality_owners(base: tuple[int, ...], x: Q) -> tuple[tuple[int, ...], tuple[int, ...]]:
    """Runners entering danger immediately to the left and to the right."""
    left = tuple(v for v in base if (v * x) % 1 == H)
    right = tuple(v for v in base if (v * x) % 1 == 1 - H)
    return left, right


def wall_graph(base: tuple[int, ...]):
    walls = sorted({
        ((Q(k) + sign * H) / v) % 1
        for v in base for k in range(v) for sign in (-1, 1)
    })
    n = len(walls)
    cell_safe = []
    for i, left in enumerate(walls):
        right = walls[(i + 1) % n]
        right_lift = right if right > left else right + 1
        midpoint = ((left + right_lift) / 2) % 1
        require(
            safe_closed(base, midpoint) == safe_strict(base, midpoint),
            f"cell convention {i}",
        )
        cell_safe.append(safe_strict(base, midpoint))
    wall_safe = [safe_closed(base, x) for x in walls]

    vertices = [("w", i) for i in range(n) if wall_safe[i]]
    vertices += [("c", i) for i in range(n) if cell_safe[i]]
    adjacency = {vertex: set() for vertex in vertices}
    for i in range(n):
        if not cell_safe[i]:
            continue
        for wall_id in (i, (i + 1) % n):
            require(wall_safe[wall_id], f"safe cell endpoint {i}:{wall_id}")
            adjacency[("c", i)].add(("w", wall_id))
            adjacency[("w", wall_id)].add(("c", i))

    components = []
    unseen = set(vertices)
    while unseen:
        seed = min(unseen)
        unseen.remove(seed)
        stack = [seed]
        component = set()
        while stack:
            vertex = stack.pop()
            component.add(vertex)
            for neighbor in adjacency[vertex]:
                if neighbor in unseen:
                    unseen.remove(neighbor)
                    stack.append(neighbor)
        components.append(component)

    isolated = sorted(
        walls[i] for i in range(n)
        if wall_safe[i] and not adjacency[("w", i)]
    )
    arcs = []
    for component in components:
        cells = sorted(i for kind, i in component if kind == "c")
        if not cells:
            continue
        length = sum(
            ((walls[(i + 1) % n] - walls[i]) % 1 for i in cells), Q(0)
        )
        arcs.append((walls[cells[0]], walls[(cells[-1] + 1) % n], length))
    return walls, cell_safe, wall_safe, components, isolated, arcs


def raw_danger(w: int) -> list[tuple[Q, Q]]:
    radius = Q(1, 14 * w)
    pieces = [(Q(0), radius), (Q(1) - radius, Q(1))]
    pieces.extend(
        (Q(k, w) - radius, Q(k, w) + radius)
        for k in range(1, w)
    )
    return pieces


def circle_open_components(frequencies: tuple[int, ...]) -> tuple[tuple[Q, Q], ...]:
    """Lifted open components; equality contacts remain separate."""
    pieces = sorted(sum((raw_danger(w) for w in frequencies), []))
    merged: list[list[Q]] = []
    for left, right in pieces:
        if not merged or left >= merged[-1][1]:
            merged.append([left, right])
        elif right > merged[-1][1]:
            merged[-1][1] = right
    if len(merged) >= 2 and merged[0][0] == 0 and merged[-1][1] == 1:
        wrap = (merged[-1][0], merged[0][1] + 1)
        return (wrap,) + tuple((left, right) for left, right in merged[1:-1])
    return tuple((left, right) for left, right in merged)


def danger(x: Q, p: int, q: int) -> bool:
    return circle_distance(p * x) < H or circle_distance(q * x) < H


def scale2_quotient_components(p: int, q: int) -> tuple[tuple[Q, Q], ...]:
    """Build C_pq directly on w=2z, without importing THM-3878 code."""
    z_walls = {
        ((Q(k) + sign * H) / v) % 1
        for v in (p, q) for k in range(v) for sign in (-1, 1)
    }
    walls = sorted({Q(0)} | {(2 * x) % 1 for x in z_walls})
    active = []
    for i, left in enumerate(walls):
        right = walls[(i + 1) % len(walls)]
        right_lift = right if right > left else right + 1
        w_mid = ((left + right_lift) / 2) % 1
        z = w_mid / 2
        active.append(danger(z, p, q) and danger(z + Q(1, 2), p, q))

    pieces = []
    for i, is_active in enumerate(active):
        if not is_active:
            continue
        left = walls[i]
        right = walls[(i + 1) % len(walls)]
        pieces.append((left, right if right > left else right + 1))

    # Every nonzero listed wall is excluded by at least one strict danger
    # equality.  Only the artificial cut at zero can join the end cells.
    if active[0] and active[-1]:
        first = pieces.pop(0)
        last = pieces.pop(-1)
        pieces.insert(0, (last[0], first[1] + 1))
    return tuple(pieces)


def beta_staircase() -> tuple[tuple[Q, int, int], ...]:
    counts = Counter(
        max(right - left for left, right in circle_open_components(pair))
        for pair in SCALE1
    )
    cumulative = 0
    result = []
    for beta, multiplicity in sorted(counts.items()):
        cumulative += multiplicity
        result.append((beta, multiplicity, cumulative))
    return tuple(result)


def clearance(speeds: tuple[int, ...], x: Q) -> Q:
    return min(circle_distance(v * x) for v in speeds)


def least_rational_witness(speeds: tuple[int, ...], max_denominator: int):
    for denominator in range(2, max_denominator + 1):
        for numerator in range(1, denominator):
            x = Q(numerator, denominator)
            if clearance(speeds, x) >= H:
                return denominator, numerator, x
    return None


def staircase_text(staircase: tuple[tuple[Q, int, int], ...]) -> str:
    return ",".join(f"{fmt(beta)}:{multiplicity}:{cumulative}" for beta, multiplicity, cumulative in staircase)


def main() -> None:
    require(len(SCALE1) == 57 and len(set(SCALE1)) == 57, "57 scale-one rows")
    require(all(math.gcd(p, q) == 1 and p < q for p, q in SCALE1), "primitive pairs")

    walls, cells, safe_walls, components, isolated, arcs = wall_graph(AP11)
    lengths = sorted(length for _, _, length in arcs)
    require(len(walls) == 128, "AP11 equality walls")
    require(sum(cells) == 14 and len(arcs) == 14, "AP11 positive arcs")
    require(sum(safe_walls) == 32 and len(components) == 18, "AP11 closed topology")
    require(isolated == [Q(3, 14), Q(5, 14), Q(9, 14), Q(11, 14)], "AP11 isolated")
    require(sum(lengths, Q(0)) == Q(10931, 194040), "AP11 measure")
    require(max(lengths) == Q(1, 77), "AP11 lambda plus")
    require(Counter(lengths) == Counter({
        Q(1, 588): 2, Q(1, 560): 2, Q(1, 504): 2,
        Q(1, 420): 2, Q(1, 308): 2, Q(1, 245): 2,
        Q(1, 77): 2,
    }), "AP11 component spectrum")

    J = (Q(1, 14), Q(13, 154))
    require(J[1] - J[0] == Q(1, 77), "J length")
    require(any(left == J[0] and right == J[1] for left, right, _ in arcs), "J is closed arc")
    require(equality_owners(AP11, J[0]) == ((1,), ()), "J left owner")
    require(equality_owners(AP11, J[1]) == ((), (11,)), "J right owner")
    require(safe_strict(AP11, sum(J, Q(0)) / 2), "J interior strict")
    isolated_owners = tuple((x, *equality_owners(AP11, x)) for x in isolated)
    require(isolated_owners == (
        (Q(3, 14), (5,), (9,)), (Q(5, 14), (3,), (11,)),
        (Q(9, 14), (11,), (3,)), (Q(11, 14), (9,), (5,)),
    ), "AP11 isolated owner switches")

    staircase = beta_staircase()
    require(staircase == EXPECTED_BETA_STAIRCASE, "full beta staircase")
    betas = {
        pair: max(right - left for left, right in circle_open_components(pair))
        for pair in SCALE1
    }
    beta_max_pairs = tuple(pair for pair in SCALE1 if betas[pair] == Q(1, 7))
    require(beta_max_pairs == ((1, 3), (1, 4), (1, 9), (1, 10)), "beta max equality rows")
    require(min(betas.values()) == Q(19, 784) > Q(1, 42), "generic LRC13 closes zero")
    selected = tuple(
        (threshold, sum(beta <= threshold for beta in betas.values()))
        for threshold in (Q(1, 35), Q(1, 28), Q(1, 21), Q(1, 14), Q(1, 7))
    )
    require(selected == (
        (Q(1, 35), 15), (Q(1, 28), 34), (Q(1, 21), 46),
        (Q(1, 14), 52), (Q(1, 7), 57),
    ), "selected response staircase")

    scale2 = scale2_quotient_components(1, 9)
    require(scale2 == ((Q(2, 21), Q(8, 63)), (Q(55, 63), Q(19, 21))), "scale-two components")
    beta2 = max(right - left for left, right in scale2)
    require(beta2 == Q(2, 63), "scale-two beta")

    # AP11 has U=11 and lambda_+=1/77, so U*lambda=1/7.
    # Compact-to-open equality closes even the four beta=1/7 rows.
    ap_response = 11 * Q(1, 77)
    require(ap_response == Q(1, 7), "AP11 response")
    require(all(ap_response >= beta for beta in betas.values()), "AP11 closes 57")
    require(ap_response > beta2, "AP11 closes scale two")

    # Canonical tight AP13: only isolated unit walls, with opposing owners.
    ap13 = tuple(range(1, 14))
    ap13_points = tuple(Q(k, 14) for k in (1, 3, 5, 9, 11, 13))
    require(all(safe_closed(ap13, x) for x in ap13_points), "AP13 unit points safe")
    require(all(not safe_strict(ap13, x) for x in ap13_points), "AP13 unit points equality")
    ap13_owners = tuple((x, *equality_owners(ap13, x)) for x in ap13_points)
    require(all(len(left) == len(right) == 1 for _, left, right in ap13_owners), "AP13 opposing owners")

    # q<=25 hostile controls: endpoint response is information that a bounded
    # rational scan misses.  The displayed witnesses are not claimed least;
    # the independently recovered least denominators are 27 and 26.
    coherent_least = least_rational_witness(COHERENT_Q25, 27)
    incoherent_least = least_rational_witness(INCOHERENT_Q25, 26)
    require(coherent_least == (27, 2, Q(2, 27)), "coherent first denominator")
    require(incoherent_least == (26, 3, Q(3, 26)), "incoherent first denominator")
    coherent_response = Q(181, 2548)
    incoherent_response = Q(589, 2744)
    require(clearance(COHERENT_Q25, coherent_response) == H, "coherent response clearance")
    require(clearance(INCOHERENT_Q25, incoherent_response) == H, "incoherent response clearance")
    coherent_eq = tuple(v for v in COHERENT_Q25 if circle_distance(v * coherent_response) == H)
    incoherent_eq = tuple(v for v in INCOHERENT_Q25 if circle_distance(v * incoherent_response) == H)
    require(coherent_eq == (182,), "coherent response owner")
    require(incoherent_eq == (196,), "incoherent response owner")
    require(coherent_response == Q(1, 14) - Q(1, 2548), "coherent response location")
    require(incoherent_response == Q(3, 14) + Q(1, 2744), "incoherent response location")

    semantic = {
        "ap11": {
            "closed_components": len(components),
            "positive_arcs": len(arcs),
            "isolated": [fmt(x) for x in isolated],
            "lambda_plus": fmt(max(lengths)),
            "J": [fmt(x) for x in J],
            "isolated_owners": [
                [fmt(x), list(left), list(right)] for x, left, right in isolated_owners
            ],
        },
        "beta_staircase": [[fmt(x), n, cumulative] for x, n, cumulative in staircase],
        "beta_max_pairs": [list(pair) for pair in beta_max_pairs],
        "scale2": [[fmt(left), fmt(right)] for left, right in scale2],
        "q25_controls": {
            "coherent": [fmt(coherent_response), list(coherent_eq)],
            "incoherent": [fmt(incoherent_response), list(incoherent_eq)],
        },
    }
    digest = sha256(
        json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()

    print("LRC14_LAMBDA_COMPONENT_RESPONSE_INDEPENDENT_AUDIT_20260823")
    print("scope=THM3878_t>=U;necessary_filter_not_counterexample;LRC14=OPEN")
    print("theorem=counterexample_requires_t*lambda_plus<beta;closed_source_to_open_target;length_equality_closes")
    print("isolated_closed_equality_walls_have_lambda=0_and_do_not_trigger_filter")
    print("AP11_topology=14_positive_arcs+4_isolated=18_closed_components;strict_components=14")
    print("AP11_isolated=3/14[L5,R9],5/14[L3,R11],9/14[L11,R3],11/14[L9,R5]")
    print("AP11_J=[1/14,13/154];length=1/77;left_equality_owner=1;right_equality_owner=11;interior=strict")
    print("AP11_Ulambda=11/77=1/7;scale1_closed=57/57;beta_equal_rows=((1,3),(1,4),(1,9),(1,10))")
    print("scale2_pair=(1,9);components=(2/21,8/63),(55/63,19/21);beta2=2/63;AP11_closed=yes")
    print("generic_LRC13_Ulambda=1/42;scale1_closed=0/57")
    print("selected_staircase=1/35:15,1/28:34,1/21:46,1/14:52,1/7:57")
    print("full_beta_staircase=beta:multiplicity:cumulative=" + staircase_text(staircase))
    print("AP13_tight_control=6_isolated_unit_walls;each_has_opposing_single_owners")
    print("q25_coherent=no_witness_q<=25;first=2/27;response=181/2548=1/14-1/2548;owner=182")
    print("q25_incoherent=no_witness_q<=25;first=3/26;response=589/2744=3/14+1/2744;owner=196")
    print("semantic_sha256=" + digest)
    print(f"checks={CHECKS}")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
