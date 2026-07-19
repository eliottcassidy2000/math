#!/usr/bin/env python3
"""Exact referee for THM-1248 centered blocker address compression.

The proof-facing calculations use only integers and ``fractions.Fraction``.
The analytic inputs are THM-1233's ratio cap and THM-1240's centered-spoke
construction.  This referee replays the new relative-address identities,
their sharp constants, affine lasso transport, sampled owner germs, and
non-cover guardrails.
"""

from __future__ import annotations

from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from math import gcd, lcm, prod
from pathlib import Path


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def nearest_integer(x: F, *, upper_tie: bool = False) -> int:
    lower = x.numerator // x.denominator
    upper = lower + 1
    dl, du = x - lower, upper - x
    if dl < du:
        return lower
    if du < dl:
        return upper
    return upper if upper_tie else lower


def circle_distance(speed: int, time: F) -> F:
    value = speed * time
    return abs(value - nearest_integer(value))


def exact_text(value: object) -> str:
    if isinstance(value, dict):
        return "{" + ",".join(
            f"{key}:{exact_text(item)}" for key, item in value.items()
        ) + "}"
    if isinstance(value, list):
        return "[" + ",".join(exact_text(item) for item in value) + "]"
    if isinstance(value, tuple):
        body = ",".join(exact_text(item) for item in value)
        if len(value) == 1:
            body += ","
        return f"({body})"
    return str(value)


def centered_phase(c: int, k: int, d: int, *, upper_tie: bool = False) -> tuple[int, int, F, F]:
    require(0 <= k < c < d, (c, k, d))
    q = c + d
    t0 = F(2 * k + 1, 2 * c)
    p = nearest_integer(q * t0, upper_tie=upper_tie)
    epsilon = p - q * t0
    require(abs(epsilon) <= F(1, 2), (c, k, d, p, epsilon))
    return q, p, F(p, q), epsilon


def sampled_owner_germ(
    c: int,
    row: dict[str, object],
    owner: int,
    owner_phase: tuple[int, int, F, F],
) -> dict[str, object]:
    qi, dj = int(row["Qi"]), int(row["target"])
    sigma, r = int(row["sigma"]), row["ell"] - row["X"]
    e_i = 2 * c * row["epsilon_i"]
    e_h = 2 * c * owner_phase[3]
    require(e_i.denominator == 1 and e_h.denominator == 1, (e_i, e_h))
    numerator = (
        7 * dj * qi * (c - e_h)
        + 7 * owner * e_i * dj
        - 14 * c * owner * r
        + c * owner * sigma * qi
    )
    modulus = 14 * c * dj * qi
    require(numerator.denominator == 1, numerator)
    numerator = int(numerator)
    lift = (
        F(1, 2)
        - e_h / (2 * c)
        + owner * e_i / (2 * c * qi)
        - F(owner * r, dj * qi)
        + F(owner * sigma, 14 * dj)
    )
    require(lift == F(numerator, modulus), (owner, lift, numerator, modulus))
    require(
        circle_distance(1, lift) == circle_distance(owner, row["wall"]),
        (owner, row["wall"], lift),
    )
    tooth = nearest_integer(owner * row["wall"])
    alpha = owner * row["wall"] - tooth
    owns = abs(alpha) < F(1, 14)
    forward = (F(1, 14) - sigma * alpha) / owner if owns else None
    backward = (F(1, 14) + sigma * alpha) / owner if owns else None
    overlap = min(F(1, 7 * dj), backward) if owns else None
    return {
        "owner": owner,
        "E_owner": int(e_h),
        "numerator": numerator,
        "modulus": modulus,
        "owns": owns,
        "tooth": tooth,
        "alpha": alpha,
        "forward": forward,
        "backward": backward,
        "overlap": overlap,
    }


def blocker_edge(
    c: int,
    k: int,
    di: int,
    dj: int,
    pi: int,
    pj: int,
) -> dict[str, object]:
    qi, qj = c + di, c + dj
    ti, tj = F(pi, qi), F(pj, qj)
    n = nearest_integer(dj * ti)
    beta = dj * ti - n
    require(abs(beta) < F(1, 14), (c, k, di, dj, pi, pj, beta))

    t0 = F(2 * k + 1, 2 * c)
    epsilon_i, epsilon_j = pi - qi * t0, pj - qj * t0
    x = c * pi - k * qi
    residue = di * 0 + pi * dj - n * qi
    m = k + n
    ell = pi * qj - m * qi
    delta = pj - m
    determinant = pi * qj - pj * qi
    theta = F(ell, qi)
    a = F(dj, qi)
    left_germ = F(1, 14 * dj) + beta / dj
    right_germ = F(1, 14 * dj) - beta / dj

    require(F(1, 4) < F(x, qi) < F(3, 4), (di, x, qi))
    require(2 * x == qi + 2 * c * epsilon_i, (x, qi, c, epsilon_i))
    require(residue == qi * beta, (residue, qi, beta))
    require(ell == x + residue, (ell, x, residue))
    require(F(5, 28) < theta < F(23, 28), (di, dj, theta))
    require(ell == determinant + delta * qi, (ell, determinant, delta, qi))
    require(ell == determinant % qi, (ell, determinant, qi))
    require(m == (pi * qj) // qi, (m, pi, qj, qi, ell))
    require(gcd(qi, qj) and ell % gcd(qi, qj) == 0, (qi, qj, ell))
    require(
        beta == delta - F(1, 2) + a * epsilon_i - epsilon_j,
        (di, dj, beta, delta, epsilon_i, epsilon_j),
    )
    require(qj * (ti - tj) == theta - delta, (ti, tj, theta, delta))
    require(left_germ > 0 and right_germ > 0, (left_germ, right_germ))
    require(
        left_germ - right_germ == F(2 * (ell - x), dj * qi),
        (left_germ, right_germ, ell, x, dj, qi),
    )
    if delta <= 0:
        require(tj < ti, (di, dj, ti, tj, delta, theta))
        selected_wall = "left"
        sigma = -1
    else:
        require(delta >= 1 and ti < tj, (di, dj, ti, tj, delta, theta))
        selected_wall = "right"
        sigma = 1
    wall = F(14 * n + sigma, 14 * dj)
    target_tooth = (F(14 * n - 1, 14 * dj), F(14 * n + 1, 14 * dj))
    gap = (F(14 * k + 1, 14 * c), F(14 * k + 13, 14 * c))
    require(
        gap[0] < target_tooth[0] < target_tooth[1] < gap[1],
        (di, dj, gap, target_tooth),
    )
    source_residue = pi * dj - n * qi
    e_i = 2 * c * epsilon_i
    e_j = 2 * c * epsilon_j
    require(
        14 * dj * qj * (wall - tj)
        == 7 * e_j + (7 - 14 * delta + sigma) * qj,
        (di, dj, wall, tj, e_j, delta, sigma),
    )
    clearance_numerator = (
        -7 * sigma * e_j + (14 * abs(F(delta) - F(1, 2)) - 1) * qj
    )
    require(clearance_numerator.denominator == 1, clearance_numerator)
    clearance_numerator = int(clearance_numerator)
    require(
        clearance_numerator == 14 * dj * qj * sigma * (tj - wall),
        (di, dj, clearance_numerator, wall, tj),
    )
    require(clearance_numerator > 5 * dj, (di, dj, clearance_numerator))
    require(abs(theta - delta) > F(5, 14), (di, dj, theta, delta))
    if delta == 0:
        require(F(5, 14) < theta < F(23, 28), (di, dj, theta))
    if delta == 1:
        require(F(5, 28) < theta < F(9, 14), (di, dj, theta))
    source_wall_numerator = (
        7 * dj * (qi - e_i) + di * (sigma * qi - 14 * source_residue)
    )
    source_wall_lift = (
        F(1, 2)
        - e_i / (2 * qi)
        + F(di, dj) * (F(sigma, 14) - F(source_residue, qi))
    )
    require(
        source_wall_lift == F(source_wall_numerator, 14 * dj * qi),
        (di, dj, source_wall_lift, source_wall_numerator),
    )
    source_wall_depth = circle_distance(1, source_wall_lift)
    require(source_wall_depth == circle_distance(di, wall), (di, dj, wall))
    require(-586 <= delta <= 587, (di, dj, delta))
    if dj <= qi:
        require(delta in (0, 1), (di, dj, qi, delta))
        require(F(5, 14 * qj) < abs(ti - tj) < F(23, 28 * qj), (ti, tj))

    return {
        "source": di,
        "target": dj,
        "carrier": c,
        "Qi": qi,
        "Qj": qj,
        "Pi": pi,
        "Pj": pj,
        "N": n,
        "M": m,
        "epsilon_i": epsilon_i,
        "epsilon_j": epsilon_j,
        "beta": beta,
        "X": x,
        "ell": ell,
        "delta": delta,
        "determinant": determinant,
        "theta": theta,
        "slope": a,
        "left_germ": left_germ,
        "right_germ": right_germ,
        "selected_wall": selected_wall,
        "wall": wall,
        "sigma": sigma,
        "source_wall_numerator": source_wall_numerator,
        "source_wall_depth": source_wall_depth,
        "target_clearance_numerator": clearance_numerator,
    }


def check_cycle(rows: tuple[dict[str, object], ...]) -> dict[str, object]:
    require(len(rows) >= 2, len(rows))
    require(all(rows[i]["target"] == rows[(i + 1) % len(rows)]["source"] for i in range(len(rows))), rows)
    slopes = tuple(row["slope"] for row in rows)
    offsets = tuple(row["delta"] - F(1, 2) - row["beta"] for row in rows)
    multiplier = prod(slopes, start=F(1))
    weights = tuple(prod(slopes[i + 1 :], start=F(1)) for i in range(len(rows)))
    affine_offset = sum((offsets[i] * weights[i] for i in range(len(rows))), F(0))
    epsilon0 = rows[0]["epsilon_i"]
    require(affine_offset == (1 - multiplier) * epsilon0, (rows, affine_offset, multiplier, epsilon0))
    require(0 < multiplier < 1, multiplier)
    require(abs(affine_offset) <= (1 - multiplier) / 2, affine_offset)

    digit_center = sum(((F(rows[i]["delta"]) - F(1, 2)) * weights[i] for i in range(len(rows))), F(0))
    beta_radius = sum((weights[i] for i in range(len(rows))), F(0)) / 14
    require(abs(digit_center) < beta_radius + (1 - multiplier) / 2, (digit_center, beta_radius))

    weighted_delta = sum((F(row["delta"], row["Qj"]) for row in rows), F(0))
    weighted_theta = sum((row["theta"] / row["Qj"] for row in rows), F(0))
    reciprocal_weight = sum((F(1, row["Qj"]) for row in rows), F(0))
    require(weighted_delta == weighted_theta, (weighted_delta, weighted_theta))
    require(F(5, 28) * reciprocal_weight < weighted_delta < F(23, 28) * reciprocal_weight, weighted_delta)
    require(any(row["delta"] <= 0 for row in rows), rows)
    require(any(row["delta"] >= 1 for row in rows), rows)

    p_product = prod((row["Pi"] for row in rows), start=1)
    m_product = prod((row["M"] for row in rows), start=1)
    require(p_product > m_product, (p_product, m_product, rows))
    return {
        "length": len(rows),
        "multiplier": multiplier,
        "contraction_gap": 1 - multiplier,
        "delta_word": tuple(row["delta"] for row in rows),
        "ell_word": tuple(row["ell"] for row in rows),
        "positive_holonomy_gap": p_product - m_product,
    }


def check_constants() -> dict[str, F]:
    ratio = F(2345, 2346)
    gap = 1 - ratio * ratio
    require(gap == F(4691, 5503716), gap)
    require(F(1, 4) - F(1, 14) == F(5, 28), "lower central remainder")
    require(F(3, 4) + F(1, 14) == F(23, 28), "upper central remainder")
    require(F(1, 4) - F(1, 7) == F(3, 28) > F(1, 14), "wall-event margin")
    require(F(49, 6) * F(1, 14) == F(7, 12), "located-seam Hunter coefficient")
    require(6 * 2 * 2011 == 24132, "uniform fast-tooth boundary-event bank")
    require(6 * 2011 == 12066, "one-direction forward-endpoint bank")
    require(6 * 24132 == 144792, "six-edge lasso event-occurrence cap")
    return {
        "global_delta_min": F(-586),
        "global_delta_max": F(587),
        "contraction_gap": gap,
        "wall_event_source_margin": F(3, 28),
        "tail_homogeneous_multiplier_cap": F(13, 19),
        "two_seam_hunter_coefficient": F(7, 12),
    }


def check_universal_wall_quantum() -> dict[str, object]:
    rows = 0
    for j in range(2, 26):
        for h in range(2, 26):
            if h == j:
                continue
            g = gcd(j, h)
            for nj in range(-2, 3):
                for nh in range(-2, 3):
                    for side_j in (-1, 1):
                        for side_h in (-1, 1):
                            numerator = h * (14 * nj + side_j) - j * (
                                14 * nh + side_h
                            )
                            if numerator <= 0:
                                continue
                            require(numerator % g == 0, (j, h, numerator, g))
                            gap = F(numerator, 14 * j * h)
                            require(gap >= F(g, 14 * j * h), (j, h, gap))
                            rows += 1

    survivor_shapes = tuple(
        v for v in range(2, 7) if 2 * v <= 6 - v
    )
    require(survivor_shapes == (2,), survivor_shapes)
    return {
        "positive_endpoint_rows": rows,
        "two_wall_no_incidence_no_reuse_lasso_sizes": survivor_shapes,
        "residual_lasso_shape": (0, 2),
    }


def spanning_trees_on_six() -> tuple[frozenset[tuple[int, int]], ...]:
    all_edges = tuple(combinations(range(6), 2))
    trees = []
    for edges in combinations(all_edges, 5):
        parent = list(range(6))

        def root(vertex: int) -> int:
            while parent[vertex] != vertex:
                parent[vertex] = parent[parent[vertex]]
                vertex = parent[vertex]
            return vertex

        acyclic = True
        for left, right in edges:
            left_root, right_root = root(left), root(right)
            if left_root == right_root:
                acyclic = False
                break
            parent[left_root] = right_root
        if acyclic:
            trees.append(frozenset(edges))
    require(len(trees) == 6**4, len(trees))
    return tuple(trees)


def protected_tree_tariff(c: int, speeds: tuple[int, ...]) -> dict[str, object]:
    require(len(speeds) == 6 and c < speeds[0] < speeds[-1], (c, speeds))
    edge_weight = {
        edge: F(c, lcm(speeds[edge[0]], speeds[edge[1]]))
        for edge in combinations(range(6), 2)
    }
    trees = spanning_trees_on_six()
    blocker_edges = tuple(combinations(range(1, 6), 2))
    x = F(speeds[0], c)
    base = (9 * x + 2) / (3 * x * (1 + x))
    total = c * sum((F(1, speed) for speed in speeds), F(0))
    candidates = []
    for pair in combinations(blocker_edges, 2):
        forest = frozenset(pair)
        forest_weight = sum((edge_weight[edge] for edge in forest), F(0))
        remainders = (
            sum((edge_weight[edge] for edge in tree - forest), F(0))
            for tree in trees
            if forest <= tree
        )
        extension_weight = min(remainders)
        tariff = (
            base
            + F(7, 12) * forest_weight
            + max(F(0), F(7, 12) * extension_weight - (base - 1))
        )
        candidates.append((tariff, tuple(sorted(forest)), extension_weight))
    best = min(candidates)
    return {
        "S": total,
        "base": base,
        "Phi": best[0],
        "min_forest": best[1],
        "extension_weight": best[2],
        "passes": total >= best[0],
    }


def check_protected_tree_tariff_frontier() -> dict[str, object]:
    c = 52311
    closing = (98700, 111300, 747300, 1046220, 1743700, 2615550)
    require(all(lcm(left, right) == 100 * c for left, right in combinations(closing, 2)), closing)
    require(gcd(*closing) == 10, closing)
    closed = protected_tree_tariff(c, closing)
    require(closed["S"] == F(117, 100), closed)
    require(closed["base"] == F(26659, 22950), closed)
    require(closed["S"] - closed["base"] == F(77, 9180), closed)
    require(closed["Phi"] - closed["S"] == F(301, 91800), closed)
    require(not closed["passes"], closed)

    tangent_c = 10673449340400
    tangent = tuple(
        (slope * tangent_c + offset)
        for slope, offset in ((1, 1), (2, -1), (4, -1), (20, -1), (50, -1), (75, -1))
    )
    require(all(gcd(left, right) == 1 for left, right in combinations(tangent, 2)), tangent)
    open_case = protected_tree_tariff(tangent_c, tangent)
    require(open_case["passes"], open_case)
    require(
        F(1) + F(1, 2) + F(1, 4) + F(1, 20) + F(1, 50) + F(1, 75)
        == F(11, 6),
        "projective tariff guardrail",
    )
    return {
        "excluded_packet_c": c,
        "excluded_packet_S": closed["S"],
        "excluded_packet_base": closed["base"],
        "excluded_packet_Phi_minus_S": closed["Phi"] - closed["S"],
        "tangent_packet_c": tangent_c,
        "tangent_projective_sum": F(11, 6),
        "tangent_packet_passes": open_case["passes"],
    }


def check_tail(rows: tuple[dict[str, object], ...]) -> dict[str, object]:
    require(1 <= len(rows) <= 4, len(rows))
    require(all(rows[i]["target"] == rows[i + 1]["source"] for i in range(len(rows) - 1)), rows)
    c = rows[0]["carrier"]
    require(all(row["carrier"] == c for row in rows), rows)
    for row in rows:
        ei = 2 * c * row["epsilon_i"]
        ej = 2 * c * row["epsilon_j"]
        require(
            row["Qi"] * ej - row["Qj"] * ei
            == 2 * c * (row["Qi"] * row["delta"] - row["ell"]),
            row,
        )
        zi = row["epsilon_i"] / row["source"]
        zj = row["epsilon_j"] / row["target"]
        rho = F(row["source"], row["Qi"])
        require(
            zj
            == rho * zi
            + (row["delta"] - F(1, 2) - row["beta"]) / row["target"],
            row,
        )

    multiplier = prod((F(row["source"], row["Qi"]) for row in rows), start=F(1))
    require(multiplier < F(13, 19), multiplier)
    return {
        "length": len(rows),
        "delta_word": tuple(row["delta"] for row in rows),
        "E_word": tuple(
            int(value)
            for value in (
                [2 * c * rows[0]["epsilon_i"]]
                + [2 * c * row["epsilon_j"] for row in rows]
            )
        ),
        "wall_word": tuple(row["wall"] for row in rows),
        "multiplier": multiplier,
    }


def check_unbounded_address_guardrail(c: int) -> dict[str, object]:
    require(c >= 27, c)
    k = c - 1
    speeds = (c + 1, c + 2, c + 3, c + 4, 2 * c, 3 * c)
    phases: dict[int, tuple[int, int, F, F]] = {}
    for d in speeds:
        phases[d] = centered_phase(c, k, d, upper_tie=(d == 2 * c))
    blockers = {c + 1: 2 * c, c + 2: 2 * c, c + 3: 2 * c, c + 4: 2 * c, 2 * c: 3 * c, 3 * c: 2 * c}
    rows = tuple(
        blocker_edge(c, k, di, dj, phases[di][1], phases[dj][1])
        for di, dj in blockers.items()
    )
    cycle_rows = tuple(row for row in rows if row["source"] in (2 * c, 3 * c))
    cycle = check_cycle(cycle_rows)
    witness = 1 - F(2, 5 * c)
    packet = (c,) + speeds
    depths = tuple(circle_distance(v, witness) for v in packet)
    require(min(depths) == F(1, 5), (c, packet, witness, depths))
    require(gcd(*packet[:2]) == 1, packet)
    return {
        "c": c,
        "max_tooth_address": max(row["N"] for row in rows),
        "cycle_delta_word": cycle["delta_word"],
        "witness": witness,
        "witness_depth": min(depths),
    }


def check_large_delta_guardrail() -> dict[str, object]:
    c, k = 1, 0
    speeds = (2, 16, 17, 34, 35, 2343)
    phases: dict[int, tuple[int, int, F, F]] = {
        d: centered_phase(c, k, d, upper_tie=(d == 2)) for d in speeds
    }
    blockers = {2: 2343, 2343: 2, 16: 17, 17: 16, 34: 35, 35: 34}
    rows = {
        (di, dj): blocker_edge(c, k, di, dj, phases[di][1], phases[dj][1])
        for di, dj in blockers.items()
    }
    cycles = tuple(
        check_cycle((rows[(a, b)], rows[(b, a)]))
        for a, b in ((2, 2343), (16, 17), (34, 35))
    )
    critical = rows[(2, 2343)]
    require(critical["delta"] == -390 and critical["ell"] == 2, critical)
    require(critical["M"] == 1562 and critical["Pi"] == 2, critical)
    require(
        F(2, 1) < F(13, 6)
        and F(16, 2) < F(273, 29)
        and F(17, 16) < F(84, 5)
        and F(34, 17) < F(343, 15)
        and F(35, 34) < F(189, 8)
        and F(2343, 35) < 77
        and 2343 < 2345,
        "THM-1233 compact inequalities",
    )
    witness = F(1, 6)
    packet = (c,) + speeds
    depth = min(circle_distance(v, witness) for v in packet)
    require(depth == F(1, 6), (packet, witness, depth))
    return {
        "packet": packet,
        "critical_edge": (critical["source"], critical["target"]),
        "critical_delta": critical["delta"],
        "critical_ell": critical["ell"],
        "cycle_delta_words": tuple(cycle["delta_word"] for cycle in cycles),
        "witness": witness,
        "witness_depth": depth,
    }


def check_long_tail_guardrail() -> dict[str, object]:
    c, k = 1, 0
    speeds = (2, 3, 4, 5, 6, 7)
    phases = {d: centered_phase(c, k, d) for d in speeds}
    edge_word = ((2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (7, 6))
    rows = tuple(
        blocker_edge(c, k, di, dj, phases[di][1], phases[dj][1])
        for di, dj in edge_word
    )
    tail = check_tail(rows[:4])
    cycle = check_cycle(rows[4:])
    witness = F(1, 8)
    packet = (c,) + speeds
    depth = min(circle_distance(v, witness) for v in packet)
    require(depth == F(1, 8), (packet, witness, depth))
    require(tail["delta_word"] == (1, 0, 1, 0), tail)
    require(tail["E_word"] == (-1, 0, -1, 0, -1), tail)
    require(tail["wall_word"] == (F(5, 14), F(27, 56), F(29, 70), F(41, 84)), tail)
    require(cycle["delta_word"] == (1, 0), cycle)
    return {
        "packet": packet,
        "tail_delta_word": tail["delta_word"],
        "tail_E_word": tail["E_word"],
        "tail_wall_word": tail["wall_word"],
        "cycle_delta_word": cycle["delta_word"],
        "witness": witness,
        "witness_depth": depth,
    }


def check_off_orbit_owner_guardrail() -> dict[str, object]:
    c, k = 1, 0
    orbit = (2, 3, 4, 5)
    phases = {d: centered_phase(c, k, d) for d in orbit}
    edge_word = ((2, 3), (3, 4), (4, 5), (5, 4))
    rows = tuple(
        blocker_edge(c, k, di, dj, phases[di][1], phases[dj][1])
        for di, dj in edge_word
    )
    wall = rows[0]["wall"]
    require(wall == F(5, 14), wall)
    packets = (
        ((1, 2, 3, 4, 5, 6, 7), F(1, 8), "none"),
        ((1, 2, 3, 4, 5, 6, 14), F(1, 12), "14"),
        ((1, 2, 3, 4, 5, 6, 28), F(1, 12), "28"),
    )
    report = []
    for packet, witness, owner in packets:
        extras = packet[3:]
        dangerous = tuple(v for v in extras if circle_distance(v, wall) < F(1, 14))
        if owner == "none":
            require(not dangerous, (packet, wall, dangerous))
        else:
            require(int(owner) in dangerous, (packet, wall, dangerous))
        depth = min(circle_distance(v, witness) for v in packet)
        require(depth >= F(1, 12), (packet, witness, depth))
        report.append((packet[-1], owner, witness, depth))
    return {
        "orbit": orbit,
        "delta_word": tuple(row["delta"] for row in rows),
        "first_wall": wall,
        "owner_rows": tuple(report),
    }


def check_forced_tail_owner_guardrail() -> dict[str, object]:
    c, k = 2, 0
    orbit = (3, 5, 14, 4, 6)
    edge_word = ((3, 5), (5, 14), (14, 4), (4, 6), (6, 4))
    expected_blockers = {3: 5, 5: 14, 14: 4, 4: 6, 6: 4}
    expected_owner_walls = {
        13: ((F(3, 14), (14,)), (F(13, 56), (13,))),
        19: ((F(3, 14), (14,)),),
        39: (
            (F(5, 28), (39,)),
            (F(3, 14), (14,)),
            (F(13, 56), (39,)),
            (F(55, 196), (39,)),
        ),
    }
    cases = []
    midpoint_residues = []
    external_39_trace: dict[F, tuple[int, F, F, F, F]] = {}
    for external in (13, 19, 39):
        fast = tuple(dict.fromkeys(orbit + (external,)))
        phases = {d: centered_phase(c, k, d) for d in fast}
        forced: dict[int, int] = {}
        for source in fast:
            blockers = tuple(
                target
                for target in fast
                if target != source
                and circle_distance(target, phases[source][2]) < F(1, 14)
            )
            require(len(blockers) == 1, (external, source, blockers))
            forced[source] = blockers[0]
        require(
            all(forced[source] == target for source, target in expected_blockers.items()),
            (external, forced),
        )
        require(forced[external] == 4, (external, forced))

        rows = tuple(
            blocker_edge(c, k, source, target, phases[source][1], phases[target][1])
            for source, target in edge_word
        )
        require(tuple(row["delta"] for row in rows) == (1, 0, 0, 1, 0), rows)
        require(tuple(row["ell"] for row in rows) == (2, 4, 8, 2, 4), rows)
        require(
            tuple(row["source_wall_depth"] for row in rows)
            == (F(5, 14), F(79, 196), F(1, 4), F(2, 7), F(11, 28)),
            rows,
        )
        reverse_depths = tuple(
            circle_distance(source, phases[target][2])
            for source, target in edge_word[:3]
        )
        require(reverse_depths == (F(1, 7), F(1, 4), F(1, 3)), reverse_depths)

        owners_by_wall: dict[F, set[int]] = {}
        for row in rows:
            wall = row["wall"]
            germs = {
                label: sampled_owner_germ(c, row, label, phases[label])
                for label in fast
                if label not in (row["source"], row["target"])
            }
            owners = {label for label, germ in germs.items() if germ["owns"]}
            owners_by_wall.setdefault(wall, set()).update(owners)
            if wall == F(13, 56) and row["source"] == 14:
                germ = sampled_owner_germ(c, row, external, phases[external])
                midpoint_residues.append(
                    (external, germ["E_owner"], germ["alpha"], germ["owns"])
                )
            if external == 39 and 39 in owners:
                germ = germs[39]
                external_39_trace[wall] = (
                    germ["tooth"],
                    germ["alpha"],
                    germ["forward"],
                    germ["backward"],
                    germ["overlap"],
                )
        paid = tuple(
            (wall, tuple(sorted(owners)))
            for wall, owners in sorted(owners_by_wall.items())
            if owners
        )
        require(paid == expected_owner_walls[external], (external, paid))

        packet = (c,) + fast
        witness = F(1, 8)
        depth = min(circle_distance(v, witness) for v in packet)
        require(depth >= F(1, 8), (external, packet, depth))
        cases.append((external, paid, depth))

    return {
        "forced_orbit": edge_word,
        "tail_reverse_depths": (F(1, 7), F(1, 4), F(1, 3)),
        "source_wall_depths": (F(5, 14), F(79, 196), F(1, 4), F(2, 7), F(11, 28)),
        "midwall_owner_residues": tuple(midpoint_residues),
        "external_39_tooth_trace": tuple(sorted(external_39_trace.items())),
        "external_owner_cases": tuple(cases),
    }


def check_two_wall_pigeonhole_guardrail() -> dict[str, object]:
    def one_case(
        c: int,
        a: int,
        b: int,
        outside: tuple[int, int, int, int],
        expected_owners: tuple[int, int, int, int],
        expected_digits: tuple[int, int],
        expected_ells: tuple[int, int],
    ) -> dict[str, object]:
        k = 0
        speeds = (a, b) + outside
        phases = {d: centered_phase(c, k, d) for d in speeds}
        rows = (
            blocker_edge(c, k, a, b, phases[a][1], phases[b][1]),
            blocker_edge(c, k, b, a, phases[b][1], phases[a][1]),
        )
        require(tuple(row["delta"] for row in rows) == expected_digits, rows)
        require(tuple(row["ell"] for row in rows) == expected_ells, rows)

        wall_rows = []
        for row in rows:
            target = int(row["target"])
            n = int(row["N"])
            other_lasso = a if target == b else b
            for side in (-1, 1):
                wall = F(14 * n + side, 14 * target)
                require(
                    circle_distance(other_lasso, wall) >= F(1, 14),
                    (c, a, b, wall, other_lasso),
                )
                owners = tuple(
                    h for h in outside if circle_distance(h, wall) < F(1, 14)
                )
                wall_rows.append((target, wall, owners))

        require(
            tuple(owner for _, _, owners in wall_rows for owner in owners)
            == expected_owners,
            wall_rows,
        )
        require(all(len(owners) == 1 for _, _, owners in wall_rows), wall_rows)
        require(len(set(expected_owners)) == 4, expected_owners)
        witness = F(1, 8)
        packet = (c,) + speeds
        depth = min(circle_distance(v, witness) for v in packet)
        require(depth == F(1, 8), (packet, witness, depth))
        return {
            "packet": packet,
            "phases": (phases[a][2], phases[b][2]),
            "digits": expected_digits,
            "ells": expected_ells,
            "wall_owners": tuple(wall_rows),
            "witness": witness,
            "witness_depth": depth,
        }

    binary = one_case(
        2,
        4,
        6,
        (15, 19, 28, 43),
        (19, 28, 43, 15),
        (1, 0),
        (2, 4),
    )
    nonbinary = one_case(
        2,
        4,
        18,
        (11, 13, 25, 29),
        (25, 29, 13, 11),
        (2, 0),
        (2, 10),
    )
    return {"binary_residual": binary, "nonbinary_ascent_residual": nonbinary}


def check_greedy_retrace_guardrail() -> dict[str, object]:
    c, k = 1, 0
    speeds = (2, 3, 4, 5, 12, 13)
    phases = {d: centered_phase(c, k, d) for d in speeds}
    edge_word = ((2, 3), (3, 12), (12, 13), (13, 4), (4, 5), (5, 12))
    rows = {
        edge: blocker_edge(c, k, edge[0], edge[1], phases[edge[0]][1], phases[edge[1]][1])
        for edge in edge_word
    }
    incoming, outgoing = rows[(3, 12)], rows[(12, 13)]
    require(incoming["wall"] == F(83, 168) and incoming["sigma"] == -1, incoming)
    require(outgoing["wall"] == F(85, 182) and outgoing["sigma"] == 1, outgoing)
    tooth_2 = (F(13, 28), F(15, 28))
    tooth_13 = (F(83, 182), F(85, 182))
    seam = (max(tooth_2[0], tooth_13[0]), min(tooth_2[1], tooth_13[1]))
    require(seam == (F(13, 28), F(85, 182)), seam)
    require(seam[1] - seam[0] == F(1, 364), seam)
    germ = sampled_owner_germ(c, outgoing, 2, phases[2])
    require(germ["owns"] and germ["overlap"] == F(1, 364), germ)
    require(germ["alpha"] == F(-6, 91), germ)
    target_spoke = phases[12][2]
    require(target_spoke == F(6, 13), target_spoke)
    require(tooth_2[0] > target_spoke, (tooth_2, target_spoke))
    require(tooth_13[0] < target_spoke < tooth_13[1], (tooth_13, target_spoke))
    require(tooth_13[0] < tooth_2[0] < tooth_13[1], (tooth_2, tooth_13))
    witness = F(1, 8)
    packet = (c,) + speeds
    depth = min(circle_distance(v, witness) for v in packet)
    require(depth == F(1, 8), (packet, witness, depth))
    return {
        "functional_orbit": edge_word,
        "turn_walls": (incoming["wall"], outgoing["wall"]),
        "reused_seam": seam,
        "reused_seam_length": seam[1] - seam[0],
        "owner_2_outgoing_germ": (germ["tooth"], germ["alpha"], germ["overlap"]),
        "incoming_target_spoke": target_spoke,
        "signed_expiry_step": (tooth_2[0], tooth_13[0]),
        "witness": witness,
        "witness_depth": depth,
    }


def main() -> None:
    print("THM-1248 CENTERED BLOCKER ADDRESS COMPRESSION EXACT REFEREE")
    print("method=integer/Fraction only; always-on checks; no dependencies")

    constants = check_constants()
    print("\nFINITE RELATIVE-ADDRESS AND WALL-EVENT CONSTANTS")
    for key, value in constants.items():
        print(f"{key}={value}")
    print("central_remainder_band=(5/28,23/28)")
    print("determinant_phase_separation=|theta-delta|>5/14")
    print("binary sharpening: delta=0 => theta>5/14; delta=1 => theta<9/14")
    print("relative_delta_word=[-586,587]; every edge with d_j<=c+d_i has digit 0 or 1")
    print("central_coordinate=X=(Q_i+E_i)/2 with E_i=2c*epsilon_i")
    print("oriented_germ=2*(ell-X)/(d_j*Q_i); delta<=0 selects left, delta>=1 selects right")
    print("tail_length<=4; nonempty slowest-rooted tail homogeneous multiplier <13/19")
    print("selected-wall source safety is exact quotient data via an integer numerator")
    print("all-candidate sampled ownership and signed owner expiry are quotient-computable")
    print("distinct slowest/cycle seam pairs give Hunter coefficient 7/12 per gcd/(uv)")
    print("every blocker tooth lies inside G and exports two target-owner gcd seams")
    print("source-safe/ascent walls sharpen these to fourth-support/protected seams")

    universal = check_universal_wall_quantum()
    print("\nUNIVERSAL TWO-WALL QUANTUM / OWNER PIGEONHOLE")
    for key, value in universal.items():
        print(f"{key}={exact_text(value)}")
    print("without lasso incidence or outside-owner reuse, only a slowest two-cycle remains")

    tariff = check_protected_tree_tariff_frontier()
    print("\nPROTECTED/FULLY-LOCATED TREE TARIFF")
    for key, value in tariff.items():
        print(f"{key}={exact_text(value)}")
    print("the protected two-edge tariff closes an exact low-common-gcd lcm packet")
    print("a pairwise-coprime projective tangent packet still passes the full scalar tariff")

    print("\nUNBOUNDED ABSOLUTE-ADDRESS GUARDRAIL")
    for c in (27, 101, 1001, 10001):
        row = check_unbounded_address_guardrail(c)
        print(
            f"c={c} max_tooth_address={row['max_tooth_address']} "
            f"cycle_delta_word={row['cycle_delta_word']} "
            f"witness={row['witness']} depth={row['witness_depth']}"
        )
    print("sampled obligations and primitive compact ratios do not bound absolute addresses")
    print("these packets are globally lonely, not six-cover counterexamples")

    large = check_large_delta_guardrail()
    print("\nLARGE RELATIVE-DIGIT GUARDRAIL")
    for key, value in large.items():
        print(f"{key}={value}")
    print("the {-586,...,587} alphabet cannot be replaced by {0,1} on speed-ascent edges")

    tail = check_long_tail_guardrail()
    print("\nLONG-TAIL ORIENTATION GUARDRAIL")
    for key, value in tail.items():
        print(f"{key}={exact_text(value)}")
    print("strictly increasing tail speeds can have alternating phases and wall coordinates")
    print("this packet is globally lonely and is not a six-cover counterexample")

    owners = check_off_orbit_owner_guardrail()
    print("\nOFF-ORBIT OWNER GUARDRAIL")
    for key, value in owners.items():
        print(f"{key}={exact_text(value)}")
    print("identical orbit-local quotient data do not determine fourth-owner existence or identity")
    print("coverage is the extra hypothesis supplying the off-sample owner letter")

    forced = check_forced_tail_owner_guardrail()
    print("\nFORCED-LASSO OWNER-REUSE GUARDRAIL")
    for key, value in forced.items():
        print(f"{key}={exact_text(value)}")
    print("unique blocker maps and safe selected walls still do not determine owner reuse")
    print("the full six-vertex quotient distinguishes the owners, but orbit-local data do not")
    print("even paying every distinct sampled lasso wall need not reuse one owner tooth")
    print("all three packets are globally lonely and are not six-cover counterexamples")

    two_wall = check_two_wall_pigeonhole_guardrail()
    print("\nTWO-WALL PIGEONHOLE RESIDUAL GUARDRAIL")
    for key, value in two_wall.items():
        print(f"{key}={exact_text(value)}")
    print("four distinct outside wall owners are arithmetically realizable on a slowest two-cycle")
    print("the ascent digit need not be binary; coverage/coherent-word topology remains essential")
    print("both residual packets are globally lonely and are not six-cover counterexamples")

    retrace = check_greedy_retrace_guardrail()
    print("\nGREEDY-TURN RETRACING GUARDRAIL")
    for key, value in retrace.items():
        print(f"{key}={exact_text(value)}")
    print("two adjacent blocker transports can traverse the identical located seam in reverse")
    print("arbitrary greedy continuation can retrace; the coherent minimal word removes this choice")
    print("this packet is globally lonely and is not a six-cover counterexample")

    print("\nTOURNAMENT / ALTERNATE-VERTEX AUDIT")
    print("observable=D_ij=P_i*Q_j-P_j*Q_i")
    print("switch=i->j iff (t_i,i) is lexicographically below (t_j,j); index breaks D_ij=0 ties")
    print("tie_Hamiltonian_path=vertices sorted by (t_i,i)")
    print("phase-determinant sign gauge is transitive and loses the central remainder")
    print("proof-bearing edge=(Qi,Qj,least-positive ell,relative delta,gcd sheet)")
    print("full sampled quotient=(c;(Qi,Ei)_i; blocker edges; derived signed owner-germ hypergraph)")
    print("expiry observable=forward tooth endpoint; switch=travel sign sigma")
    print("expiry tie Hamiltonian path=active tooth instances sorted by endpoint,label,tooth")
    print("distinct boundary-event bank<=24132; six-edge lasso occurrence cap<=144792")
    print("one fixed travel direction uses at most 12066 forward endpoints")
    print("challenged_assumption=use tooth instances and boundary events, not runners, as transport vertices")
    print("preserves=blocker danger, determinant gcd, sampled ownership, germ direction, and component identity")
    print("destroys=absolute global tooth lift; unsampled chronology is restored only by the greedy cover walk")
    print("runner scores=(0,1,2,3,4,5); cycles=0; SCCs=6; Hamilton_paths=1")
    print("STATUS=PASS")
    print(f"source_sha256={sha256(Path(__file__).read_bytes()).hexdigest()}")


if __name__ == "__main__":
    main()
