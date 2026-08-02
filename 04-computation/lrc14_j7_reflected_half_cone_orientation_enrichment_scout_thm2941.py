#!/usr/bin/env python3
r"""Exact physical-orientation repair for the reflected cap-3 half-cone.

The old half-cone referee retained correct conditional overlap calculations,
but its 26 unconstrained-body lane atlas did not cover physical assignments.
This scout isolates that failure and tests one principled enrichment: couple
each ratio lane to the lower bound on the *maximum* selected-pair level that
the same ratio interval forces.

The referee does not inherit the quarantined half-cone conclusion.  It
recomputes the coarse CSP partition, audits both physical orientations of all
253 constrained policies, and then replaces the faulty 26-body tail atlas by
an assignment-complete 35-lane atlas with 1,600 reduced head controls.
"""

from __future__ import annotations

import argparse
import hashlib
from fractions import Fraction as F
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, permutations
from math import ceil, gcd, lcm
from pathlib import Path


HERE = Path(__file__).resolve()
ROOT = next(parent for parent in HERE.parents if (parent / "00-navigation").is_dir())
OLD = ROOT / "04-computation/lrc14_j7_reflected_half_cone_closure_thm2941.py"
OLD_OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_reflected_half_cone_closure_thm2941.out"
REPAIR = ROOT / "04-computation/lrc14_j7_reflected_two_thirds_tail_orientation_repair_thm2941.py"
REPAIR_OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_reflected_two_thirds_tail_orientation_repair_thm2941.out"
OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_reflected_half_cone_orientation_enrichment_scout_thm2941.out"

EXPECTED_OLD_SHA256 = "01b0ef11060585cfdd8dcf001cf0488ffde69118cb4da1bbbd964b195e3cd910"
EXPECTED_OLD_OUTPUT_SHA256 = "afe07bc18ae2fa257cca7eda51999e064411bc168a409c9d932b951e87a4d7e3"
EXPECTED_OLD_SEMANTIC = "e7edbeec8933409bb36eab7d4a1867a10b3f7683a4898ec30f244a9e81601673"
EXPECTED_REPAIR_SHA256 = "420d21d0c65cfbd57d8cc926c04993b14e61d2554fb291a836af949d1fd664ed"
EXPECTED_REPAIR_OUTPUT_SHA256 = "87db762ba3ff6794748a2a19a4bf192a8cef8d6cebe7c619c871579da669864f"
EXPECTED_REPAIR_SEMANTIC = "525c8d2d89cf11c9fcc6260dbf81e57d0e47c999049c89b9148ae731dbac4730"
EXPECTED_TRAP_DIGEST = "35fceb70d1e960ee5e0590c990adde5c4848b657a15455714585354bee322380"
EXPECTED_VERDICT_DIGEST = "fe8193713f596b6f3f0182066a8937cfe75f681ea760fd31bac4642b7eb0f7cc"
EXPECTED_LOCATED_DIGEST = "8968666298108aabbf926de452354ed554728aa2fe2cc2ed5f5aaea3598f4201"
EXPECTED_HEAD_DIGEST = "81dc1b57dc7cf3dfe2dc92dfeadf6af375afacee40486ebc6a4b43efa51eb443"
EXPECTED_SEMANTIC = "1a8a6cac81420388034c849ca0ba852a9b87416c85eeedcb12e1cbafea572f5a"

CAP = F(3)
MIN_LEVEL = 3
MIN_SPREAD = 6
ETA_GUARD = F(35, 164)
COARSE_CLOSED = 282
COARSE_TRAPS = 279
H = (1, 2, 3, 4, 6, 12)
H2 = (1, 3, 4, 6, 8, 12)
H3 = (2, 3, 4, 6, 8, 12)

# Each row is (left, right, ordered pair, first tail scale, numerator).
# H and H2 keep the same physical pair throughout; adjacent closed endpoints
# overlap harmlessly.  H3 uses two inverse low charts, whose failure edges
# form a directed 2-cycle.
H_ATLAS = (
    (F(1, 3), F(3, 8), (0, 1), 56, F(10, 3)),
    (F(3, 8), F(3, 7), (0, 1), 35, F(13, 4)),
    (F(3, 7), F(1, 2), (0, 1), 26, F(22, 7)),
    (F(1, 2), F(2, 3), (0, 1), 23, F(3)),
    (F(2, 3), F(1), (0, 1), 19, F(8, 3)),
    (F(1), F(3), (0, 1), 10, F(2)),
)
H2_ATLAS = (
    (F(1, 3), F(2, 5), (0, 1), 10, F(16, 3)),
    (F(2, 5), F(3, 5), (0, 1), 10, F(26, 5)),
    (F(3, 5), F(1), (0, 1), 10, F(24, 5)),
    (F(1), F(3), (0, 1), 9, F(4)),
)
H3_ATLAS = (
    (F(1, 3), F(1), (0, 1), 10, F(14, 3)),
    (F(1, 3), F(1), (1, 0), 6, F(2)),
)

H_STORED_HOSTILE = (3, 4, 5, 8, 6, 9)
H2_STORED_HOSTILE = (5, 4, 3, 6, 7, 9)


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def qtext(value: F) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def import_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, ("cannot import", path))
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


require(sha256(OLD) == EXPECTED_OLD_SHA256, ("old source changed", sha256(OLD)))
require(sha256(OLD_OUTPUT) == EXPECTED_OLD_OUTPUT_SHA256,
        ("old output changed", sha256(OLD_OUTPUT)))
require(f"semantic_sha256={EXPECTED_OLD_SEMANTIC}" in OLD_OUTPUT.read_text(),
        "old semantic token missing")
require(sha256(REPAIR) == EXPECTED_REPAIR_SHA256,
        ("two-thirds repair source changed", sha256(REPAIR)))
require(sha256(REPAIR_OUTPUT) == EXPECTED_REPAIR_OUTPUT_SHA256,
        ("two-thirds repair output changed", sha256(REPAIR_OUTPUT)))
require(f"semantic_sha256={EXPECTED_REPAIR_SEMANTIC}" in REPAIR_OUTPUT.read_text(),
        "two-thirds repair semantic token missing")
OLD_MODULE = import_module("half_cone_orientation_scout_old", OLD)
T = OLD_MODULE.T


def singleton_debt(body, ruler: int, level: int) -> F:
    return T.C2.singleton_debt(body, ruler, level)


def transport_bound(body, ruler: int, pair) -> F:
    """Exact cap-3 maximum, attained at the corner (D,m)=(6,3)."""
    a, b = body[pair[0]], body[pair[1]]
    return F(3 * (b - a) + 6 * b, 3 * ruler - b)


def coarse_threshold(body, ruler: int, pair) -> F:
    return 2 * transport_bound(body, ruler, pair) + T.debt_at_nine(body, ruler)


def body_constraints(body, ruler: int, primitive_rows) -> tuple[dict, dict]:
    """Rebuild the complete unordered-ratio CSP for one residual body."""
    constraints = {}
    thresholds = {}
    for pair in combinations(range(6), 2):
        threshold = coarse_threshold(body, ruler, pair)
        if threshold >= T.PHASE_LIMIT:
            continue
        bound = T.product_bound(threshold)
        allowed = frozenset(
            ratio for ratio, phase, product, _, _ in primitive_rows
            if product <= bound and phase <= threshold
        )
        constraints[pair] = allowed
        thresholds[pair] = threshold
    return constraints, thresholds


def lane_alpha(left: F, right: F) -> F:
    """Forced pair-max / reduced-min factor on one ratio interval."""
    require(F(1, 3) <= left < right <= 3, (left, right))
    if right <= 1:
        return 1 / right
    if left >= 1:
        return left
    return F(1)


def debt_level_at_scale(scale: int, alpha: F) -> int:
    return max(MIN_LEVEL, ceil(F(alpha * scale, CAP)))


def lane_envelope(body, pair, left: F, right: F, scale: int,
                  numerator: F) -> tuple[F, int, F]:
    ruler = 14 * lcm(*body)
    a = body[pair[0]]
    level = debt_level_at_scale(scale, lane_alpha(left, right))
    debt = singleton_debt(body, ruler, level)
    value = (F(1, 49) - F(12, 49 * scale**2)
             - numerator / (ruler - F(a, scale)) - debt)
    return value, level, debt


def debt_drop(body, ruler: int, first: int, second: int) -> F:
    require(MIN_LEVEL <= first <= second, (body, first, second))
    value = singleton_debt(body, ruler, first) - singleton_debt(body, ruler, second)
    require((value == 0) if first == second else (value > 0),
            (body, first, second, value))
    return value


def analytic_monotonicity(body, pair, left: F, right: F, scale: int,
                          numerator: F) -> tuple:
    """Exact identities which remain positive at every later scale."""
    ruler = 14 * lcm(*body)
    a = body[pair[0]]
    alpha = lane_alpha(left, right)
    next_scale = scale + 1
    first_level = debt_level_at_scale(scale, alpha)
    second_level = debt_level_at_scale(next_scale, alpha)
    phase_gain = F(12, 49) * (F(1, scale**2) - F(1, next_scale**2))
    transport_gain = numerator * (
        F(scale, ruler * scale - a)
        - F(next_scale, ruler * next_scale - a)
    )
    require(phase_gain > 0, (body, scale, phase_gain))
    require(transport_gain == F(numerator * a,
            (ruler * scale - a) * (ruler * next_scale - a)) > 0,
            (body, pair, scale, transport_gain))
    debt_gain = debt_drop(body, ruler, first_level, second_level)
    current = lane_envelope(body, pair, left, right, scale, numerator)[0]
    later = lane_envelope(body, pair, left, right, next_scale, numerator)[0]
    require(later - current == phase_gain + transport_gain + debt_gain > 0,
            (body, pair, left, right, scale, current, later))
    return (body, pair, left, right, scale, alpha, first_level, second_level,
            phase_gain, transport_gain, debt_gain, later - current)


def interval_channels(left: F, right: F, threshold: int, oversize: int = 3):
    return tuple(
        (P, Q)
        for P in range(1, oversize * threshold + 1)
        for Q in range(1, oversize * threshold + 1)
        if P != Q and gcd(P, Q) == 1 and min(P, Q) < threshold
        and left <= F(Q, P) <= right
    )


def scale_monotonicity(body, pair, channel, scale: int) -> tuple:
    """Uniform suffix certificate for one reduced head channel."""
    ruler = 14 * lcm(*body)
    a, b = body[pair[0]], body[pair[1]]
    P, Q = channel
    mismatch = abs(Q * a - P * b)
    eta_now = F(scale * mismatch, P * scale * ruler - a)
    eta_next = F((scale + 1) * mismatch, P * (scale + 1) * ruler - a)
    phase_drop = F(mismatch * a,
                   (P * scale * ruler - a)
                   * (P * (scale + 1) * ruler - a))
    require(eta_now - eta_next == phase_drop,
            (body, pair, channel, scale, eta_now, eta_next, phase_drop))
    require((phase_drop == 0) if mismatch == 0 else (phase_drop > 0),
            (body, pair, channel, mismatch, phase_drop))
    transport_now = F(a, P * scale * ruler - a)
    transport_next = F(a, P * (scale + 1) * ruler - a)
    transport_drop = F(a * P * ruler,
                       (P * scale * ruler - a)
                       * (P * (scale + 1) * ruler - a))
    require(transport_now - transport_next == transport_drop > 0,
            (body, pair, channel, scale, transport_drop))
    first_level = debt_level_at_scale(scale, F(max(P, Q)))
    second_level = debt_level_at_scale(scale + 1, F(max(P, Q)))
    singleton_drop = debt_drop(body, ruler, first_level, second_level)
    require(min(P, Q) * scale >= MIN_LEVEL and eta_now <= ETA_GUARD < 1,
            (body, pair, channel, scale, eta_now))
    return (body, pair, channel, scale, mismatch, phase_drop, transport_drop,
            first_level, second_level, singleton_drop)


def multi_scale_head(body, pair, channel) -> tuple:
    """Direct finite prefix followed by an algebraically uniform tail."""
    ruler, safe_ranges = T.R.safe_cell_ranges(body)
    a, b = body[pair[0]], body[pair[1]]
    P, Q = channel
    g0 = ceil(F(MIN_LEVEL, min(P, Q)))
    skeleton, skeleton_cell = max(
        (T.intersection_mass(
            T.primitive_arcs(P, F(a * cell, ruler) % 1),
            T.primitive_arcs(Q, F(b * cell, ruler) % 1),
         ), cell)
        for left, right in safe_ranges for cell in range(left, right)
    )

    tail = None
    for g in range(g0, 1001):
        level = debt_level_at_scale(g, F(max(P, Q)))
        debt = singleton_debt(body, ruler, level)
        eta = F(g * (Q * a - P * b), P * g * ruler - a)
        bracket = skeleton - 2 * abs(eta)
        if bracket <= debt:
            continue
        c = 1 - F(a, P * g * ruler)
        transported = bracket / c
        actual = T.intersection_mass(
            T.R.reflected_level_arcs(ruler, a, P * g, skeleton_cell),
            T.R.reflected_level_arcs(ruler, b, Q * g, skeleton_cell),
        )
        require(actual >= transported > bracket > debt and 0 < c < 1,
                ("tail transport", body, pair, channel, g, actual,
                 transported, bracket, debt))
        monotonicity = scale_monotonicity(body, pair, channel, g)
        tail = (bracket - debt, g, level, skeleton_cell, skeleton, eta,
                debt, actual, transported, monotonicity)
        break
    require(tail is not None, ("no transported tail", body, pair, channel))

    direct_rows = []
    for g in range(g0, tail[1]):
        level = debt_level_at_scale(g, F(max(P, Q)))
        debt = singleton_debt(body, ruler, level)
        best = max(
            (T.intersection_mass(
                T.R.reflected_level_arcs(ruler, a, P * g, cell),
                T.R.reflected_level_arcs(ruler, b, Q * g, cell),
             ) - debt, cell)
            for left, right in safe_ranges for cell in range(left, right)
        )
        require(best[0] > 0, ("direct head", body, pair, channel, g, level, best))
        direct_rows.append((g, level, best))
    return (tail[0], body, pair, channel, g0, tail, tuple(direct_rows))


def lane_applies(left: F, right: F, pair, levels) -> bool:
    return left <= F(levels[pair[1]], levels[pair[0]]) <= right


def require_closed_cover(intervals, left: F, right: F, name: str) -> None:
    """Audit an ordered chain of closed intervals without a hidden gap."""
    require(intervals and intervals[0][0] == left and intervals[-1][1] == right,
            (name, "endpoints", intervals))
    require(all(first[1] == second[0]
                for first, second in zip(intervals, intervals[1:])),
            (name, "adjacency", intervals))


def asymptotic_limit(body, pair, left: F, right: F) -> F:
    ruler = 14 * lcm(*body)
    a, b = body[pair[0]], body[pair[1]]
    numerator = 2 * max(abs(left * a - b), abs(right * a - b))
    return F(1, 49) - numerator / ruler - singleton_debt(body, ruler, MIN_LEVEL)


def old_viable_pairs(body, left: F, right: F) -> tuple:
    return tuple(
        pair for pair in permutations(range(6), 2)
        if asymptotic_limit(body, pair, left, right) > 0
    )


def failure_cycle(edges) -> tuple[int, ...] | None:
    adjacency = {vertex: [] for vertex in range(6)}
    for left, right in edges:
        adjacency[left].append(right)

    def search(origin, vertex, path):
        for target in adjacency[vertex]:
            if target == origin:
                return tuple(path + [target])
            if target not in path:
                result = search(origin, target, path + [target])
                if result is not None:
                    return result
        return None

    for origin in range(6):
        result = search(origin, origin, [origin])
        if result is not None:
            return result
    return None


def direct_boundary_bank(body) -> tuple:
    """All-distinct exact corner (m,D)=(3,6) certificate census."""
    ruler, safe_ranges = T.R.safe_cell_ranges(body)
    cells = tuple(cell for left, right in safe_ranges for cell in range(left, right))
    values = tuple(range(3, 10))
    pair_mass = {}
    for i, j in combinations(range(6), 2):
        for qi in values:
            first_arcs = tuple(
                T.R.reflected_level_arcs(ruler, body[i], qi, cell) for cell in cells
            )
            for qj in values:
                pair_mass[i, j, qi, qj] = max(
                    T.intersection_mass(
                        first_arcs[index],
                        T.R.reflected_level_arcs(ruler, body[j], qj, cell),
                    )
                    for index, cell in enumerate(cells)
                )
    rows = []
    for levels in permutations(values, 6):
        if min(levels) != 3 or max(levels) != 9:
            continue
        debt = T.C2.C5.singleton_debt(body, ruler, levels)
        margin, pair = max(
            (pair_mass[i, j, levels[i], levels[j]] - debt, (i, j))
            for i, j in combinations(range(6), 2)
        )
        require(margin > 0, ("corner certificate failure", body, levels, margin, pair))
        rows.append((margin, levels, pair))
    require(len(rows) == 3600, (body, len(rows)))
    return tuple(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    # Rebuild the cap-3 transport maximum rather than inheriting the old split
    # conclusion.  At fixed m the worst D is 2m.  The resulting A*m/(mL-b)
    # strictly decreases with m, so D>=6 places the maximum at (D,m)=(6,3).
    bodies = T.residual_bodies()
    require(len(bodies) == 561, len(bodies))
    repeated_exceptions = {row[0] for row in T.C2.UNIVERSAL.EXPECTED_EXCEPTIONS}
    require(not (repeated_exceptions & set(bodies)),
            ("repeated-level gate", repeated_exceptions & set(bodies)))
    for body in combinations(range(1, 15), 6):
        ruler = 14 * lcm(*body)
        for pair in combinations(range(6), 2):
            a, b = body[pair[0]], body[pair[1]]
            coefficient = (b - a) + 2 * b
            corner = transport_bound(body, ruler, pair)
            limit = F(coefficient, ruler)
            require(corner - limit
                    == F(coefficient * b, ruler * (3 * ruler - b)) > 0,
                    ("cap-3 corner limit", body, pair))
            for level in (3, 4):
                now = F(level * coefficient, level * ruler - b)
                later = F((level + 1) * coefficient,
                          (level + 1) * ruler - b)
                require(now - later == F(
                    coefficient * b,
                    (level * ruler - b) * ((level + 1) * ruler - b),
                ) > 0, ("cap-3 corner monotonicity", body, pair, level))
    eta_guard = max(
        (transport_bound(body, 14 * lcm(*body), pair), body, pair,
         14 * lcm(*body))
        for body in combinations(range(1, 15), 6)
        for pair in combinations(range(6), 2)
    )
    require(eta_guard == (ETA_GUARD, H, (0, 5), 168), eta_guard)
    require(MIN_LEVEL - ETA_GUARD == F(457, 164) > 1, eta_guard)

    threshold_rows = []
    for body in bodies:
        ruler = 14 * lcm(*body)
        for pair in combinations(range(6), 2):
            threshold = coarse_threshold(body, ruler, pair)
            if threshold < T.PHASE_LIMIT:
                threshold_rows.append((threshold, body, pair))
    maximum_bound = max(
        (T.product_bound(threshold), body, pair, threshold)
        for threshold, body, pair in threshold_rows
    )
    require(len(threshold_rows) == 5666, len(threshold_rows))
    require(maximum_bound[:3] == (
        F(597757150909236102135, 22679209066623812),
        (1, 2, 6, 7, 10, 14), (0, 4)), maximum_bound)
    primitive_rows = OLD_MODULE.C.THREE.ONE.C4.primitive_universe(
        maximum_bound[0], CAP
    )
    require(len(primitive_rows) == 8796
            and max(row[0] for row in primitive_rows) == CAP,
            (len(primitive_rows), max(row[0] for row in primitive_rows)))

    closed = []
    traps = []
    witnesses = []
    reasons = {}
    constraints_by_body = {}
    verdict_digest = hashlib.sha256()
    for body in bodies:
        ruler = 14 * lcm(*body)
        constraints, _ = body_constraints(body, ruler, primitive_rows)
        witness, reason = OLD_MODULE.C.solve_constraints(body, constraints)
        reverse, reverse_reason = OLD_MODULE.C.solve_constraints(
            body, constraints, reverse=True
        )
        require((witness is None) == (reverse is None),
                ("search-order mismatch", body, reason, reverse_reason))
        verdict_digest.update(
            repr((body, witness is not None, reason, reverse_reason)).encode()
        )
        if witness is None:
            closed.append(body)
        else:
            traps.append(body)
            witnesses.append((body, witness))
            reasons[body] = reason
            constraints_by_body[body] = constraints
    trap_digest = hashlib.sha256(repr(tuple(traps)).encode()).hexdigest()
    require((len(closed), len(traps), trap_digest)
            == (COARSE_CLOSED, COARSE_TRAPS, EXPECTED_TRAP_DIGEST),
            (len(closed), len(traps), trap_digest))
    require(verdict_digest.hexdigest() == EXPECTED_VERDICT_DIGEST,
            verdict_digest.hexdigest())

    tail_bodies = tuple(body for body in traps if not constraints_by_body[body])
    require(tail_bodies == tuple(OLD_MODULE.EXPECTED_TAIL_BODIES), tail_bodies)
    require(len(tail_bodies) == 26 and {H, H2, H3} <= set(tail_bodies),
            tail_bodies)
    tail_set = set(tail_bodies)

    # A constrained survivor has an unordered gain in ``allowed`` on the
    # chosen slot pair.  Audit both ordered channels explicitly, so whichever
    # physical level is larger is covered without relabeling the body.
    located = []
    policies = []
    orientation_rows = []
    for body in traps:
        if body in tail_set:
            continue
        constraints = constraints_by_body[body]
        pair, allowed = min(constraints.items(),
                            key=lambda row: (len(row[1]), row[0]))
        channels = tuple(
            channel for ratio in sorted(allowed)
            for channel in ((ratio.denominator, ratio.numerator),
                            (ratio.numerator, ratio.denominator))
        )
        require(len(channels) == 2 * len(allowed)
                and len(set(channels)) == len(channels),
                ("oriented policy channels", body, pair, allowed, channels))
        require(all(F(max(P, Q), min(P, Q)) in allowed for P, Q in channels),
                ("unordered-to-oriented map", body, pair, allowed, channels))
        rows = tuple(OLD_MODULE.C.located_best(body, pair, channel)
                     for channel in channels)
        require(all(row[0] > 0 for row in rows),
                ("located physical policy", body, pair, rows))
        located.extend(rows)
        policies.append((body, pair, len(allowed), len(rows), min(rows)))
        orientation_rows.append((body, pair, tuple(sorted(allowed)), channels))
    require(len(policies) == 253 and len(located) == 3062,
            (len(policies), len(located)))
    located_digest = hashlib.sha256(repr(tuple(located)).encode()).hexdigest()
    orientation_digest = hashlib.sha256(
        repr(tuple(orientation_rows)).encode()
    ).hexdigest()
    require(located_digest == EXPECTED_LOCATED_DIGEST, located_digest)

    # Exhaust the old base-debt single-pair analytic family.  An asymptotic
    # limit <=0 is a proof of impossibility because its envelope is increasing.
    whole_pairs = tuple(
        (body, old_viable_pairs(body, F(1, 3), F(3))) for body in tail_bodies
    )
    no_whole = tuple(body for body, pairs in whole_pairs if not pairs)
    require(no_whole == (H, H2, H3), no_whole)
    require(sum(bool(pairs) for _, pairs in whole_pairs) == 23, whole_pairs)

    side_rows = []
    for body in no_whole:
        low = old_viable_pairs(body, F(1, 3), F(1))
        high = old_viable_pairs(body, F(1), F(3))
        # Failure of low(i,j) is q_i<q_j; failure of high(i,j) is q_j<q_i.
        edges = tuple((i, j) for i, j in low) + tuple((j, i) for i, j in high)
        side_rows.append((body, low, high, edges, failure_cycle(edges)))
    require(side_rows[0][4] is None and side_rows[1][4] is None,
            side_rows[:2])
    require(side_rows[2][4] in {(0, 1, 0), (1, 0, 1)}, side_rows[2])

    # Literal stored-atlas hostile controls, all distinct and on the integer
    # corner.  They survive the repeated-level theorem.
    h_truth = tuple(lane_applies(left, right, pair, H_STORED_HOSTILE)
                    for left, right, pair, _, _ in OLD_MODULE.H_FAREY_LANES)
    h2_low = OLD_MODULE.exact_lane(H2, F(1, 3), F(1))
    h2_high = OLD_MODULE.exact_lane(H2, F(1), F(3))
    h2_old_lanes = ((F(1, 3), F(1), h2_low[3]),
                    (F(1), F(3), h2_high[3]))
    h2_truth = tuple(lane_applies(left, right, pair, H2_STORED_HOSTILE)
                     for left, right, pair in h2_old_lanes)
    require(h_truth == (False,) * 7, h_truth)
    require(h2_truth == (False, False), h2_truth)
    corner_universe = tuple(levels for levels in permutations(range(3, 10), 6)
                            if min(levels) == 3 and max(levels) == 9)
    first_h_hostile = next(levels for levels in corner_universe
                           if not any(lane_applies(left, right, pair, levels)
                                      for left, right, pair, _, _
                                      in OLD_MODULE.H_FAREY_LANES))
    first_h2_hostile = next(levels for levels in corner_universe
                            if not any(lane_applies(left, right, pair, levels)
                                       for left, right, pair in h2_old_lanes))
    require((first_h_hostile, first_h2_hostile) ==
            (H_STORED_HOSTILE, H2_STORED_HOSTILE),
            (first_h_hostile, first_h2_hostile))

    h_boundary = direct_boundary_bank(H)
    h2_boundary = direct_boundary_bank(H2)
    boundary_digest = hashlib.sha256(repr((h_boundary, h2_boundary)).encode()).hexdigest()

    atlas = []
    for body in tail_bodies:
        if body == H:
            atlas.extend((body, *row, "H-fixed-pair") for row in H_ATLAS)
        elif body == H2:
            atlas.extend((body, *row, "H2-fixed-pair") for row in H2_ATLAS)
        elif body == H3:
            atlas.extend((body, *row, "H3-two-cycle") for row in H3_ATLAS)
        else:
            lane = OLD_MODULE.exact_lane(body, F(1, 3), F(3))
            atlas.append((body, F(1, 3), F(3), lane[3], lane[1], lane[2],
                          "whole"))
    require(len(atlas) == 35, len(atlas))
    h_intervals = tuple((row[1], row[2]) for row in atlas if row[0] == H)
    h2_intervals = tuple((row[1], row[2]) for row in atlas if row[0] == H2)
    require_closed_cover(h_intervals, F(1, 3), F(3), "H interval chain")
    require_closed_cover(h2_intervals, F(1, 3), F(3), "H2 interval chain")
    require(all(row[3] == (0, 1) for row in atlas if row[0] in {H, H2}),
            "fixed physical pair changed")
    h3_rows = tuple(row for row in atlas if row[0] == H3)
    require(tuple((row[1], row[2], row[3]) for row in h3_rows)
            == ((F(1, 3), F(1), (0, 1)),
                (F(1, 3), F(1), (1, 0))),
            ("H3 inverse-lane disjunction", h3_rows))

    lane_rows = []
    lane_monotonicity = []
    lane_channels = []
    for body, left, right, pair, start, numerator, kind in atlas:
        a, b = body[pair[0]], body[pair[1]]
        exact_numerator = 2 * max(abs(left * a - b), abs(right * a - b))
        require(numerator == exact_numerator,
                ("endpoint numerator", body, left, right, pair, numerator,
                 exact_numerator))
        margin, level, debt = lane_envelope(
            body, pair, left, right, start, numerator
        )
        previous = lane_envelope(
            body, pair, left, right, start - 1, numerator
        )[0]
        require(margin > 0 >= previous,
                ("tail threshold", body, left, right, pair, start, previous, margin))
        monotonicity = analytic_monotonicity(
            body, pair, left, right, start, numerator
        )
        channels = interval_channels(left, right, start)
        require(channels == interval_channels(left, right, start, oversize=5),
                ("head exhaustion", body, left, right, start))
        lane_rows.append((body, left, right, pair, start, numerator,
                          lane_alpha(left, right), level, debt, margin,
                          previous, len(channels), kind))
        lane_monotonicity.append(monotonicity)
        lane_channels.extend((body, pair, channel, kind) for channel in channels)
    require(len(lane_channels) == 1600, len(lane_channels))

    head_rows = tuple(
        (kind, multi_scale_head(body, pair, channel))
        for body, pair, channel, kind in lane_channels
    )
    require(max(len(row[1][6]) for row in head_rows) == 6,
            max(len(row[1][6]) for row in head_rows))
    head_digest = hashlib.sha256(repr(head_rows).encode()).hexdigest()
    lane_digest = hashlib.sha256(repr(tuple(lane_rows)).encode()).hexdigest()
    semantic_image = (
        EXPECTED_OLD_SHA256, EXPECTED_OLD_OUTPUT_SHA256, EXPECTED_OLD_SEMANTIC,
        EXPECTED_REPAIR_SHA256, EXPECTED_REPAIR_OUTPUT_SHA256,
        EXPECTED_REPAIR_SEMANTIC, tuple(bodies), eta_guard, maximum_bound,
        tuple(primitive_rows), tuple(closed), tuple(traps), tuple(witnesses),
        trap_digest, verdict_digest.hexdigest(), tuple(policies),
        tuple(orientation_rows), orientation_digest, tuple(located),
        located_digest, tail_bodies, tuple(whole_pairs),
        tuple(side_rows), H_STORED_HOSTILE, h_truth, H2_STORED_HOSTILE,
        h2_old_lanes, h2_truth, boundary_digest, min(h_boundary),
        min(h2_boundary), tuple(lane_rows), tuple(lane_monotonicity),
        lane_digest, head_rows, head_digest,
    )
    semantic = hashlib.sha256(repr(semantic_image).encode()).hexdigest()
    require((head_digest, semantic) == (EXPECTED_HEAD_DIGEST, EXPECTED_SEMANTIC),
            ("frozen digests changed", head_digest, semantic,
             EXPECTED_HEAD_DIGEST, EXPECTED_SEMANTIC))

    weakest_head = min(row[1] for row in head_rows)
    weakest_located = min(located)
    longest = max(head_rows, key=lambda row: len(row[1][6]))
    lines = [
        "LRC14 reflected cap-3 physical-orientation half-cone repair",
        "predecessor_half_cone_conclusion=RETRACTED_AND_NOT_INHERITED;partition_and_policies_recomputed_from_primitive_predicates",
        f"universe=bodies:3003;inherited_arbitrary_level_bank:2442;recomputed_residual:{len(bodies)};spread:D>=6;cone:2m>=D;ratio_cap:3",
        f"transport_corner=(D,m)=(6,3);max_abs_eta={qtext(eta_guard[0])}@body={eta_guard[1]},pair={eta_guard[2]};constrained_edges={len(threshold_rows)};primitive_bank={len(primitive_rows)}",
        f"coarse_csp=closed:{len(closed)};traps:{len(traps)};constrained_policies:{len(policies)};unconstrained_tails:{len(tail_bodies)};trap_digest={trap_digest};verdict_digest={verdict_digest.hexdigest()}",
        f"physical_policy_audit=both_orientations_per_unordered_gain;controls:{len(located)};weakest={qtext(weakest_located[0])}@body={weakest_located[1]},pair={weakest_located[2]},channel={weakest_located[3]};orientation_digest={orientation_digest};located_digest={located_digest}",
        f"old_single_pair_family=whole_interval:{len(tail_bodies)-len(no_whole)};no_whole:{no_whole}",
        f"old_full_side_failure_graphs={tuple(side_rows)}",
        f"stored_H_hostile={H_STORED_HOSTILE};lane_truth={h_truth};stored_H2_hostile={H2_STORED_HOSTILE};lane_truth={h2_truth};all_distinct;corner:m=3,D=6",
        f"corner_direct_certificate=H:3600/3600,H2:3600/3600;weakest_H={qtext(min(h_boundary)[0])}@{min(h_boundary)[1]},{min(h_boundary)[2]};weakest_H2={qtext(min(h2_boundary)[0])}@{min(h2_boundary)[1]},{min(h2_boundary)[2]};digest={boundary_digest}",
        "enrichment=on ratio lane [u,v],pair_max>=alpha*reduced_min with alpha=1/v below 1,u above 1,1 across 1;q_max<=3m gives m>=ceil(alpha*s/3)",
        f"new_atlas=assignment_complete_lanes:{len(lane_rows)};head_channels:{len(head_rows)};whole_bodies:23;H3_failure_two_cycle;H_fixed_pair_lanes:{len(H_ATLAS)};H2_fixed_pair_lanes:{len(H2_ATLAS)}",
        *(f"lane_{index+1}=body:{row[0]};interval:[{qtext(row[1])},{qtext(row[2])}];pair:{row[3]};start:{row[4]};N:{qtext(row[5])};alpha:{qtext(row[6])};debt_level:{row[7]};margin:{qtext(row[9])};previous:{qtext(row[10])};heads:{row[11]};kind:{row[12]}" for index, row in enumerate(lane_rows)),
        f"head_audit=PASS:{len(head_rows)};weakest_tail_margin={qtext(weakest_head[0])}@body={weakest_head[1]},pair={weakest_head[2]},channel={weakest_head[3]},g_tail={weakest_head[5][1]};maximum_direct_prefix={len(longest[1][6])}@body={longest[1][1]},channel={longest[1][3]},g0={longest[1][4]},g_tail={longest[1][5][1]}",
        "infinite_suffix=primitive skeleton fixed;phase displacement nonincreasing;transport displacement strictly decreasing;ceiling debt nonincreasing;first-tail transport and every earlier g checked exactly",
        f"lane_digest={lane_digest};head_digest={head_digest}",
        "proved_statement=within_the_pinned_reflected_THM-2941_k=1_sufficient_family,every_packet_with_D>=6_and_2m>=D_closes_on_all_3003_bodies;the_old_failure_was_a_tail_selector_obstruction_not_a_boundary_runner_obstruction",
        "corollary=the_assembled_reflected_certificate-failure_wedge_is_confined_to_561_bodies_with_D>=6_and_1<=m<D/2",
        "not_proved=physical_LRC14_or_absence_of_survivors_outside_the_reflected_THM-2941_sufficient_family",
        "normal_vs_python_O=BYTE_IDENTICAL",
        f"old_source_sha256={EXPECTED_OLD_SHA256}",
        f"old_output_sha256={EXPECTED_OLD_OUTPUT_SHA256}",
        f"old_semantic_sha256={EXPECTED_OLD_SEMANTIC}",
        f"two_thirds_repair_source_sha256={EXPECTED_REPAIR_SHA256}",
        f"two_thirds_repair_output_sha256={EXPECTED_REPAIR_OUTPUT_SHA256}",
        f"two_thirds_repair_semantic_sha256={EXPECTED_REPAIR_SEMANTIC}",
        f"source_sha256={sha256(HERE)}",
        f"semantic_sha256={semantic}",
        "all_exact_controls=PASS",
    ]
    text = "\n".join(lines) + "\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(text, encoding="utf-8")
    print(text, end="")


if __name__ == "__main__":
    main()
