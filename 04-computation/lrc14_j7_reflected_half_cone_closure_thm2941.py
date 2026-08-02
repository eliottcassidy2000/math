#!/usr/bin/env python3
r"""Exact proof candidate for the reflected ``2m >= D`` cone.

Inside the sufficient reflected ``k=1`` family of THM-2941, this referee
closes every residual packet with

    D >= 6,                 2m >= D.

At projective cap 3 the new phase-zero gain 3 is exponent-dependent:
``(3/2)*2=3``.  Every forced-full component contains a unique zero-gain K4
on relative scales ``{1,3/2,2,3}``, hence two copies of both the old and new
balanced triangles.  Complete located policies and exact analytic tails close
every coarse trap.

This is not a physical-survivor census and not a proof of LRC(14).
"""

from __future__ import annotations

import argparse
import hashlib
from collections import Counter
from fractions import Fraction as F
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, permutations
from math import gcd, lcm
from pathlib import Path


HERE = Path(__file__).resolve()
ROOT = next(parent for parent in HERE.parents if (parent / "00-navigation").is_dir())
TWO_THIRDS = ROOT / "04-computation/lrc14_j7_reflected_two_thirds_cone_closure_thm2941.py"
TWO_THIRDS_OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_reflected_two_thirds_cone_closure_thm2941.out"
OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_reflected_half_cone_closure_thm2941.out"

EXPECTED_TWO_THIRDS_SHA256 = "e6e64c909c6bfcc776eaf6bf2ad210f75675a6a32a46ead70d8d29a37f607eb3"
EXPECTED_TWO_THIRDS_OUTPUT_SHA256 = "e77929c87f9d9f8fb7ce3a347c48522e87f63fa7c7085eb9e2cd8fe0bb4e4a90"
EXPECTED_TWO_THIRDS_SEMANTIC = "b34db3e25b9e4c81c1549f1ec7c7ab78e935ec778610daf617934b66b3a47304"
EXPECTED_TRAP_DIGEST = "35fceb70d1e960ee5e0590c990adde5c4848b657a15455714585354bee322380"
EXPECTED_VERDICT_DIGEST = "fe8193713f596b6f3f0182066a8937cfe75f681ea760fd31bac4642b7eb0f7cc"
EXPECTED_FORCED_DIGEST = "70eab59654ebeee0fbcde8b196a2f06d10c13369e8a20ea7fadcd5a0f566eaad"
EXPECTED_K4_DIGEST = "2a90a1d9ff6b3dcdada4189eaef83d5534e8f53d2f71461dae3cb917d938d8fc"
EXPECTED_LOCATED_DIGEST = "8968666298108aabbf926de452354ed554728aa2fe2cc2ed5f5aaea3598f4201"
EXPECTED_LANE_DIGEST = "fa98d31b75ce2682d09e097cf4bcfcca4b00e0c01cac1b5a4b7f10b7676a857b"
EXPECTED_HEAD_DIGEST = "e89a7f20d72228870e488094333cf9f84d45f9f455d1bd5448a03eccca14642a"
EXPECTED_SEMANTIC_SHA256 = "e7edbeec8933409bb36eab7d4a1867a10b3f7683a4898ec30f244a9e81601673"

MIN_SPREAD = 6
MIN_LEVEL = 3
RATIO_CAP = F(3)
ETA_GUARD = F(35, 164)
COARSE_CLOSED = 282
COARSE_TRAPS = 279

H = (1, 2, 3, 4, 6, 12)
EXPECTED_TAIL_BODIES = (
    (1, 2, 3, 4, 6, 8),
    (1, 2, 3, 4, 6, 9),
    H,
    (1, 2, 3, 4, 8, 12),
    (1, 2, 3, 4, 9, 12),
    (1, 2, 3, 5, 6, 10),
    (1, 2, 3, 6, 8, 12),
    (1, 2, 3, 6, 9, 12),
    (1, 2, 4, 5, 8, 10),
    (1, 2, 4, 6, 8, 12),
    (1, 2, 4, 6, 9, 12),
    (1, 3, 4, 5, 6, 10),
    (1, 3, 4, 5, 6, 12),
    (1, 3, 4, 5, 10, 12),
    (1, 3, 4, 6, 8, 12),
    (1, 3, 4, 6, 9, 12),
    (1, 3, 4, 6, 10, 12),
    (1, 3, 5, 6, 10, 12),
    (1, 4, 5, 6, 10, 12),
    (1, 4, 6, 7, 12, 14),
    (1, 4, 6, 8, 9, 12),
    (2, 3, 4, 6, 8, 12),
    (2, 3, 4, 6, 9, 12),
    (2, 4, 5, 6, 10, 12),
    (2, 4, 6, 8, 9, 12),
    (3, 4, 5, 6, 10, 12),
)

# The only difficult tail body uses seven exact Farey charts.  The tuple is
# (left endpoint, right endpoint, ordered slot pair, first tail scale,
# exact endpoint numerator bound).
H_FAREY_LANES = (
    (F(1, 3), F(8, 11), (1, 0), 6, F(10, 11)),
    (F(8, 11), F(11, 12), (3, 2), 8, F(4, 3)),
    (F(11, 12), F(1), (1, 0), 24, F(2)),
    (F(1), F(13, 12), (0, 1), 24, F(2)),
    (F(13, 12), F(17, 11), (2, 3), 9, F(3, 2)),
    (F(17, 11), F(29, 12), (0, 1), 6, F(10, 11)),
    (F(29, 12), F(3), (0, 2), 7, F(7, 6)),
)

ZERO_VECTORS = {
    F(4, 3): (2, -1, 0),
    F(3, 2): (-1, 1, 0),
    F(2): (1, 0, 0),
    F(5, 2): (-1, 0, 1),
    F(3): (0, 1, 0),
}
FORCED_PROFILE = {
    (4, 6, 3): 3,
    (5, 6, 3): 1,
    (5, 7, 3): 30,
    (6, 7, 3): 7,
    (6, 8, 3): 12,
    (6, 8, 4): 2,
}


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


require(sha256(TWO_THIRDS) == EXPECTED_TWO_THIRDS_SHA256,
        ("two-thirds source changed", sha256(TWO_THIRDS)))
require(sha256(TWO_THIRDS_OUTPUT) == EXPECTED_TWO_THIRDS_OUTPUT_SHA256,
        ("two-thirds output changed", sha256(TWO_THIRDS_OUTPUT)))
require(f"semantic_sha256={EXPECTED_TWO_THIRDS_SEMANTIC}" in
        TWO_THIRDS_OUTPUT.read_text(), "two-thirds semantic token missing")
C = import_module("half_cone_two_thirds", TWO_THIRDS)
T = C.T
C.MIN_LEVEL = MIN_LEVEL
C.RATIO_CAP = RATIO_CAP
C.ETA_GUARD = ETA_GUARD
T.MIN_LEVEL = MIN_LEVEL
T.MIN_SPREAD = MIN_SPREAD
T.RATIO_CAP = RATIO_CAP


def transport_bound(body, ruler, pair) -> F:
    """Exact cone maximum, attained at ``(D,m)=(6,3)``."""
    a, b = body[pair[0]], body[pair[1]]
    return F(3 * (b - a) + 6 * b, 3 * ruler - b)


def coarse_threshold(body, ruler, pair) -> F:
    return 2 * transport_bound(body, ruler, pair) + T.debt_at_nine(body, ruler)


def body_constraints(body, ruler, primitive_rows):
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


def interval_channels(left: F, right: F, threshold: int, oversize: int = 3):
    return tuple(
        (P, Q)
        for P in range(1, oversize * threshold + 1)
        for Q in range(1, oversize * threshold + 1)
        if P != Q and gcd(P, Q) == 1 and min(P, Q) < threshold
        and left <= F(Q, P) <= right
    )


def exact_lane(body, left: F, right: F, forced_pair=None,
               forced_start=None, forced_numerator=None):
    """Choose or verify one uniform analytic tail lane."""
    candidates = []
    pairs = (forced_pair,) if forced_pair is not None else tuple(permutations(range(6), 2))
    for pair in pairs:
        a, b = body[pair[0]], body[pair[1]]
        numerator = 2 * max(abs(left * a - b), abs(right * a - b))
        if forced_numerator is not None:
            require(numerator == forced_numerator,
                    ("lane numerator", body, left, right, pair, numerator, forced_numerator))
        starts = (forced_start,) if forced_start is not None else range(2, 61)
        for start in starts:
            margin = T.tail_envelope(body, pair, start, numerator)
            if margin <= 0:
                continue
            previous = T.tail_envelope(body, pair, start - 1, numerator)
            step = T.tail_step_gain(body, pair, start, numerator)
            ruler = 14 * lcm(*body)
            eta_bound = numerator / 2 / (ruler - F(a, start))
            require(previous <= 0 and step > 0 and eta_bound < 1,
                    ("lane threshold", body, left, right, pair, start,
                     previous, margin, step, eta_bound))
            channels = interval_channels(left, right, start)
            require(channels == interval_channels(left, right, start, oversize=5),
                    ("head exhaustion", body, left, right, start))
            candidates.append((len(channels), start, numerator, pair, margin,
                               previous, step, eta_bound, channels))
            break
    require(candidates, ("no analytic lane", body, left, right))
    return min(candidates, key=lambda row: (row[0], row[1], row[2], row[3]))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    bodies = T.residual_bodies()
    universal_exceptions = {row[0] for row in T.C2.UNIVERSAL.EXPECTED_EXCEPTIONS}
    require(not (universal_exceptions & set(bodies)),
            ("repeated-level discharge", universal_exceptions & set(bodies)))

    # On 2m>=D the largest D at fixed m is 2m.  With A=delta+2b,
    # B=mA/(mL-b) strictly decreases in m and D>=6 forces m>=3.
    for body in combinations(range(1, 15), 6):
        ruler = 14 * lcm(*body)
        for pair in combinations(range(6), 2):
            a, b = body[pair[0]], body[pair[1]]
            A = b - a + 2 * b
            corner = transport_bound(body, ruler, pair)
            limit = F(A, ruler)
            require(corner - limit == F(A * b, ruler * (3 * ruler - b)) > 0,
                    ("corner limit", body, pair))
            for m in (3, 4):
                now = F(m * A, m * ruler - b)
                later = F((m + 1) * A, (m + 1) * ruler - b)
                require(now - later == F(A * b,
                        (m * ruler - b) * ((m + 1) * ruler - b)) > 0,
                        ("corner monotonicity", body, pair, m))

    eta_guard = max(
        (transport_bound(body, 14 * lcm(*body), pair), body, pair, 14 * lcm(*body))
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
    primitive_rows = C.THREE.ONE.C4.primitive_universe(maximum_bound[0], RATIO_CAP)
    require(len(primitive_rows) == 8796 and max(row[0] for row in primitive_rows) == RATIO_CAP,
            (len(primitive_rows), max(row[0] for row in primitive_rows)))
    phase_zero = tuple((P, Q) for _, phase, _, P, Q in primitive_rows if phase == 0)
    require(phase_zero == ((1, 2), (1, 3), (2, 3), (2, 5), (3, 4)), phase_zero)

    old_gains = (F(4, 3), F(3, 2), F(2), F(5, 2))
    old_vectors = tuple(ZERO_VECTORS[gain] for gain in old_gains)
    all_vectors = old_vectors + (ZERO_VECTORS[F(3)],)
    old_rank = C.rational_rank(old_vectors)
    full_rank = C.rational_rank(all_vectors)
    require((old_rank, len(old_vectors) - old_rank) == (3, 1), old_vectors)
    require((full_rank, len(all_vectors) - full_rank) == (3, 2), all_vectors)
    require(tuple(x + y for x, y in zip(ZERO_VECTORS[F(3, 2)], ZERO_VECTORS[F(2)]))
            == ZERO_VECTORS[F(3)], "(3/2)*2=3 exponent relation")

    closed = []
    traps = []
    witnesses = []
    reasons = {}
    constraints_by_body = {}
    verdict_digest = hashlib.sha256()
    for body in bodies:
        ruler = 14 * lcm(*body)
        constraints, _ = body_constraints(body, ruler, primitive_rows)
        witness, reason = C.solve_constraints(body, constraints)
        reverse, reverse_reason = C.solve_constraints(body, constraints, reverse=True)
        require((witness is None) == (reverse is None),
                ("search-order mismatch", body, reason, reverse_reason))
        verdict_digest.update(repr((body, witness is not None, reason, reverse_reason)).encode())
        if witness is None:
            closed.append(body)
        else:
            traps.append(body)
            witnesses.append((body, witness))
            reasons[body] = reason
            constraints_by_body[body] = constraints
    trap_digest = hashlib.sha256(repr(tuple(traps)).encode()).hexdigest()
    require((len(closed), len(traps), trap_digest) ==
            (COARSE_CLOSED, COARSE_TRAPS, EXPECTED_TRAP_DIGEST),
            (len(closed), len(traps), trap_digest))
    require(verdict_digest.hexdigest() == EXPECTED_VERDICT_DIGEST,
            verdict_digest.hexdigest())
    witness_spans = tuple((body, F(max(witness), min(witness)))
                          for body, witness in witnesses)
    require(max(span for _, span in witness_spans) == RATIO_CAP, witness_spans)
    require(sum(span == RATIO_CAP for _, span in witness_spans) == 55, witness_spans)

    zero_gains = set(ZERO_VECTORS)
    profile = Counter()
    triangle_profile = Counter()
    gainset_profile = Counter()
    forced_rows = []
    k4_rows = []
    for body, witness in witnesses:
        forced = tuple(choice[0] for choice in reasons[body][0] if not choice[1])
        require(len(forced) <= 1, ("multiple forced components", body, forced))
        if not forced:
            continue
        component = forced[0]
        zero_edges = []
        gains_by_edge = {}
        gain_counts = Counter()
        for edge in combinations(component, 2):
            if edge not in constraints_by_body[body]:
                continue
            gain = F(max(witness[edge[0]], witness[edge[1]]),
                     min(witness[edge[0]], witness[edge[1]]))
            if gain in zero_gains:
                zero_edges.append(edge)
                gains_by_edge[edge] = gain
                gain_counts[gain] += 1

        old_triangles = []
        new_triangles = []
        for triangle in combinations(component, 3):
            edges = tuple(combinations(triangle, 2))
            if not all(edge in gains_by_edge for edge in edges):
                continue
            gains = {gains_by_edge[edge] for edge in edges}
            if gains == {F(4, 3), F(3, 2), F(2)}:
                old_triangles.append(triangle)
            if gains == {F(3, 2), F(2), F(3)}:
                new_triangles.append(triangle)

        k4_cores = []
        expected_k4_multiset = Counter({F(4, 3): 1, F(3, 2): 2,
                                        F(2): 2, F(3): 1})
        for quad in combinations(component, 4):
            edges = tuple(combinations(quad, 2))
            if (all(edge in gains_by_edge for edge in edges)
                    and Counter(gains_by_edge[edge] for edge in edges) == expected_k4_multiset):
                k4_cores.append(quad)

        used = {vertex for edge in zero_edges for vertex in edge}
        components = C.component_count(used, zero_edges)
        cycle_rank = len(zero_edges) - len(used) + components
        gain_three_edges = tuple(edge for edge in zero_edges if gains_by_edge[edge] == 3)
        require(len(gain_three_edges) == 1, ("gain-three multiplicity", body, gain_three_edges))
        gain_three_edge = gain_three_edges[0]
        is_bridge = (C.component_count(
            used, tuple(edge for edge in zero_edges if edge != gain_three_edge)) == components + 1)
        require(not is_bridge, ("gain three unexpectedly bridges", body, gain_three_edge))
        require(len(k4_cores) == 1 and all(v in k4_cores[0] for v in gain_three_edge),
                ("unique doubled core", body, k4_cores, gain_three_edge))
        require(len(old_triangles) >= 2 and len(new_triangles) == 2,
                ("triangle multiplicity", body, old_triangles, new_triangles))
        require({F(4, 3), F(3, 2), F(2), F(3)} <= set(gain_counts),
                ("missing doubled-core gain", body, gain_counts))

        profile[(len(component), len(zero_edges), cycle_rank)] += 1
        triangle_profile[(len(old_triangles), len(new_triangles))] += 1
        gainset_profile[tuple(sorted(gain_counts))] += 1
        forced_rows.append((body, component, witness, tuple(zero_edges),
                            tuple(old_triangles), tuple(new_triangles),
                            gain_three_edges, (is_bridge,), cycle_rank,
                            tuple(sorted(gain_counts.items()))))
        k4_rows.append((body, component, tuple(zero_edges), tuple(k4_cores), gain_three_edge))

    require(len(forced_rows) == 55 and dict(profile) == FORCED_PROFILE,
            (len(forced_rows), profile))
    require(dict(triangle_profile) == {(2, 2): 53, (3, 2): 2}, triangle_profile)
    require(gainset_profile[
        (F(4, 3), F(3, 2), F(2), F(3))] == 6, gainset_profile)
    require(gainset_profile[
        (F(4, 3), F(3, 2), F(2), F(5, 2), F(3))] == 49, gainset_profile)
    forced_digest = hashlib.sha256(repr(tuple(forced_rows)).encode()).hexdigest()
    k4_digest = hashlib.sha256(repr(tuple(k4_rows)).encode()).hexdigest()
    require(forced_digest == EXPECTED_FORCED_DIGEST, forced_digest)
    require(k4_digest == EXPECTED_K4_DIGEST, k4_digest)

    tail_bodies = tuple(body for body in traps if not constraints_by_body[body])
    require(tail_bodies == EXPECTED_TAIL_BODIES, tail_bodies)
    tail_set = set(tail_bodies)

    located = []
    policies = []
    for body in traps:
        if body in tail_set:
            continue
        constraints = constraints_by_body[body]
        pair, allowed = min(constraints.items(), key=lambda row: (len(row[1]), row[0]))
        channels = tuple(
            channel for ratio in sorted(allowed)
            for channel in ((ratio.denominator, ratio.numerator),
                            (ratio.numerator, ratio.denominator))
        )
        rows = tuple(C.located_best(body, pair, channel) for channel in channels)
        require(all(row[0] > 0 for row in rows), ("located policy", body, pair, rows))
        located.extend(rows)
        policies.append((body, pair, len(rows), min(rows)))
    require(len(policies) == 253 and len(located) == 3062,
            (len(policies), len(located)))
    located_digest = hashlib.sha256(repr(tuple(located)).encode()).hexdigest()
    require(located_digest == EXPECTED_LOCATED_DIGEST, located_digest)

    lanes = []
    for body in tail_bodies:
        if body == H:
            for left, right, pair, start, numerator in H_FAREY_LANES:
                lane = exact_lane(body, left, right, pair, start, numerator)
                lanes.append((body, left, right, lane))
        else:
            lanes.append((body, F(1, 3), F(1), exact_lane(body, F(1, 3), F(1))))
            lanes.append((body, F(1), F(3), exact_lane(body, F(1), F(3))))
    require(len(lanes) == 57 and sum(lane[3][0] for lane in lanes) == 792,
            (len(lanes), sum(lane[3][0] for lane in lanes)))

    all_heads = []
    lane_rows = []
    for body, left, right, lane in lanes:
        count, start, numerator, pair, margin, previous, step, eta_bound, channels = lane
        rows = tuple(C.located_best(body, pair, channel) for channel in channels)
        require(len(rows) == count and all(row[0] > 0 for row in rows),
                ("head failure", body, left, right, pair, rows))
        all_heads.extend(rows)
        lane_rows.append((body, left, right, pair, start, numerator, margin,
                          previous, step, eta_bound, count, min(rows)))
    require(len(all_heads) == 792, len(all_heads))
    head_digest = hashlib.sha256(repr(tuple(all_heads)).encode()).hexdigest()
    lane_digest = hashlib.sha256(repr(tuple(lane_rows)).encode()).hexdigest()
    require(head_digest == EXPECTED_HEAD_DIGEST, head_digest)
    require(lane_digest == EXPECTED_LANE_DIGEST, lane_digest)

    semantic_image = (
        tuple(bodies), eta_guard, maximum_bound, tuple(primitive_rows), phase_zero,
        tuple(closed), tuple(traps), tuple(witnesses), witness_spans,
        trap_digest, verdict_digest.hexdigest(), tuple(forced_rows), forced_digest,
        tuple(k4_rows), k4_digest, tuple(sorted(profile.items())),
        tuple(policies), tuple(located), located_digest,
        tuple(lane_rows), lane_digest, tuple(all_heads), head_digest,
    )
    semantic_sha = hashlib.sha256(repr(semantic_image).encode()).hexdigest()
    require(semantic_sha == EXPECTED_SEMANTIC_SHA256,
            ("semantic image changed", semantic_sha, EXPECTED_SEMANTIC_SHA256))

    weakest_located = min(located)
    weakest_head = min(all_heads)
    lines = [
        "LRC14 reflected half-cone exact proof candidate",
        "universe=bodies:3003;arbitrary_level_bank:2442;residual:561;spread:D>=6;cone:2m>=D",
        f"corner=integer boundary D=6,m=3;ratio_cap=3;max_abs_eta={qtext(eta_guard[0])};intermediate_slope_floor={qtext(MIN_LEVEL-eta_guard[0])}",
        f"phase=1/49+c/(PQ),c>=-12/49;constrained_edges={len(threshold_rows)};max_product_bound={qtext(maximum_bound[0])}@{maximum_bound[1]},{maximum_bound[2]};primitive_bank={len(primitive_rows)};cap_endpoint=(1,3)",
        f"phase_zero_channels={phase_zero};old_exponent_rank_nullity=3,1;full_exponent_rank_nullity=3,2",
        f"coarse_csp=closed:{len(closed)};traps:{len(traps)};reverse_order_agrees;largest_witness_span=3;forced_full_components={len(forced_rows)}",
        f"coarse_traps={tuple(traps)}",
        f"trap_digest={trap_digest};verdict_digest={verdict_digest.hexdigest()}",
        f"forced_profile={tuple(sorted(profile.items()))};triangle_profile={tuple(sorted(triangle_profile.items()))};unique_doubled_K4=55;gain_3_cycle_edge_in_all=55;forced_digest={forced_digest};k4_digest={k4_digest}",
        f"located_policies={len(policies)};oriented_controls={len(located)};weakest={qtext(weakest_located[0])}@body={weakest_located[1]},pair={weakest_located[2]},channel={weakest_located[3]},cell={weakest_located[5]};digest={located_digest}",
        *(f"lane_{index+1}=body:{row[0]};interval:[{qtext(row[1])},{qtext(row[2])}];ordered_pair:{row[3]};start:{row[4]};numerator:{qtext(row[5])};margin:{qtext(row[6])};previous:{qtext(row[7])};first_step:{qtext(row[8])};eta_bound:{qtext(row[9])};head_channels:{row[10]};head_weakest:{qtext(row[11][0])}@{row[11][3]},cell={row[11][5]}" for index, row in enumerate(lane_rows)),
        f"tail_lanes={len(lane_rows)};head_controls={len(all_heads)};weakest={qtext(weakest_head[0])}@body={weakest_head[1]},pair={weakest_head[2]},channel={weakest_head[3]},cell={weakest_head[5]};lane_digest={lane_digest};head_digest={head_digest}",
        "transport=exact integer corner;oriented eta_g retained;g/(PgL-a) strictly decreases;c_inverse dropped only after positive bracket",
        "tail=PQ>=s^2,Pg>=s;endpoint numerator bounds exact on each closed Farey lane;losses decrease strictly;finite heads 3s-vs-5s audited",
        "conclusion=all reflected residual packets with D>=6 and 2m>=D close on all 3003 bodies",
        "corollary=assembled reflected certificate-failure wedge is confined to 561 bodies,D>=6,1<=m<D/2",
        "scope=reflected THM-2941 sufficient family only;physical LRC14 remains open",
        "normal_vs_python_O=BYTE_IDENTICAL",
        f"two_thirds_source_sha256={EXPECTED_TWO_THIRDS_SHA256}",
        f"two_thirds_output_sha256={EXPECTED_TWO_THIRDS_OUTPUT_SHA256}",
        f"two_thirds_semantic_sha256={EXPECTED_TWO_THIRDS_SEMANTIC}",
        f"source_sha256={sha256(HERE)}",
        f"semantic_sha256={semantic_sha}",
        "all_exact_controls=PASS",
    ]
    text = "\n".join(lines) + "\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(text, encoding="utf-8")
    print(text, end="")


if __name__ == "__main__":
    main()
