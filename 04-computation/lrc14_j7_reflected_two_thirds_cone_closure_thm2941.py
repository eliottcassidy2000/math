#!/usr/bin/env python3
r"""Exact proof candidate for the reflected ``3m >= 2D`` cone.

Inside the sufficient reflected ``k=1`` family of THM-2941, this referee
closes every residual packet with

    D >= 6,                 3m >= 2D.

At projective cap 5/2 the new phase-zero gain 5/2 is prime-exponent
independent of the old balanced triangle.  The exact CSP leaves 172 traps;
65 forced-full components consist of an old zero-holonomy core plus a 5/2
bridge.  Complete located policies and analytic tails close every trap.

This is not a physical-survivor census and not a proof of LRC(14).
"""

from __future__ import annotations

import argparse
import hashlib
from collections import Counter
from fractions import Fraction as F
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from math import gcd, lcm
from pathlib import Path


HERE = Path(__file__).resolve()
ROOT = next(parent for parent in HERE.parents if (parent / "00-navigation").is_dir())
THREE_QUARTER = ROOT / "04-computation/lrc14_j7_reflected_three_quarter_cone_closure_thm2941.py"
THREE_QUARTER_OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_reflected_three_quarter_cone_closure_thm2941.out"
OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_reflected_two_thirds_cone_closure_thm2941.out"

EXPECTED_THREE_QUARTER_SHA256 = "85bb9bd1613abd5cd7a877958b5c89a10172a5034c37cae5cc0467bd8ba4c0d3"
EXPECTED_THREE_QUARTER_OUTPUT_SHA256 = "451823bd9bb11ae4af5c1fb4675fba2f163d172c178eae26c4cd115eb945ba7c"
EXPECTED_THREE_QUARTER_SEMANTIC = "05a900c80283d5bc9a7b01c1c4ad045b889aebe2e4d5798eaa66cb8907a7ae9f"
EXPECTED_TRAP_DIGEST = "b8e1f964d610af641a229d06951c4805167af4310429740501809034e9b2a716"
EXPECTED_SEMANTIC_SHA256 = "b34db3e25b9e4c81c1549f1ec7c7ab78e935ec778610daf617934b66b3a47304"

MIN_SPREAD = 6
MIN_LEVEL = 4
RATIO_CAP = F(5, 2)
ETA_GUARD = F(29, 165)
COARSE_CLOSED = 389
COARSE_TRAPS = 172

H = (1, 2, 3, 4, 6, 12)
H2 = (1, 3, 4, 6, 8, 12)

TAIL_POLICIES = (
    ((1, 2, 3, 4, 6, 8), "all", (0, 1), 6, F(16, 5)),
    (H, "low", (1, 0), 10, F(2)),
    (H, "high", (0, 1), 10, F(2)),
    ((1, 2, 3, 4, 8, 12), "all", (0, 1), 6, F(16, 5)),
    ((1, 2, 3, 5, 6, 10), "all", (0, 1), 5, F(16, 5)),
    ((1, 2, 3, 6, 8, 12), "all", (0, 1), 6, F(16, 5)),
    ((1, 2, 4, 6, 8, 12), "all", (0, 1), 6, F(16, 5)),
    (H2, "low", (2, 1), 6, F(14, 5)),
    (H2, "high", (0, 1), 8, F(4)),
    ((1, 3, 4, 6, 9, 12), "all", (0, 1), 6, F(26, 5)),
    ((1, 4, 5, 6, 10, 12), "all", (0, 1), 5, F(36, 5)),
    ((2, 3, 4, 6, 8, 12), "all", (0, 1), 9, F(22, 5)),
    ((2, 3, 4, 6, 9, 12), "all", (0, 1), 6, F(22, 5)),
    ((2, 4, 5, 6, 10, 12), "all", (0, 1), 5, F(32, 5)),
)

ZERO_VECTORS = {
    F(4, 3): (2, -1, 0),
    F(3, 2): (-1, 1, 0),
    F(2): (1, 0, 0),
    F(5, 2): (-1, 0, 1),
}
FORCED_PROFILE = {
    (4, 4, 1): 8,
    (5, 4, 1): 24,
    (5, 5, 1): 4,
    (6, 4, 1): 20,
    (6, 5, 1): 2,
    (6, 7, 2): 7,
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


def rational_rank(rows) -> int:
    matrix = [[F(value) for value in row] for row in rows]
    rank = 0
    for column in range(3):
        pivot = next((r for r in range(rank, len(matrix)) if matrix[r][column]), None)
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        scale = matrix[rank][column]
        matrix[rank] = [value / scale for value in matrix[rank]]
        for r in range(len(matrix)):
            if r == rank or not matrix[r][column]:
                continue
            factor = matrix[r][column]
            matrix[r] = [x - factor * y for x, y in zip(matrix[r], matrix[rank])]
        rank += 1
    return rank


def component_count(vertices, edges) -> int:
    unseen = set(vertices)
    count = 0
    while unseen:
        count += 1
        stack = [unseen.pop()]
        while stack:
            vertex = stack.pop()
            for left, right in edges:
                other = right if left == vertex else left if right == vertex else None
                if other in unseen:
                    unseen.remove(other)
                    stack.append(other)
    return count


require(sha256(THREE_QUARTER) == EXPECTED_THREE_QUARTER_SHA256,
        ("three-quarter source changed", sha256(THREE_QUARTER)))
require(sha256(THREE_QUARTER_OUTPUT) == EXPECTED_THREE_QUARTER_OUTPUT_SHA256,
        ("three-quarter output changed", sha256(THREE_QUARTER_OUTPUT)))
require(f"semantic_sha256={EXPECTED_THREE_QUARTER_SEMANTIC}" in
        THREE_QUARTER_OUTPUT.read_text(), "three-quarter semantic token missing")
THREE = import_module("two_thirds_three_quarter", THREE_QUARTER)
T = THREE.T
T.MIN_LEVEL = MIN_LEVEL
T.MIN_SPREAD = MIN_SPREAD
T.RATIO_CAP = RATIO_CAP


def transport_bound(body, ruler, pair) -> F:
    a, b = body[pair[0]], body[pair[1]]
    return F(4 * (b - a) + 6 * b, 4 * ruler - b)


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


def solve_constraints(body, constraints, reverse=False):
    """Solve one already-generated CSP in either deterministic order."""
    if any(not allowed for allowed in constraints.values()):
        return None, ("empty edge",)
    choices = []
    component_rows = []
    nodes = 0
    forced_full = 0
    for component in T.connected_components(constraints):
        short, full, count = T.component_types(component, constraints, reverse=reverse)
        nodes += count
        if short is None and full is None:
            return None, ("empty component", component, nodes)
        if short is None:
            forced_full += 1
            component_rows.append(full)
        else:
            component_rows.append(short)
        choices.append((component, short is not None, full is not None, count))
    if forced_full > 1:
        return None, ("two forced full components", tuple(choices), nodes)
    witness = T.scale_components(tuple(component_rows))
    for pair, allowed in constraints.items():
        i, j = pair
        ratio = F(max(witness[i], witness[j]), min(witness[i], witness[j]))
        require(ratio in allowed, (body, witness, pair, ratio, allowed))
    return witness, (tuple(choices), nodes)


def located_best(body, pair, channel):
    ruler, safe_ranges = T.R.safe_cell_ranges(body)
    a, b = body[pair[0]], body[pair[1]]
    P, Q = channel
    g0 = (MIN_LEVEL + min(P, Q) - 1) // min(P, Q)
    debt = T.debt_at_nine(body, ruler)
    eta = F(g0 * (Q * a - P * b), P * g0 * ruler - a)
    rows = []
    for left, right in safe_ranges:
        for cell in range(left, right):
            skeleton = T.intersection_mass(
                T.primitive_arcs(P, F(a * cell, ruler) % 1),
                T.primitive_arcs(Q, F(b * cell, ruler) % 1),
            )
            rows.append((skeleton - 2 * abs(eta) - debt, cell, skeleton))
    margin, cell, skeleton = max(rows)
    bracket = skeleton - 2 * abs(eta)
    c = 1 - F(a, P * g0 * ruler)
    transported = bracket / c
    actual = T.intersection_mass(
        T.R.reflected_level_arcs(ruler, a, P * g0, cell),
        T.R.reflected_level_arcs(ruler, b, Q * g0, cell),
    )
    require(min(P, Q) * g0 >= MIN_LEVEL and abs(eta) <= ETA_GUARD < 1,
            ("homotopy gate", body, pair, channel, g0, eta))
    require(bracket > debt > 0 and 0 < c < 1 and transported >= bracket > debt,
            ("positive transport", body, pair, channel, c, bracket, debt))
    require(actual >= transported > debt,
            ("direct transport", body, pair, channel, actual, transported, debt))
    A = P * ruler
    require(F(g0, A * g0 - a) - F(g0 + 1, A * (g0 + 1) - a)
            == F(a, (A * g0 - a) * (A * (g0 + 1) - a)) > 0,
            ("scale monotonicity", body, pair, channel))
    return margin, body, pair, channel, g0, cell, skeleton, eta, actual, debt


def oriented_channels_below(side: str, threshold: int, oversize: int = 3):
    def in_side(P: int, Q: int) -> bool:
        return side == "all" or (side == "low" and Q < P) or (side == "high" and Q > P)

    return tuple(
        (P, Q)
        for P in range(1, oversize * threshold + 1)
        for Q in range(1, oversize * threshold + 1)
        if P != Q and gcd(P, Q) == 1 and min(P, Q) < threshold
        and F(2, 5) <= F(Q, P) <= F(5, 2) and in_side(P, Q)
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    bodies = T.residual_bodies()
    universal_exceptions = {row[0] for row in T.C2.UNIVERSAL.EXPECTED_EXCEPTIONS}
    require(not (universal_exceptions & set(bodies)),
            ("repeated-level discharge", universal_exceptions & set(bodies)))

    # If m=2k, the boundary D=3k gives kA/(2kL-b), where
    # A=2(b-a)+3b; it decreases in k and begins at (D,m)=(6,4).
    # The odd residue class lies strictly below the common limit A/(2L).
    for body in bodies:
        ruler = 14 * lcm(*body)
        for pair in combinations(range(6), 2):
            a, b = body[pair[0]], body[pair[1]]
            delta = b - a
            A = 2 * delta + 3 * b
            B64 = transport_bound(body, ruler, pair)
            limit = F(A, 2 * ruler)
            require(B64 - limit == F(A * b, 2 * ruler * (4 * ruler - b)) > 0,
                    ("boundary limit", body, pair))
            require(ruler - A > 0, ("odd residue sign", body, pair))
            for k in (2, 3):
                odd = F(k * A + delta + b, 2 * k * ruler + ruler - b)
                require(limit - odd == F(b * (ruler - A),
                        2 * ruler * (2 * k * ruler + ruler - b)) > 0,
                        ("odd residue", body, pair, k))

    eta_guard = max(
        (transport_bound(body, 14 * lcm(*body), pair), body, pair, 14 * lcm(*body))
        for body in combinations(range(1, 15), 6)
        for pair in combinations(range(6), 2)
    )
    require(eta_guard == (ETA_GUARD, H, (0, 5), 168), eta_guard)
    require(MIN_LEVEL - ETA_GUARD == F(631, 165) > 1, eta_guard)

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
    primitive_rows = THREE.ONE.C4.primitive_universe(maximum_bound[0], RATIO_CAP)
    require(len(primitive_rows) == 15995 and max(row[0] for row in primitive_rows) == RATIO_CAP,
            (len(primitive_rows), max(row[0] for row in primitive_rows)))
    phase_zero = tuple((P, Q) for _, phase, _, P, Q in primitive_rows if phase == 0)
    require(phase_zero == ((1, 2), (2, 3), (2, 5), (3, 4)), phase_zero)
    old_vectors = tuple(ZERO_VECTORS[gain] for gain in (F(4, 3), F(3, 2), F(2)))
    all_vectors = old_vectors + (ZERO_VECTORS[F(5, 2)],)
    require((rational_rank(old_vectors), len(old_vectors) - rational_rank(old_vectors)) == (2, 1),
            old_vectors)
    require((rational_rank(all_vectors), len(all_vectors) - rational_rank(all_vectors)) == (3, 1),
            all_vectors)

    closed = []
    traps = []
    witnesses = []
    reasons = {}
    constraints_by_body = {}
    verdict_digest = hashlib.sha256()
    for body in bodies:
        ruler = 14 * lcm(*body)
        constraints, _ = body_constraints(body, ruler, primitive_rows)
        witness, reason = solve_constraints(body, constraints)
        reverse, reverse_reason = solve_constraints(body, constraints, reverse=True)
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
    witness_spans = tuple((body, F(max(witness), min(witness)))
                          for body, witness in witnesses)
    require(max(witness_spans, key=lambda row: row[1]) ==
            ((1, 2, 3, 4, 6, 7), RATIO_CAP), witness_spans)
    require(sum(span == RATIO_CAP for _, span in witness_spans) == 65, witness_spans)

    zero_gains = set(ZERO_VECTORS)
    forced_rows = []
    profile = Counter()
    for body, witness in witnesses:
        forced = tuple(choice[0] for choice in reasons[body][0] if not choice[1])
        require(len(forced) <= 1, ("multiple forced components", body, forced))
        if not forced:
            continue
        component = forced[0]
        zero_edges = []
        new_edges = []
        gains = set()
        old_triangles = []
        for edge in combinations(component, 2):
            if edge not in constraints_by_body[body]:
                continue
            gain = F(max(witness[edge[0]], witness[edge[1]]),
                     min(witness[edge[0]], witness[edge[1]]))
            if gain in zero_gains:
                zero_edges.append(edge)
                gains.add(gain)
            if gain == F(5, 2):
                new_edges.append(edge)
        for triangle in combinations(component, 3):
            edge_set = tuple(combinations(triangle, 2))
            if all(edge in zero_edges for edge in edge_set):
                triangle_gains = {
                    F(max(witness[i], witness[j]), min(witness[i], witness[j]))
                    for i, j in edge_set
                }
                if triangle_gains == {F(4, 3), F(3, 2), F(2)}:
                    old_triangles.append(triangle)
        used = {vertex for edge in zero_edges for vertex in edge}
        components = component_count(used, zero_edges)
        cycle_rank = len(zero_edges) - len(used) + components
        require(gains == zero_gains and rational_rank(tuple(ZERO_VECTORS[g] for g in gains)) == 3,
                ("gain rank", body, gains))
        require(old_triangles and len(new_edges) == 1,
                ("core/bridge multiplicity", body, old_triangles, new_edges))
        bridge = new_edges[0]
        require(component_count(used, tuple(edge for edge in zero_edges if edge != bridge))
                == components + 1, ("new gain is not a bridge", body, bridge, zero_edges))
        profile[(len(component), len(zero_edges), cycle_rank)] += 1
        forced_rows.append((body, component, tuple(zero_edges), tuple(old_triangles), bridge,
                            cycle_rank))
    require(len(forced_rows) == 65 and dict(profile) == FORCED_PROFILE,
            (len(forced_rows), profile))
    forced_digest = hashlib.sha256(repr(tuple(forced_rows)).encode()).hexdigest()

    tail_bodies = {row[0] for row in TAIL_POLICIES}
    require(len(tail_bodies) == 12, tail_bodies)
    located = []
    policies = []
    for body in traps:
        if body in tail_bodies:
            continue
        constraints = constraints_by_body[body]
        require(constraints, ("unexpected empty policy graph", body))
        pair, allowed = min(constraints.items(), key=lambda row: (len(row[1]), row[0]))
        channels = tuple(
            channel for ratio in sorted(allowed)
            for channel in ((ratio.denominator, ratio.numerator),
                            (ratio.numerator, ratio.denominator))
        )
        rows = tuple(located_best(body, pair, channel) for channel in channels)
        require(all(row[0] > 0 for row in rows), ("located policy", body, pair, rows))
        located.extend(rows)
        policies.append((body, pair, len(rows), min(rows)))
    require(len(policies) == 160 and len(located) == 1566,
            (len(policies), len(located)))
    located_digest = hashlib.sha256(repr(tuple(located)).encode()).hexdigest()

    intervals = {"low": (F(2, 5), F(1)), "high": (F(1), F(5, 2)),
                 "all": (F(2, 5), F(5, 2))}
    tail_rows = []
    all_heads = []
    for body, side, pair, start, numerator_bound in TAIL_POLICIES:
        require(not constraints_by_body[body],
                ("tail body unexpectedly constrained", body, constraints_by_body[body]))
        a, b = body[pair[0]], body[pair[1]]
        lo, hi = intervals[side]
        exact_numerator = 2 * max(abs(lo * a - b), abs(hi * a - b))
        require(numerator_bound == exact_numerator,
                ("endpoint numerator", body, side, pair, numerator_bound, exact_numerator))
        margin = T.tail_envelope(body, pair, start, numerator_bound)
        previous = T.tail_envelope(body, pair, start - 1, numerator_bound)
        step = T.tail_step_gain(body, pair, start, numerator_bound)
        ruler = 14 * lcm(*body)
        tail_eta = numerator_bound / 2 / (ruler - F(a, start))
        require(margin > 0 >= previous and step > 0 and tail_eta < 1,
                ("tail threshold", body, side, pair, previous, margin, step, tail_eta))
        channels = oriented_channels_below(side, start)
        oversized = oriented_channels_below(side, start, oversize=5)
        require(channels == oversized,
                ("head exhaustion", body, side, len(channels), len(oversized)))
        rows = tuple(located_best(body, pair, channel) for channel in channels)
        require(all(row[0] > 0 for row in rows),
                ("head failure", body, side, min(rows)))
        all_heads.extend(rows)
        tail_rows.append((body, side, pair, start, numerator_bound, margin, step,
                          tail_eta, len(channels), min(rows)))
    require(tuple(row[8] for row in tail_rows) ==
            (30, 42, 42, 30, 18, 30, 30, 15, 27, 30, 18, 66, 30, 18), tail_rows)
    require(len(all_heads) == 426, len(all_heads))
    head_digest = hashlib.sha256(repr(tuple(all_heads)).encode()).hexdigest()

    semantic_image = (
        tuple(bodies), maximum_bound, tuple(primitive_rows), phase_zero,
        tuple(closed), tuple(traps), tuple(witnesses), witness_spans,
        trap_digest, verdict_digest.hexdigest(), tuple(forced_rows), forced_digest,
        tuple(sorted(profile.items())), tuple(policies), tuple(located), located_digest,
        tuple(tail_rows), tuple(all_heads), head_digest, eta_guard,
    )
    semantic_sha = hashlib.sha256(repr(semantic_image).encode()).hexdigest()
    require(semantic_sha == EXPECTED_SEMANTIC_SHA256,
            ("semantic image changed", semantic_sha, EXPECTED_SEMANTIC_SHA256))

    weakest_located = min(located)
    weakest_head = min(all_heads)
    lines = [
        "LRC14 reflected two-thirds cone exact proof candidate",
        "universe=bodies:3003;arbitrary_level_bank:2442;residual:561;spread:D>=6;cone:3m>=2D",
        f"corner=integer boundary D=6,m=4;ratio_cap=5/2;max_abs_eta={qtext(eta_guard[0])};intermediate_slope_floor={qtext(MIN_LEVEL-eta_guard[0])}",
        f"phase=1/49+c/(PQ),c>=-12/49;constrained_edges={len(threshold_rows)};max_product_bound={qtext(maximum_bound[0])}@{maximum_bound[1]},{maximum_bound[2]};primitive_bank={len(primitive_rows)};cap_endpoint=(2,5)",
        f"phase_zero_channels={phase_zero};old_exponent_rank_nullity=2,1;full_exponent_rank_nullity=3,1",
        f"coarse_csp=closed:{len(closed)};traps:{len(traps)};reverse_order_agrees;largest_witness_span=5/2;forced_full_components={len(forced_rows)}",
        f"coarse_traps={tuple(traps)}",
        f"trap_digest={trap_digest};verdict_digest={verdict_digest.hexdigest()}",
        f"forced_profile={tuple(sorted(profile.items()))};new_5/2_edge_is_bridge_in_all=65;forced_digest={forced_digest}",
        f"located_policies={len(policies)};oriented_controls={len(located)};weakest={qtext(weakest_located[0])}@body={weakest_located[1]},pair={weakest_located[2]},channel={weakest_located[3]},cell={weakest_located[5]};digest={located_digest}",
        *(f"tail_{index+1}=body:{row[0]};side:{row[1]};ordered_pair:{row[2]};start:{row[3]};numerator:{qtext(row[4])};margin:{qtext(row[5])};first_step:{qtext(row[6])};eta_bound:{qtext(row[7])};head_channels:{row[8]};head_weakest:{qtext(row[9][0])}@{row[9][3]},cell={row[9][5]}" for index, row in enumerate(tail_rows)),
        f"head_controls={len(all_heads)};weakest={qtext(weakest_head[0])}@body={weakest_head[1]},pair={weakest_head[2]},channel={weakest_head[3]},cell={weakest_head[5]};head_digest={head_digest}",
        "transport=exact integer corner;oriented eta_g retained;g/(PgL-a) strictly decreases;c_inverse dropped only after positive bracket",
        "tail=PQ>=s^2,Pg>=s;endpoint numerator bounds exact on declared intervals;losses decrease strictly;finite heads 3s-vs-5s audited",
        "conclusion=all reflected residual packets with D>=6 and 3m>=2D close on all 3003 bodies",
        "corollary=assembled reflected certificate-failure wedge is confined to 561 bodies,D>=6,1<=m<2D/3",
        "scope=reflected THM-2941 sufficient family only;physical LRC14 remains open",
        "normal_vs_python_O=BYTE_IDENTICAL",
        f"three_quarter_source_sha256={EXPECTED_THREE_QUARTER_SHA256}",
        f"three_quarter_output_sha256={EXPECTED_THREE_QUARTER_OUTPUT_SHA256}",
        f"three_quarter_semantic_sha256={EXPECTED_THREE_QUARTER_SEMANTIC}",
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
