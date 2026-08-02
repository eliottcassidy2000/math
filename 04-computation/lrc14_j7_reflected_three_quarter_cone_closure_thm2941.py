#!/usr/bin/env python3
r"""Exact proof candidate for the reflected ``4m >= 3D`` cone.

Inside the sufficient reflected ``k=1`` family of THM-2941, this referee
closes every residual packet with

    D >= 6,                 4m >= 3D.

The exact integer transport corner is (D,m)=(8,6), while singleton debt is
bounded at the independent minimum level five.  A cap-aware CSP closes 459
bodies; 806 complete located controls and 276 finite tail-head controls close
the remaining 102 coarse traps.

This is not a physical-survivor census and not a proof of LRC(14).
"""

from __future__ import annotations

import argparse
import hashlib
from fractions import Fraction as F
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from math import gcd, lcm
from pathlib import Path


HERE = Path(__file__).resolve()
ROOT = next(parent for parent in HERE.parents if (parent / "00-navigation").is_dir())
ONE_CONE = ROOT / "04-computation/lrc14_j7_reflected_one_cone_closure_thm2941.py"
ONE_CONE_OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_reflected_one_cone_closure_thm2941.out"
OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_reflected_three_quarter_cone_closure_thm2941.out"

EXPECTED_ONE_CONE_SHA256 = "800395ae242860094fed3db9638a93ebc2faba7973558a3eaa51f3af62145200"
EXPECTED_ONE_CONE_OUTPUT_SHA256 = "e63ee74f42d8f213196cf2907bbc25dc20b93a3e64b8a40aa85d5d81c8b11ee6"
EXPECTED_ONE_CONE_SEMANTIC = "dc84fbea8c3e951aa510f02776e5ad300e7d50843da6137fd967f14faff7d2d9"
EXPECTED_TRAP_DIGEST = "cff2ad5603eb88999cc9a5bebb8cdc4a6e3a2cb29d52fab433a935ed4453a971"
EXPECTED_SEMANTIC_SHA256 = "05a900c80283d5bc9a7b01c1c4ad045b889aebe2e4d5798eaa66cb8907a7ae9f"

MIN_SPREAD = 6
MIN_LEVEL = 5
RATIO_CAP = F(7, 3)
ETA_GUARD = F(27, 166)
COARSE_CLOSED = 459
COARSE_TRAPS = 102

H = (1, 2, 3, 4, 6, 12)
H2 = (1, 3, 4, 6, 8, 12)

TAIL_POLICIES = (
    ((1, 2, 3, 4, 6, 8), "all", (0, 1), 6, F(22, 7)),
    (H, "low", (1, 0), 9, F(2)),
    (H, "high", (0, 1), 9, F(2)),
    ((1, 2, 3, 4, 8, 12), "all", (0, 1), 6, F(22, 7)),
    ((1, 2, 3, 6, 8, 12), "all", (0, 1), 6, F(22, 7)),
    ((1, 2, 4, 6, 8, 12), "all", (0, 1), 6, F(22, 7)),
    (H2, "low", (2, 1), 5, F(18, 7)),
    (H2, "high", (0, 1), 7, F(4)),
    ((1, 3, 4, 6, 9, 12), "all", (0, 1), 6, F(36, 7)),
    ((2, 3, 4, 6, 8, 12), "all", (0, 1), 8, F(30, 7)),
    ((2, 3, 4, 6, 9, 12), "all", (0, 1), 5, F(30, 7)),
)


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


require(sha256(ONE_CONE) == EXPECTED_ONE_CONE_SHA256,
        ("one-cone source changed", sha256(ONE_CONE)))
require(sha256(ONE_CONE_OUTPUT) == EXPECTED_ONE_CONE_OUTPUT_SHA256,
        ("one-cone output changed", sha256(ONE_CONE_OUTPUT)))
require(f"semantic_sha256={EXPECTED_ONE_CONE_SEMANTIC}" in ONE_CONE_OUTPUT.read_text(),
        "one-cone semantic token missing")
ONE = import_module("three_quarter_one_cone", ONE_CONE)
T = ONE.T
T.MIN_LEVEL = MIN_LEVEL
T.MIN_SPREAD = MIN_SPREAD
T.RATIO_CAP = RATIO_CAP


def transport_bound(body, ruler, pair) -> F:
    """Exact integer-cone maximum, attained at ``(D,m)=(8,6)``."""
    a, b = body[pair[0]], body[pair[1]]
    return F(6 * (b - a) + 8 * b, 6 * ruler - b)


def coarse_threshold(body, ruler, pair) -> F:
    # Combining the exact transport corner with the independent level-five
    # debt corner is an over-bound, hence safe even though no packet need
    # realize both corners simultaneously.
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


def solve_body(body, ruler, primitive_rows, reverse=False):
    """Cap-parametric CSP solver without monkeypatching a threshold global."""
    constraints, thresholds = body_constraints(body, ruler, primitive_rows)
    if any(not allowed for allowed in constraints.values()):
        return None, constraints, thresholds, ("empty edge",)
    choices = []
    component_rows = []
    nodes = 0
    forced_full = 0
    for component in T.connected_components(constraints):
        short, full, count = T.component_types(component, constraints, reverse=reverse)
        nodes += count
        if short is None and full is None:
            return None, constraints, thresholds, ("empty component", component, nodes)
        if short is None:
            forced_full += 1
            component_rows.append(full)
        else:
            component_rows.append(short)
        choices.append((component, short is not None, full is not None, count))
    if forced_full > 1:
        return None, constraints, thresholds, ("two forced full components", tuple(choices), nodes)
    witness = T.scale_components(tuple(component_rows))
    for pair, allowed in constraints.items():
        i, j = pair
        ratio = F(max(witness[i], witness[j]), min(witness[i], witness[j]))
        require(ratio in allowed, (body, witness, pair, ratio, allowed))
    return witness, constraints, thresholds, (tuple(choices), nodes)


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
    require(
        F(g0, A * g0 - a) - F(g0 + 1, A * (g0 + 1) - a)
        == F(a, (A * g0 - a) * (A * (g0 + 1) - a)) > 0,
        ("scale monotonicity", body, pair, channel),
    )
    return margin, body, pair, channel, g0, cell, skeleton, eta, actual, debt


def oriented_channels_below(side: str, threshold: int, oversize: int = 3):
    def in_side(P: int, Q: int) -> bool:
        return side == "all" or (side == "low" and Q < P) or (side == "high" and Q > P)

    return tuple(
        (P, Q)
        for P in range(1, oversize * threshold + 1)
        for Q in range(1, oversize * threshold + 1)
        if P != Q and gcd(P, Q) == 1 and min(P, Q) < threshold
        and F(3, 7) <= F(Q, P) <= F(7, 3) and in_side(P, Q)
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    bodies = T.residual_bodies()
    universal_exceptions = {row[0] for row in T.C2.UNIVERSAL.EXPECTED_EXCEPTIONS}
    require(not (universal_exceptions & set(bodies)),
            ("repeated-level discharge", universal_exceptions & set(bodies)))

    # Write A=3(b-a)+4b.  For m=3k the boundary D=4k gives
    # kA/(3kL-b), decreasing in k.  The other two residue classes lie below
    # the common limit A/(3L).  The first exact-boundary point (8,6) exceeds
    # that limit and also the two smaller admissible corners (6,5),(9,7).
    for body in bodies:
        ruler = 14 * lcm(*body)
        for pair in combinations(range(6), 2):
            a, b = body[pair[0]], body[pair[1]]
            delta = b - a
            A = 3 * delta + 4 * b
            B86 = transport_bound(body, ruler, pair)
            B65 = F(5 * delta + 6 * b, 5 * ruler - b)
            B97 = F(7 * delta + 9 * b, 7 * ruler - b)
            limit = F(A, 3 * ruler)
            require(B86 - B65 == F(b * (4 * ruler + a - 3 * b),
                    (6 * ruler - b) * (5 * ruler - b)) > 0,
                    ("corner 6,5", body, pair))
            require(B86 - B97 == F(b * (2 * ruler + delta + b),
                    (6 * ruler - b) * (7 * ruler - b)) > 0,
                    ("corner 9,7", body, pair))
            require(B86 - limit == F(A * b, 3 * ruler * (6 * ruler - b)) > 0,
                    ("boundary limit", body, pair))
            require(ruler - 7 * b + 3 * a > 0 and 2 * ruler - A > 0,
                    ("off-residue signs", body, pair))
            for k in (2, 3):
                residue_one = F(k * A + delta + b, 3 * k * ruler + ruler - b)
                residue_two = F(k * A + 2 * (delta + b),
                                3 * k * ruler + 2 * ruler - b)
                require(limit - residue_one == F(b * (ruler - 7 * b + 3 * a),
                        3 * ruler * (3 * k * ruler + ruler - b)) > 0,
                        ("residue one", body, pair, k))
                require(limit - residue_two == F(b * (2 * ruler - A),
                        3 * ruler * (3 * k * ruler + 2 * ruler - b)) > 0,
                        ("residue two", body, pair, k))

    eta_guard = max(
        (transport_bound(body, 14 * lcm(*body), pair), body, pair, 14 * lcm(*body))
        for body in combinations(range(1, 15), 6)
        for pair in combinations(range(6), 2)
    )
    require(eta_guard == (ETA_GUARD, H, (0, 5), 168), eta_guard)
    require(MIN_LEVEL - ETA_GUARD == F(803, 166) > 1, eta_guard)

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
    primitive_rows = ONE.C4.primitive_universe(maximum_bound[0], RATIO_CAP)
    require(len(primitive_rows) == 9857 and max(row[0] for row in primitive_rows) == RATIO_CAP,
            (len(primitive_rows), max(row[0] for row in primitive_rows)))
    phase_zero = tuple((P, Q) for _, phase, _, P, Q in primitive_rows if phase == 0)
    require(phase_zero == ((1, 2), (2, 3), (3, 4)), phase_zero)

    closed = []
    traps = []
    witnesses = []
    constraints_by_body = {}
    verdict_digest = hashlib.sha256()
    for body in bodies:
        ruler = 14 * lcm(*body)
        witness, constraints, _, reason = solve_body(body, ruler, primitive_rows)
        reverse, reverse_constraints, _, reverse_reason = solve_body(
            body, ruler, primitive_rows, reverse=True
        )
        require((witness is None) == (reverse is None),
                ("search-order mismatch", body, reason, reverse_reason))
        require(constraints == reverse_constraints, ("constraint mismatch", body))
        verdict_digest.update(repr((body, witness is not None, reason, reverse_reason)).encode())
        if witness is None:
            closed.append(body)
        else:
            traps.append(body)
            witnesses.append((body, witness))
            constraints_by_body[body] = constraints
    trap_digest = hashlib.sha256(repr(tuple(traps)).encode()).hexdigest()
    require((len(closed), len(traps)) == (COARSE_CLOSED, COARSE_TRAPS),
            (len(closed), len(traps)))
    require(trap_digest == EXPECTED_TRAP_DIGEST,
            ("trap census changed", trap_digest, EXPECTED_TRAP_DIGEST, traps))
    witness_spans = tuple((body, F(max(witness), min(witness)))
                          for body, witness in witnesses)
    require(max(witness_spans, key=lambda row: row[1]) ==
            ((2, 3, 4, 5, 10, 12), F(128, 55)), witness_spans)
    require(all(span < RATIO_CAP for _, span in witness_spans), witness_spans)

    tail_bodies = {row[0] for row in TAIL_POLICIES}
    require(len(tail_bodies) == 9, tail_bodies)
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
    require(len(policies) == 93 and len(located) == 806,
            (len(policies), len(located)))
    located_digest = hashlib.sha256(repr(tuple(located)).encode()).hexdigest()

    intervals = {"low": (F(3, 7), F(1)), "high": (F(1), F(7, 3)),
                 "all": (F(3, 7), F(7, 3))}
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
            (26, 29, 29, 26, 26, 26, 8, 16, 26, 48, 16), tail_rows)
    require(len(all_heads) == 276, len(all_heads))
    head_digest = hashlib.sha256(repr(tuple(all_heads)).encode()).hexdigest()

    semantic_image = (
        tuple(bodies), maximum_bound, tuple(primitive_rows), phase_zero,
        tuple(closed), tuple(traps), tuple(witnesses), witness_spans,
        trap_digest, verdict_digest.hexdigest(), tuple(policies), tuple(located),
        located_digest, tuple(tail_rows), tuple(all_heads), head_digest, eta_guard,
    )
    semantic_sha = hashlib.sha256(repr(semantic_image).encode()).hexdigest()
    require(semantic_sha == EXPECTED_SEMANTIC_SHA256,
            ("semantic image changed", semantic_sha, EXPECTED_SEMANTIC_SHA256))

    weakest_located = min(located)
    weakest_head = min(all_heads)
    lines = [
        "LRC14 reflected three-quarter cone exact proof candidate",
        "universe=bodies:3003;arbitrary_level_bank:2442;residual:561;spread:D>=6;cone:4m>=3D",
        f"corner=integer boundary D=8,m=6;debt_level=5;ratio_cap=7/3;max_abs_eta={qtext(eta_guard[0])};intermediate_slope_floor={qtext(MIN_LEVEL-eta_guard[0])}",
        f"phase=1/49+c/(PQ),c>=-12/49;constrained_edges={len(threshold_rows)};max_product_bound={qtext(maximum_bound[0])}@{maximum_bound[1]},{maximum_bound[2]};primitive_bank={len(primitive_rows)};cap_endpoint=(3,7)",
        f"phase_zero_channels={phase_zero};no_new_zero_gain_beyond_cap2",
        f"coarse_csp=closed:{len(closed)};traps:{len(traps)};reverse_order_agrees;largest_witness_span={qtext(max(span for _,span in witness_spans))}<7/3;forced_full_components=0",
        f"coarse_traps={tuple(traps)}",
        f"trap_digest={trap_digest};verdict_digest={verdict_digest.hexdigest()}",
        f"located_policies={len(policies)};oriented_controls={len(located)};weakest={qtext(weakest_located[0])}@body={weakest_located[1]},pair={weakest_located[2]},channel={weakest_located[3]},cell={weakest_located[5]};digest={located_digest}",
        *(f"tail_{index+1}=body:{row[0]};side:{row[1]};ordered_pair:{row[2]};start:{row[3]};numerator:{qtext(row[4])};margin:{qtext(row[5])};first_step:{qtext(row[6])};eta_bound:{qtext(row[7])};head_channels:{row[8]};head_weakest:{qtext(row[9][0])}@{row[9][3]},cell={row[9][5]}" for index, row in enumerate(tail_rows)),
        f"head_controls={len(all_heads)};weakest={qtext(weakest_head[0])}@body={weakest_head[1]},pair={weakest_head[2]},channel={weakest_head[3]},cell={weakest_head[5]};head_digest={head_digest}",
        "transport=exact integer corner separated from debt corner;oriented eta_g retained;g/(PgL-a) strictly decreases;c_inverse dropped only after positive bracket",
        "tail=PQ>=s^2,Pg>=s;endpoint numerator bounds exact on declared intervals;losses decrease strictly;finite heads 3s-vs-5s audited",
        "conclusion=all reflected residual packets with D>=6 and 4m>=3D close on all 3003 bodies",
        "corollary=assembled reflected certificate-failure wedge is confined to 561 bodies,D>=6,1<=m<3D/4",
        "scope=reflected THM-2941 sufficient family only;physical LRC14 remains open",
        "normal_vs_python_O=BYTE_IDENTICAL",
        f"one_cone_source_sha256={EXPECTED_ONE_CONE_SHA256}",
        f"one_cone_output_sha256={EXPECTED_ONE_CONE_OUTPUT_SHA256}",
        f"one_cone_semantic_sha256={EXPECTED_ONE_CONE_SEMANTIC}",
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
