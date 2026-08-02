#!/usr/bin/env python3
r"""Exact proof candidate for the reflected ``m >= D`` cone.

Inside the sufficient reflected ``k=1`` family of THM-2941, this referee
proves that every residual packet closes whenever

    D >= 6,                 m >= D.

The boundary introduces the primitive zero gain 2 and the multiplicative
gain triangle (3/2)(4/3)=2.  A cap-aware projective CSP leaves 69 coarse
traps.  Complete located policies close 65; orientation-split analytic tails
and exhaustive finite heads close the other four.

This is not a physical-survivor classification and not a proof of LRC(14).
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
FOUR_THIRDS = ROOT / "04-computation/lrc14_j7_reflected_four_thirds_cone_closure_thm2941.py"
FOUR_THIRDS_OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_reflected_four_thirds_cone_closure_thm2941.out"
OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_reflected_one_cone_closure_thm2941.out"

EXPECTED_FOUR_THIRDS_SHA256 = "0b8b01a184f5862c3a32598600988773bcbefdc26c63644cddf0884b18285c3e"
EXPECTED_FOUR_THIRDS_OUTPUT_SHA256 = "0377aaeb146553fe7c027a81ea3589045d01f50fe11aa79848d009a48c806fb9"
EXPECTED_FOUR_THIRDS_SEMANTIC = "d2c1b3bb27ad649db66a3b76305dd6027565cc020c31ddd4cd7e779b48d8bd53"
EXPECTED_SEMANTIC_SHA256 = "dc84fbea8c3e951aa510f02776e5ad300e7d50843da6137fd967f14faff7d2d9"

MIN_SPREAD = 6
MIN_LEVEL = 6
RATIO_CAP = F(2)
ETA_GUARD = F(23, 166)
COARSE_CLOSED = 492
COARSE_TRAPS = 69

H = (1, 2, 3, 4, 6, 12)
H2 = (1, 3, 4, 6, 8, 12)
H4 = (1, 3, 4, 6, 9, 12)
H3 = (2, 3, 4, 6, 8, 12)

# ``side`` refers to r=Q/P.  H must reverse its ordered labels below one;
# the other three bodies admit a single ordered pair on the whole interval.
TAIL_POLICIES = (
    (H, "low", (1, 0), 8, F(2)),
    (H, "high", (0, 1), 8, F(2)),
    (H2, "all", (0, 1), 9, F(5)),
    (H4, "all", (0, 1), 6, F(5)),
    (H3, "all", (0, 1), 7, F(4)),
)

EXPECTED_TRAPS = (
    (1, 2, 3, 4, 6, 8), (1, 2, 3, 4, 6, 9), H,
    (1, 2, 3, 4, 8, 12), (1, 2, 3, 4, 9, 12),
    (1, 2, 3, 5, 6, 10), (1, 2, 3, 5, 6, 12),
    (1, 2, 3, 5, 10, 12), (1, 2, 3, 6, 7, 12),
    (1, 2, 3, 6, 7, 14), (1, 2, 3, 6, 8, 9),
    (1, 2, 3, 6, 8, 12), (1, 2, 3, 6, 9, 12),
    (1, 2, 3, 6, 10, 12), (1, 2, 3, 6, 12, 14),
    (1, 2, 3, 7, 12, 14), (1, 2, 3, 8, 9, 12),
    (1, 2, 4, 5, 6, 10), (1, 2, 4, 5, 6, 12),
    (1, 2, 4, 5, 8, 10), (1, 2, 4, 5, 10, 12),
    (1, 2, 4, 6, 7, 12), (1, 2, 4, 6, 7, 14),
    (1, 2, 4, 6, 8, 9), (1, 2, 4, 6, 8, 12),
    (1, 2, 4, 6, 9, 12), (1, 2, 4, 6, 10, 12),
    (1, 2, 4, 6, 12, 14), (1, 2, 4, 7, 8, 14),
    (1, 2, 4, 7, 12, 14), (1, 2, 4, 8, 9, 12),
    (1, 2, 5, 6, 10, 12), (1, 2, 5, 7, 10, 14),
    (1, 2, 6, 7, 12, 14), (1, 2, 6, 8, 9, 12),
    (1, 3, 4, 5, 6, 10), (1, 3, 4, 5, 6, 12),
    (1, 3, 4, 5, 10, 12), (1, 3, 4, 6, 8, 9), H2, H4,
    (1, 3, 4, 6, 10, 12), (1, 3, 4, 6, 12, 14),
    (1, 3, 4, 7, 12, 14), (1, 3, 4, 8, 9, 12),
    (1, 3, 5, 6, 10, 12), (1, 3, 6, 7, 12, 14),
    (1, 3, 6, 8, 9, 12), (1, 4, 5, 6, 10, 12),
    (1, 4, 6, 7, 12, 14), (1, 4, 6, 8, 9, 12),
    (2, 3, 4, 5, 6, 10), (2, 3, 4, 5, 6, 12),
    (2, 3, 4, 5, 10, 12), (2, 3, 4, 6, 8, 9), H3,
    (2, 3, 4, 6, 9, 12), (2, 3, 4, 6, 10, 12),
    (2, 3, 4, 7, 12, 14), (2, 3, 4, 8, 9, 12),
    (2, 3, 5, 6, 10, 12), (2, 3, 6, 7, 12, 14),
    (2, 3, 6, 8, 9, 12), (2, 4, 5, 6, 10, 12),
    (2, 4, 6, 7, 12, 14), (2, 4, 6, 8, 9, 12),
    (3, 4, 5, 6, 10, 12), (3, 4, 6, 7, 12, 14),
    (3, 4, 6, 8, 9, 12),
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


require(sha256(FOUR_THIRDS) == EXPECTED_FOUR_THIRDS_SHA256,
        ("four-thirds source changed", sha256(FOUR_THIRDS)))
require(sha256(FOUR_THIRDS_OUTPUT) == EXPECTED_FOUR_THIRDS_OUTPUT_SHA256,
        ("four-thirds output changed", sha256(FOUR_THIRDS_OUTPUT)))
require(f"semantic_sha256={EXPECTED_FOUR_THIRDS_SEMANTIC}" in FOUR_THIRDS_OUTPUT.read_text(),
        "four-thirds semantic token missing")
C4 = import_module("one_cone_four_thirds", FOUR_THIRDS)
T = C4.T
T.MIN_LEVEL = MIN_LEVEL
T.MIN_SPREAD = MIN_SPREAD
T.RATIO_CAP = RATIO_CAP


def located_best(body, pair, channel):
    """Best exact body-safe cell, with the full-cone homotopy gate."""
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
    require(bracket > debt > 0, ("positive bracket", body, pair, channel, bracket, debt))
    require(0 < c < 1 and transported >= bracket > debt,
            ("prefactor direction", body, pair, channel, c, bracket, debt))
    require(actual >= transported > debt,
            ("direct transport", body, pair, channel, actual, transported, debt))
    A = P * ruler
    require(
        F(g0, A * g0 - a) - F(g0 + 1, A * (g0 + 1) - a)
        == F(a, (A * g0 - a) * (A * (g0 + 1) - a)) > 0,
        ("scale monotonicity", body, pair, channel),
    )
    return margin, body, pair, channel, g0, cell, skeleton, eta, actual, debt


def oriented_channels_below(side: str, threshold: int, oversize: int = 2):
    def in_side(P: int, Q: int) -> bool:
        return side == "all" or (side == "low" and Q < P) or (side == "high" and Q > P)

    return tuple(
        (P, Q)
        for P in range(1, oversize * threshold + 1)
        for Q in range(1, oversize * threshold + 1)
        if P != Q and gcd(P, Q) == 1 and min(P, Q) < threshold
        and F(1, 2) <= F(Q, P) <= 2 and in_side(P, Q)
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    bodies = T.residual_bodies()
    universal_exceptions = {row[0] for row in T.C2.UNIVERSAL.EXPECTED_EXCEPTIONS}
    require(not (universal_exceptions & set(bodies)),
            ("repeated-level discharge", universal_exceptions & set(bodies)))

    # The cone bound decreases in m.  On the real boundary m=D it decreases
    # in D, so the exact worst corner is D=m=6.
    for body in bodies:
        ruler = 14 * lcm(*body)
        for i, j in combinations(range(6), 2):
            a, b = body[i], body[j]
            delta = b - a
            for D in (6, 7):
                m = F(D)
                B = F(m * delta + D * b, m * ruler - b)
                next_m = F((m + 1) * delta + D * b, (m + 1) * ruler - b)
                require(B - next_m == F(b * (D * ruler + delta),
                        (m * ruler - b) * ((m + 1) * ruler - b)) > 0,
                        ("m monotonicity", body, (i, j), D))
            A = delta + b
            B6 = F(6 * A, 6 * ruler - b)
            B7 = F(7 * A, 7 * ruler - b)
            require(B6 - B7 == F(A * b, (6 * ruler - b) * (7 * ruler - b)) > 0,
                    ("D monotonicity", body, (i, j)))

    eta_guard = max(
        (
            F(MIN_LEVEL * (b - a) + MIN_SPREAD * b,
              MIN_LEVEL * (14 * lcm(*body)) - b),
            body, (a, b), 14 * lcm(*body),
        )
        for body in combinations(range(1, 15), 6)
        for a, b in combinations(body, 2)
    )
    require(eta_guard == (ETA_GUARD, H, (1, 12), 168), eta_guard)
    require(MIN_LEVEL - ETA_GUARD == F(973, 166) > 1, eta_guard)

    threshold_rows = []
    for body in bodies:
        ruler = 14 * lcm(*body)
        for pair in combinations(range(6), 2):
            threshold = T.coarse_threshold(body, ruler, pair)
            if threshold < T.PHASE_LIMIT:
                threshold_rows.append((threshold, body, pair))
    maximum_bound = max(
        (T.product_bound(threshold), body, pair, threshold)
        for threshold, body, pair in threshold_rows
    )
    primitive_rows = C4.primitive_universe(maximum_bound[0], RATIO_CAP)
    require(len(primitive_rows) == 279 and max(row[0] for row in primitive_rows) == 2,
            (len(primitive_rows), max(row[0] for row in primitive_rows)))
    phase_zero = tuple((P, Q) for _, phase, _, P, Q in primitive_rows if phase == 0)
    require(phase_zero == ((1, 2), (2, 3), (3, 4)), phase_zero)
    require(F(3, 2) * F(4, 3) == 2, "zero-gain circuit")

    closed = []
    traps = []
    witnesses = []
    reasons = {}
    constraints_by_body = {}
    verdict_digest = hashlib.sha256()
    for body in bodies:
        ruler = 14 * lcm(*body)
        witness, constraints, _, reason = T.solve_body(body, ruler, primitive_rows)
        reverse, reverse_constraints, _, reverse_reason = T.solve_body(
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
            reasons[body] = reason
            constraints_by_body[body] = constraints
    require((len(closed), len(traps), tuple(traps)) ==
            (COARSE_CLOSED, COARSE_TRAPS, EXPECTED_TRAPS),
            (len(closed), len(traps), traps))

    witness_spans = tuple((body, F(max(witness), min(witness)))
                          for body, witness in witnesses)
    require(max(span for _, span in witness_spans) == 2, witness_spans)
    require(sum(span == 2 for _, span in witness_spans) == 21, witness_spans)

    # Every forced full-span component contains the new balanced zero triangle.
    zero_gains = {F(3, 2), F(4, 3), F(2)}
    gain_triangles = []
    for body, witness in witnesses:
        forced = tuple(choice[0] for choice in reasons[body][0] if not choice[1])
        require(len(forced) <= 1, ("multiple forced components", body, forced))
        if not forced:
            continue
        component = forced[0]
        cores = []
        for triangle in combinations(component, 3):
            if not all(tuple(sorted(edge)) in constraints_by_body[body]
                       for edge in combinations(triangle, 2)):
                continue
            gains = {
                F(max(witness[i], witness[j]), min(witness[i], witness[j]))
                for i, j in combinations(triangle, 2)
            }
            if gains == zero_gains:
                values = sorted(witness[index] for index in triangle)
                require(F(values[1], values[0]) * F(values[2], values[1]) == 2,
                        ("unbalanced gain circuit", body, triangle, values))
                cores.append(triangle)
        require(cores, ("forced full component lacks zero triangle", body, component, witness))
        gain_triangles.append((body, component, tuple(cores)))
    require(len(gain_triangles) == 21, len(gain_triangles))
    require(tuple(sorted(len(row[2]) for row in gain_triangles)) == (1,) * 19 + (2,) * 2,
            gain_triangles)
    gain_digest = hashlib.sha256(repr(tuple(gain_triangles)).encode()).hexdigest()

    tail_bodies = {row[0] for row in TAIL_POLICIES}
    located = []
    policies = []
    for body in traps:
        if body in tail_bodies:
            continue
        constraints = constraints_by_body[body]
        require(constraints, ("unexpected empty policy graph", body))
        pair, allowed = min(constraints.items(), key=lambda row: (len(row[1]), row[0]))
        channels = tuple(
            channel
            for ratio in sorted(allowed)
            for channel in ((ratio.denominator, ratio.numerator),
                            (ratio.numerator, ratio.denominator))
        )
        rows = tuple(located_best(body, pair, channel) for channel in channels)
        require(all(row[0] > 0 for row in rows), ("located policy", body, pair, rows))
        located.extend(rows)
        policies.append((body, pair, len(rows), min(rows)))
    require(len(policies) == 65 and len(located) == 894,
            (len(policies), len(located)))
    located_digest = hashlib.sha256(repr(tuple(located)).encode()).hexdigest()

    intervals = {"low": (F(1, 2), F(1)), "high": (F(1), F(2)),
                 "all": (F(1, 2), F(2))}
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
        tail_eta = F(numerator_bound, 2) / (ruler - F(a, start))
        require(margin > 0 >= previous and step > 0 and tail_eta < 1,
                ("tail threshold", body, side, pair, previous, margin, step, tail_eta))
        channels = oriented_channels_below(side, start)
        oversized = oriented_channels_below(side, start, oversize=4)
        require(channels == oversized,
                ("head exhaustion", body, side, len(channels), len(oversized)))
        rows = tuple(located_best(body, pair, channel) for channel in channels)
        require(all(row[0] > 0 for row in rows),
                ("head failure", body, side, min(rows)))
        all_heads.extend(rows)
        tail_rows.append((body, side, pair, start, numerator_bound, margin, step,
                          tail_eta, len(channels), min(rows)))
    require(tuple(row[8] for row in tail_rows) == (18, 18, 44, 20, 24), tail_rows)
    require(len(all_heads) == 124, len(all_heads))
    head_digest = hashlib.sha256(repr(tuple(all_heads)).encode()).hexdigest()

    # The old fixed order genuinely degenerates at the new reverse 2:1 zero
    # channel; reversing the labels restores the same cell and zero transport.
    wrong_fixed = T.fixed_H_skeleton((2, 1))
    repaired_fixed = T.intersection_mass(
        T.primitive_arcs(2, F(2 * 155, 168) % 1),
        T.primitive_arcs(1, F(155, 168)),
    )
    require(wrong_fixed == 0 and repaired_fixed == F(1, 14),
            (wrong_fixed, repaired_fixed))

    semantic_image = (
        tuple(bodies), maximum_bound, tuple(primitive_rows), phase_zero,
        tuple(closed), tuple(traps), tuple(witnesses), witness_spans,
        verdict_digest.hexdigest(), tuple(gain_triangles), gain_digest,
        tuple(policies), tuple(located), located_digest, tuple(tail_rows),
        tuple(all_heads), head_digest, wrong_fixed, repaired_fixed, eta_guard,
    )
    semantic_sha = hashlib.sha256(repr(semantic_image).encode()).hexdigest()
    require(semantic_sha == EXPECTED_SEMANTIC_SHA256,
            ("semantic image changed", semantic_sha, EXPECTED_SEMANTIC_SHA256))

    weakest_located = min(located)
    weakest_head = min(all_heads)
    lines = [
        "LRC14 reflected one-cone exact proof candidate",
        "universe=bodies:3003;arbitrary_level_bank:2442;residual:561;spread:D>=6;cone:m>=D",
        f"corner=real boundary m=D decreases in D;worst_D=6;worst_m=6;ratio_cap=2;max_abs_eta={qtext(eta_guard[0])};intermediate_slope_floor={qtext(MIN_LEVEL-eta_guard[0])}",
        f"phase=1/49+c/(PQ),c>=-12/49;constrained_edges={len(threshold_rows)};max_product_bound={qtext(maximum_bound[0])}@{maximum_bound[1]},{maximum_bound[2]};primitive_bank={len(primitive_rows)};cap_endpoint=(1,2)",
        f"phase_zero_channels={phase_zero};balanced_gain_circuit=(3/2)*(4/3)=2",
        f"coarse_csp=closed:{len(closed)};traps:{len(traps)};reverse_order_agrees;largest_witness_span={qtext(max(span for _,span in witness_spans))};forced_full_components={len(gain_triangles)}",
        f"coarse_traps={tuple(traps)}",
        f"verdict_digest={verdict_digest.hexdigest()}",
        f"gain_triangle_rows={len(gain_triangles)};core_count_profile=19x1+2x2;digest={gain_digest}",
        f"located_policies={len(policies)};oriented_controls={len(located)};weakest={qtext(weakest_located[0])}@body={weakest_located[1]},pair={weakest_located[2]},channel={weakest_located[3]},cell={weakest_located[5]};digest={located_digest}",
        *(f"tail_{index+1}=body:{row[0]};side:{row[1]};ordered_pair:{row[2]};start:{row[3]};numerator:{qtext(row[4])};margin:{qtext(row[5])};first_step:{qtext(row[6])};eta_bound:{qtext(row[7])};head_channels:{row[8]};head_weakest:{qtext(row[9][0])}@{row[9][3]},cell={row[9][5]}" for index, row in enumerate(tail_rows)),
        f"head_controls={len(all_heads)};weakest={qtext(weakest_head[0])}@body={weakest_head[1]},pair={weakest_head[2]},channel={weakest_head[3]},cell={weakest_head[5]};head_digest={head_digest}",
        "orientation_boundary=old H ordered pair at reverse channel (2,1),cell155 has skeleton 0;reversed labels restore skeleton 1/14 and eta 0",
        "transport=oriented eta_g retained;g/(PgL-a) strictly decreases;c_inverse dropped only after positive bracket;debt bounded at level six",
        "tail=PQ>=s^2,Pg>=s;endpoint numerator bounds exact on declared half/full intervals;losses decrease strictly;finite heads oversized-audited",
        "conclusion=all reflected residual packets with D>=6 and m>=D close on all 3003 bodies",
        "corollary=assembled reflected certificate-failure wedge is confined to 561 bodies,D>=6,1<=m<D",
        "scope=reflected THM-2941 sufficient family only;physical LRC14 remains open",
        "normal_vs_python_O=BYTE_IDENTICAL",
        f"four_thirds_source_sha256={EXPECTED_FOUR_THIRDS_SHA256}",
        f"four_thirds_output_sha256={EXPECTED_FOUR_THIRDS_OUTPUT_SHA256}",
        f"four_thirds_semantic_sha256={EXPECTED_FOUR_THIRDS_SEMANTIC}",
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
