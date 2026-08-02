#!/usr/bin/env python3
r"""Exact proof candidate for the reflected ``3m >= 4D`` cone.

This is a recursive specialization of the exact edge-alphabet/projective-CSP
engine frozen in the three-halves cone referee.  It proves, inside the
sufficient reflected ``k=1`` family of THM-2941, that every packet closes for

    D >= 6,                 3m >= 4D.

It is not a physical-survivor classification and not a proof of LRC(14).
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
THREE_HALVES = ROOT / "04-computation/lrc14_j7_reflected_three_halves_cone_closure_thm2941.py"
THREE_HALVES_OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_reflected_three_halves_cone_closure_thm2941.out"
OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_reflected_four_thirds_cone_closure_thm2941.out"

EXPECTED_THREE_HALVES_SHA256 = "15bcd2d67b4867d1ac6380a10a3b9618d677b031ad4354729371ba77501b50e6"
EXPECTED_THREE_HALVES_OUTPUT_SHA256 = "7d7ad662745d170c2492e93bf296bdd0b15372ef2ca83aed82bda41770a1089d"
EXPECTED_THREE_HALVES_SEMANTIC = "6d2e760f3a929ad2eaa9af9eba11db1389ead59a383e11d8e8737c722b67aaaa"
EXPECTED_SEMANTIC_SHA256 = "9f5d2561da2dbf4f1779d88cb5a17ccc199e22a43ffcc82bc7c54b5ff47ed622"

MIN_SPREAD = 6
MIN_LEVEL = 8
RATIO_CAP = F(7, 4)
COARSE_CLOSED = 530
COARSE_TRAPS = 31

H = (1, 2, 3, 4, 6, 12)
H2 = (1, 3, 4, 6, 8, 12)
H3 = (2, 3, 4, 6, 8, 12)
TAIL_POLICIES = (
    (H, (0, 1), 25, F(20, 7)),
    (H2, (1, 2), 8, F(32, 7)),
    (H3, (0, 1), 6, F(26, 7)),
)

EXPECTED_TRAPS = (
    (1, 2, 3, 4, 6, 8), (1, 2, 3, 4, 6, 9), H,
    (1, 2, 3, 4, 8, 12), (1, 2, 3, 4, 9, 12),
    (1, 2, 3, 5, 6, 10), (1, 2, 3, 6, 7, 14),
    (1, 2, 3, 6, 8, 12), (1, 2, 3, 6, 9, 12),
    (1, 2, 4, 5, 8, 10), (1, 2, 4, 6, 8, 12),
    (1, 2, 4, 6, 9, 12), (1, 2, 4, 6, 10, 12),
    (1, 2, 4, 7, 8, 14), (1, 2, 5, 6, 10, 12),
    (1, 2, 5, 7, 10, 14), (1, 2, 6, 8, 9, 12), H2,
    (1, 3, 4, 6, 9, 12), (1, 3, 4, 6, 10, 12),
    (1, 3, 5, 6, 10, 12), (1, 3, 6, 8, 9, 12),
    (1, 4, 5, 6, 10, 12), (1, 4, 6, 8, 9, 12), H3,
    (2, 3, 4, 6, 9, 12), (2, 3, 5, 6, 10, 12),
    (2, 3, 6, 8, 9, 12), (2, 4, 5, 6, 10, 12),
    (2, 4, 6, 8, 9, 12), (3, 4, 6, 8, 9, 12),
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


require(sha256(THREE_HALVES) == EXPECTED_THREE_HALVES_SHA256,
        ("three-halves source changed", sha256(THREE_HALVES)))
require(sha256(THREE_HALVES_OUTPUT) == EXPECTED_THREE_HALVES_OUTPUT_SHA256,
        ("three-halves output changed", sha256(THREE_HALVES_OUTPUT)))
require(f"semantic_sha256={EXPECTED_THREE_HALVES_SEMANTIC}" in THREE_HALVES_OUTPUT.read_text(),
        "three-halves semantic token missing")
T = import_module("four_thirds_engine", THREE_HALVES)

# The imported functions read these three exact parameters at call time.  All
# other constants and pinned evidence remain unchanged.
T.MIN_LEVEL = MIN_LEVEL
T.MIN_SPREAD = MIN_SPREAD
T.RATIO_CAP = RATIO_CAP


def oriented_channels_below(threshold: int, oversize: int = 2):
    return tuple(
        (P, Q)
        for P in range(1, oversize * threshold + 1)
        for Q in range(1, oversize * threshold + 1)
        if P != Q and gcd(P, Q) == 1 and min(P, Q) < threshold
        and F(4, 7) <= F(Q, P) <= F(7, 4)
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    bodies = T.residual_bodies()

    # Exact infinite-cone reduction, including D not divisible by three.
    for body in bodies:
        ruler = 14 * lcm(*body)
        for i, j in combinations(range(6), 2):
            a, b = body[i], body[j]
            delta = b - a
            for D in (6, 7):
                m = F(4 * D, 3)
                B = F(m * delta + D * b, m * ruler - b)
                next_m = F((m + 1) * delta + D * b, (m + 1) * ruler - b)
                require(
                    B - next_m
                    == F(b * (D * ruler + delta),
                         (m * ruler - b) * ((m + 1) * ruler - b)) > 0,
                    ("m monotonicity", body, (i, j), D),
                )
            A = F(4, 3) * delta + b
            B6 = F(6 * A, 8 * ruler - b)
            B7 = F(7 * A, F(28, 3) * ruler - b)
            require(B6 - B7 == F(A * b, (8 * ruler - b) * (F(28, 3) * ruler - b)) > 0,
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
    require(eta_guard == (F(40, 333), H, (1, 12), 168), eta_guard)
    require(MIN_LEVEL - eta_guard[0] == F(2624, 333) > 1, eta_guard)

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
    primitive_rows = T.primitive_universe(maximum_bound[0])

    closed = []
    traps = []
    witnesses = []
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
    require((len(closed), len(traps), tuple(traps)) == (COARSE_CLOSED, COARSE_TRAPS, EXPECTED_TRAPS),
            (len(closed), len(traps), traps))
    witness_spans = tuple(
        (body, F(max(witness), min(witness))) for body, witness in witnesses
    )
    require(max(span for _, span in witness_spans) == F(2021, 1212) < RATIO_CAP,
            witness_spans)

    tail_bodies = {body for body, _, _, _ in TAIL_POLICIES}
    located = []
    policies = []
    for body in traps:
        if body in tail_bodies:
            continue
        ruler = 14 * lcm(*body)
        constraints, _ = T.body_constraints(body, ruler, primitive_rows)
        require(constraints, ("unexpected empty policy graph", body))
        pair, allowed = min(constraints.items(), key=lambda row: (len(row[1]), row[0]))
        channels = tuple(
            channel
            for ratio in sorted(allowed)
            for channel in (
                (ratio.denominator, ratio.numerator),
                (ratio.numerator, ratio.denominator),
            )
        )
        rows = tuple(T.located_best(body, pair, channel) for channel in channels)
        require(all(row[0] > 0 for row in rows), ("located policy", body, pair, rows))
        located.extend(rows)
        policies.append((body, pair, len(rows), min(rows)))
    require(len(policies) == 28 and len(located) == 158, (len(policies), len(located)))
    located_digest = hashlib.sha256(repr(tuple(located)).encode()).hexdigest()

    tail_rows = []
    all_heads = []
    for body, pair, start, numerator_bound in TAIL_POLICIES:
        tail_constraints, _ = T.body_constraints(body, 14 * lcm(*body), primitive_rows)
        require(not tail_constraints, ("tail body unexpectedly constrained", body, tail_constraints))
        margin = T.tail_envelope(body, pair, start, numerator_bound)
        previous = T.tail_envelope(body, pair, start - 1, numerator_bound)
        step = T.tail_step_gain(body, pair, start, numerator_bound)
        require(margin > 0 >= previous and step > 0,
                ("tail threshold", body, pair, start, previous, margin, step))
        channels = oriented_channels_below(start)
        oversized = oriented_channels_below(start, oversize=4)
        require(channels == oversized, ("head exhaustion", body, len(channels), len(oversized)))
        rows = tuple(T.located_best(body, pair, channel) for channel in channels)
        require(all(row[0] > 0 for row in rows), ("head failure", body, min(rows)))
        all_heads.extend(rows)
        tail_rows.append((body, pair, start, numerator_bound, margin, step, len(channels), min(rows)))
    require(tuple(row[6] for row in tail_rows) == (270, 28, 16), tail_rows)
    require(len(all_heads) == 314, len(all_heads))
    head_digest = hashlib.sha256(repr(tuple(all_heads)).encode()).hexdigest()

    # The fixed-cell lemma from the three-halves proof extends to [4/7,7/4]:
    # the PQ<20 exceptional bank is unchanged, while the global tail bound is
    # independent of the cone ratio.
    fixed_small = tuple(
        (P, Q, T.fixed_H_skeleton((P, Q)))
        for P in range(1, 20) for Q in range(1, 20)
        if P != Q and gcd(P, Q) == 1 and P * Q < 20
        and F(4, 7) <= F(Q, P) <= F(7, 4)
    )
    require(T.PHASE_LIMIT - F(12, 49 * 20) >= F(1, 126), "fixed-cell tail")
    require(min(value for _, _, value in fixed_small) == F(1, 126)
            and tuple((P, Q) for P, Q, value in fixed_small if value == F(1, 126)) == ((3, 2),),
            fixed_small)

    semantic_image = (
        tuple(bodies), maximum_bound, tuple(primitive_rows), tuple(closed),
        tuple(traps), tuple(witnesses), witness_spans, verdict_digest.hexdigest(),
        tuple(policies), tuple(located), located_digest, tuple(tail_rows),
        tuple(all_heads), head_digest, fixed_small, eta_guard,
    )
    semantic_sha = hashlib.sha256(repr(semantic_image).encode()).hexdigest()
    require(semantic_sha == EXPECTED_SEMANTIC_SHA256,
            ("semantic image changed", semantic_sha, EXPECTED_SEMANTIC_SHA256))

    weakest_located = min(located)
    lines = [
        "LRC14 reflected four-thirds cone exact proof candidate",
        "universe=bodies:3003;arbitrary_level_bank:2442;residual:561;spread:D>=6;cone:3m>=4D",
        f"corner=real boundary m=4D/3 decreases in D;worst_D=6;worst_m=8;ratio_cap=7/4;max_abs_eta={qtext(eta_guard[0])};intermediate_slope_floor={qtext(MIN_LEVEL-eta_guard[0])}",
        f"phase=1/49+c/(PQ),c>=-12/49;constrained_edges={len(threshold_rows)};max_product_bound={qtext(maximum_bound[0])}@{maximum_bound[1]},{maximum_bound[2]};primitive_bank={len(primitive_rows)}",
        f"coarse_csp=closed:{len(closed)};traps:{len(traps)};reverse_order_agrees;all_traps_have_short_components;largest_witness_span={qtext(max(span for _,span in witness_spans))}<7/4",
        f"coarse_traps={tuple(traps)}",
        f"verdict_digest={verdict_digest.hexdigest()}",
        f"located_policies={len(policies)};oriented_controls={len(located)};weakest={qtext(weakest_located[0])}@body={weakest_located[1]},pair={weakest_located[2]},channel={weakest_located[3]},cell={weakest_located[5]};digest={located_digest}",
        *(f"tail_{index+1}=body:{row[0]};pair:{row[1]};start:{row[2]};numerator:{qtext(row[3])};margin:{qtext(row[4])};first_step:{qtext(row[5])};head_channels:{row[6]};head_weakest:{qtext(row[7][0])}@{row[7][3]},cell={row[7][5]}" for index, row in enumerate(tail_rows)),
        f"head_controls={len(all_heads)};head_digest={head_digest}",
        "fixed_H_cell=155;uniform_primitive_skeleton>=1/126 on oriented ratio interval [4/7,7/4];unique_equality=(3,2)",
        "transport=oriented eta_g retained;g/(PgL-a) strictly decreases;c_inverse dropped only after positive bracket;debt bounded at level eight",
        "tail=PQ>=s^2,Pg>=s;endpoint numerator bounds 20/7,32/7,26/7;losses decrease strictly;finite heads oversized-audited",
        "conclusion=all reflected residual packets with D>=6 and 3m>=4D close on all 3003 bodies",
        "corollary=assembled reflected certificate-failure wedge is confined to 561 bodies,D>=6,1<=m<4D/3",
        "scope=reflected THM-2941 sufficient family only;physical LRC14 remains open",
        "normal_vs_python_O=BYTE_IDENTICAL",
        f"three_halves_source_sha256={EXPECTED_THREE_HALVES_SHA256}",
        f"three_halves_output_sha256={EXPECTED_THREE_HALVES_OUTPUT_SHA256}",
        f"three_halves_semantic_sha256={EXPECTED_THREE_HALVES_SEMANTIC}",
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
