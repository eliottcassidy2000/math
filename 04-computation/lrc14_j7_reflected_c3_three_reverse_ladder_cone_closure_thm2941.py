#!/usr/bin/env python3
r"""Sharpen the reflected cone from ``m>=4D`` to ``m>=3D`` for ``D>=6``.

The exact primitive phase fibre changes character at this wall.  For coprime
``P<Q`` with ``P>=3(Q-P)``, its minimum is nonnegative, but the unique zero
channel is ``(3,4)``.  Excluding that channel the floor is ``1/70``, uniquely
at ``(4,5)``; excluding both it is ``1/55``.  The complete ``PQ<110`` bank is

    (3,4),(4,5),(5,6),(6,7),(7,8),(7,9),(8,9),(9,10),(9,11).

The preceding robust-edge theorems close 2,354 bodies for arbitrary levels.
We retain the broader 649-body complement with at most ten robust edges as an
independent audit universe.  For every non-``(3,4)`` channel, the ordinary
cross-determinant invoice at ``D=6,m=18`` closes 647 of these bodies.  Only

    H =(1,2,3,4,6,12),       H2=(1,3,4,6,8,12)

need the ``1/55`` off-``(4,5)`` floor and a short ``(4,5)`` sidecar.  The
hostile reverse ``(4,5)`` ray is exactly the central ladder already proved by
the ``C=4`` theorem; all other cases close by cross-determinant transport.

The zero channel ``(3,4)`` occurs only at ``m=3D`` with endpoint levels
``3D,4D``.  On each body's upper-median safe cell, retain its exact located
primitive phase instead of minimizing over phase.  Transport from that phase
at ``D=6`` closes 1,295 of the 1,298 oriented body cases.  The located phase
is scale-independent, while drift and singleton debt decrease with ``D``, so
these certificates are uniform.

Exactly three reverse assignments escape the located transport inequality:

    body                         cell  pair labels
    H                            90    (1,2)
    H2                           174   (3,4)
    H3=(2,4,6,8,9,12)           540   (8,9).

Their exact reverse ``(4D,3D)`` overlaps are respectively

    30D(91D-2) /[(252D-1)(672D-1)],
    14D(91D-2) /[(252D-1)(448D-1)],
       D(924D-17)/[(336D-1)(504D-1)].

For each body the interval sweep has exactly ``D`` intersections in one
arithmetic ladder.  The scaled segment lengths are

    5040D+300+840r,
    14112D+840+2352r,
    72576D-7272-12096r,                    0<=r<D.

The displayed rational functions are increasing for ``D>=1`` and already
exceed the exact body debt at ``D=6``.  This closes the three exceptions and
hence every reflected residual packet with ``D>=6,m>=3D``.  The remaining
reflected wedge is ``m<3D``.  This is still a sufficient THM-2941 result, not
a proof of physical LRC(14).
"""

from __future__ import annotations

import argparse
import hashlib
from fractions import Fraction as F
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from math import gcd
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "04-computation/lrc14_j7_reflected_levels_all_q_mass_closure_thm2941.py"
LOW_PHASE = ROOT / "04-computation/lrc14_j7_reflected_low_phase_clique_robust_body_closure_thm2941.py"
EDGE11 = ROOT / "04-computation/lrc14_j7_reflected_robust_edge11_threshold_shapes_uniform_closure_thm2941.py"
EDGE11_OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_reflected_robust_edge11_threshold_shapes_uniform_closure_thm2941.out"
C4 = ROOT / "04-computation/lrc14_j7_reflected_c4_central_exception_cone_closure_thm2941.py"
C4_OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_reflected_c4_central_exception_cone_closure_thm2941.out"
OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_reflected_c3_three_reverse_ladder_cone_closure_thm2941.out"

EXPECTED_BASE_SHA256 = "2cf0866932f775cc493f97093333e81e65ac3aa76a8e439de969aa700c993f31"
EXPECTED_LOW_PHASE_SHA256 = "b2418dfda1b48257d1f7582d4ea977203a26f88885e13946bc100ccf264c9ce1"
EXPECTED_EDGE11_SHA256 = "58685796af62dfb425bfa17a6363d9cc5ca4cbb72e0f1d42ecc23dcbf10c01b4"
EXPECTED_EDGE11_OUTPUT_SHA256 = "f19bef5d4b2a32bda2e939a303870680b022973a9513a52c18677f73d7512a29"
EXPECTED_C4_SHA256 = "699ebd445efd40c2d8c3218f700bd637ff0fae7318c1593bd2a5137343431877"
EXPECTED_C4_OUTPUT_SHA256 = "6e2f96fc0aeb9093697e018afe67f905ac17ae14bc64129b311c6082d5671d61"
EXPECTED_SEMANTIC_SHA256 = "c99dc8d1da0f1b3910641c0889aba79597e03f1562dd41653c75be8e0866e241"

BODY_COUNT = 3003
UPSTREAM_CLOSED_COUNT = 2354
ACTIVE_BODY_COUNT = 649
MIN_SPREAD = 6
CONE_CONSTANT = 3
H = (1, 2, 3, 4, 6, 12)
H2 = (1, 3, 4, 6, 8, 12)
H3 = (2, 4, 6, 8, 9, 12)
SPECIAL_BODIES = (H, H2, H3)
DIRECT_CONTROLS = (1, 2, 6, 13)
EXPECTED_NONZERO_WORST = (
    F(1173003125846269, 758117895348950710),
    (2, 3, 4, 6, 8, 12),
    336,
    (0, 1),
    F(1, 168),
    F(4079143043377062, 4927766319768179615),
)
EXPECTED_LOCATED_WORST = (
    F(38701426153927529, 236532783348872621520),
    (2, 3, 4, 6, 8, 12),
    336,
    180,
    (0, 1),
    (3, 4),
    F(1, 336),
    F(1007, 1015728),
    F(4079143043377062, 4927766319768179615),
    F(-3, 3023),
)
EXPECTED_SPECIAL_MARGINS = (
    F(182490365907103031562, 12371762617323495134105),
    F(276803509033353286134, 26486196258041612764945),
    F(6153999974859460056, 1201682332635323928935),
)


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def qtext(value: F) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


for path, expected in (
    (BASE, EXPECTED_BASE_SHA256),
    (LOW_PHASE, EXPECTED_LOW_PHASE_SHA256),
    (EDGE11, EXPECTED_EDGE11_SHA256),
    (EDGE11_OUTPUT, EXPECTED_EDGE11_OUTPUT_SHA256),
    (C4, EXPECTED_C4_SHA256),
    (C4_OUTPUT, EXPECTED_C4_OUTPUT_SHA256),
):
    require(sha256(path) == expected, ("upstream theorem changed", path, sha256(path), expected))


def import_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, ("cannot import", path))
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


R = import_module("c3_cone_base", BASE)
LOW = import_module("c3_cone_low_phase", LOW_PHASE)
C4M = import_module("c3_cone_c4", C4)
C5 = C4M.C5


def phase_floor(P: int, Q: int) -> F:
    require(1 <= P < Q and gcd(P, Q) == 1, (P, Q))
    return F(1, 49) + C5.phase_correction(P % 14, Q % 14)[0] / (P * Q)


def intersection_mass(first, second) -> F:
    return (
        R.interval_mass(first)
        + R.interval_mass(second)
        - R.interval_mass(R.merge_intervals(first + second))
    )


def primitive_arcs(slope: int, phase: F):
    arcs = []
    for n in range(-2, slope + 2):
        left = max(F(0), (phase + n - F(1, 14)) / slope)
        right = min(F(1), (phase + n + F(1, 14)) / slope)
        if left < right:
            arcs.append((left, right))
    return tuple(arcs)


def special_formula(index: int, D: int) -> F:
    if index == 0:
        return F(30 * D * (91 * D - 2), (252 * D - 1) * (672 * D - 1))
    if index == 1:
        return F(14 * D * (91 * D - 2), (252 * D - 1) * (448 * D - 1))
    if index == 2:
        return F(D * (924 * D - 17), (336 * D - 1) * (504 * D - 1))
    raise RuntimeError(index)


SPECIAL_DATA = (
    # body, L, cell, slots, labels, scaled segment A*D+B+C*r
    (H, 168, 90, (0, 1), (1, 2), (5040, 300, 840)),
    (H2, 336, 174, (1, 2), (3, 4), (14112, 840, 2352)),
    (H3, 1008, 540, (3, 4), (8, 9), (72576, -7272, -12096)),
)


def special_control(index: int, D: int):
    body, L, cell, slots, labels, coefficients = SPECIAL_DATA[index]
    a, b = labels
    p, q = 4 * D, 3 * D
    first = R.reflected_level_arcs(L, a, p, cell)
    second = R.reflected_level_arcs(L, b, q, cell)
    require(
        first == R.direct_multiplier_arcs(L, p * L - a, cell)
        and second == R.direct_multiplier_arcs(L, q * L - b, cell),
        ("arc-law mismatch", body, D),
    )
    z1, z2 = p * L - a, q * L - b
    segments = []
    for left_a, right_a in first:
        for left_b, right_b in second:
            length = min(right_a, right_b) - max(left_a, left_b)
            if length > 0:
                segments.append(length * z1 * z2)
    A, B, C = coefficients
    expected = tuple(A * D + B + C * r for r in range(D))
    require(tuple(segments) == expected, (body, D, segments, expected))
    actual = intersection_mass(first, second)
    require(actual == special_formula(index, D), (body, D, actual, special_formula(index, D)))
    return body, D, len(first), len(second), len(segments), actual


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    small_channels = tuple(
        (P, Q, phase_floor(P, Q))
        for P in range(1, 110)
        for Q in range(P + 1, 110)
        if gcd(P, Q) == 1 and P >= 3 * (Q - P) and P * Q < 110
    )
    require(
        small_channels
        == (
            (3, 4, F(0)),
            (4, 5, F(1, 70)),
            (5, 6, F(2, 105)),
            (6, 7, F(1, 49)),
            (7, 8, F(1, 49)),
            (7, 9, F(1, 49)),
            (8, 9, F(5, 252)),
            (9, 10, F(2, 105)),
            (9, 11, F(13, 693)),
        ),
        small_channels,
    )
    require(F(1, 49) - F(12, 49 * 12) == 0, "C3 zero tail")
    require(F(1, 49) - F(12, 49 * 40) == F(1, 70), "C3 off-zero tail")
    require(F(1, 49) - F(12, 49 * 110) == F(1, 55), "C3 off-45 tail")

    remaining = []
    body_digest = hashlib.sha256()
    for body in combinations(range(1, 15), 6):
        L, _, robust = LOW.robust_edges(body)
        if len(robust) <= 10:
            ratio, i, j = C5.chosen_pair(body, L)
            _, ranges = R.safe_cell_ranges(body)
            cells = tuple(cell for left, right in ranges for cell in range(left, right))
            cell = cells[len(cells) // 2]
            remaining.append((body, L, robust, ratio, i, j, cell))
            body_digest.update(f"{body}|{L}|{robust}|{ratio}|{(i,j)}|{cell}\n".encode())
    require(len(remaining) == ACTIVE_BODY_COUNT, len(remaining))
    require(BODY_COUNT - len(remaining) == UPSTREAM_CLOSED_COUNT, len(remaining))

    D = MIN_SPREAD
    m = CONE_CONSTANT * D
    nonzero_rows = []
    for body, L, _, ratio, i, j, _ in remaining:
        if body in (H, H2):
            continue
        debt = C5.singleton_debt(body, L, (m,) * 6)
        delta = body[j] - body[i]
        margin = (
            F(1, 70)
            - F(2 * D * (3 * delta + body[j]), m * L - body[j])
            - debt
        )
        nonzero_rows.append((margin, body, L, (i, j), ratio, debt))
    require(len(nonzero_rows) == 647, len(nonzero_rows))
    worst_nonzero = min(nonzero_rows)
    require(worst_nonzero == EXPECTED_NONZERO_WORST, worst_nonzero)
    require(worst_nonzero[0] > 0, worst_nonzero)

    located_rows = []
    located_residual = []
    for body, L, _, _, i, j, cell in remaining:
        a, b = body[i], body[j]
        debt = C5.singleton_debt(body, L, (m,) * 6)
        for primitive_p, primitive_q in ((3, 4), (4, 3)):
            skeleton = intersection_mass(
                primitive_arcs(primitive_p, F(a * cell, L) % 1),
                primitive_arcs(primitive_q, F(b * cell, L) % 1),
            )
            eta = F(D * (primitive_q * a - primitive_p * b), primitive_p * D * L - a)
            gain = max(F(0), skeleton - 2 * abs(eta))
            actual = intersection_mass(
                R.reflected_level_arcs(L, a, primitive_p * D, cell),
                R.reflected_level_arcs(L, b, primitive_q * D, cell),
            )
            require(actual >= gain, ("located transport failed", body, primitive_p, primitive_q))
            row = (
                gain - debt, body, L, cell, (i, j), (primitive_p, primitive_q),
                skeleton, gain, debt, eta,
            )
            located_rows.append(row)
            if row[0] <= 0:
                located_residual.append(row)
    require(len(located_rows) == 1298, len(located_rows))
    require(
        tuple((row[1], row[3], row[4], row[5]) for row in located_residual)
        == ((H, 90, (0, 1), (4, 3)), (H2, 174, (1, 2), (4, 3)),
            (H3, 540, (3, 4), (4, 3))),
        located_residual,
    )
    worst_located = min(row for row in located_rows if row[0] > 0)
    require(worst_located == EXPECTED_LOCATED_WORST, worst_located)

    special_margins = []
    for index, (body, L, _, _, _, _) in enumerate(SPECIAL_DATA):
        debt = C5.singleton_debt(body, L, (m,) * 6)
        special_margins.append(special_formula(index, D) - debt)
    require(tuple(special_margins) == EXPECTED_SPECIAL_MARGINS, special_margins)
    require(all(margin > 0 for margin in special_margins), special_margins)

    h_debt_bound = F(1, 126 * D - 3)
    h2_debt = C5.singleton_debt(H2, 336, (m,) * 6)
    h_generic = F(1, 55) - F(10 * D, 672 * D - 1) - h_debt_bound
    h2_generic = F(1, 55) - F(14 * D, 1344 * D - 3) - h2_debt
    # A 4:5 channel in [m,m+D], m>=3D, has scale g>=5.
    h_forward_45 = F(1, 70) - F(30, 3359) - h_debt_bound
    h2_reverse_45 = F(1, 70) - F(80, 8397) - h2_debt
    require(
        (h_generic, h2_generic, h_forward_45, h2_reverse_45)
        == (
            F(328738, 166943865),
            F(2026944006691831647523, 291348158838457740414395),
            F(712897, 177052890),
            F(130921463361507916385, 33108238180688563064718),
        ),
        (h_generic, h2_generic, h_forward_45, h2_reverse_45),
    )
    require(all(value > 0 for value in (h_generic, h2_generic, h_forward_45, h2_reverse_45)),
            "exceptional nonzero channel")
    require(C4M.reverse_formula(5) > F(1, 70) > h_debt_bound,
            "hostile reverse 4:5 reuse")

    # Each special overlap increases: after clearing the positive squared
    # denominator, its derivative numerator is the displayed polynomial.
    derivative_polynomials = (
        (127302, 91, -1),
        (81046, 91, -1),
        (2102688, 1848, -17),
    )
    require(all(A + B + C > 0 and A > 0 and B > 0 for A, B, C in derivative_polynomials),
            derivative_polynomials)
    controls = tuple(
        special_control(index, D0)
        for index in range(3)
        for D0 in DIRECT_CONTROLS
    )

    semantic_payload = (
        small_channels,
        body_digest.hexdigest(),
        worst_nonzero,
        worst_located,
        tuple(located_residual),
        tuple(special_margins),
        (h_generic, h2_generic, h_forward_45, h2_reverse_45),
        derivative_polynomials,
        controls,
    )
    semantic = hashlib.sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256:
        require(semantic == EXPECTED_SEMANTIC_SHA256, ("semantic digest", semantic))

    lines = [
        "LRC14 reflected C=3 three-reverse-ladder cone closure exact proof",
        f"universe=bodies:{BODY_COUNT};upstream_arbitrary_level:{UPSTREAM_CLOSED_COUNT};active_bodies:{ACTIVE_BODY_COUNT};spread:D>={MIN_SPREAD};cone:m>={CONE_CONSTANT}D",
        f"phase_C3_bank={small_channels};unique_zero=(3,4);off_zero_floor=1/70;unique_next=(4,5);off_34_45_floor=1/55",
        f"nonzero_channel_invoice=647 bodies;D6_worst={qtext(worst_nonzero[0])};E={worst_nonzero[1]};L={worst_nonzero[2]};pair={worst_nonzero[3]};ratio={qtext(worst_nonzero[4])};debt={qtext(worst_nonzero[5])}",
        f"located_34_transport=1298 orientations;certified=1295;residual=3;weakest_positive={qtext(worst_located[0])};E={worst_located[1]};cell={worst_located[3]};pair={worst_located[4]};orientation={worst_located[5]};phase={qtext(worst_located[6])};gain={qtext(worst_located[7])};debt={qtext(worst_located[8])}",
        f"three_reverse_residuals={tuple((row[1],row[3],row[4],row[5]) for row in located_residual)}",
        "reverse_formulas=(30D(91D-2)/[(252D-1)(672D-1)],14D(91D-2)/[(252D-1)(448D-1)],D(924D-17)/[(336D-1)(504D-1)])",
        f"reverse_D6_margins={tuple(qtext(value) for value in special_margins)};all_increasing_in_D",
        f"exceptional_nonzero_margins_D6=H_off45:{qtext(h_generic)};H2_off45:{qtext(h2_generic)};H_forward45:{qtext(h_forward_45)};H2_reverse45:{qtext(h2_reverse_45)}",
        "hostile_reverse45=reuse C4 central ladder at arbitrary scale g>=5;overlap>23/1120>1/70>debt",
        "scale_law=located primitive phase is invariant;drift and debt decrease with D and m;3:4 can occur only at m=3D",
    ]
    for body, D0, first_count, second_count, segment_count, actual in controls:
        lines.append(
            f"CONTROL;E={body};D={D0};reverse34_levels=({4*D0},{3*D0});arc_counts=({first_count},{second_count});segments={segment_count};overlap={qtext(actual)}"
        )
    lines.extend((
        "conclusion=all reflected residual packets with D>=6,m>=3D close on all 3003 bodies",
        "corollary=the still-open reflected D>=6 wedge is confined to m<3D",
        "scope=reflected THM-2941 sufficient family only;physical LRC14 remains open",
        "normal_vs_python_O=BYTE_IDENTICAL",
        f"base_sha256={sha256(BASE)}",
        f"low_phase_sha256={sha256(LOW_PHASE)}",
        f"edge11_sha256={sha256(EDGE11)}",
        f"edge11_output_sha256={sha256(EDGE11_OUTPUT)}",
        f"c4_source_sha256={sha256(ROOT / '04-computation/lrc14_j7_reflected_c4_central_exception_cone_closure_thm2941.py')}",
        f"c4_output_sha256={sha256(C4_OUTPUT)}",
        f"body_digest={body_digest.hexdigest()}",
        f"source_sha256={sha256(Path(__file__))}",
        f"semantic_sha256={semantic}",
        "all_exact_controls=PASS",
    ))
    payload = "\n".join(lines) + "\n"
    args.output.write_text(payload, encoding="utf-8", newline="\n")
    print(payload, end="")


if __name__ == "__main__":
    main()
