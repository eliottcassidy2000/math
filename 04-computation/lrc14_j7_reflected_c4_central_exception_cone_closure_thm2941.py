#!/usr/bin/env python3
r"""Sharpen the reflected all-spread cone from ``m>=5D`` to ``m>=4D``.

Let ``E`` be a six-label body, ``L=14*lcm(E)``, and suppose its reflected
levels lie in ``[m,m+D]``.  On each of the 3,001 ordinary bodies, an
uncertified packet has six distinct levels.  On either universal exception,
the fixed slots ``0,1`` are same-level-good and therefore distinct.  Thus the
closest normalized label pair from the preceding cone theorem always has
distinct levels.

If coprime ``P<Q`` satisfy ``P>=4(Q-P)``, the exact phase fibre obeys

    F_PQ,min >= 1/70.

The only equality channel is ``(P,Q)=(4,5)``.  Away from it the stronger
floor ``1/55`` holds.  Indeed the residue correction is at least ``-12/49``;
this proves ``1/55`` once ``PQ>=110``, while the complete smaller bank is

    (4,5),(5,6),(6,7),(7,8),(8,9),(9,10),(9,11).

For every body except ``H=(1,2,3,4,6,12)``, the closest-pair ratio is at most
``1/168``.  At ``D=6,m=24`` the exact 3,002-body invoice

    1/70 - 2D[4(b-a)+b]/[4DL-b] - debt_E(4D)

is strictly positive.  Both subtracted terms decrease with ``D`` and with
``m>=4D``, so this closes all nonhostile bodies uniformly.

For ``H``, use labels ``a=1,b=2``.  Every non-``(4,5)`` reduced channel has
phase floor ``1/55`` and closes by the sharper exact drift bound

    2|qa-pb|/(pL-a) <= 12D/(840D-1).

The channel ``(4,5)`` can occur inside ``[m,m+D]`` only when ``m=4D`` and
the two levels are ``4D,5D``.  The forward orientation closes by the ordinary
cross-determinant bound.  The reverse orientation, label 1 at ``5D`` and
label 2 at ``4D``, is the sole cross-determinant exception.

It has an exact located certificate on the central safe cell ``j=90``.
Writing ``z_1=840D-1,z_2=672D-2``, its two arc families are

    A_k=[78+168k,102+168k]/z_1,       0<=k<5D,
    B_l=[168l,24+168l]/z_2,           0<=l<4D.

The endpoint inequalities say that ``A_k`` and ``B_l`` overlap exactly when
``4k-5l=-2``, hence at ``(k,l)=(5r+2,4r+2)``, ``0<=r<D``.  The scaled length
of this intersection is ``12096D-540-1008r``.  Summing gives

    I_D = 36D(322D-1)/[(840D-1)(672D-2)] > 23/1120 > 1/70.

Meanwhile the total singleton debt is at most
``1/(168D-3)<1/70``.  Bonferroni therefore closes the exceptional packet.
This proves every reflected residual packet with ``D>=6,m>=4D`` closes.
The surviving reflected wedge is ``m<4D``.  The theorem remains a sufficient
THM-2941 certificate, not a proof of physical LRC(14).
"""

from __future__ import annotations

import argparse
import hashlib
from fractions import Fraction as F
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from math import gcd, lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "04-computation/lrc14_j7_reflected_levels_all_q_mass_closure_thm2941.py"
UNIVERSAL = ROOT / "04-computation/lrc14_j7_reflected_universal_pair_chromatic_closure_thm2941.py"
UNIVERSAL_OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_reflected_universal_pair_chromatic_closure_thm2941.out"
C5 = ROOT / "04-computation/lrc14_j7_reflected_all_spread_conical_tail_closure_thm2941.py"
C5_OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_reflected_all_spread_conical_tail_closure_thm2941.out"
OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_reflected_c4_central_exception_cone_closure_thm2941.out"

EXPECTED_BASE_SHA256 = "2cf0866932f775cc493f97093333e81e65ac3aa76a8e439de969aa700c993f31"
EXPECTED_UNIVERSAL_SHA256 = "a6f58c1a52dfc1fca61a239068dbe0b216bac41f1622b98748bc4a6d213fb6e8"
EXPECTED_UNIVERSAL_OUTPUT_SHA256 = "7364d5866171405fa90539a9ad76727c0c52f020ac1a104a1ab4f0276aedd115"
EXPECTED_C5_SHA256 = "d2c8b4a26e87bce98c8ea94e82c38e81eae95504a6a7ef70db53670f197e8ec8"
EXPECTED_C5_OUTPUT_SHA256 = "08553b2e00589e65fe3127d7ec5c9d90257978bd03d6c2f11d953527ff4675a7"
EXPECTED_SEMANTIC_SHA256 = "fd355370e96f889667ea1c5b391eccb1144135ea0f73a24b5a73527536268c97"

BODY_COUNT = 3003
HOSTILE = (1, 2, 3, 4, 6, 12)
MIN_SPREAD = 6
CONE_CONSTANT = 4
PHASE_FLOOR = F(1, 70)
NEXT_PHASE_FLOOR = F(1, 55)
HOSTILE_CELL = 90
DIRECT_CONTROLS = (1, 2, 6, 13)
EXPECTED_NONHOSTILE_WORST = (
    F(1374624312892139, 775692109936614794),
    (1, 3, 4, 6, 8, 12),
    336,
    (1, 2),
    F(1, 168),
    F(15199561634868606, 25209993572939980805),
)
EXPECTED_GENERIC_HOSTILE_MARGIN = F(32290, 11141229)
EXPECTED_FORWARD_HOSTILE_MARGIN = F(247277, 56716170)


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def qtext(value: F) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


for path, expected in (
    (BASE, EXPECTED_BASE_SHA256),
    (UNIVERSAL, EXPECTED_UNIVERSAL_SHA256),
    (UNIVERSAL_OUTPUT, EXPECTED_UNIVERSAL_OUTPUT_SHA256),
    (C5, EXPECTED_C5_SHA256),
    (C5_OUTPUT, EXPECTED_C5_OUTPUT_SHA256),
):
    require(sha256(path) == expected, ("upstream theorem changed", path, sha256(path), expected))


def import_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, ("cannot import", path))
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


R = import_module("c4_cone_base", BASE)
U = import_module("c4_cone_universal", UNIVERSAL)
C5 = import_module("c4_cone_predecessor", C5)


def phase_floor(P: int, Q: int) -> F:
    require(1 <= P < Q and gcd(P, Q) == 1, (P, Q))
    return F(1, 49) + C5.phase_correction(P % 14, Q % 14)[0] / (P * Q)


def intersection_mass(first, second) -> F:
    i = j = 0
    total = F(0)
    while i < len(first) and j < len(second):
        total += max(
            F(0), min(first[i][1], second[j][1]) - max(first[i][0], second[j][0])
        )
        if first[i][1] < second[j][1]:
            i += 1
        else:
            j += 1
    return total


def reverse_formula(D: int) -> F:
    return F(36 * D * (322 * D - 1), (840 * D - 1) * (672 * D - 2))


def reverse_direct_control(D: int):
    L, safe_ranges = R.safe_cell_ranges(HOSTILE)
    require(L == 168 and any(left <= HOSTILE_CELL < right for left, right in safe_ranges),
            (L, safe_ranges, HOSTILE_CELL))
    p, q = 5 * D, 4 * D
    first = R.reflected_level_arcs(L, 1, p, HOSTILE_CELL)
    second = R.reflected_level_arcs(L, 2, q, HOSTILE_CELL)
    require(
        first == R.direct_multiplier_arcs(L, p * L - 1, HOSTILE_CELL)
        and second == R.direct_multiplier_arcs(L, q * L - 2, HOSTILE_CELL),
        ("arc-law mismatch", D),
    )
    z1, z2 = 840 * D - 1, 672 * D - 2
    segments = []
    for left_a, right_a in first:
        for left_b, right_b in second:
            length = min(right_a, right_b) - max(left_a, left_b)
            if length > 0:
                segments.append(length * z1 * z2)
    expected_segments = tuple(12096 * D - 540 - 1008 * r for r in range(D))
    require(tuple(segments) == expected_segments, (D, segments, expected_segments))
    require(
        sum(expected_segments) == 36 * D * (322 * D - 1),
        (D, sum(expected_segments)),
    )
    actual = intersection_mass(first, second)
    require(actual == reverse_formula(D), (D, actual, reverse_formula(D)))
    require(
        actual - F(23, 1120)
        == F(6888 * D - 23, 1120 * (336 * D - 1) * (840 * D - 1)) > 0,
        (D, actual),
    )
    return D, len(first), len(second), len(segments), actual


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    require(tuple((row[0], row[1]) for row in U.EXPECTED_EXCEPTIONS) == C5.EXCEPTIONS,
            "universal exception mismatch")

    small_channels = tuple(
        (P, Q, phase_floor(P, Q))
        for P in range(1, 110)
        for Q in range(P + 1, 110)
        if gcd(P, Q) == 1 and P >= 4 * (Q - P) and P * Q < 110
    )
    require(
        small_channels
        == (
            (4, 5, F(1, 70)),
            (5, 6, F(2, 105)),
            (6, 7, F(1, 49)),
            (7, 8, F(1, 49)),
            (8, 9, F(5, 252)),
            (9, 10, F(2, 105)),
            (9, 11, F(13, 693)),
        ),
        small_channels,
    )
    require(F(1, 49) - F(12, 49 * 40) == PHASE_FLOOR, "C4 phase tail")
    require(F(1, 49) - F(12, 49 * 110) == NEXT_PHASE_FLOOR, "next phase tail")
    require(
        all(value > NEXT_PHASE_FLOOR for P, Q, value in small_channels if (P, Q) != (4, 5)),
        small_channels,
    )

    rows = []
    pair_rows = []
    body_digest = hashlib.sha256()
    D = MIN_SPREAD
    m = CONE_CONSTANT * D
    for body in combinations(range(1, 15), 6):
        L = 14 * lcm(*body)
        ratio, i, j = C5.chosen_pair(body, L)
        pair_rows.append((ratio, body, L, (i, j)))
        body_digest.update(f"{body}|{L}|{ratio}|{(i,j)}\n".encode())
        if body == HOSTILE:
            continue
        debt = C5.singleton_debt(body, L, (m,) * 6)
        delta = body[j] - body[i]
        margin = (
            PHASE_FLOOR
            - F(2 * D * (CONE_CONSTANT * delta + body[j]), m * L - body[j])
            - debt
        )
        rows.append((margin, body, L, (i, j), ratio, debt))

    require(len(pair_rows) == BODY_COUNT and len(rows) == BODY_COUNT - 1,
            (len(pair_rows), len(rows)))
    require(
        max(pair_rows)
        == (F(1, 84), HOSTILE, 168, (0, 1)),
        max(pair_rows),
    )
    require(max(row[0] for row in pair_rows if row[1] != HOSTILE) == F(1, 168),
            "nonhostile ratio boundary")
    worst_nonhostile = min(rows)
    require(worst_nonhostile == EXPECTED_NONHOSTILE_WORST, worst_nonhostile)
    require(worst_nonhostile[0] > 0, worst_nonhostile)

    debt_bound = F(1, 168 * D - 3)
    generic_hostile_margin = NEXT_PHASE_FLOOR - F(12 * D, 840 * D - 1) - debt_bound
    forward_hostile_margin = PHASE_FLOOR - F(6 * D, 672 * D - 1) - debt_bound
    require(generic_hostile_margin == EXPECTED_GENERIC_HOSTILE_MARGIN > 0,
            generic_hostile_margin)
    require(forward_hostile_margin == EXPECTED_FORWARD_HOSTILE_MARGIN > 0,
            forward_hostile_margin)
    require(F(23, 1120) - PHASE_FLOOR == F(1, 160), "reverse overlap gap")
    require(debt_bound < PHASE_FLOOR, debt_bound)
    require(F(6 * D, 840 * D - 1) < 1 and 4 * D - F(6 * D, 840 * D - 1) > 1,
            "homotopy target slope gate")

    # Algebraic monotonicity controls: x/(Ax-b) decreases by exactly b over
    # the positive product of consecutive denominators.  The all-D argument
    # in the docstring applies this identity termwise.
    for x in range(MIN_SPREAD, 40):
        require(
            F(x, 840 * x - 1) - F(x + 1, 840 * (x + 1) - 1)
            == F(1, (840 * x - 1) * (840 * (x + 1) - 1)) > 0,
            ("drift monotonicity", x),
        )
        require(F(1, 168 * (x + 1) - 3) < F(1, 168 * x - 3),
                ("debt monotonicity", x))

    controls = tuple(reverse_direct_control(D) for D in DIRECT_CONTROLS)
    semantic_payload = (
        small_channels,
        tuple(sorted(pair_rows)),
        body_digest.hexdigest(),
        worst_nonhostile,
        generic_hostile_margin,
        forward_hostile_margin,
        controls,
    )
    semantic = hashlib.sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256:
        require(semantic == EXPECTED_SEMANTIC_SHA256, ("semantic digest", semantic))

    lines = [
        "LRC14 reflected C=4 central-exception cone closure exact proof",
        f"universe=bodies:{BODY_COUNT};spread:D>={MIN_SPREAD};cone:m>={CONE_CONSTANT}D",
        f"phase_C4_bank={small_channels};global_floor={qtext(PHASE_FLOOR)};unique_equality_channel=(4,5);off_channel_floor={qtext(NEXT_PHASE_FLOOR)}",
        "pair_selection=ordinary complete bodies use closest normalized label pair;exception slots (0,1) are same-level-good",
        "pair_ratio_dichotomy=hostile E=(1,2,3,4,6,12) unique at 1/84;all other 3002 bodies <=1/168",
        f"nonhostile_D6_invoice_worst={qtext(worst_nonhostile[0])};E={worst_nonhostile[1]};L={worst_nonhostile[2]};pair={worst_nonhostile[3]};ratio={qtext(worst_nonhostile[4])};debt={qtext(worst_nonhostile[5])}",
        f"hostile_off_45_D6_margin={qtext(generic_hostile_margin)};phase=1/55;drift<=12D/(840D-1);debt<=1/(168D-3)",
        f"hostile_forward_45_D6_margin={qtext(forward_hostile_margin)};levels=(4D,5D);phase=1/70;drift=6D/(672D-1)",
        "hostile_reverse_45=sole crossdet exception;levels=(5D,4D);central_cell=90;overlap=36D(322D-1)/[(840D-1)(672D-2)]>23/1120>1/70>debt",
        "reverse_arc_proof=overlaps iff 4k-5l=-2;(k,l)=(5r+2,4r+2);scaled_lengths=12096D-540-1008r",
        "homotopy_gate=|eta|<=6D/(840D-1)<1 and all transported target slopes remain >1",
    ]
    for D, first_count, second_count, segment_count, actual in controls:
        lines.append(
            f"CONTROL;D={D};reverse_levels=({5*D},{4*D});arc_counts=({first_count},{second_count});segments={segment_count};overlap={qtext(actual)}"
        )
    lines.extend((
        "monotonicity=transport and singleton debt decrease with m and D on m>=4D;D=6 invoice is uniform",
        "conclusion=all reflected residual packets with D>=6,m>=4D close on all 3003 bodies",
        "corollary=the still-open reflected D>=6 wedge is confined to m<4D",
        "scope=reflected THM-2941 sufficient family only;physical LRC14 remains open",
        "normal_vs_python_O=BYTE_IDENTICAL",
        f"base_sha256={sha256(BASE)}",
        f"universal_sha256={sha256(UNIVERSAL)}",
        f"universal_output_sha256={sha256(UNIVERSAL_OUTPUT)}",
        f"c5_source_sha256={sha256(ROOT / '04-computation/lrc14_j7_reflected_all_spread_conical_tail_closure_thm2941.py')}",
        f"c5_output_sha256={sha256(C5_OUTPUT)}",
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
