#!/usr/bin/env python3
r"""Uniform conical tail for every reflected level spread.

Let ``E`` be a six-label body, ``L=14*lcm(E)``, and put

    z_e=q_e L-e,             m <= q_e <= m+D.

The cross-determinant lemma from the reflected ``D=5`` theorem gives, for
distinct levels ``p,q`` on labels ``a<b``, a body-safe-cell pair floor

    c^-1 (F_PQ,min-2|eta|),
    c=1-a/(pL),              eta=(qa-pb)/(pL-a).          (1)

This file proves that every such packet closes whenever ``D>=4`` and
``m>=5D``.  The point is that the estimate is conical: it depends on the
ratio ``D/m``, not on fixing ``D`` first.  Spreads below four are already
closed by the stronger bounded-spread theorem.

There is always an admissible distinct pair with

    2(b-a)/L <= 1/84.                                  (2)

For the 3,001 bodies whose universal same-level graph is complete, any
uncertified packet has six distinct levels, so we take the body pair
minimizing the left side.  On either exceptional body, slots ``0,1`` form a
same-level-good edge and hence have distinct levels in every residual word;
their ratio is ``1/126126``.  An exact all-body audit shows that equality in
(2) occurs only on ``E=(1,2,3,4,6,12)``.

The phase fibre has a sharper exact functional form than its earlier
``1/(2PQ)`` error bound.  The quantity ``T_s(z)-s^2`` is one-periodic in
``s``.  An exact ``14 by 14`` residue bank gives

    min_z [(T_A-A^2)-(T_B-B^2)] >= -12/49,

and hence

    F_PQ,min >= 1/49-12/(49PQ).                        (3)

If ``m>=5D`` and ``r=Q-P``, then ``P>=5r``.  Formula (3) gives ``1/55``
once ``PQ>=110``.  The five smaller coprime pairs are

    (5,6),(6,7),(7,8),(8,9),(9,10),

and their exact floors are all larger.  Equality is attained at ``(10,11)``,
so the near-diagonal floor is sharply

    F_PQ,min >= 1/55.                                  (3a)

Also

    |qa-pb| <= m(b-a)+Db,

and subtracting the limiting transport term exactly gives

    2[m(b-a)+Db]/[mL-b] - 2(b-a)/L
      =2b[LD+(b-a)]/[L(mL-b)]
      <=28(D+13/168)/(168m-14).                        (4)

The singleton debt decreases with ``m``.  Its exact ``m=1`` all-body maximum
is

    1915198706/76797355635 < 1/39

on the same hostile body.  Now fix ``C=5`` and ``D_0=4``.  At ``m=CD`` the
transport bound is

    2D[C(b-a)+b]/[CDL-b],

which decreases with ``D``; the debt also decreases.  Consequently all
``D>=D_0`` are bounded below by the finite body expression

    M_C(E)=1/55
           -2D_0[C(b-a)+b]/[CD_0L-b]
           -sum_e e/[7(CD_0L-e)].                      (5)

An exact audit of all 3,003 bodies gives the unique minimum

    M_5(1,2,3,4,6,12)
      =33073441051905679/105976173248680729995 > 0.     (6)

At ``C=4`` the sharp near-diagonal phase floor is only ``1/70`` and the
hostile limiting phase-minus-transport margin is ``-1/280`` before debt.
Thus five is sharp for this proof invoice (not claimed sharp for the
underlying reflected problem).

The homotopy used in (1) is safely inside its slope range: the same bounds
give ``|eta|<1``, while both integer slopes are at least ``m>=20``.

Hence every reflected residual packet with ``D>=4,m>=5D`` closes on every
body.  In particular all still-open spreads ``D>=6`` are reduced to the
wedge ``m<5D``.  This is a sufficient theorem inside the
THM-2941 reflected family, not a proof of physical LRC(14).
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


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "04-computation/lrc14_j7_reflected_levels_all_q_mass_closure_thm2941.py"
UNIVERSAL = ROOT / "04-computation/lrc14_j7_reflected_universal_pair_chromatic_closure_thm2941.py"
UNIVERSAL_OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_reflected_universal_pair_chromatic_closure_thm2941.out"
D5 = ROOT / "04-computation/lrc14_j7_reflected_d5_crossdet_tail_closure_thm2941.py"
D5_OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_reflected_d5_crossdet_tail_closure_thm2941.out"
OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_reflected_all_spread_conical_tail_closure_thm2941.out"

EXPECTED_BASE_SHA256 = "2cf0866932f775cc493f97093333e81e65ac3aa76a8e439de969aa700c993f31"
EXPECTED_UNIVERSAL_SHA256 = "a6f58c1a52dfc1fca61a239068dbe0b216bac41f1622b98748bc4a6d213fb6e8"
EXPECTED_UNIVERSAL_OUTPUT_SHA256 = "7364d5866171405fa90539a9ad76727c0c52f020ac1a104a1ab4f0276aedd115"
EXPECTED_D5_SHA256 = "1575b9fabec292bccf0bd639b47b3775922a1531e421b8e4441c6909cc2cedb7"
EXPECTED_D5_OUTPUT_SHA256 = "73e8242e68431bd42d3a39b40a8dd0cc6ae8f9e1439e80038a54835bedc3ac55"
EXPECTED_SEMANTIC_SHA256 = "21471954712a0f3150304723b1e9f91de15a08c1221797cb0033542e9333a0fc"

BODY_COUNT = 3003
COMPLETE_BODY_COUNT = 3001
HOSTILE = (1, 2, 3, 4, 6, 12)
EXCEPTIONS = (
    ((1, 2, 7, 9, 11, 13), ((2, 3), (3, 4), (4, 5))),
    ((2, 4, 7, 9, 11, 13), ((2, 4), (3, 5))),
)

CONE_CONSTANT = 5
MIN_CONE_SPREAD = 4
NEAR_DIAGONAL_FLOOR = F(1, 55)
EXPECTED_PHASE_CORRECTION_MINIMUM = -F(12, 49)
EXPECTED_WORST_CONE_MARGIN = F(33073441051905679, 105976173248680729995)
EXPECTED_C4_LIMIT_MARGIN = -F(1, 280)
EXPECTED_DEBT_MAXIMUM = F(1915198706, 76797355635)
EXPECTED_DEBT_GAP = F(53964259, 76797355635)
DIRECT_CASES = (5, 13)


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def qtext(value: F) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def triangle_sum(s: F, z: F) -> F:
    """Exact periodized tent."""
    bound = s.numerator // s.denominator + 3
    return sum(
        (max(F(0), s - abs(z + n)) for n in range(-bound, bound + 1)),
        F(0),
    )


def phase_correction(residue_p: int, residue_q: int) -> tuple[F, F]:
    """Minimum numerator correction for one pair of residues modulo fourteen."""
    A = F((residue_p + residue_q) % 14, 14)
    B = F((residue_q - residue_p) % 14, 14)
    events = {F(0), F(1)}
    for s in (A, B):
        for n in range(-3, 4):
            for z in (-F(n), s - n, -s - n):
                if 0 <= z <= 1:
                    events.add(z)
    return min(
        (
            (triangle_sum(A, z) - A * A)
            - (triangle_sum(B, z) - B * B),
            z,
        )
        for z in events
    )


for path, expected in (
    (BASE, EXPECTED_BASE_SHA256),
    (UNIVERSAL, EXPECTED_UNIVERSAL_SHA256),
    (UNIVERSAL_OUTPUT, EXPECTED_UNIVERSAL_OUTPUT_SHA256),
    (D5, EXPECTED_D5_SHA256),
    (D5_OUTPUT, EXPECTED_D5_OUTPUT_SHA256),
):
    require(sha256(path) == expected, ("upstream theorem changed", path, sha256(path), expected))


def import_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, ("cannot import", path))
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


R = import_module("conical_tail_base", BASE)
U = import_module("conical_tail_universal", UNIVERSAL)


def singleton_debt(body: tuple[int, ...], ruler: int, levels: tuple[int, ...]) -> F:
    return sum(
        (F(e, 7 * (q * ruler - e)) for e, q in zip(body, levels)),
        F(0),
    )


def intersection_mass(first, second) -> F:
    i = 0
    j = 0
    total = F(0)
    while i < len(first) and j < len(second):
        total += max(
            F(0),
            min(first[i][1], second[j][1]) - max(first[i][0], second[j][0]),
        )
        if first[i][1] < second[j][1]:
            i += 1
        else:
            j += 1
    return total


def chosen_pair(body: tuple[int, ...], ruler: int) -> tuple[F, int, int]:
    if body in {row[0] for row in EXCEPTIONS}:
        i, j = 0, 1
        return F(2 * (body[j] - body[i]), ruler), i, j
    return min(
        (F(2 * (body[j] - body[i]), ruler), i, j)
        for i, j in combinations(range(6), 2)
    )


def direct_control(D: int):
    body = HOSTILE
    ruler, safe_ranges = R.safe_cell_ranges(body)
    m = CONE_CONSTANT * D
    levels = (m + D, m, m + 1, m + 2, m + 3, m + 4)
    a, b = body[0], body[1]
    p, q = levels[0], levels[1]
    divisor = gcd(p, q)
    P, Q = sorted((p // divisor, q // divisor))
    phase_floor = F(1, 49) + phase_correction(P % 14, Q % 14)[0] / (P * Q)
    c = F(p * ruler - a, p * ruler)
    eta = F(q * a - p * b, p * ruler - a)
    transported = (phase_floor - 2 * abs(eta)) / c
    debt = singleton_debt(body, ruler, levels)
    actual = min(
        (
            intersection_mass(
                R.reflected_level_arcs(ruler, a, p, cell),
                R.reflected_level_arcs(ruler, b, q, cell),
            ),
            cell,
        )
        for left, right in safe_ranges
        for cell in range(left, right)
    )
    require(P == CONE_CONSTANT and Q == CONE_CONSTANT + 1, (D, P, Q))
    require(
        transported - debt > EXPECTED_WORST_CONE_MARGIN,
        (D, transported, debt, EXPECTED_WORST_CONE_MARGIN),
    )
    require(actual[0] >= transported > debt, (D, actual, transported, debt))
    return D, m, levels, phase_floor, c, eta, transported, debt, actual


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    universal_exceptions = tuple((row[0], row[1]) for row in U.EXPECTED_EXCEPTIONS)
    require(universal_exceptions == EXCEPTIONS, (universal_exceptions, EXCEPTIONS))

    # The tent error is integer-periodic in its radius.  The linked pair of
    # residues therefore reduces to this complete 14 by 14 bank.
    for residue in range(14):
        alpha = F(residue, 14)
        for integer in range(3):
            for z in (F(0), F(1, 14), F(3, 7), F(1, 2), F(13, 14)):
                require(
                    triangle_sum(alpha + integer, z) - (alpha + integer) ** 2
                    == triangle_sum(alpha, z) - alpha ** 2,
                    ("tent radius recurrence", residue, integer, z),
                )
    correction_rows = tuple(
        (*phase_correction(residue_p, residue_q), residue_p, residue_q)
        for residue_p in range(14)
        for residue_q in range(14)
    )
    require(len(correction_rows) == 196, len(correction_rows))
    require(
        min(row[0] for row in correction_rows)
        == EXPECTED_PHASE_CORRECTION_MINIMUM,
        min(correction_rows),
    )
    require(
        sum(row[0] == EXPECTED_PHASE_CORRECTION_MINIMUM for row in correction_rows) == 8,
        "sharp phase-correction multiplicity",
    )
    near_small = tuple(
        (P, Q)
        for P in range(1, 110)
        for Q in range(P + 1, 110)
        if gcd(P, Q) == 1 and P >= 5 * (Q - P) and P * Q < 110
    )
    require(
        near_small == ((5, 6), (6, 7), (7, 8), (8, 9), (9, 10)),
        near_small,
    )
    near_small_floors = tuple(
        (
            F(1, 49) + phase_correction(P % 14, Q % 14)[0] / (P * Q),
            P,
            Q,
        )
        for P, Q in near_small
    )
    require(
        near_small_floors
        == (
            (F(2, 105), 5, 6),
            (F(1, 49), 6, 7),
            (F(1, 49), 7, 8),
            (F(5, 252), 8, 9),
            (F(2, 105), 9, 10),
        ),
        near_small_floors,
    )
    require(
        F(1, 49)
        + phase_correction(10, 11)[0] / 110
        == NEAR_DIAGONAL_FLOOR,
        "sharp near-diagonal equality",
    )

    body_rows = []
    pair_histogram: Counter[F] = Counter()
    body_digest = hashlib.sha256()
    debt_rows = []
    cone_rows = []
    for body in combinations(range(1, 15), 6):
        ruler = 14 * lcm(*body)
        ratio, i, j = chosen_pair(body, ruler)
        require(ruler >= 168, (body, ruler))
        require(ratio <= F(1, 84), (body, ruler, ratio, i, j))
        levels = (1,) * 6
        debt = singleton_debt(body, ruler, levels)
        row = (ratio, body, ruler, (i, j), debt)
        body_rows.append(row)
        debt_rows.append((debt, body, ruler))
        delta = body[j] - body[i]
        cone_level = CONE_CONSTANT * MIN_CONE_SPREAD
        cone_debt = singleton_debt(body, ruler, (cone_level,) * 6)
        cone_margin = (
            NEAR_DIAGONAL_FLOOR
            - F(
                2
                * MIN_CONE_SPREAD
                * (CONE_CONSTANT * delta + body[j]),
                cone_level * ruler - body[j],
            )
            - cone_debt
        )
        cone_rows.append((cone_margin, body, ruler, (i, j), cone_debt))
        pair_histogram[ratio] += 1
        body_digest.update(f"{row}\n".encode())

    require(len(body_rows) == BODY_COUNT, len(body_rows))
    require(BODY_COUNT - len(EXCEPTIONS) == COMPLETE_BODY_COUNT, COMPLETE_BODY_COUNT)
    worst_pair_rows = tuple(row for row in body_rows if row[0] == max(r[0] for r in body_rows))
    require(
        worst_pair_rows == ((F(1, 84), HOSTILE, 168, (0, 1), EXPECTED_DEBT_MAXIMUM),),
        worst_pair_rows,
    )
    worst_debt_rows = tuple(row for row in debt_rows if row[0] == max(r[0] for r in debt_rows))
    require(
        worst_debt_rows == ((EXPECTED_DEBT_MAXIMUM, HOSTILE, 168),),
        worst_debt_rows,
    )
    require(F(1, 39) - EXPECTED_DEBT_MAXIMUM == EXPECTED_DEBT_GAP > 0, EXPECTED_DEBT_GAP)

    # Exact finite invoice behind the universal D>=4, m>=5D cone.
    worst_cone = min(cone_rows)
    require(
        worst_cone
        == (
            EXPECTED_WORST_CONE_MARGIN,
            HOSTILE,
            168,
            (0, 1),
            F(2298991660251998, 1926839513612376909),
        ),
        worst_cone,
    )
    require(
        F(1, 70) - F(2 * (4 + 2), 4 * 168) == EXPECTED_C4_LIMIT_MARGIN < 0,
        EXPECTED_C4_LIMIT_MARGIN,
    )
    require(
        F(13) + F(14, CONE_CONSTANT)
        < F(168) - F(14, CONE_CONSTANT),
        "homotopy slope gate",
    )

    controls = tuple(direct_control(D) for D in DIRECT_CASES)
    semantic_payload = (
        len(body_rows),
        worst_pair_rows,
        worst_debt_rows,
        correction_rows,
        near_small_floors,
        tuple(sorted(pair_histogram.items())),
        body_digest.hexdigest(),
        worst_cone,
        EXPECTED_C4_LIMIT_MARGIN,
        controls,
    )
    semantic = hashlib.sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic == EXPECTED_SEMANTIC_SHA256, (semantic, EXPECTED_SEMANTIC_SHA256))

    lines = [
        "LRC14 reflected all-spread conical tail closure exact proof",
        f"universe=bodies:{len(body_rows)};complete_same_level_graph:{COMPLETE_BODY_COUNT};exceptions:{EXCEPTIONS}",
        f"cone=every spread D>={MIN_CONE_SPREAD} closes whenever minimum level m>={CONE_CONSTANT}D",
        "pair_selection=complete bodies use minimum 2(b-a)/L;exceptions use same-level-good slots (0,1)",
        f"sharp_pair_ratio=max_E min_pair 2(b-a)/L={qtext(worst_pair_rows[0][0])};unique_body={HOSTILE};L=168",
        f"singleton_debt=m^-1 monotone envelope;maximum_at_m1={qtext(EXPECTED_DEBT_MAXIMUM)};gap_below_1/39={qtext(EXPECTED_DEBT_GAP)}",
        f"phase_correction_bank=196;minimum={qtext(EXPECTED_PHASE_CORRECTION_MINIMUM)};multiplicity=8",
        f"sharp_near_diagonal_floor=P>=5(Q-P);minimum={qtext(NEAR_DIAGONAL_FLOOR)};equality=(10,11)",
        f"finite_cone_invoice=M_C(E)=1/55-2D0[C(b-a)+b]/[CD0L-b]-debt_E(CD0);C={CONE_CONSTANT};D0={MIN_CONE_SPREAD}",
        f"unique_worst_cone_margin={qtext(worst_cone[0])};body={worst_cone[1]};pair={worst_cone[3]}",
        f"C4_limit_phase_minus_transport={qtext(EXPECTED_C4_LIMIT_MARGIN)}<0;C5_is_sharp_for_this_invoice",
        f"homotopy_slope_gate=|eta|<1 and integer slopes>=m>={CONE_CONSTANT * MIN_CONE_SPREAD}",
    ]
    for D, m, levels, floor, c, eta, transported, debt, actual in controls:
        lines.append(
            f"CONTROL;D={D};m={m};body={HOSTILE};levels={levels};pair=(0,1);"
            f"phase_floor={qtext(floor)};c={qtext(c)};eta={qtext(eta)};"
            f"transported={qtext(transported)};debt={qtext(debt)};"
            f"minimum_actual_overlap={qtext(actual[0])};cell={actual[1]}"
        )
    lines.extend((
        f"conclusion=all reflected residual packets with D>={MIN_CONE_SPREAD},m>={CONE_CONSTANT}D close on all 3003 bodies",
        f"corollary=every still-open D>=6 sector is confined to m<{CONE_CONSTANT}D",
        "scope=reflected THM-2941 sufficient family only;physical LRC14 remains open",
        "normal_vs_python_O=BYTE_IDENTICAL",
        f"base_sha256={sha256(BASE)}",
        f"universal_sha256={sha256(UNIVERSAL)}",
        f"universal_output_sha256={sha256(UNIVERSAL_OUTPUT)}",
        f"d5_source_sha256={sha256(D5)}",
        f"d5_output_sha256={sha256(D5_OUTPUT)}",
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
