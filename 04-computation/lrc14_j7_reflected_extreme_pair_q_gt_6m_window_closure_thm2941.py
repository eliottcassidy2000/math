#!/usr/bin/env python3
r"""Close every reflected physical extreme pair with ``Q>6m``.

This verifier works inside THM-2941's sufficient reflected ``k=1`` family.
On a body-safe integer cell, let ``a`` and ``b`` be the labels at the
physical levels ``m`` and ``Q`` and put

    L=14*lcm(E),  z1=mL-a,  z2=QL-b,  alpha=z2/z1.

The low channel has exactly ``m`` full teeth, each of length
``L/(7z1)``.  In units of the high period ``L/z2``, one low tooth is a
window of length ``W=alpha/7``.  A period-one comb with one interval of
length ``1/7`` has exact phase-free mass floor

    J(W)=floor(W)/7 + max(0, frac(W)-6/7).

Thus the physical pair overlap is at least ``mL*J(W)/z2``.  This is an
exact circle-window discrepancy formula, replacing the older loss of two
whole endpoint periods per low tooth.

For ``alpha>=7``, the piecewise formula gives

    J(alpha/7)/alpha >= 1/91,

with equality at ``alpha=13``.  Every six-label body satisfies the sharp
elementary inequality

    L >= 6*sum(E),

with equality only at ``(1,2,3,4,6,12)``.  Since every label is at most
``L/14``, the complete six-distinct-level singleton debt for ``m>=2`` is
at most

    (1/7)[14a/(27L)+14(sum(E)-a)/(41L)] <= 11/1107 < 1/91.

For ``6<alpha<7``, the invoice increases strictly with ``Q`` and reduces
to ``Q=6m+1``.  At that boundary, subtracting the two endpoint debts from
the window lower bound cancels to

    (L-2b)/(7((6m+1)L-b)).

If the four interior labels are ``e_r`` at levels ``m+r``, then

    6((m+r)L-e_r)-((6m+1)L-b)
      =(6r-1)L-6e_r+b > 0,

while

    L-2b-6*sum(e_r)=(L-6*sum(E))+6a+4b > 0.

Consequently the endpoint remainder strictly dominates all four interior
debts.  This closes every integer ``m>=2,Q>6m``, without a gcd restriction.
Together with the already proved rays 4, 5, and 6, the reflected residual is

    m>=2, 3<Q/m<6, Q/m not in {4,5}, gcd(m,Q)<=47.

The abstract ``alpha=6`` window is a hostile sharp control: its length is
exactly the ``6/7`` complement of the high tooth and its phase-free floor is
zero.  The script also checks every breakpoint through alpha 21, the equality
at alpha 13, all 3,003 body inequalities, exact boundary cancellations, and
an independent physical full-tooth engine on a finite control family.
"""

from __future__ import annotations

import argparse
import hashlib
from fractions import Fraction as F
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from math import lcm
from pathlib import Path


HERE = Path(__file__).resolve()
ROOT = next(parent for parent in HERE.parents if (parent / "00-navigation").is_dir())
UPSTREAM_SOURCE = ROOT / "04-computation/lrc14_j7_reflected_extreme_pair_resonance_g48_closure_thm2941.py"
UPSTREAM_OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_reflected_extreme_pair_resonance_g48_closure_thm2941.out"
OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_reflected_extreme_pair_q_gt_6m_window_closure_thm2941.out"

EXPECTED_UPSTREAM_SOURCE_SHA256 = "5ed1a4235466a801022b67342f6ded07cdea8e1449889ac3cc4b68ed5eb1d87f"
EXPECTED_UPSTREAM_OUTPUT_SHA256 = "817047f36413e0ee707bcdb92dd20475b554da740651f387da86587ec281663c"
EXPECTED_UPSTREAM_SEMANTIC = "1fb9c1ee2f0c02e2033b940001b58a76232478a2ae94166e8d1f8cb3d54ae6ac"

ALL_BODY_COUNT = 3003
RESIDUAL_BODY_COUNT = 561
ORIENTATION_COUNT = 30
H = (1, 2, 3, 4, 6, 12)
WINDOW_BOUNDARIES = (6, 7, 13, 14, 20, 21)
ENGINE_LEVEL_CONTROLS = (
    (2, 13), (2, 14), (2, 25), (2, 26), (2, 39), (2, 40), (3, 19),
)

# Frozen after the first complete exact replay.
EXPECTED_BODY_DIGEST = "595ddd1f476c03ff66aeb662efd38b0ad20aded26ad5859b4b4f6667d75e042a"
EXPECTED_BOUNDARY_DIGEST = "93a00c461e7eec97e2044362f6a217515de1ba563b2d0e4e7610ede367641a2d"
EXPECTED_ENGINE_DIGEST = "368fe909b514e1ae1ba79c37c34cd276a55dc08e9f384268ba36c6fbac6ef4a2"
EXPECTED_SEMANTIC = "cac060914df0e16366945480c259e920d6afc8b41038ebcdc4ed19549b72c094"


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


require(sha256(UPSTREAM_SOURCE) == EXPECTED_UPSTREAM_SOURCE_SHA256,
        ("upstream source changed", sha256(UPSTREAM_SOURCE), EXPECTED_UPSTREAM_SOURCE_SHA256))
require(sha256(UPSTREAM_OUTPUT) == EXPECTED_UPSTREAM_OUTPUT_SHA256,
        ("upstream output changed", sha256(UPSTREAM_OUTPUT), EXPECTED_UPSTREAM_OUTPUT_SHA256))
require(f"semantic_sha256={EXPECTED_UPSTREAM_SEMANTIC}" in UPSTREAM_OUTPUT.read_text(),
        "upstream semantic token missing")

U = import_module("reflected_extreme_resonance_base", UPSTREAM_SOURCE)
T = U.T
B = U.B


def circle_window_floor(window: F) -> F:
    """Exact minimum comb mass in a real window, in period-one units."""
    require(window >= 0, window)
    whole = window.numerator // window.denominator
    remainder = window - whole
    return F(whole, 7) + max(F(0), remainder - F(6, 7))


def window_overlap_bound(body: tuple[int, ...], pair: tuple[int, int],
                         minimum: int, maximum: int):
    """Phase-free physical overlap lower bound and its normalized data."""
    require(1 <= minimum < maximum, (minimum, maximum))
    ruler = 14 * lcm(*body)
    a, b = body[pair[0]], body[pair[1]]
    z1 = minimum * ruler - a
    z2 = maximum * ruler - b
    alpha = F(z2, z1)
    window = alpha / 7
    lower = F(minimum * ruler, z2) * circle_window_floor(window)
    return lower, alpha, window, z1, z2, ruler


def coarse_debt_bound(body: tuple[int, ...], pair: tuple[int, int]) -> F:
    """Uniform m>=2 bound: low slot at level 2, all others at level 3."""
    ruler = 14 * lcm(*body)
    a = body[pair[0]]
    return F(1, 7) * (
        F(14 * a, 27 * ruler)
        + F(14 * (sum(body) - a), 41 * ruler)
    )


def boundary_row(body: tuple[int, ...], pair: tuple[int, int], minimum: int):
    """Exact Q=6m+1 cancellation and strict interior comparison."""
    maximum = 6 * minimum + 1
    debt, levels, ruler = U.selected_debt(body, pair, minimum, maximum)
    lower, alpha, window, z1, z2, bound_ruler = window_overlap_bound(
        body, pair, minimum, maximum
    )
    require(bound_ruler == ruler, (body, pair, ruler, bound_ruler))
    a, b = body[pair[0]], body[pair[1]]
    require(F(6, 1) < alpha < F(7, 1), (body, pair, minimum, alpha))
    endpoint_debt = F(a, 7 * z1) + F(b, 7 * z2)
    endpoint_remainder = lower - endpoint_debt
    require(endpoint_remainder == F(ruler - 2 * b, 7 * z2),
            (body, pair, minimum, endpoint_remainder))

    remaining = tuple(k for k in range(6) if k not in set(pair))
    descending = tuple(sorted(remaining, key=lambda k: body[k], reverse=True))
    interiors = tuple(body[k] for k in descending)
    require(tuple(levels[k] for k in descending) == tuple(range(minimum + 1, minimum + 5)),
            (body, pair, minimum, levels, descending))
    comparison_gaps = tuple(
        6 * ((minimum + offset) * ruler - label) - z2
        for offset, label in enumerate(interiors, 1)
    )
    require(min(comparison_gaps) > 0,
            (body, pair, minimum, comparison_gaps))
    numerator_slack = ruler - 2 * b - 6 * sum(interiors)
    require(numerator_slack == ruler - 6 * sum(body) + 6 * a + 4 * b > 0,
            (body, pair, minimum, numerator_slack))
    interior_debt = debt - endpoint_debt
    require(interior_debt < F(6 * sum(interiors), 7 * z2),
            (body, pair, minimum, interior_debt, interiors, z2))
    margin = lower - debt
    certified_slack = F(numerator_slack, 7 * z2)
    require(margin > certified_slack > 0,
            (body, pair, minimum, margin, certified_slack))
    return (
        margin, certified_slack, numerator_slack, min(comparison_gaps),
        body, pair, minimum, maximum, alpha, lower, debt, levels,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    all_bodies = tuple(combinations(range(1, 15), 6))
    require(len(all_bodies) == ALL_BODY_COUNT, len(all_bodies))
    body_digest = hashlib.sha256()
    body_rows = []
    for body in all_bodies:
        ruler = 14 * lcm(*body)
        row = (ruler - 6 * sum(body), body, ruler, sum(body))
        require(row[0] >= 0, ("body inequality failure", row))
        body_rows.append(row)
        body_digest.update(f"{row}\n".encode())
    body_minimum = min(body_rows)
    body_equalities = tuple(row for row in body_rows if row[0] == 0)
    require(body_minimum == (0, H, 168, 28), body_minimum)
    require(body_equalities == (body_minimum,), body_equalities)

    small_lcm = {}
    for body in all_bodies:
        primitive_ruler = lcm(*body)
        if primitive_ruler < 30:
            previous = small_lcm.get(primitive_ruler)
            row = (sum(body), body)
            if previous is None or row > previous:
                small_lcm[primitive_ruler] = row
    require(small_lcm == {
        12: (28, H),
        24: (35, (2, 3, 4, 6, 8, 12)),
    }, small_lcm)

    boundary_values = tuple(
        (alpha, circle_window_floor(F(alpha, 7)))
        for alpha in WINDOW_BOUNDARIES
    )
    require(boundary_values == (
        (6, F(0)), (7, F(1, 7)), (13, F(1, 7)),
        (14, F(2, 7)), (20, F(2, 7)), (21, F(3, 7)),
    ), boundary_values)
    require(circle_window_floor(F(6, 7)) == 0,
            "alpha=6 hostile complement must have zero floor")
    require(circle_window_floor(F(13, 7)) / 13 == F(1, 91),
            "alpha=13 must attain the alpha>=7 normalized floor")
    coarse_gap = F(1, 91) - F(11, 1107)
    require(coarse_gap == F(106, 100737) > 0, coarse_gap)

    residual_bodies = T.residual_bodies()
    require(len(residual_bodies) == RESIDUAL_BODY_COUNT, len(residual_bodies))
    debt_surrogates = []
    for body in residual_bodies:
        for i in range(6):
            for j in range(6):
                if i == j:
                    continue
                pair = (i, j)
                surrogate = coarse_debt_bound(body, pair)
                require(surrogate <= F(11, 1107) < F(1, 91),
                        (body, pair, surrogate))
                debt, _, _ = U.selected_debt(body, pair, 2, 13)
                require(debt <= surrogate,
                        ("debt surrogate failure", body, pair, debt, surrogate))
                debt_surrogates.append((surrogate, debt, body, pair))
    maximum_surrogate = max(debt_surrogates)
    maximum_base_debt = max(
        (debt, body, pair) for _, debt, body, pair in debt_surrogates
    )

    boundary_digest = hashlib.sha256()
    boundary_rows = []
    for minimum in (2, 3, 47):
        for body in residual_bodies:
            for i in range(6):
                for j in range(6):
                    if i == j:
                        continue
                    row = boundary_row(body, (i, j), minimum)
                    boundary_rows.append(row)
                    boundary_digest.update(f"{row}\n".encode())
    boundary_minimum_margin = min(boundary_rows)
    boundary_minimum_numerator = min(
        (row[2], row[4], row[5]) for row in boundary_rows
    )
    boundary_minimum_comparison = min(
        (row[3], row[4], row[5], row[6]) for row in boundary_rows
    )
    require(boundary_minimum_numerator == (14, H, (0, 1)),
            boundary_minimum_numerator)

    engine_digest = hashlib.sha256()
    engine_count = 0
    engine_weakest_gap = None
    for minimum, maximum in ENGINE_LEVEL_CONTROLS:
        for body in residual_bodies:
            ruler, safe_ranges = T.R.safe_cell_ranges(body)
            cell = safe_ranges[0][0]
            for i in range(6):
                for j in range(6):
                    if i == j:
                        continue
                    pair = (i, j)
                    lower, alpha, window, _, _, bound_ruler = window_overlap_bound(
                        body, pair, minimum, maximum
                    )
                    require(bound_ruler == ruler, (body, pair, ruler, bound_ruler))
                    independent = U.full_tooth_overlap(
                        ruler, body[i], minimum, body[j], maximum, cell
                    )
                    imported = B.intersection_mass(
                        T.R.reflected_level_arcs(ruler, body[i], minimum, cell),
                        T.R.reflected_level_arcs(ruler, body[j], maximum, cell),
                    )
                    require(independent == imported >= lower,
                            ("engine/window disagreement", body, pair, minimum,
                             maximum, cell, independent, imported, lower))
                    row = (
                        independent - lower, body, pair, minimum, maximum,
                        cell, independent, lower, alpha, window,
                    )
                    if engine_weakest_gap is None or row < engine_weakest_gap:
                        engine_weakest_gap = row
                    engine_digest.update(f"{row}\n".encode())
                    engine_count += 1
    require(engine_count == (
        RESIDUAL_BODY_COUNT * ORIENTATION_COUNT * len(ENGINE_LEVEL_CONTROLS)
    ), engine_count)

    digests = {
        "body": body_digest.hexdigest(),
        "boundary": boundary_digest.hexdigest(),
        "engine": engine_digest.hexdigest(),
    }
    for label, expected in (
        ("body", EXPECTED_BODY_DIGEST),
        ("boundary", EXPECTED_BOUNDARY_DIGEST),
        ("engine", EXPECTED_ENGINE_DIGEST),
    ):
        if expected is not None:
            require(digests[label] == expected,
                    (label, digests[label], expected))

    semantic_payload = (
        tuple(all_bodies), body_minimum, body_equalities, small_lcm,
        boundary_values, coarse_gap, tuple(residual_bodies),
        maximum_surrogate, maximum_base_debt,
        boundary_minimum_margin, boundary_minimum_numerator,
        boundary_minimum_comparison, ENGINE_LEVEL_CONTROLS,
        engine_count, engine_weakest_gap, digests,
        "m>=2;3<Q/m<6;Q/m not in {4,5};gcd(m,Q)<=47;ray6 boundary closed",
    )
    semantic = hashlib.sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC is not None:
        require(semantic == EXPECTED_SEMANTIC, (semantic, EXPECTED_SEMANTIC))

    lines = [
        "LRC14 reflected extreme-pair Q>6m exact window closure",
        f"all_bodies={len(all_bodies)};sharp_body_inequality=L-6sum(E)>=0;unique_equality:{body_minimum}",
        f"small_lcm_classification={small_lcm}",
        "window_lemma=J(W)=floor(W)/7+max(0,frac(W)-6/7);overlap>=mL*J((QL-b)/(7(mL-a)))/(QL-b)",
        f"window_boundaries={boundary_values};alpha6_hostile_floor:0;alpha13_normalized_equality:1/91",
        f"alpha_ge_7=overlap>1/91;debt<=11/1107;uniform_gap:{qtext(coarse_gap)}",
        f"maximum_body_surrogate={maximum_surrogate}",
        f"maximum_m2_q13_debt={maximum_base_debt}",
        "six_lt_alpha_lt_7=invoice strictly increases in Q;reduce to Q=6m+1",
        "boundary_identity=overlap-endpoint_debt=(L-2b)/(7((6m+1)L-b))",
        "interior_comparison=6((m+r)L-e_r)-((6m+1)L-b)=(6r-1)L-6e_r+b>0",
        "boundary_slack=L-2b-6sum(e_r)=(L-6sum(E))+6a+4b>0",
        f"boundary_control_rows={len(boundary_rows)};minimum_checked_margin:{boundary_minimum_margin};minimum_numerator:{boundary_minimum_numerator};minimum_denominator_gap:{boundary_minimum_comparison}",
        f"independent_engine_controls=levels:{ENGINE_LEVEL_CONTROLS};rows:{engine_count};weakest_actual_minus_window:{engine_weakest_gap}",
        "corollary=every reflected physical extreme pair with integer m>=2,Q>6m closes;no gcd restriction",
        "residual=inside inherited D>=6 stage:561 bodies;m>=2;3<Q/m<6;Q/m not in {4,5};gcd(m,Q)<=47;ray6 is the closed boundary",
        "scope=reflected THM-2941 sufficient family only;not a physical-survivor census and not LRC14",
        f"body_digest={digests['body']}",
        f"boundary_digest={digests['boundary']}",
        f"engine_digest={digests['engine']}",
        "normal_vs_python_O=BYTE_IDENTICAL",
        f"upstream_source_sha256={EXPECTED_UPSTREAM_SOURCE_SHA256}",
        f"upstream_output_sha256={EXPECTED_UPSTREAM_OUTPUT_SHA256}",
        f"upstream_semantic_sha256={EXPECTED_UPSTREAM_SEMANTIC}",
        f"source_sha256={sha256(HERE)}",
        f"semantic_sha256={semantic}",
        "all_exact_controls=PASS",
    ]
    text = "\n".join(lines) + "\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(text, encoding="utf-8", newline="\n")
    print(text, end="")


if __name__ == "__main__":
    main()
