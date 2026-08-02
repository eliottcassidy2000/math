#!/usr/bin/env python3
r"""Complete all-``m`` closure of reflected ``D=1`` level packets.

Let ``E`` be a six-element subset of ``{1,...,14}``, put
``L=14*lcm(E)``, and assign each label one of the two levels ``m,m+1``.
There are exactly

    binom(14,6)*(2^6-2)=3003*62=186186

labelled nonconstant packets.  The adjacent-pair referee pinned below closes
the 133200 packets whose two level classes contain a consecutive-label edge.
This referee closes every one of the remaining 52986 packets, uniformly for
all integers ``m>=1``.

The missing object is the symmetric *difference graph*, not a tournament.
For same-level labels ``a<b``, put ``d=b-a``.  On a body-safe cell ``j``, the
limiting overlap is the tent transform

    phi(d t),  phi(x)=max(0,1/7-||x||).

It is enough to retain a cell on which the relative phase stays wholly in
one linear side of this tent.  Write ``r=dj mod L`` and ``h=floor(aj/L)``.
There are two signed charts.

Positive chart (``r+d<=L/7``), with ``H=floor(dj/L)``:

    A=L/7-r-d/2,
    C=d/2+d*h-a*H-(a+b)/14.

Negative chart (``d<=L-r<=L/7``), with ``H=ceil(dj/L)``:

    A=L/7-(L-r)+d/2,
    C=-d/2-d*h+a*H-(a+b)/14.

Body safety means neither clause has a clipped boundary tooth.  Solving the
two affine endpoint inequalities on each of the ``Q`` teeth and summing gives
the exact same-level overlap

    P_Q = LQ(QA+C)/[(LQ-a)(LQ-b)].                       (1)

If ``Q>=Q0`` and ``Q0*A+C>0``, the denominator in (1) is below ``(LQ)^2``.
Therefore

    P_Q > A/L+C/(LQ)
        >= A/L+min(0,C)/(LQ0) =: F.                     (2)

For a binary packet, ``Q0`` is 1 for a minimum-level pair and 2 for a
higher-level pair.  Its singleton excess

    epsilon_m=sum_e e/[7(q_e L-e)]

decreases term by term with ``m``.  Thus one exact test ``F>epsilon_1``
proves ``P_Q>epsilon_m`` for every ``m>=1``; one-edge Bonferroni then makes
the six-clause union strictly smaller than ``6/7``.

The finite search is complete at a small breakpoint set.  On a fixed safe
range with ``H=floor/ceil(dj/L)`` and ``h=floor(aj/L)`` fixed, ``F`` is affine
decreasing in ``j`` on a positive chart and affine increasing on a negative
chart.  Hence a maximum occurs at a safe-range endpoint or immediately beside
``j=kL/d`` or ``j=kL/a``.  The script includes both range endpoints and a
redundant ``+-2`` integer collar around every such rational breakpoint.  The
``breakpoint_piece_audit`` verifies that the first and last integer in every
monotonicity piece are retained.  The remote tent-support boundary is a
minimum on its signed piece, so it creates no additional maximizer.

An exact census finds 750910 signed chart candidates.  For every available
best labelled-pair chart at ``Q0=1,2``, direct reflected-interval pullback
agrees with (1) and exceeds (2) at ``Q0,Q0+1,Q0+4``: 270240 independent
controls.
Every residual packet has a strict certificate.  The weakest has margin

    765484673861/12987405485040 > 0

at ``E=(1,2,3,4,6,8)``, mask 58, pair ``(6,8)``, ``Q0=1``, ``j=180``.
Together with the pinned adjacent theorem this proves that the reflected
``D=1`` sector is empty for all ``m``.  This is a scoped certificate theorem:
it does not assert that arbitrary spread ``D>=2`` packets close, and failure
of a signed-pair certificate elsewhere would not be a physical survivor.
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
BASE_RELATIVE = Path("04-computation") / "lrc14_j7_reflected_levels_all_q_mass_closure_thm2941.py"
ROOT = next(parent for parent in HERE.parents if (parent / BASE_RELATIVE).is_file())
BASE = ROOT / BASE_RELATIVE
ADJACENT = ROOT / "04-computation" / "lrc14_j7_reflected_adjacent_pair_all_m_closure_thm2941.py"
OUTPUT = ROOT / "05-knowledge" / "results" / "lrc14_j7_reflected_d1_signed_pair_complete_closure_thm2941.out"
EXPECTED_BASE_SHA256 = "2cf0866932f775cc493f97093333e81e65ac3aa76a8e439de969aa700c993f31"
EXPECTED_ADJACENT_SHA256 = "e66cf9712535077a338c2ddb7970aa2e9cf1e59e135413cffcd34416eef61de8"
EXPECTED_ADJACENT_COVERAGE = 133200
EXPECTED_RESIDUAL_COVERAGE = 52986
EXPECTED_TOTAL_COVERAGE = 186186
EXPECTED_SIGNED_CHART_CANDIDATES = 750910
EXPECTED_BREAKPOINT_PIECE_ENDPOINT_CHECKS = 2148769
EXPECTED_DIRECT_FORMULA_FLOOR_CHECKS = 270240
EXPECTED_WEAKEST_MARGIN = F(765484673861, 12987405485040)
EXPECTED_SEMANTIC_SHA256 = "6aedecdb0b663e357ee573b2937bea4d12797a0d9020215684270b9c1b20577e"


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


require(sha256(BASE) == EXPECTED_BASE_SHA256, "reflected interval engine changed")
require(sha256(ADJACENT) == EXPECTED_ADJACENT_SHA256, "adjacent-pair referee changed")
SPEC = spec_from_file_location("reflected_base_for_signed_pair_atlas", BASE)
require(SPEC is not None and SPEC.loader is not None, "cannot load reflected base")
R = module_from_spec(SPEC)
SPEC.loader.exec_module(R)


def qtext(value: F) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def intersection_mass(first, second) -> F:
    i = k = 0
    total = F(0)
    while i < len(first) and k < len(second):
        total += max(F(0), min(first[i][1], second[k][1]) - max(first[i][0], second[k][0]))
        if first[i][1] <= second[k][1]:
            i += 1
        else:
            k += 1
    return total


def candidate_cells(ruler: int, ranges: tuple[tuple[int, int], ...], a: int, b: int) -> tuple[int, ...]:
    """All monotonicity breakpoints relevant to the signed chart floor."""
    d = b - a
    candidates = set()
    for left, right in ranges:
        candidates.update((left, right - 1))
    for denominator, count in ((d, d), (a, a)):
        for k in range(count + 1):
            numerator = k * ruler
            base = numerator // denominator
            for offset in (-2, -1, 0, 1, 2):
                candidates.add(base + offset)
    return tuple(sorted(
        j for j in candidates
        if 0 <= j < ruler and any(left <= j < right for left, right in ranges)
    ))


def ceiling(value: F) -> int:
    return -((-value.numerator) // value.denominator)


def breakpoint_piece_audit(
    ruler: int,
    ranges: tuple[tuple[int, int], ...],
    a: int,
    b: int,
    candidates: tuple[int, ...],
) -> int:
    """Check both integer ends of every ``H,h`` monotonicity piece.

    Support clipping can only shorten such a piece from the remote side of a
    tent.  Positive floors decrease and negative floors increase, so the
    retained near end remains the only possible maximum.
    """
    d = b - a
    candidate_set = set(candidates)
    checks = 0
    rational_breaks = {F(k * ruler, d) for k in range(d + 1)}
    rational_breaks.update(F(k * ruler, a) for k in range(a + 1))
    for left, right in ranges:
        cuts = sorted({F(left), F(right), *(x for x in rational_breaks if left < x < right)})
        for low, high in zip(cuts, cuts[1:]):
            first = max(left, ceiling(low))
            last = min(right - 1, ceiling(high) - 1)
            if first > last:
                continue
            require(first in candidate_set and last in candidate_set,
                    ("breakpoint collar incomplete", ruler, a, b, left, right, low, high,
                     first, last))
            checks += 2 if first != last else 1
    return checks


def signed_chart(ruler: int, a: int, b: int, cell: int):
    """Return sign,A,C,H,h when the whole relative cell lies in the tent."""
    d = b - a
    quotient, remainder = divmod(d * cell, ruler)
    h = (a * cell) // ruler
    if remainder + d <= ruler // 7:
        # theta=(remainder+d*u)/L in [0,1/7].
        H = quotient
        A = F(ruler, 7) - remainder - F(d, 2)
        C = F(d, 2) + d * h - a * H - F(a + b, 14)
        return 1, A, C, H, h
    if remainder and d <= ruler - remainder <= ruler // 7:
        # theta=-(eta-d*u)/L in [-1/7,0].
        eta = ruler - remainder
        H = quotient + 1
        A = F(ruler, 7) - eta + F(d, 2)
        C = -F(d, 2) - d * h + a * H - F(a + b, 14)
        return -1, A, C, H, h
    return None


def overlap_formula(ruler: int, a: int, b: int, level: int, chart) -> F:
    _sign, A, C, _H, _h = chart
    numerator = ruler * level * (level * A + C)
    denominator = (level * ruler - a) * (level * ruler - b)
    require(numerator > 0, ("nonpositive pair numerator", ruler, a, b, level, chart))
    return numerator / denominator


def uniform_floor(ruler: int, first_level: int, chart) -> F:
    _sign, A, C, _H, _h = chart
    # Positive numerator and denominator <(LQ)^2 give this for Q>=Q0.
    return A / ruler + min(F(0), C) / (ruler * first_level)


def epsilon(body: tuple[int, ...], ruler: int, mask: int, m: int = 1) -> F:
    return sum(
        (F(e, 7 * ((m if mask & (1 << i) else m + 1) * ruler - e)) for i, e in enumerate(body)),
        F(0),
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    rows = []
    pair_profiles = {}
    direct_checks = 0
    breakpoint_piece_checks = 0
    candidate_count = 0
    body_count = 0
    adjacent_coverage = 0
    for body in combinations(range(1, 15), 6):
        ruler, ranges = R.safe_cell_ranges(body)
        require(ruler == 14 * lcm(*body), ("ruler mismatch", body, ruler))
        body_count += 1
        for i, k in combinations(range(6), 2):
            a, b = body[i], body[k]
            pair_candidate_cells = candidate_cells(ruler, ranges, a, b)
            breakpoint_piece_checks += breakpoint_piece_audit(
                ruler, ranges, a, b, pair_candidate_cells
            )
            charts = []
            for cell in pair_candidate_cells:
                chart = signed_chart(ruler, a, b, cell)
                if chart is None:
                    continue
                require(R.body_cell_is_safe(ruler, body, cell), ("candidate cell unsafe", body, cell))
                candidate_count += 1
                charts.append((cell, chart))
            for first_level in (1, 2):
                ranked = []
                for cell, chart in charts:
                    floor = uniform_floor(ruler, first_level, chart)
                    if floor <= 0 or first_level * chart[1] + chart[2] <= 0:
                        continue
                    ranked.append((floor, cell, chart))
                best = max(ranked, default=None)
                pair_profiles[(body, i, k, first_level)] = best
                if best is not None:
                    floor, cell, chart = best
                    for level in (first_level, first_level + 1, first_level + 4):
                        actual = intersection_mass(
                            R.reflected_level_arcs(ruler, a, level, cell),
                            R.reflected_level_arcs(ruler, b, level, cell),
                        )
                        formula = overlap_formula(ruler, a, b, level, chart)
                        require(actual == formula,
                                ("signed pair formula mismatch", body, i, k, first_level, cell, chart,
                                 level, actual, formula))
                        require(actual >= floor, ("uniform floor failed", body, i, k, level, actual, floor))
                        direct_checks += 1

        for mask in range(1, 63):
            # The consecutive theorem already handles these; retain them as a
            # positive control, but report the genuinely new no-mono sector.
            mono_consecutive = any(
                a in body and a + 1 in body
                and bool(mask & (1 << body.index(a))) == bool(mask & (1 << body.index(a + 1)))
                for a in range(1, 14)
            )
            if mono_consecutive:
                adjacent_coverage += 1
                continue
            debt = epsilon(body, ruler, mask)
            candidates = []
            for i, k in combinations(range(6), 2):
                if bool(mask & (1 << i)) != bool(mask & (1 << k)):
                    continue
                first_level = 1 if mask & (1 << i) else 2
                profile = pair_profiles[(body, i, k, first_level)]
                if profile is None:
                    continue
                floor, cell, chart = profile
                candidates.append((floor - debt, floor, i, k, first_level, cell, chart))
            best = max(candidates, default=None)
            rows.append((best is not None and best[0] > 0, body, mask, debt, best))

    require(body_count == 3003, ("body census changed", body_count))
    require(adjacent_coverage == EXPECTED_ADJACENT_COVERAGE,
            ("adjacent composition census changed", adjacent_coverage))
    require(len(rows) == EXPECTED_RESIDUAL_COVERAGE,
            ("no-mono mask census changed", len(rows)))
    require(adjacent_coverage + len(rows) == EXPECTED_TOTAL_COVERAGE,
            ("D=1 universe composition changed", adjacent_coverage, len(rows)))
    require(candidate_count == EXPECTED_SIGNED_CHART_CANDIDATES,
            ("signed chart candidate census changed", candidate_count))
    require(breakpoint_piece_checks == EXPECTED_BREAKPOINT_PIECE_ENDPOINT_CHECKS,
            ("breakpoint piece audit changed", breakpoint_piece_checks))
    require(direct_checks == EXPECTED_DIRECT_FORMULA_FLOOR_CHECKS,
            ("direct formula/floor control census changed", direct_checks))
    closed = tuple(row for row in rows if row[0])
    open_rows = tuple(row for row in rows if not row[0])
    weakest = min(closed, key=lambda row: row[4][0]) if closed else None
    require(weakest is not None and weakest[4][0] == EXPECTED_WEAKEST_MARGIN,
            ("weakest signed-pair margin changed", weakest))
    reason_counts = {"NO_SIGNED_CHART": 0, "FLOOR_NOT_ABOVE_EPSILON": 0}
    for _ok, _body, _mask, _debt, best in open_rows:
        reason_counts["NO_SIGNED_CHART" if best is None else "FLOOR_NOT_ABOVE_EPSILON"] += 1

    payload = (
        EXPECTED_BASE_SHA256,
        EXPECTED_ADJACENT_SHA256,
        adjacent_coverage,
        rows,
        candidate_count,
        breakpoint_piece_checks,
        direct_checks,
    )
    semantic = hashlib.sha256(repr(payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256:
        require(semantic == EXPECTED_SEMANTIC_SHA256,
                ("semantic payload changed", semantic, EXPECTED_SEMANTIC_SHA256))
    source_sha = sha256(Path(__file__))
    lines = [
        "LRC14 reflected D=1 signed equal-level pair complete closure exact referee",
        f"upstream_base_sha256={EXPECTED_BASE_SHA256};adjacent_referee_sha256={EXPECTED_ADJACENT_SHA256}",
        f"universe=bodies:{body_count};masks_per_body:62;total:{EXPECTED_TOTAL_COVERAGE};"
        f"adjacent_closed:{adjacent_coverage};signed_residual:{len(rows)}",
        f"signed_chart_candidates={candidate_count};breakpoint_piece_endpoint_checks={breakpoint_piece_checks};"
        f"selected_direct_formula_floor_checks={direct_checks}",
        f"closed={len(closed)};open={len(open_rows)};reason_counts={tuple(sorted(reason_counts.items()))}",
        "positive_chart: A=L/7-r-d/2,C=d/2+d*h-a*H-(a+b)/14",
        "negative_chart: A=L/7-(L-r)+d/2,C=-d/2-d*h+a*(H+1)-(a+b)/14",
        "pair_formula=LQ(QA+C)/((LQ-a)(LQ-b));uniform_floor=A/L+min(0,C)/(LQ0)",
        "all_m_reduction=epsilon_m_termwise_decreasing;denominator<(LQ)^2;one_edge_Bonferroni",
        "D1_COMPLETE_CLOSURE=adjacent:133200+signed_difference:52986=186186;ALL integers m>=1",
    ]
    if weakest is not None:
        _ok, body, mask, debt, best = weakest
        margin, floor, i, k, first_level, cell, chart = best
        lines.append(
            f"weakest_closed_margin={qtext(margin)};E={body};mask={mask};epsilon={qtext(debt)};"
            f"floor={qtext(floor)};pair={(body[i],body[k])};Q0={first_level};j={cell};chart={chart}"
        )
    ranked_open = sorted(open_rows, key=lambda row: row[4][0] if row[4] is not None else -row[3], reverse=True)
    for rank, row in enumerate(ranked_open[:50], 1):
        _ok, body, mask, debt, best = row
        low = tuple(body[i] for i in range(6) if mask & (1 << i))
        high = tuple(e for e in body if e not in low)
        if best is None:
            detail = "reason=NO_SIGNED_CHART"
        else:
            margin, floor, i, k, first_level, cell, chart = best
            detail = (f"reason=FLOOR_NOT_ABOVE_EPSILON;margin={qtext(margin)};floor={qtext(floor)};"
                      f"pair={(body[i],body[k])};Q0={first_level};j={cell};chart={chart}")
        lines.append(
            f"OPEN_RANK={rank};E={body};mask={mask};low={low};high={high};epsilon={qtext(debt)};{detail}"
        )
    lines.extend((f"source_sha256={source_sha}", f"semantic_sha256={semantic}", "all_exact_controls=PASS"))
    output = "\n".join(lines) + "\n"
    args.output.write_text(output, encoding="utf-8", newline="\n")
    print(output, end="")


if __name__ == "__main__":
    main()
