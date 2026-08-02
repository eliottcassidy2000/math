#!/usr/bin/env python3
r"""Exact closure of the complete reflected ``m=1`` branch.

This verifier works inside THM-2941's sufficient reflected ``k=1`` family.
The arbitrary-level bank already closes 2,442 bodies.  On each of the 561
residual bodies the same-level graph is K6, so every packet not already
closed by a repeated-level pair has six distinct levels.

For minimum level ``m=1`` and spread ``D``, the physical global-extreme pair
has

    q_i=1,             q_j=D+1=Q,

hence its reduced channel is exactly ``(P,Q,g)=(1,D+1,1)``.  The singleton
debt is maximized by assigning the four smallest available intermediate
levels ``2,3,4,5`` to the remaining labels in decreasing-label order; the
exact rearrangement proof is inherited from the pinned ``(D,m)=(6,1)``
verifier.

Finite head
-----------
The pinned D=6 theorem supplies one maximizing body-safe cell for every one
of the 561*30 physical endpoint orientations.  This verifier reuses that
cell at D=5 (a bounded-spread hostile control) and D=7..44.  Exactly three
rows fail to persist, all on H=(1,2,3,4,6,12) at D=9.  Exhaustive reselection
over H's 88 safe cells repairs them.  Every resulting selected-cell margin
is strictly positive.  D=6 is imported, not recomputed or silently claimed
by this head.

Analytic tail
-------------
On any body-safe integer cell the reflected q=1 arc for label a is one full
interval of length

    ell = L/[7(L-a)].

The q=Q arcs for label b are the intersection with [0,1] of a periodic comb
of period ``h=L/(QL-b)`` and duty cycle 1/7.  Partition the line into period
cells centered on the comb teeth.  An interval I contains every full period
cell except possibly its two endpoint cells, so

    mu(I intersect comb) >= (|I|-2h)/7
                            = |I|/7-2h/7.                 (1)

For the extreme pair, (1) minus the maximal compatible distinct-level debt
is

    L/[49(L-a)] - a/[7(L-a)] - interior_debt
      - (2L+b)/[7(QL-b)].                                (2)

At Q=46 the minimum of (2) over all 561*30 rows is positive.  Moreover

    ell-2h = L((Q-14)L+14a-b)/[7(L-a)(QL-b)] > 0,

and ``h(Q+1)<h(Q)``; hence the endpoint-period argument remains in its
positive-length domain for every later Q.  For fixed body and orientation,
the increment of (2) from Q to Q+1 is exactly

    (2L+b)L / [7(QL-b)((Q+1)L-b)] > 0.                   (3)

Thus Q>=46, equivalently D>=45, is an analytic tail.  The same coarse invoice
is negative at Q=45 on the hostile row, so 46 is sharp for this invoice (not
for the underlying extreme-pair certificate).

Together with the inherited D<=5 and exact D=6 closures, the finite head and
analytic tail close every reflected packet with m=1.  This is still only a
sufficient-family theorem, not a physical-survivor census or LRC(14).
"""

from __future__ import annotations

import argparse
import hashlib
from fractions import Fraction as F
from importlib.util import module_from_spec, spec_from_file_location
from math import lcm
from multiprocessing import get_context
from pathlib import Path


HERE = Path(__file__).resolve()
ROOT = next(parent for parent in HERE.parents if (parent / "00-navigation").is_dir())
BASE = ROOT / "04-computation/lrc14_j7_reflected_extreme_pair_d6_m1_distinct_debt_closure_thm2941.py"
BASE_OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_reflected_extreme_pair_d6_m1_distinct_debt_closure_thm2941.out"
OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_reflected_extreme_pair_m1_complete_closure_thm2941.out"

EXPECTED_BASE_SHA256 = "e605c2e02ba253c687c7fb8c082d543fde9b59ff4212cb4107b85705238a89ed"
EXPECTED_BASE_OUTPUT_SHA256 = "c5fe094fcdf3092f96957c9d80ac2589460f6050072cb490270bca45f2d2c16c"
EXPECTED_BASE_SEMANTIC = "a74dfc2377029e62aa549e7ccfcf84b14348c56b507074efdb92e97fa9ca5dfd"

BODY_COUNT = 561
ORIENTATION_COUNT = 30
HEAD_SPREADS = (5,) + tuple(range(7, 45))
TAIL_Q = 46
TAIL_D = TAIL_Q - 1

EXPECTED_HEAD_DIGEST = "250ef090ff29cc5bffacf205369ae1115aaa1ebb88bba5f14fd62df286a3310a"
EXPECTED_TAIL_DIGEST = "aeaf679cc34288b637658c3c6dd09d9c72090836fd6e08e53a97497cfff25d21"
EXPECTED_TAIL_INVOICE = F(3316914368, 90338710837375)
EXPECTED_TAIL_ACTUAL = F(543172541118, 90338710837375)
EXPECTED_TAIL_DOMAIN = (
    F(16257024, 163556065), (1, 4, 8, 9, 12, 14), (0, 5)
)
EXPECTED_Q45_INVOICE = -F(8992715944, 88374571660375)
EXPECTED_SEMANTIC = "1b19c0356bb6ac3b88ba93c9a775ff2d8db02d4353789d14d93b0441ef5ad2cc"

H = (1, 2, 3, 4, 6, 12)
EXPECTED_FIXED_FAILURES = (
    (-F(32954646, 23410512125), H, (5, 1), 9, 146,
     F(12, 839), F(367789146, 23410512125), (5, 10, 4, 3, 2, 1)),
    (-F(384189062, 302447019875), H, (5, 2), 9, 142,
     F(8, 559), F(4712590062, 302447019875), (5, 4, 10, 3, 2, 1)),
    (-F(3708507072, 10990170183265), H, (5, 4), 9, 151,
     F(4, 279), F(1451463038908, 98911531649385), (5, 4, 3, 2, 10, 1)),
)
EXPECTED_REPAIRS = (
    (F(301879854, 23410512125), H, (5, 1), 9, 152,
     F(24, 839), F(367789146, 23410512125), (5, 10, 4, 3, 2, 1), "reselected"),
    (F(3944211938, 302447019875), H, (5, 2), 9, 152,
     F(16, 559), F(4712590062, 302447019875), (5, 4, 10, 3, 2, 1), "reselected"),
    (F(1384709911612, 98911531649385), H, (5, 4), 9, 150,
     F(8, 279), F(1451463038908, 98911531649385), (5, 4, 3, 2, 10, 1), "reselected"),
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


require(sha256(BASE) == EXPECTED_BASE_SHA256,
        ("D6 source changed", sha256(BASE), EXPECTED_BASE_SHA256))
require(sha256(BASE_OUTPUT) == EXPECTED_BASE_OUTPUT_SHA256,
        ("D6 output changed", sha256(BASE_OUTPUT), EXPECTED_BASE_OUTPUT_SHA256))
require(f"semantic_sha256={EXPECTED_BASE_SEMANTIC}" in BASE_OUTPUT.read_text(),
        "D6 semantic token missing")
B = import_module("extreme_pair_m1_base", BASE)
T = B.T


def compatible_debt(body: tuple[int, ...], pair: tuple[int, int],
                    maximum: int):
    """Worst debt for six distinct levels with endpoints 1 and maximum."""
    require(maximum >= 6, maximum)
    ruler = 14 * lcm(*body)
    i, j = pair
    remaining = tuple(k for k in range(6) if k not in {i, j})
    descending = tuple(sorted(remaining, key=lambda k: body[k], reverse=True))
    levels = [None] * 6
    levels[i] = 1
    levels[j] = maximum
    for slot, level in zip(descending, (2, 3, 4, 5)):
        levels[slot] = level
    levels = tuple(levels)
    return T.C2.C5.singleton_debt(body, ruler, levels), levels, ruler


def finite_head_body(body: tuple[int, ...]):
    """Reuse D6 cells and reselect exactly when their invoice is nonpositive."""
    base_rows = B.audit_body(body)
    ruler = base_rows[0][7]
    require(all(row[7] == ruler for row in base_rows), (body, ruler))
    safe_cells = None
    minima = []
    failures = []
    repairs = []
    row_digest = hashlib.sha256()
    for spread in HEAD_SPREADS:
        maximum = spread + 1  # m=1 pins Q=D+1 exactly.
        selected_rows = []
        for base in base_rows:
            i, j = base[1]
            cell = base[5]
            debt, levels, debt_ruler = compatible_debt(body, (i, j), maximum)
            require(debt_ruler == ruler, (body, ruler, debt_ruler))
            overlap = B.intersection_mass(
                T.R.reflected_level_arcs(ruler, body[i], 1, cell),
                T.R.reflected_level_arcs(ruler, body[j], maximum, cell),
            )
            margin = overlap - debt
            selected = (margin, body, (i, j), spread, cell, overlap,
                        debt, levels, "base-D6-cell")
            if margin <= 0:
                if safe_cells is None:
                    _, safe_ranges = T.R.safe_cell_ranges(body)
                    safe_cells = tuple(
                        c for left, right in safe_ranges for c in range(left, right)
                    )
                best_overlap, best_cell = max(
                    (B.intersection_mass(
                        T.R.reflected_level_arcs(ruler, body[i], 1, c),
                        T.R.reflected_level_arcs(ruler, body[j], maximum, c),
                     ), c)
                    for c in safe_cells
                )
                failures.append((margin, body, (i, j), spread, cell,
                                 overlap, debt, levels))
                selected = (best_overlap - debt, body, (i, j), spread,
                            best_cell, best_overlap, debt, levels, "reselected")
                require(selected[0] > 0, ("unrepaired finite row", selected))
                repairs.append(selected)
            selected_rows.append(selected)
            row_digest.update(repr(selected).encode())
        minima.append((spread, min(selected_rows)))
    return tuple(minima), tuple(failures), tuple(repairs), row_digest.hexdigest()


def periodic_comb_lower(length: F, period: F) -> F:
    """Universal interval/comb lower bound from complete period cells.

    Period cells meeting the two endpoints of the interval have total length
    below ``2*period``.  Every complete intervening cell contributes exactly
    one seventh of its length.  No phase or residue assumption is used.
    """
    require(length > 2 * period > 0, (length, period))
    return length / 7 - 2 * period / 7


def tail_domain_persistence(body: tuple[int, ...], pair: tuple[int, int],
                            q: int):
    """Pin the endpoint-period domain and its all-later-Q persistence.

    The displayed closed form is positive at the audited threshold.  Since
    the period denominator gains ``L`` at every step, the comb period then
    strictly decreases, so ``ell>2h`` holds for every integer ``Q>=q``.
    """
    ruler = 14 * lcm(*body)
    i, j = pair
    a, b = body[i], body[j]
    length = F(ruler, 7 * (ruler - a))
    period = F(ruler, q * ruler - b)
    later_period = F(ruler, (q + 1) * ruler - b)
    domain = length - 2 * period
    closed = F(
        ruler * ((q - 14) * ruler + 14 * a - b),
        7 * (ruler - a) * (q * ruler - b),
    )
    require(domain == closed > 0, (body, pair, q, domain, closed))
    require(0 < later_period < period, (body, pair, q, period, later_period))
    return domain, body, pair, length, period, later_period


def tail_row(body: tuple[int, ...], pair: tuple[int, int], q: int):
    ruler, safe_ranges = T.R.safe_cell_ranges(body)
    cell = safe_ranges[0][0]
    i, j = pair
    a, b = body[i], body[j]
    require(q >= 6, q)

    # Body safety gives L/14 <= (a*cell mod L) <= 13L/14-a.  These are
    # precisely the two non-clipping inequalities for the q=1 reflected arc.
    residue = (a * cell) % ruler
    require(ruler // 14 <= residue <= 13 * ruler // 14 - a,
            (body, pair, cell, ruler, residue))
    low_arcs = T.R.reflected_level_arcs(ruler, a, 1, cell)
    low_length = F(ruler, 7 * (ruler - a))
    require(len(low_arcs) == 1 and T.R.interval_mass(low_arcs) == low_length,
            (body, pair, cell, low_arcs, low_length))

    period = F(ruler, q * ruler - b)
    lower = periodic_comb_lower(low_length, period)
    debt, levels, debt_ruler = compatible_debt(body, pair, q)
    require(debt_ruler == ruler, (body, ruler, debt_ruler))
    actual = B.intersection_mass(
        low_arcs, T.R.reflected_level_arcs(ruler, b, q, cell)
    )
    require(actual >= lower, ("comb discrepancy", body, pair, q,
                              cell, actual, lower, period))

    interior = debt - F(a, 7 * (ruler - a)) - F(b, 7 * (q * ruler - b))
    closed = (
        F(ruler, 49 * (ruler - a))
        - F(a, 7 * (ruler - a))
        - interior
        - F(2 * ruler + b, 7 * (q * ruler - b))
    )
    require(lower - debt == closed,
            (body, pair, q, lower, debt, closed, levels))
    return (actual - debt, closed, body, pair, cell, actual, lower, debt)


def tail_monotonicity(body: tuple[int, ...], pair: tuple[int, int], q: int):
    """Exact positive increment of the analytic invoice for every later q."""
    i, j = pair
    ruler = 14 * lcm(*body)
    b = body[j]
    now = tail_row(body, pair, q)[1]
    later_debt, _, _ = compatible_debt(body, pair, q + 1)

    # Rebuild the next invoice algebraically without requiring its direct
    # overlap control; the periodic-comb lemma already holds at every q.
    a = body[i]
    remaining = sorted(
        (body[k] for k in range(6) if k not in {i, j}), reverse=True
    )
    interior = sum(
        (F(e, 7 * (level * ruler - e))
         for e, level in zip(remaining, (2, 3, 4, 5))),
        F(0),
    )
    later = (
        F(ruler, 49 * (ruler - a))
        - F(a, 7 * (ruler - a))
        - interior
        - F(2 * ruler + b, 7 * ((q + 1) * ruler - b))
    )
    expected = F(
        (2 * ruler + b) * ruler,
        7 * (q * ruler - b) * ((q + 1) * ruler - b),
    )
    later_from_lemma = (
        periodic_comb_lower(
            F(ruler, 7 * (ruler - a)),
            F(ruler, (q + 1) * ruler - b),
        )
        - later_debt
    )
    require(later == later_from_lemma,
            (body, pair, q, later, later_from_lemma, later_debt))
    require(later - now == expected > 0,
            (body, pair, q, now, later, expected, later_debt))
    return now, later, expected


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--processes", type=int, default=8)
    args = parser.parse_args()
    require(args.processes >= 1, args.processes)

    bodies = T.residual_bodies()
    require(len(bodies) == BODY_COUNT, len(bodies))
    repeated_exceptions = {row[0] for row in T.C2.UNIVERSAL.EXPECTED_EXCEPTIONS}
    require(not (repeated_exceptions & set(bodies)),
            ("same-level exceptions in residual", repeated_exceptions & set(bodies)))

    if args.processes == 1:
        body_rows = tuple(finite_head_body(body) for body in bodies)
    else:
        with get_context("fork").Pool(args.processes) as pool:
            body_rows = tuple(pool.map(finite_head_body, bodies, chunksize=1))

    head_minima = tuple(
        (spread, min(dict(packet[0])[spread] for packet in body_rows))
        for spread in HEAD_SPREADS
    )
    failures = tuple(sorted(row for packet in body_rows for row in packet[1]))
    repairs = tuple(sorted(row for packet in body_rows for row in packet[2]))
    body_digest_image = tuple(
        (body, packet[3]) for body, packet in zip(bodies, body_rows)
    )
    head_digest = hashlib.sha256(repr((
        head_minima, failures, repairs, body_digest_image,
    )).encode()).hexdigest()
    require(failures == EXPECTED_FIXED_FAILURES, (failures, EXPECTED_FIXED_FAILURES))
    require(repairs == EXPECTED_REPAIRS, (repairs, EXPECTED_REPAIRS))
    if EXPECTED_HEAD_DIGEST is not None:
        require(head_digest == EXPECTED_HEAD_DIGEST,
                (head_digest, EXPECTED_HEAD_DIGEST))
    require(all(row[1][0] > 0 for row in head_minima), head_minima)

    tail_rows = []
    tail_domains = []
    tail_digest_state = hashlib.sha256()
    for body in bodies:
        for i in range(6):
            for j in range(6):
                if i == j:
                    continue
                row = tail_row(body, (i, j), TAIL_Q)
                tail_rows.append(row)
                tail_domains.append(tail_domain_persistence(
                    body, (i, j), TAIL_Q
                ))
                tail_digest_state.update(repr(row).encode())
                tail_monotonicity(body, (i, j), TAIL_Q)
    require(len(tail_rows) == BODY_COUNT * ORIENTATION_COUNT, len(tail_rows))
    weakest_invoice = min(tail_rows, key=lambda row: (row[1], row[2], row[3]))
    weakest_actual = min(tail_rows, key=lambda row: (row[0], row[2], row[3]))
    weakest_domain = min(tail_domains)
    require(weakest_invoice[1] == EXPECTED_TAIL_INVOICE,
            (weakest_invoice, EXPECTED_TAIL_INVOICE))
    require(weakest_actual[0] == EXPECTED_TAIL_ACTUAL,
            (weakest_actual, EXPECTED_TAIL_ACTUAL))
    require(weakest_domain[:3] == EXPECTED_TAIL_DOMAIN,
            (weakest_domain, EXPECTED_TAIL_DOMAIN))
    tail_digest = tail_digest_state.hexdigest()
    require(tail_digest == EXPECTED_TAIL_DIGEST,
            (tail_digest, EXPECTED_TAIL_DIGEST))

    hostile_q45 = min(
        tail_row(body, (i, j), TAIL_Q - 1)[1]
        for body in bodies for i in range(6) for j in range(6) if i != j
    )
    require(hostile_q45 == EXPECTED_Q45_INVOICE < 0,
            (hostile_q45, EXPECTED_Q45_INVOICE))

    semantic_image = (
        tuple(bodies), tuple(sorted(repeated_exceptions)), HEAD_SPREADS,
        head_minima, failures, repairs, body_digest_image, head_digest,
        TAIL_Q, TAIL_D,
        tuple(tail_rows), weakest_invoice, weakest_actual, tail_digest,
        tuple(tail_domains), weakest_domain, hostile_q45,
    )
    semantic = hashlib.sha256(repr(semantic_image).encode()).hexdigest()
    if EXPECTED_SEMANTIC is not None:
        require(semantic == EXPECTED_SEMANTIC,
                (semantic, EXPECTED_SEMANTIC))

    lines = [
        "LRC14 reflected m=1 complete extreme-pair closure",
        f"universe=residual_bodies:{len(bodies)};physical_orientations:{ORIENTATION_COUNT};same_level_exceptions_disjoint:{tuple(sorted(repeated_exceptions))}",
        "parameter_translation=m=1;q_min=1;q_max=1+D;reduced_channel=(P,Q,g)=(1,D+1,1);therefore Q=D+1 exactly",
        f"finite_head=D=5 hostile control and D=7..44;rows:{len(bodies) * ORIENTATION_COUNT * len(HEAD_SPREADS)};D=6 imported from pinned predecessor",
        f"fixed_cell_failures={failures}",
        f"exact_reselections={repairs}",
    ]
    lines.extend(
        f"head_D{spread}=margin:{qtext(row[0])};body:{row[1]};pair:{row[2]};cell:{row[4]};mode:{row[8]};overlap:{qtext(row[5])};debt:{qtext(row[6])};levels:{row[7]}"
        for spread, row in head_minima
    )
    lines.extend((
        f"head_digest={head_digest}",
        "head_digest_binding=ordered_(body,per_body_digest);body_count=561",
        "periodic_comb_lemma=for period h and duty 1/7,an interval loses at most its two partial endpoint periods,so overlap>=length/7-2h/7;phase-free",
        f"analytic_tail=Q>={TAIL_Q},equivalently D>={TAIL_D};q=1 arc length L/[7(L-a)];h=L/(QL-b)",
        f"tail_domain=ell-2h=L((Q-14)L+14a-b)/[7(L-a)(QL-b)]>0;weakest_at_Q46:{qtext(weakest_domain[0])}@body={weakest_domain[1]},pair={weakest_domain[2]};h(Q+1)<h(Q),so_domain_persists_for_all_Q>=46",
        f"tail_weakest_invoice={qtext(weakest_invoice[1])}@body={weakest_invoice[2]},pair={weakest_invoice[3]},cell={weakest_invoice[4]};actual_margin={qtext(weakest_invoice[0])};overlap={qtext(weakest_invoice[5])};lower={qtext(weakest_invoice[6])};debt={qtext(weakest_invoice[7])}",
        f"tail_weakest_actual={qtext(weakest_actual[0])}@body={weakest_actual[2]},pair={weakest_actual[3]},cell={weakest_actual[4]}",
        f"tail_monotonicity=M(Q+1)-M(Q)=(2L+b)L/[7(QL-b)((Q+1)L-b)]>0;Q45_hostile_invoice={qtext(hostile_q45)};Q46_is_sharp_for_this_invoice",
        f"tail_rows={len(tail_rows)};tail_digest={tail_digest}",
        "assembly=D<=5 inherited bounded-spread closure;D5 hostile control positive;D6 pinned exact extreme-pair theorem;D7..44 finite head;D>=45 analytic tail",
        "conclusion=within the pinned reflected THM-2941 sufficient family,the complete m=1 branch closes on all 3003 bodies",
        "corollary=the reflected certificate-failure wedge is confined to the 561 residual bodies with D>=6 and 2<=m<D/2",
        "scope=not a physical-survivor census and not LRC14",
        "normal_vs_python_O=BYTE_IDENTICAL",
        f"base_source_sha256={EXPECTED_BASE_SHA256}",
        f"base_output_sha256={EXPECTED_BASE_OUTPUT_SHA256}",
        f"base_semantic_sha256={EXPECTED_BASE_SEMANTIC}",
        f"source_sha256={sha256(HERE)}",
        f"semantic_sha256={semantic}",
        "all_exact_controls=PASS",
    ))
    text = "\n".join(lines) + "\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(text, encoding="utf-8", newline="\n")
    print(text, end="")


if __name__ == "__main__":
    main()
