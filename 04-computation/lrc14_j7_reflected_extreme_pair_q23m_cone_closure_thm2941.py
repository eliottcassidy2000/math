#!/usr/bin/env python3
r"""Exact closure of the reflected global-extreme cone ``Q>=23m``.

This verifier works inside THM-2941's sufficient reflected ``k=1`` family.
The arbitrary-level bank already closes 2,442 bodies.  On each of the 561
residual bodies the same-level graph is K6, so a packet not already closed by
a repeated-level pair has six distinct levels.

Fix the physical global-min/global-max slots at levels ``m`` and ``Q=m+D``.
If the endpoint labels are ``a,b`` and ``L=14*lcm(E)``, body safety on an
integer cell gives

    L/14 <= r=(a*j mod L) <= 13L/14-a.

For ``z=mL-a``, the reflected low channel consequently consists of exactly
the ``m`` full teeth

    [(r+kL-L/14)/z, (r+kL+L/14)/z],  k=0,...,m-1.

The ``k=-1`` tooth lies wholly left of zero and the ``k=m`` tooth wholly
right of one.  Each retained tooth has length

    ell=L/[7(mL-a)].

The high channel is a periodic comb of period ``h=L/(QL-b)`` and duty 1/7.
On each low tooth the two endpoint-period cells are the only possible partial
cells.  Therefore, whenever ``ell>2h``, summing the phase-free interval/comb
bound over the ``m`` disjoint low teeth gives

    overlap >= mL/[49(mL-a)] - 2mL/[7(QL-b)].              (1)

Six distinct levels make the worst compatible singleton debt use the four
smallest interior levels ``m+1,...,m+4``, assigned in decreasing-label order.
Subtracting this debt from (1) gives

    M(m,Q) = mL/[49(mL-a)] - a/[7(mL-a)]
             - sum_r e_r/[7((m+r)L-e_r)]
             - (2mL+b)/[7(QL-b)],                         (2)

where ``r=1,...,4`` and the ``e_r`` are decreasing.

The finite base ``(m,Q)=(2,46)`` is positive on all ``561*30`` physical
orientations.  The same invoice is negative at ``(m,Q)=(2,45)`` on its
hostile row.  The all-scale propagation is algebraic, not a finite
extrapolation.  At ``Q=23m``,

    M(m+1,23(m+1))-M(m,23m)

is the sum of the following positive terms:

    6aL/[49(mL-a)((m+1)L-a)],
    e_r L/[7((m+r)L-e_r)((m+r+1)L-e_r)]  (r=1,...,4),
    25bL/[7(23mL-b)(23(m+1)L-b)].

For fixed ``m``, increasing ``Q`` by one increases (2) by

    (2mL+b)L/[7(QL-b)((Q+1)L-b)] > 0.

Finally, at ``Q>=23m``,

    ell-2h = L((Q-14m)L+14a-b)
               / [7(mL-a)(QL-b)] > 0,

and the domain persists because ``h`` decreases with ``Q``.  Thus the exact
base proves the full cone ``m>=2,Q>=23m``, equivalently ``D>=22m``.

Together with the complete ``m=1`` theorem and the corrected cap-three cone,
this confines reflected certificate failure, within the inherited ``D>=6``
stage, to the 561 residual bodies with ``m>=2`` and ``2m<D<22m``.  This is a
sufficient-family theorem, not a physical-survivor census or LRC(14).
"""

from __future__ import annotations

import argparse
import hashlib
from fractions import Fraction as F
from importlib.util import module_from_spec, spec_from_file_location
from math import lcm
from pathlib import Path


HERE = Path(__file__).resolve()
ROOT = next(parent for parent in HERE.parents if (parent / "00-navigation").is_dir())
M1_SOURCE = ROOT / "04-computation/lrc14_j7_reflected_extreme_pair_m1_complete_closure_thm2941.py"
M1_OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_reflected_extreme_pair_m1_complete_closure_thm2941.out"
OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_reflected_extreme_pair_q23m_cone_closure_thm2941.out"

EXPECTED_M1_SOURCE_SHA256 = "691cdad79c7d719cb7c1268c867b919919c46b45ee632d423ac8f8b7f957bcfd"
EXPECTED_M1_OUTPUT_SHA256 = "af1e7868e9e22dc56165ff1f14141c12ff6b7fee43934bdf5a6cdab7de76706d"
EXPECTED_M1_SEMANTIC = "1b19c0356bb6ac3b88ba93c9a775ff2d8db02d4353789d14d93b0441ef5ad2cc"

BODY_COUNT = 561
ORIENTATION_COUNT = 30
BASE_M = 2
RAY_RATIO = 23
BASE_Q = RAY_RATIO * BASE_M
HOSTILE_Q = BASE_Q - 1
TOOTH_CONTROL_LEVELS = (1, 2, 3, 5, 11, 23)
MONOTONICITY_CONTROLS = (2, 3, 4, 5, 11)

H = (1, 2, 3, 4, 6, 12)
EXPECTED_WEAKEST_INVOICE = (
    F(2031150202, 35071508756991), H, (5, 0)
)
EXPECTED_WEAKEST_ACTUAL = (
    F(5826258981442, 569255309820855), H, (5, 2)
)
EXPECTED_HOSTILE = (-F(7500666974, 34308985983447), H, (5, 0))

# Frozen after the first exact replay with body-bound digests.
EXPECTED_TOOTH_DIGEST = "2019f1a46362859bb1bcb5b67cf003ebf01644380255af578fac2d547447d610"
EXPECTED_BASE_DIGEST = "2a5586443d08ef5248ec6300291b32a61b55b0746959f0ea604e37dfa065b730"
EXPECTED_MONOTONICITY_DIGEST = "48a1b0aa0e626ebce69b2bb9f7fe9fbb4350a60782cf01f94ee983d2e3f2efd7"
EXPECTED_SEMANTIC = "88644e8527610dcfd9c1108889a533ceee5883aeb56a84ed0da5edb6891c7e5e"


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


require(sha256(M1_SOURCE) == EXPECTED_M1_SOURCE_SHA256,
        ("m1 source changed", sha256(M1_SOURCE), EXPECTED_M1_SOURCE_SHA256))
require(sha256(M1_OUTPUT) == EXPECTED_M1_OUTPUT_SHA256,
        ("m1 output changed", sha256(M1_OUTPUT), EXPECTED_M1_OUTPUT_SHA256))
require(f"semantic_sha256={EXPECTED_M1_SEMANTIC}" in M1_OUTPUT.read_text(),
        "m1 semantic token missing")
M1 = import_module("extreme_pair_m1_complete", M1_SOURCE)
T = M1.T
B = M1.B


def debt_term(label: int, level: int, ruler: int) -> F:
    require(0 < label < level * ruler, (label, level, ruler))
    return F(label, 7 * (level * ruler - label))


def level_increment_gap(label: int, level: int, ruler: int) -> F:
    """Exact debt decrease when a label moves from level q to q+1."""
    direct = debt_term(label, level, ruler) - debt_term(label, level + 1, ruler)
    closed = F(
        label * ruler,
        7 * (level * ruler - label) * ((level + 1) * ruler - label),
    )
    require(direct == closed > 0, (label, level, ruler, direct, closed))
    return direct


def rearrangement_gap(low_label: int, high_label: int,
                      low_level: int, high_level: int, ruler: int) -> F:
    """Gain from assigning the larger label to the smaller level."""
    require(0 < low_label < high_label < ruler,
            (low_label, high_label, ruler))
    require(0 < low_level < high_level, (low_level, high_level))
    direct = (
        debt_term(high_label, low_level, ruler)
        + debt_term(low_label, high_level, ruler)
        - debt_term(high_label, high_level, ruler)
        - debt_term(low_label, low_level, ruler)
    )
    numerator = (
        ruler * (high_level - low_level) * (high_label - low_label)
        * (low_level * high_level * ruler * ruler - low_label * high_label)
    )
    denominator = 7
    for level, label in (
        (low_level, low_label), (low_level, high_label),
        (high_level, low_label), (high_level, high_label),
    ):
        denominator *= level * ruler - label
    closed = F(numerator, denominator)
    require(direct == closed > 0,
            (low_label, high_label, low_level, high_level, direct, closed))
    return direct


def compatible_debt(body: tuple[int, ...], pair: tuple[int, int],
                    minimum: int, maximum: int):
    """Worst six-distinct-level debt with fixed global extreme slots."""
    require(minimum >= 1 and maximum >= minimum + 5, (minimum, maximum))
    ruler = 14 * lcm(*body)
    i, j = pair
    remaining = tuple(k for k in range(6) if k not in {i, j})
    descending = tuple(sorted(remaining, key=lambda k: body[k], reverse=True))
    levels = [None] * 6
    levels[i] = minimum
    levels[j] = maximum
    for offset, slot in enumerate(descending, 1):
        levels[slot] = minimum + offset
    levels = tuple(levels)
    debt = T.C2.C5.singleton_debt(body, ruler, levels)

    # Exact algebraic controls for both parts of the maximization proof.
    for slot in remaining:
        level_increment_gap(body[slot], minimum + 4, ruler)
    for left in range(len(descending)):
        for right in range(left + 1, len(descending)):
            high_label = body[descending[left]]
            low_label = body[descending[right]]
            rearrangement_gap(
                low_label, high_label,
                minimum + left + 1, minimum + right + 1, ruler,
            )
    return debt, levels, ruler


def full_low_teeth(ruler: int, label: int, cell: int, minimum: int):
    """Recover the exact m full teeth on a body-safe cell."""
    residue = (label * cell) % ruler
    require(ruler // 14 <= residue <= 13 * ruler // 14 - label,
            (ruler, label, cell, residue))
    z = minimum * ruler - label
    radius = F(ruler, 14 * z)
    expected = tuple(
        (F(residue + k * ruler, z) - radius,
         F(residue + k * ruler, z) + radius)
        for k in range(minimum)
    )
    require(expected[0][0] >= 0 and expected[-1][1] <= 1,
            (ruler, label, cell, minimum, expected[0], expected[-1]))

    # The two neighboring teeth are strictly outside [0,1].
    require(F(residue - ruler, z) + radius < 0,
            ("left neighbor", ruler, label, cell, minimum, residue))
    require(F(residue + minimum * ruler, z) - radius > 1,
            ("right neighbor", ruler, label, cell, minimum, residue))
    imported = T.R.reflected_level_arcs(ruler, label, minimum, cell)
    require(imported == expected, (ruler, label, cell, minimum, imported, expected))
    tooth_length = F(ruler, 7 * z)
    require(len(imported) == minimum
            and all(right - left == tooth_length for left, right in imported),
            (ruler, label, cell, minimum, imported, tooth_length))
    return imported, tooth_length, residue


def interval_comb_lower(length: F, period: F) -> F:
    require(length > 2 * period > 0, (length, period))
    return length / 7 - 2 * period / 7


def domain_identity(ruler: int, a: int, b: int,
                    minimum: int, maximum: int) -> F:
    length = F(ruler, 7 * (minimum * ruler - a))
    period = F(ruler, maximum * ruler - b)
    direct = length - 2 * period
    closed = F(
        ruler * ((maximum - 14 * minimum) * ruler + 14 * a - b),
        7 * (minimum * ruler - a) * (maximum * ruler - b),
    )
    require(direct == closed, (ruler, a, b, minimum, maximum, direct, closed))
    return direct


def invoice_data(body: tuple[int, ...], pair: tuple[int, int],
                 minimum: int, maximum: int):
    debt, levels, ruler = compatible_debt(body, pair, minimum, maximum)
    i, j = pair
    a, b = body[i], body[j]
    low_tooth_length = F(ruler, 7 * (minimum * ruler - a))
    high_period = F(ruler, maximum * ruler - b)
    require(domain_identity(ruler, a, b, minimum, maximum) > 0,
            (body, pair, minimum, maximum))
    lower = minimum * interval_comb_lower(low_tooth_length, high_period)
    closed_lower = (
        F(minimum * ruler, 49 * (minimum * ruler - a))
        - F(2 * minimum * ruler, 7 * (maximum * ruler - b))
    )
    require(lower == closed_lower,
            (body, pair, minimum, maximum, lower, closed_lower))
    return lower - debt, lower, debt, levels, ruler


def base_row(body: tuple[int, ...], pair: tuple[int, int]):
    ruler, safe_ranges = T.R.safe_cell_ranges(body)
    cell = safe_ranges[0][0]
    i, j = pair
    low, _, _ = full_low_teeth(ruler, body[i], cell, BASE_M)
    high = T.R.reflected_level_arcs(ruler, body[j], BASE_Q, cell)
    actual = B.intersection_mass(low, high)
    invoice, lower, debt, levels, debt_ruler = invoice_data(
        body, pair, BASE_M, BASE_Q
    )
    require(debt_ruler == ruler, (body, pair, ruler, debt_ruler))
    require(actual >= lower and invoice > 0,
            (body, pair, cell, actual, lower, debt, invoice, levels))
    return (invoice, body, pair, cell, actual - debt,
            actual, lower, debt, levels)


def q_increment(body: tuple[int, ...], pair: tuple[int, int],
                minimum: int, maximum: int):
    now = invoice_data(body, pair, minimum, maximum)[0]
    later = invoice_data(body, pair, minimum, maximum + 1)[0]
    ruler = 14 * lcm(*body)
    b = body[pair[1]]
    closed = F(
        (2 * minimum * ruler + b) * ruler,
        7 * (maximum * ruler - b) * ((maximum + 1) * ruler - b),
    )
    require(later - now == closed > 0,
            (body, pair, minimum, maximum, now, later, closed))
    return closed, body, pair, minimum, maximum


def ray_m_increment(body: tuple[int, ...], pair: tuple[int, int], minimum: int):
    """Exact increment from (m,23m) to (m+1,23(m+1))."""
    maximum = RAY_RATIO * minimum
    now = invoice_data(body, pair, minimum, maximum)[0]
    later = invoice_data(
        body, pair, minimum + 1, RAY_RATIO * (minimum + 1)
    )[0]
    ruler = 14 * lcm(*body)
    i, j = pair
    a, b = body[i], body[j]
    remaining = tuple(sorted(
        (body[k] for k in range(6) if k not in {i, j}), reverse=True
    ))
    low = F(
        6 * a * ruler,
        49 * (minimum * ruler - a) * ((minimum + 1) * ruler - a),
    )
    interior = sum((
        F(
            e * ruler,
            7 * ((minimum + offset) * ruler - e)
            * ((minimum + offset + 1) * ruler - e),
        )
        for offset, e in enumerate(remaining, 1)
    ), F(0))
    high = F(
        25 * b * ruler,
        7 * (RAY_RATIO * minimum * ruler - b)
        * (RAY_RATIO * (minimum + 1) * ruler - b),
    )
    closed = low + interior + high
    require(later - now == closed > 0,
            (body, pair, minimum, now, later, low, interior, high, closed))
    return closed, body, pair, minimum, low, interior, high


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    bodies = T.residual_bodies()
    require(len(bodies) == BODY_COUNT, len(bodies))
    repeated_exceptions = {row[0] for row in T.C2.UNIVERSAL.EXPECTED_EXCEPTIONS}
    require(not (repeated_exceptions & set(bodies)),
            ("same-level exceptions in residual", repeated_exceptions & set(bodies)))

    tooth_body_digests = []
    base_body_digests = []
    monotonicity_body_digests = []
    base_rows = []
    m_increment_minima = {m: None for m in MONOTONICITY_CONTROLS}
    q_increment_minima = {m: None for m in MONOTONICITY_CONTROLS}

    for body in bodies:
        ruler, safe_ranges = T.R.safe_cell_ranges(body)
        cell = safe_ranges[0][0]
        tooth_state = hashlib.sha256()
        for slot, label in enumerate(body):
            for minimum in TOOTH_CONTROL_LEVELS:
                row = full_low_teeth(ruler, label, cell, minimum)
                tooth_state.update(repr((body, slot, minimum, row)).encode())
        tooth_body_digests.append((body, tooth_state.hexdigest()))

        base_state = hashlib.sha256()
        monotonicity_state = hashlib.sha256()
        for i in range(6):
            for j in range(6):
                if i == j:
                    continue
                pair = (i, j)
                row = base_row(body, pair)
                base_rows.append(row)
                base_state.update(repr(row).encode())
                for minimum in MONOTONICITY_CONTROLS:
                    m_row = ray_m_increment(body, pair, minimum)
                    q_row = q_increment(
                        body, pair, minimum, RAY_RATIO * minimum
                    )
                    monotonicity_state.update(repr((m_row, q_row)).encode())
                    if (m_increment_minima[minimum] is None
                            or m_row < m_increment_minima[minimum]):
                        m_increment_minima[minimum] = m_row
                    if (q_increment_minima[minimum] is None
                            or q_row < q_increment_minima[minimum]):
                        q_increment_minima[minimum] = q_row
        base_body_digests.append((body, base_state.hexdigest()))
        monotonicity_body_digests.append((body, monotonicity_state.hexdigest()))

    require(len(base_rows) == BODY_COUNT * ORIENTATION_COUNT, len(base_rows))
    weakest_invoice = min(base_rows, key=lambda row: (row[0], row[1], row[2]))
    weakest_actual = min(base_rows, key=lambda row: (row[4], row[1], row[2]))
    require(weakest_invoice[:3] == EXPECTED_WEAKEST_INVOICE,
            (weakest_invoice, EXPECTED_WEAKEST_INVOICE))
    require((weakest_actual[4], weakest_actual[1], weakest_actual[2])
            == EXPECTED_WEAKEST_ACTUAL,
            (weakest_actual, EXPECTED_WEAKEST_ACTUAL))

    hostile = min(
        (invoice_data(body, (i, j), BASE_M, HOSTILE_Q)[0], body, (i, j))
        for body in bodies for i in range(6) for j in range(6) if i != j
    )
    require(hostile == EXPECTED_HOSTILE and hostile[0] < 0,
            (hostile, EXPECTED_HOSTILE))

    tooth_digest = hashlib.sha256(repr(tuple(tooth_body_digests)).encode()).hexdigest()
    base_digest = hashlib.sha256(repr((
        tuple(base_body_digests), weakest_invoice, weakest_actual,
    )).encode()).hexdigest()
    monotonicity_digest = hashlib.sha256(repr((
        tuple(monotonicity_body_digests),
        tuple(sorted(m_increment_minima.items())),
        tuple(sorted(q_increment_minima.items())),
    )).encode()).hexdigest()
    if EXPECTED_TOOTH_DIGEST is not None:
        require(tooth_digest == EXPECTED_TOOTH_DIGEST,
                (tooth_digest, EXPECTED_TOOTH_DIGEST))
    if EXPECTED_BASE_DIGEST is not None:
        require(base_digest == EXPECTED_BASE_DIGEST,
                (base_digest, EXPECTED_BASE_DIGEST))
    if EXPECTED_MONOTONICITY_DIGEST is not None:
        require(monotonicity_digest == EXPECTED_MONOTONICITY_DIGEST,
                (monotonicity_digest, EXPECTED_MONOTONICITY_DIGEST))

    semantic_image = (
        tuple(bodies), tuple(sorted(repeated_exceptions)),
        BASE_M, BASE_Q, HOSTILE_Q, RAY_RATIO,
        TOOTH_CONTROL_LEVELS, tuple(tooth_body_digests), tooth_digest,
        tuple(base_rows), tuple(base_body_digests), weakest_invoice,
        weakest_actual, base_digest, hostile,
        MONOTONICITY_CONTROLS, tuple(monotonicity_body_digests),
        tuple(sorted(m_increment_minima.items())),
        tuple(sorted(q_increment_minima.items())), monotonicity_digest,
    )
    semantic = hashlib.sha256(repr(semantic_image).encode()).hexdigest()
    if EXPECTED_SEMANTIC is not None:
        require(semantic == EXPECTED_SEMANTIC,
                (semantic, EXPECTED_SEMANTIC))

    lines = [
        "LRC14 reflected physical-extreme Q>=23m cone closure",
        f"universe=residual_bodies:{len(bodies)};physical_orientations:{ORIENTATION_COUNT};same_level_exceptions_disjoint:{tuple(sorted(repeated_exceptions))}",
        "parameters=physical_levels:q_min=m,q_max=Q=m+D;Q_here_is_not_the_reduced_primitive_channel_coordinate",
        "full_tooth_lemma=safe residue L/14<=r<=13L/14-a gives exactly k=0..m-1 full teeth of length L/[7(mL-a)];k=-1 left and k=m right",
        f"full_tooth_controls=levels:{TOOTH_CONTROL_LEVELS};rows:{len(bodies) * 6 * len(TOOTH_CONTROL_LEVELS)};body_bound_digest:{tooth_digest}",
        "interval_comb_lemma=each low tooth loses at most two partial Q-period cells;sum gives mL/[49(mL-a)]-2mL/[7(QL-b)]",
        "domain=ell-2h=L((Q-14m)L+14a-b)/[7(mL-a)(QL-b)]>0 at Q>=23m and persists as Q increases",
        "debt=interior levels m+1..m+4 with decreasing-label assignment;strict level decrease plus exact rearrangement gap",
        f"base=m:{BASE_M};Q:{BASE_Q};rows:{len(base_rows)};weakest_invoice:{qtext(weakest_invoice[0])}@body={weakest_invoice[1]},pair={weakest_invoice[2]},cell={weakest_invoice[3]};actual_margin:{qtext(weakest_invoice[4])};overlap:{qtext(weakest_invoice[5])};lower:{qtext(weakest_invoice[6])};debt:{qtext(weakest_invoice[7])};levels:{weakest_invoice[8]}",
        f"base_weakest_actual={qtext(weakest_actual[4])}@body={weakest_actual[1]},pair={weakest_actual[2]},cell={weakest_actual[3]};invoice:{qtext(weakest_actual[0])}",
        f"base_body_bound_digest={base_digest}",
        f"hostile=m:{BASE_M};Q:{HOSTILE_Q};invoice:{qtext(hostile[0])}@body={hostile[1]},pair={hostile[2]};Q46_is_sharp_for_this_invoice_at_m2",
        "ray_m_increment=6aL/[49(mL-a)((m+1)L-a)]+sum_r e_rL/[7((m+r)L-e_r)((m+r+1)L-e_r)]+25bL/[7(23mL-b)(23(m+1)L-b)]>0",
        "Q_increment=(2mL+b)L/[7(QL-b)((Q+1)L-b)]>0",
        f"monotonicity_controls={MONOTONICITY_CONTROLS};body_bound_digest:{monotonicity_digest}",
    ]
    lines.extend(
        f"m{minimum}_ray_increment_min={qtext(row[0])}@body={row[1]},pair={row[2]}"
        for minimum, row in sorted(m_increment_minima.items())
    )
    lines.extend((
        "analytic_quantifier=the displayed rational identities and positive factors prove every integer m>=2 and Q>=23m;finite controls are hostile replays,not extrapolation",
        "conclusion=within reflected THM-2941 sufficient family,Q>=23m is closed for every m>=2;equivalently D>=22m",
        "corollary=inside the inherited D>=6 stage,certificate failure is confined to 561 bodies with m>=2 and 2m<D<22m",
        "scope=not a physical-survivor census and not LRC14",
        "normal_vs_python_O=BYTE_IDENTICAL",
        f"m1_source_sha256={EXPECTED_M1_SOURCE_SHA256}",
        f"m1_output_sha256={EXPECTED_M1_OUTPUT_SHA256}",
        f"m1_semantic_sha256={EXPECTED_M1_SEMANTIC}",
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
