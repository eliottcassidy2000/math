#!/usr/bin/env python3
"""Close the sole projected ``k=3,z1=324`` state by an antipodal H-drift pair.

The exact ray/status frontier leaves one denominator state on

    E=(2,8,10,11,12,14),  L=129360,
    (d1,d2,d3,d4)=(3920,4620,10780,10780).

Exact scalar thresholds force the denominator-4620 and later
denominator-10780 labels to be 364 and 492.  They also rule out every packet
whose denominator-3920 label lies below the projected high wall 12741.  The
remaining label has no scalar cutoff: it may be any high label of exact
denominator 3920.  This is a genuine H-drift tail, not a finite-horizon row.

Two body-safe cells close that whole tail.  The fixed labels 324,364,492 miss
cells j=5880 and k=19600 throughout.  If z has exact denominator 3920, write

    z = 33u + mL,  gcd(u,3920)=1.

Then u is odd and k-j=13720=3*3920+3920/2, so the z-phases on the two cells
differ by u/2=1/2 modulo one, independently of m.  Two open danger arcs of
radius 1/14 whose phases differ by 1/2 are disjoint.  Hence the two local
unions have empty intersection, the all-cell common danger is empty, and the
lossless projected residual has mass one.  This is strictly above the
three-aligned open-union cap 36/91.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
from fractions import Fraction as F
from math import gcd
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FRONTIER_PATH = (
    ROOT / "04-computation" / "lrc14_j7_k3_z324_ray_status_frontier_thm2941.py"
)
FRONTIER_OUTPUT_PATH = (
    ROOT / "05-knowledge" / "results" /
    "lrc14_j7_k3_z324_ray_status_frontier_thm2941.out"
)
PROJECTED_PATH = (
    ROOT / "04-computation" /
    "lrc14_j7_five_aligned_two_drift_projected_closure_thm2941.py"
)
DEFAULT_OUTPUT = (
    ROOT / "05-knowledge" / "results" /
    "lrc14_j7_k3_z324_antipodal_h_drift_closure_thm2941.out"
)
EXPECTED_FRONTIER_SHA256 = (
    "f934aba3b305aa928cdf1099e23ae0a982c8cf6e3003ce923b66b40365d8cb11"
)
EXPECTED_FRONTIER_OUTPUT_SHA256 = (
    "db3c5c68c4aa2f61584ef91dd2171901888270edbf17e860d40f16a64d3a9242"
)
EXPECTED_PROJECTED_SHA256 = (
    "76f891edfcc029a08202481304a809e03e8bd81f247afaeabab685825c4d3662"
)
EXPECTED_SEMANTIC_SHA256 = (
    "613d6daf11b6cd5b4292d43820badbb40a0951bf9e58a38108ceb68cf8e360f9"
)

BODY = (2, 8, 10, 11, 12, 14)
FIRST = 324
DS = (3920, 4620, 10780, 10780)
TAIL_DS = (3920, 4620, 10780)
ALIGNED_THREE_UNION_CAP = F(36, 91)
CELL_J = 5880
CELL_K = 19600
FIXED_LABELS = (324, 364, 492)
REPRESENTATIVE_H = 12771


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def file_sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


if EXPECTED_FRONTIER_SHA256 is not None:
    require(file_sha256(FRONTIER_PATH) == EXPECTED_FRONTIER_SHA256, "frontier changed")
if EXPECTED_FRONTIER_OUTPUT_SHA256 is not None:
    require(
        file_sha256(FRONTIER_OUTPUT_PATH) == EXPECTED_FRONTIER_OUTPUT_SHA256,
        "frontier transcript changed",
    )
require(file_sha256(PROJECTED_PATH) == EXPECTED_PROJECTED_SHA256, "projection engine changed")
frontier = load_module("z324_exact_frontier", FRONTIER_PATH)
P = load_module("z324_lossless_projection", PROJECTED_PATH)
ray = frontier.ray


def first_on_ray(residue, modulus, threshold):
    if residue >= threshold:
        return residue
    return residue + ((threshold - residue + modulus - 1) // modulus) * modulus


def denominator(L, label):
    return L // gcd(L, label)


def unit_rows(stream, d):
    """Exact amplitude and first arbitrary/high point on every unit ray."""
    rows = []
    for unit in range(1, d):
        if gcd(unit, d) != 1:
            continue
        residue = (stream.L // d) * unit
        amplitude = residue * ray.local.delta(stream.carrier, stream.h, residue)
        require(
            (residue + stream.L)
            * ray.local.delta(stream.carrier, stream.h, residue + stream.L)
            == amplitude,
            ("ray recurrence failed", d, unit),
        )
        arbitrary = first_on_ray(residue, stream.L, FIRST + 1)
        high = first_on_ray(residue, stream.L, stream.high_floor)
        rows.append((unit, residue, amplitude, arbitrary, high))
    require(len(rows) > 0, ("empty denominator class", d))
    return tuple(rows)


def best_row(rows, high):
    candidates = []
    for unit, residue, amplitude, arbitrary, high_label in rows:
        label = high_label if high else arbitrary
        candidates.append((amplitude / label, label, unit, residue))
    return max(candidates, key=lambda row: (row[0], -row[1], -row[2]))


def labels_meeting(stream, rows, threshold, *, below_high=False):
    """Enumerate all labels on a denominator class above a positive threshold."""
    require(threshold > 0, ("nonpositive finite threshold", threshold))
    labels = set()
    for _unit, residue, amplitude, arbitrary, _high in rows:
        if amplitude <= 0:
            continue
        label = arbitrary
        while amplitude / label >= threshold:
            if not below_high or label < stream.high_floor:
                labels.add(label)
            if below_high and label >= stream.high_floor:
                break
            label += stream.L
    return tuple(sorted(labels))


def projected_full_mass(cells, L, labels):
    common = ((F(0), F(1)),)
    for cell in cells:
        local_union = P.merge_fraction(
            [
                interval
                for label in labels
                for interval in P.phase_danger(cell, label, L)
            ]
        )
        common = P.intersect_fraction(common, local_union)
        if not common:
            break
    return 1 - P.interval_mass(common), common


def direct_projected_mass(body, L, labels):
    """Independent carrier subtraction followed by direct phi_L splitting."""
    carrier = tuple(
        (F(left, P.A.RULER), F(right, P.A.RULER))
        for left, right in P.A.carrier_for(body)
    )
    removed = P.merge_fraction(
        [
            interval
            for label in labels
            for interval in P.danger_fraction(label)
        ]
    )
    residual = P.subtract_fraction(carrier, removed)
    projected = []
    for left, right in residual:
        scaled_left = L * left
        scaled_right = L * right
        for integer in range(P.floor_fraction(scaled_left), P.ceil_fraction(scaled_right)):
            piece_left = max(scaled_left, F(integer)) - integer
            piece_right = min(scaled_right, F(integer + 1)) - integer
            if piece_left < piece_right:
                projected.append((piece_left, piece_right))
    return P.interval_mass(P.merge_fraction(projected))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    stream = ray.Stream(BODY)
    trials, states, recurrence_checks, signs = ray.ray_quotient_states(stream)
    require(
        (stream.L, stream.high_floor, stream.first_d)
        == (129360, 12741, 10780),
        "body constants changed",
    )
    require(DS in states, "sole denominator state disappeared")
    state = states[DS]
    require(
        state
        == {
            "labels": (324, 12771, 364, 492),
            "total": F(58218833, 5387957190),
            "excess": F(22661393, 280173773880),
        },
        ("maximizing state changed", state),
    )

    rows = {d: unit_rows(stream, d) for d in TAIL_DS}
    best_any = {d: best_row(rows[d], False) for d in TAIL_DS}
    best_high = {d: best_row(rows[d], True) for d in TAIL_DS}
    require(
        best_any
        == {
            3920: (F(8147, 5885880), 429, 13, 429),
            4620: (F(2059, 452760), 364, 13, 364),
            10780: (F(56489, 18563160), 492, 41, 492),
        },
        ("arbitrary ray heads changed", best_any),
    )
    require(
        best_high
        == {
            3920: (F(2449, 35043624), 12771, 387, 12771),
            4620: (F(31639, 272108760), 16828, 601, 16828),
            10780: (F(58687, 481283880), 12756, 1063, 12756),
        },
        ("high ray heads changed", best_high),
    )

    required = stream.lower - stream.first_delta
    low_3920_threshold = required - max(
        best_high[4620][0] + best_any[10780][0],
        best_any[4620][0] + best_high[10780][0],
    )
    label_4620_threshold = required - best_high[3920][0] - best_any[10780][0]
    label_10780_threshold = required - best_high[3920][0] - best_any[4620][0]
    thresholds = (
        low_3920_threshold,
        label_4620_threshold,
        label_10780_threshold,
    )
    require(
        thresholds
        == (
            F(35115043, 12066474420),
            F(89391041, 20012412420),
            F(722933, 244053810),
        ),
        ("scalar thresholds changed", thresholds),
    )
    low_3920 = labels_meeting(
        stream, rows[3920], low_3920_threshold, below_high=True
    )
    eligible_4620 = labels_meeting(stream, rows[4620], label_4620_threshold)
    eligible_10780 = labels_meeting(stream, rows[10780], label_10780_threshold)
    require(
        low_3920 == () and eligible_4620 == (364,) and eligible_10780 == (492,),
        ("literal packet classification changed", low_3920, eligible_4620, eligible_10780),
    )
    base_surplus = (
        stream.first_delta
        + ray.local.delta(stream.carrier, stream.h, 364)
        + ray.local.delta(stream.carrier, stream.h, 492)
        - stream.lower
    )
    require(base_surplus == F(5119, 465404940) > 0, "base surplus changed")

    # Literal denominator and distinctness gates.
    require(
        tuple(denominator(stream.L, label) for label in FIXED_LABELS)
        == (10780, 4620, 10780),
        "fixed denominator ledger changed",
    )
    require(
        denominator(stream.L, REPRESENTATIVE_H) == 3920
        and REPRESENTATIVE_H >= stream.high_floor
        and len(set((*FIXED_LABELS, REPRESENTATIVE_H))) == 4,
        "representative H label changed",
    )

    cells = P.body_cells(P.A.carrier_for(BODY), stream.L)
    require(len(cells) == 42082, ("body-cell count changed", len(cells)))
    require(CELL_J in cells and CELL_K in cells, "antipodal cells left the carrier")
    fixed_miss_ledger = tuple(
        (cell, label, P.phase_danger(cell, label, stream.L))
        for cell in (CELL_J, CELL_K)
        for label in FIXED_LABELS
    )
    require(
        all(intervals == () for _cell, _label, intervals in fixed_miss_ledger),
        ("fixed drift entered an antipodal cell", fixed_miss_ledger),
    )

    difference = CELL_K - CELL_J
    require(difference == 13720 and difference % 3920 == 1960, "cell shift changed")
    units = tuple(unit for unit in range(1, 3920) if gcd(unit, 3920) == 1)
    require(len(units) == 1344 and all(unit % 2 for unit in units), "unit parity changed")
    phase_shifts = tuple((unit * difference) % 3920 for unit in units)
    require(set(phase_shifts) == {1960}, "H phase is not uniformly antipodal")

    # Exact representative instance of the symbolic two-cell argument.
    packet = (*FIXED_LABELS, REPRESENTATIVE_H)
    local_j = P.merge_fraction(
        [interval for label in packet for interval in P.phase_danger(CELL_J, label, stream.L)]
    )
    local_k = P.merge_fraction(
        [interval for label in packet for interval in P.phase_danger(CELL_K, label, stream.L)]
    )
    require(P.intersect_fraction(local_j, local_k) == (), "half-turn sections overlap")
    full_mass, common = projected_full_mass(cells, stream.L, packet)
    direct_mass = direct_projected_mass(BODY, stream.L, packet)
    require(
        common == () and full_mass == direct_mass == 1,
        ("independent projected control failed", common, full_mass, direct_mass),
    )
    require(full_mass > ALIGNED_THREE_UNION_CAP, "projected cap not beaten")

    semantic_payload = (
        BODY,
        FIRST,
        DS,
        stream.L,
        stream.h,
        stream.lower,
        stream.first_delta,
        stream.high_floor,
        state,
        trials,
        recurrence_checks,
        signs,
        best_any,
        best_high,
        thresholds,
        low_3920,
        eligible_4620,
        eligible_10780,
        base_surplus,
        len(cells),
        fixed_miss_ledger,
        difference,
        len(units),
        set(phase_shifts),
        packet,
        full_mass,
        direct_mass,
    )
    semantic_hash = hashlib.sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256, "semantic digest changed")

    lines = [
        "LRC14 projected k=3 z1=324 antipodal H-drift closure",
        f"frontier_source_sha256={file_sha256(FRONTIER_PATH)}",
        f"frontier_output_sha256={file_sha256(FRONTIER_OUTPUT_PATH)}",
        f"projection_source_sha256={file_sha256(PROJECTED_PATH)}",
        f"body={BODY};L={stream.L};h={stream.h};lower={stream.lower};"
        f"first_delta={stream.first_delta};high_floor={stream.high_floor}",
        f"sole_state={DS};maximizing_packet={state['labels']};state_excess={state['excess']}",
        f"best_any={best_any}",
        f"best_high={best_high}",
        f"positive_scalar_thresholds={thresholds}",
        (
            f"literal_classification=low_d3920:{low_3920};"
            f"d4620:{eligible_4620};d10780:{eligible_10780};"
            "remaining_H=every high exact-denominator-3920 label"
        ),
        f"fixed_base_surplus={base_surplus}",
        (
            f"two_cell_certificate=j:{CELL_J};k:{CELL_K};difference:{difference};"
            "fixed_phase_dangers:empty;unit_count:1344;H_phase_shift:1/2"
        ),
        "symbolic_identity=z=33u+mL => phase_k-phase_j=u(k-j)/3920=1/2 mod 1",
        (
            f"projected_residual_mass=1;three_aligned_open_union_cap="
            f"{ALIGNED_THREE_UNION_CAP};margin={1-ALIGNED_THREE_UNION_CAP}"
        ),
        f"representative_packet={packet};full_cell_mass={full_mass};direct_subtraction_mass={direct_mass}",
        "independent_control=full-cell De Morgan agrees with direct carrier subtraction/projection",
        "conclusion=the sole z1=324 status residual is empty uniformly over all H heights",
        f"semantic_sha256={semantic_hash}",
        "all_exact_controls=PASS",
    ]
    payload = "\n".join(lines) + "\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(payload, encoding="utf-8", newline="\n")
    print(payload, end="")


if __name__ == "__main__":
    main()
