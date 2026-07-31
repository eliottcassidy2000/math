#!/usr/bin/env python3
"""Close projected k=3 heights 275 through 272 by exact septimal descent.

The pinned cap-275 package leaves ten atlas bodies at z1=275.  Exact periodic
rays and the common 16-status table leave 65 states on five bodies.  On every
residual body the projected high gate is essential, a duplicate-permitting
two-high upper has a positive exact gap, and literal low enumeration leaves
one fixed low pair and exactly one case per state.

There is a uniform geometric reason all 65 cases close.  The number C of
complete cells missed by the first and two low labels is greater than L/7,
and every residual high denominator d is divisible by seven.  Each residue
class modulo d has at most L/d representatives among the L body cells, so
the set S of occupied residues has |S|>d/7.  Two distinct residues share a
coset modulo d/7; because seven is prime, their difference has exact order
seven.  A primitive high-ray unit preserves that order.  Its two phases are
therefore separated by at least 1/7, and their strict-open radius-1/14 danger
arcs are disjoint, including the endpoint-touching boundary.

The subsequent occupied levels 274, 273, and 272 have no terminal residuals.
Together the 39 atlas rows contain 2,279 states: 370 crude kills, 1,844 exact
common-status/Farkas kills, and the 65 septimal terminal kills.  Height 271
is empty and the next occupied atlas height is 270.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import math
import multiprocessing as mp
import re
from collections import Counter, defaultdict
from fractions import Fraction as Q
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
UPSTREAM_PATH = (
    ROOT / "04-computation/"
    "lrc14_j7_k3_z294_to_z276_ray_status_torsion_descent_thm2941.py"
)
UPSTREAM_OUTPUT_PATH = (
    ROOT / "05-knowledge/results/"
    "lrc14_j7_k3_z294_to_z276_ray_status_torsion_descent_thm2941.out"
)
DEFAULT_OUTPUT = (
    ROOT / "05-knowledge/results/"
    "lrc14_j7_k3_z275_to_z272_septimal_torsion_descent_thm2941.out"
)
EXPECTED_UPSTREAM_SHA256 = (
    "61868bc605b10bf3649dbed0b2323133ff2ed5df8d0c93587f3645975a6c274a"
)
EXPECTED_UPSTREAM_OUTPUT_SHA256 = (
    "a9db04e9bf53709199a8fb18f618ae941e58afb21c7ade2078769ce8313305f4"
)
EXPECTED_PROFILE_SHA256 = None
EXPECTED_SEMANTIC_SHA256 = None
EXPECTED_M_HISTOGRAM = (
    (2, 1), (3, 41), (4, 74), (5, 346), (6, 223), (7, 487),
    (9, 16), (10, 3), (11, 118), (13, 474), (14, 16), (18, 8),
    (22, 4), (24, 8), (39, 3), (55, 7), (110, 9), (117, 2), (165, 4),
)

LEVELS = (275, 274, 273, 272)
EXPECTED_LEVEL_TOTALS = {
    275: (2033, 273, 1695, 65),
    274: (67, 1, 66, 0),
    273: (175, 96, 79, 0),
    272: (4, 0, 4, 0),
}
EXPECTED_TOTALS = (2279, 370, 1844, 65)
INHERITED_LEDGER = 375_713
FINAL_LEDGER = 375_674

EXPECTED_COUNTS = {
    (275, (1, 5, 7, 9, 11, 13)): (426, 0, 407, 19),
    (275, (1, 5, 8, 9, 11, 13)): (477, 239, 238, 0),
    (275, (1, 5, 9, 11, 12, 13)): (34, 34, 0, 0),
    (275, (1, 5, 9, 11, 13, 14)): (215, 0, 211, 4),
    (275, (1, 8, 9, 10, 11, 12)): (0, 0, 0, 0),
    (275, (1, 8, 10, 11, 12, 14)): (119, 0, 113, 6),
    (275, (2, 5, 7, 9, 11, 13)): (429, 0, 401, 28),
    (275, (2, 5, 9, 11, 13, 14)): (215, 0, 207, 8),
    (275, (2, 8, 10, 11, 12, 14)): (118, 0, 118, 0),
    (275, (2, 9, 11, 12, 13, 14)): (0, 0, 0, 0),
    (274, (1, 2, 6, 8, 12, 14)): (3, 0, 3, 0),
    (274, (1, 6, 8, 9, 12, 14)): (1, 0, 1, 0),
    (274, (1, 8, 9, 12, 13, 14)): (0, 0, 0, 0),
    (274, (2, 4, 6, 8, 12, 14)): (2, 0, 2, 0),
    (274, (2, 4, 6, 9, 12, 14)): (2, 0, 2, 0),
    (274, (2, 6, 8, 9, 12, 14)): (34, 1, 33, 0),
    (274, (2, 6, 8, 10, 12, 14)): (25, 0, 25, 0),
    (273, (1, 2, 6, 8, 12, 14)): (1, 1, 0, 0),
    (273, (1, 2, 6, 9, 12, 14)): (3, 1, 2, 0),
    (273, (1, 2, 7, 10, 12, 14)): (2, 1, 1, 0),
    (273, (1, 3, 6, 7, 9, 12)): (1, 1, 0, 0),
    (273, (1, 3, 6, 9, 12, 14)): (2, 2, 0, 0),
    (273, (1, 4, 7, 10, 12, 14)): (1, 1, 0, 0),
    (273, (1, 6, 7, 8, 12, 14)): (3, 1, 2, 0),
    (273, (1, 6, 7, 9, 12, 14)): (33, 22, 11, 0),
    (273, (1, 6, 7, 10, 12, 14)): (40, 22, 18, 0),
    (273, (1, 6, 8, 9, 12, 14)): (1, 1, 0, 0),
    (273, (1, 6, 9, 10, 12, 14)): (4, 3, 1, 0),
    (273, (1, 7, 8, 10, 12, 14)): (35, 16, 19, 0),
    (273, (2, 4, 6, 8, 12, 14)): (2, 1, 1, 0),
    (273, (2, 4, 6, 9, 12, 14)): (1, 0, 1, 0),
    (273, (2, 6, 7, 9, 12, 14)): (18, 10, 8, 0),
    (273, (2, 6, 7, 10, 12, 14)): (22, 10, 12, 0),
    (273, (2, 6, 8, 9, 12, 14)): (4, 2, 2, 0),
    (273, (2, 6, 9, 10, 12, 14)): (1, 0, 1, 0),
    (273, (3, 6, 7, 9, 12, 14)): (1, 1, 0, 0),
    (272, (1, 4, 6, 8, 12, 14)): (1, 0, 1, 0),
    (272, (2, 4, 6, 8, 12, 14)): (2, 0, 2, 0),
    (272, (2, 6, 7, 9, 12, 14)): (1, 0, 1, 0),
}

SURVIVOR_BODIES = (
    (1, 5, 7, 9, 11, 13),
    (1, 5, 9, 11, 13, 14),
    (1, 8, 10, 11, 12, 14),
    (2, 5, 7, 9, 11, 13),
    (2, 5, 9, 11, 13, 14),
)
EXPECTED_GAPS = {
    SURVIVOR_BODIES[0]: Q(55771569587684, 24705788305673475),
    SURVIVOR_BODIES[1]: Q(7112494172759, 2381741545283100),
    SURVIVOR_BODIES[2]: Q(1433340859, 380271992100),
    SURVIVOR_BODIES[3]: Q(3605935026371, 1661178218629930),
    SURVIVOR_BODIES[4]: Q(3723587148503, 1291438601683230),
}
EXPECTED_TERMINAL = {
    SURVIVOR_BODIES[0]: (19, 19, 1),
    SURVIVOR_BODIES[1]: (4, 4, 1),
    SURVIVOR_BODIES[2]: (6, 6, 1),
    SURVIVOR_BODIES[3]: (28, 28, 1),
    SURVIVOR_BODIES[4]: (8, 8, 1),
}
EXPECTED_LOW_PAIRS = {
    SURVIVOR_BODIES[0]: (287, 351),
    SURVIVOR_BODIES[1]: (351, 378),
    SURVIVOR_BODIES[2]: (312, 364),
    SURVIVOR_BODIES[3]: (287, 351),
    SURVIVOR_BODIES[4]: (351, 378),
}
EXPECTED_CLEAN_COUNTS = {
    SURVIVOR_BODIES[0]: 115_050,
    SURVIVOR_BODIES[1]: 240_580,
    SURVIVOR_BODIES[2]: 24_978,
    SURVIVOR_BODIES[3]: 232_068,
    SURVIVOR_BODIES[4]: 244_176,
}
EXPECTED_HIGH_RAY_CHECKS = 2_585_952
EXPECTED_DIRECT_QUALIFYING = ((2, 33), (3, 20), (4, 1), (5, 6), (6, 3), (7, 2))
EXPECTED_DIRECT_EFFECTIVE = ((2, 33), (3, 22), (4, 1), (5, 6), (6, 1), (7, 2))


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def file_sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


require(file_sha256(UPSTREAM_PATH) == EXPECTED_UPSTREAM_SHA256, "cap275 source changed")
require(
    file_sha256(UPSTREAM_OUTPUT_PATH) == EXPECTED_UPSTREAM_OUTPUT_SHA256,
    "cap275 output changed",
)
upstream = load("z275_to_z272_upstream", UPSTREAM_PATH)
z297 = upstream.z297
ATLAS_PATH = upstream.ATLAS_PATH


def atlas_data():
    counts = Counter()
    rows = {level: [] for level in LEVELS}
    high_gate_bodies = []
    pattern = re.compile(r"^row=E=([0-9,]+);.*;z1=([0-9]+);")
    for line in ATLAS_PATH.read_text().splitlines():
        match = pattern.match(line)
        if not match:
            continue
        body = tuple(map(int, match.group(1).split(",")))
        first = int(match.group(2))
        counts[first] += 1
        if first in rows:
            rows[first].append(body)
        if first == 275 and body in SURVIVOR_BODIES:
            require("suffix=HIGH-TAIL:" in line, (body, "high gate disappeared"))
            high_gate_bodies.append(body)
    tasks = tuple((first, body) for first in LEVELS for body in rows[first])
    require(tasks == tuple(EXPECTED_COUNTS), ("atlas body order changed", tasks))
    require(tuple(high_gate_bodies) == SURVIVOR_BODIES, "survivor high-gate order changed")
    require(counts[271] == 0 and counts[270] == 26, ("next frontier changed", counts[271], counts[270]))
    return tasks, tuple(sorted(counts.items()))


def evaluate_task(task):
    first, body = task
    z297.FIRST = first
    z297.ray.FIRST = first
    z297.EXPECTED_COUNTS = {body: EXPECTED_COUNTS[task]}
    return first, z297.evaluate(body)


def analyze_terminal_task(task):
    body, residuals = task
    z297.FIRST = 275
    z297.ray.FIRST = 275
    z297.EXPECTED_GAPS = EXPECTED_GAPS
    z297.EXPECTED_TERMINAL = EXPECTED_TERMINAL
    return z297.analyze_residual_body(body, residuals)


def septimal_certificates(terminal_records):
    """Rebuild the stronger body-mass-to-order-seven implication."""
    z297.FIRST = 275
    z297.ray.FIRST = 275
    certificates = []
    for terminal in terminal_records:
        body, L = terminal[0], terminal[1]
        low_pair = EXPECTED_LOW_PAIRS[body]
        stream = z297.ray.Stream(body)
        clean = z297.fixed_safe_cells(stream, low_pair)
        require(len(clean) == EXPECTED_CLEAN_COUNTS[body], (body, len(clean)))
        require(7 * len(clean) > L, (body, "clean mass not above one seventh"))
        for case in terminal[-1]:
            ds, high_d, low_rows, scalar_excess, _direct_witness = case
            labels = tuple(sorted(label for _d, label in low_rows))
            require(labels == low_pair and scalar_excess >= 0, (body, case))
            require(L % high_d == 0 and high_d % 7 == 0, (body, high_d))
            cell_for_residue = {}
            for cell in clean:
                cell_for_residue.setdefault(cell % high_d, cell)
            residues = tuple(sorted(cell_for_residue))
            fibre_size = L // high_d
            require(len(clean) <= len(residues) * fibre_size, (body, high_d, "fibre bound"))
            require(7 * len(residues) > high_d, (body, high_d, "distinct density"))
            quotient = high_d // 7
            cosets = defaultdict(list)
            for residue in residues:
                cosets[residue % quotient].append(residue)
            crowded = next((values for values in cosets.values() if len(values) >= 2), None)
            require(crowded is not None, (body, high_d, "septimal pigeonhole"))
            first_residue, second_residue = crowded[:2]
            shift = (second_residue - first_residue) % high_d
            order = high_d // math.gcd(high_d, shift)
            require(order == 7, (body, high_d, shift, order))
            require(all(min(unit, 7 - unit) >= 1 for unit in range(1, 7)), "unit gap")
            certificates.append(
                (
                    body,
                    ds,
                    high_d,
                    low_pair,
                    len(clean),
                    len(residues),
                    fibre_size,
                    cell_for_residue[first_residue],
                    cell_for_residue[second_residue],
                    first_residue,
                    second_residue,
                    shift,
                    order,
                )
            )
    require(len(certificates) == 65, ("septimal certificate count", len(certificates)))
    return tuple(certificates)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--processes", type=int, default=8)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    require(args.processes >= 1, "process count must be positive")
    require("projected k3 cap<=275" in UPSTREAM_OUTPUT_PATH.read_text(), "upstream cap changed")
    tasks, atlas_counts = atlas_data()
    controls = z297.ray.local.controls()
    if args.processes == 1:
        records = tuple(evaluate_task(task) for task in tasks)
    else:
        with mp.get_context("spawn").Pool(min(args.processes, len(tasks))) as pool:
            records = tuple(pool.map(evaluate_task, tasks))
    totals = tuple(sum(len(row[index]) for _first, row in records) for index in (7, 8, 9, 10))
    require(totals == EXPECTED_TOTALS, ("global counts changed", totals))
    level_totals = {
        level: tuple(
            sum(len(row[index]) for first, row in records if first == level)
            for index in (7, 8, 9, 10)
        )
        for level in LEVELS
    }
    require(level_totals == EXPECTED_LEVEL_TOTALS, ("level totals changed", level_totals))
    m_histogram = Counter()
    for _first, row in records:
        m_histogram.update(dict(row[11]))
    m_histogram = tuple(sorted(m_histogram.items()))
    if EXPECTED_M_HISTOGRAM is not None:
        require(m_histogram == EXPECTED_M_HISTOGRAM, ("M histogram changed", m_histogram))

    residual_ledger = tuple(
        (row[0], tuple(packet[0] for packet in row[10]))
        for first, row in records
        if first == 275 and row[10]
    )
    require(tuple(body for body, _residuals in residual_ledger) == SURVIVOR_BODIES, "residual bodies changed")
    if args.processes == 1:
        terminal_records = tuple(analyze_terminal_task(task) for task in residual_ledger)
    else:
        with mp.get_context("spawn").Pool(min(args.processes, len(residual_ledger))) as pool:
            terminal_records = tuple(pool.map(analyze_terminal_task, residual_ledger))
    require(sum(row[7] for row in terminal_records) == 65, "zero-high ledger changed")
    require(sum(row[8] for row in terminal_records) == 65, "one-high ledger changed")
    require(sum(row[9] for row in terminal_records) == 5, "low-pair ledger changed")
    require(sum(row[11] for row in terminal_records) == EXPECTED_HIGH_RAY_CHECKS, "high-ray checks changed")
    direct_qualifying = Counter()
    direct_effective = Counter()
    for row in terminal_records:
        direct_qualifying.update(dict(row[12]))
        direct_effective.update(dict(row[13]))
        require(dict(row[10]).get(-1, 0) > 0, (row[0], "negative lows lost"))
    require(tuple(sorted(direct_qualifying.items())) == EXPECTED_DIRECT_QUALIFYING, direct_qualifying)
    require(tuple(sorted(direct_effective.items())) == EXPECTED_DIRECT_EFFECTIVE, direct_effective)
    septimal = septimal_certificates(terminal_records)
    require(INHERITED_LEDGER - len(tasks) == FINAL_LEDGER, "ledger arithmetic changed")

    profile_payload = (
        tasks,
        tuple(
            (
                first,
                row[0],
                row[1],
                row[2],
                row[3],
                row[4],
                row[5],
                row[6],
                tuple(len(row[index]) for index in (7, 8, 9, 10)),
                row[11],
                row[12],
                row[13],
            )
            for first, row in records
        ),
        tuple(row[:-1] for row in terminal_records),
        totals,
        tuple(sorted(level_totals.items())),
        m_histogram,
        atlas_counts,
        controls,
        INHERITED_LEDGER,
        FINAL_LEDGER,
    )
    profile = hashlib.sha256(repr(profile_payload).encode()).hexdigest()
    if EXPECTED_PROFILE_SHA256 is not None:
        require(profile == EXPECTED_PROFILE_SHA256, "profile digest changed")

    # ``records`` come from the repaired MISTAKE-333 engine.  Every rejected
    # common-status instance has its exact dual checked there, but its
    # solver-selected dual basis is deliberately excluded: the stored record
    # is (denominator state, certificate-free problem data).  Consequently
    # this semantic payload is invariant under an equally valid LP basis.
    semantic_payload = (
        tasks,
        records,
        terminal_records,
        septimal,
        totals,
        tuple(sorted(level_totals.items())),
        m_histogram,
        atlas_counts,
        controls,
        INHERITED_LEDGER,
        FINAL_LEDGER,
        EXPECTED_UPSTREAM_OUTPUT_SHA256,
    )
    semantic = hashlib.sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic == EXPECTED_SEMANTIC_SHA256, "semantic digest changed")

    lines = [
        "LRC14 projected k=3 z1=275 through z1=272 septimal torsion/status descent",
        f"upstream_cap275_source_sha256={file_sha256(UPSTREAM_PATH)}",
        f"upstream_cap275_output_sha256={file_sha256(UPSTREAM_OUTPUT_PATH)}",
        f"atlas_output_sha256={file_sha256(ATLAS_PATH)}",
        "scope=all39 exact atlas body rows at occupied heights275,274,273,272;no finite label horizon",
        f"level_totals={level_totals}",
        f"totals=states:{totals[0]};crude_kills:{totals[1]};status_kills:{totals[2]};septimal_residuals:{totals[3]}",
        f"status_M_histogram={dict(m_histogram)}",
        f"canonical_profile_sha256={profile}",
        "independent_exact_farkas_checks=1844/1844:PASS",
        f"pair_overlap_exhaustive_controls={controls[0]}",
        "common_table_controls=coverable positive instance accepted;incompatible hostile instance rejected",
        "z275_scalar_split=zero-high hostile relaxation passes65/65;five positive exact duplicate-permitting >=2-high gaps;therefore exactly one high",
        "z275_finite_reduction=65 states;65 one-high cases;one fixed low pair per survivor body;negative low amplitudes retained",
        f"direct_least_r_histogram={dict(sorted(direct_qualifying.items()))};direct_effective_order_histogram={dict(sorted(direct_effective.items()))}",
        f"high_unit_ray_recurrence_checks={EXPECTED_HIGH_RAY_CHECKS}",
        "uniform_septimal_law=clean-cell count C>L/7 and each residue fibre has at most L/d cells imply |S|>d/7;with7|d two distinct residues share a mod-(d/7) coset and have exact order7",
        "phase_law=primitive high-ray units preserve order7;strict-open radius-1/14 dangers at circular separation1/7 cannot both hold,including endpoint equality",
        "projected_consequence=all65 terminal cases retain full drift-safe mass1>three-aligned cap36/91",
    ]
    for first, row in records:
        (
            body, L, high, first_d, trials, checks, signs, states, crude,
            status, packets, mhist, verified_instance_digest,
            verified_instance_count, representative_instance,
        ) = row
        require(
            verified_instance_count == len(status),
            (body, verified_instance_count, len(status)),
        )
        lines.append(
            f"BODY;z1={first};E={body};L={L};high={high};d1={first_d};trials={trials};checks={checks};"
            f"signs={dict(signs)};states={len(states)};crude={len(crude)};status={len(status)};"
            f"residual={len(packets)};M={dict(mhist)};"
            f"verified_farkas_instance_sha256={verified_instance_digest};"
            f"verified_farkas_checks={verified_instance_count};all_negative=1;"
            f"representative_infeasible_instance={representative_instance};"
            "solver_basis_not_frozen;contradiction_magnitudes_not_frozen"
        )
    for terminal in terminal_records:
        body, L, high, required, residual_count, gap, gap_witness, zero_count, case_count, pair_count, low_signs, high_checks, qualifying, effective, _cases = terminal
        lines.append(
            f"TERMINAL;z1=275;E={body};L={L};high={high};required_suffix={z297.ftext(required)};"
            f"residual_states={residual_count};two_high_gap={z297.ftext(gap)};gap_witness={gap_witness};"
            f"zero_high_passes={zero_count};one_high_cases={case_count};low_pairs={pair_count};"
            f"fixed_low_pair={EXPECTED_LOW_PAIRS[body]};clean_cells={EXPECTED_CLEAN_COUNTS[body]};"
            f"low_sign_census={dict(low_signs)};high_ray_checks={high_checks};"
            f"direct_least_r={dict(qualifying)};direct_effective_orders={dict(effective)}"
        )
    for certificate in septimal:
        lines.append(f"SEPTIMAL;row={certificate}")
    lines.extend(
        (
            "atlas_exact_check=z275:10,z274:7,z273:19,z272:3 body rows;z271 empty;next occupied z270 with26 bodies",
            f"ledger_decrement={INHERITED_LEDGER}-{len(tasks)}={FINAL_LEDGER};decrement counts atlas body rows,not quotient states",
            "conclusion=all projected k3 rows at occupied heights275..272 are empty;with pinned cap275 package,projected k3 cap<=270;next exact frontier=270",
            f"semantic_sha256={semantic}",
            "all_exact_controls=PASS",
        )
    )
    payload = "\n".join(lines) + "\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(payload, encoding="utf-8", newline="\n")
    print(payload, end="")


if __name__ == "__main__":
    main()
