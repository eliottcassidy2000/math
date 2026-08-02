#!/usr/bin/env python3
"""Exact z1=234 boundary compositor with a local Farkas normalization repair.

The inherited status generator has a sound early certificate when a threshold
row has no good status.  The legacy independent auditor rebuilds every earlier
nontrivial threshold and therefore rejects that smaller, already sufficient
certificate.  This script repairs only that verifier interface, without
modifying any promoted dependency.  It then uses a memory-bounded exact
cardinality implementation for the terminal translated-band check.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import multiprocessing as mp
import os
import pickle
import re
import sys
import tempfile
from collections import defaultdict
from fractions import Fraction as F
from math import gcd, lcm
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
BASE_SOURCE = (
    ROOT
    / "04-computation/lrc14_j7_k3_z239_gap238_z237_compositional_descent_thm3061.py"
)
BASE_OUTPUT = (
    ROOT
    / "05-knowledge/results/lrc14_j7_k3_z239_gap238_z237_compositional_descent_thm3061.out"
)
ATLAS = (
    ROOT
    / "05-knowledge/results/lrc14_j7_k3_projected_scalar_body_atlas_thm2941.out"
)
OUTPUT = (
    ROOT
    / "05-knowledge/results/lrc14_j7_k3_z234_direct_farkas_four_two_high_boundary_thm3078.out"
)
BASE_SOURCE_SHA256 = "93a90c48ebed4bc782bca31f378ebb0d4f7ee19ef471a60220fcda3e8927e2fb"
BASE_OUTPUT_SHA256 = "b30961378985a543a86a82681f8556353243703aca9bacf09bb6b4f61648274c"
BASE_SEMANTIC_SHA256 = "2b825c5f12a59048b497b73dba234b79301d78e2db65325c63a580870eaccc88"
ATLAS_SHA256 = "cee82237ce1f51729813b9c916edd3353204c18172abe1d71278dee2c5562eda"
LEVEL = 234
NEIGHBOR_COUNTS = {
    236: (1, 1, 0),
    235: (11, 5, 6),
    234: (381, 330, 51),
    233: (62, 45, 17),
}
EXPECTED_SCREEN = (356804, 155330, 191734, 9740)
EXPECTED_ORDER = (1587, 946, 641, 0)
EXPECTED_TERMINAL = (74, 70, 70, 9028, 1245093, 1245093, 0, 0, 0)
EXPECTED_ROW_ORDER_SHA256 = (
    "47fd0c39966bb04799f458ac6492471f80189b783086ead1fbc435587e14d2ee"
)
EXPECTED_SCREEN_RECORD_SHA256 = (
    "7aa128c9558052010093a36d33ca666b7eaf27aab50033b643a29e6a9a5c00f2"
)
EXPECTED_DIRECT_AUDIT = (16, 191718, 1, 380)
EXPECTED_CARDINALITY_METHODS = (1162904, 82189)
EXPECTED_TERMINAL_RECORD_SHA256 = (
    "cfd70cf409c6938b1af8a431527f8a8abd9c6a8b3e3e50cafe03eec1d81d130f"
)
EXPECTED_SEMANTIC_SHA256 = (
    "749a2292da58925c9653ce919d5bab6b374f64f8b713ac52efa1f18fba74c918"
)
CACHE_SCHEMA = "lrc14-k3-z234-direct-farkas-four-two-high-v1"
RULER = 5_045_040

EXPECTED_SURVIVORS = {
    (1, 5, 6, 9, 12, 14): (
        F(-26554, 71461845),
        (
            (2, 490, 980, 2520),
            5,
            (
                (2, 8820, True, F(0)),
                (490, 324, False, F(1181, 277830)),
                (2520, 1939, True, F(653, 1710198)),
            ),
        ),
    ),
    (1, 5, 9, 11, 12, 14): (
        F(-4325478623, 10049021194212),
        (
            (1260, 10780, 16170, 21560),
            3,
            (
                (1260, 20174, True, F(6091, 88967340)),
                (16170, 19308, True, F(19286, 273184065)),
                (21560, 243, False, F(20378, 4584195)),
            ),
        ),
    ),
    (1, 9, 10, 11, 12, 14): (
        F(-413298076, 464380481565),
        (
            (924, 9702, 10780, 21560),
            5,
            (
                (924, 20370, True, F(157, 1996260)),
                (9702, 260, False, F(2809, 490490)),
                (21560, 19521, True, F(18329, 294610932)),
            ),
        ),
    ),
    (2, 5, 9, 11, 12, 14): (
        F(-138056623, 3349673731404),
        (
            (1260, 10780, 16170, 21560),
            3,
            (
                (1260, 20174, True, F(1789, 22241835)),
                (16170, 19308, True, F(70379, 1092736260)),
                (21560, 243, False, F(6646, 1528065)),
            ),
        ),
    ),
}
EXPECTED_LITERAL_SAFE_MASS = {
    (1, 5, 6, 9, 12, 14): F(191990179, 1000465830),
    (1, 5, 9, 11, 12, 14): F(7803086659, 46609560270),
    (1, 9, 10, 11, 12, 14): F(7485939823, 44226712530),
    (2, 5, 9, 11, 12, 14): F(302562628919, 1794468070395),
}
EXPECTED_SURVIVOR_CASE_SHA256 = {
    (1, 5, 6, 9, 12, 14):
        "e9cacadc16e735389099bf3e9a4d88347d20f4410da45095ee2fbe5eebf23594",
    (1, 5, 9, 11, 12, 14):
        "2eeb9ddacc01c89e8bd896f09eb54d86c2ef54b36b86b58d6890684c15d05868",
    (1, 9, 10, 11, 12, 14):
        "dcc8dd2b7f7fa7d5819e5c319c101e45c0cff02cd3637f56d63b56bd5921f445",
    (2, 5, 9, 11, 12, 14):
        "0569ae0be4016b522b7686302c15817bab59c7568352a3fcb095a79992015bd1",
}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def sha(path):
    payload = path.read_bytes()
    require(b"\r" not in payload.replace(b"\r\n", b""), f"bare CR in {path}")
    return hashlib.sha256(payload.replace(b"\r\n", b"\n")).hexdigest()


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


require(sha(BASE_SOURCE) == BASE_SOURCE_SHA256, "THM-3061 source changed")
require(sha(BASE_OUTPUT) == BASE_OUTPUT_SHA256, "THM-3061 output changed")
require(
    f"semantic_sha256={BASE_SEMANTIC_SHA256}" in BASE_OUTPUT.read_text(),
    "THM-3061 semantic changed",
)
require(sha(ATLAS) == ATLAS_SHA256, "THM-2941 projected body atlas changed")
base = load("thm3061_z234_direct_farkas_boundary_base", BASE_SOURCE)
screen_engine = base.base.thm
eng = screen_engine.eng
legacy_farkas_check = eng.independent_farkas_check

_audit_direct = 0
_audit_legacy = 0
_audit_digest = None


def reset_audit():
    global _audit_direct, _audit_legacy, _audit_digest
    _audit_direct = 0
    _audit_legacy = 0
    _audit_digest = hashlib.sha256()


def direct_aware_farkas_check(q, marginals, capacities, histogram, certificate):
    """Audit the generator's exact early zero-good row, else use legacy logic."""
    global _audit_direct, _audit_legacy
    first_zero_good = None
    for threshold, _count in histogram:
        if threshold <= 0:
            continue
        demand = sum(count for load, count in histogram if load >= threshold)
        good = tuple(int(capacity >= threshold) for capacity in capacities)
        if all(good):
            continue
        if not any(good):
            first_zero_good = (threshold, demand, good)
            break
    if first_zero_good is not None:
        threshold, demand, good = first_zero_good
        expected = ((threshold,), (F(1),), (F(0),) * 5)
        require(certificate == expected, ("direct certificate shape", certificate, expected))
        require(demand > 0, ("direct demand", threshold, demand))
        equality_rows = [
            (1,) * 16,
            *[
                tuple((pattern >> index) & 1 for pattern in range(16))
                for index in range(4)
            ],
        ]
        alpha = certificate[1]
        z = certificate[2]
        slacks = tuple(
            sum(z[row] * equality_rows[row][pattern] for row in range(5))
            - alpha[0] * good[pattern]
            for pattern in range(16)
        )
        contradiction = sum(z[row] * (q, *marginals)[row] for row in range(5))
        contradiction -= alpha[0] * demand
        require(slacks == (F(0),) * 16, ("direct slacks", slacks))
        require(contradiction == -demand < 0, ("direct contradiction", contradiction))
        _audit_direct += 1
        branch = ("direct-zero-good", threshold, demand, slacks, contradiction)
    else:
        contradiction = legacy_farkas_check(
            q, marginals, capacities, histogram, certificate
        )
        _audit_legacy += 1
        branch = ("legacy-full-table", contradiction)
    require(_audit_digest is not None, "audit digest not initialized")
    _audit_digest.update(
        repr((q, marginals, capacities, histogram, certificate, branch)).encode()
        + b"\n"
    )
    return contradiction


def require_reject(function, expected_fragment):
    try:
        function()
    except RuntimeError as exc:
        require(expected_fragment in str(exc), (expected_fragment, str(exc)))
        return str(exc)
    raise RuntimeError(("mutation unexpectedly accepted", expected_fragment))


def farkas_controls():
    """One repaired direct row and one untouched nontrivial legacy row."""
    direct_ds = (72, 882, 980, 8820)
    direct_q = 1470
    direct_marginals, direct_capacities = eng.ray.local.hunter_status_data(
        17640, direct_ds, direct_q
    )
    require(direct_marginals == (1225, 630, 420, 1260), direct_marginals)
    direct_histogram = (
        (0, 210),
        (3, 96),
        (4, 624),
        (5, 144),
        (6, 156),
        (8, 198),
        (10, 12),
        (12, 30),
    )
    direct_certificate = ((12,), (F(1),), (F(0),) * 5)
    reset_audit()
    require(
        direct_aware_farkas_check(
            direct_q,
            direct_marginals,
            direct_capacities,
            direct_histogram,
            direct_certificate,
        )
        == -30,
        "direct control contradiction",
    )
    require((_audit_direct, _audit_legacy) == (1, 0), "direct branch control")
    direct_mutations = (
        require_reject(
            lambda: direct_aware_farkas_check(
                direct_q,
                direct_marginals,
                direct_capacities,
                direct_histogram,
                ((10,), (F(1),), (F(0),) * 5),
            ),
            "direct certificate shape",
        ),
        require_reject(
            lambda: direct_aware_farkas_check(
                direct_q,
                direct_marginals,
                direct_capacities,
                direct_histogram,
                ((12,), (F(2),), (F(0),) * 5),
            ),
            "direct certificate shape",
        ),
    )

    legacy_ds = (2, 140, 882, 980)
    legacy_q = 1260
    legacy_marginals, legacy_capacities = eng.ray.local.hunter_status_data(
        8820, legacy_ds, legacy_q
    )
    require(legacy_marginals == (630, 180, 0, 0), legacy_marginals)
    require(
        legacy_capacities
        == (2, 7, 7, 7, 3, 7, 7, 7, 3, 7, 7, 7, 4, 7, 7, 7),
        legacy_capacities,
    )
    legacy_histogram = ((0, 180), (3, 420), (4, 588), (5, 42), (6, 30))
    legacy_certificate = (
        (3, 4, 5, 6),
        (F(1, 270), F(0), F(0), F(0)),
        (F(0), F(1, 270), F(1, 270), F(1, 270), F(1, 270)),
    )
    reset_audit()
    require(
        direct_aware_farkas_check(
            legacy_q,
            legacy_marginals,
            legacy_capacities,
            legacy_histogram,
            legacy_certificate,
        )
        == -1,
        "legacy control contradiction",
    )
    require((_audit_direct, _audit_legacy) == (0, 1), "legacy branch control")
    legacy_mutations = (
        require_reject(
            lambda: direct_aware_farkas_check(
                legacy_q,
                legacy_marginals,
                legacy_capacities,
                legacy_histogram,
                ((3, 4, 5), legacy_certificate[1], legacy_certificate[2]),
            ),
            "threshold mismatch",
        ),
        require_reject(
            lambda: direct_aware_farkas_check(
                legacy_q,
                legacy_marginals,
                legacy_capacities,
                legacy_histogram,
                (
                    legacy_certificate[0],
                    (-F(1, 270), F(0), F(0), F(0)),
                    legacy_certificate[2],
                ),
            ),
            "negative Farkas alpha",
        ),
        require_reject(
            lambda: direct_aware_farkas_check(
                legacy_q,
                legacy_marginals,
                legacy_capacities,
                legacy_histogram,
                (
                    legacy_certificate[0],
                    legacy_certificate[1],
                    (F(0), F(0), F(1, 270), F(1, 270), F(1, 270)),
                ),
            ),
            "negative Farkas column",
        ),
    )
    return (
        direct_ds,
        direct_q,
        direct_marginals,
        direct_histogram,
        direct_certificate,
        -30,
        direct_mutations,
        legacy_ds,
        legacy_q,
        legacy_marginals,
        legacy_capacities,
        legacy_histogram,
        legacy_certificate,
        -1,
        legacy_mutations,
    )


def atlas_tasks():
    pattern = re.compile(
        r"^row=E=([0-9,]+);.*;L=([0-9]+);high=([0-9]+);z1=([0-9]+);"
    )
    row_lines = tuple(
        line for line in ATLAS.read_text().splitlines() if line.startswith("row=")
    )
    require(len(row_lines) == 6060, ("atlas row-line count", len(row_lines)))
    rows = defaultdict(list)
    tasks = []
    for line in row_lines:
        match = pattern.match(line)
        require(match is not None, ("unparsed atlas row", line))
        body = tuple(map(int, match.group(1).split(",")))
        L, high, first = map(int, match.group(2, 3, 4))
        require(
            high
            == screen_engine.base.WALL.numerator
            * L
            // screen_engine.base.WALL.denominator
            + 1,
            (body, L, high, "high floor"),
        )
        rows[first].append((body, L, high, first < high))
        if first == LEVEL:
            tasks.append((first, body, L, high, first < high))
    derived = {
        z: (
            len(rows[z]),
            sum(row[3] for row in rows[z]),
            sum(not row[3] for row in rows[z]),
        )
        for z in NEIGHBOR_COUNTS
    }
    require(derived == NEIGHBOR_COUNTS, ("neighbor census", derived))
    require(len(tasks) == 381, ("z234 rows", len(tasks)))
    row_order = tuple((task[1], task[2], task[3], task[4]) for task in tasks)
    require(
        hashlib.sha256(repr(row_order).encode()).hexdigest()
        == EXPECTED_ROW_ORDER_SHA256,
        "z234 row-order digest",
    )
    return tuple(tasks), derived


def worker_screen(task):
    global _audit_digest
    reset_audit()
    eng.independent_farkas_check = direct_aware_farkas_check
    try:
        row = screen_engine.evaluate(task)
    finally:
        eng.independent_farkas_check = legacy_farkas_check
    require(_audit_direct + _audit_legacy == row[16], (task, "audit count"))
    return task, (*row, _audit_direct, _audit_legacy, _audit_digest.hexdigest())


def vector_fixed_safe_cells(stream, low_labels):
    chunks = [
        np.arange(left, right, dtype=np.int64)
        for left, right in stream.ranges
        if right > left
    ]
    cells = np.concatenate(chunks) if chunks else np.empty(0, dtype=np.int64)
    clean = np.ones(cells.shape, dtype=np.bool_)
    for label in (stream.first, *low_labels):
        residues = (cells * int(label)) % int(stream.L)
        clean &= 14 * residues >= stream.L
        clean &= 14 * (residues + int(label)) <= 13 * stream.L
    return cells[clean]


def hybrid_translated_certificate(clean_cells, L, d):
    """Exact cardinality certificate, computing residues only near the wall."""
    require(L % d == 0, (L, d, "nondivisor translated modulus"))
    clean_count = int(clean_cells.size)
    per_residue_capacity = L // d
    lower_bound = (clean_count + per_residue_capacity - 1) // per_residue_capacity
    kappa = (d + 6) // 7
    if lower_bound > kappa:
        return (
            d,
            clean_count,
            per_residue_capacity,
            lower_bound,
            None,
            kappa,
            lower_bound - kappa,
            "coarse-cardinality",
            True,
            0,
            None,
            None,
        )
    residues_array = np.unique(clean_cells % d)
    residues = tuple(int(value) for value in residues_array.tolist())
    exact_count = len(residues)
    residue_sha = hashlib.sha256(repr(residues).encode()).hexdigest()
    slack = exact_count - kappa
    if slack > 0:
        return (
            d,
            clean_count,
            per_residue_capacity,
            lower_bound,
            exact_count,
            kappa,
            slack,
            "exact-cardinality",
            True,
            0,
            residue_sha,
            None,
        )
    inherited = screen_engine.translated_certificate(residues, d)
    require(inherited[0] == d and inherited[1] == exact_count, inherited)
    return (
        d,
        clean_count,
        per_residue_capacity,
        lower_bound,
        exact_count,
        kappa,
        slack,
        "max-gap",
        inherited[5],
        inherited[6],
        residue_sha,
        inherited,
    )


def merge_intervals(rows):
    merged = []
    for left, right in sorted(rows):
        if not merged or left > merged[-1][1]:
            merged.append([left, right])
        else:
            merged[-1][1] = max(merged[-1][1], right)
    return tuple((left, right) for left, right in merged)


def danger_intervals(label, ruler):
    require(ruler % (14 * label) == 0, (label, ruler, "inexact danger ruler"))
    step = ruler // label
    radius = ruler // (14 * label)
    return (
        ((0, radius),)
        + tuple(
            (index * step - radius, index * step + radius)
            for index in range(1, label)
        )
        + ((ruler - radius, ruler),)
    )


def body_carrier(body):
    require(all(RULER % (14 * label) == 0 for label in body), body)
    danger = merge_intervals(
        interval
        for label in body
        for interval in danger_intervals(label, RULER)
    )
    carrier = []
    cursor = 0
    for left, right in danger:
        if cursor < left:
            carrier.append((cursor, left))
        cursor = max(cursor, right)
    if cursor < RULER:
        carrier.append((cursor, RULER))
    return tuple(carrier)


def intersect_intervals(first, second):
    rows = []
    left_index = 0
    right_index = 0
    while left_index < len(first) and right_index < len(second):
        left = max(first[left_index][0], second[right_index][0])
        right = min(first[left_index][1], second[right_index][1])
        if left < right:
            rows.append((left, right))
        if first[left_index][1] <= second[right_index][1]:
            left_index += 1
        else:
            right_index += 1
    return tuple(rows)


def exact_maximizer_safe_mass(body, labels):
    ruler = lcm(RULER, *(14 * label for label in labels))
    scale = ruler // RULER
    carrier = tuple(
        (scale * left, scale * right) for left, right in body_carrier(body)
    )
    danger = merge_intervals(
        interval
        for label in labels
        for interval in danger_intervals(label, ruler)
    )
    overlap = intersect_intervals(carrier, danger)
    safe_numerator = sum(right - left for left, right in carrier)
    safe_numerator -= sum(right - left for left, right in overlap)
    return F(safe_numerator, ruler), ruler, len(carrier), len(danger)


def literal_boundary_controls():
    rows = []
    for body, (_gap, witness) in EXPECTED_SURVIVORS.items():
        labels = (LEVEL, *tuple(row[1] for row in witness[2]))
        mass, ruler, carrier_components, danger_components = (
            exact_maximizer_safe_mass(body, labels)
        )
        require(mass == EXPECTED_LITERAL_SAFE_MASS[body] > 0, (body, mass))
        rows.append(
            (body, labels, mass, ruler, carrier_components, danger_components)
        )
    return tuple(rows)


def terminal_probe(task):
    first, body, residual = task
    eng.FIRST = first
    eng.ray.FIRST = first
    stream = eng.ray.Stream(body)
    require(first < stream.high_floor, (body, first, stream.high_floor))
    needed = {d for ds in residual for d in eng.suffix_slots(ds, stream.first_d)}
    low, high, low_signs, ray_checks = eng.build_literal_tables(stream, needed)
    gap, gap_witness = eng.duplicate_two_high_gap(stream, residual, low, high)
    zero = eng.zero_high_scalar_passes(stream, residual, low)
    cases = eng.one_high_cases(stream, residual, low, high)

    needed_by_labels = defaultdict(set)
    for _ds, high_d, low_rows, _excess in cases:
        labels = tuple(sorted(label for _d, label in low_rows))
        needed_by_labels[labels].add(high_d)
    certificate_cache = {}
    for labels in sorted(needed_by_labels):
        cells = vector_fixed_safe_cells(stream, labels)
        require(
            all(
                eng.cell_clean(int(cell), label, stream.L)
                for cell in cells[: min(64, len(cells))]
                for label in (stream.first, *labels)
            ),
            (body, labels, "vector clean sample"),
        )
        for high_d in sorted(needed_by_labels[labels]):
            certificate_cache[(labels, high_d)] = hybrid_translated_certificate(
                cells, stream.L, high_d
            )

    coarse_cases = 0
    exact_cases = 0
    maxgap_cases = 0
    cardinality_cases = 0
    failed_cases = []
    unit_checks = 0
    minimum_slack = None
    digest = hashlib.sha256()
    for ds, high_d, low_rows, excess in cases:
        labels = tuple(sorted(label for _d, label in low_rows))
        certificate = certificate_cache[(labels, high_d)]
        method = certificate[7]
        passed = certificate[8]
        unit_checks += certificate[9]
        minimum_slack = (
            certificate[6]
            if minimum_slack is None
            else min(minimum_slack, certificate[6])
        )
        if method == "coarse-cardinality":
            coarse_cases += 1
            cardinality_cases += 1
        elif method == "exact-cardinality":
            exact_cases += 1
            cardinality_cases += 1
        else:
            maxgap_cases += 1
        if not passed:
            failed_cases.append((ds, high_d, low_rows, excess, certificate))
        digest.update(
            repr((ds, high_d, low_rows, excess, certificate)).encode() + b"\n"
        )
    require(
        coarse_cases + exact_cases + maxgap_cases == len(cases),
        (body, "terminal partition"),
    )
    return (
        first,
        body,
        stream.L,
        stream.high_floor,
        len(residual),
        gap,
        gap_witness,
        len(zero),
        len(cases),
        len(needed_by_labels),
        coarse_cases,
        exact_cases,
        maxgap_cases,
        cardinality_cases,
        len(failed_cases),
        unit_checks,
        minimum_slack,
        digest.hexdigest(),
        tuple(failed_cases),
        gap > 0,
        gap > 0 and not failed_cases,
        low_signs,
        ray_checks,
    )


def worker_terminal(task):
    return task, terminal_probe(task)


def checkpoint_fingerprint():
    return hashlib.sha256(
        repr(
            (
                CACHE_SCHEMA,
                BASE_SOURCE_SHA256,
                BASE_OUTPUT_SHA256,
                BASE_SEMANTIC_SHA256,
                ATLAS_SHA256,
                LEVEL,
                NEIGHBOR_COUNTS,
                EXPECTED_SCREEN,
                EXPECTED_ORDER,
                EXPECTED_TERMINAL,
                EXPECTED_ROW_ORDER_SHA256,
                "direct-first-zero-good-exact-v1",
                "vector-cardinality-lower-bound-v1",
            )
        ).encode()
    ).hexdigest()


CACHE_FINGERPRINT = checkpoint_fingerprint()


def checkpoint_path(directory, label, task):
    task_digest = hashlib.sha256(repr(task).encode()).hexdigest()
    return directory / f"{label}-{CACHE_FINGERPRINT[:16]}-{task_digest}.pickle"


def run_pool(function, tasks, processes, label, directory):
    if not tasks:
        return ()
    directory.mkdir(parents=True, exist_ok=True)
    rows = []
    pending = []
    for task in tasks:
        path = checkpoint_path(directory, label, task)
        if path.exists():
            saved_fingerprint, saved_task, row = pickle.loads(path.read_bytes())
            require(saved_fingerprint == CACHE_FINGERPRINT, (path, "fingerprint mismatch"))
            require(saved_task == task, (path, "task mismatch"))
            rows.append(row)
        else:
            pending.append(task)
    print(
        f"{label}_recovered={len(rows)};pending={len(pending)}",
        file=sys.stderr,
        flush=True,
    )

    def retain(envelope):
        task, row = envelope
        path = checkpoint_path(directory, label, task)
        temporary = path.with_name(f"{path.name}.tmp-{os.getpid()}")
        temporary.write_bytes(
            pickle.dumps((CACHE_FINGERPRINT, task, row), protocol=5)
        )
        os.replace(temporary, path)
        rows.append(row)
        print(
            f"{label}_progress={len(rows)}/{len(tasks)};z1={task[0]};body={task[1]}",
            file=sys.stderr,
            flush=True,
        )

    if processes == 1:
        for task in pending:
            retain(function(task))
    elif pending:
        with mp.get_context("spawn").Pool(min(processes, len(pending))) as pool:
            for envelope in pool.imap_unordered(function, pending, chunksize=1):
                retain(envelope)
    return tuple(rows)


def four_totals(rows):
    return tuple(sum(row[index] for row in rows) for index in (9, 10, 11, 12))


def terminal_totals(rows):
    return (
        len(rows),
        sum(row[19] for row in rows),
        sum(row[20] for row in rows),
        sum(row[7] for row in rows),
        sum(row[8] for row in rows),
        sum(row[13] for row in rows),
        sum(row[12] for row in rows),
        sum(row[14] for row in rows),
        sum(row[15] for row in rows),
    )


def ftext(value):
    return f"{value.numerator}/{value.denominator}"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--processes", type=int, default=4)
    parser.add_argument("--terminal-processes", type=int, default=2)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument(
        "--checkpoint-dir",
        type=Path,
        default=Path(tempfile.gettempdir()) / "lrc14-k3-z234-thm3078-v1",
    )
    args = parser.parse_args()

    controls = farkas_controls()
    literal_controls = literal_boundary_controls()
    tasks, neighbor_census = atlas_tasks()
    scheduled = tuple(sorted(tasks, key=lambda task: (-task[2], task[1])))
    records = run_pool(
        worker_screen, scheduled, args.processes, "screen", args.checkpoint_dir
    )
    by_key = {(row[0], row[1]): row for row in records}
    require(len(by_key) == len(tasks), "lost screen row")
    records = tuple(by_key[(task[0], task[1])] for task in tasks)
    base_records = tuple(row[:19] for row in records)
    require(
        hashlib.sha256(repr(base_records).encode()).hexdigest()
        == EXPECTED_SCREEN_RECORD_SHA256,
        "screen record digest",
    )
    totals = four_totals(records)
    order = four_totals(tuple(row for row in records if not row[5]))
    require(totals == EXPECTED_SCREEN, ("screen", totals))
    require(order == EXPECTED_ORDER, ("order screen", order))

    direct_certificates = sum(row[19] for row in records)
    legacy_certificates = sum(row[20] for row in records)
    direct_rows = sum(row[19] > 0 for row in records)
    legacy_only_rows = len(records) - direct_rows
    direct_audit = (
        direct_certificates,
        legacy_certificates,
        direct_rows,
        legacy_only_rows,
    )
    require(
        direct_certificates + legacy_certificates == totals[2],
        ("status audit total", direct_audit),
    )
    require(direct_rows == 1 and legacy_only_rows == 380, direct_audit)
    if EXPECTED_DIRECT_AUDIT is not None:
        require(direct_audit == EXPECTED_DIRECT_AUDIT, direct_audit)
    direct_record = next(row for row in records if row[19])
    require(
        (direct_record[1], direct_record[19], direct_record[20])
        == ((2, 4, 9, 10, 12, 14), 16, 1161),
        ("direct row", direct_record[1], direct_record[19:21]),
    )

    terminal_tasks = tuple(
        (row[0], row[1], row[13]) for row in records if row[5] and row[12]
    )
    require(len(terminal_tasks) == 74, ("residual bodies", len(terminal_tasks)))
    terminals = run_pool(
        worker_terminal,
        terminal_tasks,
        args.terminal_processes,
        "terminal",
        args.checkpoint_dir,
    )
    terminal_by_key = {(row[0], row[1]): row for row in terminals}
    require(len(terminal_by_key) == len(terminal_tasks), "lost terminal row")
    terminals = tuple(
        terminal_by_key[(task[0], task[1])] for task in terminal_tasks
    )
    terminal_summary = terminal_totals(terminals)
    require(terminal_summary == EXPECTED_TERMINAL, terminal_summary)
    cardinality_methods = (
        sum(row[10] for row in terminals),
        sum(row[11] for row in terminals),
    )
    require(cardinality_methods == EXPECTED_CARDINALITY_METHODS,
            ("cardinality methods", cardinality_methods))
    require(sum(cardinality_methods) == 1245093, "cardinality method total")

    survivors = tuple(row for row in terminals if not row[20])
    require(tuple(row[1] for row in survivors) == tuple(EXPECTED_SURVIVORS), survivors)
    for row in survivors:
        expected_gap, expected_witness = EXPECTED_SURVIVORS[row[1]]
        require(row[5] == expected_gap, (row[1], row[5], expected_gap))
        require(row[6] == expected_witness, (row[1], row[6], expected_witness))
        require(
            row[17] == EXPECTED_SURVIVOR_CASE_SHA256[row[1]],
            (row[1], "survivor case digest", row[17]),
        )
        require(row[14] == 0 and row[13] == row[8], (row[1], "one-high closure"))
    positive_gaps = tuple(row[5] for row in terminals if row[19])
    require(len(positive_gaps) == 70, "positive-gap census")
    minimum_positive_gap = min(positive_gaps)
    minimum_row = next(row for row in terminals if row[5] == minimum_positive_gap)
    require(
        (minimum_positive_gap, minimum_row[1])
        == (F(13057, 124176780), (1, 6, 9, 10, 12, 14)),
        (minimum_positive_gap, minimum_row[1]),
    )

    terminal_record_sha = hashlib.sha256(repr(terminals).encode()).hexdigest()
    require(
        terminal_record_sha == EXPECTED_TERMINAL_RECORD_SHA256,
        ("terminal record digest", terminal_record_sha),
    )
    semantic_packet = (
        tasks,
        neighbor_census,
        records,
        terminals,
        totals,
        order,
        direct_audit,
        terminal_summary,
        survivors,
        minimum_positive_gap,
        controls,
        literal_controls,
        CACHE_FINGERPRINT,
    )
    semantic = hashlib.sha256(repr(semantic_packet).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic == EXPECTED_SEMANTIC_SHA256, ("semantic", semantic))

    lines = [
        "LRC14 projected k3 z1=234 direct-Farkas and four-two-high exact boundary",
        f"dependency=THM3061_source:{BASE_SOURCE_SHA256};output:{BASE_OUTPUT_SHA256};semantic:{BASE_SEMANTIC_SHA256};atlas:{ATLAS_SHA256}",
        f"checkpoint_fingerprint={CACHE_FINGERPRINT}",
        "neighbor_census="
        + ";".join(
            f"z{z}:rows{values[0]}/wall{values[1]}/order{values[2]}"
            for z, values in sorted(neighbor_census.items(), reverse=True)
        ),
        f"universe=occupied_rows:{len(tasks)};wall:{sum(t[4] for t in tasks)};order:{sum(not t[4] for t in tasks)};row_order_sha256:{EXPECTED_ROW_ORDER_SHA256}",
        f"screen=states:{totals[0]};crude:{totals[1]};status:{totals[2]};residual:{totals[3]};residual_bodies:{len(terminal_tasks)};screen_record_sha256:{EXPECTED_SCREEN_RECORD_SHA256}",
        f"order_screen=states:{order[0]};crude:{order[1]};status:{order[2]};residual:{order[3]}",
        f"farkas_audit=direct_certificates:{direct_audit[0]};legacy_certificates:{direct_audit[1]};direct_rows:{direct_audit[2]};legacy_only_rows:{direct_audit[3]}",
        "direct_control=body:(2,4,9,10,12,14);ds:(72,882,980,8820);q:1470;M:12;threshold:12;demand:30;slacks:0x16;contradiction:-30;mutations:2",
        "legacy_control=ds:(2,140,882,980);q:1260;M:7;thresholds:(3,4,5,6);contradiction:-1;mutations:3",
        "literal_maximizer_controls="
        + ";".join(
            f"E={row[0]}/labels={row[1]}/safe_mass={ftext(row[2])}/ruler={row[3]}/carrier_components={row[4]}/danger_components={row[5]}"
            for row in literal_controls
        ),
        f"terminal=residual_bodies:{terminal_summary[0]};two_high_positive:{terminal_summary[1]};closed:{terminal_summary[2]};zero_high:{terminal_summary[3]};one_high:{terminal_summary[4]};cardinality:{terminal_summary[5]};maxgap:{terminal_summary[6]};failed:{terminal_summary[7]};unit_checks:{terminal_summary[8]};coarse_cardinality:{sum(row[10] for row in terminals)};exact_cardinality:{sum(row[11] for row in terminals)};terminal_record_sha256:{terminal_record_sha}",
        f"minimum_positive_gap={ftext(minimum_positive_gap)};body={minimum_row[1]}",
    ]
    for row in records:
        lines.append(
            f"BODY;E={row[1]};L={row[2]};high={row[3]};wall={row[5]};"
            f"states={row[9]};crude={row[10]};status={row[11]};residual={row[12]};"
            f"direct={row[19]};legacy={row[20]};stage_sha256={row[18]};"
            f"residual_sha256={hashlib.sha256(repr(row[13]).encode()).hexdigest()};audit_sha256={row[21]}"
        )
    for row in terminals:
        lines.append(
            f"TERMINAL;E={row[1]};residual={row[4]};gap={ftext(row[5])};"
            f"two_high={row[19]};zero_high={row[7]};one_high={row[8]};"
            f"coarse={row[10]};exact={row[11]};maxgap={row[12]};"
            f"cardinality={row[13]};failed={row[14]};slack={row[16]};"
            f"closed={row[20]};case_sha256={row[17]}"
        )
    for row in survivors:
        lines.append(
            f"SURVIVOR;E={row[1]};gap={ftext(row[5])};witness={row[6]};"
            f"one_high={row[8]};cardinality={row[13]};"
            f"literal_maximizer_safe_mass={ftext(EXPECTED_LITERAL_SAFE_MASS[row[1]])};"
            f"case_sha256={row[17]}"
        )
    lines += [
        "first_failed_implication=duplicate_two_high_gap > 0",
        "strongest_survivor=all 1,245,093 one-high terminal cases close by exact complete-cell cardinality; only four two-high scalar gaps are nonpositive",
        "repair=accept only the generator's exact first zero-good one-row Farkas certificate locally; every other certificate is delegated unchanged to the legacy full-table auditor",
        "memory_safety=each fixed low-label pair is vectorized and discarded; if C clean cells lie over d residues, support is at least ceil(C/(L/d)); exact residues and inherited max-gap are used only when this lower bound is insufficient",
        "two_high_boundary=the four scalar maximizing literal packets all retain positive exact interval mass on ruler lcm(5045040,14*labels); low-cell cardinality cannot be reused for high labels and no physical cover is asserted",
        "scope=projected necessary sector only;candidate boundary;no ledger decrement, projected-cap change, navigation mutation, or LRC consequence",
        f"semantic_sha256={semantic}",
        "all_exact_controls=PASS",
    ]
    payload = "\n".join(lines) + "\n"
    args.output.write_text(payload, encoding="utf-8", newline="\n")
    print(payload, end="")


if __name__ == "__main__":
    main()
