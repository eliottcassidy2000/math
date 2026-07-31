#!/usr/bin/env python3
"""Uniformly close all fifteen projected k=2 rows at z1=1736.

This is a self-contained hybrid replay against the current canonical
THM-2941 uniform-ray and projected-cell engines.  It does not import the
older k=2 descent wrappers, whose guards still pin a superseded byte hash.

The exact all-body atlas partitions the height into nine ordinary rows and
six explicit HIGH-TAIL rows.  The ordinary rows use

  ray maxima -> denominator states -> crude capacities -> common K5 status
  -> finite scalar-slack packets -> lossless projected residual.

Five HIGH rows have a negative exact forced-high scalar maximum.  The sixth
has one positive scalar packet; a positive all-label cutoff makes its
universe finite and its projected residual closes it.  The three routes are
disjoint and exhaust the atlas height, uniformly over all distinct later
nonaligned labels and without a label horizon.
"""

from __future__ import annotations

import argparse
import multiprocessing as mp
import os
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, combinations_with_replacement, product
from math import comb, gcd, lcm
from pathlib import Path

import numpy as np
from scipy.optimize import linprog


ROOT = Path(__file__).resolve().parents[1]
UNIFORM = ROOT / "04-computation/lrc14_j7_k3_uniform_ray_status_closure_thm2941.py"
PROJECTED = ROOT / "04-computation/lrc14_j7_five_aligned_two_drift_projected_closure_thm2941.py"
ATLAS_SOURCE = ROOT / "04-computation/lrc14_j7_k2_scalar_band_1680_1742_thm2941.py"
ATLAS_OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_k2_scalar_band_1680_1742_thm2941.out"
CAP_SOURCE = ROOT / "04-computation/lrc14_j7_k2_cap1736_scalar_splice_thm2941.py"
CAP_OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_k2_cap1736_scalar_splice_thm2941.out"
OUTPUT_PATH = ROOT / "05-knowledge/results/lrc14_j7_k2_z1736_hybrid_closure_thm2941.out"

EXPECTED_HASHES = (
    "34ab29162ed33d90093e6d2bf781def36c420a1cd6596158b5d6579a3a8f3f46",
    "76f891edfcc029a08202481304a809e03e8bd81f247afaeabab685825c4d3662",
    "89016f939c961fa979ec5b30812981456df5bfb2af3066f1f1b38e5a83f1a412",
    "4a36611b26585964e185bbaa3d583be3f1c67a7b608cca785920266bc217a779",
    "89c09d1199c94cdd95b79589aebc7b3f9a3ee3ef1ef63dd8a2c52e6a19ea22d7",
    "084ba0464f16b45c05ccf90e828d9be9d5cb55d1181e8dd60e6ac663d830827c",
)
EXPECTED_PROFILE_SHA256 = (
    "5ff5a8584cdd29f4dbc042a1436aa9b4e8f3b15107d475f2bf365cb8c4b12fd4"
)
EXPECTED_SEMANTIC_SHA256 = (
    "4eb8cf508b9b2b80f7f69cf411117740a7ad3ed93527f41dd1a565077f290e0a"
)

QUANTIFIER = "distinct later nonaligned labels"
FIRST = 1736
SUFFIX_SLOTS = 4
HIGH_WALL_RATIO = F(13, 150)
SCALAR_ETA = F(1, 91)
ALIGNED_TWO_UNION_CAP = F(25, 91)

ORDINARY_CASES = tuple(
    (FIRST, body)
    for body in (
        (1, 2, 8, 10, 12, 14),
        (1, 4, 8, 10, 12, 14),
        (1, 6, 8, 10, 12, 14),
        (2, 3, 8, 10, 12, 14),
        (2, 4, 6, 8, 10, 14),
        (2, 4, 8, 10, 12, 14),
        (2, 6, 8, 10, 12, 14),
        (3, 4, 8, 10, 12, 14),
        (4, 6, 8, 10, 12, 14),
    )
)
HIGH_NEGATIVE_CASES = tuple(
    (FIRST, body)
    for body in (
        (1, 2, 10, 11, 12, 14),
        (1, 2, 10, 12, 13, 14),
        (1, 4, 10, 12, 13, 14),
        (2, 3, 10, 11, 12, 14),
        (2, 8, 9, 10, 12, 14),
    )
)
EXCEPTIONAL_CASE = (FIRST, (1, 8, 10, 12, 13, 14))
ALL_CASES = tuple(sorted((*ORDINARY_CASES, *HIGH_NEGATIVE_CASES, EXCEPTIONAL_CASE)))

EXPECTED_HEIGHTS = (
    (1683, 1),
    (1694, 10),
    (1702, 3),
    (1708, 14),
    (1722, 11),
    (1724, 2),
    (1732, 2),
    (1736, 15),
)
# scalar, crude kills, common-status kills, status survivors, literal packets
EXPECTED_COUNTS = {
    ORDINARY_CASES[0]: (46, 0, 16, 30, 43),
    ORDINARY_CASES[1]: (617, 0, 611, 6, 8),
    ORDINARY_CASES[2]: (30, 0, 23, 7, 10),
    ORDINARY_CASES[3]: (3, 0, 1, 2, 2),
    ORDINARY_CASES[4]: (1, 0, 1, 0, 0),
    ORDINARY_CASES[5]: (198, 0, 109, 89, 194),
    ORDINARY_CASES[6]: (30, 0, 9, 21, 29),
    ORDINARY_CASES[7]: (3, 0, 3, 0, 0),
    ORDINARY_CASES[8]: (1, 0, 1, 0, 0),
}
EMPTY_SHA256 = "2e38e77b22c314a449e91fafed92a43826ac6aa403ae6a8acb6cf58239fbaf5d"
EXPECTED_EXACT_AUDIT = {
    ORDINARY_CASES[0]: (
        (
            "d7f2e2b289c739bb667af041a6617fc915ce88f0ddd47d6e56837057cea99cb7",
            EMPTY_SHA256,
            "9e7a9bb79b7261393314d6ed223c4374c7ed6aa9226aa036f80e140f0fe46732",
            "ac04d1fbda030f66542a66e06989e6abbe92b2593534ef3fc1bd4b0d465a1131",
        ),
        F(23897, 1434069), 4,
        "77491cb5de95588fee59c489e0f8427b9365799aa40d02350f8c68765e35e5e3",
    ),
    ORDINARY_CASES[1]: (
        (
            "45d8ffc36a4488ae7e95a7dffd7f7f5840c4d00c198696c11f5fbc8f72b58ebe",
            EMPTY_SHA256,
            "31efc3c1d5482d4297d9bcd92a0728c6b0363455a20492a4b9809e0a5ffc326b",
            "97997123d9471015e1971c89dbf2a4e74eaf2c0e21af90603532affe64207b6e",
        ),
        F(161107, 6585579), 3,
        "9cf36665aa9b9aa3b5d1845bb8abc3d8b23d1c917610ba951cb2c1cd651cc490",
    ),
    ORDINARY_CASES[2]: (
        (
            "60fbbf3a845c112400e2d50b49879b4edfb7dd09ffccc9d6444daccdee212a3b",
            EMPTY_SHA256,
            "b696cbb6f897347c675a8945c079958f2365ff010c40dbb4f389cd8de5a05902",
            "c6a2125986b0135103f443cd8b8cd5b3873d15477a735baad86b02e4ec02781a",
        ),
        F(590, 2821), 2,
        "ae15f835044caddbfe3a81d637ce360b8d7d3a331c7e9906d2bbbbf86405a2d1",
    ),
    ORDINARY_CASES[3]: (
        (
            "cebb57972815261050413c464f8e1a83754632091322901fcc57f738d3551692",
            EMPTY_SHA256,
            "cfd55da8473cec11821f112d0628d75c72a2ead8cd1c0b93bd2a7d8162964a05",
            "650d42500026f5fb381e6a6fd204a0555b014dabd4b094e90e434f145e7e94df",
        ),
        F(681, 2821), 1,
        "bb6d1f3147315ebda37641915eb8a385283bfad504e12c5e6fd4e27162958439",
    ),
    ORDINARY_CASES[4]: (
        (
            "d6891c20082b2dd5793b175241f569fc3227a8b44e7626b1ae136018c3707d8a",
            EMPTY_SHA256,
            "34a53d625d8af279bc28ca392bdb1e9791c417918b815e62c2753fa2bc05561d",
            EMPTY_SHA256,
        ),
        None, 0, EMPTY_SHA256,
    ),
    ORDINARY_CASES[5]: (
        (
            "31dfd0f4e4ed72f597d4906f4267d649bd5095833b23ead5352175dbc5b4bf9f",
            EMPTY_SHA256,
            "1a745151627d0e4944a283ee72c54d77de1abefad03f0ab6d4333bc0bfb6fada",
            "7b2555bbbd300cd6b33fa4505587bfeb2991837a2e10625575acd9a7d48bd454",
        ),
        F(681, 2821), 2,
        "d3484bf4398240d2a64b78dd5f4075037bfde530a599d8f9fea77f9123140748",
    ),
    ORDINARY_CASES[6]: (
        (
            "1261a50f3b553e4b6b96d573fb5c3207ec5679687d47247db37dbea669d8c825",
            EMPTY_SHA256,
            "e95586b6978c432a610350d11ad411148c8f01c9369e3810bf9214b8e90ce8c2",
            "de7635dde176a1f80604e781a7102d5869880883e668773c78c1b4511821233b",
        ),
        F(49915, 510601), 2,
        "34c8a735ae7198acdef49ff96ceff9b151226e97cf0880c98abbba038f4b6e1c",
    ),
    ORDINARY_CASES[7]: (
        (
            "c659774c4dbcdce1977a3dd711a88e955c7d10772d170bca6485ca8039de9a7d",
            EMPTY_SHA256,
            "3bd88628d5ce3fd6bf0b0bb0ab5af1aba403fe51da1c559bf8cdb66e4a2ede1f",
            EMPTY_SHA256,
        ),
        None, 0, EMPTY_SHA256,
    ),
    ORDINARY_CASES[8]: (
        (
            "e8f6ee0c0fdfe775f96b8041771f97cfc7d2d615177d10e9f62adbe99d367f1c",
            EMPTY_SHA256,
            "156856e01f7110d608b617993336df659a9bfe64e83ee7118a2c5eb50299e162",
            EMPTY_SHA256,
        ),
        None, 0, EMPTY_SHA256,
    ),
}
EXPECTED_HIGH_GAPS = {
    HIGH_NEGATIVE_CASES[0]: F(-17704278379, 151099294846500),
    HIGH_NEGATIVE_CASES[1]: F(-2122333151, 17733467524500),
    HIGH_NEGATIVE_CASES[2]: F(-7353136579, 118459563063660),
    HIGH_NEGATIVE_CASES[3]: F(-47509753883, 178751453380500),
    HIGH_NEGATIVE_CASES[4]: F(-7891025249, 164964861300785),
}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def file_sha256(path):
    return sha256(path.read_bytes()).hexdigest()


def load_module(name, path):
    spec = spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, f"cannot load {path}")
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def ftext(value):
    if value is None:
        return "NONE"
    return f"{value.numerator}/{value.denominator}"


PATHS = (UNIFORM, PROJECTED, ATLAS_SOURCE, ATLAS_OUTPUT, CAP_SOURCE, CAP_OUTPUT)
require(tuple(file_sha256(path) for path in PATHS) == EXPECTED_HASHES, "dependency changed")
U = load_module("k2_z1736_uniform", UNIFORM)
P = load_module("k2_z1736_projected", PROJECTED)
require(P.A.RULER == U.suffix.A.RULER, "master rulers disagree")


def atlas_partition():
    keys = []
    high = []
    heights = Counter()
    for line in ATLAS_OUTPUT.read_text(encoding="utf-8").splitlines():
        if not line.startswith("SURVIVOR;"):
            continue
        fields = dict(field.split("=", 1) for field in line.split(";")[1:] if "=" in field)
        key = (int(fields["z1"]), tuple(map(int, fields["E"].split(","))))
        keys.append(key)
        heights[key[0]] += 1
        if "HIGH-TAIL" in line:
            high.append(key)
    return tuple(sorted(keys)), tuple(sorted(high)), tuple(sorted(heights.items()))


ATLAS_KEYS, ATLAS_HIGH, ATLAS_HEIGHTS = atlas_partition()
require(ATLAS_HEIGHTS == EXPECTED_HEIGHTS, "atlas height profile changed")
require(tuple(key for key in ATLAS_KEYS if key[0] == FIRST) == ALL_CASES, "z1736 cases changed")
require(
    tuple(key for key in ATLAS_HIGH if key[0] == FIRST)
    == tuple(sorted((*HIGH_NEGATIVE_CASES, EXCEPTIONAL_CASE))),
    "z1736 HIGH partition changed",
)
require(
    not any(key[0] in (*range(1733, 1736), *range(1737, 1743)) for key in ATLAS_KEYS),
    "intervening height occupied",
)
require(sum(1 for key in ATLAS_KEYS if key[0] == 1732) == 2, "next height changed")
require("consequence=projected k=2 first drift label z1<=1736" in CAP_OUTPUT.read_text(), "cap input changed")


def spanning_trees(vertex_count):
    edges = tuple(combinations(range(vertex_count), 2))
    result = []
    for chosen in combinations(edges, vertex_count - 1):
        seen = {0}
        changed = True
        while changed:
            changed = False
            for left, right in chosen:
                if left in seen and right not in seen:
                    seen.add(right)
                    changed = True
                elif right in seen and left not in seen:
                    seen.add(left)
                    changed = True
        if len(seen) == vertex_count:
            result.append(tuple(sorted(chosen)))
    require(len(result) == vertex_count ** (vertex_count - 2), "Cayley count changed")
    return tuple(sorted(result))


def prufer_trees(vertex_count):
    result = set()
    for word in product(range(vertex_count), repeat=vertex_count - 2):
        degree = [1] * vertex_count
        for vertex in word:
            degree[vertex] += 1
        edges = []
        for vertex in word:
            leaf = next(index for index, value in enumerate(degree) if value == 1)
            edges.append(tuple(sorted((leaf, vertex))))
            degree[leaf] -= 1
            degree[vertex] -= 1
        leaves = tuple(index for index, value in enumerate(degree) if value == 1)
        require(len(leaves) == 2, "Pruefer terminal changed")
        edges.append(tuple(sorted(leaves)))
        result.add(tuple(sorted(edges)))
    return tuple(sorted(result))


TREES5 = spanning_trees(5)
require(TREES5 == prufer_trees(5), "independent K5 tree atlases disagree")


@lru_cache(maxsize=None)
def hunter_cap5(M, ds, es):
    sizes = tuple((M // d) * e for d, e in zip(ds, es))
    overlaps = {
        (left, right): U.pair_lower(M, ds[left], es[left], ds[right], es[right])
        for left, right in combinations(range(5), 2)
    }
    invoice = max(sum(overlaps[edge] for edge in tree) for tree in TREES5)
    return min(M, sum(sizes) - invoice)


def hunter_status_data5(D, ds, q):
    M = D // q
    inner = []
    lows = []
    marginals = []
    for d in ds:
        common = gcd(d, q)
        low, remainder = divmod((d + 6) // 7, common)
        inner.append(d // common)
        lows.append(low)
        marginals.append((q // common) * remainder)
    capacities = tuple(
        hunter_cap5(
            M,
            tuple(inner),
            tuple(lows[index] + ((pattern >> index) & 1) for index in range(5)),
        )
        for pattern in range(32)
    )
    return tuple(marginals), capacities


@lru_cache(maxsize=None)
def common_status_feasible5(q, marginals, capacities, histogram):
    tail_rows = []
    tail_rhs = []
    thresholds = []
    for threshold, _count in histogram:
        if threshold <= 0:
            continue
        demand = sum(count for load, count in histogram if load >= threshold)
        good = tuple(int(capacity >= threshold) for capacity in capacities)
        if all(good):
            continue
        thresholds.append(threshold)
        if not any(good):
            return False, ((threshold,), (F(1),), (F(0),) * 6)
        tail_rows.append(good)
        tail_rhs.append(demand)
    equality_rows = [
        (1,) * 32,
        *[tuple((pattern >> index) & 1 for pattern in range(32)) for index in range(5)],
    ]
    equality_rhs = (q, *marginals)
    if not tail_rows:
        return True, None
    primal = linprog(
        np.zeros(32),
        A_ub=-np.asarray(tail_rows, dtype=float),
        b_ub=-np.asarray(tail_rhs, dtype=float),
        A_eq=np.asarray(equality_rows, dtype=float),
        b_eq=np.asarray(equality_rhs, dtype=float),
        bounds=(0, None),
        method="highs",
    )
    if primal.success:
        return True, None
    require(primal.status == 2, ("unexpected primal status", primal.status))
    tail_count = len(tail_rows)
    dual_rows = []
    for pattern in range(32):
        dual_rows.append(
            tuple(tail_rows[row][pattern] for row in range(tail_count))
            + tuple(-equality_rows[row][pattern] for row in range(6))
        )
    dual_rows.append(tuple(-value for value in tail_rhs) + tuple(equality_rhs))
    dual = linprog(
        np.zeros(tail_count + 6),
        A_ub=np.asarray(dual_rows, dtype=float),
        b_ub=np.asarray((0,) * 32 + (-1,), dtype=float),
        bounds=[(0, None)] * tail_count + [(None, None)] * 6,
        method="highs",
    )
    if not dual.success:
        return True, None
    alpha = tuple(F(float(value)).limit_denominator(1_000_000) for value in dual.x[:tail_count])
    z = tuple(F(float(value)).limit_denominator(1_000_000) for value in dual.x[tail_count:])
    slacks = tuple(
        sum(z[row] * equality_rows[row][pattern] for row in range(6))
        - sum(alpha[row] * tail_rows[row][pattern] for row in range(tail_count))
        for pattern in range(32)
    )
    contradiction = sum(z[row] * equality_rhs[row] for row in range(6)) - sum(
        alpha[row] * tail_rhs[row] for row in range(tail_count)
    )
    if all(value >= 0 for value in alpha) and all(value >= 0 for value in slacks) and contradiction < 0:
        return False, (tuple(thresholds), alpha, z)
    return True, None


def projected_safe_lower_bound(cells, L, labels, stop_at_cap=True):
    common_danger = ((F(0), F(1)),)
    cells_used = 0
    for cell in cells:
        local_union = P.merge_fraction(
            [interval for label in labels for interval in P.phase_danger(cell, label, L)]
        )
        common_danger = P.intersect_fraction(common_danger, local_union)
        cells_used += 1
        if stop_at_cap and P.interval_mass(common_danger) <= 1 - ALIGNED_TWO_UNION_CAP:
            break
    return 1 - P.interval_mass(common_danger), cells_used, common_danger


def direct_projected_mass(body, L, labels):
    carrier = tuple((F(left, P.A.RULER), F(right, P.A.RULER)) for left, right in P.A.carrier_for(body))
    removed = P.merge_fraction(
        [interval for label in labels for interval in P.danger_fraction(label)]
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


def ray_and_status(first, body, carrier, h, lower, L):
    first_delta = U.delta(carrier, h, first)
    first_d = L // gcd(L, first)
    amplitudes = [F(0)]
    signs = Counter()
    for residue in range(1, L):
        amplitude = residue * U.delta(carrier, h, residue)
        require((residue + L) * U.delta(carrier, h, residue + L) == amplitude, (body, "ray"))
        amplitudes.append(amplitude)
        signs[(amplitude > 0) - (amplitude < 0)] += 1
    require(all(amplitudes[L - r] == -amplitudes[r] for r in range(1, L)), (body, "antipode"))
    require(L * U.delta(carrier, h, L) == 0, (body, "zero ray"))
    require(signs[1] == signs[-1] and sum(signs.values()) == L - 1, (body, "signs"))
    ray_digest = sha256(repr(tuple(amplitudes)).encode()).hexdigest()
    divisors = tuple(d for d in U.support.divisors(L) if d > 1)

    @lru_cache(maxsize=None)
    def top_values(d, multiplicity):
        candidates = []
        for direction in range(1, d):
            if gcd(direction, d) != 1:
                continue
            residue = (L // d) * direction
            amplitude = amplitudes[residue]
            if amplitude < 0:
                continue
            label = residue
            if label <= first:
                label += ((first + 1 - label + L - 1) // L) * L
            candidates.extend(
                (amplitude / (label + offset * L), label + offset * L, residue)
                for offset in range(multiplicity)
            )
        candidates.sort(key=lambda row: (-row[0], row[1], row[2]))
        require(len(candidates) >= multiplicity, (body, "missing ray", d))
        return tuple(candidates[:multiplicity])

    trials = comb(len(divisors) + SUFFIX_SLOTS - 1, SUFFIX_SLOTS)
    scalar = []
    for tail in combinations_with_replacement(divisors, SUFFIX_SLOTS):
        upper = first_delta
        labels = []
        for d, multiplicity in Counter(tail).items():
            chosen = top_values(d, multiplicity)
            upper += sum((value for value, _label, _residue in chosen), F())
            labels.extend(label for _value, label, _residue in chosen)
        if upper >= lower:
            scalar.append((tuple(sorted((first_d, *tail))), upper, (first, *sorted(labels))))
    require(len({row[0] for row in scalar}) == len(scalar), (body, "duplicate state"))

    actual_L, ranges = U.support.safe_cell_ranges(body)
    require(actual_L == L, (body, "safe-cell ruler"))
    arcs_cache = {}
    crude_kills = []
    crude_survivors = []
    for ds, upper, labels in scalar:
        D = lcm(*ds)
        if D not in arcs_cache:
            arcs_cache[D] = U.fibre.projected_support_arcs(D, ranges)
        arcs = arcs_cache[D]
        witness = None
        for q in U.support.divisors(D):
            histogram = U.fibre.residue_load_histogram(arcs, q)
            target = max(load for load, count in histogram if count)
            capacity = sum(U.fibre_cap(D, d, q) for d in ds)
            if target > capacity:
                witness = (q, D // q, target, capacity)
                break
        row = (ds, upper, labels, witness)
        (crude_survivors if witness is None else crude_kills).append(row)

    status_kills = []
    status_survivors = []
    for ds, upper, labels, _witness in crude_survivors:
        D = lcm(*ds)
        arcs = arcs_cache[D]
        witness = None
        for M in U.support.divisors(D):
            q = D // M
            marginals, capacities = hunter_status_data5(D, ds, q)
            histogram = U.fibre.residue_load_histogram(arcs, q)
            feasible, certificate = common_status_feasible5(q, marginals, capacities, histogram)
            if not feasible:
                require(certificate is not None, (body, ds, "missing Farkas witness"))
                witness = (q, M, marginals, tuple(sorted(set(capacities))), histogram, certificate)
                break
        row = (ds, upper, labels, witness)
        (status_survivors if witness is None else status_kills).append(row)
    stages = (scalar, crude_kills, status_kills, status_survivors)
    stage_digests = tuple(sha256(repr(tuple(rows)).encode()).hexdigest() for rows in stages)
    return (
        tuple(amplitudes), ray_digest, len(divisors), trials, first_delta, first_d,
        tuple(scalar), tuple(crude_kills), tuple(status_kills), tuple(status_survivors), stage_digests,
    )


def projected_packets(first, body, carrier, h, lower, L, amplitudes, first_delta, states):
    cells = P.body_cells(P.A.carrier_for(body), L)

    @lru_cache(maxsize=None)
    def ray_heads(d, multiplicity):
        candidates = []
        for direction in range(1, d):
            if gcd(direction, d) != 1:
                continue
            residue = (L // d) * direction
            amplitude = amplitudes[residue]
            if amplitude < 0:
                continue
            label = residue
            if label <= first:
                label += ((first + 1 - label + L - 1) // L) * L
            candidates.extend(
                (amplitude / (label + offset * L), label + offset * L, residue)
                for offset in range(multiplicity)
            )
        candidates.sort(key=lambda row: (-row[0], row[1], row[2]))
        return tuple(candidates[:multiplicity])

    @lru_cache(maxsize=None)
    def eligible_group(d, multiplicity, slack):
        heads = ray_heads(d, multiplicity)
        top_sum = sum((row[0] for row in heads), F())
        threshold = heads[-1][0] - slack
        require(threshold > 0, (body, "nonfinite ray cutoff", d, multiplicity))
        candidates = []
        for direction in range(1, d):
            if gcd(direction, d) != 1:
                continue
            residue = (L // d) * direction
            amplitude = amplitudes[residue]
            if amplitude <= 0:
                continue
            label = residue
            if label <= first:
                label += ((first + 1 - label + L - 1) // L) * L
            while amplitude / label >= threshold:
                candidates.append((amplitude / label, label, residue))
                label += L
        candidates.sort(key=lambda row: (-row[0], row[1], row[2]))
        choices = []
        for chosen in combinations(candidates, multiplicity):
            labels = tuple(sorted(row[1] for row in chosen))
            if len(set(labels)) != multiplicity:
                continue
            value = sum((row[0] for row in chosen), F())
            deficit = top_sum - value
            if deficit <= slack:
                choices.append((deficit, labels, value))
        choices.sort()
        require(choices and choices[0][0] == 0, (body, "maximizer missing", d))
        return top_sum, threshold, tuple(candidates), tuple(choices)

    packet_count = 0
    killed_count = 0
    survivors = []
    minimum_margin = None
    maximum_cells = 0
    minimum_row = None
    state_audit = []
    first_d = L // gcd(L, first)
    for state_index, (ds, upper, maximizing_labels, witness) in enumerate(states, 1):
        require(witness is None, (body, "status witness on survivor"))
        slack = upper - lower
        tail = list(ds)
        tail.remove(first_d)
        groups = []
        for d, multiplicity in sorted(Counter(tail).items()):
            top_sum, threshold, candidates, choices = eligible_group(d, multiplicity, slack)
            groups.append((d, multiplicity, top_sum, threshold, len(candidates), choices))
        require(first_delta + sum((group[2] for group in groups), F()) == upper, (body, "upper"))
        state_packets = 0
        state_kills = 0
        state_minimum = None
        for selection in product(*(group[5] for group in groups)):
            deficit = sum((choice[0] for choice in selection), F())
            if deficit > slack:
                continue
            labels = (first, *sorted(label for choice in selection for label in choice[1]))
            require(len(labels) == 5 and len(set(labels)) == 5, (body, "literal labels"))
            actual_upper = first_delta + sum((choice[2] for choice in selection), F())
            require(actual_upper == upper - deficit and actual_upper >= lower, (body, "scalar"))
            projected, cells_used, _common = projected_safe_lower_bound(cells, L, labels)
            margin = projected - ALIGNED_TWO_UNION_CAP
            packet_count += 1
            state_packets += 1
            maximum_cells = max(maximum_cells, cells_used)
            minimum_margin = margin if minimum_margin is None else min(minimum_margin, margin)
            state_minimum = margin if state_minimum is None else min(state_minimum, margin)
            if minimum_row is None or (margin, labels) < (minimum_row[0], minimum_row[1]):
                minimum_row = (margin, labels, projected, cells_used)
            if margin > 0:
                killed_count += 1
                state_kills += 1
            else:
                survivors.append((state_index, ds, labels, actual_upper, projected))
        require(state_packets > 0, (body, "empty literal state"))
        state_audit.append(
            (
                state_index, ds, upper, maximizing_labels, slack,
                tuple((d, m, top, threshold, count, len(choices)) for d, m, top, threshold, count, choices in groups),
                state_packets, state_kills, state_minimum,
            )
        )
    require(packet_count == killed_count and not survivors, (body, "projected survivor", survivors))
    direct = None
    if minimum_row is not None:
        _margin, labels, prefix_mass, _cells = minimum_row
        full_mass, full_cells, _common = projected_safe_lower_bound(cells, L, labels, False)
        direct = direct_projected_mass(body, L, labels)
        require(full_cells == len(cells) and direct == full_mass and direct >= prefix_mass, (body, "direct"))
    return (
        len(cells), packet_count, killed_count, minimum_margin, maximum_cells,
        minimum_row, direct, sha256(repr(tuple(state_audit)).encode()).hexdigest(), tuple(state_audit),
    )


def ordinary_profile(case):
    first, body = case
    carrier = U.suffix.A.carrier_for(body)
    require(P.A.carrier_for(body) == carrier, (body, "carrier disagreement"))
    h = F(sum(right - left for left, right in carrier), U.suffix.A.RULER)
    lower = h * SCALAR_ETA
    L = 14 * lcm(*body)
    require(U.suffix.A.RULER % L == 0 and first >= HIGH_WALL_RATIO * L, (body, "ruler/wall"))
    data = ray_and_status(first, body, carrier, h, lower, L)
    amplitudes, ray_digest, divisor_count, trials, first_delta, first_d, scalar, crude, status, states, stage = data
    projected = projected_packets(first, body, carrier, h, lower, L, amplitudes, first_delta, states) if states else (
        len(P.body_cells(carrier, L)), 0, 0, None, 0, None, None, sha256(b"()").hexdigest(), (),
    )
    counts = (len(scalar), len(crude), len(status), len(states), projected[1])
    require(counts == EXPECTED_COUNTS[case], (case, "counts", counts))
    require(projected[1] == projected[2], (case, "packet survived"))
    audit = (stage, projected[3], projected[4], projected[7])
    if EXPECTED_EXACT_AUDIT[case] is not None:
        require(audit == EXPECTED_EXACT_AUDIT[case], (case, "audit", audit))
    return (
        "ORDINARY", first, body, h, len(carrier), L, lower, first_delta, first_d,
        ray_digest, divisor_count, trials, counts, stage, *projected[:-1],
    )


def forced_high_profile(case):
    first, body = case
    carrier = U.suffix.A.carrier_for(body)
    h = F(sum(right - left for left, right in carrier), U.suffix.A.RULER)
    lower = h * SCALAR_ETA
    L = 14 * lcm(*body)
    high_floor = max(15, (HIGH_WALL_RATIO * L).numerator // (HIGH_WALL_RATIO * L).denominator + 1)
    require(first < high_floor, (case, "not HIGH"))
    amplitudes = [F(0)]
    signs = Counter()
    for residue in range(1, L):
        amplitude = residue * U.delta(carrier, h, residue)
        require((residue + L) * U.delta(carrier, h, residue + L) == amplitude, (case, "ray"))
        amplitudes.append(amplitude)
        signs[(amplitude > 0) - (amplitude < 0)] += 1
    require(all(amplitudes[L-r] == -amplitudes[r] for r in range(1, L)), (case, "antipode"))
    arbitrary = []
    high = []
    omitted = []
    for residue in range(1, L):
        amplitude = amplitudes[residue]
        if amplitude <= 0:
            continue
        first_label = residue
        if first_label <= first:
            first_label += ((first + 1 - first_label + L - 1) // L) * L
        for offset in range(SUFFIX_SLOTS):
            label = first_label + offset * L
            arbitrary.append((amplitude / label, label, residue, offset))
        fifth = first_label + SUFFIX_SLOTS * L
        omitted.append((amplitude / fifth, fifth, residue))
        high_label = residue
        if high_label < high_floor:
            high_label += ((high_floor - high_label + L - 1) // L) * L
        high.append((amplitude / high_label, high_label, residue))
    rank4 = tuple(sorted(arbitrary, key=lambda row: (-row[0], row[1:]))[:4])
    omitted_max = min(omitted, key=lambda row: (-row[0], row[1:]))
    require(omitted_max[0] <= rank4[-1][0], (case, "ray truncation"))
    best_high = min(high, key=lambda row: (-row[0], row[1:]))
    if any(label >= high_floor for _value, label, _residue, _offset in rank4):
        constrained = rank4
        branch = "UNRESTRICTED_TOP4_ALREADY_HIGH"
    else:
        constrained = (*rank4[:3], (best_high[0], best_high[1], best_high[2], 0))
        branch = "GLOBAL_TOP3_PLUS_BEST_HIGH"
    require(len({row[1] for row in constrained}) == 4 and any(row[1] >= high_floor for row in constrained), (case, "constraint"))
    first_delta = U.delta(carrier, h, first)
    upper = first_delta + sum((row[0] for row in constrained), F())
    gap = upper - lower
    require(gap == EXPECTED_HIGH_GAPS[case] < 0, (case, "gap", gap))
    return (
        "HIGH-NEGATIVE", first, body, h, len(carrier), L, high_floor, first_delta, lower,
        signs[1], signs[-1], signs[0], sha256(repr(tuple(amplitudes)).encode()).hexdigest(),
        rank4, omitted_max, best_high, branch, tuple(constrained), upper, gap,
    )


def exceptional_profile(case):
    first, body = case
    carrier = U.suffix.A.carrier_for(body)
    h = F(sum(right - left for left, right in carrier), U.suffix.A.RULER)
    lower = h * SCALAR_ETA
    L = 14 * lcm(*body)
    high_floor = max(15, (HIGH_WALL_RATIO * L).numerator // (HIGH_WALL_RATIO * L).denominator + 1)
    require((h, len(carrier), L, high_floor) == (F(811,2548), 38, 152880, 13250), (case, "body"))
    amplitudes = [F(0)]
    signs = Counter()
    ray_heads = []
    for residue in range(1, L):
        amplitude = residue * U.delta(carrier, h, residue)
        require((residue + L) * U.delta(carrier, h, residue + L) == amplitude, (case, "ray"))
        amplitudes.append(amplitude)
        signs[(amplitude > 0) - (amplitude < 0)] += 1
        if amplitude > 0:
            label = residue
            if label <= first:
                label += ((first + 1 - label + L - 1) // L) * L
            ray_heads.append((amplitude / label, label, residue, amplitude))
    require(all(amplitudes[L-r] == -amplitudes[r] for r in range(1,L)), (case, "antipode"))
    ray_heads.sort(key=lambda row: (-row[0], row[1:]))
    first_delta = U.delta(carrier, h, first)
    target = lower - first_delta
    top3 = tuple(ray_heads[:3])
    threshold = target - sum((row[0] for row in top3), F())
    require(threshold == F(4876247,30609706218) > 0, (case, "cutoff", threshold))
    candidates = []
    for _value, label, residue, amplitude in ray_heads:
        while amplitude / label >= threshold:
            candidates.append((amplitude / label, label, residue, amplitude))
            label += L
    candidates.sort(key=lambda row: (-row[0], row[1:]))
    require(len(candidates) == 501 and len({row[2] for row in candidates}) == 501, (case, "candidates"))
    high = tuple(row for row in candidates if row[1] >= high_floor)
    require(high == ((F(149,909636),13260,13260,F(745,343)),), (case, "high"))
    forced = high[0]
    others = tuple(row for row in candidates if row != forced)
    values = tuple(row[0] for row in others)
    packets = []
    nonempty_pairs = 0
    for left in range(len(others)-2):
        for middle in range(left+1,len(others)-1):
            needed = target - forced[0] - values[left] - values[middle]
            low = middle + 1
            upper_index = len(others)
            while low < upper_index:
                pivot = (low + upper_index) // 2
                if values[pivot] >= needed:
                    low = pivot + 1
                else:
                    upper_index = pivot
            if low > middle + 1:
                nonempty_pairs += 1
            for right in range(middle+1,low):
                chosen = (forced, others[left], others[middle], others[right])
                labels = tuple(sorted(row[1] for row in chosen))
                gap = sum((row[0] for row in chosen), F()) - target
                require(gap >= 0 and len(set(labels)) == 4, (case, "packet"))
                packets.append((gap, labels, chosen))
    packets.sort()
    require(nonempty_pairs == 1 and len(packets) == 1, (case, "packet count"))
    gap, suffix_labels, chosen = packets[0]
    require((gap,suffix_labels) == (F(91785,20406470812),(1836,2004,2340,13260)), (case, "unique"))
    labels = (first,*suffix_labels)
    cells = P.body_cells(carrier,L)
    hostile, hostile_cells, hostile_common = projected_prefix(cells,L,labels,13)
    projected, cells_used, common = projected_safe_lower_bound(cells,L,labels)
    margin = projected - ALIGNED_TWO_UNION_CAP
    direct = direct_projected_mass(body,L,labels)
    require((hostile,hostile_cells,hostile_common) == (F(0),13,((F(0),F(1)),)), (case,"hostile"))
    require((projected,cells_used,common,margin,direct) == (F(14,17),14,((F(0),F(3,17)),),F(849,1547),F(1)), (case,"projection"))
    return (
        "HIGH-PACKET", first, body, h, len(carrier), L, high_floor, first_delta, lower,
        target, tuple(signs[index] for index in (1,-1,0)), sha256(repr(tuple(amplitudes)).encode()).hexdigest(),
        top3, threshold, len(candidates), sha256(repr(tuple(candidates)).encode()).hexdigest(), high,
        nonempty_pairs, len(packets), sha256(repr(tuple(packets)).encode()).hexdigest(), gap, labels,
        len(cells), hostile, projected, cells_used, common, margin, direct,
    )


def projected_prefix(cells, L, labels, count):
    common = ((F(0),F(1)),)
    for index, cell in enumerate(cells[:count],1):
        local = P.merge_fraction([interval for label in labels for interval in P.phase_danger(cell,label,L)])
        common = P.intersect_fraction(common,local)
    return 1-P.interval_mass(common), index, common


TASKS = tuple(
    [("ORDINARY",case) for case in ORDINARY_CASES]
    + [("HIGH-NEGATIVE",case) for case in HIGH_NEGATIVE_CASES]
    + [("HIGH-PACKET",EXCEPTIONAL_CASE)]
)


def profile(task):
    route, case = task
    if route == "ORDINARY":
        return ordinary_profile(case)
    if route == "HIGH-NEGATIVE":
        return forced_high_profile(case)
    require(route == "HIGH-PACKET", (task,"route"))
    return exceptional_profile(case)


def render(profiles):
    require(len(profiles) == 15, "profile count changed")
    ordinary = tuple(row for row in profiles if row[0] == "ORDINARY")
    negative = tuple(row for row in profiles if row[0] == "HIGH-NEGATIVE")
    exceptional = tuple(row for row in profiles if row[0] == "HIGH-PACKET")
    require((len(ordinary),len(negative),len(exceptional)) == (9,5,1), "route partition changed")
    scalar = sum(row[12][0] for row in ordinary)
    crude = sum(row[12][1] for row in ordinary)
    status = sum(row[12][2] for row in ordinary)
    states = sum(row[12][3] for row in ordinary)
    packets = sum(row[15] for row in ordinary)
    kills = sum(row[16] for row in ordinary)
    require((scalar,crude,status,states,packets,kills) == (929,0,774,155,286,286), "ordinary ledger changed")
    require(all(row[-1] < 0 for row in negative), "negative HIGH row survived")
    require(exceptional[0][-2] == F(849,1547) and exceptional[0][-1] == F(1), "exception changed")
    prior_ledger = 2_239_789
    final_ledger = prior_ledger - 15
    profile_digest = sha256(repr(tuple(profiles)).encode()).hexdigest()
    if EXPECTED_PROFILE_SHA256 is not None:
        require(profile_digest == EXPECTED_PROFILE_SHA256, "profile digest changed")
    semantic_payload = (
        EXPECTED_HASHES, QUANTIFIER, HIGH_WALL_RATIO, SUFFIX_SLOTS, EXPECTED_HEIGHTS,
        ORDINARY_CASES, HIGH_NEGATIVE_CASES, EXCEPTIONAL_CASE,
        (scalar,crude,status,states,packets,kills), tuple(row[-1] for row in negative),
        exceptional[0][-2:], prior_ledger, final_ledger, profile_digest,
    )
    semantic_digest = sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_digest == EXPECTED_SEMANTIC_SHA256, "semantic digest changed")
    lines = [
        "LRC14 projected k=2 hybrid closure at z1=1736",
        *(f"dependency_sha256={path.name}:{file_sha256(path)}" for path in PATHS),
        f"scope=all fifteen atlas rows;all {QUANTIFIER};no finite label horizon",
        f"atlas=height_counts:{EXPECTED_HEIGHTS};top_rows:15;ordinary:9;HIGH-TAIL:6;empty:1733..1735,1737..1742;next:1732",
        f"high_tail=projected_wall_ratio:{ftext(HIGH_WALL_RATIO)};suffix_slots:{SUFFIX_SLOTS};classification:pinned",
        "routes=9 all-label ray/status/projected;5 exact forced-high scalar;1 finite-cutoff high packet",
        f"ordinary_counts=scalar:{scalar};crude:{crude};status:{status};status_survivors:{states};literal_packets:{packets};projected_kills:{kills};survivors:0",
    ]
    for row in profiles:
        route = row[0]
        if route == "ORDINARY":
            (
                _route,first,body,h,components,L,lower,first_delta,first_d,ray_digest,divisors,trials,counts,stage,
                body_cells,packet_count,killed_count,min_margin,max_cells,min_row,direct,state_digest,
            ) = row
            lines.append(
                f"CASE;route=ORDINARY;z1={first};E={','.join(map(str,body))};h={ftext(h)};r={components};L={L};"
                f"lower={ftext(lower)};delta1={ftext(first_delta)};first_d={first_d};ray_sha256={ray_digest};"
                f"denominators={divisors};trials={trials};counts={counts};stage_sha256={stage};body_cells={body_cells};"
                f"packets={packet_count};kills={killed_count};min_margin={ftext(min_margin)};max_prefix_cells={max_cells};"
                f"minimum_row={min_row};direct_control_mass={ftext(direct)};state_sha256={state_digest};conclusion=EMPTY"
            )
        elif route == "HIGH-NEGATIVE":
            (
                _route,first,body,h,components,L,high_floor,first_delta,lower,positive,negative_count,zero,ray_digest,
                rank4,omitted,best_high,branch,constrained,upper,gap,
            ) = row
            lines.append(
                f"CASE;route=HIGH-NEGATIVE;z1={first};E={','.join(map(str,body))};h={ftext(h)};r={components};L={L};"
                f"high_floor={high_floor};delta1={ftext(first_delta)};lower={ftext(lower)};ray_signs=+{positive}/-{negative_count}/0:{zero};"
                f"ray_sha256={ray_digest};rank4={rank4};omitted={omitted};best_high={best_high};branch={branch};"
                f"constrained={constrained};upper={ftext(upper)};gap={ftext(gap)};conclusion=SCALAR-EMPTY"
            )
        else:
            lines.append(
                f"CASE;route=HIGH-PACKET;z1={row[1]};E={','.join(map(str,row[2]))};h={ftext(row[3])};r={row[4]};L={row[5]};"
                f"high_floor={row[6]};delta1={ftext(row[7])};lower={ftext(row[8])};target={ftext(row[9])};"
                f"ray_signs={row[10]};ray_sha256={row[11]};top3={row[12]};threshold={ftext(row[13])};"
                f"candidates={row[14]};candidate_sha256={row[15]};high_candidates={row[16]};nonempty_pairs={row[17]};"
                f"packets={row[18]};packet_sha256={row[19]};scalar_gap={ftext(row[20])};labels={row[21]};body_cells={row[22]};"
                f"hostile_13_cell_lower={ftext(row[23])};projected_lower={ftext(row[24])};cells_used={row[25]};"
                f"common_danger={row[26]};margin={ftext(row[27])};direct_mass={ftext(row[28])};conclusion=EMPTY"
            )
    lines.extend(
        (
            "global_z1736_rows=15;empty=15;survivors=0;partition_complete=yes",
            f"necessary_ledger={prior_ledger}-15={final_ledger}",
            "splice=no occupied height at 1733..1735;no interpolation",
            "consequence=projected k=2 first drift label z1<=1732",
            f"profile_sha256={profile_digest}",
            f"semantic_sha256={semantic_digest}",
            "all_exact_controls=PASS",
        )
    )
    return "\n".join(lines) + "\n"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers",type=int,default=min(8,mp.cpu_count() or 1))
    parser.add_argument("--hash-seed",type=int,default=0)
    parser.add_argument("--output",type=Path,default=OUTPUT_PATH)
    args = parser.parse_args()
    require(args.workers >= 1 and args.hash_seed >= 0, "bad replay options")
    os.environ["PYTHONHASHSEED"] = str(args.hash_seed)
    if args.workers == 1:
        profiles = [profile(task) for task in TASKS]
    else:
        with mp.get_context("spawn").Pool(min(args.workers,len(TASKS))) as pool:
            profiles = list(pool.imap(profile,TASKS,chunksize=1))
    order = {task:index for index,task in enumerate(TASKS)}
    profiles.sort(key=lambda row: order[(row[0],(row[1],row[2]))])
    output = render(tuple(profiles))
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(output,encoding="utf-8",newline="\n")
    print(output,end="")


if __name__ == "__main__":
    main()
