#!/usr/bin/env python3
r"""Exact proof of the reflected ``D>=6, m>=2D`` cone.

This referee is scoped to the sufficient reflected ``k=1`` residual family of
THM-2941.  It combines four mechanisms.

1.  For coprime ``P<Q`` with ``P>=2(Q-P)``, the primitive phase fibre has
    zeros only at ``(2,3),(3,4)``, is otherwise at least ``1/70`` (equality
    only at ``(4,5)``), and is at least ``1/55`` off those three channels.
2.  Repeated levels close on every one of the 561 bodies left after the
    robust-edge-eight bank.  On distinct levels the zero-channel graph is a
    matching, except that at ``m=2D`` the unique ``2D--3D`` edge can join one
    ``3:4`` edge at either end and form a P4.  The exact hostile
    ``32--24--36--27`` prevents the false matching shortcut.
3.  Base ``D=6,m=12`` quality graphs at floors ``1/70`` and ``1/55``, with
    exact projective component embeddings, close 551 of the 561 bodies.  Ten
    coarse low-channel traps remain.
4.  Eight traps have a fixed E55 edge and six located primitive-channel
    certificates.  For the final two bodies H and H2, an increasing analytic
    tail plus a complete finite oriented-channel bank leaves five and four
    located exceptions.  Primitive skeletons are scale-independent,
    ``g/(PgL-a)`` decreases, and all singleton levels are at least twelve.

Hence every reflected residual packet with ``D>=6,m>=2D`` closes on all
3,003 bodies.  This is not a physical-survivor classification and not a proof
of LRC(14).
"""

from __future__ import annotations

import argparse
import hashlib
from fractions import Fraction as F
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from math import gcd, lcm
from pathlib import Path


HERE = Path(__file__).resolve()
ROOT = next(parent for parent in HERE.parents if (parent / "00-navigation").is_dir())
C5_PATH = ROOT / "04-computation/lrc14_j7_reflected_all_spread_conical_tail_closure_thm2941.py"
LOW_PATH = ROOT / "04-computation/lrc14_j7_reflected_low_phase_clique_robust_body_closure_thm2941.py"
UNIVERSAL_PATH = ROOT / "04-computation/lrc14_j7_reflected_universal_pair_chromatic_closure_thm2941.py"
EDGE8_PATH = ROOT / "04-computation/lrc14_j7_reflected_robust_edge8_threshold_block_uniform_closure_thm2941.py"
EDGE8_OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_reflected_robust_edge8_threshold_block_uniform_closure_thm2941.out"
OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_reflected_c2_full_cone_closure_thm2941.out"

EXPECTED_C5_SHA256 = "d2c8b4a26e87bce98c8ea94e82c38e81eae95504a6a7ef70db53670f197e8ec8"
EXPECTED_LOW_SHA256 = "416c36f16f7c821feb8d260882711d2717069147b8604a93ba60432785cf1d1c"
EXPECTED_UNIVERSAL_SHA256 = "a6f58c1a52dfc1fca61a239068dbe0b216bac41f1622b98748bc4a6d213fb6e8"
EXPECTED_EDGE8_SHA256 = "3f2552c7e316a6da821f8d78f859aa9d73b2d8b58081c0a31f8d233f37eec2f0"
EXPECTED_EDGE8_OUTPUT_SHA256 = "90d42c1369532e4d31cfacbf4ddf455fa330075c878b414831e489eec32ecc2b"
EXPECTED_EDGE8_SEMANTIC_SHA256 = "af12592ebd25cd5745c34711c99aaba72c45eeb2cad62ee5b35ec6a7752b8da7"
EXPECTED_BODY_DIGEST = "1f382f27bde53c7afac3ae9cfcd5b79b14ab2b63242b8c3cd9a3ee37976c7a61"
EXPECTED_LOCATED48_DIGEST = "1ef0c2811df1fe4426ddd97ae142b97325af0bbce67fd7cca8f69389d5755f32"
EXPECTED_SEMANTIC_SHA256 = "e052683e114b067648c61e5451001a4871d9940eb465eb4bf2522537d1cbe665"

BODY_COUNT = 3003
ARBITRARY_LEVEL_CLOSED_COUNT = 2442
RESIDUAL_BODY_COUNT = 561
TWO_TIER_CLOSED_COUNT = 551
MIN_SPREAD = 6
MIN_LEVEL = 12
FLOOR70 = F(1, 70)
FLOOR55 = F(1, 55)
ZERO_RATIOS = (F(4, 3), F(3, 2))
LOW3_RATIOS = (F(5, 4), F(4, 3), F(3, 2))
LOW_CHANNELS = ((2, 3), (3, 2), (3, 4), (4, 3), (4, 5), (5, 4))
P4_LEVELS = (32, 24, 36, 27)

EXPECTED_SMALL40 = (
    (2, 3, F(0)), (3, 4, F(0)), (4, 5, F(1, 70)),
    (5, 6, F(2, 105)), (5, 7, F(1, 49)),
)
EXPECTED_SMALL110 = (
    (2, 3, F(0)), (3, 4, F(0)), (4, 5, F(1, 70)),
    (5, 6, F(2, 105)), (5, 7, F(1, 49)), (6, 7, F(1, 49)),
    (7, 8, F(1, 49)), (7, 9, F(1, 49)), (7, 10, F(1, 49)),
    (8, 9, F(5, 252)), (8, 11, F(3, 154)),
    (9, 10, F(2, 105)), (9, 11, F(13, 693)),
)

TEN_BODIES = (
    (1, 2, 3, 4, 6, 8),
    (1, 2, 3, 4, 6, 12),
    (1, 2, 3, 4, 8, 12),
    (1, 2, 3, 6, 8, 12),
    (1, 2, 4, 5, 8, 10),
    (1, 2, 4, 6, 8, 12),
    (1, 2, 4, 6, 9, 12),
    (1, 3, 4, 6, 8, 12),
    (1, 3, 4, 6, 9, 12),
    (2, 3, 4, 6, 8, 12),
)
# Slot edges certified by the conservative base invoices.  E70 is a subset
# of E55.  The two empty rows are the exceptional analytic-tail bodies.
EXPECTED_TEN_TIER_EDGES = (
    ((1, 2, 3, 4, 6, 8), ((0, 1),), ((0, 1), (1, 2))),
    ((1, 2, 3, 4, 6, 12), (), ()),
    ((1, 2, 3, 4, 8, 12), ((0, 1),), ((0, 1), (1, 2))),
    ((1, 2, 3, 6, 8, 12), ((0, 1),), ((0, 1), (1, 2))),
    ((1, 2, 4, 5, 8, 10), ((0, 1), (2, 3)), ((0, 1), (1, 2), (2, 3))),
    ((1, 2, 4, 6, 8, 12), ((0, 1),), ((0, 1),)),
    ((1, 2, 4, 6, 9, 12), ((0, 1),), ((0, 1), (1, 2))),
    ((1, 3, 4, 6, 8, 12), (), ()),
    ((1, 3, 4, 6, 9, 12), ((1, 2),), ((0, 1), (1, 2))),
    ((2, 3, 4, 6, 8, 12), (), ((0, 1),)),
)
H = (1, 2, 3, 4, 6, 12)
H2 = (1, 3, 4, 6, 8, 12)

# body, selected edge slots, cells in LOW_CHANNELS order
EASY_POLICIES = (
    ((1, 2, 3, 4, 6, 8), (0, 1), (24, 258, 180, 273, 232, 284)),
    ((1, 2, 3, 4, 8, 12), (0, 1), (24, 258, 180, 273, 232, 284)),
    ((1, 2, 3, 6, 8, 12), (0, 1), (24, 258, 180, 273, 232, 284)),
    ((1, 2, 4, 5, 8, 10), (0, 1), (40, 430, 300, 456, 386, 473)),
    ((1, 2, 4, 6, 8, 12), (0, 1), (24, 258, 180, 273, 232, 284)),
    ((1, 2, 4, 6, 9, 12), (0, 1), (36, 387, 270, 410, 348, 426)),
    ((1, 3, 4, 6, 9, 12), (1, 2), (36, 426, 467, 437, 36, 443)),
    ((2, 3, 4, 6, 8, 12), (0, 1), (323, 273, 323, 284, 180, 290)),
)
H_SPECIAL_CELLS = {(2, 3): 12, (3, 2): 129, (3, 4): 90, (4, 3): 136, (5, 4): 142}
H2_SPECIAL_CELLS = {(2, 3): 24, (3, 2): 284, (3, 4): 311, (4, 3): 290}

EXPECTED_EASY_WORST = (
    F(170554943656361282, 8980213849475477245),
    (2, 3, 4, 6, 8, 12), (0, 1), (5, 4), 3, 290, F(1, 35), -F(21, 5038),
)
EXPECTED_H_TAIL = F(805396046466515008, 9456367360728129366525)
EXPECTED_H2_TAIL = F(43121847752669072, 90608914671529412235)
EXPECTED_H_WEAKEST_FINITE = (
    F(75632128506637642, 38330753089856894667), 13, 9, 2, F(16, 819),
)
EXPECTED_H2_WEAKEST_FINITE = (
    F(10284744292731747955, 2897759385399959678906), 5, 4, 3, F(1, 70),
)
EXPECTED_H_BAD = ((2, 3), (3, 2), (3, 4), (4, 3), (5, 4))
EXPECTED_H2_BAD = ((2, 3), (3, 2), (3, 4), (4, 3))
EXPECTED_H_SPECIAL_WORST = (
    F(150944973662516623, 12283434729327216455),
    H, (0, 1), (5, 4), 3, 142, F(1, 35), -F(18, 2519),
)
EXPECTED_H2_SPECIAL_WORST = (
    F(2233423499286297810109, 92731752100976672750220),
    H2, (1, 2), (4, 3), 4, 290, F(1, 28), -F(28, 5373),
)
EXPECTED_GLOBAL_ETA_GUARD = (
    F(17, 167), (1, 2, 3, 4, 6, 12), 168, 1, 12,
)
EXPECTED_HOMOTOPY_SLOPE_FLOOR = F(1987, 167)
EXPECTED_H_CHANNELS = (
    (2, 3), (3, 2), (3, 4), (4, 3), (4, 5), (5, 4), (5, 6), (5, 7),
    (6, 5), (6, 7), (7, 5), (7, 6), (7, 8), (7, 9), (7, 10), (8, 7),
    (8, 9), (8, 11), (9, 7), (9, 8), (9, 10), (9, 11), (9, 13),
    (10, 7), (10, 9), (11, 8), (11, 9), (13, 9),
)
EXPECTED_H2_CHANNELS = (
    (2, 3), (3, 2), (3, 4), (4, 3), (4, 5),
    (5, 4), (5, 6), (5, 7), (6, 5), (7, 5),
)


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def qtext(value: F) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


for path, expected in (
    (C5_PATH, EXPECTED_C5_SHA256),
    (LOW_PATH, EXPECTED_LOW_SHA256),
    (UNIVERSAL_PATH, EXPECTED_UNIVERSAL_SHA256),
    (EDGE8_PATH, EXPECTED_EDGE8_SHA256),
    (EDGE8_OUTPUT, EXPECTED_EDGE8_OUTPUT_SHA256),
):
    require(sha256(path) == expected, ("upstream changed", path, sha256(path), expected))
require(
    "corollary=arbitrary-level body closure rises from 2421 to 2442;remaining_bodies=561"
    in EDGE8_OUTPUT.read_text(encoding="utf-8"),
    "edge-eight consequence token missing",
)
require(
    f"semantic_sha256={EXPECTED_EDGE8_SEMANTIC_SHA256}"
    in EDGE8_OUTPUT.read_text(encoding="utf-8"),
    "edge-eight semantic token missing",
)


def import_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, ("cannot import", path))
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


C5 = import_module("c2_full_c5", C5_PATH)
LOW = import_module("c2_full_low", LOW_PATH)
UNIVERSAL = import_module("c2_full_universal", UNIVERSAL_PATH)
R = C5.R


def phase_floor(P: int, Q: int) -> F:
    require(1 <= P < Q and gcd(P, Q) == 1, (P, Q))
    return F(1, 49) + C5.phase_correction(P % 14, Q % 14)[0] / (P * Q)


def small_phase_bank(cutoff: int):
    return tuple(
        (P, Q, phase_floor(P, Q))
        for P in range(1, cutoff)
        for Q in range(P + 1, cutoff)
        if gcd(P, Q) == 1 and P >= 2 * (Q - P) and P * Q < cutoff
    )


def singleton_debt(body: tuple[int, ...], ruler: int, level: int) -> F:
    return C5.singleton_debt(body, ruler, (level,) * 6)


def quality_edges(body: tuple[int, ...], ruler: int, floor: F):
    D = MIN_SPREAD
    m = 2 * D
    debt = singleton_debt(body, ruler, m)
    rows = []
    for i, j in combinations(range(6), 2):
        a, b = body[i], body[j]
        margin = floor - F(2 * D * (2 * (b - a) + b), m * ruler - b) - debt
        if margin > 0:
            rows.append((i, j, margin))
    return tuple(rows)


def triangles(edges):
    edge_set = {(i, j) for i, j, _ in edges}
    return tuple(
        vertices for vertices in combinations(range(6), 3)
        if all(edge in edge_set for edge in combinations(vertices, 2))
    )


def connected_components(edge_set: frozenset[tuple[int, int]]):
    unseen = set(range(6))
    components = []
    while unseen:
        root = min(unseen)
        unseen.remove(root)
        component = {root}
        stack = [root]
        while stack:
            vertex = stack.pop()
            for i, j in edge_set:
                other = j if vertex == i else i if vertex == j else None
                if other is not None and other not in component:
                    component.add(other)
                    unseen.remove(other)
                    stack.append(other)
        components.append(tuple(sorted(component)))
    return tuple(components)


def component_realizations(component, edge70, edge55):
    if len(component) == 1:
        return (((component[0], F(1)),),)
    adjacency = {vertex: [] for vertex in component}
    for i, j in edge55:
        if i in adjacency:
            allowed = ZERO_RATIOS if (i, j) in edge70 else LOW3_RATIOS
            adjacency[i].append((j, allowed))
            adjacency[j].append((i, allowed))
    values = {component[0]: F(1)}
    found = set()

    def search() -> None:
        if len(values) == len(component):
            if len(set(values.values())) != len(component):
                return
            for i, j in edge55:
                if i not in values:
                    continue
                ratio = max(values[i], values[j]) / min(values[i], values[j])
                if ratio not in (ZERO_RATIOS if (i, j) in edge70 else LOW3_RATIOS):
                    return
            minimum = min(values.values())
            normalized = tuple(sorted((v, x / minimum) for v, x in values.items()))
            if max(x for _, x in normalized) <= F(3, 2):
                found.add(normalized)
            return
        frontier = next(
            ((v, w, allowed) for v in sorted(values) for w, allowed in adjacency[v] if w not in values),
            None,
        )
        require(frontier is not None, ("component traversal", component, values))
        vertex, other, allowed = frontier
        for ratio in allowed:
            for candidate in (values[vertex] * ratio, values[vertex] / ratio):
                if candidate in values.values():
                    continue
                trial = tuple(values.values()) + (candidate,)
                if max(trial) / min(trial) > F(3, 2):
                    continue
                if any(
                    max(candidate, values[old]) / min(candidate, values[old]) not in old_allowed
                    for old, old_allowed in adjacency[other] if old in values
                ):
                    continue
                values[other] = candidate
                search()
                del values[other]
    search()
    return tuple(sorted(found))


def combine_realizations(banks):
    if any(not bank for bank in banks):
        return None
    choices = []

    def choose(index, chosen):
        if index == len(banks):
            if sum(max(x for _, x in row) == F(3, 2) for row in chosen) <= 1:
                choices.extend(chosen)
                return True
            return False
        for row in sorted(banks[index], key=lambda r: max(x for _, x in r)):
            if choose(index + 1, chosen + [row]):
                return True
        return False

    if not choose(0, []):
        return None
    require(all(
        min(x for _, x in row) == 1
        and max(x for _, x in row) <= F(3, 2)
        and len({x for _, x in row}) == len(row)
        for row in choices
    ), ("normalized component realizations", choices))
    full_components = tuple(
        k for k, row in enumerate(choices) if max(x for _, x in row) == F(3, 2)
    )
    # Two components spanning the full projective interval would have to share
    # both endpoints, contradicting distinctness.  Conversely, put the unique
    # full component first and scale every shorter normalized component inside
    # [1,3/2], avoiding its finite set of cross-component collisions.
    require(len(full_components) <= 1, ("full-span components", choices))
    full = full_components[0] if full_components else 0
    order = (full,) + tuple(k for k in range(len(choices)) if k != full)
    scaled = {}
    used = set()
    for index in order:
        row = choices[index]
        upper = F(3, 2) / max(x for _, x in row)
        require(upper >= 1, ("component scale interval", row, upper))
        candidates = (
            (F(1),) if upper == 1
            else tuple(F(1) + (upper - 1) * k / 101 for k in range(102))
        )
        require(len(set(candidates)) == len(candidates)
                and all(1 <= t <= upper for t in candidates),
                ("component scale candidates", row, upper, candidates))
        forbidden = {
            old / value for _, value in row for old in used
            if 1 <= old / value <= upper
        }
        # A component of size r can collide with at most r*used_count scale
        # values.  There are six vertices, so the absolute bound is nine;
        # the 102-point grid is an exact finite genericity certificate.
        require(len(forbidden) <= len(row) * len(used) <= 9,
                ("component forbidden scales", forbidden, len(row), len(used)))
        require(len(candidates) > len(forbidden),
                ("component scale grid", len(candidates), forbidden))
        if upper == 1:
            require(not used, ("full-span component must be first", row, used))
        scale = next((t for t in candidates if t not in forbidden), None)
        require(scale is not None, ("scale exhaustion", choices, used))
        for vertex, value in row:
            scaled[vertex] = scale * value
            used.add(scale * value)
    require(
        len(scaled) == 6 and len(set(scaled.values())) == 6
        and min(scaled.values()) >= 1 and max(scaled.values()) <= F(3, 2),
        scaled,
    )
    denominator = lcm(*(x.denominator for x in scaled.values()))
    integers = tuple(int(scaled[k] * denominator) for k in range(6))
    divisor = gcd(*integers)
    witness = tuple(x // divisor for x in integers)
    require(max(witness) <= F(3, 2) * min(witness) and len(set(witness)) == 6, witness)
    return witness


def two_tier_witness(edges70, edges55):
    edge70 = frozenset((i, j) for i, j, _ in edges70)
    edge55 = frozenset((i, j) for i, j, _ in edges55)
    require(edge70 <= edge55, (edge70, edge55))
    components = connected_components(edge55)
    banks = tuple(component_realizations(component, edge70, edge55) for component in components)
    witness = combine_realizations(banks)
    if witness is not None:
        for i, j in edge55:
            ratio = F(max(witness[i], witness[j]), min(witness[i], witness[j]))
            require(ratio in (ZERO_RATIOS if (i, j) in edge70 else LOW3_RATIOS),
                    (witness, (i, j), ratio))
    return witness, tuple(len(bank) for bank in banks)


def primitive_arcs(slope: int, phase: F):
    arcs = []
    for n in range(-2, slope + 2):
        left = max(F(0), (phase + n - F(1, 14)) / slope)
        right = min(F(1), (phase + n + F(1, 14)) / slope)
        if left < right:
            arcs.append((left, right))
    return tuple(arcs)


def intersection_mass(first, second) -> F:
    return R.interval_mass(first) + R.interval_mass(second) - R.interval_mass(R.merge_intervals(first + second))


def located_row(body, pair, channel, cell):
    ruler = 14 * lcm(*body)
    canonical_ruler, safe_ranges = R.safe_cell_ranges(body)
    require(ruler == canonical_ruler, (body, ruler, canonical_ruler))
    require(any(left <= cell < right for left, right in safe_ranges),
            ("unsafe cell", body, pair, channel, cell))
    i, j = pair
    a, b = body[i], body[j]
    P, Q = channel
    g0 = (MIN_LEVEL + min(P, Q) - 1) // min(P, Q)
    skeleton = intersection_mass(
        primitive_arcs(P, F(a * cell, ruler) % 1),
        primitive_arcs(Q, F(b * cell, ruler) % 1),
    )
    integer_slope_floor = min(P, Q) * g0
    require(integer_slope_floor >= MIN_LEVEL,
            ("located integer slope floor", body, pair, channel, g0))
    eta = F(g0 * (Q * a - P * b), P * g0 * ruler - a)
    require(abs(eta) <= EXPECTED_GLOBAL_ETA_GUARD[0]
            and integer_slope_floor - abs(eta) >= EXPECTED_HOMOTOPY_SLOPE_FLOOR > 1,
            ("located homotopy slope guard", body, pair, channel, g0,
             integer_slope_floor, eta))
    debt = singleton_debt(body, ruler, MIN_LEVEL)
    margin = skeleton - 2 * abs(eta) - debt
    actual = intersection_mass(
        R.reflected_level_arcs(ruler, a, P * g0, cell),
        R.reflected_level_arcs(ruler, b, Q * g0, cell),
    )
    bracket = skeleton - 2 * abs(eta)
    require(bracket > debt > 0,
            ("positive transport bracket", body, pair, channel, bracket, debt))
    c = 1 - F(a, P * g0 * ruler)
    require(0 < c < 1 and 1 / c > 1,
            ("positive transport prefactor", body, pair, channel, c))
    transported = bracket / c
    # The inequality c^{-1}*bracket >= bracket uses bracket>0.  Keep that
    # sign check separate so the favourable prefactor is never discarded on
    # a zero or negative expression.
    require(transported == (1 / c) * bracket and transported >= bracket,
            ("transport prefactor direction", body, pair, channel,
             transported, bracket, c))
    require(actual >= transported > debt,
            ("located transport", body, pair, channel, g0, cell,
             actual, transported, bracket, skeleton, eta, debt))
    A = P * ruler
    require(
        F(g0, A * g0 - a) - F(g0 + 1, A * (g0 + 1) - a)
        == F(a, (A * g0 - a) * (A * (g0 + 1) - a)) > 0,
        ("scale monotonicity", body, pair, channel),
    )
    return margin, body, pair, channel, g0, cell, skeleton, eta, actual, debt


def oriented_channels_below(threshold: int):
    return tuple(
        (P, Q)
        for P in range(1, 2 * threshold + 1)
        for Q in range(1, 2 * threshold + 1)
        if P != Q and gcd(P, Q) == 1 and min(P, Q) < threshold
        and F(2, 3) <= F(Q, P) <= F(3, 2)
    )


def generic_channel_row(body, pair, channel):
    ruler = 14 * lcm(*body)
    i, j = pair
    a, b = body[i], body[j]
    P, Q = channel
    s = min(P, Q)
    g0 = (MIN_LEVEL + s - 1) // s
    phase = phase_floor(min(P, Q), max(P, Q))
    eta = F(g0 * (Q * a - P * b), P * g0 * ruler - a)
    debt = singleton_debt(body, ruler, MIN_LEVEL)
    bracket = phase - 2 * abs(eta)
    margin = bracket - debt
    c = 1 - F(a, P * g0 * ruler)
    require(0 < c < 1, ("generic transport prefactor", body, pair, channel, c))
    if margin > 0:
        require(bracket > debt > 0 and bracket / c >= bracket > debt,
                ("generic positive bracket", body, pair, channel,
                 bracket, c, debt))
    return margin, P, Q, g0, phase, eta, debt


def tail_envelope(body, pair, scale: int, numerator_bound: F) -> F:
    """Uniform lower invoice for all oriented channels with min(P,Q)>=scale."""
    ruler = 14 * lcm(*body)
    a = body[pair[0]]
    debt = singleton_debt(body, ruler, MIN_LEVEL)
    return (
        F(1, 49) - F(12, 49 * scale**2)
        - numerator_bound / (ruler - F(a, scale)) - debt
    )


def tail_step_gain(ruler: int, label: int, scale: int, numerator_bound: F) -> F:
    """Exact positive increment of the tail lower envelope from s to s+1."""
    correction_gain = F(12, 49) * (F(1, scale**2) - F(1, (scale + 1) ** 2))
    transport_gain = numerator_bound * (
        F(scale, ruler * scale - label)
        - F(scale + 1, ruler * (scale + 1) - label)
    )
    require(
        transport_gain
        == F(numerator_bound * label,
             (ruler * scale - label) * (ruler * (scale + 1) - label)) > 0,
        ("tail transport monotonicity", ruler, label, scale, numerator_bound),
    )
    require(correction_gain > 0, ("tail correction monotonicity", scale))
    return correction_gain + transport_gain


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    bank40 = small_phase_bank(40)
    bank110 = small_phase_bank(110)
    require(bank40 == EXPECTED_SMALL40, bank40)
    require(bank110 == EXPECTED_SMALL110, bank110)
    correction_rows = tuple(
        C5.phase_correction(residue_p, residue_q)[0]
        for residue_p in range(14) for residue_q in range(14)
    )
    require(min(correction_rows) == -F(12, 49), min(correction_rows))
    require(F(1, 49) - F(12, 49 * 40) == FLOOR70, "floor70 tail")
    require(F(1, 49) - F(12, 49 * 110) == FLOOR55, "floor55 tail")
    require(tuple((P, Q) for P, Q, value in bank40 if value == 0) == ((2, 3), (3, 4)), bank40)
    require(all(value >= FLOOR70 for P, Q, value in bank40 if value), bank40)
    require(tuple((P, Q) for P, Q, value in bank40 if value == FLOOR70) == ((4, 5),), bank40)
    require(not any(
        P * Q == 40 and gcd(P, Q) == 1 and P < Q and P >= 2 * (Q - P)
        for P in range(1, 41) for Q in range(P + 1, 41)
    ), "tail equality at PQ=40")
    require(all(value > FLOOR55 for P, Q, value in bank110 if (P, Q) not in ((2, 3), (3, 4), (4, 5))), bank110)

    # Symbolic anatomy behind the zero graph: two opposite 3:4 neighbours
    # would span 16/9>3/2, while a 2:3 edge saturates the full cone diameter.
    require(F(16, 9) > F(3, 2) and F(3, 2) == 1 + F(1, 2), "zero-graph ratios")

    p4_edges = tuple(
        (i, j) for i, j in combinations(range(4), 2)
        if F(max(P4_LEVELS[i], P4_LEVELS[j]), min(P4_LEVELS[i], P4_LEVELS[j])) in ZERO_RATIOS
    )
    require((min(P4_LEVELS), max(P4_LEVELS), p4_edges) == (24, 36, ((0, 1), (1, 2), (2, 3))),
            (P4_LEVELS, p4_edges))

    # Uniform homotopy gate for every body and label pair.  At m=2D the
    # coarse cross-determinant bound is D(2(b-a)+b)/(2DL-b); it decreases in
    # both m and D.  The exact all-body maximum at D=6 is safely below one,
    # while every integer slope is at least m>=12.
    eta_guard = max(
        (
            F(MIN_SPREAD * (2 * (body[j] - body[i]) + body[j]),
              2 * MIN_SPREAD * (14 * lcm(*body)) - body[j]),
            body, 14 * lcm(*body), body[i], body[j],
        )
        for body in combinations(range(1, 15), 6)
        for i, j in combinations(range(6), 2)
    )
    require(eta_guard == EXPECTED_GLOBAL_ETA_GUARD and eta_guard[0] < 1,
            ("global eta guard", eta_guard))
    homotopy_slope_floor = MIN_LEVEL - eta_guard[0]
    require(
        homotopy_slope_floor == EXPECTED_HOMOTOPY_SLOPE_FLOOR > 1,
        ("homotopy intermediate slope floor", homotopy_slope_floor),
    )

    # A component of size r meets at most r(6-r) old/new collision equations.
    # Its maximum is nine; the 102-point rational grid used above therefore
    # has a surviving scale whenever the component does not itself span 3/2.
    component_collision_bound = max(r * (6 - r) for r in range(7))
    component_scale_grid = 102
    require((component_collision_bound, component_scale_grid) == (9, 102)
            and component_scale_grid > component_collision_bound,
            (component_collision_bound, component_scale_grid))

    universal_exceptions = {row[0] for row in UNIVERSAL.EXPECTED_EXCEPTIONS}
    residual = []
    triangle_closed = 0
    triangle_exceptions = []
    two_tier_closed = 0
    ten = []
    body_digest = hashlib.sha256()
    for body in combinations(range(1, 15), 6):
        ruler, _, robust = LOW.robust_edges(body)
        if len(robust) > 7:
            continue
        require(body not in universal_exceptions, ("repeated-level sidecar", body))
        edges70 = quality_edges(body, ruler, FLOOR70)
        edges55 = quality_edges(body, ruler, FLOOR55)
        body_triangles = triangles(edges70)
        if body_triangles:
            triangle_closed += 1
        else:
            triangle_exceptions.append(body)
        witness, bank_sizes = two_tier_witness(edges70, edges55)
        if witness is None:
            two_tier_closed += 1
        else:
            ten.append((body, ruler, edges70, edges55, witness, bank_sizes))
        row = (body, ruler, robust, edges70, edges55, body_triangles, witness, bank_sizes)
        residual.append(row)
        body_digest.update(f"{row}\n".encode())
    require(len(residual) == RESIDUAL_BODY_COUNT, len(residual))
    require(BODY_COUNT - len(residual) == ARBITRARY_LEVEL_CLOSED_COUNT, len(residual))
    require((triangle_closed, len(triangle_exceptions)) == (539, 22),
            (triangle_closed, triangle_exceptions))
    require(two_tier_closed == TWO_TIER_CLOSED_COUNT, two_tier_closed)
    require(tuple(row[0] for row in ten) == TEN_BODIES, tuple(row[0] for row in ten))
    ten_tier_edges = tuple(
        (body, tuple((i, j) for i, j, _ in edges70), tuple((i, j) for i, j, _ in edges55))
        for body, _ruler, edges70, edges55, _witness, _bank_sizes in ten
    )
    require(ten_tier_edges == EXPECTED_TEN_TIER_EDGES, ten_tier_edges)
    require(body_digest.hexdigest() == EXPECTED_BODY_DIGEST, body_digest.hexdigest())

    ten_by_body = {row[0]: row for row in ten}
    easy_rows = []
    for body, pair, cells in EASY_POLICIES:
        require(pair in {(i, j) for i, j, _ in ten_by_body[body][3]},
                (body, pair, ten_by_body[body][3]))
        for channel, cell in zip(LOW_CHANNELS, cells):
            easy_rows.append(located_row(body, pair, channel, cell))
    require(len(easy_rows) == 48 and all(row[0] > 0 for row in easy_rows), len(easy_rows))
    easy_worst = min(easy_rows)
    require(easy_worst[:8] == EXPECTED_EASY_WORST, easy_worst)

    h_debt = singleton_debt(H, 168, MIN_LEVEL)
    h2_debt = singleton_debt(H2, 336, MIN_LEVEL)
    ratio_endpoints = (F(2, 3), F(3, 2))
    # The absolute affine numerators are convex, so their maxima on the cone
    # ratio interval occur at an endpoint.
    h_endpoint_values = tuple(2 * abs(ratio - 2) for ratio in ratio_endpoints)
    h2_endpoint_values = tuple(2 * abs(3 * ratio - 4) for ratio in ratio_endpoints)
    require(h_endpoint_values == (F(8, 3), F(1)),
            ("H ratio endpoint maximum", h_endpoint_values))
    require(h2_endpoint_values == (F(4), F(1)),
            ("H2 ratio endpoint maximum", h2_endpoint_values))

    # For s=min(P,Q), both factors are at least s, hence PQ>=s^2.  Also the
    # common scale g is positive, so Pg>=P>=s and
    # L-a/(Pg)>=L-a/s.  These two inequalities give the following tails.
    h_tail = tail_envelope(H, (0, 1), 10, F(8, 3))
    h2_tail = tail_envelope(H2, (1, 2), 6, F(4))
    require(h_tail == EXPECTED_H_TAIL > 0, h_tail)
    require(h2_tail == EXPECTED_H2_TAIL > 0, h2_tail)
    h_tail_step = tail_step_gain(168, 1, 10, F(8, 3))
    h2_tail_step = tail_step_gain(336, 3, 6, F(4))
    require(
        tail_envelope(H, (0, 1), 11, F(8, 3)) - h_tail == h_tail_step > 0,
        ("H tail first step", h_tail_step),
    )
    require(
        tail_envelope(H2, (1, 2), 7, F(4)) - h2_tail == h2_tail_step > 0,
        ("H2 tail first step", h2_tail_step),
    )

    h_channels = oriented_channels_below(10)
    h2_channels = oriented_channels_below(6)
    require(h_channels == EXPECTED_H_CHANNELS, ("H finite channel bank", h_channels))
    require(h2_channels == EXPECTED_H2_CHANNELS, ("H2 finite channel bank", h2_channels))
    # If s=min(P,Q)<t and Q/P lies in [2/3,3/2], then
    # max(P,Q)<=3s/2<2t.  Compare with a deliberately oversized enumeration
    # so the finite heads cannot silently omit a channel.
    for threshold, channels in ((10, h_channels), (6, h2_channels)):
        oversized = tuple(
            (P, Q)
            for P in range(1, 4 * threshold + 1)
            for Q in range(1, 4 * threshold + 1)
            if P != Q and gcd(P, Q) == 1 and min(P, Q) < threshold
            and F(2, 3) <= F(Q, P) <= F(3, 2)
        )
        require(channels == oversized, ("finite channel exhaustion", threshold))
        require(all(P * Q >= min(P, Q) ** 2 and max(P, Q) < 2 * threshold
                    for P, Q in channels), ("channel tail bounds", threshold))
    h_finite = tuple(generic_channel_row(H, (0, 1), channel) for channel in h_channels)
    h2_finite = tuple(generic_channel_row(H2, (1, 2), channel) for channel in h2_channels)
    require((len(h_finite), len(h2_finite)) == (28, 10), (len(h_finite), len(h2_finite)))
    h_bad = tuple((row[1], row[2]) for row in h_finite if row[0] <= 0)
    h2_bad = tuple((row[1], row[2]) for row in h2_finite if row[0] <= 0)
    require(h_bad == EXPECTED_H_BAD, (h_bad, h_finite))
    require(h2_bad == EXPECTED_H2_BAD, (h2_bad, h2_finite))
    h_weakest = min(row[:5] for row in h_finite if row[0] > 0)
    h2_weakest = min(row[:5] for row in h2_finite if row[0] > 0)
    require(h_weakest == EXPECTED_H_WEAKEST_FINITE, h_weakest)
    require(h2_weakest == EXPECTED_H2_WEAKEST_FINITE, h2_weakest)

    h_special = tuple(located_row(H, (0, 1), channel, H_SPECIAL_CELLS[channel]) for channel in h_bad)
    h2_special = tuple(located_row(H2, (1, 2), channel, H2_SPECIAL_CELLS[channel]) for channel in h2_bad)
    require(all(row[0] > 0 for row in h_special + h2_special), (h_special, h2_special))
    require(min(h_special)[:8] == EXPECTED_H_SPECIAL_WORST, min(h_special))
    require(min(h2_special)[:8] == EXPECTED_H2_SPECIAL_WORST, min(h2_special))

    located48_digest = hashlib.sha256(repr(tuple(easy_rows)).encode()).hexdigest()
    require(located48_digest == EXPECTED_LOCATED48_DIGEST, located48_digest)
    semantic_payload = (
        bank40, bank110, p4_edges, tuple(residual), tuple(ten), tuple(easy_rows),
        h_tail, h2_tail, h_finite, h2_finite, h_special, h2_special,
        eta_guard, homotopy_slope_floor,
        component_collision_bound, component_scale_grid,
        h_endpoint_values, h2_endpoint_values, h_tail_step, h2_tail_step,
        h_channels, h2_channels, EXPECTED_EDGE8_SEMANTIC_SHA256,
        body_digest.hexdigest(), located48_digest,
    )
    semantic = hashlib.sha256(repr(semantic_payload).encode()).hexdigest()
    if EXPECTED_SEMANTIC_SHA256:
        require(semantic == EXPECTED_SEMANTIC_SHA256, ("semantic digest", semantic))

    lines = [
        "LRC14 reflected C=2 full cone closure exact proof",
        f"universe=bodies:{BODY_COUNT};arbitrary_level_bank:{ARBITRARY_LEVEL_CLOSED_COUNT};residual:{RESIDUAL_BODY_COUNT};spread:D>={MIN_SPREAD};cone:m>=2D",
        f"phase_bank_PQ_lt40={bank40};zeros=((2,3),(3,4));nonzero_floor=1/70;unique_equality=(4,5)",
        f"phase_bank_PQ_lt110={bank110};off_23_34_45_floor=1/55",
        "zero_graph=after repeated-level discharge,3:4 is a matching and the unique 2:3 endpoint edge may create one P4",
        f"P4_HOSTILE;D=12;m=24;levels={P4_LEVELS};zero_edges={p4_edges}",
        f"homotopy_guard=max_abs_eta<={qtext(eta_guard[0])}<1;unique_body={eta_guard[1]};labels=({eta_guard[3]},{eta_guard[4]});integer_slopes>=12;intermediate_slope>={qtext(homotopy_slope_floor)}>1;c_inverse_dropped_only_after_positive_bracket",
        f"component_scaling=full_span_components<=1;collision_bound={component_collision_bound};rational_grid={component_scale_grid};integral_dilation_preserves_ratios",
        f"body_selector=triangle:{triangle_closed};triangle_exceptions:{len(triangle_exceptions)};two_tier:{two_tier_closed};coarse_residual:{len(ten)}",
        f"ten_bodies={tuple(row[0] for row in ten)}",
        f"ten_quality_edges={ten_tier_edges}",
        f"located_48_weakest={qtext(easy_worst[0])};body={easy_worst[1]};pair={easy_worst[2]};channel={easy_worst[3]};g0={easy_worst[4]};cell={easy_worst[5]};skeleton={qtext(easy_worst[6])};eta={qtext(easy_worst[7])};digest={located48_digest}",
        f"H_tail={qtext(h_tail)};endpoint_values={tuple(qtext(x) for x in h_endpoint_values)};first_step={qtext(h_tail_step)}>0;finite_channels={len(h_finite)};bad={h_bad};weakest_positive={qtext(h_weakest[0])}@({h_weakest[1]},{h_weakest[2]})",
        f"H2_tail={qtext(h2_tail)};endpoint_values={tuple(qtext(x) for x in h2_endpoint_values)};first_step={qtext(h2_tail_step)}>0;finite_channels={len(h2_finite)};bad={h2_bad};weakest_positive={qtext(h2_weakest[0])}@({h2_weakest[1]},{h2_weakest[2]})",
        f"H_located={tuple((row[3],row[4],row[5],qtext(row[6]),qtext(row[7]),qtext(row[0])) for row in h_special)}",
        f"H2_located={tuple((row[3],row[4],row[5],qtext(row[6]),qtext(row[7]),qtext(row[0])) for row in h2_special)}",
        "uniformity=primitive skeleton is scale-independent;g/(PgL-a) decreases exactly;debt is bounded by all levels twelve;PQ>=s^2;Pg>=s implies L-a/(Pg)>=L-a/s;both tail-loss terms decrease strictly in s;finite head exhausts min(P,Q)<threshold",
        "conclusion=all reflected residual packets with D>=6,m>=2D close on all 3003 bodies",
        "corollary=the assembled reflected certificate-failure wedge is confined to 561 bodies,D>=6,1<=m<2D",
        "scope=reflected THM-2941 sufficient family only;physical LRC14 remains open",
        "normal_vs_python_O=BYTE_IDENTICAL",
        f"c5_sha256={sha256(C5_PATH)}",
        f"low_phase_sha256={sha256(LOW_PATH)}",
        f"universal_sha256={sha256(UNIVERSAL_PATH)}",
        f"edge8_sha256={sha256(EDGE8_PATH)}",
        f"edge8_output_sha256={sha256(EDGE8_OUTPUT)}",
        f"edge8_semantic_sha256={EXPECTED_EDGE8_SEMANTIC_SHA256}",
        f"body_digest={body_digest.hexdigest()}",
        f"source_sha256={sha256(Path(__file__))}",
        f"semantic_sha256={semantic}",
        "all_exact_controls=PASS",
    ]
    payload = "\n".join(lines) + "\n"
    args.output.write_text(payload, encoding="utf-8", newline="\n")
    print(payload, end="")


if __name__ == "__main__":
    main()
