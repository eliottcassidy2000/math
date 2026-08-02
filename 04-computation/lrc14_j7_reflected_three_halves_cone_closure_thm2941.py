#!/usr/bin/env python3
r"""Exact proof candidate for the reflected ``2m >= 3D`` cone.

This referee is scoped to the sufficient reflected ``k=1`` family assembled
in THM-2941.  It starts with the 561 bodies left by the arbitrary-level
robust-edge-eight bank and proves that all of them close whenever

    D >= 6,                 2m >= 3D.

The proof has three layers.  First, the exact primitive phase formula and a
body-specific loss threshold give a finite rational constraint problem on
the six distinct levels; it closes 540 bodies.  Second, 19 of the 21 coarse
traps have a selected label pair whose complete finite set of oriented low
channels admits an exact located cell.  Finally, two analytic tails and
complete finite oriented heads close the bodies H and H2.

This is a certificate theorem inside THM-2941, not a physical-survivor census
and not a proof of LRC(14).
"""

from __future__ import annotations

import argparse
import hashlib
from fractions import Fraction as F
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from math import gcd, isqrt, lcm
from pathlib import Path


HERE = Path(__file__).resolve()
ROOT = next(parent for parent in HERE.parents if (parent / "00-navigation").is_dir())
C2_PATH = ROOT / "04-computation/lrc14_j7_reflected_c2_full_cone_closure_thm2941.py"
C2_OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_reflected_c2_full_cone_closure_thm2941.out"
OUTPUT = ROOT / "05-knowledge/results/lrc14_j7_reflected_three_halves_cone_closure_thm2941.out"

EXPECTED_C2_SHA256 = "f7eb06c50b11fd810320146a2cded2590790fdbcf6db58cff2a358470fcfd0b6"
EXPECTED_C2_OUTPUT_SHA256 = "7449cb0ab6b969fa4653d98582423df262329e087cc73d8ea33499cb23ec8422"
EXPECTED_C2_SEMANTIC_SHA256 = "e052683e114b067648c61e5451001a4871d9940eb465eb4bf2522537d1cbe665"
EXPECTED_SEMANTIC_SHA256 = "6d2e760f3a929ad2eaa9af9eba11db1389ead59a383e11d8e8737c722b67aaaa"

BODY_COUNT = 3003
ARBITRARY_LEVEL_CLOSED_COUNT = 2442
RESIDUAL_BODY_COUNT = 561
COARSE_CLOSED_COUNT = 540
COARSE_TRAP_COUNT = 21
MIN_SPREAD = 6
MIN_LEVEL = 9
RATIO_CAP = F(5, 3)
PHASE_LIMIT = F(1, 49)
PHASE_CORRECTION_FLOOR = -F(12, 49)
COMPONENT_SCALE_GRID = 102

H = (1, 2, 3, 4, 6, 12)
H2 = (1, 3, 4, 6, 8, 12)
H_PAIR = (0, 1)
H2_PAIR = (1, 2)
H_TAIL_START = 16
H2_TAIL_START = 7
H_NUMERATOR_BOUND = F(14, 5)
H2_NUMERATOR_BOUND = F(22, 5)

EXPECTED_TRAPS = (
    (1, 2, 3, 4, 6, 8),
    (1, 2, 3, 4, 6, 9),
    H,
    (1, 2, 3, 4, 8, 12),
    (1, 2, 3, 4, 9, 12),
    (1, 2, 3, 5, 6, 10),
    (1, 2, 3, 6, 7, 14),
    (1, 2, 3, 6, 8, 12),
    (1, 2, 3, 6, 9, 12),
    (1, 2, 4, 5, 8, 10),
    (1, 2, 4, 6, 8, 12),
    (1, 2, 4, 6, 9, 12),
    (1, 2, 4, 7, 8, 14),
    (1, 2, 5, 6, 10, 12),
    (1, 2, 5, 7, 10, 14),
    H2,
    (1, 3, 4, 6, 9, 12),
    (1, 3, 5, 6, 10, 12),
    (1, 4, 5, 6, 10, 12),
    (2, 3, 4, 6, 8, 12),
    (2, 3, 4, 6, 9, 12),
)

# A single selected pair suffices on every coarse trap except H and H2.  The
# channel bank is generated from the exact body-specific threshold; it is not
# hard-coded or sampled.
SELECTED_PAIRS = {
    body: (0, 1) for body in EXPECTED_TRAPS if body not in (H, H2)
}


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


require(sha256(C2_PATH) == EXPECTED_C2_SHA256, ("C2 source changed", sha256(C2_PATH)))
require(sha256(C2_OUTPUT) == EXPECTED_C2_OUTPUT_SHA256,
        ("C2 output changed", sha256(C2_OUTPUT)))
require(f"semantic_sha256={EXPECTED_C2_SEMANTIC_SHA256}" in C2_OUTPUT.read_text(),
        "C2 semantic token missing")
C2 = import_module("three_halves_c2", C2_PATH)
R = C2.R


def phase_floor(P: int, Q: int) -> F:
    require(1 <= P < Q and gcd(P, Q) == 1, (P, Q))
    return F(1, 49) + C2.C5.phase_correction(P % 14, Q % 14)[0] / (P * Q)


def debt_at_nine(body: tuple[int, ...], ruler: int) -> F:
    return C2.singleton_debt(body, ruler, MIN_LEVEL)


def coarse_threshold(body: tuple[int, ...], ruler: int, pair: tuple[int, int]) -> F:
    """Worst ``2|eta| + debt`` over D>=6 and 2m>=3D."""
    i, j = pair
    a, b = body[i], body[j]
    delta = b - a
    coarse_eta = F(MIN_LEVEL * delta + MIN_SPREAD * b, MIN_LEVEL * ruler - b)
    return 2 * coarse_eta + debt_at_nine(body, ruler)


def product_bound(threshold: F) -> F:
    """Complete bound for a channel whose phase is at most threshold."""
    require(threshold < PHASE_LIMIT, threshold)
    return F(12, 1 - 49 * threshold)


def residual_bodies() -> tuple[tuple[int, ...], ...]:
    rows = []
    for body in combinations(range(1, 15), 6):
        _, _, robust = C2.LOW.robust_edges(body)
        if len(robust) <= 7:
            rows.append(body)
    require(len(rows) == RESIDUAL_BODY_COUNT, len(rows))
    return tuple(rows)


def primitive_universe(bound: F):
    """Every reduced unordered channel under the rational product bound."""
    rows = []
    P = 1
    while P * (P + 1) <= bound:
        for Q in range(P + 1, bound.numerator // (bound.denominator * P) + 1):
            if P * Q <= bound and gcd(P, Q) == 1 and 3 * Q <= 5 * P:
                rows.append((F(Q, P), phase_floor(P, Q), P * Q, P, Q))
        P += 1
    return tuple(rows)


def body_constraints(body, ruler, primitive_rows):
    constraints = {}
    thresholds = {}
    for pair in combinations(range(6), 2):
        threshold = coarse_threshold(body, ruler, pair)
        if threshold >= PHASE_LIMIT:
            continue
        bound = product_bound(threshold)
        allowed = frozenset(
            ratio for ratio, phase, product, _, _ in primitive_rows
            if product <= bound and phase <= threshold
        )
        constraints[pair] = allowed
        thresholds[pair] = threshold
    return constraints, thresholds


def connected_components(constraints):
    adjacency = {vertex: set() for vertex in range(6)}
    for i, j in constraints:
        adjacency[i].add(j)
        adjacency[j].add(i)
    unseen = set(range(6))
    components = []
    while unseen:
        root = min(unseen)
        unseen.remove(root)
        stack = [root]
        component = {root}
        while stack:
            vertex = stack.pop()
            for other in adjacency[vertex]:
                if other not in component:
                    component.add(other)
                    unseen.remove(other)
                    stack.append(other)
        components.append(tuple(sorted(component)))
    return tuple(components)


def component_types(component, constraints, reverse=False):
    """Exhaustively find a short and/or full-span projective realization."""
    if len(component) == 1:
        return ((component[0], F(1)),), None, 1
    adjacency = {vertex: [] for vertex in component}
    for (i, j), allowed in constraints.items():
        if i in adjacency and j in adjacency:
            adjacency[i].append((j, allowed))
            adjacency[j].append((i, allowed))
    root = min(
        component,
        key=lambda vertex: (
            min(len(allowed) for _, allowed in adjacency[vertex]),
            -len(adjacency[vertex]),
            -vertex if reverse else vertex,
        ),
    )
    values = {root: F(1)}
    found = [None, None]  # short, full
    nodes = 0

    def search() -> None:
        nonlocal nodes
        nodes += 1
        if found[0] is not None and found[1] is not None:
            return
        if len(values) == len(component):
            if len(set(values.values())) != len(values):
                return
            minimum = min(values.values())
            normalized = tuple(sorted((v, x / minimum) for v, x in values.items()))
            span = max(x for _, x in normalized)
            if span > RATIO_CAP:
                return
            for (i, j), allowed in constraints.items():
                if i not in values or j not in values:
                    continue
                ratio = max(values[i], values[j]) / min(values[i], values[j])
                if ratio not in allowed:
                    return
            found[span == RATIO_CAP] = normalized
            return

        options = []
        for vertex in component:
            if vertex in values:
                continue
            links = [(old, allowed) for old, allowed in adjacency[vertex] if old in values]
            if not links:
                continue
            old, allowed = min(links, key=lambda item: len(item[1]))
            candidates = {
                candidate
                for ratio in allowed
                for candidate in (values[old] * ratio, values[old] / ratio)
            }
            options.append((len(candidates), -len(links), vertex, links, candidates))
        require(options, ("component traversal", component, values))
        _, _, vertex, links, candidates = min(
            options,
            key=lambda row: (row[0], row[1], -row[2] if reverse else row[2]),
        )
        for candidate in sorted(candidates, reverse=reverse):
            if candidate in values.values():
                continue
            trial = tuple(values.values()) + (candidate,)
            if max(trial) / min(trial) > RATIO_CAP:
                continue
            if any(
                max(candidate, values[old]) / min(candidate, values[old]) not in allowed
                for old, allowed in links
            ):
                continue
            values[vertex] = candidate
            search()
            del values[vertex]

    search()
    return found[0], found[1], nodes


def scale_components(rows):
    """Place projective components distinctly in [1,5/3], then integralize."""
    full = tuple(k for k, row in enumerate(rows) if max(x for _, x in row) == RATIO_CAP)
    require(len(full) <= 1, ("multiple full components", rows))
    first = full[0] if full else 0
    order = (first,) + tuple(k for k in range(len(rows)) if k != first)
    scaled = {}
    used = set()
    for index in order:
        row = rows[index]
        span = max(x for _, x in row)
        upper = RATIO_CAP / span
        candidates = (
            (F(1),) if upper == 1 else
            tuple(F(1) + (upper - 1) * k / (COMPONENT_SCALE_GRID - 1)
                  for k in range(COMPONENT_SCALE_GRID))
        )
        forbidden = {
            old / value for _, value in row for old in used
            if 1 <= old / value <= upper
        }
        require(len(forbidden) <= len(row) * len(used) <= 9,
                ("collision bound", row, forbidden, used))
        require(len(candidates) > len(forbidden), (candidates, forbidden))
        if upper == 1:
            require(not used, ("full component not first", row, used))
        scale = next((value for value in candidates if value not in forbidden), None)
        require(scale is not None, ("scale exhaustion", row, forbidden))
        for vertex, value in row:
            scaled[vertex] = scale * value
            used.add(scale * value)
    require(len(scaled) == 6 and len(set(scaled.values())) == 6, scaled)
    denominator = lcm(*(value.denominator for value in scaled.values()))
    integers = tuple(int(scaled[k] * denominator) for k in range(6))
    divisor = gcd(*integers)
    witness = tuple(value // divisor for value in integers)
    require(len(set(witness)) == 6 and max(witness) <= RATIO_CAP * min(witness), witness)
    return witness


def solve_body(body, ruler, primitive_rows, reverse=False):
    constraints, thresholds = body_constraints(body, ruler, primitive_rows)
    if any(not allowed for allowed in constraints.values()):
        return None, constraints, thresholds, ("empty edge",)
    choices = []
    component_rows = []
    nodes = 0
    forced_full = 0
    for component in connected_components(constraints):
        short, full, count = component_types(component, constraints, reverse=reverse)
        nodes += count
        if short is None and full is None:
            return None, constraints, thresholds, ("empty component", component, nodes)
        if short is None:
            forced_full += 1
            component_rows.append(full)
        else:
            component_rows.append(short)
        choices.append((component, short is not None, full is not None, count))
    if forced_full > 1:
        return None, constraints, thresholds, ("two forced full components", tuple(choices), nodes)
    witness = scale_components(tuple(component_rows))
    for pair, allowed in constraints.items():
        i, j = pair
        ratio = F(max(witness[i], witness[j]), min(witness[i], witness[j]))
        require(ratio in allowed, (body, witness, pair, ratio, allowed))
    return witness, constraints, thresholds, (tuple(choices), nodes)


def primitive_arcs(slope: int, phase: F):
    return C2.primitive_arcs(slope, phase)


def intersection_mass(first, second) -> F:
    return C2.intersection_mass(first, second)


def located_best(body, pair, channel):
    """Best exact body-safe cell at the least possible common channel scale."""
    ruler, safe_ranges = R.safe_cell_ranges(body)
    i, j = pair
    a, b = body[i], body[j]
    P, Q = channel
    g0 = (MIN_LEVEL + min(P, Q) - 1) // min(P, Q)
    debt = debt_at_nine(body, ruler)
    eta = F(g0 * (Q * a - P * b), P * g0 * ruler - a)
    rows = []
    for left, right in safe_ranges:
        for cell in range(left, right):
            skeleton = intersection_mass(
                primitive_arcs(P, F(a * cell, ruler) % 1),
                primitive_arcs(Q, F(b * cell, ruler) % 1),
            )
            rows.append((skeleton - 2 * abs(eta) - debt, cell, skeleton))
    margin, cell, skeleton = max(rows)
    bracket = skeleton - 2 * abs(eta)
    c = 1 - F(a, P * g0 * ruler)
    transported = bracket / c
    actual = intersection_mass(
        R.reflected_level_arcs(ruler, a, P * g0, cell),
        R.reflected_level_arcs(ruler, b, Q * g0, cell),
    )
    require(min(P, Q) * g0 >= MIN_LEVEL and abs(eta) <= F(57, 500) < 1,
            ("homotopy gate", body, pair, channel, g0, eta))
    require(bracket > debt > 0, ("nonpositive bracket", body, pair, channel, bracket, debt))
    require(0 < c < 1 and 1 / c > 1 and transported >= bracket > debt,
            ("prefactor direction", body, pair, channel, c, bracket, debt))
    require(actual >= transported > debt,
            ("direct transport control", body, pair, channel, cell, actual, transported, debt))
    A = P * ruler
    require(
        F(g0, A * g0 - a) - F(g0 + 1, A * (g0 + 1) - a)
        == F(a, (A * g0 - a) * (A * (g0 + 1) - a)) > 0,
        ("scale monotonicity", body, pair, channel),
    )
    return margin, body, pair, channel, g0, cell, skeleton, eta, actual, debt


def oriented_channels_below(threshold: int, oversize: int = 2):
    rows = []
    for P in range(1, oversize * threshold + 1):
        for Q in range(1, oversize * threshold + 1):
            if (
                P != Q and gcd(P, Q) == 1 and min(P, Q) < threshold
                and F(3, 5) <= F(Q, P) <= F(5, 3)
            ):
                rows.append((P, Q))
    return tuple(rows)


def tail_envelope(body, pair, scale: int, numerator_bound: F) -> F:
    ruler = 14 * lcm(*body)
    a = body[pair[0]]
    return (
        F(1, 49) - F(12, 49 * scale**2)
        - numerator_bound / (ruler - F(a, scale))
        - debt_at_nine(body, ruler)
    )


def tail_step_gain(body, pair, scale: int, numerator_bound: F) -> F:
    ruler = 14 * lcm(*body)
    a = body[pair[0]]
    correction = F(12, 49) * (F(1, scale**2) - F(1, (scale + 1)**2))
    transport = numerator_bound * (
        F(scale, ruler * scale - a)
        - F(scale + 1, ruler * (scale + 1) - a)
    )
    require(
        transport == F(numerator_bound * a,
                       (ruler * scale - a) * (ruler * (scale + 1) - a)) > 0,
        ("tail transport step", body, scale),
    )
    require(correction > 0, ("tail phase step", body, scale))
    return correction + transport


def fixed_H_skeleton(channel):
    P, Q = channel
    return intersection_mass(
        primitive_arcs(P, F(155, 168)),
        primitive_arcs(Q, F(310, 168) % 1),
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()

    bodies = residual_bodies()
    universal_exceptions = {row[0] for row in C2.UNIVERSAL.EXPECTED_EXCEPTIONS}
    require(all(body not in universal_exceptions for body in bodies),
            ("repeated-level discharge failed", universal_exceptions & set(bodies)))

    # The coarse eta bound B(D,m) decreases with m.  On the real boundary
    # m=3D/2 it also decreases with D, so this handles odd D without rounding
    # a fictitious integer corner.  The actual integer condition is 2m>=3D.
    for body in bodies:
        ruler = 14 * lcm(*body)
        for i, j in combinations(range(6), 2):
            a, b = body[i], body[j]
            delta = b - a
            for D in (6, 7):
                m = F(3 * D, 2)
                B = F(m * delta + D * b, m * ruler - b)
                B_next_m = F((m + 1) * delta + D * b, (m + 1) * ruler - b)
                require(
                    B - B_next_m
                    == F(b * (D * ruler + delta),
                         (m * ruler - b) * ((m + 1) * ruler - b)) > 0,
                    ("m monotonicity", body, (i, j), D),
                )
            A = F(3, 2) * delta + b
            B6 = F(6 * A, 9 * ruler - b)
            B7 = F(7 * A, F(21, 2) * ruler - b)
            require(B6 - B7 == F(A * b, (9 * ruler - b) * (F(21, 2) * ruler - b)) > 0,
                    ("D monotonicity", body, (i, j)))

    eta_guard = max(
        (
            F(MIN_LEVEL * (b - a) + MIN_SPREAD * b, MIN_LEVEL * (14 * lcm(*body)) - b),
            body, (a, b), 14 * lcm(*body),
        )
        for body in combinations(range(1, 15), 6)
        for a, b in combinations(body, 2)
    )
    require(eta_guard == (F(57, 500), H, (1, 12), 168), eta_guard)
    require(MIN_LEVEL - eta_guard[0] == F(4443, 500) > 1, eta_guard)

    corrections = tuple(
        C2.C5.phase_correction(p, q)[0]
        for p in range(14) for q in range(14)
    )
    require(min(corrections) == PHASE_CORRECTION_FLOOR, min(corrections))

    thresholds = []
    for body in bodies:
        ruler = 14 * lcm(*body)
        for pair in combinations(range(6), 2):
            threshold = coarse_threshold(body, ruler, pair)
            if threshold < PHASE_LIMIT:
                thresholds.append((threshold, body, pair))
    maximum_bound = max((product_bound(t), body, pair, t) for t, body, pair in thresholds)
    primitive_rows = primitive_universe(maximum_bound[0])
    require(all(
        phase >= PHASE_LIMIT + PHASE_CORRECTION_FLOOR / product
        for _, phase, product, _, _ in primitive_rows
    ), "phase tail inequality")

    closed = []
    traps = []
    witness_digest = hashlib.sha256()
    verdict_digest = hashlib.sha256()
    for body in bodies:
        ruler = 14 * lcm(*body)
        witness, constraints, _, reason = solve_body(body, ruler, primitive_rows)
        reverse_witness, reverse_constraints, _, reverse_reason = solve_body(
            body, ruler, primitive_rows, reverse=True
        )
        require((witness is None) == (reverse_witness is None),
                ("search-order mismatch", body, reason, reverse_reason))
        require(constraints == reverse_constraints, ("constraint mismatch", body))
        verdict_digest.update(repr((body, witness is not None, reason, reverse_reason)).encode())
        if witness is None:
            closed.append(body)
        else:
            traps.append(body)
            witness_digest.update(repr((body, witness)).encode())
    require((len(closed), len(traps), tuple(traps))
            == (COARSE_CLOSED_COUNT, COARSE_TRAP_COUNT, EXPECTED_TRAPS),
            (len(closed), len(traps), traps))

    located = []
    policy_rows = []
    for body, pair in SELECTED_PAIRS.items():
        ruler = 14 * lcm(*body)
        constraints, _ = body_constraints(body, ruler, primitive_rows)
        require(pair in constraints and constraints[pair], ("missing policy edge", body, pair))
        unordered = tuple(sorted(
            (ratio.denominator, ratio.numerator) for ratio in constraints[pair]
        ))
        oriented = tuple(
            channel for P, Q in unordered for channel in ((P, Q), (Q, P))
        )
        rows = tuple(located_best(body, pair, channel) for channel in oriented)
        require(all(row[0] > 0 for row in rows), ("located policy failure", body, pair, rows))
        located.extend(rows)
        policy_rows.append((body, pair, len(rows), min(rows)))
    require(len(located) == 138, len(located))
    located_digest = hashlib.sha256(repr(tuple(located)).encode()).hexdigest()
    located_worst = min(located)

    H_tail = tail_envelope(H, H_PAIR, H_TAIL_START, H_NUMERATOR_BOUND)
    H2_tail = tail_envelope(H2, H2_PAIR, H2_TAIL_START, H2_NUMERATOR_BOUND)
    require(H_tail > 0 >= tail_envelope(H, H_PAIR, H_TAIL_START - 1, H_NUMERATOR_BOUND),
            ("H tail threshold", H_tail))
    require(H2_tail > 0 >= tail_envelope(H2, H2_PAIR, H2_TAIL_START - 1, H2_NUMERATOR_BOUND),
            ("H2 tail threshold", H2_tail))
    H_step = tail_step_gain(H, H_PAIR, H_TAIL_START, H_NUMERATOR_BOUND)
    H2_step = tail_step_gain(H2, H2_PAIR, H2_TAIL_START, H2_NUMERATOR_BOUND)

    H_channels = oriented_channels_below(H_TAIL_START)
    H_channels_large = oriented_channels_below(H_TAIL_START, oversize=4)
    H2_channels = oriented_channels_below(H2_TAIL_START)
    H2_channels_large = oriented_channels_below(H2_TAIL_START, oversize=4)
    require(H_channels == H_channels_large and len(H_channels) == 96,
            (len(H_channels), len(H_channels_large)))
    require(H2_channels == H2_channels_large and len(H2_channels) == 16,
            (len(H2_channels), len(H2_channels_large)))
    H_head = tuple(located_best(H, H_PAIR, channel) for channel in H_channels)
    H2_head = tuple(located_best(H2, H2_PAIR, channel) for channel in H2_channels)
    require(all(row[0] > 0 for row in H_head), ("H head", min(H_head)))
    require(all(row[0] > 0 for row in H2_head), ("H2 head", min(H2_head)))
    head_digest = hashlib.sha256(repr((H_head, H2_head)).encode()).hexdigest()

    # Fixed-cell sidecar requested by the hostile audit.  The global phase
    # bound gives >=1/126 once PQ>=20.  The complete smaller bank proves the
    # same bound at cell 155, sharply at the reverse 3:2 channel.
    fixed_small = tuple(
        (P, Q, fixed_H_skeleton((P, Q)))
        for P in range(1, 20) for Q in range(1, 20)
        if P != Q and gcd(P, Q) == 1 and P * Q < 20
        and F(3, 5) <= F(Q, P) <= F(5, 3)
    )
    _, H_safe_ranges = R.safe_cell_ranges(H)
    require(any(left <= 155 < right for left, right in H_safe_ranges),
            ("fixed H cell is not body-safe", H_safe_ranges))
    require(F(1, 49) - F(12, 49 * 20) >= F(1, 126), "fixed-cell phase tail")
    require(min(value for _, _, value in fixed_small) == F(1, 126)
            and tuple((P, Q) for P, Q, value in fixed_small if value == F(1, 126)) == ((3, 2),),
            fixed_small)
    fixed_transport_bad = []
    fixed_direct_controls = []
    H_debt = debt_at_nine(H, 168)
    for P, Q in H_channels:
        g0 = (MIN_LEVEL + min(P, Q) - 1) // min(P, Q)
        skeleton = fixed_H_skeleton((P, Q))
        eta = F(g0 * (Q - 2 * P), P * g0 * 168 - 1)
        if skeleton - 2 * abs(eta) <= H_debt:
            fixed_transport_bad.append((P, Q))
            actual = intersection_mass(
                R.reflected_level_arcs(168, 1, P * g0, 155),
                R.reflected_level_arcs(168, 2, Q * g0, 155),
            )
            require(actual > H_debt, ("fixed direct hostile", P, Q, actual, H_debt))
            fixed_direct_controls.append((P, Q, g0, actual - H_debt))
    require(tuple(fixed_transport_bad) == ((3, 2), (4, 3), (5, 3), (5, 4)),
            fixed_transport_bad)

    semantic_image = (
        tuple(bodies), maximum_bound, tuple(primitive_rows),
        tuple(closed), tuple(traps), witness_digest.hexdigest(), verdict_digest.hexdigest(),
        tuple(policy_rows), tuple(located), located_digest,
        H_tail, H2_tail, H_step, H2_step, H_channels, H2_channels,
        H_head, H2_head, head_digest, fixed_small,
        tuple(fixed_transport_bad), tuple(fixed_direct_controls), eta_guard,
    )
    semantic_sha = hashlib.sha256(repr(semantic_image).encode()).hexdigest()
    require(semantic_sha == EXPECTED_SEMANTIC_SHA256,
            ("semantic image changed", semantic_sha, EXPECTED_SEMANTIC_SHA256))

    lines = [
        "LRC14 reflected three-halves cone exact proof candidate",
        "universe=bodies:3003;arbitrary_level_bank:2442;residual:561;spread:D>=6;cone:2m>=3D",
        f"corner=real boundary m=3D/2 decreases in D;worst_D=6;worst_m=9;debt_levels=all_nine;max_abs_eta={qtext(eta_guard[0])};intermediate_slope_floor={qtext(MIN_LEVEL-eta_guard[0])}",
        f"phase=1/49+c/(PQ);c_min={qtext(PHASE_CORRECTION_FLOOR)};ratio_cap=5/3;constrained_edges={len(thresholds)};max_product_bound={qtext(maximum_bound[0])}@{maximum_bound[1]},{maximum_bound[2]};primitive_bank={len(primitive_rows)}",
        f"coarse_csp=closed:{len(closed)};traps:{len(traps)};reverse_order_agrees;component_full_span_at_most_one;collision_bound=9;scale_grid={COMPONENT_SCALE_GRID}",
        f"coarse_traps={tuple(traps)}",
        f"witness_digest={witness_digest.hexdigest()};verdict_digest={verdict_digest.hexdigest()}",
        f"located_policies={len(policy_rows)};oriented_controls={len(located)};weakest={qtext(located_worst[0])}@body={located_worst[1]},pair={located_worst[2]},channel={located_worst[3]},g0={located_worst[4]},cell={located_worst[5]};digest={located_digest}",
        f"H_tail=start:{H_TAIL_START};margin={qtext(H_tail)};first_step={qtext(H_step)};head_channels={len(H_channels)};head_weakest={qtext(min(H_head)[0])}@{min(H_head)[3]},cell={min(H_head)[5]}",
        f"H2_tail=start:{H2_TAIL_START};margin={qtext(H2_tail)};first_step={qtext(H2_step)};head_channels={len(H2_channels)};head_weakest={qtext(min(H2_head)[0])}@{min(H2_head)[3]},cell={min(H2_head)[5]};head_digest={head_digest}",
        f"fixed_H_cell=155;uniform_primitive_skeleton>=1/126;equality=(3,2);transport_bound_exceptions={tuple(fixed_transport_bad)};direct_g0_controls={tuple((P,Q,g,qtext(margin)) for P,Q,g,margin in fixed_direct_controls)}",
        "transport=oriented eta_g=g(Qa-Pb)/(PgL-a);g/(PgL-a) strictly decreases;c_inverse_dropped_only_after_positive_bracket;debt decreases from all-nine bound",
        "tail=PQ>=s^2 and Pg>=s;endpoint_numerators H:2|max(r-2)|<=14/5,H2:2|3r-4|<=22/5;both losses decrease strictly in s;heads are oversized-audited",
        "conclusion=all reflected residual packets with D>=6 and 2m>=3D close on all 3003 bodies",
        "corollary=assembled reflected certificate-failure wedge is confined to 561 bodies,D>=6,1<=m<3D/2",
        "scope=reflected THM-2941 sufficient family only;physical LRC14 remains open",
        "normal_vs_python_O=BYTE_IDENTICAL",
        f"c2_source_sha256={EXPECTED_C2_SHA256}",
        f"c2_output_sha256={EXPECTED_C2_OUTPUT_SHA256}",
        f"c2_semantic_sha256={EXPECTED_C2_SEMANTIC_SHA256}",
        f"source_sha256={sha256(HERE)}",
        f"semantic_sha256={semantic_sha}",
        "all_exact_controls=PASS",
    ]
    text = "\n".join(lines) + "\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(text, encoding="utf-8")
    print(text, end="")


if __name__ == "__main__":
    main()
