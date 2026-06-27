#!/usr/bin/env python3
"""HYP-3111 scout: Minkowski, circuit complexity, Ising, and De Moivre.

This script joins four requested outside carriers to the current LRC14 proof
frontier.  It is a route-selection scout, not a proof of LRC14.

Measured carriers:
  * a q-vector lattice sidecar with a Minkowski-ball forcing threshold;
  * a monotone proof-state circuit with size/depth/minterm data;
  * finite ferromagnetic Ising partition-polynomial zeros;
  * De Moivre's quintic fold as an exact Laurent-polynomial identity;
  * Tournament Analysis on proof-carrier sidecars rather than runners.
"""

from __future__ import annotations

import itertools
import math
from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction as F

SECTORS = 7
EPS = 1e-10


def fstr(x: F) -> str:
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def sector_of(x: F) -> int:
    return int((x % 1) * SECTORS)


def missdist(row: tuple[int, ...]) -> tuple[F, ...]:
    """Distribution of the number of empty inner sectors for an LRC row."""
    row = tuple(sorted(set(row)))
    breaks = {F(0), F(1)}
    for speed in row:
        if speed == 0:
            continue
        for m in range(SECTORS * speed + 1):
            breaks.add(F(m, SECTORS * speed))
    cuts = sorted(breaks)
    q = [F(0) for _ in range(SECTORS)]
    for lo, hi in zip(cuts, cuts[1:]):
        if hi <= lo:
            continue
        mid = (lo + hi) / 2
        occupied = {sector_of(speed * mid) for speed in row}
        empty = SECTORS - len(occupied)
        if 0 <= empty < SECTORS:
            q[empty] += hi - lo
    return tuple(q)


def q_drop_last(q: tuple[F, ...]) -> tuple[F, ...]:
    """Affine face coordinates; q_6 is determined by sum q_t = 1."""
    return tuple(q[:-1])


def vsub(a: tuple[F, ...], b: tuple[F, ...]) -> tuple[F, ...]:
    return tuple(x - y for x, y in zip(a, b))


def dot(a: tuple[F, ...], b: tuple[F, ...]) -> F:
    return sum((x * y for x, y in zip(a, b)), F(0))


def norm(v: tuple[F, ...]) -> float:
    return math.sqrt(float(dot(v, v)))


def rank_fraction(rows: list[tuple[F, ...]]) -> int:
    if not rows:
        return 0
    mat = [list(row) for row in rows if any(x != 0 for x in row)]
    if not mat:
        return 0
    m = len(mat)
    n = len(mat[0])
    rank = 0
    col = 0
    while rank < m and col < n:
        pivot = None
        for r in range(rank, m):
            if mat[r][col] != 0:
                pivot = r
                break
        if pivot is None:
            col += 1
            continue
        mat[rank], mat[pivot] = mat[pivot], mat[rank]
        scale = mat[rank][col]
        mat[rank] = [x / scale for x in mat[rank]]
        for r in range(m):
            if r == rank:
                continue
            factor = mat[r][col]
            if factor:
                mat[r] = [x - factor * y for x, y in zip(mat[r], mat[rank])]
        rank += 1
        col += 1
    return rank


def independent_basis(vectors: list[tuple[str, tuple[F, ...]]]) -> list[tuple[str, tuple[F, ...]]]:
    basis: list[tuple[str, tuple[F, ...]]] = []
    for name, vec in vectors:
        if rank_fraction([v for _, v in basis] + [vec]) > len(basis):
            basis.append((name, vec))
    return basis


def determinant_fraction(matrix: list[list[F]]) -> F:
    n = len(matrix)
    mat = [row[:] for row in matrix]
    det = F(1)
    sign = 1
    for col in range(n):
        pivot = None
        for r in range(col, n):
            if mat[r][col] != 0:
                pivot = r
                break
        if pivot is None:
            return F(0)
        if pivot != col:
            mat[col], mat[pivot] = mat[pivot], mat[col]
            sign *= -1
        pivot_val = mat[col][col]
        det *= pivot_val
        for r in range(col + 1, n):
            factor = mat[r][col] / pivot_val
            for c in range(col, n):
                mat[r][c] -= factor * mat[col][c]
    return det if sign > 0 else -det


def poly_eval(coeffs_high_to_low: list[float], z: complex) -> complex:
    out = 0j
    for coeff in coeffs_high_to_low:
        out = out * z + coeff
    return out


def durand_kerner(coeffs_high_to_low: list[float]) -> list[complex]:
    coeffs = coeffs_high_to_low[:]
    while len(coeffs) > 1 and abs(coeffs[0]) < 1e-14:
        coeffs = coeffs[1:]
    degree = len(coeffs) - 1
    if degree <= 0:
        return []
    lead = coeffs[0]
    coeffs = [c / lead for c in coeffs]
    radius = 1.0 + max(abs(c) for c in coeffs[1:])
    roots = [
        radius
        * complex(
            math.cos(2.0 * math.pi * (i + 0.31) / degree),
            math.sin(2.0 * math.pi * (i + 0.31) / degree),
        )
        for i in range(degree)
    ]
    for _ in range(500):
        max_delta = 0.0
        new_roots = roots[:]
        for i, root in enumerate(roots):
            denom = 1.0 + 0j
            for j, other in enumerate(roots):
                if i != j:
                    denom *= root - other
            if abs(denom) < 1e-24:
                denom += complex(1e-12, -1e-12)
            delta = poly_eval(coeffs, root) / denom
            new_roots[i] = root - delta
            max_delta = max(max_delta, abs(delta))
        roots = new_roots
        if max_delta < 1e-13:
            break
    return sorted(roots, key=lambda z: (abs(z), math.atan2(z.imag, z.real)))


@dataclass(frozen=True)
class Packet:
    name: str
    row: tuple[int, ...]
    note: str


@dataclass(frozen=True)
class BridgeRow:
    carrier: str
    source_hook: str
    hyp3113_role: str
    keeps: str
    destroys: str
    cut_payload: str
    legal_gate: str
    next_measure: str


PACKETS = [
    Packet("consec_8", tuple(range(8)), "HYP-3109 root-locus leader"),
    Packet("even_AP_8", tuple(2 * i for i in range(8)), "dilation-tied AP row"),
    Packet("top_cluster", (0, 7, 8, 9, 10, 11, 12, 13), "high root-clearance control"),
    Packet("break_spread", (0, 1, 5, 7, 8, 9, 11, 13), "root-collision control"),
    Packet("covering_probe", (0, 1, 2, 3, 5, 8, 13, 14), "covering probe"),
    Packet("random_spread", (0, 1, 3, 7, 12, 20, 33, 54), "spread control"),
    Packet("cap_j2", (0, 1, 13), "HYP-3104 cap atom"),
    Packet("cap_j4_dip", (0, 1, 11, 12, 13), "1/4004 dip sidecar"),
    Packet("cap_j5_dip", (0, 1, 5, 7, 8, 9), "1081/76440 dip sidecar"),
    Packet("gw10", (0, 1, 3, 5, 7, 9, 11, 13, 15, 17), "wide binding bank"),
    Packet("gw11", (0, 1, 2, 3, 4, 5, 6, 7, 8, 21, 22), "wide binding bank"),
    Packet("gw12", (0, 2, 4, 6, 7, 8, 10, 11, 12, 14, 18, 20), "wide binding bank"),
    Packet("E_star_k12", (0, 2, 4, 6, 8, 9, 10, 11, 12, 14, 16, 18), "mac-mini k=12 breaker"),
]


BRIDGE_ROWS = [
    BridgeRow(
        carrier="minkowski_lattice_body",
        source_hook="symmetric convex body volume > 2^rank * covolume forces a nonzero lattice point",
        hyp3113_role="memory-lattice-ear map: Bravais_shape_tensor / successive_minima_covolume",
        keeps="rank, covolume, declared q-body inequalities, and forced relation vector",
        destroys="witness time, endpoint owner, root stratum, and route status if the body is only Euclidean",
        cut_payload="q-body inequality word plus forced-vector-to-packet decoding map",
        legal_gate="body must be symmetric, convex, LRC-native, and tied to witness/lift/observer output",
        next_measure="reciprocal-flat q-body using root stratum, segment clearance, entropy, finite-address, observer debt",
    ),
    BridgeRow(
        carrier="proof_circuit_complexity",
        source_hook="Boolean circuit DAG; size is gate count and depth is longest input-output path",
        hyp3113_role="memory-lattice-ear map: packet_sheaf_legal_exit / Savitch_midpoint_certificate",
        keeps="proof obligation DAG, essential inputs, minterm exits, and shortcut regression tests",
        destroys="metric geometry and analytic root motion when a gate is treated as a theorem",
        cut_payload="missing-input bit vector for every proposed shortcut or quotient",
        legal_gate="shortcut must supply a minimal certificate or prove the removed input is reconstructible",
        next_measure="lower-bound style tests for quotients that delete finite_address or observer_gluing fields",
    ),
    BridgeRow(
        carrier="ising_partition_zero",
        source_hook="spin interactions define a Hamiltonian and partition function over configurations",
        hyp3113_role="root-curve map: Lee_Yang_zero_free_region / tournament_Iomega_root_spectrum",
        keeps="interaction graph, coupling sign, partition coefficients, and whole zero multiset",
        destroys="endpoint-owner and branch-route semantics if only a moment or nearest zero is kept",
        cut_payload="spin/proof-obligation incidence word plus partition-zero arc signature",
        legal_gate="couplings must map to named proof obligations and zeros must be compared as a curve",
        next_measure="replace toy carrier graph by HYP-3109 zero-real one-swap component before/after root wall",
    ),
    BridgeRow(
        carrier="demoivre_quintic_fold",
        source_hook="x^5 + 5*a*x^3 + 5*a^2*x + b folds through y^2 + b*y - a^5 = 0",
        hyp3113_role="memory-lattice-ear map: first_obstruction_cocycle / finite branch address",
        keeps="exact fifth-root branch orbit, algebraic cancellation depth, and fold normal form",
        destroys="which LRC packet branch is legal unless finite-address data is attached",
        cut_payload="five-branch orbit word plus finite-address sidecar for the selected branch",
        legal_gate="fold must feed an observer-gluing or finite-address certificate, not a naked algebraic analogy",
        next_measure="formalize Laurent identity and test HYP-3110 branch packets for fifth-root orbit debt",
    ),
]


def minkowski_lattice_report() -> dict[str, object]:
    qmap = {packet.name: missdist(packet.row) for packet in PACKETS}
    coords = {name: q_drop_last(q) for name, q in qmap.items()}
    base_name = "consec_8"
    base = coords[base_name]
    diffs = [
        (name, vsub(vec, base))
        for name, vec in coords.items()
        if name != base_name and vsub(vec, base) != tuple(F(0) for _ in base)
    ]
    diffs.sort(key=lambda item: (norm(item[1]), item[0]))
    basis = independent_basis(diffs)
    rank = len(basis)
    gram = [[dot(a, b) for _, a in basis] for _, b in basis]
    gram_det = determinant_fraction(gram) if rank else F(0)
    covolume = math.sqrt(max(0.0, float(gram_det))) if rank else 0.0
    unit_ball = math.pi ** (rank / 2.0) / math.gamma(rank / 2.0 + 1.0) if rank else 0.0
    force_radius = ((2.0**rank * covolume) / unit_ball) ** (1.0 / rank) if rank else 0.0
    all_diffs: list[tuple[str, str, tuple[F, ...]]] = []
    for a, b in itertools.combinations(sorted(coords), 2):
        diff = vsub(coords[a], coords[b])
        if diff != tuple(F(0) for _ in diff):
            all_diffs.append((a, b, diff))
    shortest = min(all_diffs, key=lambda item: norm(item[2]))
    longest = max(all_diffs, key=lambda item: norm(item[2]))
    radii = {
        "shortest_named_diff": norm(shortest[2]),
        "half_named_diameter": 0.5 * norm(longest[2]),
        "longest_named_diff": norm(longest[2]),
    }
    ratios = {
        key: unit_ball * (value**rank) / ((2.0**rank) * covolume)
        for key, value in radii.items()
    }
    return {
        "qmap": qmap,
        "rank": rank,
        "basis": basis,
        "gram_det": gram_det,
        "covolume": covolume,
        "force_radius": force_radius,
        "shortest": shortest,
        "longest": longest,
        "radii": radii,
        "ratios": ratios,
    }


CIRCUIT_INPUTS = [
    "q_witness_gate",
    "level7_sieve",
    "dyadic_lift",
    "pgf_root_curve",
    "zero_real_ear_map",
    "minkowski_lattice_body",
    "ising_partition_zero",
    "demoivre_fold",
    "observer_gluing_certificate",
    "finite_address_packet",
]

CIRCUIT_GATES = {
    "composite_lift": ("AND", ("level7_sieve", "dyadic_lift")),
    "early_exit": ("OR", ("q_witness_gate", "composite_lift")),
    "analytic_certificate": (
        "AND",
        ("pgf_root_curve", "zero_real_ear_map", "minkowski_lattice_body"),
    ),
    "stat_mech_certificate": (
        "AND",
        ("ising_partition_zero", "pgf_root_curve", "zero_real_ear_map"),
    ),
    "algebraic_certificate": ("AND", ("demoivre_fold", "finite_address_packet")),
    "carrier_certificate": (
        "OR",
        ("analytic_certificate", "stat_mech_certificate", "algebraic_certificate"),
    ),
    "observer_exit": (
        "AND",
        ("observer_gluing_certificate", "finite_address_packet", "carrier_certificate"),
    ),
    "LRC14Statement": ("OR", ("early_exit", "observer_exit")),
}


def circuit_eval(node: str, assignment: dict[str, bool], memo: dict[str, bool]) -> bool:
    if node in memo:
        return memo[node]
    if node in assignment:
        memo[node] = assignment[node]
        return memo[node]
    op, children = CIRCUIT_GATES[node]
    vals = [circuit_eval(child, assignment, memo) for child in children]
    if op == "AND":
        memo[node] = all(vals)
    elif op == "OR":
        memo[node] = any(vals)
    else:
        raise ValueError(op)
    return memo[node]


def circuit_depth(node: str, memo: dict[str, int]) -> int:
    if node in memo:
        return memo[node]
    if node in CIRCUIT_INPUTS:
        memo[node] = 0
    else:
        _, children = CIRCUIT_GATES[node]
        memo[node] = 1 + max(circuit_depth(child, memo) for child in children)
    return memo[node]


def proof_circuit_report() -> dict[str, object]:
    depth = circuit_depth("LRC14Statement", {})
    size = len(CIRCUIT_GATES)
    max_fanin = max(len(children) for _, children in CIRCUIT_GATES.values())
    satisfying_masks: list[int] = []
    for mask in range(1 << len(CIRCUIT_INPUTS)):
        assignment = {
            name: bool(mask & (1 << i)) for i, name in enumerate(CIRCUIT_INPUTS)
        }
        if circuit_eval("LRC14Statement", assignment, {}):
            satisfying_masks.append(mask)
    minterms: list[int] = []
    sat_set = set(satisfying_masks)
    for mask in satisfying_masks:
        minimal = True
        sub = (mask - 1) & mask
        while sub:
            if sub in sat_set:
                minimal = False
                break
            sub = (sub - 1) & mask
        if minimal:
            minterms.append(mask)
    essential = []
    for i, name in enumerate(CIRCUIT_INPUTS):
        bit = 1 << i
        found = False
        for mask in range(1 << len(CIRCUIT_INPUTS)):
            a = {n: bool(mask & (1 << j)) for j, n in enumerate(CIRCUIT_INPUTS)}
            bmask = mask ^ bit
            b = {n: bool(bmask & (1 << j)) for j, n in enumerate(CIRCUIT_INPUTS)}
            if circuit_eval("LRC14Statement", a, {}) != circuit_eval("LRC14Statement", b, {}):
                found = True
                break
        if found:
            essential.append(name)
    minterm_names = [
        tuple(CIRCUIT_INPUTS[i] for i in range(len(CIRCUIT_INPUTS)) if mask & (1 << i))
        for mask in minterms
    ]
    return {
        "size": size,
        "depth": depth,
        "max_fanin": max_fanin,
        "satisfying_count": len(satisfying_masks),
        "minterms": minterm_names,
        "minterm_size_hist": Counter(len(m) for m in minterm_names),
        "essential_inputs": essential,
    }


def ising_coefficients(n: int, edges: tuple[tuple[int, int], ...], r: int = 2) -> list[int]:
    coeffs = [0 for _ in range(n + 1)]
    for mask in range(1 << n):
        down = mask.bit_count()
        weight = 1
        for a, b in edges:
            same = bool(mask & (1 << a)) == bool(mask & (1 << b))
            if same:
                weight *= r
        coeffs[down] += weight
    return coeffs


ISING_GRAPHS = {
    "path6": (6, tuple((i, i + 1) for i in range(5))),
    "cycle6": (6, tuple((i, (i + 1) % 6) for i in range(6))),
    "carrier7": (
        7,
        (
            (0, 1),  # pgf/root
            (0, 2),  # pgf/minkowski
            (0, 3),  # pgf/ising
            (0, 5),  # pgf/observer
            (1, 3),  # root/ising
            (1, 5),  # root/observer
            (2, 6),  # minkowski/finite address
            (3, 5),  # ising/observer
            (4, 5),  # demoivre/observer
            (4, 6),  # demoivre/finite address
            (5, 6),  # observer/finite address
        ),
    ),
}


def ising_report() -> list[dict[str, object]]:
    rows = []
    for name, (n, edges) in ISING_GRAPHS.items():
        coeffs = ising_coefficients(n, edges, r=2)
        roots = durand_kerner([float(c) for c in reversed(coeffs)])
        max_radial_error = max(abs(abs(z) - 1.0) for z in roots) if roots else 0.0
        angles = sorted(math.degrees(math.atan2(z.imag, z.real)) for z in roots)
        nearest_to_apex = min(abs(abs(a) - 360.0 / 7.0) for a in angles) if angles else 0.0
        rows.append(
            {
                "name": name,
                "n": n,
                "edges": len(edges),
                "coeffs": coeffs,
                "max_radial_error": max_radial_error,
                "angles": angles,
                "nearest_to_apex": nearest_to_apex,
            }
        )
    return rows


Laurent = dict[tuple[int, int], int]  # (u_power, a_power) -> coefficient


def l_add(a: Laurent, b: Laurent) -> Laurent:
    out = defaultdict(int)
    for key, value in itertools.chain(a.items(), b.items()):
        out[key] += value
    return {key: value for key, value in out.items() if value}


def l_mul(a: Laurent, b: Laurent) -> Laurent:
    out = defaultdict(int)
    for (up1, ap1), c1 in a.items():
        for (up2, ap2), c2 in b.items():
            out[(up1 + up2, ap1 + ap2)] += c1 * c2
    return {key: value for key, value in out.items() if value}


def l_scale_shift(a: Laurent, coeff: int, a_power: int) -> Laurent:
    return {
        (up, ap + a_power): coeff * value
        for (up, ap), value in a.items()
        if coeff * value
    }


def l_power(a: Laurent, k: int) -> Laurent:
    out: Laurent = {(0, 0): 1}
    for _ in range(k):
        out = l_mul(out, a)
    return out


def demoivre_report() -> dict[str, object]:
    # x = u - a/u.  Verify:
    # x^5 + 5*a*x^3 + 5*a^2*x = u^5 - a^5/u^5.
    x: Laurent = {(1, 0): 1, (-1, 1): -1}
    lhs = l_add(
        l_add(l_power(x, 5), l_scale_shift(l_power(x, 3), 5, 1)),
        l_scale_shift(x, 5, 2),
    )
    expected: Laurent = {(5, 0): 1, (-5, 5): -1}
    residual = l_add(lhs, {key: -value for key, value in expected.items()})
    return {
        "identity_ok": residual == {},
        "lhs": lhs,
        "expected": expected,
        "residual": residual,
        "fold_depth": "quintic -> auxiliary quadratic in y=u^5 -> five root-of-unity branches",
        "branch_orbit_size": 5,
    }


TOURNAMENT_VERTICES = [
    "observer_gluing_certificate",
    "finite_address_packet",
    "lee_yang_root_curve",
    "zero_real_ear_map",
    "proof_circuit_complexity",
    "minkowski_lattice_body",
    "ising_partition_zero",
    "demoivre_quintic_fold",
    "raw_scalar_p0",
]

TOURNAMENT_AXES = {
    "preserves_lrc_predicate": {
        "observer_gluing_certificate": 5,
        "finite_address_packet": 5,
        "lee_yang_root_curve": 4,
        "zero_real_ear_map": 4,
        "proof_circuit_complexity": 3,
        "minkowski_lattice_body": 3,
        "ising_partition_zero": 3,
        "demoivre_quintic_fold": 3,
        "raw_scalar_p0": 1,
    },
    "controls_destroyed_coordinate": {
        "finite_address_packet": 5,
        "observer_gluing_certificate": 5,
        "proof_circuit_complexity": 4,
        "zero_real_ear_map": 4,
        "demoivre_quintic_fold": 3,
        "lee_yang_root_curve": 3,
        "minkowski_lattice_body": 3,
        "ising_partition_zero": 2,
        "raw_scalar_p0": 0,
    },
    "finite_checkable": {
        "demoivre_quintic_fold": 5,
        "proof_circuit_complexity": 5,
        "ising_partition_zero": 5,
        "minkowski_lattice_body": 4,
        "zero_real_ear_map": 4,
        "finite_address_packet": 4,
        "lee_yang_root_curve": 3,
        "observer_gluing_certificate": 3,
        "raw_scalar_p0": 5,
    },
    "lean_facing": {
        "demoivre_quintic_fold": 5,
        "proof_circuit_complexity": 4,
        "finite_address_packet": 4,
        "observer_gluing_certificate": 4,
        "minkowski_lattice_body": 3,
        "zero_real_ear_map": 3,
        "lee_yang_root_curve": 2,
        "ising_partition_zero": 2,
        "raw_scalar_p0": 2,
    },
    "analytic_strength": {
        "lee_yang_root_curve": 5,
        "ising_partition_zero": 4,
        "minkowski_lattice_body": 4,
        "zero_real_ear_map": 4,
        "observer_gluing_certificate": 3,
        "finite_address_packet": 2,
        "proof_circuit_complexity": 2,
        "demoivre_quintic_fold": 1,
        "raw_scalar_p0": 1,
    },
    "algebraic_fold_exactness": {
        "demoivre_quintic_fold": 5,
        "finite_address_packet": 4,
        "proof_circuit_complexity": 4,
        "minkowski_lattice_body": 3,
        "ising_partition_zero": 3,
        "lee_yang_root_curve": 3,
        "zero_real_ear_map": 3,
        "observer_gluing_certificate": 2,
        "raw_scalar_p0": 1,
    },
    "handoff_legality": {
        "observer_gluing_certificate": 5,
        "finite_address_packet": 5,
        "proof_circuit_complexity": 4,
        "zero_real_ear_map": 3,
        "lee_yang_root_curve": 3,
        "minkowski_lattice_body": 3,
        "ising_partition_zero": 2,
        "demoivre_quintic_fold": 2,
        "raw_scalar_p0": 0,
    },
}

TIE_ORDER = {
    name: i
    for i, name in enumerate(
        [
            "observer_gluing_certificate",
            "finite_address_packet",
            "lee_yang_root_curve",
            "zero_real_ear_map",
            "proof_circuit_complexity",
            "minkowski_lattice_body",
            "ising_partition_zero",
            "demoivre_quintic_fold",
            "raw_scalar_p0",
        ]
    )
}


def orient(a: str, b: str) -> tuple[str, int, int]:
    awins = 0
    bwins = 0
    for scores in TOURNAMENT_AXES.values():
        if scores[a] > scores[b]:
            awins += 1
        elif scores[b] > scores[a]:
            bwins += 1
    if awins > bwins:
        return a, awins, bwins
    if bwins > awins:
        return b, bwins, awins
    winner = a if TIE_ORDER[a] < TIE_ORDER[b] else b
    return winner, awins, bwins


def sccs(vertices: list[str], edges: dict[str, set[str]]) -> list[list[str]]:
    index = 0
    stack: list[str] = []
    onstack: set[str] = set()
    indices: dict[str, int] = {}
    low: dict[str, int] = {}
    comps: list[list[str]] = []

    def visit(v: str) -> None:
        nonlocal index
        indices[v] = index
        low[v] = index
        index += 1
        stack.append(v)
        onstack.add(v)
        for w in edges[v]:
            if w not in indices:
                visit(w)
                low[v] = min(low[v], low[w])
            elif w in onstack:
                low[v] = min(low[v], indices[w])
        if low[v] == indices[v]:
            comp = []
            while True:
                w = stack.pop()
                onstack.remove(w)
                comp.append(w)
                if w == v:
                    break
            comps.append(sorted(comp, key=TIE_ORDER.get))

    for vertex in vertices:
        if vertex not in indices:
            visit(vertex)
    return comps


def hamiltonian_path_count(vertices: list[str], edges: dict[str, set[str]]) -> int:
    n = len(vertices)
    idx = {v: i for i, v in enumerate(vertices)}
    dp = [[0 for _ in vertices] for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for mask in range(1 << n):
        for last in range(n):
            count = dp[mask][last]
            if not count:
                continue
            last_v = vertices[last]
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                if vertices[nxt] in edges[last_v]:
                    dp[mask | (1 << nxt)][nxt] += count
    full = (1 << n) - 1
    return sum(dp[full][i] for i in range(n))


def tournament_report() -> dict[str, object]:
    edges = {v: set() for v in TOURNAMENT_VERTICES}
    margins = []
    for a, b in itertools.combinations(TOURNAMENT_VERTICES, 2):
        winner, win_axes, lose_axes = orient(a, b)
        loser = b if winner == a else a
        edges[winner].add(loser)
        margins.append((win_axes - lose_axes, winner, loser, win_axes, lose_axes))
    score_hist = Counter(len(edges[v]) for v in TOURNAMENT_VERTICES)
    tri_cycles = []
    for triple in itertools.combinations(TOURNAMENT_VERTICES, 3):
        a, b, c = triple
        if (
            (b in edges[a] and c in edges[b] and a in edges[c])
            or (c in edges[a] and b in edges[c] and a in edges[b])
        ):
            tri_cycles.append(triple)
    return {
        "scores": {v: len(edges[v]) for v in TOURNAMENT_VERTICES},
        "score_hist": dict(sorted(score_hist.items())),
        "directed_3cycles": len(tri_cycles),
        "sample_3cycles": tri_cycles[:6],
        "sccs": sccs(TOURNAMENT_VERTICES, edges),
        "hamiltonian_path_count": hamiltonian_path_count(TOURNAMENT_VERTICES, edges),
        "edge_flip_risks": sorted(margins)[:12],
        "priority_path": sorted(TOURNAMENT_VERTICES, key=lambda v: -len(edges[v])),
    }


def print_minkowski(report: dict[str, object]) -> None:
    print("=" * 96)
    print("MAP 1: MINKOWSKI q-LATTICE BODY SIDEcar")
    print("=" * 96)
    print("q-vector rows used:")
    qmap: dict[str, tuple[F, ...]] = report["qmap"]  # type: ignore[assignment]
    for packet in PACKETS:
        q = qmap[packet.name]
        print(f"  {packet.name:16s} q0={float(q[0]):.5f} q={tuple(fstr(x) for x in q)}")
    print()
    print(f"affine q-rank from named packets: {report['rank']} (max 6)")
    print(f"basis from consec_8 differences:")
    for name, vec in report["basis"]:  # type: ignore[assignment]
        print(f"  {name:16s} norm={norm(vec):.6f}")
    print(f"Gram determinant: {report['gram_det']}")
    print(f"covolume proxy sqrt(det Gram): {report['covolume']:.10g}")
    print(f"Euclidean ball radius with volume = 2^rank*covolume: {report['force_radius']:.6f}")
    shortest = report["shortest"]  # type: ignore[assignment]
    longest = report["longest"]  # type: ignore[assignment]
    print(
        "shortest nonzero named q-difference: "
        f"{shortest[0]} <-> {shortest[1]} norm={norm(shortest[2]):.6f}"
    )
    print(
        "longest named q-difference: "
        f"{longest[0]} <-> {longest[1]} norm={norm(longest[2]):.6f}"
    )
    print("Minkowski ratios vol(B_R)/(2^rank*covolume):")
    for key, value in report["ratios"].items():  # type: ignore[union-attr]
        print(f"  {key:22s} R={report['radii'][key]:.6f} ratio={value:.6f}")  # type: ignore[index]
    print(
        "readout: Minkowski supplies a forcing threshold for a declared symmetric "
        "convex q-body; it does not by itself identify an LRC witness time."
    )


def print_circuit(report: dict[str, object]) -> None:
    print()
    print("=" * 96)
    print("MAP 2: PROOF-STATE CIRCUIT COMPLEXITY SIDEcar")
    print("=" * 96)
    print(f"inputs={len(CIRCUIT_INPUTS)} gates(size)={report['size']} depth={report['depth']} max_fanin={report['max_fanin']}")
    print(f"truth-table satisfying assignments: {report['satisfying_count']} / {1 << len(CIRCUIT_INPUTS)}")
    print(f"essential inputs ({len(report['essential_inputs'])}): {', '.join(report['essential_inputs'])}")
    print(f"minimal certificate size histogram: {dict(sorted(report['minterm_size_hist'].items()))}")
    print("minimal monotone certificates:")
    for term in report["minterms"]:  # type: ignore[assignment]
        print("  " + " AND ".join(term))
    print(
        "readout: the current proof frontier has a small monotone OR of exits, "
        "but observer_exit still needs finite_address_packet plus one carrier certificate."
    )


def print_ising(rows: list[dict[str, object]]) -> None:
    print()
    print("=" * 96)
    print("MAP 3: FINITE ISING PARTITION-ZERO SIDEcar")
    print("=" * 96)
    for row in rows:
        angles = row["angles"]  # type: ignore[assignment]
        upper = [a for a in angles if a >= 0]
        print(f"{row['name']:10s} spins={row['n']} edges={row['edges']} coeffs={row['coeffs']}")
        print(
            f"  max ||z|-1|={row['max_radial_error']:.3e} "
            f"nearest |angle|-360/7={row['nearest_to_apex']:.3f} deg "
            f"upper_angles={[round(a, 2) for a in upper]}"
        )
    print(
        "readout: these finite ferromagnetic packets keep zeros numerically on "
        "the unit circle; the useful LRC sidecar is the whole zero set, not one moment."
    )


def monomial_text(key: tuple[int, int], coeff: int) -> str:
    up, ap = key
    factors = []
    if ap:
        factors.append("a" if ap == 1 else f"a^{ap}")
    if up:
        if up == 1:
            factors.append("u")
        elif up == -1:
            factors.append("u^-1")
        else:
            factors.append(f"u^{up}")
    if not factors:
        factors.append("1")
    body = "*".join(factors)
    if coeff == 1:
        return body
    if coeff == -1:
        return "-" + body
    return f"{coeff}*{body}"


def print_demoivre(report: dict[str, object]) -> None:
    print()
    print("=" * 96)
    print("MAP 4: DE MOIVRE QUINTIC FOLD SIDEcar")
    print("=" * 96)
    lhs = report["lhs"]  # type: ignore[assignment]
    ordered = sorted(lhs.items(), key=lambda item: (-item[0][0], item[0][1]))
    print(f"identity_ok={report['identity_ok']}")
    print("x = u - a/u")
    print("x^5 + 5*a*x^3 + 5*a^2*x = " + " + ".join(monomial_text(k, c) for k, c in ordered).replace("+ -", "- "))
    print("If y=u^5 and y^2 + b*y - a^5 = 0, then y - a^5/y = -b, so")
    print("  x^5 + 5*a*x^3 + 5*a^2*x + b = 0.")
    print(f"fold_depth={report['fold_depth']}")
    print(f"branch_orbit_size={report['branch_orbit_size']}")
    print(
        "readout: De Moivre is the cleanest Lean-facing algebraic fold, but it "
        "preserves only a cancellation normal form unless paired with finite-address data."
    )


def print_tournament(report: dict[str, object]) -> None:
    print()
    print("=" * 96)
    print("MAP 5: TOURNAMENT ANALYSIS ON PROOF-CARRIER SIDECARS")
    print("=" * 96)
    print("vertices are proof carriers/obligations, not runners or arcs.")
    print("axes=" + ", ".join(TOURNAMENT_AXES))
    print(f"score_hist={report['score_hist']}")
    print(f"scores={report['scores']}")
    print(f"directed_3cycles={report['directed_3cycles']}")
    print(f"sample_3cycles={report['sample_3cycles']}")
    print(f"SCCs={report['sccs']}")
    print(f"Hamiltonian_path_count={report['hamiltonian_path_count']}")
    print("low-margin edge-flip risks:")
    for margin, winner, loser, win_axes, lose_axes in report["edge_flip_risks"]:  # type: ignore[assignment]
        print(f"  margin={margin}: {winner} -> {loser} ({win_axes}-{lose_axes} axes)")
    print("one priority path: " + " -> ".join(report["priority_path"]))  # type: ignore[arg-type]


def print_bridge_rows() -> None:
    print()
    print("=" * 96)
    print("MAP 6: HYP-3113 TWO-MAP CUT-PAYLOAD BRIDGE")
    print("=" * 96)
    print(
        "The four imported carriers are retention tests for HYP-3113's "
        "root-curve portfolio and memory-lattice-ear certificate ladder."
    )
    print(
        "duodecimal audit: 4 carriers x 3 legal cells "
        "(predicate, destroyed coordinate, handoff payload) = 12 named cells; "
        "proof-closed carriers = 0 / 4."
    )
    for row in BRIDGE_ROWS:
        print()
        print(f"- {row.carrier}")
        print(f"  source_hook: {row.source_hook}")
        print(f"  HYP-3113 role: {row.hyp3113_role}")
        print(f"  keeps: {row.keeps}")
        print(f"  destroys: {row.destroys}")
        print(f"  cut_payload: {row.cut_payload}")
        print(f"  legal_gate: {row.legal_gate}")
        print(f"  next_measure: {row.next_measure}")


def main() -> None:
    print("HYP-3111 Minkowski/circuit/Ising/De Moivre carrier scout -- codex S264")
    print("This scout integrates HYP-3108/HYP-3109/HYP-3110 with the user-requested sources.")
    print()
    print("External source hooks used:")
    print("  Minkowski theorem: symmetric convex body volume threshold 2^n * covolume.")
    print("  Circuit complexity: circuit size = gate count; depth = longest input-output path.")
    print("  Lee-Yang/Ising: finite ferromagnetic partition zeros are the whole zero-set sidecar.")
    print("  De Moivre quintic: x^5 + 5*a*x^3 + 5*a^2*x + b = 0 folds through y^2+b*y-a^5.")
    print()
    minkowski = minkowski_lattice_report()
    circuit = proof_circuit_report()
    ising = ising_report()
    demoivre = demoivre_report()
    tournament = tournament_report()
    print_minkowski(minkowski)
    print_circuit(circuit)
    print_ising(ising)
    print_demoivre(demoivre)
    print_tournament(tournament)
    print_bridge_rows()
    print()
    print("=" * 96)
    print("SYNTHESIS")
    print("=" * 96)
    print("strengthened:")
    print("  Minkowski converts the Bravais/q-rank story into an explicit convex-body threshold.")
    print("  Circuit complexity exposes the current proof frontier as a shallow monotone exit DAG.")
    print("  The Ising finite-zero check supports keeping the entire Lee-Yang zero geometry.")
    print("  De Moivre gives an exact algebraic fold that is plausibly formalizable.")
    print("still open:")
    print("  None of the four carriers proves the residual LRC14 observer certificate alone.")
    print("  The live bridge remains: finite-address packet + observer gluing + root/ear sidecar.")
    print("  HYP-3113 reframes this as packet-sheaf legal exit after a root-curve or memory-lattice-ear sidecar.")
    print("challenged assumption:")
    print("  The tournament vertices need not be runners, gaps, arcs, or roots; in this scout")
    print("  they are proof carriers.  This preserves route legality and destroys raw labels/time.")
    print("DONE.")


if __name__ == "__main__":
    main()
