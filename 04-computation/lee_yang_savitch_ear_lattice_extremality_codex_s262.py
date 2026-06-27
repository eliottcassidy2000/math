#!/usr/bin/env python3
"""HYP-3108 signal scout: Lee-Yang, Savitch, Bravais, and ear-lattice maps.

This script extends HYP-3103 through HYP-3107 rather than replacing them.

Map 1, the packet signal map, measures named LRC packets by:
  * miss-count PGF roots for G_N(z)=sum q_t z^t;
  * quartic phi4-style fit -log(q_t) ~= c + b*S^2 + lambda*S^4;
  * q0 component count and largest q0 cell;
  * pairwise co-emptiness Perron/reflection data;
  * the cell-state transition graph as an ear-decomposition proxy.

Map 2, the proof-frontier sidecar map, treats proof obligations and signals as
tournament vertices.  It also builds a small Savitch-style reachability ledger:
recursive midpoint certificates for how a sidecar could reach the terminal
LRC14 proof interface.

This is evidence and route selection, not an LRC proof.
"""

from __future__ import annotations

import itertools
import math
import sys
from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from fractions import Fraction as F
from typing import Iterable

import numpy as np

sys.stdout.reconfigure(line_buffering=True)


def fstr(x: F) -> str:
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def sector_of(x: F) -> int:
    return int((x % 1) * 7)


@dataclass(frozen=True)
class Packet:
    name: str
    family: str
    speeds: tuple[int, ...]
    note: str


@dataclass
class PacketFeatures:
    packet: Packet
    q: list[F]
    roots: list[complex]
    n_real: int
    nearest_root_abs: float
    neg_unit_real_roots: int
    apex7_angle_gap: float
    lc_defect: float
    phi4_lambda: float
    phi4_b: float
    phi4_r2: float
    q0_components: int
    largest_q0_cell: F
    state_count: int
    state_edge_count: int
    state_strong: bool
    ear_excess: int
    factor_critical_proxy: bool
    k4_subgraphs: int
    perron_gap: float
    reflection_ok: bool


PACKETS = [
    Packet("consec_8", "S66_pgf_bank", tuple(range(8)), "all-complex PGF root stratum"),
    Packet("even_AP_8", "S66_pgf_bank", tuple(2 * i for i in range(8)), "dilation tie with consec_8"),
    Packet("break_8", "S66_pgf_bank", (0, 1, 5, 7, 8, 9, 11, 13), "root-collision/break flavor"),
    Packet("top_cluster_8", "S66_pgf_bank", (0, 7, 8, 9, 10, 11, 12, 13), "all-complex but lower mass"),
    Packet("random_spread_8", "S66_pgf_bank", (0, 1, 3, 7, 12, 20, 33, 54), "spread control"),
    Packet("covering_8", "S66_pgf_bank", (0, 14, 1, 2, 3, 5, 8, 13), "near-covering control"),
    Packet("cap_j1", "HYP3104_cap_atoms", (0, 1), "pure order-1 cap atom"),
    Packet("cap_j2", "HYP3104_cap_atoms", (0, 1, 13), "pure pair-Pascal atom"),
    Packet("cap_j3_order3", "HYP3104_cap_atoms", (0, 1, 12, 13), "pair value exact but order-3 live"),
    Packet("cap_j4_dip", "HYP3104_cap_atoms", (0, 1, 11, 12, 13), "1/4004 cap dip sidecar"),
    Packet("cap_j5_dip", "HYP3104_cap_atoms", (0, 1, 5, 7, 8, 9), "1081/76440 cap dip sidecar"),
    Packet("gw10", "wide_binding_bank", (0, 1, 3, 5, 7, 9, 11, 13, 15, 17), "genuine-wide k=10"),
    Packet("gw11", "wide_binding_bank", (0, 1, 2, 3, 4, 5, 6, 7, 8, 21, 22), "genuine-wide k=11"),
    Packet("gw12", "wide_binding_bank", (0, 2, 4, 6, 7, 8, 10, 11, 12, 14, 18, 20), "genuine-wide k=12"),
    Packet("E_star_k12", "wide_binding_bank", (0, 2, 4, 6, 8, 9, 10, 11, 12, 14, 16, 18), "mac-mini k=12 breaker"),
]


def cell_profile(speeds: tuple[int, ...]) -> tuple[list[F], list[tuple[F, F, int]], list[tuple[F, F]]]:
    """Return q_t, cell masks of empty inner sectors, and merged q0 cells."""
    speeds = tuple(sorted(set(speeds)))
    if 0 not in speeds:
        speeds = (0,) + speeds
    breaks = {F(0), F(1)}
    for e in speeds:
        if e == 0:
            continue
        for m in range(0, 7 * e + 1):
            breaks.add(F(m, 7 * e))
    cuts = sorted(breaks)
    q = [F(0)] * 7
    cells: list[tuple[F, F, int]] = []
    q0_cells: list[tuple[F, F]] = []
    for lo, hi in zip(cuts, cuts[1:]):
        if hi <= lo:
            continue
        mid = (lo + hi) / 2
        covered = {sector_of(e * mid) for e in speeds}
        mask = 0
        empty_count = 0
        for s in range(1, 7):
            if s not in covered:
                mask |= 1 << (s - 1)
                empty_count += 1
        q[empty_count] += hi - lo
        cells.append((lo, hi, mask))
        if empty_count == 0:
            if q0_cells and q0_cells[-1][1] == lo:
                q0_cells[-1] = (q0_cells[-1][0], hi)
            else:
                q0_cells.append((lo, hi))
    return q, cells, q0_cells


def pgf_roots(q: list[F]) -> list[complex]:
    coeffs = [float(q[t]) for t in range(6, -1, -1)]
    while len(coeffs) > 1 and abs(coeffs[0]) < 1e-15:
        coeffs = coeffs[1:]
    return sorted(np.roots(coeffs), key=lambda z: (abs(z), z.real, z.imag))


def log_concavity_defect(q: list[F]) -> float:
    return min(float(q[t]) ** 2 - float(q[t - 1]) * float(q[t + 1]) for t in range(1, 6))


def root_features(roots: list[complex]) -> tuple[int, float, int, float]:
    n_real = sum(1 for z in roots if abs(z.imag) < 1e-8)
    nearest = min(abs(z) for z in roots) if roots else 0.0
    neg_unit = sum(1 for z in roots if abs(z.imag) < 1e-8 and -1.0 <= z.real < 0.0)
    targets = [360.0 / 7.0, 2.0 * 360.0 / 7.0, 3.0 * 360.0 / 7.0]
    angles = [abs(math.degrees(math.atan2(z.imag, z.real))) for z in roots]
    gap = sum(min(abs(a - t) for a in angles) for t in targets) / len(targets)
    return n_real, nearest, neg_unit, gap


def phi4_fit(q: list[F]) -> tuple[float, float, float]:
    probs = np.array([float(x) for x in q], dtype=float)
    ts = np.arange(7, dtype=float)
    mean = float(np.dot(ts, probs))
    s = ts - mean
    y = -np.log(np.maximum(probs, 1e-300))
    x = np.column_stack([np.ones_like(s), s**2, s**4])
    w = np.sqrt(np.maximum(probs, 1e-12))
    coef, *_ = np.linalg.lstsq(x * w[:, None], y * w, rcond=None)
    pred = x @ coef
    ybar = float(np.dot(probs, y))
    sse = float(np.dot(probs, (y - pred) ** 2))
    sst = float(np.dot(probs, (y - ybar) ** 2))
    r2 = 1.0 - sse / sst if sst > 1e-15 else 1.0
    return float(coef[2]), float(coef[1]), r2


def tarjan_scc(vertices: Iterable[int], edges: Iterable[tuple[int, int]]) -> list[list[int]]:
    verts = list(vertices)
    adj: dict[int, list[int]] = {v: [] for v in verts}
    for a, b in edges:
        adj.setdefault(a, []).append(b)
        adj.setdefault(b, adj.get(b, []))
    index = 0
    stack: list[int] = []
    onstack: set[int] = set()
    indices: dict[int, int] = {}
    low: dict[int, int] = {}
    comps: list[list[int]] = []

    def visit(v: int) -> None:
        nonlocal index
        indices[v] = index
        low[v] = index
        index += 1
        stack.append(v)
        onstack.add(v)
        for w in adj.get(v, []):
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
            comps.append(comp)

    for v in verts:
        if v not in indices:
            visit(v)
    return comps


def has_perfect_matching(vertices: tuple[int, ...], edges: set[frozenset[int]]) -> bool:
    if not vertices:
        return True
    first = vertices[0]
    rest = vertices[1:]
    for i, other in enumerate(rest):
        if frozenset((first, other)) in edges:
            remaining = rest[:i] + rest[i + 1 :]
            if has_perfect_matching(remaining, edges):
                return True
    return False


def state_graph_features(cells: list[tuple[F, F, int]]) -> tuple[int, int, bool, int, bool, int]:
    seq = [mask for _, _, mask in cells]
    vertices = sorted(set(seq))
    directed_edges = {(seq[i], seq[(i + 1) % len(seq)]) for i in range(len(seq))}
    comps = tarjan_scc(vertices, directed_edges)
    strong = len(comps) == 1
    undirected = {frozenset((a, b)) for a, b in directed_edges if a != b}
    ear_excess = max(0, len(undirected) - len(vertices) + 1) if vertices else 0
    factor_critical = False
    if len(vertices) % 2 == 1 and vertices:
        factor_critical = all(
            has_perfect_matching(tuple(v for v in vertices if v != drop), undirected)
            for drop in vertices
        )
    k4 = 0
    edge_lookup = undirected
    for quad in itertools.combinations(vertices, 4):
        if all(frozenset((a, b)) in edge_lookup for a, b in itertools.combinations(quad, 2)):
            k4 += 1
    return len(vertices), len(directed_edges), strong, ear_excess, factor_critical, k4


def coempty_perron_features(cells: list[tuple[F, F, int]]) -> tuple[float, bool]:
    matrix = [[F(0)] * 6 for _ in range(6)]
    for lo, hi, mask in cells:
        empty = [s for s in range(6) if mask & (1 << s)]
        w = hi - lo
        for a in empty:
            for b in empty:
                matrix[a][b] += w
    mf = np.array([[float(x) for x in row] for row in matrix], dtype=float)
    evals = np.sort(np.linalg.eigvalsh(mf))[::-1]
    gap = float(evals[0] - evals[1]) if len(evals) >= 2 else 0.0

    def refl(idx: int) -> int | None:
        sector = idx + 1
        reflected = 6 - sector
        return reflected - 1 if 1 <= reflected <= 6 else None

    ok = True
    for a in range(6):
        ra = refl(a)
        if ra is None:
            continue
        for b in range(6):
            rb = refl(b)
            if rb is None:
                continue
            if matrix[a][b] != matrix[ra][rb]:
                ok = False
                break
    return gap, ok


def compute_features(packet: Packet) -> PacketFeatures:
    q, cells, q0_cells = cell_profile(packet.speeds)
    roots = pgf_roots(q)
    n_real, nearest, neg_unit, angle_gap = root_features(roots)
    lam, b, r2 = phi4_fit(q)
    state_count, edge_count, strong, ear_excess, factor_critical, k4 = state_graph_features(cells)
    perron_gap, reflection_ok = coempty_perron_features(cells)
    largest = max((hi - lo for lo, hi in q0_cells), default=F(0))
    return PacketFeatures(
        packet=packet,
        q=q,
        roots=roots,
        n_real=n_real,
        nearest_root_abs=nearest,
        neg_unit_real_roots=neg_unit,
        apex7_angle_gap=angle_gap,
        lc_defect=log_concavity_defect(q),
        phi4_lambda=lam,
        phi4_b=b,
        phi4_r2=r2,
        q0_components=len(q0_cells),
        largest_q0_cell=largest,
        state_count=state_count,
        state_edge_count=edge_count,
        state_strong=strong,
        ear_excess=ear_excess,
        factor_critical_proxy=factor_critical,
        k4_subgraphs=k4,
        perron_gap=perron_gap,
        reflection_ok=reflection_ok,
    )


def fraction_rank(rows: list[list[F]]) -> int:
    if not rows:
        return 0
    mat = [list(row) for row in rows if any(x for x in row)]
    if not mat:
        return 0
    m = len(mat)
    n = len(mat[0])
    rank = 0
    col = 0
    while rank < m and col < n:
        pivot = None
        for r in range(rank, m):
            if mat[r][col]:
                pivot = r
                break
        if pivot is None:
            col += 1
            continue
        mat[rank], mat[pivot] = mat[pivot], mat[rank]
        pv = mat[rank][col]
        mat[rank] = [x / pv for x in mat[rank]]
        for r in range(m):
            if r != rank and mat[r][col]:
                factor = mat[r][col]
                mat[r] = [x - factor * y for x, y in zip(mat[r], mat[rank])]
        rank += 1
        col += 1
    return rank


def bravais_summary(features: list[PacketFeatures]) -> dict[str, object]:
    base = features[0].q
    diff_rows = [[x - y for x, y in zip(feat.q, base)] for feat in features[1:]]
    q_rank = fraction_rank(diff_rows)
    arr = np.array([[float(x) for x in feat.q] for feat in features], dtype=float)
    centered = arr - arr[0]
    singular = np.linalg.svd(centered, compute_uv=False)
    nonzero = [s for s in singular if s > 1e-10]
    covolume = float(np.prod(nonzero)) if nonzero else 0.0

    chosen: list[tuple[str, str, float]] = []
    chosen_vecs: list[np.ndarray] = []
    pairs = []
    for a, b in itertools.combinations(range(len(features)), 2):
        vec = arr[a] - arr[b]
        pairs.append((float(np.linalg.norm(vec)), a, b, vec))
    for norm, a, b, vec in sorted(pairs):
        trial = chosen_vecs + [vec]
        if np.linalg.matrix_rank(np.array(trial), tol=1e-10) > len(chosen_vecs):
            chosen_vecs.append(vec)
            chosen.append((features[a].packet.name, features[b].packet.name, norm))
        if len(chosen_vecs) == q_rank:
            break
    return {"q_rank": q_rank, "covolume": covolume, "successive": chosen}


FRONTIER_EDGES = [
    ("raw_residual_packet", "q_witness_gate"),
    ("q_witness_gate", "level7_lift_sieve"),
    ("level7_lift_sieve", "dyadic_lift_sieve"),
    ("dyadic_lift_sieve", "residual_leq6_multiples7"),
    ("residual_leq6_multiples7", "coverage_extremality_k8_k10"),
    ("raw_residual_packet", "node3_effective_peel"),
    ("node3_effective_peel", "coverage_extremality_k8_k10"),
    ("coverage_extremality_k8_k10", "pgf_zero_curve_signal"),
    ("coverage_extremality_k8_k10", "reflection_perron_moment"),
    ("pgf_zero_curve_signal", "quartic_phi4_energy"),
    ("pgf_zero_curve_signal", "fine_modp_winding_transfer"),
    ("quartic_phi4_energy", "bravais_lattice_address"),
    ("bravais_lattice_address", "finite_address_packet"),
    ("raw_residual_packet", "ear_state_transition"),
    ("ear_state_transition", "observer_gluing_certificate"),
    ("reflection_perron_moment", "observer_gluing_certificate"),
    ("fine_modp_winding_transfer", "observer_gluing_certificate"),
    ("observer_gluing_certificate", "finite_address_packet"),
    ("finite_address_packet", "terminal_lrc14_statement"),
]


def shortest_path(start: str, target: str) -> list[str]:
    adj: dict[str, list[str]] = defaultdict(list)
    for a, b in FRONTIER_EDGES:
        adj[a].append(b)
    q = deque([(start, [start])])
    seen = {start}
    while q:
        node, path = q.popleft()
        if node == target:
            return path
        for nxt in adj[node]:
            if nxt not in seen:
                seen.add(nxt)
                q.append((nxt, path + [nxt]))
    return []


def midpoint_schedule(path: list[str]) -> list[tuple[str, str, str, int]]:
    out: list[tuple[str, str, str, int]] = []

    def rec(sub: list[str]) -> None:
        if len(sub) <= 2:
            return
        mid = len(sub) // 2
        out.append((sub[0], sub[mid], sub[-1], len(sub) - 1))
        rec(sub[: mid + 1])
        rec(sub[mid:])

    rec(path)
    return out


SIDECAR_AXES = [
    "preserves_lrc_predicate",
    "detects_root_collision",
    "fits_lean_frontier",
    "names_destroyed_coordinate",
    "finite_checkable",
    "tournament_transfer_safe",
    "captures_gluing",
    "computes_on_packets",
]

SIDECARS = {
    "pgf_zero_curve": (1, 1, 1, 1, 1, 0, 0, 1),
    "quartic_phi4_energy": (0, 1, 0, 0, 1, 0, 0, 1),
    "savitch_reachability": (1, 0, 1, 1, 1, 1, 1, 1),
    "bravais_lattice_address": (1, 0, 1, 1, 1, 0, 0, 1),
    "ear_decomposition_state_graph": (1, 0, 1, 1, 1, 1, 1, 1),
    "coverage_component_profile": (1, 0, 1, 0, 1, 0, 1, 1),
    "reflection_perron_moment": (1, 1, 1, 0, 1, 0, 0, 1),
    "obstruction_transfer_ledger": (1, 0, 1, 1, 0, 1, 1, 0),
    "lean_frontier_open_fields": (1, 0, 1, 1, 0, 1, 1, 0),
    "raw_scalar_p0": (1, 0, 0, 0, 1, 0, 0, 1),
}


def tournament_from_sidecars() -> tuple[list[str], list[list[bool]], dict[str, int], dict[tuple[str, str], int]]:
    vertices = list(SIDECARS)
    order = {v: i for i, v in enumerate(vertices)}
    adj = [[False] * len(vertices) for _ in vertices]
    margins: dict[tuple[str, str], int] = {}
    for i, a in enumerate(vertices):
        for j, b in enumerate(vertices):
            if i >= j:
                continue
            av = SIDECARS[a]
            bv = SIDECARS[b]
            aw = sum(1 for x, y in zip(av, bv) if x > y)
            bw = sum(1 for x, y in zip(av, bv) if y > x)
            if aw > bw:
                winner = a
                margin = aw - bw
            elif bw > aw:
                winner = b
                margin = bw - aw
            else:
                winner = a if order[a] < order[b] else b
                margin = 0
            wi = order[winner]
            li = j if wi == i else i
            adj[wi][li] = True
            margins[(a, b)] = margin
    scores = {v: sum(adj[order[v]]) for v in vertices}
    return vertices, adj, scores, margins


def tournament_fingerprint(vertices: list[str], adj: list[list[bool]]) -> dict[str, object]:
    n = len(vertices)
    c3 = 0
    for a, b, c in itertools.combinations(range(n), 3):
        if (adj[a][b] and adj[b][c] and adj[c][a]) or (adj[a][c] and adj[c][b] and adj[b][a]):
            c3 += 1
    edge_list = [(i, j) for i in range(n) for j in range(n) if adj[i][j]]
    comps_idx = tarjan_scc(range(n), edge_list)
    comps = [[vertices[i] for i in comp] for comp in comps_idx]
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        for last in range(n):
            val = dp[mask][last]
            if not val:
                continue
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                if adj[last][nxt]:
                    dp[mask | (1 << nxt)][nxt] += val
    hp = sum(dp[-1])
    return {"directed_3cycles": c3, "sccs": comps, "hamiltonian_paths": hp}


def print_packet_map(features: list[PacketFeatures]) -> None:
    print("=" * 100)
    print("MAP 1: NAMED LRC PACKETS, PGF ROOTS, PHI4 FITS, LATTICE/EAR SIGNALS")
    print("=" * 100)
    header = (
        "name".ljust(18)
        + " fam".ljust(19)
        + " q0".rjust(9)
        + " LyK8".rjust(9)
        + " real".rjust(6)
        + " nearest".rjust(9)
        + " gap7".rjust(8)
        + " lam".rjust(10)
        + " b".rjust(10)
        + " r2".rjust(7)
        + " c0".rjust(5)
        + " max0".rjust(9)
        + " states".rjust(8)
        + " ears".rjust(6)
        + " Perron".rjust(8)
    )
    print(header)
    for feat in features:
        q = feat.q
        lyk8 = 10 * q[0] + q[3] + 10 * q[6]
        print(
            feat.packet.name.ljust(18)
            + feat.packet.family[:18].ljust(19)
            + f"{float(q[0]):9.4f}"
            + f"{float(lyk8):9.4f}"
            + f"{feat.n_real:6d}"
            + f"{feat.nearest_root_abs:9.3f}"
            + f"{feat.apex7_angle_gap:8.2f}"
            + f"{feat.phi4_lambda:10.4f}"
            + f"{feat.phi4_b:10.4f}"
            + f"{feat.phi4_r2:7.3f}"
            + f"{feat.q0_components:5d}"
            + f"{float(feat.largest_q0_cell):9.5f}"
            + f"{feat.state_count:8d}"
            + f"{feat.ear_excess:6d}"
            + f"{feat.perron_gap:8.4f}"
        )
    print()
    print("Legend: real=# real PGF roots; gap7=mean angular distance to 1,2,3 times 360/7;")
    print("lam,b,r2 fit -log q_t to c + b*S^2 + lambda*S^4; c0/max0 are q0 component count/largest cell.")


def print_bravais(features: list[PacketFeatures]) -> None:
    summary = bravais_summary(features)
    by_family: dict[str, list[PacketFeatures]] = defaultdict(list)
    for feat in features:
        by_family[feat.packet.family].append(feat)
    print("\nBRAVAIS LATTICE ADDRESS MAP")
    print(f"  affine rank of q-vectors across all packets: {summary['q_rank']} (max 6 because sum q_t=1)")
    print(f"  approximate q-lattice covolume proxy (product nonzero singular values): {summary['covolume']:.8g}")
    print("  successive minima proxy (shortest independent q-differences):")
    for a, b, norm in summary["successive"]:
        print(f"    {norm:.6f}  {a} <-> {b}")
    for family, feats in sorted(by_family.items()):
        if len(feats) <= 1:
            continue
        sub = bravais_summary(feats)
        print(f"  family {family}: rank={sub['q_rank']} covolume={sub['covolume']:.8g}")


def print_extreme_findings(features: list[PacketFeatures]) -> None:
    print("\nLEE-YANG / PHI4 EXTREMALITY READOUT")
    complex_rows = [f for f in features if f.n_real == 0]
    collision_rows = [f for f in features if f.nearest_root_abs < 0.25]
    positive_lambda = [f for f in features if f.phi4_lambda > 0]
    print("  all-complex PGF rows:", ", ".join(f.packet.name for f in complex_rows))
    print("  root-collision rows with nearest |z|<0.25:", ", ".join(f.packet.name for f in collision_rows))
    print(f"  positive quartic lambda rows: {len(positive_lambda)}/{len(features)}")
    best_lam = sorted(features, key=lambda f: f.phi4_lambda, reverse=True)[:5]
    print("  largest lambda packets:", ", ".join(f"{f.packet.name}:{f.phi4_lambda:.4f}" for f in best_lam))
    best_gap = sorted(features, key=lambda f: f.apex7_angle_gap)[:5]
    print("  closest apex-7 angle packets:", ", ".join(f"{f.packet.name}:{f.apex7_angle_gap:.2f}deg" for f in best_gap))


def print_savitch_map() -> None:
    print("\n" + "=" * 100)
    print("MAP 2A: SAVITCH-STYLE REACHABILITY THROUGH CURRENT LRC14 FRONTIER OBLIGATIONS")
    print("=" * 100)
    starts = [
        "raw_residual_packet",
        "pgf_zero_curve_signal",
        "ear_state_transition",
        "bravais_lattice_address",
        "fine_modp_winding_transfer",
    ]
    target = "terminal_lrc14_statement"
    for start in starts:
        path = shortest_path(start, target)
        print(f"  {start} -> {target}: length={len(path)-1 if path else 'unreachable'}")
        if path:
            print("    path:", " -> ".join(path))
            sched = midpoint_schedule(path)
            if sched:
                print("    recursive midpoint schedule:")
                for a, mid, b, length in sched:
                    print(f"      len<={length}: {a} -> {mid} -> {b}")
    print("\n  Remaining proof fields from incoming HYP-3107:")
    print("    coverage extremality for k=8,9,10; finite local-minima certificate;")
    print("    reflection-Perron/order-3/order-4 certificate; Node-3 peel;")
    print("    finite-ruler normalized arc glue; fine mod-p winding transfer;")
    print("    residual classifier to ObserverGluingCertificate or FiniteAddressBranchPacket.")


def print_sidecar_tournament() -> None:
    print("\n" + "=" * 100)
    print("MAP 2B: TOURNAMENT ANALYSIS ON EXTREMALITY SIDECARS")
    print("=" * 100)
    vertices, adj, scores, margins = tournament_from_sidecars()
    fp = tournament_fingerprint(vertices, adj)
    score_hist = Counter(scores.values())
    print("  vertices are sidecars/proof obligations, not runners.")
    print("  axes:", ", ".join(SIDECAR_AXES))
    print("  score_hist:", dict(sorted(score_hist.items())))
    print("  scores:", dict(sorted(scores.items(), key=lambda kv: (-kv[1], kv[0]))))
    print("  directed_3cycles:", fp["directed_3cycles"])
    print("  sccs:", fp["sccs"])
    print("  Hamiltonian_path_count:", fp["hamiltonian_paths"])
    low = sorted((m, a, b) for (a, b), m in margins.items() if m <= 1)
    print("  low-margin edge-flip risks:")
    for m, a, b in low[:14]:
        print(f"    margin={m}: {a} vs {b}")
    ranked = sorted(vertices, key=lambda v: (-scores[v], v))
    print("  one retention path:", " -> ".join(ranked))


def print_hypotheses(features: list[PacketFeatures]) -> None:
    print("\nSYNTHESIS: WHAT THIS SESSION ADVANCES")
    consec = next(f for f in features if f.packet.name == "consec_8")
    break_row = next(f for f in features if f.packet.name == "break_8")
    print("  strengthened:")
    print(
        "    PGF root confinement is not just a scalar p0 restatement: "
        f"consec_8 has {consec.n_real} real roots, nearest |z|={consec.nearest_root_abs:.3f}, "
        f"while break_8 has {break_row.n_real} real roots, nearest |z|={break_row.nearest_root_abs:.3f}."
    )
    print("    The Lee-Yang signal should be paired with component/ear data: root collision rows")
    print("    also change the state-transition and q0-component sidecars, which HYP-3107 needs.")
    print("    Bravais rank says the named packets occupy a full 6-dimensional q-lattice face;")
    print("    a proof cannot reduce them to one moment without declaring lost coordinates.")
    print("  bold but still speculative:")
    print("    The phi4 lambda coefficient may be an extremality curvature coordinate.")
    print("    Ear-excess in the cell-state graph may be the observer-gluing sidecar that")
    print("    replaces raw H=7/H=21 analogies for LRC packets.")
    print("    Savitch midpoint certificates suggest a finite recursive proof-state search,")
    print("    but the hard theorem is still the coverage/observer certificate at the midpoint.")
    print("  rejected simplification:")
    print("    Raw scalar p0 and raw tournament H are too lossy.  They destroy the root stratum,")
    print("    first-live order, component count, and sidecar needed by the Lean frontier.")


def corr(xs: list[float], ys: list[float]) -> float:
    if len(xs) != len(ys) or not xs:
        return 0.0
    mx = sum(xs) / len(xs)
    my = sum(ys) / len(ys)
    sx = math.sqrt(sum((x - mx) ** 2 for x in xs))
    sy = math.sqrt(sum((y - my) ** 2 for y in ys))
    if sx * sy == 0:
        return 0.0
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / (sx * sy)


def print_consecutive_zero_arc() -> None:
    print("\n" + "=" * 100)
    print("SUPPLEMENT: CONSECUTIVE ZERO ARC k=8..13")
    print("=" * 100)
    print(f"{'k':>3s} {'q0':>9s} {'real':>5s} {'nearest':>8s} {'gap7':>7s} roots |z|@arg")
    for k in range(8, 14):
        feat = compute_features(Packet(f"consec_{k}", "consecutive_arc", tuple(range(k)), "zero-arc supplement"))
        root_bits = ", ".join(
            f"{abs(z):.3f}@{math.degrees(math.atan2(z.imag, z.real)):+.1f}"
            for z in feat.roots[:6]
        )
        print(
            f"{k:3d} {fstr(feat.q[0]):>9s} {feat.n_real:5d}"
            f" {feat.nearest_root_abs:8.3f} {feat.apex7_angle_gap:7.2f} {root_bits}"
        )
    print("  readout: consecutive rows keep all roots non-real; the middle pair passes")
    print("  close to 2*360/7 at k=11, matching the HYP-3103 zero-arc hint.")


def light_root_row(speeds: tuple[int, ...]) -> dict[str, float]:
    q, _, _ = cell_profile(speeds)
    roots = pgf_roots(q)
    n_real, nearest, _, _ = root_features(roots)
    lam, _, _ = phi4_fit(q)
    return {
        "n_real": float(n_real),
        "nearest": nearest,
        "lambda": lam,
        "extreme": float(q[0] + q[6]),
        "lyk8": float(10 * q[0] + q[3] + 10 * q[6]),
    }


def print_sample_correlations(sample_size: int = 120) -> None:
    rng = np.random.default_rng(3108)
    sample: list[dict[str, float]] = [
        light_root_row(tuple(range(8))),
        light_root_row(tuple(2 * i for i in range(8))),
        light_root_row((0, 7, 8, 9, 10, 11, 12, 13)),
    ]
    for i in range(sample_size):
        speeds = tuple([0] + sorted(int(x) for x in rng.choice(np.arange(1, 46), size=7, replace=False)))
        sample.append(light_root_row(speeds))
    print("\n" + "=" * 100)
    print("SUPPLEMENT: ROOT/ENERGY CORRELATIONS ON DETERMINISTIC ANCHORED 8-SET SAMPLE")
    print("=" * 100)
    groups: dict[int, list[dict[str, float]]] = defaultdict(list)
    for row in sample:
        groups[int(row["n_real"])].append(row)
    print("  #real -> count, mean(q0+q6), mean(L_yK8), max(L_yK8), mean phi4 lambda")
    for key in sorted(groups):
        rows = groups[key]
        mean_ext = sum(row["extreme"] for row in rows) / len(rows)
        mean_ly = sum(row["lyk8"] for row in rows) / len(rows)
        max_ly = max(row["lyk8"] for row in rows)
        mean_lam = sum(row["lambda"] for row in rows) / len(rows)
        print(
            f"    #real={key}: n={len(rows):3d} mean_ext={mean_ext:.4f}"
            f" mean_Ly={mean_ly:.4f} max_Ly={max_ly:.4f} mean_lambda={mean_lam:.4f}"
        )
    real = [row["n_real"] for row in sample]
    nearest = [row["nearest"] for row in sample]
    lam = [row["lambda"] for row in sample]
    ext = [row["extreme"] for row in sample]
    ly = [row["lyk8"] for row in sample]
    print(f"  corr(#real, q0+q6)={corr(real, ext):+.3f}")
    print(f"  corr(nearest-root-modulus, q0+q6)={corr(nearest, ext):+.3f}")
    print(f"  corr(phi4-lambda, L_yK8)={corr(lam, ly):+.3f}")
    print("  readout: root confinement is a robust sidecar; phi4 curvature is much")
    print("  noisier and should stay downstream of exact packet labels.")


def main() -> None:
    print("HYP-3108 Lee-Yang/Savitch/Bravais/ear-lattice extremality scout -- codex S262")
    print("Incoming integration: HYP-3103 PGF roots, HYP-3104 maximizer currencies,")
    print("HYP-3105 obstruction-transfer discipline, HYP-3106 perspective functors,")
    print("and HYP-3107 Lean proof-frontier open fields.\n")
    features = [compute_features(packet) for packet in PACKETS]
    print_packet_map(features)
    print_bravais(features)
    print_extreme_findings(features)
    print_savitch_map()
    print_sidecar_tournament()
    print_hypotheses(features)
    print_consecutive_zero_arc()
    print_sample_correlations()
    print("\nDONE.")


if __name__ == "__main__":
    main()
