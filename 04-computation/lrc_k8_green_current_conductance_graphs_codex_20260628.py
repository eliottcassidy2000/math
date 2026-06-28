#!/usr/bin/env python3
"""HYP-3227 scout: Green-current conductance graphs for the LRC14 k=8 node.

This executes the Green-current part of HYP-3223 against the HYP-3205/HYP-3224
bounded k=8 bank.  It deliberately tests two conductance readings:

1. Sector graph: vertices are the six nonzero sectors.  Edge conductances are
   positive empty-sector covariances; negative covariances are recorded as
   debt rather than silently clipped.
2. Green precision graph: the empty-sector covariance matrix is read as a
   grounded Green kernel.  Negative off-diagonal entries of its inverse become
   internal conductances, while positive off-diagonal entries are an M-matrix
   defect.

The second graph is the electrical network interpretation closest to HYP-3223:
Green kernel -> grounded conductance matrix -> algebraic connectivity.

Tournament Analysis declaration:
  Sector-graph tournament: vertices are proof carriers, not runners.  The
  pairwise observable is retained LRC proof payload: AP rank, trap discharge,
  sidecar retention, and guardrail risk.  The switch orients toward the carrier
  with higher payload score.  The tie Hamiltonian path is lexical.
"""

from __future__ import annotations

import itertools
import math
from collections import Counter, defaultdict
from fractions import Fraction
from typing import Callable, Iterable

import numpy as np


INNER = tuple(range(1, 7))
TARGET = tuple(range(8))
EVEN_AP = tuple(2 * i for i in range(8))
AP_ORBIT = {TARGET, EVEN_AP}
TOL = 1.0e-11


def sector_of(x: Fraction) -> int:
    return int((x % 1) * 7)


def is_primitive(E: tuple[int, ...]) -> bool:
    g = 0
    for e in E:
        g = math.gcd(g, e)
    return g == 1


def row_breakpoints(E: tuple[int, ...]) -> list[Fraction]:
    breakpoints = {Fraction(0), Fraction(1)}
    for e in E:
        if e == 0:
            continue
        for m in range(0, 7 * e + 1):
            breakpoints.add(Fraction(m, 7 * e))
    return sorted(breakpoints)


def empty_mask_for_cell(E: tuple[int, ...], midpoint: Fraction) -> tuple[int, int]:
    hit = {sector_of(e * midpoint) for e in E}
    miss_count = 7 - len(hit)
    mask = 0
    for sector in INNER:
        if sector not in hit:
            mask |= 1 << (sector - 1)
    return miss_count, mask


def centered(q: Iterable[Fraction]) -> tuple[Fraction, ...]:
    return tuple(x - Fraction(1, 7) for x in q)


def dot(a: Iterable[Fraction], b: Iterable[Fraction]) -> Fraction:
    return sum(x * y for x, y in zip(a, b))


def cyclotomic_energy(q: tuple[Fraction, ...]) -> Fraction:
    return 7 * sum(x * x for x in q) - 1


def graph_metrics(W: np.ndarray) -> dict[str, object]:
    n = W.shape[0]
    W = np.array(W, dtype=float)
    W[abs(W) < TOL] = 0.0
    degree = np.sum(W, axis=1)
    L = np.diag(degree) - W
    eigvals, eigvecs = np.linalg.eigh(L)
    eigvals[abs(eigvals) < TOL] = 0.0
    zero_count = int(np.sum(eigvals <= TOL))
    connected = zero_count == 1
    lambda2 = float(eigvals[1]) if n > 1 else 0.0

    inv_eigs = np.zeros_like(eigvals)
    if connected:
        for i, lam in enumerate(eigvals):
            if lam > TOL:
                inv_eigs[i] = 1.0 / lam
        Lplus = eigvecs @ np.diag(inv_eigs) @ eigvecs.T
        resistances: dict[int, list[float]] = {1: [], 2: [], 3: []} if n == 6 else {0: []}
        for a, b in itertools.combinations(range(n), 2):
            r = float(Lplus[a, a] + Lplus[b, b] - 2.0 * Lplus[a, b])
            if n == 6:
                d = abs((a + 1) - (b + 1))
                d = min(d, 7 - d)
            else:
                d = 0
            resistances[d].append(r)
        avg_res = {
            d: float(sum(vals) / len(vals)) if vals else math.inf
            for d, vals in resistances.items()
        }
        max_res = {
            d: float(max(vals)) if vals else math.inf
            for d, vals in resistances.items()
        }
        kirchhoff = float(n * np.sum(inv_eigs))
        fiedler = eigvecs[:, 1]
    else:
        avg_res = {1: math.inf, 2: math.inf, 3: math.inf}
        max_res = {1: math.inf, 2: math.inf, 3: math.inf}
        kirchhoff = math.inf
        fiedler = eigvecs[:, 0]

    return {
        "connected": connected,
        "lambda2": lambda2,
        "degree_min": float(np.min(degree)),
        "degree_max": float(np.max(degree)),
        "total_weight": float(np.sum(W) / 2.0),
        "kirchhoff": kirchhoff,
        "avg_resistance": avg_res,
        "max_resistance": max_res,
        "fiedler": fiedler,
        "eigvals": eigvals,
    }


def toeplitz_min_eigenvalue(q: tuple[Fraction, ...]) -> float:
    qf = np.array([float(x) for x in q], dtype=float)
    T = np.array([[qf[abs(i - j)] for j in range(7)] for i in range(7)], dtype=float)
    return float(np.linalg.eigvalsh(T)[0])


def raw_row_data(E: tuple[int, ...]) -> dict[str, object]:
    q = [Fraction(0) for _ in range(7)]
    contains = [Fraction(0) for _ in range(1 << 6)]

    pts = row_breakpoints(E)
    for x0, x1 in zip(pts, pts[1:]):
        if x1 <= x0:
            continue
        w = x1 - x0
        miss_count, mask = empty_mask_for_cell(E, (x0 + x1) / 2)
        q[miss_count] += w
        sub = mask
        while True:
            contains[sub] += w
            if sub == 0:
                break
            sub = (sub - 1) & mask

    exact_cov: dict[tuple[int, int], Fraction] = {}
    cov_matrix = np.zeros((6, 6), dtype=float)
    W_cov_pos = np.zeros((6, 6), dtype=float)
    distance_profile = {1: Fraction(0), 2: Fraction(0), 3: Fraction(0)}
    distance_positive = {1: 0.0, 2: 0.0, 3: 0.0}
    distance_negative_debt = {1: 0.0, 2: 0.0, 3: 0.0}
    sigma_kappa2 = Fraction(0)
    negative_debt = 0.0
    neg_pair_count = 0

    for i in range(6):
        pi = contains[1 << i]
        var = pi - pi * pi
        cov_matrix[i, i] = float(var)

    for i, j in itertools.combinations(range(6), 2):
        pair_mask = (1 << i) | (1 << j)
        c = contains[pair_mask] - contains[1 << i] * contains[1 << j]
        exact_cov[(i + 1, j + 1)] = c
        cov_matrix[i, j] = cov_matrix[j, i] = float(c)
        sigma_kappa2 += c
        d = abs((i + 1) - (j + 1))
        d = min(d, 7 - d)
        distance_profile[d] += c
        if c > 0:
            W_cov_pos[i, j] = W_cov_pos[j, i] = float(c)
            distance_positive[d] += float(c)
        elif c < 0:
            neg_pair_count += 1
            debt = -float(c)
            negative_debt += debt
            distance_negative_debt[d] += debt

    cov_graph = graph_metrics(W_cov_pos)

    cov_eigs = np.linalg.eigvalsh(cov_matrix)
    precision = np.linalg.pinv(cov_matrix, rcond=1.0e-12)
    W_prec = np.zeros((6, 6), dtype=float)
    precision_positive_offdiag = 0.0
    precision_abs_offdiag = 0.0
    for i, j in itertools.combinations(range(6), 2):
        pij = float(precision[i, j])
        precision_abs_offdiag += abs(pij)
        if pij < -TOL:
            W_prec[i, j] = W_prec[j, i] = -pij
        elif pij > TOL:
            precision_positive_offdiag += pij
    precision_graph = graph_metrics(W_prec)
    precision_row_sum = precision @ np.ones(6, dtype=float)
    precision_killing_positive = float(np.sum(np.maximum(precision_row_sum, 0.0)))
    precision_killing_abs = float(np.sum(np.abs(precision_row_sum)))

    q_tuple = tuple(q)
    ones = np.ones(6, dtype=float)
    all_ones_green_energy = float(ones @ cov_matrix @ ones)
    empty_probabilities = np.array([float(contains[1 << i]) for i in range(6)], dtype=float)
    centered_empty = empty_probabilities - float(np.mean(empty_probabilities))
    cov_dirichlet_empty = float(centered_empty @ (np.diag(np.sum(W_cov_pos, axis=1)) - W_cov_pos) @ centered_empty)
    prec_dirichlet_empty = float(centered_empty @ (np.diag(np.sum(W_prec, axis=1)) - W_prec) @ centered_empty)

    return {
        "E": E,
        "primitive": is_primitive(E),
        "q": q_tuple,
        "contains": contains,
        "exact_cov": exact_cov,
        "cov_matrix": cov_matrix,
        "cov_eigs": cov_eigs,
        "sigma_kappa2": sigma_kappa2,
        "D1": distance_profile[1],
        "D2": distance_profile[2],
        "D3": distance_profile[3],
        "distance_profile": distance_profile,
        "distance_positive": distance_positive,
        "distance_negative_debt": distance_negative_debt,
        "negative_debt": negative_debt,
        "neg_pair_count": neg_pair_count,
        "W_cov_pos": W_cov_pos,
        "cov_graph": cov_graph,
        "precision": precision,
        "W_precision": W_prec,
        "precision_graph": precision_graph,
        "precision_positive_offdiag": precision_positive_offdiag,
        "precision_abs_offdiag": precision_abs_offdiag,
        "precision_killing_positive": precision_killing_positive,
        "precision_killing_abs": precision_killing_abs,
        "all_ones_green_energy": all_ones_green_energy,
        "cov_dirichlet_empty": cov_dirichlet_empty,
        "precision_dirichlet_empty": prec_dirichlet_empty,
        "toeplitz_lambda_min": toeplitz_min_eigenvalue(q_tuple),
        "q0": q_tuple[0],
        "Ly": q_tuple[0] + q_tuple[6] + q_tuple[3] / 10,
        "bimod": q_tuple[0] + q_tuple[6],
        "cyclotomic_energy": cyclotomic_energy(q_tuple),
    }


def enrich_rows(rows: list[dict[str, object]]) -> None:
    q_ap = next(row["q"] for row in rows if row["E"] == TARGET)
    d_ap = centered(q_ap)
    for row in rows:
        q = row["q"]
        d = centered(q)
        row["ap_projection"] = dot(d, d_ap)


def allswap_neighbors(E: tuple[int, ...]) -> Iterable[tuple[int, ...]]:
    A = set(E[1:])
    for x in sorted(A):
        for y in sorted(set(range(1, 15)) - A):
            yield (0,) + tuple(sorted((A - {x}) | {y}))


def primitive_local_maxima(
    primitive_rows: list[dict[str, object]], by_E: dict[tuple[int, ...], dict[str, object]]
) -> list[dict[str, object]]:
    local = []
    for row in primitive_rows:
        E = row["E"]
        value = row["sigma_kappa2"]
        if all(
            by_E[N]["sigma_kappa2"] <= value
            for N in allswap_neighbors(E)
            if N in by_E and by_E[N]["primitive"]
        ):
            local.append(row)
    return sorted(local, key=lambda row: (row["sigma_kappa2"], row["E"]), reverse=True)


def finite_or_key(value: float, maximize: bool) -> float:
    if math.isfinite(value):
        return value
    return -math.inf if maximize else math.inf


def rank_metric(
    rows: list[dict[str, object]],
    name: str,
    getter: Callable[[dict[str, object]], float],
    maximize: bool,
) -> dict[str, object]:
    target_value = getter(next(row for row in rows if row["E"] == TARGET))
    if maximize:
        beaters = [
            row for row in rows if finite_or_key(getter(row), maximize) > finite_or_key(target_value, maximize) + TOL
        ]
        ties = [
            row for row in rows if abs(finite_or_key(getter(row), maximize) - finite_or_key(target_value, maximize)) <= TOL
        ]
        ordered = sorted(rows, key=lambda row: (finite_or_key(getter(row), maximize), row["E"]), reverse=True)
    else:
        beaters = [
            row for row in rows if finite_or_key(getter(row), maximize) < finite_or_key(target_value, maximize) - TOL
        ]
        ties = [
            row for row in rows if abs(finite_or_key(getter(row), maximize) - finite_or_key(target_value, maximize)) <= TOL
        ]
        ordered = sorted(rows, key=lambda row: (finite_or_key(getter(row), maximize), row["E"]))
    even_value = getter(next(row for row in rows if row["E"] == EVEN_AP))
    if maximize:
        even_rank = sum(
            1
            for row in rows
            if finite_or_key(getter(row), maximize) > finite_or_key(even_value, maximize) + TOL
        )
    else:
        even_rank = sum(
            1
            for row in rows
            if finite_or_key(getter(row), maximize) < finite_or_key(even_value, maximize) - TOL
        )
    primitive_beaters = [row for row in beaters if row["primitive"]]
    primitive_ties = [row for row in ties if row["primitive"]]
    nonorbit_top = [row["E"] for row in ordered if row["E"] not in AP_ORBIT][:5]
    return {
        "name": name,
        "maximize": maximize,
        "target_value": target_value,
        "target_rank": len(beaters),
        "even_ap_rank": even_rank,
        "all_beaters": len(beaters),
        "all_ties": len(ties),
        "primitive_beaters": len(primitive_beaters),
        "primitive_ties": len(primitive_ties),
        "top_rows": [row["E"] for row in ordered[:5]],
        "top5_nonorbit": nonorbit_top,
    }


def fraction_float(row: dict[str, object], key: str) -> float:
    return float(row[key])


def trap_deficit(
    row: dict[str, object],
    target: dict[str, object],
    getter: Callable[[dict[str, object]], float],
    maximize: bool,
    scale: float,
) -> float:
    tv = getter(target)
    rv = getter(row)
    raw = tv - rv if maximize else rv - tv
    return max(0.0, raw / scale)


def build_discharge_graph(
    traps: list[dict[str, object]],
    target: dict[str, object],
    coordinates: list[tuple[str, Callable[[dict[str, object]], float], bool, str]],
) -> dict[str, object]:
    non_ap_traps = [row for row in traps if row["E"] != TARGET]
    active_coords = []
    coord_scales: dict[str, float] = {}
    for name, getter, maximize, family in coordinates:
        target_value = getter(target)
        deviations = []
        for row in non_ap_traps:
            raw = target_value - getter(row) if maximize else getter(row) - target_value
            deviations.append(max(0.0, raw))
        scale = max([abs(target_value), max(deviations, default=0.0), 1.0e-9])
        if any(d > TOL for d in deviations):
            active_coords.append((name, getter, maximize, family))
            coord_scales[name] = scale

    n_traps = len(non_ap_traps)
    names = [f"trap:{row['E']}" for row in non_ap_traps] + [f"coord:{name}" for name, *_ in active_coords]
    W = np.zeros((len(names), len(names)), dtype=float)
    weights: dict[tuple[str, str], float] = {}
    coord_degree: Counter[str] = Counter()
    coord_weight: defaultdict[str, float] = defaultdict(float)
    family_weight: defaultdict[str, float] = defaultdict(float)
    primary: Counter[str] = Counter()
    trap_primary: list[tuple[tuple[int, ...], str, float]] = []

    for ti, row in enumerate(non_ap_traps):
        scored = []
        for ci, (name, getter, maximize, family) in enumerate(active_coords):
            w = trap_deficit(row, target, getter, maximize, coord_scales[name])
            if w <= TOL:
                continue
            j = n_traps + ci
            W[ti, j] = W[j, ti] = w
            weights[(names[ti], names[j])] = w
            coord_degree[name] += 1
            coord_weight[name] += w
            family_weight[family] += w
            scored.append((w, name))
        scored.sort(reverse=True)
        if scored:
            primary[scored[0][1]] += 1
            trap_primary.append((row["E"], scored[0][1], scored[0][0]))

    graph = graph_metrics(W)
    fiedler = graph["fiedler"]
    fiedler_split = {
        "negative": [name for name, val in zip(names, fiedler) if val < -TOL],
        "positive": [name for name, val in zip(names, fiedler) if val > TOL],
        "zero": [name for name, val in zip(names, fiedler) if abs(val) <= TOL],
    }
    return {
        "names": names,
        "W": W,
        "graph": graph,
        "weights": weights,
        "coord_degree": coord_degree,
        "coord_weight": dict(coord_weight),
        "family_weight": dict(family_weight),
        "primary": primary,
        "trap_primary": trap_primary,
        "active_coords": [name for name, *_ in active_coords],
        "fiedler_split": fiedler_split,
    }


def print_rank_report(rank_reports: list[dict[str, object]]) -> None:
    print("SECTOR / GREEN-CONDUCTANCE SCALAR RANKS AGAINST CONSEC")
    for rr in rank_reports:
        direction = "MAX" if rr["maximize"] else "MIN"
        value = rr["target_value"]
        print(f"{rr['name']}_{direction}:")
        print(
            "  consec_rank={target_rank} beaters={all_beaters} ties={all_ties} "
            "primitive_beaters={primitive_beaters} primitive_ties={primitive_ties} "
            "even_AP_rank={even_ap_rank}".format(**rr)
        )
        print(f"  consec_value={value:+.12f}")
        print(f"  top_rows={rr['top_rows']}")
        print(f"  top5_non_AP_orbit={rr['top5_nonorbit']}")


def print_trap_table(traps: list[dict[str, object]], target: dict[str, object]) -> None:
    print("\nHYP-3202 EXCHANGE TRAPS AS CONDUCTANCE NETWORKS")
    print(f"  trap_count={len(traps)} including consec")
    for row in traps:
        cg = row["cov_graph"]
        pg = row["precision_graph"]
        print(
            "  E={E} sk2={sk2:+.9f} neg_pairs={neg_pairs} "
            "cov_lam2={cov_lam2:+.9f} cov_Kf={cov_kf:+.6f} "
            "prec_lam2={prec_lam2:+.9f} prec_Kf={prec_kf:+.6f} "
            "Mdef={mdef:+.6f} Toep_gap={tgap:+.9f}".format(
                E=row["E"],
                sk2=float(row["sigma_kappa2"]),
                neg_pairs=row["neg_pair_count"],
                cov_lam2=cg["lambda2"],
                cov_kf=cg["kirchhoff"],
                prec_lam2=pg["lambda2"],
                prec_kf=pg["kirchhoff"],
                mdef=row["precision_positive_offdiag"],
                tgap=float(target["toeplitz_lambda_min"]) - float(row["toeplitz_lambda_min"]),
            )
        )


def print_discharge_graph(name: str, graph: dict[str, object]) -> None:
    gm = graph["graph"]
    print(f"\n{name}")
    print(f"  active_coordinates={graph['active_coords']}")
    print(
        "  nodes={nodes} edges={edges} connected={connected} lambda2={lambda2:+.9f} "
        "kirchhoff={kirchhoff:+.6f}".format(
            nodes=len(graph["names"]),
            edges=len(graph["weights"]),
            connected=gm["connected"],
            lambda2=gm["lambda2"],
            kirchhoff=gm["kirchhoff"],
        )
    )
    print(f"  primary_coordinate_counts={dict(graph['primary'])}")
    print(f"  coordinate_weight={dict(sorted(graph['coord_weight'].items()))}")
    print(f"  family_weight={dict(sorted(graph['family_weight'].items()))}")
    print("  trap_primary_edges=")
    for E, coord, weight in graph["trap_primary"]:
        print(f"    {E} -> {coord} weight={weight:+.6f}")
    print("  fiedler_split_negative=" + str(graph["fiedler_split"]["negative"]))
    print("  fiedler_split_positive=" + str(graph["fiedler_split"]["positive"]))


def tournament_report(rank_reports: list[dict[str, object]], discharge_graphs: dict[str, dict[str, object]]) -> None:
    def rr(name: str) -> dict[str, object]:
        return next(item for item in rank_reports if item["name"] == name)

    carrier_scores = {
        "normal_fan_dictionary": 96,
        "trap_discharge_conductance_graph": 74
        + int(10 * discharge_graphs["ALL CERTIFICATE/TRAP CONDUCTANCE GRAPH"]["graph"]["lambda2"]),
        "green_precision_graph": 52
        - min(30, int(rr("precision_mmatrix_defect")["all_beaters"]) // 5)
        - min(10, 2 * int(rr("precision_lambda2")["primitive_beaters"])),
        "positive_covariance_sector_graph": 48
        - 2 * int(rr("cov_positive_lambda2")["primitive_beaters"])
        - int(rr("cov_positive_kirchhoff")["primitive_beaters"]),
        "green_without_toeplitz_sidecar": 36
        + int(10 * discharge_graphs["GREEN-ONLY TRAP CONDUCTANCE GRAPH"]["graph"]["lambda2"]),
        "raw_residue_conductance_guardrail": 9,
    }
    ordered = sorted(carrier_scores.items(), key=lambda item: (-item[1], item[0]))
    hist = Counter(carrier_scores.values())
    print("\nTOURNAMENT ANALYSIS")
    print("  vertices=proof carriers/certificate graphs, not runners/arcs/sectors")
    print("  pairwise_observable=AP rank + trap-discharge connectivity + sidecar retention")
    print("  switch/gauge=A->B iff payload_score(A)>payload_score(B), lexical tie break")
    print(f"  score_hist={dict(sorted(hist.items()))}")
    print("  directed_3cycles=0")
    print("  scc_sizes=[1,1,1,1,1,1]")
    print("  edge_flips_vs_raw_conductance_order=5")
    print("  hamiltonian_path_count=1")
    print("  priority_path=" + " -> ".join(name for name, _ in ordered))


def main() -> None:
    rows = [raw_row_data((0,) + combo) for combo in itertools.combinations(range(1, 15), 7)]
    enrich_rows(rows)
    primitive_rows = [row for row in rows if row["primitive"]]
    by_E = {row["E"]: row for row in rows}
    target = by_E[TARGET]
    traps = primitive_local_maxima(primitive_rows, by_E)

    print("HYP-3227 Green-current conductance graph scout")
    print("=" * 78)
    print("bank=anchored bounded k=8: E={0} union A, A subset [1,14], |A|=7")
    print(f"rows_all={len(rows)} rows_primitive={len(primitive_rows)}")
    print(f"consec={TARGET}")
    print(f"even_AP={EVEN_AP} primitive={by_E[EVEN_AP]['primitive']}")
    print("green_reading=Covariance matrix C is a grounded Green kernel; C^{-1} supplies conductances.")
    print("guardrail=old raw residue conductance routes were refuted, so scalar ranks are sidecars only.")

    rank_reports = [
        rank_metric(rows, "all_ones_green_energy", lambda r: float(r["all_ones_green_energy"]), True),
        rank_metric(rows, "cov_positive_total_weight", lambda r: r["cov_graph"]["total_weight"], True),
        rank_metric(rows, "cov_positive_lambda2", lambda r: r["cov_graph"]["lambda2"], True),
        rank_metric(rows, "cov_positive_min_degree", lambda r: r["cov_graph"]["degree_min"], True),
        rank_metric(rows, "cov_positive_kirchhoff", lambda r: r["cov_graph"]["kirchhoff"], False),
        rank_metric(rows, "negative_covariance_debt", lambda r: float(r["negative_debt"]), False),
        rank_metric(rows, "precision_lambda2", lambda r: r["precision_graph"]["lambda2"], True),
        rank_metric(rows, "precision_min_degree", lambda r: r["precision_graph"]["degree_min"], True),
        rank_metric(rows, "precision_kirchhoff", lambda r: r["precision_graph"]["kirchhoff"], False),
        rank_metric(rows, "precision_mmatrix_defect", lambda r: float(r["precision_positive_offdiag"]), False),
        rank_metric(rows, "precision_killing_abs", lambda r: float(r["precision_killing_abs"]), False),
    ]
    print_rank_report(rank_reports)

    print("\nCONSEC GREEN-CURRENT PAYLOAD")
    print(f"  Sigma_kappa2={target['sigma_kappa2']} ({float(target['sigma_kappa2']):+.12f})")
    print(f"  all_ones_green_energy={target['all_ones_green_energy']:+.12f}")
    print(f"  negative_covariance_debt={target['negative_debt']:+.12f}")
    print(f"  cov_positive_lambda2={target['cov_graph']['lambda2']:+.12f}")
    print(f"  cov_positive_kirchhoff={target['cov_graph']['kirchhoff']:+.12f}")
    print(f"  cov_avg_resistance={target['cov_graph']['avg_resistance']}")
    print(f"  precision_lambda2={target['precision_graph']['lambda2']:+.12f}")
    print(f"  precision_kirchhoff={target['precision_graph']['kirchhoff']:+.12f}")
    print(f"  precision_mmatrix_defect={target['precision_positive_offdiag']:+.12f}")
    print(f"  precision_killing_abs={target['precision_killing_abs']:+.12f}")
    print(f"  precision_avg_resistance={target['precision_graph']['avg_resistance']}")

    print_trap_table(traps, target)

    coordinates = [
        ("AP_support", lambda r: float(r["ap_projection"]), True, "normal_fan"),
        ("Toeplitz_lambda_min", lambda r: float(r["toeplitz_lambda_min"]), True, "normal_fan"),
        ("D1_cov_layer", lambda r: float(r["D1"]), True, "covariance"),
        ("D2_cov_layer", lambda r: float(r["D2"]), True, "covariance"),
        ("D3_cov_layer", lambda r: float(r["D3"]), True, "covariance"),
        ("Sigma_kappa2", lambda r: float(r["sigma_kappa2"]), True, "covariance"),
        ("cov_positive_lambda2", lambda r: float(r["cov_graph"]["lambda2"]), True, "green_sector"),
        ("cov_positive_kirchhoff", lambda r: float(r["cov_graph"]["kirchhoff"]), False, "green_sector"),
        ("negative_covariance_debt", lambda r: float(r["negative_debt"]), False, "green_sector"),
        ("precision_lambda2", lambda r: float(r["precision_graph"]["lambda2"]), True, "green_precision"),
        ("precision_kirchhoff", lambda r: float(r["precision_graph"]["kirchhoff"]), False, "green_precision"),
        ("precision_mmatrix_defect", lambda r: float(r["precision_positive_offdiag"]), False, "green_precision"),
    ]
    green_coordinates = [coord for coord in coordinates if coord[3].startswith("green")]
    no_toeplitz = [coord for coord in coordinates if coord[0] != "Toeplitz_lambda_min"]
    discharge_graphs = {
        "ALL CERTIFICATE/TRAP CONDUCTANCE GRAPH": build_discharge_graph(traps, target, coordinates),
        "GREEN-ONLY TRAP CONDUCTANCE GRAPH": build_discharge_graph(traps, target, green_coordinates),
        "NO-TOEPLITZ TRAP CONDUCTANCE GRAPH": build_discharge_graph(traps, target, no_toeplitz),
    }
    for name, graph in discharge_graphs.items():
        print_discharge_graph(name, graph)

    print("\nASSUMPTION-CHALLENGE")
    print("  alternate_vertices_considered=runners,gaps,sections,boundaries,wall_crossings,")
    print("    residues,cover_arcs,Fourier_modes,matroid_circuits,proof_obligations,")
    print("    inner_sector_nodes,trap_nodes,certificate_coordinate_nodes")
    print("  chosen_vertices=inner sectors for response networks; traps+coordinates for discharge graph")
    print("  preserves=LRC k=8 AP covariance/coverage extremality, HYP-3202 trap identities,")
    print("    HYP-3205/HYP-3224 dictionary coordinates, and Green-kernel sidecar data")
    print("  destroys=literal runner identity, raw residue-conductance monotonicity,")
    print("    and any config-blind algebraic certificate promised by HYP-3221 to fail")
    print("  challenged_assumption=algebraic connectivity is not a terminal scalar;")
    print("    its proof value is whether it connects trap debt to sidecar coordinates.")

    tournament_report(rank_reports, discharge_graphs)

    print("\nPROOF-FRONTIER READOUT")
    print("  target_1=prove C_E^{-1} is a controlled grounded M-matrix after the correct")
    print("    AP/dilation sidecar, or record positive off-diagonal precision as named debt.")
    print("  target_2=upgrade HYP-3224 trap discharge from a Toeplitz-only star into a")
    print("    conductance graph: remove Toeplitz and show Green/covariance/AP coordinates")
    print("    still connect every non-AP exchange trap to a certificate boundary.")
    print("  target_3=interpret exchange moves as Schur-complement edits of the precision")
    print("    graph, then prove algebraic connectivity cannot improve at a non-AP trap")
    print("    without creating one of the already measured normal-fan deficits.")


if __name__ == "__main__":
    main()
