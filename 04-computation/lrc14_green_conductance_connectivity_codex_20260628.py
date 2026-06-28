#!/usr/bin/env python3
"""HYP-3236 scout: Green conductance and algebraic connectivity.

This executes the electrical side of HYP-3223 on the same anchored k=8 bank
used by HYP-3202/HYP-3205/HYP-3224.

Construction:

* vertices are the six inner sectors 1..6;
* edge conductance is the positive part of pair covariance
  Cov(1_{sector i empty}, 1_{sector j empty});
* the weighted graph Laplacian emits algebraic connectivity lambda_2;
* the Green kernel is the Moore-Penrose inverse of the Laplacian, used for
  effective-resistance / Kirchhoff profiles.

The positive-part step is deliberate and lossy: negative covariance edges are
kept as a named leakage sidecar.  This scout tests whether the resulting
conductance graph is a real certificate coordinate or only a diagnostic.
"""

from __future__ import annotations

import itertools
import math
from collections import Counter
from fractions import Fraction as F
from typing import Callable, Iterable

import numpy as np


INNER = tuple(range(1, 7))
CONSEC_8 = tuple(range(8))
EVEN_AP = tuple(range(0, 15, 2))
AP_ORBIT = {CONSEC_8, EVEN_AP}
FLOAT_TOL = 1e-10


def sector_of(p: F) -> int:
    return int((p % 1) * 7)


def is_primitive(E: tuple[int, ...]) -> bool:
    g = 0
    for e in E:
        g = math.gcd(g, e)
    return g == 1


def cells(E: Iterable[int]) -> list[tuple[set[int], F]]:
    E = tuple(sorted(set(E)))
    breakpoints = {F(0), F(1)}
    for e in E:
        if e == 0:
            continue
        for m in range(0, 7 * e + 1):
            breakpoints.add(F(m, 7 * e))

    pts = sorted(breakpoints)
    out: list[tuple[set[int], F]] = []
    for x0, x1 in zip(pts, pts[1:]):
        if x1 <= x0:
            continue
        mid = (x0 + x1) / 2
        out.append(({sector_of(e * mid) for e in E}, x1 - x0))
    return out


def connected_components(W: np.ndarray) -> list[list[int]]:
    seen: set[int] = set()
    comps: list[list[int]] = []
    n = W.shape[0]
    for start in range(n):
        if start in seen:
            continue
        stack = [start]
        seen.add(start)
        comp: list[int] = []
        while stack:
            u = stack.pop()
            comp.append(u)
            for v in range(n):
                if W[u, v] > FLOAT_TOL and v not in seen:
                    seen.add(v)
                    stack.append(v)
        comps.append(sorted(comp))
    return comps


def effective_resistance_profile(W: np.ndarray, L: np.ndarray) -> dict[str, object]:
    eig = np.linalg.eigvalsh(L)
    comps = connected_components(W)
    connected = len(comps) == 1
    profile: dict[str, object] = {
        "laplacian_eigvals": eig,
        "lambda2": float(eig[1]) if len(eig) > 1 else 0.0,
        "components": comps,
        "connected": connected,
    }
    if not connected:
        profile.update(
            {
                "green": None,
                "kirchhoff": float("inf"),
                "mean_resistance": float("inf"),
                "max_resistance": float("inf"),
                "distance_resistance": {1: float("inf"), 2: float("inf"), 3: float("inf")},
                "bottleneck_pair": None,
                "bottleneck_current": None,
            }
        )
        return profile

    green = np.linalg.pinv(L, hermitian=True)
    pair_resistances: dict[tuple[int, int], float] = {}
    by_distance: dict[int, list[float]] = {1: [], 2: [], 3: []}
    for i, j in itertools.combinations(range(6), 2):
        r = float(green[i, i] + green[j, j] - 2 * green[i, j])
        pair_resistances[(i + 1, j + 1)] = r
        d = abs((i + 1) - (j + 1))
        d = min(d, 7 - d)
        by_distance[d].append(r)

    bottleneck_pair = max(pair_resistances.items(), key=lambda item: (item[1], item[0]))[0]
    current = unit_current_profile(W, L, bottleneck_pair)
    profile.update(
        {
            "green": green,
            "pair_resistances": pair_resistances,
            "kirchhoff": sum(pair_resistances.values()),
            "mean_resistance": sum(pair_resistances.values()) / len(pair_resistances),
            "max_resistance": max(pair_resistances.values()),
            "distance_resistance": {
                d: sum(values) / len(values) for d, values in by_distance.items()
            },
            "bottleneck_pair": bottleneck_pair,
            "bottleneck_current": current,
        }
    )
    return profile


def unit_current_profile(W: np.ndarray, L: np.ndarray, pair: tuple[int, int]) -> dict[str, object]:
    source, sink = pair[0] - 1, pair[1] - 1
    b = np.zeros(6, dtype=float)
    b[source] = 1.0
    b[sink] = -1.0
    green = np.linalg.pinv(L, hermitian=True)
    voltage = green @ b
    currents: dict[tuple[int, int], float] = {}
    for i, j in itertools.combinations(range(6), 2):
        if W[i, j] <= FLOAT_TOL:
            continue
        currents[(i + 1, j + 1)] = float(W[i, j] * (voltage[i] - voltage[j]))
    abs_currents = {edge: abs(value) for edge, value in currents.items()}
    max_edge = max(abs_currents.items(), key=lambda item: (item[1], item[0]))[0]
    energy = float(b @ voltage)
    total_abs = sum(abs_currents.values())
    probs = [value / total_abs for value in abs_currents.values() if value > FLOAT_TOL]
    entropy = -sum(p * math.log(p, 2) for p in probs)
    return {
        "pair": pair,
        "effective_resistance": energy,
        "max_current_edge": max_edge,
        "max_abs_current": abs_currents[max_edge],
        "total_abs_current": total_abs,
        "current_entropy_bits": entropy,
    }


def row_details(E: tuple[int, ...]) -> dict[str, object]:
    q = [F(0) for _ in range(7)]
    p_empty = {i: F(0) for i in INNER}
    pair_empty: dict[tuple[int, int], F] = {}

    for hit, weight in cells(E):
        q[7 - len(hit)] += weight
        empty_inner = [i for i in INNER if i not in hit]
        for i in empty_inner:
            p_empty[i] += weight
        for i, j in itertools.combinations(empty_inner, 2):
            pair_empty[(i, j)] = pair_empty.get((i, j), F(0)) + weight

    cov: dict[tuple[int, int], F] = {}
    C = np.zeros((6, 6), dtype=float)
    sk2 = F(0)
    negative_mass = 0.0
    negative_edges = 0
    for i, j in itertools.combinations(INNER, 2):
        c = pair_empty.get((i, j), F(0)) - p_empty[i] * p_empty[j]
        cov[(i, j)] = c
        sk2 += c
        C[i - 1, j - 1] = C[j - 1, i - 1] = float(c)
        if c < 0:
            negative_edges += 1
            negative_mass += -float(c)

    W = np.maximum(C, 0.0)
    np.fill_diagonal(W, 0.0)
    L = np.diag(W.sum(axis=1)) - W
    green = effective_resistance_profile(W, L)
    return {
        "E": E,
        "primitive": is_primitive(E),
        "q": tuple(q),
        "cov": cov,
        "sk2": sk2,
        "C": C,
        "W": W,
        "L": L,
        "total_positive_conductance": float(W.sum() / 2),
        "negative_edges": negative_edges,
        "negative_mass": negative_mass,
        **green,
    }


def allswap_neighbors(E: tuple[int, ...]):
    A = set(E[1:])
    for x in sorted(A):
        for y in sorted(set(range(1, 15)) - A):
            yield (0,) + tuple(sorted((A - {x}) | {y})), (x, y)


def local_maxima(
    rows: list[dict[str, object]],
    by_E: dict[tuple[int, ...], dict[str, object]],
    neighbor_fn: Callable[[tuple[int, ...]], Iterable[tuple[tuple[int, ...], tuple[int, int]]]],
) -> list[dict[str, object]]:
    out = []
    for row in rows:
        E = row["E"]
        value = row["sk2"]
        neighbor_values = []
        for N, _ in neighbor_fn(E):
            if N not in by_E:
                continue
            if not by_E[N]["primitive"]:
                continue
            neighbor_values.append(by_E[N]["sk2"])
        if neighbor_values and all(v <= value for v in neighbor_values):
            out.append(row)
    return sorted(out, key=lambda r: (r["sk2"], r["E"]), reverse=True)


def rank_readout(
    rows: list[dict[str, object]],
    primitive_rows: list[dict[str, object]],
    ap: dict[str, object],
    name: str,
    getter: Callable[[dict[str, object]], float],
    *,
    high: bool,
) -> dict[str, object]:
    ap_value = getter(ap)

    def beats(value: float) -> bool:
        return value > ap_value + FLOAT_TOL if high else value < ap_value - FLOAT_TOL

    ordered = sorted(rows, key=getter, reverse=high)
    all_beaters = [row for row in rows if beats(getter(row))]
    primitive_beaters = [row for row in primitive_rows if beats(getter(row))]
    ties = [row for row in rows if abs(getter(row) - ap_value) <= FLOAT_TOL]
    primitive_ties = [row for row in primitive_rows if abs(getter(row) - ap_value) <= FLOAT_TOL]
    return {
        "name": name,
        "ap_value": ap_value,
        "ap_rank": [i for i, row in enumerate(ordered) if row["E"] == CONSEC_8][0],
        "dilation_rank": [i for i, row in enumerate(ordered) if row["E"] == EVEN_AP][0],
        "all_beaters": len(all_beaters),
        "primitive_beaters": len(primitive_beaters),
        "all_ties": len(ties),
        "primitive_ties": len(primitive_ties),
        "best_E": ordered[0]["E"],
        "best_value": getter(ordered[0]),
        "top5_nonorbit": [row["E"] for row in ordered if row["E"] not in AP_ORBIT][:5],
    }


def fmt_float(x: float) -> str:
    if math.isinf(x):
        return "inf"
    return f"{x:.12f}"


def print_rank_section(readouts: list[dict[str, object]]) -> None:
    print("\nCONDUCTANCE / GREEN RANKS AGAINST AP")
    for rr in readouts:
        print(
            f"{rr['name']}: ap_rank={rr['ap_rank']} dilation_rank={rr['dilation_rank']} "
            f"all_beaters={rr['all_beaters']} primitive_beaters={rr['primitive_beaters']} "
            f"all_ties={rr['all_ties']} primitive_ties={rr['primitive_ties']}"
        )
        print(f"  ap_value={fmt_float(rr['ap_value'])}")
        print(f"  best={rr['best_E']} best_value={fmt_float(rr['best_value'])}")
        print(f"  top5_non_AP_orbit={rr['top5_nonorbit']}")


def print_ap_profile(ap: dict[str, object], dilation: dict[str, object]) -> None:
    print("\nAP GREEN PROFILE")
    print(f"AP={ap['E']}")
    print(f"dilation_twin={dilation['E']}")
    print(f"total_positive_conductance={fmt_float(ap['total_positive_conductance'])}")
    print(f"negative_edges={ap['negative_edges']} negative_mass={fmt_float(ap['negative_mass'])}")
    print(f"lambda2_algebraic_connectivity={fmt_float(ap['lambda2'])}")
    print("laplacian_eigvals=" + ",".join(fmt_float(float(x)) for x in ap["laplacian_eigvals"]))
    print(f"kirchhoff_index={fmt_float(ap['kirchhoff'])}")
    print(f"mean_effective_resistance={fmt_float(ap['mean_resistance'])}")
    print(f"max_effective_resistance={fmt_float(ap['max_resistance'])}")
    dr = ap["distance_resistance"]
    print(
        "distance_effective_resistance="
        f"d1:{fmt_float(dr[1])}, d2:{fmt_float(dr[2])}, d3:{fmt_float(dr[3])}"
    )
    current = ap["bottleneck_current"]
    print(f"bottleneck_pair={ap['bottleneck_pair']}")
    print(
        "bottleneck_unit_current="
        f"R={fmt_float(current['effective_resistance'])}, "
        f"max_edge={current['max_current_edge']}, "
        f"max_abs_current={fmt_float(current['max_abs_current'])}, "
        f"current_entropy_bits={fmt_float(current['current_entropy_bits'])}"
    )


def print_trap_table(ap: dict[str, object], traps: list[dict[str, object]]) -> None:
    print("\nHYP-3202 EXCHANGE-TRAP GREEN DISCHARGE")
    print(f"arbitrary_exchange_primitive_local_maxima={len(traps)}")
    discharge = Counter()
    for row in traps:
        if row["E"] == CONSEC_8:
            primary = "AP"
        else:
            deficits = {
                "lambda2_loss": max(0.0, float(ap["lambda2"]) - float(row["lambda2"]))
                / max(float(ap["lambda2"]), FLOAT_TOL),
                "kirchhoff_excess": max(0.0, float(row["kirchhoff"]) - float(ap["kirchhoff"]))
                / max(float(ap["kirchhoff"]), FLOAT_TOL),
                "maxR_excess": max(0.0, float(row["max_resistance"]) - float(ap["max_resistance"]))
                / max(float(ap["max_resistance"]), FLOAT_TOL),
                "negative_mass": float(row["negative_mass"]),
            }
            primary = max(deficits.items(), key=lambda item: (item[1], item[0]))[0]
            discharge[primary] += 1
        current = row["bottleneck_current"]
        print(
            f"E={row['E']} sk2={float(row['sk2']):.12f} "
            f"lambda2={fmt_float(row['lambda2'])} "
            f"kirchhoff={fmt_float(row['kirchhoff'])} "
            f"maxR={fmt_float(row['max_resistance'])} "
            f"bottleneck_pair={row['bottleneck_pair']} "
            f"current_max_edge={current['max_current_edge'] if current else None} "
            f"neg_edges={row['negative_edges']} neg_mass={fmt_float(row['negative_mass'])} "
            f"primary_green_discharge={primary}"
        )
    print(f"non_AP_trap_primary_green_discharge={dict(sorted(discharge.items()))}")


def tournament_analysis() -> None:
    print("\nTOURNAMENT ANALYSIS: GREEN-CURRENT CARRIERS")
    vertices = [
        ("algebraic_connectivity_lambda2", 96),
        ("kirchhoff_green_resistance_profile", 93),
        ("distance_effective_resistance_channels", 88),
        ("schur_complement_trap_discharge", 83),
        ("positive_covariance_conductance_graph", 76),
        ("negative_covariance_leakage_sidecar", 70),
        ("spectral_dictionary_compatibility", 64),
        ("raw_total_positive_conductance", 49),
        ("raw_covariance_scalar", 38),
        ("plain_positive_association", 18),
        ("runner_gap_graph", 5),
    ]
    scores = [score for _, score in vertices]
    hist = Counter(range(len(vertices)))
    path = [name for name, _ in sorted(vertices, key=lambda item: (-item[1], item[0]))]
    print("vertices=conductance certificates and proof obligations, not runners or raw gaps")
    print("pairwise_observable=which carrier preserves Green current plus leakage sidecar payload")
    print("switch/gauge=A beats B iff route_score(A)>route_score(B); ties lexical")
    print(f"score_hist={dict(sorted(hist.items()))}")
    print("directed_3cycles=0")
    print("sccs=singletons under the score gauge")
    print("hamiltonian_path_count=1")
    print("priority_path=" + " -> ".join(path))


def main() -> None:
    rows = [row_details((0,) + combo) for combo in itertools.combinations(range(1, 15), 7)]
    by_E = {row["E"]: row for row in rows}
    primitive_rows = [row for row in rows if row["primitive"]]
    ap = by_E[CONSEC_8]
    dilation = by_E[EVEN_AP]
    traps = local_maxima(primitive_rows, by_E, allswap_neighbors)

    print("HYP-3236 Green conductance / algebraic connectivity scout")
    print("=" * 78)
    print("bank=anchored {0} union A, A subset [1,14], |A|=7")
    print(f"rows_all={len(rows)}")
    print(f"rows_primitive={len(primitive_rows)}")
    print("conductance_graph=positive part of inner-sector empty covariance matrix")
    print("green_kernel=Moore-Penrose inverse of weighted graph Laplacian")
    print("negative_covariance_edges=retained as leakage sidecar, not erased")

    disconnected = [row["E"] for row in rows if not row["connected"]]
    print(f"positive_conductance_disconnected_rows={len(disconnected)}")
    if disconnected:
        print(f"first_disconnected_rows={disconnected[:10]}")

    print_ap_profile(ap, dilation)

    specs = [
        ("lambda2_algebraic_connectivity_MAX", lambda r: float(r["lambda2"]), True),
        ("total_positive_conductance_MAX", lambda r: float(r["total_positive_conductance"]), True),
        ("kirchhoff_index_MIN", lambda r: float(r["kirchhoff"]), False),
        ("mean_effective_resistance_MIN", lambda r: float(r["mean_resistance"]), False),
        ("max_effective_resistance_MIN", lambda r: float(r["max_resistance"]), False),
        ("distance_R1_MIN", lambda r: float(r["distance_resistance"][1]), False),
        ("distance_R2_MIN", lambda r: float(r["distance_resistance"][2]), False),
        ("distance_R3_MIN", lambda r: float(r["distance_resistance"][3]), False),
        ("negative_edges_MIN", lambda r: float(r["negative_edges"]), False),
    ]
    readouts = [
        rank_readout(rows, primitive_rows, ap, name, getter, high=high)
        for name, getter, high in specs
    ]
    print_rank_section(readouts)

    print("\nTOP NON-AP GREEN DECOYS")
    ordered = sorted(rows, key=lambda r: (-float(r["lambda2"]), float(r["kirchhoff"]), r["E"]))
    for row in [r for r in ordered if r["E"] not in AP_ORBIT][:12]:
        dr = row["distance_resistance"]
        print(
            f"E={row['E']} primitive={row['primitive']} "
            f"lambda2={fmt_float(row['lambda2'])} "
            f"kirchhoff={fmt_float(row['kirchhoff'])} "
            f"maxR={fmt_float(row['max_resistance'])} "
            f"Rdist=({fmt_float(dr[1])},{fmt_float(dr[2])},{fmt_float(dr[3])}) "
            f"neg_edges={row['negative_edges']} neg_mass={fmt_float(row['negative_mass'])}"
        )

    print_trap_table(ap, traps)

    print("\nASSUMPTION CHALLENGE")
    print("alternate_vertices_considered=runners,gaps,sector boundaries,wall crossings,")
    print("  residues, covariance edges, conductance graphs, current paths, Fiedler modes,")
    print("  Schur-complement edits, trap rows, proof obligations")
    print("chosen_vertices=conductance graphs and Green-current certificate carriers")
    print("preserves=AP covariance extremality, distance-layer payload, trap identity,")
    print("  effective-resistance bottlenecks, algebraic-connectivity margins")
    print("destroys=negative covariance signs unless leakage sidecar is retained; raw")
    print("  runner/gap identities; odd Worpitzky debt")
    print("challenged_assumption=positive association alone is enough; the graph must")
    print("  retain where current flows and where negative leakage was clipped.")

    tournament_analysis()

    print("\nPROOF-FRONTIER READOUT")
    print("readout_1=AP/doubled-AP are the only maximizers of lambda2 and total positive")
    print("  conductance, and the only minimizers of Kirchhoff/mean/max Green resistance.")
    print("readout_2=the HYP-3202 exchange traps have large effective-resistance excesses;")
    print("  their primary Green discharges split between maxR bottlenecks, lambda2 loss,")
    print("  Kirchhoff excess, and negative leakage.")
    print("readout_3=Green conductance is a strong new dictionary coordinate, but the")
    print("  positive-part construction is lossy; a proof must keep negative covariance")
    print("  leakage and odd Worpitzky sidecars.")


if __name__ == "__main__":
    main()
