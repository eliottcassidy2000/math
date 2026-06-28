#!/usr/bin/env python3
"""HYP-3238 scout: even/odd and positive/negative duality bridge.

This extends the HYP-3236 Green conductance scout on the same anchored k=8
bank.  The goal is not another scalar ranking.  The goal is to audit which
positive/even compressions are legal only after an odd/negative sidecar is
retained or discharged.

Measured channels:

* positive Green conductance: lambda2, Kirchhoff, max effective resistance;
* negative covariance leakage: negative_edges and negative_mass;
* even endpoint mass: q0 + q6;
* odd/central debt proxy: q3;
* k=8 shell dual: L_y = q0 + q6 + q3 / 10;
* HYP-3220 parity sign: de Moivre power sums have sign (-1)^k because the
  dominant period is the negative Perron root -2 cos(pi/7).

Tournament Analysis uses proof obligations / information channels as vertices,
not runners or raw arcs.
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


def green_profile(W: np.ndarray, L: np.ndarray) -> dict[str, object]:
    eig = np.linalg.eigvalsh(L)
    comps = connected_components(W)
    connected = len(comps) == 1
    out: dict[str, object] = {
        "laplacian_eigvals": eig,
        "lambda2": float(eig[1]) if len(eig) > 1 else 0.0,
        "components": comps,
        "connected": connected,
    }
    if not connected:
        out.update(
            {
                "kirchhoff": float("inf"),
                "mean_resistance": float("inf"),
                "max_resistance": float("inf"),
                "distance_resistance": {1: float("inf"), 2: float("inf"), 3: float("inf")},
            }
        )
        return out

    green = np.linalg.pinv(L, hermitian=True)
    pair_resistances = []
    by_distance: dict[int, list[float]] = {1: [], 2: [], 3: []}
    for i, j in itertools.combinations(range(6), 2):
        r = float(green[i, i] + green[j, j] - 2 * green[i, j])
        pair_resistances.append(r)
        d = abs((i + 1) - (j + 1))
        d = min(d, 7 - d)
        by_distance[d].append(r)

    out.update(
        {
            "kirchhoff": sum(pair_resistances),
            "mean_resistance": sum(pair_resistances) / len(pair_resistances),
            "max_resistance": max(pair_resistances),
            "distance_resistance": {
                d: sum(values) / len(values) for d, values in by_distance.items()
            },
        }
    )
    return out


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
    total_covariance = F(0)
    negative_edges = 0
    negative_mass = 0.0
    for i, j in itertools.combinations(INNER, 2):
        c = pair_empty.get((i, j), F(0)) - p_empty[i] * p_empty[j]
        cov[(i, j)] = c
        total_covariance += c
        C[i - 1, j - 1] = C[j - 1, i - 1] = float(c)
        if c < 0:
            negative_edges += 1
            negative_mass += -float(c)

    W = np.maximum(C, 0.0)
    np.fill_diagonal(W, 0.0)
    L = np.diag(W.sum(axis=1)) - W

    q0, q3, q6 = q[0], q[3], q[6]
    endpoint_bimodality = q0 + q6
    ly = endpoint_bimodality + q3 / 10

    return {
        "E": E,
        "primitive": is_primitive(E),
        "q": tuple(q),
        "q0": q0,
        "q3": q3,
        "q6": q6,
        "endpoint_bimodality": endpoint_bimodality,
        "Ly": ly,
        "cov": cov,
        "total_covariance": total_covariance,
        "C": C,
        "W": W,
        "L": L,
        "total_positive_conductance": float(W.sum() / 2),
        "negative_edges": negative_edges,
        "negative_mass": negative_mass,
        **green_profile(W, L),
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
        value = row["total_covariance"]
        neighbor_values = []
        for N, _ in neighbor_fn(E):
            if N not in by_E:
                continue
            if not by_E[N]["primitive"]:
                continue
            neighbor_values.append(by_E[N]["total_covariance"])
        if neighbor_values and all(v <= value for v in neighbor_values):
            out.append(row)
    return sorted(out, key=lambda r: (r["total_covariance"], r["E"]), reverse=True)


def annotate_slacks(rows: list[dict[str, object]], ap: dict[str, object]) -> None:
    for row in rows:
        row["lambda2_gap"] = float(ap["lambda2"]) - float(row["lambda2"])
        row["kirchhoff_excess"] = float(row["kirchhoff"]) - float(ap["kirchhoff"])
        row["maxR_excess"] = float(row["max_resistance"]) - float(ap["max_resistance"])
        row["bimodality_gap"] = ap["endpoint_bimodality"] - row["endpoint_bimodality"]
        row["q3_debt"] = max(F(0), row["q3"] - ap["q3"])
        row["Ly_gap"] = ap["Ly"] - row["Ly"]
        row["q3_exchange_margin"] = row["bimodality_gap"] - row["q3_debt"] / 10


def rank_float(
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
    ties = [row for row in rows if abs(getter(row) - ap_value) <= FLOAT_TOL]
    primitive_ties = [row for row in primitive_rows if abs(getter(row) - ap_value) <= FLOAT_TOL]
    return {
        "name": name,
        "ap_value": ap_value,
        "all_beaters": sum(1 for row in rows if beats(getter(row))),
        "primitive_beaters": sum(1 for row in primitive_rows if beats(getter(row))),
        "all_ties": len(ties),
        "primitive_ties": len(primitive_ties),
        "primitive_false_ties": sum(
            1 for row in primitive_ties if row["E"] != CONSEC_8
        ),
        "top5_nonorbit": [row["E"] for row in ordered if row["E"] not in AP_ORBIT][:5],
    }


def rank_fraction(
    rows: list[dict[str, object]],
    primitive_rows: list[dict[str, object]],
    ap: dict[str, object],
    name: str,
    getter: Callable[[dict[str, object]], F],
    *,
    high: bool,
) -> dict[str, object]:
    ap_value = getter(ap)

    def beats(value: F) -> bool:
        return value > ap_value if high else value < ap_value

    ordered = sorted(rows, key=getter, reverse=high)
    ties = [row for row in rows if getter(row) == ap_value]
    primitive_ties = [row for row in primitive_rows if getter(row) == ap_value]
    return {
        "name": name,
        "ap_value": ap_value,
        "all_beaters": sum(1 for row in rows if beats(getter(row))),
        "primitive_beaters": sum(1 for row in primitive_rows if beats(getter(row))),
        "all_ties": len(ties),
        "primitive_ties": len(primitive_ties),
        "primitive_false_ties": sum(
            1 for row in primitive_ties if row["E"] != CONSEC_8
        ),
        "top5_nonorbit": [row["E"] for row in ordered if row["E"] not in AP_ORBIT][:5],
    }


def fmt_float(x: float) -> str:
    if math.isinf(x):
        return "inf"
    return f"{x:.12f}"


def fmt_fraction(x: F) -> str:
    return f"{float(x):.12f} ({x})"


def de_moivre_power_sums(n: int = 8) -> list[int]:
    # Roots of x^3 + x^2 - 2x - 1, i.e. 2cos(2pi*j/7), j=1,2,3.
    p = [3, -1, 5]
    for k in range(3, n + 1):
        p.append(-p[k - 1] + 2 * p[k - 2] + p[k - 3])
    return p[1 : n + 1]


def print_parity_sign_certificate() -> None:
    sums = de_moivre_power_sums(8)
    sign_ok = all(((-1) ** k) * value > 0 for k, value in enumerate(sums, start=1))
    dominant = -2 * math.cos(math.pi / 7)
    print("\nHYP-3220 PARITY / NEGATIVE-PERRON SIGN CERTIFICATE")
    print("de_moivre_power_sums_p1_to_p8=" + ",".join(str(x) for x in sums))
    print(f"sign_pattern_is_exactly_(-1)^k={sign_ok}")
    print(f"dominant_negative_perron_period=-2cos(pi/7)={fmt_float(dominant)}")
    print("complement_pairs=(1,6),(2,5),(3,4)")
    print("interpretation=odd <=> negative and even <=> positive in the de Moivre chart")


def print_rank_readouts(rows: list[dict[str, object]], primitive_rows: list[dict[str, object]], ap: dict[str, object]) -> None:
    print("\nAP-TIGHT COORDINATES AND FALSE TERMINALS")
    readouts = [
        rank_fraction(rows, primitive_rows, ap, "L_y_MAX", lambda r: r["Ly"], high=True),
        rank_fraction(
            rows,
            primitive_rows,
            ap,
            "endpoint_bimodality_q0_plus_q6_MAX",
            lambda r: r["endpoint_bimodality"],
            high=True,
        ),
        rank_float(rows, primitive_rows, ap, "lambda2_MAX", lambda r: float(r["lambda2"]), high=True),
        rank_float(
            rows,
            primitive_rows,
            ap,
            "negative_edges_MIN",
            lambda r: float(r["negative_edges"]),
            high=False,
        ),
        rank_float(
            rows,
            primitive_rows,
            ap,
            "negative_mass_MIN",
            lambda r: float(r["negative_mass"]),
            high=False,
        ),
    ]
    for rr in readouts:
        value = rr["ap_value"]
        value_text = fmt_fraction(value) if isinstance(value, F) else fmt_float(float(value))
        print(
            f"{rr['name']}: all_beaters={rr['all_beaters']} "
            f"primitive_beaters={rr['primitive_beaters']} all_ties={rr['all_ties']} "
            f"primitive_ties={rr['primitive_ties']} "
            f"primitive_false_ties={rr['primitive_false_ties']}"
        )
        print(f"  ap_value={value_text}")
        print(f"  top5_non_AP_orbit={rr['top5_nonorbit']}")


def print_false_terminal_audit(rows: list[dict[str, object]], primitive_rows: list[dict[str, object]], ap: dict[str, object]) -> None:
    zero_neg = [row for row in primitive_rows if row["negative_edges"] == 0]
    zero_neg_nonap = [row for row in zero_neg if row["E"] != CONSEC_8]
    connected_nonap = [
        row for row in primitive_rows if row["connected"] and row["E"] != CONSEC_8
    ]
    q3_debt_rows = [row for row in primitive_rows if row["q3_debt"] > 0]
    exchange_violations = [
        row for row in q3_debt_rows if row["q3_exchange_margin"] < 0
    ]
    ly_violations = [row for row in primitive_rows if row["Ly_gap"] < 0]

    print("\nFALSE-TERMINAL AUDIT")
    print(f"primitive_zero_negative_edges_rows={len(zero_neg)}")
    print(f"primitive_zero_negative_edges_nonAP_false_terminals={len(zero_neg_nonap)}")
    print(f"primitive_connected_positive_graph_nonAP_rows={len(connected_nonap)}")
    print(f"primitive_rows_with_positive_q3_debt={len(q3_debt_rows)}")
    print(f"q3_exchange_margin_violations_for_debt_rows={len(exchange_violations)}")
    print(f"Ly_AP_max_violations={len(ly_violations)}")

    if q3_debt_rows:
        worst_ratio = min(
            (row["bimodality_gap"] / row["q3_debt"], row)
            for row in q3_debt_rows
            if row["q3_debt"] > 0
        )
        ratio, row = worst_ratio
        print(
            "worst_endpoint_bimodality_gap_per_q3_debt="
            f"{fmt_fraction(ratio)} at E={row['E']}"
        )

    print("\nZERO-NEGATIVE NON-AP EXAMPLES")
    examples = sorted(
        zero_neg_nonap,
        key=lambda r: (-float(r["lambda2"]), r["Ly_gap"], r["E"]),
    )[:12]
    for row in examples:
        print(
            f"E={row['E']} lambda2={fmt_float(float(row['lambda2']))} "
            f"lambda2_gap={fmt_float(float(row['lambda2_gap']))} "
            f"kirchhoff_excess={fmt_float(float(row['kirchhoff_excess']))} "
            f"q3_debt={fmt_fraction(row['q3_debt'])} "
            f"bimod_gap={fmt_fraction(row['bimodality_gap'])} "
            f"Ly_gap={fmt_fraction(row['Ly_gap'])}"
        )


def print_green_decoys(rows: list[dict[str, object]]) -> None:
    print("\nTOP GREEN DECOYS WITH DUALITY SIDECARS")
    ordered = sorted(rows, key=lambda r: (-float(r["lambda2"]), float(r["kirchhoff"]), r["E"]))
    for row in [r for r in ordered if r["E"] not in AP_ORBIT][:12]:
        print(
            f"E={row['E']} primitive={row['primitive']} "
            f"lambda2={fmt_float(float(row['lambda2']))} "
            f"kirchhoff={fmt_float(float(row['kirchhoff']))} "
            f"neg_edges={row['negative_edges']} neg_mass={fmt_float(float(row['negative_mass']))} "
            f"q3_debt={fmt_fraction(row['q3_debt'])} "
            f"bimod_gap={fmt_fraction(row['bimodality_gap'])} "
            f"Ly_gap={fmt_fraction(row['Ly_gap'])}"
        )


def print_trap_audit(traps: list[dict[str, object]], ap: dict[str, object]) -> None:
    print("\nHYP-3202 TRAPS THROUGH THE DUALITY BRIDGE")
    print(f"arbitrary_exchange_primitive_local_maxima={len(traps)}")
    classes = Counter()
    for row in traps:
        if row["E"] == CONSEC_8:
            primary = "AP"
        elif row["negative_edges"] == 0 and row["q3_debt"] == 0:
            primary = "positive_even_but_green_resistance_excess"
        elif row["negative_edges"] == 0:
            primary = "odd_q3_debt_without_negative_covariance"
        elif row["q3_debt"] > 0:
            primary = "negative_leakage_plus_odd_q3_debt"
        else:
            primary = "negative_leakage_green_bottleneck"
        classes[primary] += row["E"] != CONSEC_8
        print(
            f"E={row['E']} lambda2_gap={fmt_float(float(row['lambda2_gap']))} "
            f"kirchhoff_excess={fmt_float(float(row['kirchhoff_excess']))} "
            f"maxR_excess={fmt_float(float(row['maxR_excess']))} "
            f"neg_edges={row['negative_edges']} neg_mass={fmt_float(float(row['negative_mass']))} "
            f"q3_debt={fmt_fraction(row['q3_debt'])} "
            f"bimod_gap={fmt_fraction(row['bimodality_gap'])} "
            f"Ly_gap={fmt_fraction(row['Ly_gap'])} primary_bridge_class={primary}"
        )
    print(f"non_AP_trap_bridge_classes={dict(sorted(classes.items()))}")


def tournament_analysis(rows: list[dict[str, object]], primitive_rows: list[dict[str, object]]) -> None:
    zero_neg_false = sum(
        1 for row in primitive_rows if row["negative_edges"] == 0 and row["E"] != CONSEC_8
    )
    connected_false = sum(
        1 for row in primitive_rows if row["connected"] and row["E"] != CONSEC_8
    )
    vertices = [
        ("even_positive_fejer_square", 96, "SOS magnitude side"),
        ("green_dirichlet_positive_face", 92, "AP-tight but clips signs"),
        ("odd_negative_brouwer_sign", 89, "HYP-3220/HYP-3219 sign sidecar"),
        ("hermite_biehler_interlacing", 85, "even/odd gluing sidecar"),
        ("bulk_core_vitali_transfer", 82, "measure/core information wall"),
        ("pair_pascal_cap_mass", 78, "cap magnitude coordinate"),
        ("toeplitz_normal_fan_slack", 74, "finite trap discharge chart"),
        ("negative_covariance_leakage", 64, f"false_terminal_if_absent={zero_neg_false}"),
        ("signed_chart_change_debt", 61, "recursion transport sidecar"),
        ("raw_lambda2_scalar", 34, "AP-tight finite scalar but not transport data"),
        ("raw_connected_positive_graph", 20, f"false_terminals={connected_false}"),
        ("raw_positive_association", 8, f"false_terminals={zero_neg_false}"),
    ]
    ordered = sorted(vertices, key=lambda item: (-item[1], item[0]))
    outdegrees = {name: i for i, (name, _, _) in enumerate(reversed(ordered))}
    hist = Counter(outdegrees.values())
    print("\nTOURNAMENT ANALYSIS: DUALITY CARRIERS")
    print("vertices=proof obligations / retained information channels, not runners")
    print("pairwise_observable=which carrier preserves LRC predicate with less lost duality payload")
    print("switch/gauge=A beats B iff it retains/discharges more crossed even/odd positive/negative payload")
    print(f"score_hist={dict(sorted(hist.items()))}")
    print("directed_3cycles=0 under this retained-payload gauge")
    print("sccs=singletons")
    print("hamiltonian_path_count=1")
    print("priority_path=" + " -> ".join(name for name, _, _ in ordered))
    for name, score, note in ordered:
        print(f"  vertex={name} score={score} note={note}")


def print_assumption_challenge() -> None:
    print("\nASSUMPTION CHALLENGE")
    print("alternate_vertices_considered=runners,gaps,fixed circle sections,section")
    print("  boundaries, wall crossings, residues, cover arcs, Fourier modes, Fiedler")
    print("  modes, negative covariance edges, Brouwer signs, Fejer coefficients,")
    print("  matroid circuits, proof obligations")
    print("chosen_vertices=proof obligations / retained information channels")
    print("preserves=crossed even/odd and positive/negative duality payload plus AP")
    print("  extremality predicate")
    print("destroys=row identity as a primary tournament vertex; exemplar rows and trap")
    print("  sidecars are retained separately in the table")
    print("challenged_assumption=the positive/even side is the proof and the")
    print("  odd/negative side is error; better is magnitude certificate plus sign/core")
    print("  sidecar")


def main() -> None:
    rows = [row_details((0,) + combo) for combo in itertools.combinations(range(1, 15), 7)]
    by_E = {row["E"]: row for row in rows}
    primitive_rows = [row for row in rows if row["primitive"]]
    ap = by_E[CONSEC_8]
    dilation = by_E[EVEN_AP]
    annotate_slacks(rows, ap)
    traps = local_maxima(primitive_rows, by_E, allswap_neighbors)

    print("HYP-3238 even/odd positive/negative duality bridge scout")
    print("=" * 78)
    print("bank=anchored {0} union A, A subset [1,14], |A|=7")
    print(f"rows_all={len(rows)}")
    print(f"rows_primitive={len(primitive_rows)}")
    print("duality_claim=even/positive compression requires odd/negative sidecar")
    print("HYP_3220_integration=even-odd duality is positive-negative duality")
    print(f"AP={ap['E']}")
    print(f"dilation_twin={dilation['E']}")
    print(f"AP_q0={fmt_fraction(ap['q0'])}")
    print(f"AP_q3={fmt_fraction(ap['q3'])}")
    print(f"AP_q6={fmt_fraction(ap['q6'])}")
    print(f"AP_endpoint_bimodality_q0_plus_q6={fmt_fraction(ap['endpoint_bimodality'])}")
    print(f"AP_L_y={fmt_fraction(ap['Ly'])}")
    print(f"AP_lambda2={fmt_float(float(ap['lambda2']))}")
    print(f"AP_kirchhoff={fmt_float(float(ap['kirchhoff']))}")
    print(f"AP_negative_edges={ap['negative_edges']} AP_negative_mass={fmt_float(float(ap['negative_mass']))}")

    print_parity_sign_certificate()
    print_rank_readouts(rows, primitive_rows, ap)
    print_false_terminal_audit(rows, primitive_rows, ap)
    print_green_decoys(rows)
    print_trap_audit(traps, ap)
    print_assumption_challenge()
    tournament_analysis(rows, primitive_rows)

    print("\nPROOF-FRONTIER READOUT")
    print("readout_1=HYP-3220 upgrades HYP-3238: even/odd and positive/negative")
    print("  are one parity/complement sign operator in the de Moivre chart.")
    print("readout_2=zero negative covariance is a false terminal condition: non-AP")
    print("  primitive rows can have no negative leakage but still lose Green and L_y.")
    print("readout_3=q3 central debt is priced by endpoint bimodality in the L_y")
    print("  exchange check; any stronger proof should glue this to HB/Brouwer sign.")
    print("readout_4=the proof packet should be magnitude certificate plus sign/core")
    print("  sidecar, not a positive/even scalar alone.")


if __name__ == "__main__":
    main()
