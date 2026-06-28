#!/usr/bin/env python3
"""HYP-3205 scout: spectral dictionary compatibility for the k=8 frontier.

This is a synthesis scout rather than a standalone proof.  It places the
current k=8 LRC14 certificates in one dictionary:

* covariance and cyclic-distance covariance layers (HYP-3202),
* AP-polarized cyclotomic support and orbit-aware compression (HYP-3203),
* Caratheodory-Toeplitz moment-cone margin (KPS S31al / S73d),
* Perron uniform-mode diagnostics and the Joukowski/Hermite-Biehler bridge
  (HYP-3210).

Tournament Analysis uses proof channels as vertices, not runners, sectors, or
arcs.  The pairwise observable is certificate payload retained for the LRC
coverage/covariance extremality.  The switch orients toward the channel with
fewer AP beaters, fewer nontrivial ties, more localization of the obstruction,
and better sidecar retention.
"""

from __future__ import annotations

import itertools
import math
from collections import Counter
from fractions import Fraction as F
from typing import Callable, Iterable, Sequence

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

    out: list[tuple[set[int], F]] = []
    pts = sorted(breakpoints)
    for x0, x1 in zip(pts, pts[1:]):
        if x1 <= x0:
            continue
        mid = (x0 + x1) / 2
        out.append(({sector_of(e * mid) for e in E}, x1 - x0))
    return out


def cyclotomic_energy(q: Sequence[F]) -> F:
    # Parseval for a 7-coefficient probability vector, excluding the DC term.
    return 7 * sum(x * x for x in q) - 1


def centered(q: Sequence[F]) -> tuple[F, ...]:
    return tuple(x - F(1, 7) for x in q)


def dot(a: Sequence[F], b: Sequence[F]) -> F:
    return sum(x * y for x, y in zip(a, b))


def toeplitz_min_eigenvalue(q: Sequence[F]) -> float:
    qf = np.array([float(x) for x in q], dtype=float)
    T = np.array([[qf[abs(i - j)] for j in range(7)] for i in range(7)], dtype=float)
    return float(np.linalg.eigvalsh(T)[0])


def row_details(E: tuple[int, ...], ap_direction: Sequence[F] | None = None) -> dict[str, object]:
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
    cov_matrix = np.zeros((6, 6), dtype=float)
    distance_profile = {1: F(0), 2: F(0), 3: F(0)}
    sk2 = F(0)
    min_cov: F | None = None
    neg_pair_count = 0
    for i, j in itertools.combinations(INNER, 2):
        c = pair_empty.get((i, j), F(0)) - p_empty[i] * p_empty[j]
        cov[(i, j)] = c
        cov_matrix[i - 1, j - 1] = cov_matrix[j - 1, i - 1] = float(c)
        sk2 += c
        d = abs(i - j)
        d = min(d, 7 - d)
        distance_profile[d] += c
        if c < 0:
            neg_pair_count += 1
        if min_cov is None or c < min_cov:
            min_cov = c
    for i in INNER:
        var = p_empty[i] * (1 - p_empty[i])
        cov_matrix[i - 1, i - 1] = float(var)

    eigvals, eigvecs = np.linalg.eigh(cov_matrix)
    perron_vec = eigvecs[:, int(np.argmax(eigvals))]
    ones = np.ones(6, dtype=float)
    perron_align = abs(float(np.dot(perron_vec, ones) / (np.linalg.norm(perron_vec) * np.linalg.norm(ones))))
    q_tuple = tuple(q)
    d = centered(q_tuple)

    return {
        "E": E,
        "primitive": is_primitive(E),
        "q": q_tuple,
        "q0": q_tuple[0],
        "Ly": q_tuple[0] + q_tuple[6] + q_tuple[3] / 10,
        "cov": cov,
        "sk2": sk2,
        "distance_profile": distance_profile,
        "min_cov": min_cov,
        "neg_pair_count": neg_pair_count,
        "cyclotomic_energy": cyclotomic_energy(q_tuple),
        "ap_projection": dot(d, ap_direction) if ap_direction is not None else F(0),
        "toeplitz_min": toeplitz_min_eigenvalue(q_tuple),
        "cov_lambda_max": float(eigvals[-1]),
        "cov_lambda_min": float(eigvals[0]),
        "perron_uniform_alignment": perron_align,
    }


def allswap_neighbors(E: tuple[int, ...]):
    A = set(E[1:])
    for x in sorted(A):
        for y in sorted(set(range(1, 15)) - A):
            yield (0,) + tuple(sorted((A - {x}) | {y})), (x, y)


def left_moves(E: tuple[int, ...]):
    S = set(E)
    for e in sorted(S, reverse=True):
        if e == 0:
            continue
        for m in range(1, e):
            if m not in S:
                yield tuple(sorted((S - {e}) | {m})), (e, m)


def local_maxima(
    rows: list[dict[str, object]],
    by_E: dict[tuple[int, ...], dict[str, object]],
    neighbor_fn: Callable[[tuple[int, ...]], Iterable[tuple[tuple[int, ...], tuple[int, int]]]],
    *,
    primitive_neighbors: bool,
) -> list[dict[str, object]]:
    out = []
    for row in rows:
        E = row["E"]
        value = row["sk2"]
        neighbor_values = []
        for N, _ in neighbor_fn(E):
            if N not in by_E:
                continue
            if primitive_neighbors and not by_E[N]["primitive"]:
                continue
            neighbor_values.append(by_E[N]["sk2"])
        if not neighbor_values or all(v <= value for v in neighbor_values):
            out.append(row)
    return sorted(out, key=lambda r: (r["sk2"], r["E"]), reverse=True)


def fmt_value(value: object) -> str:
    if isinstance(value, F):
        return f"{value} ({float(value):+.8f})"
    if isinstance(value, float):
        return f"{value:+.10f}"
    return str(value)


def exact_equal(a: object, b: object) -> bool:
    if isinstance(a, float) or isinstance(b, float):
        return abs(float(a) - float(b)) <= FLOAT_TOL
    return a == b


def ap_deficit(row: dict[str, object], ap: dict[str, object], key: str) -> float:
    return max(0.0, float(ap[key]) - float(row[key]))


def layer_deficit(row: dict[str, object], ap: dict[str, object], d: int) -> float:
    return max(0.0, float(ap["distance_profile"][d]) - float(row["distance_profile"][d]))


def rank_readout(
    rows: list[dict[str, object]],
    primitive_rows: list[dict[str, object]],
    ap: dict[str, object],
    name: str,
    getter: Callable[[dict[str, object]], object],
    *,
    high: bool = True,
) -> dict[str, object]:
    ap_value = getter(ap)

    def beats(value: object) -> bool:
        if high:
            return float(value) > float(ap_value) + FLOAT_TOL
        return float(value) < float(ap_value) - FLOAT_TOL

    def ties(value: object) -> bool:
        return exact_equal(value, ap_value)

    ordered = sorted(rows, key=lambda row: float(getter(row)), reverse=high)
    ranks = {row["E"]: i for i, row in enumerate(ordered)}
    all_beaters = [row["E"] for row in rows if beats(getter(row))]
    all_ties = [row["E"] for row in rows if ties(getter(row))]
    primitive_beaters = [row["E"] for row in primitive_rows if beats(getter(row))]
    primitive_ties = [row["E"] for row in primitive_rows if ties(getter(row))]
    best = ordered[0]
    return {
        "name": name,
        "ap_value": ap_value,
        "best_E": best["E"],
        "best_value": getter(best),
        "ap_rank": ranks[CONSEC_8],
        "dilation_rank": ranks[EVEN_AP],
        "all_beaters": len(all_beaters),
        "all_ties": len(all_ties),
        "primitive_beaters": len(primitive_beaters),
        "primitive_ties": len(primitive_ties),
        "top5_nonorbit": [row["E"] for row in ordered if row["E"] not in AP_ORBIT][:5],
    }


def channel_tournament(
    exact_readouts: dict[str, dict[str, object]],
    left_trap_count: int,
    exchange_trap_count: int,
) -> None:
    def strength(name: str, route_bonus: int) -> int:
        rr = exact_readouts[name]
        return (
            route_bonus
            + (40 if rr["all_beaters"] == 0 else -3 * int(rr["all_beaters"]))
            + (20 if rr["primitive_beaters"] == 0 else -2 * int(rr["primitive_beaters"]))
            - 4 * max(0, int(rr["all_ties"]) - 2)
            - 5 * max(0, int(rr["primitive_ties"]) - 1)
        )

    scores = {
        "cyclic_distance_layer_bundle": min(
            strength("distance_1", 14), strength("distance_2", 14), strength("distance_3", 14)
        )
        + 5,
        "AP_residual_support_hyperplane": strength("AP_residual_projection", 16),
        "Toeplitz_moment_margin": strength("Toeplitz_lambda_min", 15),
        "total_covariance_scalar": strength("Sigma_kappa2", 9),
        "exchange_trap_sheaf": 48 - exchange_trap_count,
        "orbit_aware_left_compression": 38 - left_trap_count,
        "Perron_uniform_mode_diagnostic": strength("Perron_uniform_alignment", 5),
        "raw_cyclotomic_norm": strength("raw_cyclotomic_energy_MIN", -5),
        "entropy_or_min_description_scalar": -11,
        "one_seventh_associator_scalar": -19,
    }
    ordered = sorted(scores.items(), key=lambda item: (-item[1], item[0]))
    hist = Counter(scores.values())
    print("\nTOURNAMENT ANALYSIS: PROOF-CHANNEL DICTIONARY")
    print("vertices=proof channels/certificates, not runners, gaps, arcs, roots, or constants")
    print("pairwise_observable=retained certificate payload for k=8 covariance/coverage extremality")
    print("switch/gauge=A->B iff channel_score(A)>channel_score(B), lexical tie break")
    print(f"score_hist={dict(sorted(hist.items()))}")
    print("directed_3cycles=0")
    print("sccs=singletons under the score gauge")
    print("hamiltonian_path_count=1")
    print("priority_path=" + " -> ".join(name for name, _ in ordered))


def main() -> None:
    q_ap = tuple(x["q"] for x in [row_details(CONSEC_8)])[0]
    ap_direction = centered(q_ap)
    rows = [
        row_details((0,) + combo, ap_direction)
        for combo in itertools.combinations(range(1, 15), 7)
    ]
    by_E = {row["E"]: row for row in rows}
    primitive_rows = [row for row in rows if row["primitive"]]
    ap = by_E[CONSEC_8]
    dilation = by_E[EVEN_AP]

    print("HYP-3205 k=8 spectral dictionary compatibility scout")
    print("=" * 78)
    print("bank=anchored {0} union A, A subset [1,14], |A|=7")
    print(f"rows_all={len(rows)}")
    print(f"rows_primitive={len(primitive_rows)}")
    print(f"AP={CONSEC_8}")
    print(f"dilation_twin={EVEN_AP}")
    print(f"AP_q={ap['q']}")
    print(f"dilation_q={dilation['q']}")

    readout_specs = [
        ("coverage_q0", lambda r: r["q0"], True),
        ("Ly", lambda r: r["Ly"], True),
        ("Sigma_kappa2", lambda r: r["sk2"], True),
        ("distance_1", lambda r: r["distance_profile"][1], True),
        ("distance_2", lambda r: r["distance_profile"][2], True),
        ("distance_3", lambda r: r["distance_profile"][3], True),
        ("AP_residual_projection", lambda r: r["ap_projection"], True),
        ("Toeplitz_lambda_min", lambda r: r["toeplitz_min"], True),
        ("Perron_uniform_alignment", lambda r: r["perron_uniform_alignment"], True),
        ("raw_cyclotomic_energy_MIN", lambda r: r["cyclotomic_energy"], False),
    ]

    readouts = {
        name: rank_readout(rows, primitive_rows, ap, name, getter, high=high)
        for name, getter, high in readout_specs
    }

    print("\nCHANNEL RANKS AGAINST AP")
    for name, _, _ in readout_specs:
        rr = readouts[name]
        print(
            f"{name}: ap_rank={rr['ap_rank']} dilation_rank={rr['dilation_rank']} "
            f"all_beaters={rr['all_beaters']} all_ties={rr['all_ties']} "
            f"primitive_beaters={rr['primitive_beaters']} primitive_ties={rr['primitive_ties']}"
        )
        print(f"  ap_value={fmt_value(rr['ap_value'])}")
        print(f"  best={rr['best_E']} best_value={fmt_value(rr['best_value'])}")
        print(f"  top5_non_AP_orbit={rr['top5_nonorbit']}")

    exact_certificate_names = [
        "q0",
        "Ly",
        "sk2",
        "ap_projection",
    ]
    exact_layer_names = [1, 2, 3]
    simultaneous_exact = []
    simultaneous_with_toeplitz = []
    for row in rows:
        exact_ok = all(row[name] == ap[name] for name in exact_certificate_names)
        exact_ok = exact_ok and all(
            row["distance_profile"][d] == ap["distance_profile"][d] for d in exact_layer_names
        )
        if exact_ok:
            simultaneous_exact.append(row["E"])
            if abs(float(row["toeplitz_min"]) - float(ap["toeplitz_min"])) <= FLOAT_TOL:
                simultaneous_with_toeplitz.append(row["E"])

    print("\nCOMPATIBILITY READOUT")
    print(f"simultaneous_exact_maximizers={simultaneous_exact}")
    print(f"simultaneous_exact_plus_toeplitz_maximizers={simultaneous_with_toeplitz}")
    print("dictionary_claim_candidate=AP and its 2x dilation are the only rows coherent")
    print("across coverage, Ly, covariance, three distance layers, AP support, and Toeplitz margin.")

    certs = [
        ("q0", lambda r: float(ap["q0"] - r["q0"])),
        ("Ly", lambda r: float(ap["Ly"] - r["Ly"])),
        ("sk2", lambda r: float(ap["sk2"] - r["sk2"])),
        ("d1", lambda r: float(ap["distance_profile"][1] - r["distance_profile"][1])),
        ("d2", lambda r: float(ap["distance_profile"][2] - r["distance_profile"][2])),
        ("d3", lambda r: float(ap["distance_profile"][3] - r["distance_profile"][3])),
        ("ap_support", lambda r: float(ap["ap_projection"] - r["ap_projection"])),
        ("toeplitz", lambda r: float(ap["toeplitz_min"] - r["toeplitz_min"])),
    ]
    scale = {name: max(fn(row) for row in rows) or 1.0 for name, fn in certs}
    danger_rows = []
    for row in rows:
        if row["E"] in AP_ORBIT:
            continue
        deficits = {name: max(0.0, fn(row)) for name, fn in certs}
        normalized = {name: deficits[name] / scale[name] for name, _ in certs}
        danger_rows.append(
            {
                "row": row,
                "deficits": deficits,
                "normalized": normalized,
                "max_norm": max(normalized.values()),
                "sum_norm": sum(normalized.values()),
                "zeroish": sum(1 for value in deficits.values() if value <= FLOAT_TOL),
            }
        )
    danger_rows.sort(key=lambda item: (item["max_norm"], item["sum_norm"], item["row"]["E"]))

    print("\nNEAREST NON-AP DICTIONARY DECOYS")
    print("criterion=smallest maximum normalized deficit across q0,Ly,sk2,D1,D2,D3,AP-support,Toeplitz")
    for item in danger_rows[:12]:
        row = item["row"]
        deficits = item["deficits"]
        print(
            f"E={row['E']} primitive={row['primitive']} "
            f"max_norm={item['max_norm']:.6f} sum_norm={item['sum_norm']:.6f} "
            f"zeroish_channels={item['zeroish']} "
            f"def_sk2={deficits['sk2']:.10f} def_ap_support={deficits['ap_support']:.10f} "
            f"def_toeplitz={deficits['toeplitz']:.10f} "
            f"layer_defs=({deficits['d1']:.10f},{deficits['d2']:.10f},{deficits['d3']:.10f})"
        )

    swap_locals = local_maxima(primitive_rows, by_E, allswap_neighbors, primitive_neighbors=True)
    left_locals = local_maxima(rows, by_E, left_moves, primitive_neighbors=False)

    print("\nTRAP-SHEAF DISCHARGE TABLE")
    print(f"arbitrary_exchange_primitive_local_maxima={len(swap_locals)}")
    exchange_discharge = Counter()
    for row in swap_locals:
        layer_defs = tuple(layer_deficit(row, ap, d) for d in (1, 2, 3))
        if row["E"] != CONSEC_8:
            deficits = {
                name: max(0.0, fn(row)) / scale[name]
                for name, fn in certs
            }
            exchange_discharge[max(deficits.items(), key=lambda item: (item[1], item[0]))[0]] += 1
        print(
            f"E={row['E']} sk2={fmt_value(row['sk2'])} "
            f"def_ap_support={ap_deficit(row, ap, 'ap_projection'):.10f} "
            f"def_toeplitz={ap_deficit(row, ap, 'toeplitz_min'):.10f} "
            f"layer_defs=({layer_defs[0]:.10f},{layer_defs[1]:.10f},{layer_defs[2]:.10f}) "
            f"neg_pairs={row['neg_pair_count']} min_cov={row['min_cov']}"
        )
    print(f"naive_left_compression_local_maxima={len(left_locals)}")
    print(f"naive_left_compression_first_12={[row['E'] for row in left_locals[:12]]}")
    left_discharge = Counter()
    for row in left_locals:
        if row["E"] == CONSEC_8:
            continue
        deficits = {
            name: max(0.0, fn(row)) / scale[name]
            for name, fn in certs
        }
        left_discharge[max(deficits.items(), key=lambda item: (item[1], item[0]))[0]] += 1
    print(f"exchange_trap_primary_discharge={dict(sorted(exchange_discharge.items()))}")
    print(f"left_compression_non_AP_primary_discharge={dict(sorted(left_discharge.items()))}")

    print("\nASSUMPTION CHALLENGE")
    print("alternate_vertices_considered=runners,gaps,sections,boundaries,wall_crossings,residues,")
    print("  cover_arcs,Fourier_modes,matroid_circuits,proof_obligations,certificate_channels")
    print("chosen_vertices=certificate_channels plus bounded-bank rows as test objects")
    print("preserves=LRC coverage/covariance extremality, dilation exception, trap identities,")
    print("  root/moment/covariance sidecar coordinates")
    print("forgets=literal runner identities, raw arc orientations, and any scalar-only proof story")
    print("challenged_assumption=the obstruction is not a single max scalar; it is compatibility")
    print("  failure across adjoint certificate languages.")

    channel_tournament(readouts, len(left_locals) - 1, len(swap_locals) - 1)

    print("\nPROOF-FRONTIER READOUT")
    print("new_target_1=prove a certificate-Helly lemma: every primitive non-AP row violates")
    print("  at least one of the layer/AP-support/Toeplitz certificates by a usable margin.")
    print("new_target_2=treat the exchange and compression traps as sheaf stalks; discharge")
    print("  each stalk by the first dictionary coordinate where it has visible deficit.")
    print("new_target_3=transport Toeplitz moment-cone slack through the Joukowski/Hermite-Biehler")
    print("  bridge to an even/odd interlacing statement, with Worpitzky/associator debt named.")


if __name__ == "__main__":
    main()
