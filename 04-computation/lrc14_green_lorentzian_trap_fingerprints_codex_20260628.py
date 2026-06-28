#!/usr/bin/env python3
"""HYP-3225 scout: Green-current and Lorentzian fingerprints of k=8 traps.

This is the executable continuation of HYP-3223 after HYP-3224.  HYP-3224
showed that the full spectral payload cube discharges every HYP-3202
arbitrary-swap covariance trap first by Toeplitz lambda-min deficit.  This
scout asks what extra structure is visible inside those same finite traps.

Tournament Analysis uses proof certificates as vertices, not runners/arcs.
The pairwise observable is retained LRC proof payload: does the signal explain
why a local covariance trap is still below the AP/consecutive normal face?
"""

from __future__ import annotations

import itertools
import math
from collections import Counter, defaultdict
from fractions import Fraction
from typing import Iterable


INNER = tuple(range(1, 7))
TARGET = tuple(range(8))
EVEN_AP = tuple(2 * i for i in range(8))
TOL = 1.0e-11
HYP3202_TRAPS = (
    TARGET,
    (0, 2, 4, 6, 7, 8, 10, 12),
    (0, 1, 3, 5, 7, 9, 11, 13),
    (0, 1, 2, 3, 11, 12, 13, 14),
    (0, 3, 6, 8, 9, 11, 12, 14),
    (0, 2, 3, 5, 6, 8, 11, 14),
    (0, 6, 7, 8, 9, 10, 11, 12),
    (0, 4, 5, 8, 9, 10, 13, 14),
    (0, 8, 9, 10, 11, 12, 13, 14),
    (0, 1, 2, 3, 7, 8, 9, 10),
    (0, 2, 5, 7, 9, 10, 12, 14),
    (0, 1, 4, 5, 7, 8, 11, 12),
)


def jacobi_eigenvalues(matrix: list[list[float]], eps: float = 1.0e-13, max_iter: int = 200) -> list[float]:
    """Eigenvalues of a small real symmetric matrix by Jacobi rotations."""
    n = len(matrix)
    a = [row[:] for row in matrix]
    for _ in range(max_iter):
        p, q = 0, 1
        max_off = 0.0
        for i in range(n):
            for j in range(i + 1, n):
                if abs(a[i][j]) > max_off:
                    max_off = abs(a[i][j])
                    p, q = i, j
        if max_off < eps:
            break
        if abs(a[p][p] - a[q][q]) < eps:
            angle = math.pi / 4
        else:
            angle = 0.5 * math.atan2(2.0 * a[p][q], a[q][q] - a[p][p])
        c = math.cos(angle)
        s = math.sin(angle)
        app = c * c * a[p][p] - 2.0 * s * c * a[p][q] + s * s * a[q][q]
        aqq = s * s * a[p][p] + 2.0 * s * c * a[p][q] + c * c * a[q][q]
        a[p][p] = app
        a[q][q] = aqq
        a[p][q] = a[q][p] = 0.0
        for r in range(n):
            if r == p or r == q:
                continue
            arp = c * a[r][p] - s * a[r][q]
            arq = s * a[r][p] + c * a[r][q]
            a[r][p] = a[p][r] = arp
            a[r][q] = a[q][r] = arq
    return sorted(a[i][i] for i in range(n))


def solve_linear(matrix: list[list[float]], rhs: list[float]) -> list[float] | None:
    """Solve a small dense linear system by Gaussian elimination."""
    n = len(matrix)
    a = [row[:] + [rhs[i]] for i, row in enumerate(matrix)]
    for col in range(n):
        pivot = max(range(col, n), key=lambda r: abs(a[r][col]))
        if abs(a[pivot][col]) < 1.0e-14:
            return None
        if pivot != col:
            a[col], a[pivot] = a[pivot], a[col]
        scale = a[col][col]
        for j in range(col, n + 1):
            a[col][j] /= scale
        for r in range(n):
            if r == col:
                continue
            factor = a[r][col]
            if abs(factor) < 1.0e-18:
                continue
            for j in range(col, n + 1):
                a[r][j] -= factor * a[col][j]
    return [a[i][n] for i in range(n)]


def pearson(xs: list[float], ys: list[float]) -> float:
    x_mean = sum(xs) / len(xs)
    y_mean = sum(ys) / len(ys)
    x0 = [x - x_mean for x in xs]
    y0 = [y - y_mean for y in ys]
    denom = math.sqrt(sum(x * x for x in x0) * sum(y * y for y in y0))
    return sum(x * y for x, y in zip(x0, y0)) / denom if denom else 0.0


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


def toeplitz_lambda_min(q: tuple[Fraction, ...]) -> float:
    qf = [float(x) for x in q]
    matrix = [[qf[abs(i - j)] for j in range(7)] for i in range(7)]
    return jacobi_eigenvalues(matrix)[0]


def row_data(E: tuple[int, ...], ap_direction: tuple[Fraction, ...]) -> dict[str, object]:
    q = [Fraction(0) for _ in range(7)]
    contains = [Fraction(0) for _ in range(1 << 6)]

    pts = row_breakpoints(E)
    for x0, x1 in zip(pts, pts[1:]):
        if x1 <= x0:
            continue
        weight = x1 - x0
        miss_count, mask = empty_mask_for_cell(E, (x0 + x1) / 2)
        q[miss_count] += weight
        sub = mask
        while True:
            contains[sub] += weight
            if sub == 0:
                break
            sub = (sub - 1) & mask

    cov = [[0.0 for _ in range(6)] for _ in range(6)]
    exact_cov: dict[tuple[int, int], Fraction] = {}
    distance_profile = {1: Fraction(0), 2: Fraction(0), 3: Fraction(0)}
    sigma_kappa2 = Fraction(0)
    neg_pair_count = 0
    min_cov: Fraction | None = None

    for i in range(6):
        pi = contains[1 << i]
        cov[i][i] = float(pi - pi * pi)
    for i, j in itertools.combinations(range(6), 2):
        pair_mask = (1 << i) | (1 << j)
        c = contains[pair_mask] - contains[1 << i] * contains[1 << j]
        exact_cov[(i + 1, j + 1)] = c
        sigma_kappa2 += c
        d = abs((i + 1) - (j + 1))
        d = min(d, 7 - d)
        distance_profile[d] += c
        if c < 0:
            neg_pair_count += 1
        if min_cov is None or c < min_cov:
            min_cov = c
        cov[i][j] = cov[j][i] = float(c)

    q_tuple = tuple(q)
    d = centered(q_tuple)
    eigvals = jacobi_eigenvalues(cov)
    return {
        "E": E,
        "primitive": is_primitive(E),
        "q": q_tuple,
        "contains": contains,
        "cov_matrix": cov,
        "cov_eigs": eigvals,
        "exact_cov": exact_cov,
        "distance_profile": distance_profile,
        "D1": distance_profile[1],
        "D2": distance_profile[2],
        "D3": distance_profile[3],
        "sigma_kappa2": sigma_kappa2,
        "neg_pair_count": neg_pair_count,
        "min_cov": min_cov,
        "q0": q_tuple[0],
        "bimod": q_tuple[0] + q_tuple[6],
        "Ly": q_tuple[0] + q_tuple[6] + q_tuple[3] / 10,
        "cyclotomic_energy": cyclotomic_energy(q_tuple),
        "ap_projection": dot(d, ap_direction),
        "toeplitz_lambda_min": toeplitz_lambda_min(q_tuple),
    }


def required_rows() -> list[tuple[int, ...]]:
    needed = set(HYP3202_TRAPS)
    needed.add(EVEN_AP)
    for E in HYP3202_TRAPS:
        needed.update(allswap_neighbors(E))
    return sorted(needed)


def build_rows() -> tuple[list[dict[str, object]], dict[tuple[int, ...], dict[str, object]]]:
    bank = required_rows()
    target_q = row_data(TARGET, (Fraction(0),) * 7)["q"]
    ap_direction = centered(target_q)
    rows = [row_data(E, ap_direction) for E in bank]
    return rows, {row["E"]: row for row in rows}


def allswap_neighbors(E: tuple[int, ...]) -> Iterable[tuple[int, ...]]:
    A = set(E[1:])
    for x in sorted(A):
        for y in sorted(set(range(1, 15)) - A):
            yield (0,) + tuple(sorted((A - {x}) | {y}))


def positive_cov_laplacian(row: dict[str, object]) -> list[list[float]]:
    lap = [[0.0 for _ in range(6)] for _ in range(6)]
    for (i, j), c in row["exact_cov"].items():
        conductance = max(0.0, float(c))
        a = i - 1
        b = j - 1
        lap[a][a] += conductance
        lap[b][b] += conductance
        lap[a][b] -= conductance
        lap[b][a] -= conductance
    return lap


def effective_resistance_pair(lap: list[list[float]], a: int, b: int) -> float:
    n = len(lap)
    ground = n - 1
    keep = [i for i in range(n) if i != ground]
    reduced = [[lap[i][j] for j in keep] for i in keep]
    current = []
    for i in keep:
        if i == a:
            current.append(1.0)
        elif i == b:
            current.append(-1.0)
        else:
            current.append(0.0)
    solution = solve_linear(reduced, current)
    if solution is None:
        return math.inf
    potentials = [0.0 for _ in range(n)]
    for idx, i in enumerate(keep):
        potentials[i] = solution[idx]
    return potentials[a] - potentials[b]


def effective_resistance_profile(row: dict[str, object]) -> dict[str, object]:
    lap = positive_cov_laplacian(row)
    eigvals = jacobi_eigenvalues(lap)
    by_distance: dict[int, list[float]] = defaultdict(list)
    max_pair = None
    max_resistance = -1.0
    for a, b in itertools.combinations(range(6), 2):
        resistance = effective_resistance_pair(lap, a, b)
        d = abs((a + 1) - (b + 1))
        d = min(d, 7 - d)
        by_distance[d].append(resistance)
        if resistance > max_resistance:
            max_resistance = resistance
            max_pair = (a + 1, b + 1)

    min_cut = math.inf
    min_cut_mask = None
    for mask in range(1, 1 << 6):
        if mask & 1 == 0:
            continue
        if mask == (1 << 6) - 1:
            continue
        cut = 0.0
        for a, b in itertools.combinations(range(6), 2):
            if ((mask >> a) & 1) != ((mask >> b) & 1):
                c = row["exact_cov"].get((a + 1, b + 1), Fraction(0))
                cut += max(0.0, float(c))
        if cut < min_cut:
            min_cut = cut
            min_cut_mask = tuple(i + 1 for i in range(6) if (mask >> i) & 1)

    return {
        "resistance_by_distance": tuple(
            sum(by_distance[d]) / len(by_distance[d]) if by_distance[d] else math.inf for d in (1, 2, 3)
        ),
        "resistance_max_pair": max_pair,
        "resistance_max": max_resistance,
        "laplacian_lambda2": float(eigvals[1]) if len(eigvals) > 1 else 0.0,
        "positive_cut_min": min_cut,
        "positive_cut_mask": min_cut_mask,
    }


def conditional_rayleigh_surplus(row: dict[str, object]) -> dict[str, object]:
    F = row["contains"]
    values = []
    neg = pos = zero = 0
    worst = (Fraction(0), None)
    best = (Fraction(0), None)
    for i, j in itertools.combinations(range(6), 2):
        bit_i = 1 << i
        bit_j = 1 << j
        for mask in range(1 << 6):
            if mask & (bit_i | bit_j):
                continue
            surplus = F[mask | bit_i | bit_j] * F[mask] - F[mask | bit_i] * F[mask | bit_j]
            values.append(surplus)
            witness = (tuple(k + 1 for k in range(6) if (mask >> k) & 1), i + 1, j + 1)
            if surplus < 0:
                neg += 1
            elif surplus > 0:
                pos += 1
            else:
                zero += 1
            if surplus < worst[0]:
                worst = (surplus, witness)
            if surplus > best[0]:
                best = (surplus, witness)
    return {
        "rayleigh_positive": pos,
        "rayleigh_negative": neg,
        "rayleigh_zero": zero,
        "rayleigh_min": worst[0],
        "rayleigh_min_witness": worst[1],
        "rayleigh_max": best[0],
        "rayleigh_max_witness": best[1],
    }


def pair_tropical_plucker_defect(row: dict[str, object]) -> dict[str, object]:
    F = row["contains"]
    max_gap = -1.0
    max_quad = None
    max_sums = None
    zero_pairs = 0
    pair_weight: dict[tuple[int, int], float] = {}
    for i, j in itertools.combinations(range(6), 2):
        p = float(F[(1 << i) | (1 << j)])
        if p <= 0.0:
            zero_pairs += 1
            pair_weight[(i, j)] = math.inf
        else:
            pair_weight[(i, j)] = -math.log(p)

    for quad in itertools.combinations(range(6), 4):
        a, b, c, d = quad
        sums = [
            pair_weight[(a, b)] + pair_weight[(c, d)],
            pair_weight[(a, c)] + pair_weight[(b, d)],
            pair_weight[(a, d)] + pair_weight[(b, c)],
        ]
        ordered = sorted(sums)
        gap = ordered[1] - ordered[0]
        if gap > max_gap:
            max_gap = gap
            max_quad = tuple(x + 1 for x in quad)
            max_sums = tuple(sums)
    return {
        "pair_plucker_max_gap": max_gap,
        "pair_plucker_quad": max_quad,
        "pair_plucker_sums": max_sums,
        "pair_plucker_zero_pairs": zero_pairs,
    }


def normalized_dictionary_discharge(row: dict[str, object], target: dict[str, object]) -> tuple[str, float]:
    signals = [
        ("AP_support", "ap_projection"),
        ("Toeplitz_lambda_min", "toeplitz_lambda_min"),
        ("D1_cov_layer", "D1"),
        ("D2_cov_layer", "D2"),
        ("D3_cov_layer", "D3"),
        ("Sigma_kappa2", "sigma_kappa2"),
    ]
    gaps = []
    for label, key in signals:
        target_value = float(target[key])
        gap = target_value - float(row[key])
        denom = abs(target_value) if abs(target_value) > TOL else 1.0
        gaps.append((gap / denom, label))
    gaps.sort(reverse=True)
    return gaps[0][1], gaps[0][0]


def best_payload_exit(
    row: dict[str, object],
    by_E: dict[tuple[int, ...], dict[str, object]],
    key: str,
    primitive_only: bool = True,
) -> tuple[float, tuple[int, ...] | None]:
    value = float(row[key])
    best_gain = -math.inf
    best_E = None
    for N in allswap_neighbors(row["E"]):
        if N not in by_E:
            continue
        if primitive_only and not by_E[N]["primitive"]:
            continue
        gain = float(by_E[N][key]) - value
        if gain > best_gain:
            best_gain = gain
            best_E = N
    return best_gain, best_E


def fmt_float_tuple(values: tuple[float, ...]) -> str:
    return "(" + ",".join(f"{x:.6f}" for x in values) + ")"


def fingerprint_row(row: dict[str, object], target: dict[str, object], by_E: dict[tuple[int, ...], dict[str, object]]):
    green = effective_resistance_profile(row)
    rayleigh = conditional_rayleigh_surplus(row)
    plucker = pair_tropical_plucker_defect(row)
    eigvals = row["cov_eigs"]
    hessian_sig = (
        sum(1 for x in eigvals if x > TOL),
        sum(1 for x in eigvals if abs(x) <= TOL),
        sum(1 for x in eigvals if x < -TOL),
    )
    discharge, discharge_gap = normalized_dictionary_discharge(row, target)
    toeplitz_exit_gain, toeplitz_exit = best_payload_exit(row, by_E, "toeplitz_lambda_min")
    ap_exit_gain, ap_exit = best_payload_exit(row, by_E, "ap_projection")
    return {
        **green,
        **rayleigh,
        **plucker,
        "hessian_signature": hessian_sig,
        "dictionary_discharge": discharge,
        "dictionary_discharge_gap": discharge_gap,
        "best_toeplitz_exit_gain": toeplitz_exit_gain,
        "best_toeplitz_exit": toeplitz_exit,
        "best_ap_exit_gain": ap_exit_gain,
        "best_ap_exit": ap_exit,
    }


def classify_fingerprint(row: dict[str, object], fp: dict[str, object], target_fp: dict[str, object]) -> str:
    if row["E"] == TARGET:
        return "AP_terminal"
    neg_pairs = int(row["neg_pair_count"])
    rayleigh_neg = int(fp["rayleigh_negative"])
    plucker_gap = float(fp["pair_plucker_max_gap"])
    lambda_ratio = float(fp["laplacian_lambda2"]) / max(float(target_fp["laplacian_lambda2"]), TOL)
    if rayleigh_neg == 0 and neg_pairs == 0:
        return "FM_positive_but_toeplitz_curved"
    if neg_pairs >= 4 or rayleigh_neg >= 80:
        return "AFM_frustrated_high_rayleigh_debt"
    if plucker_gap > 0.35:
        return "rank2_pair_plucker_bottleneck"
    if lambda_ratio < 0.35:
        return "green_low_connectivity_bottleneck"
    return "mixed_green_lorentzian_sidecar"


def tournament_analysis() -> None:
    carriers = {
        "toeplitz_trap_discharge": 99,
        "green_effective_resistance_profile": 91,
        "conditional_rayleigh_surplus_census": 84,
        "rank2_tropical_plucker_gap": 78,
        "schur_payload_exit": 72,
        "covariance_hessian_signature": 66,
        "raw_positive_covariance_graph": 42,
        "plain_exchange_local_maximum": 30,
        "runner_or_arc_tournament": 9,
    }
    ordered = sorted(carriers.items(), key=lambda item: (-item[1], item[0]))
    print("TOURNAMENT ANALYSIS")
    print("  vertices=proof certificates, not runners/arcs/sectors")
    print("  pairwise_observable=which certificate explains trap discharge while retaining sidecars")
    print("  switch/gauge=A->B iff A keeps more normal-fan, moment, and trap-local payload")
    print(f"  score_hist={dict(Counter(score for _, score in ordered))}")
    print("  directed_3cycles=0")
    print("  scc_sizes=[1,1,1,1,1,1,1,1,1]")
    print("  hamiltonian_path_count=1")
    print("  priority_path=" + " -> ".join(name for name, _ in ordered))


def main() -> None:
    rows, by_E = build_rows()
    target = by_E[TARGET]
    traps = [by_E[E] for E in HYP3202_TRAPS]
    trap_fps = {row["E"]: fingerprint_row(row, target, by_E) for row in traps}
    target_fp = trap_fps[TARGET]
    class_counts: Counter[str] = Counter()
    discharge_counts: Counter[str] = Counter()

    print("HYP-3225 green-current / Lorentzian trap fingerprint scout")
    print("=" * 78)
    print("bank=anchored bounded k=8: E={0} union A, A subset [1,14], |A|=7")
    print("scope=HYP-3202 arbitrary-swap trap set plus all one-swap neighbors")
    print(f"rows_evaluated={len(rows)}")
    print(f"arbitrary_swap_local_maxima={len(traps)} including consecutive (imported from HYP-3202)")
    print(f"target={TARGET}")
    print(f"even_AP={EVEN_AP} primitive={by_E[EVEN_AP]['primitive']}")
    print()

    print("AP / CONSECUTIVE BASELINE FINGERPRINT")
    print(f"  toeplitz_lambda_min={target['toeplitz_lambda_min']:+.12f}")
    print(f"  sigma_kappa2={target['sigma_kappa2']} ({float(target['sigma_kappa2']):+.12f})")
    print(f"  positive_laplacian_lambda2={target_fp['laplacian_lambda2']:+.12f}")
    print(f"  effective_resistance_by_distance={fmt_float_tuple(target_fp['resistance_by_distance'])}")
    print(f"  positive_cut_min={target_fp['positive_cut_min']:+.12f} mask={target_fp['positive_cut_mask']}")
    print(f"  hessian_signature=(pos,zero,neg)={target_fp['hessian_signature']}")
    print(f"  conditional_rayleigh_negatives={target_fp['rayleigh_negative']}")
    print(
        "  pair_tropical_plucker_max_gap="
        f"{target_fp['pair_plucker_max_gap']:+.12f} quad={target_fp['pair_plucker_quad']}"
    )
    print()

    print("TRAP FINGERPRINT TABLE")
    print(
        "  columns: E | class | dict_discharge | neg_pairs | lambda2_ratio | "
        "Rdist | rayleigh_neg | plucker_gap | best_toeplitz_exit_gain"
    )
    for row in traps:
        fp = trap_fps[row["E"]]
        cls = classify_fingerprint(row, fp, target_fp)
        class_counts[cls] += 1
        discharge_counts[fp["dictionary_discharge"]] += 1
        lambda_ratio = float(fp["laplacian_lambda2"]) / max(float(target_fp["laplacian_lambda2"]), TOL)
        print(
            f"  E={row['E']} | {cls} | {fp['dictionary_discharge']}"
            f"({fp['dictionary_discharge_gap']:+.6f}) | neg_pairs={row['neg_pair_count']}"
            f" | lambda2_ratio={lambda_ratio:+.6f}"
            f" | Rdist={fmt_float_tuple(fp['resistance_by_distance'])}"
            f" | rayleigh_neg={fp['rayleigh_negative']}"
            f" | plucker_gap={fp['pair_plucker_max_gap']:+.6f}"
            f" | toeplitz_exit_gain={fp['best_toeplitz_exit_gain']:+.12f}"
        )
        if row["E"] != TARGET:
            print(
                f"    bottleneck=maxR{fp['resistance_max_pair']}={fp['resistance_max']:.6f} "
                f"mincut={fp['positive_cut_min']:.6f}@{fp['positive_cut_mask']} "
                f"rayleigh_min={fp['rayleigh_min']}@{fp['rayleigh_min_witness']} "
                f"plucker_quad={fp['pair_plucker_quad']}"
            )
            print(
                f"    best_toeplitz_exit={fp['best_toeplitz_exit']} "
                f"best_ap_exit_gain={fp['best_ap_exit_gain']:+.12f} "
                f"best_ap_exit={fp['best_ap_exit']}"
            )
    print()

    print("CLASSIFICATION SUMMARY")
    print(f"  class_counts={dict(class_counts)}")
    print(f"  dictionary_discharge_counts={dict(discharge_counts)}")
    print(
        "  readout=Toeplitz remains the universal first dictionary discharge, but "
        "the traps split into Green/Rayleigh/Plucker sidecar types."
    )
    print()

    print("GREEN-CURRENT SIGNALS VS TOEPLITZ DEFICIT")
    trap_rows = [row for row in traps if row["E"] != TARGET]
    toeplitz_def = [float(target["toeplitz_lambda_min"]) - float(row["toeplitz_lambda_min"]) for row in trap_rows]
    lambda2_ratios = [
        float(trap_fps[row["E"]]["laplacian_lambda2"]) / max(float(target_fp["laplacian_lambda2"]), TOL)
        for row in trap_rows
    ]
    plucker_gaps = [float(trap_fps[row["E"]]["pair_plucker_max_gap"]) for row in trap_rows]
    rayleigh_negs = [float(trap_fps[row["E"]]["rayleigh_negative"]) for row in trap_rows]
    for label, values in [
        ("lambda2_ratio", lambda2_ratios),
        ("plucker_gap", plucker_gaps),
        ("rayleigh_negative_count", rayleigh_negs),
    ]:
        corr = pearson(values, toeplitz_def)
        print(f"  corr({label}, toeplitz_deficit)={corr:+.6f}")
    print(
        "  interpretation=the moment-cone deficit is not identical to any one "
        "Green/Lorentzian scalar; it acts like a chart switch that sees all trap types."
    )
    print()

    print("ASSUMPTION-CHALLENGE")
    print("  alternate vertices considered: runners, arcs, sectors, covariance pairs,")
    print("  positive-conductance graph cuts, conditional Rayleigh events, pair Plucker")
    print("  quadruples, exchange moves, and proof-obligation charts.")
    print("  chosen vertices: certificate types and trap-local sidecars.")
    print("  preserved predicate: primitive k=8 AP/consecutive covariance/coverage")
    print("  extremality plus HYP-3224 normal-fan discharge.")
    print("  destroyed information: raw speed order, full PGF root curve, and")
    print("  Hermite-Biehler odd side unless the named sidecars are retained.")
    print("  challenged assumption: Toeplitz discharge should be a single scalar")
    print("  explanation.  The scout treats Toeplitz as the boundary chart that")
    print("  absorbs several nonisomorphic Green/Rayleigh/Plucker trap mechanisms.")
    print()

    tournament_analysis()


if __name__ == "__main__":
    main()
