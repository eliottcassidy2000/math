#!/usr/bin/env python3
"""Lee-Yang root-curve and ear-map scout for the LRC miss-count PGF.

This extends HYP-3103/HYP-3104/HYP-3106 from a real-root count signal to the
whole root curve of

    G_N(z) = sum_t q_t z^t,  q_t = P(# empty inner sectors = t).

New signals measured here:
  * root clearance from the forbidden real segment [-1, 0];
  * apex-7 root-angle deviation for the upper half-plane roots;
  * an even quartic potential fit V(s) ~= c + b*s^2 + lambda*s^4 to the
    symmetrized miss-count law, inspired by exp(-lambda*S^4 - b*S^2)dS;
  * the one-swap graph of the #real=0 Lee-Yang stratum, as an "ear-map"
    candidate rather than a proof.

Tournament Analysis uses proof signals as vertices, not runners.
"""

from __future__ import annotations

import itertools
import math
from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction


SECTORS = 7
BANK_N = 14
K = 8
EPS = 1e-8


def sector_of(x: Fraction) -> int:
    return int((x % 1) * SECTORS)


def missdist(row: tuple[int, ...]) -> tuple[Fraction, ...]:
    row = tuple(sorted(set(row)))
    breaks = {Fraction(0), Fraction(1)}
    for e in row:
        if e == 0:
            continue
        for m in range(SECTORS * e + 1):
            breaks.add(Fraction(m, SECTORS * e))
    cuts = sorted(breaks)
    q = [Fraction(0) for _ in range(SECTORS)]
    for left, right in zip(cuts, cuts[1:]):
        if right <= left:
            continue
        mid = (left + right) / 2
        occupied = {sector_of(e * mid) for e in row}
        empty = SECTORS - len(occupied)
        if 0 <= empty <= SECTORS - 1:
            q[empty] += right - left
    return tuple(q)


def poly_eval(coeffs: list[float], z: complex) -> complex:
    out = 0j
    for c in coeffs:
        out = out * z + c
    return out


def durand_kerner(coeffs: list[float]) -> list[complex]:
    while len(coeffs) > 1 and abs(coeffs[0]) < 1e-15:
        coeffs = coeffs[1:]
    degree = len(coeffs) - 1
    if degree <= 0:
        return []
    lead = coeffs[0]
    coeffs = [c / lead for c in coeffs]
    radius = 1.0 + max(abs(c) for c in coeffs[1:])
    roots = [
        radius * complex(math.cos(2.0 * math.pi * (i + 0.37) / degree),
                         math.sin(2.0 * math.pi * (i + 0.37) / degree))
        for i in range(degree)
    ]
    for _ in range(300):
        max_delta = 0.0
        new_roots = roots[:]
        for i, root in enumerate(roots):
            denom = 1.0 + 0j
            for j, other in enumerate(roots):
                if i != j:
                    denom *= root - other
            if abs(denom) < 1e-18:
                denom += complex(1e-12, -1e-12)
            delta = poly_eval(coeffs, root) / denom
            new_roots[i] = root - delta
            max_delta = max(max_delta, abs(delta))
        roots = new_roots
        if max_delta < 1e-12:
            break
    return roots


def pgf_roots(q: tuple[Fraction, ...]) -> list[complex]:
    coeffs = [float(q[t]) for t in range(SECTORS - 1, -1, -1)]
    return durand_kerner(coeffs)


def eval_pgf(q: tuple[Fraction, ...], z: float) -> float:
    return sum(float(q[t]) * (z**t) for t in range(SECTORS))


def distance_to_segment(z: complex, left: float = -1.0, right: float = 0.0) -> float:
    x = min(max(z.real, left), right)
    return abs(z - x)


def pearson(xs: list[float], ys: list[float]) -> float:
    mx = sum(xs) / len(xs)
    my = sum(ys) / len(ys)
    vx = sum((x - mx) ** 2 for x in xs)
    vy = sum((y - my) ** 2 for y in ys)
    if vx <= 0 or vy <= 0:
        return 0.0
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / math.sqrt(vx * vy)


@dataclass(frozen=True)
class QuarticFit:
    b: float
    lam: float
    r2: float
    well_radius: float
    phase: str


def solve_3x3(a: list[list[float]], b: list[float]) -> list[float]:
    mat = [row[:] + [rhs] for row, rhs in zip(a, b)]
    for col in range(3):
        pivot = max(range(col, 3), key=lambda r: abs(mat[r][col]))
        mat[col], mat[pivot] = mat[pivot], mat[col]
        if abs(mat[col][col]) < 1e-15:
            return [0.0, 0.0, 0.0]
        scale = mat[col][col]
        for k in range(col, 4):
            mat[col][k] /= scale
        for r in range(3):
            if r == col:
                continue
            factor = mat[r][col]
            for k in range(col, 4):
                mat[r][k] -= factor * mat[col][k]
    return [mat[i][3] for i in range(3)]


def quartic_fit(q: tuple[Fraction, ...]) -> QuarticFit:
    # Symmetrize around the middle sector count 3.  The potential fit is not a
    # probability model; it is a two-well/single-well diagnostic.
    sym = [
        float(q[3]),
        0.5 * (float(q[2]) + float(q[4])),
        0.5 * (float(q[1]) + float(q[5])),
        0.5 * (float(q[0]) + float(q[6])),
    ]
    rows = [[1.0, s * s, s**4] for s in [0.0, 1.0, 2.0, 3.0]]
    y = [-math.log(max(v, 1e-15)) for v in sym]
    ata = [[sum(row[i] * row[j] for row in rows) for j in range(3)] for i in range(3)]
    aty = [sum(row[i] * yy for row, yy in zip(rows, y)) for i in range(3)]
    coef = solve_3x3(ata, aty)
    pred = [sum(row[i] * coef[i] for i in range(3)) for row in rows]
    ss_res = sum((yy - pp) ** 2 for yy, pp in zip(y, pred))
    mean_y = sum(y) / len(y)
    ss_tot = sum((yy - mean_y) ** 2 for yy in y)
    r2 = 1.0 if ss_tot <= 1e-15 else 1.0 - ss_res / ss_tot
    b = float(coef[1])
    lam = float(coef[2])
    if lam > 0 and b < 0:
        phase = "double-well"
        well_radius = math.sqrt(max(0.0, -b / (2.0 * lam)))
    elif lam > 0:
        phase = "single-well"
        well_radius = 0.0
    else:
        phase = "unstable-fit"
        well_radius = 0.0
    return QuarticFit(b=b, lam=lam, r2=r2, well_radius=well_radius, phase=phase)


@dataclass(frozen=True)
class Metrics:
    row: tuple[int, ...]
    q: tuple[Fraction, ...]
    roots: tuple[complex, ...]
    n_real: int
    nearest_modulus: float
    segment_clearance: float
    angle_deviation_7: float
    radius_spread: float
    log_concavity_defect: float
    quartic: QuarticFit
    p0: float
    extreme_mass: float
    lyk8: float


def root_metrics(row: tuple[int, ...]) -> Metrics:
    q = missdist(row)
    roots = tuple(complex(z) for z in pgf_roots(q))
    n_real = sum(1 for z in roots if abs(z.imag) < 1e-7)
    nearest_modulus = min(abs(z) for z in roots)
    segment_clearance = min(distance_to_segment(z) for z in roots)
    radii = [abs(z) for z in roots]
    upper_angles = sorted(
        abs(math.degrees(math.atan2(z.imag, z.real)))
        for z in roots
        if z.imag >= -1e-7
    )
    upper_angles = [a for a in upper_angles if 1e-7 < a < 180.0 - 1e-7]
    targets = [360.0 / 7.0, 2.0 * 360.0 / 7.0, 3.0 * 360.0 / 7.0]
    if len(upper_angles) >= 3:
        angle_deviation_7 = sum(abs(a - b) for a, b in zip(upper_angles[:3], targets)) / 3.0
    else:
        angle_deviation_7 = 999.0
    lc = min(float(q[t]) ** 2 - float(q[t - 1]) * float(q[t + 1]) for t in range(1, 6))
    return Metrics(
        row=row,
        q=q,
        roots=roots,
        n_real=n_real,
        nearest_modulus=nearest_modulus,
        segment_clearance=segment_clearance,
        angle_deviation_7=angle_deviation_7,
        radius_spread=max(radii) - min(radii),
        log_concavity_defect=lc,
        quartic=quartic_fit(q),
        p0=float(q[0]),
        extreme_mass=float(q[0] + q[6]),
        lyk8=float(10 * q[0] + q[3] + 10 * q[6]),
    )


def row_label(row: tuple[int, ...]) -> str:
    return "(" + ",".join(str(x) for x in row) + ")"


def fmt(x: float) -> str:
    return f"{x:.5f}"


def one_swap_neighbors(row: tuple[int, ...], universe: range = range(BANK_N)) -> list[tuple[int, ...]]:
    s = set(row)
    out = []
    for old in row:
        base = s - {old}
        for new in universe:
            if new in s:
                continue
            out.append(tuple(sorted(base | {new})))
    return out


def connected_components(vertices: set[tuple[int, ...]], adj: dict[tuple[int, ...], set[tuple[int, ...]]]):
    seen: set[tuple[int, ...]] = set()
    comps = []
    for v in vertices:
        if v in seen:
            continue
        stack = [v]
        seen.add(v)
        comp = []
        while stack:
            x = stack.pop()
            comp.append(x)
            for y in adj.get(x, ()):
                if y not in seen:
                    seen.add(y)
                    stack.append(y)
        comps.append(comp)
    return comps


def articulation_points(vertices: set[tuple[int, ...]], adj: dict[tuple[int, ...], set[tuple[int, ...]]]):
    timer = 0
    disc: dict[tuple[int, ...], int] = {}
    low: dict[tuple[int, ...], int] = {}
    parent: dict[tuple[int, ...], tuple[int, ...] | None] = {}
    arts: set[tuple[int, ...]] = set()

    def dfs(u: tuple[int, ...]) -> None:
        nonlocal timer
        timer += 1
        disc[u] = low[u] = timer
        child_count = 0
        for v in adj.get(u, ()):
            if v not in disc:
                parent[v] = u
                child_count += 1
                dfs(v)
                low[u] = min(low[u], low[v])
                if parent.get(u) is None and child_count > 1:
                    arts.add(u)
                if parent.get(u) is not None and low[v] >= disc[u]:
                    arts.add(u)
            elif v != parent.get(u):
                low[u] = min(low[u], disc[v])

    for v in vertices:
        if v not in disc:
            parent[v] = None
            dfs(v)
    return arts


def hamiltonian_paths_tournament(adj: list[list[bool]]) -> int:
    n = len(adj)
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        for last in range(n):
            count = dp[mask][last]
            if not count:
                continue
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                if adj[last][nxt]:
                    dp[mask | (1 << nxt)][nxt] += count
    return sum(dp[-1])


def tournament_sccs(adj: list[list[bool]]) -> list[list[int]]:
    n = len(adj)
    order: list[int] = []
    seen = [False] * n

    def dfs(v: int) -> None:
        seen[v] = True
        for w in range(n):
            if adj[v][w] and not seen[w]:
                dfs(w)
        order.append(v)

    def rdfs(v: int, comp: list[int]) -> None:
        seen[v] = True
        comp.append(v)
        for w in range(n):
            if adj[w][v] and not seen[w]:
                rdfs(w, comp)

    for v in range(n):
        if not seen[v]:
            dfs(v)
    seen = [False] * n
    comps = []
    for v in reversed(order):
        if not seen[v]:
            comp: list[int] = []
            rdfs(v, comp)
            comps.append(comp)
    return comps


def print_named_rows() -> None:
    print("=" * 96)
    print("NAMED ROWS: whole PGF roots, Lee-Yang clearance, and quartic potential")
    print("=" * 96)
    named = {
        "consec_8": tuple(range(8)),
        "even_AP_8": tuple(2 * i for i in range(8)),
        "top_cluster": (0, 7, 8, 9, 10, 11, 12, 13),
        "break_spread": (0, 1, 5, 7, 8, 9, 11, 13),
        "covering_probe": (0, 1, 2, 3, 5, 8, 13, 14),
        "random_spread": (0, 1, 3, 7, 12, 20, 33, 54),
    }
    for name, row in named.items():
        m = root_metrics(row)
        angles = sorted(
            abs(math.degrees(math.atan2(z.imag, z.real)))
            for z in m.roots
            if z.imag > 1e-7
        )
        print(f"\n{name:15s} row={row_label(row)}")
        print(
            f"  p0={fmt(m.p0)}  extreme={fmt(m.extreme_mass)}  L_yK8={fmt(m.lyk8)} "
            f"nreal={m.n_real}  min|z|={fmt(m.nearest_modulus)}  "
            f"dist(roots,[-1,0])={fmt(m.segment_clearance)}"
        )
        print(
            "  upper root angles="
            + "[" + ", ".join(f"{a:.1f}" for a in angles) + "]"
            + f"  mean apex7 angle error={m.angle_deviation_7:.3f} deg"
        )
        print(
            f"  quartic: phase={m.quartic.phase:12s} "
            f"b={m.quartic.b:+.5f} lambda={m.quartic.lam:+.5f} "
            f"well_radius={m.quartic.well_radius:.3f} r2={m.quartic.r2:.4f}"
        )


def print_consecutive_ladder() -> None:
    print("\n" + "=" * 96)
    print("CONSECUTIVE LADDER k=8..13: apex-7 root arc and phi4 phase")
    print("=" * 96)
    print("k   p0       nreal  min|z|   seg_clear  apex7_err  b        lambda    phase")
    for k in range(8, 14):
        row = tuple(range(k))
        m = root_metrics(row)
        print(
            f"{k:<2d}  {m.p0:.5f}  {m.n_real:>2d}     {m.nearest_modulus:.4f}  "
            f"{m.segment_clearance:.4f}    {m.angle_deviation_7:8.3f}  "
            f"{m.quartic.b:+.4f}  {m.quartic.lam:+.4f}  {m.quartic.phase}"
        )


def analyze_bank() -> tuple[list[Metrics], dict[tuple[int, ...], Metrics]]:
    bank = [tuple((0,) + c) for c in itertools.combinations(range(1, BANK_N), K - 1)]
    metrics = [root_metrics(row) for row in bank]
    by_row = {m.row: m for m in metrics}
    print("\n" + "=" * 96)
    print(f"EXHAUSTIVE ANCHORED RESIDUE BANK: rows 0 union A, |A|={K - 1}, A subset {{1..13}}")
    print(f"count=C({BANK_N - 1},{K - 1})={len(bank)}")
    print("=" * 96)
    groups: dict[int, list[Metrics]] = defaultdict(list)
    for m in metrics:
        groups[m.n_real].append(m)
    print("#real -> count, max p0, mean p0, max L_yK8, phase counts")
    for nr in sorted(groups):
        rows = groups[nr]
        phases = Counter(m.quartic.phase for m in rows)
        print(
            f"  {nr:>2d}: n={len(rows):4d}  max_p0={max(m.p0 for m in rows):.5f} "
            f"mean_p0={sum(m.p0 for m in rows)/len(rows):.5f}  "
            f"max_L_yK8={max(m.lyk8 for m in rows):.5f}  phases={dict(phases)}"
        )
    top_p0 = sorted(metrics, key=lambda m: (-m.p0, m.n_real, -m.nearest_modulus))[:8]
    top_clear = sorted(metrics, key=lambda m: (-m.segment_clearance, m.n_real, -m.p0))[:8]
    top_radius = sorted(metrics, key=lambda m: (-m.nearest_modulus, m.n_real, -m.p0))[:8]
    print("\nTop p0 rows:")
    for m in top_p0:
        print(
            f"  {row_label(m.row):34s} p0={m.p0:.5f} nreal={m.n_real} "
            f"min|z|={m.nearest_modulus:.4f} phase={m.quartic.phase}"
        )
    print("\nTop Lee-Yang segment-clearance rows:")
    for m in top_clear:
        print(
            f"  {row_label(m.row):34s} clear={m.segment_clearance:.4f} "
            f"p0={m.p0:.5f} nreal={m.n_real} phase={m.quartic.phase}"
        )
    print("\nTop nearest-root-radius rows:")
    for m in top_radius:
        print(
            f"  {row_label(m.row):34s} min|z|={m.nearest_modulus:.4f} "
            f"p0={m.p0:.5f} nreal={m.n_real} phase={m.quartic.phase}"
        )
    print("\nCorrelations over the residue bank:")
    print(f"  corr(p0, -#real)               = {pearson([m.p0 for m in metrics], [-m.n_real for m in metrics]):+.3f}")
    print(f"  corr(p0, min root modulus)     = {pearson([m.p0 for m in metrics], [m.nearest_modulus for m in metrics]):+.3f}")
    print(f"  corr(p0, segment clearance)    = {pearson([m.p0 for m in metrics], [m.segment_clearance for m in metrics]):+.3f}")
    print(f"  corr(p0, -apex7 angle error)   = {pearson([m.p0 for m in metrics], [-m.angle_deviation_7 for m in metrics]):+.3f}")
    print(f"  corr(p0, -quartic b)           = {pearson([m.p0 for m in metrics], [-m.quartic.b for m in metrics]):+.3f}")
    return metrics, by_row


def print_ear_map(metrics: list[Metrics], by_row: dict[tuple[int, ...], Metrics]) -> None:
    print("\n" + "=" * 96)
    print("ROOT-STRATUM ONE-SWAP GRAPH: ear-map signals")
    print("=" * 96)
    rows = {m.row for m in metrics}
    zero_rows = {m.row for m in metrics if m.n_real == 0}
    adj_zero: dict[tuple[int, ...], set[tuple[int, ...]]] = {r: set() for r in zero_rows}
    transition = Counter()
    total_edges = 0
    zero_boundary = Counter()
    for row in rows:
        for nb in one_swap_neighbors(row):
            if nb not in rows or row >= nb:
                continue
            total_edges += 1
            a = by_row[row].n_real
            b = by_row[nb].n_real
            transition[tuple(sorted((a, b)))] += 1
            if row in zero_rows and nb in zero_rows:
                adj_zero[row].add(nb)
                adj_zero[nb].add(row)
            elif row in zero_rows or nb in zero_rows:
                other = nb if row in zero_rows else row
                zero_boundary[by_row[other].n_real] += 1
    comps = connected_components(zero_rows, adj_zero) if zero_rows else []
    comp_sizes = sorted((len(c) for c in comps), reverse=True)
    edge_zero = sum(len(v) for v in adj_zero.values()) // 2
    cycle_rank = edge_zero - len(zero_rows) + len(comps)
    arts = articulation_points(zero_rows, adj_zero) if zero_rows else set()
    consec = tuple(range(8))
    consec_comp = next((c for c in comps if consec in c), [])
    print(f"Full Johnson one-swap graph edges: {total_edges}")
    print("Root-count transition edges:")
    for key in sorted(transition):
        print(f"  #real {key[0]} <-> {key[1]}: {transition[key]}")
    print(
        f"\n#real=0 stratum: vertices={len(zero_rows)} edges={edge_zero} "
        f"components={len(comps)} largest={comp_sizes[:5]} cycle_rank={cycle_rank}"
    )
    print(f"  articulation_points={len(arts)}")
    print(f"  consec component size={len(consec_comp)}")
    print(f"  boundary edges from #real=0 into other strata={dict(sorted(zero_boundary.items()))}")
    if consec_comp:
        best_in_comp = sorted((by_row[r] for r in consec_comp), key=lambda m: (-m.p0, m.row))[:5]
        print("  leading rows inside consec #real=0 component:")
        for m in best_in_comp:
            print(
                f"    {row_label(m.row):34s} p0={m.p0:.5f} "
                f"min|z|={m.nearest_modulus:.4f} clear={m.segment_clearance:.4f}"
            )
    print("\nEar reading:")
    print("  A bidirected connected zero-real component is strongly connected, hence has")
    print("  a directed ear decomposition.  The useful next test is not this tautology,")
    print("  but whether a small seed around consec plus root-collision ears generates")
    print("  the whole zero-real component before crossing into the #real=2 wall.")


def print_signal_tournament() -> None:
    print("\n" + "=" * 96)
    print("TOURNAMENT ANALYSIS ON LEE-YANG SIGNALS")
    print("=" * 96)
    signals = [
        "whole_pgf_curve",
        "root_zero_locus",
        "quartic_phi4_fit",
        "zero_stratum_ear_graph",
        "root_collision_discriminant",
        "nearest_zero_clearance",
        "real_root_count",
        "factorial_moments",
        "single_value_p0",
    ]
    axes = [
        "predicate",
        "whole_curve",
        "root_geometry",
        "phase_transition",
        "local_graph",
        "controlled_forgetting",
        "computable",
        "hypothesis_yield",
    ]
    score = {
        "whole_pgf_curve":              [2, 2, 1, 1, 0, 2, 2, 2],
        "root_zero_locus":              [2, 2, 2, 2, 0, 2, 2, 2],
        "quartic_phi4_fit":             [1, 1, 1, 2, 0, 1, 2, 2],
        "zero_stratum_ear_graph":       [1, 1, 1, 2, 2, 2, 1, 2],
        "root_collision_discriminant":  [1, 1, 2, 2, 1, 1, 1, 2],
        "nearest_zero_clearance":       [1, 1, 1, 1, 0, 1, 2, 1],
        "real_root_count":              [1, 0, 1, 1, 0, 1, 2, 1],
        "factorial_moments":            [1, 0, 0, 0, 0, 1, 2, 0],
        "single_value_p0":              [1, 0, 0, 0, 0, 0, 2, 0],
    }
    tie_path = {name: i for i, name in enumerate(signals)}
    n = len(signals)
    adj = [[False] * n for _ in range(n)]
    votes: dict[tuple[int, int], int] = {}
    for i, a in enumerate(signals):
        for j, b in enumerate(signals):
            if i == j:
                continue
            if i > j:
                continue
            av = sum(1 for ax in range(len(axes)) if score[a][ax] > score[b][ax])
            bv = sum(1 for ax in range(len(axes)) if score[b][ax] > score[a][ax])
            if av > bv or (av == bv and tie_path[a] < tie_path[b]):
                adj[i][j] = True
                votes[(i, j)] = av - bv
            else:
                adj[j][i] = True
                votes[(j, i)] = bv - av
    out_scores = {signals[i]: sum(adj[i]) for i in range(n)}
    tri = 0
    for a, b, c in itertools.combinations(range(n), 3):
        e = [adj[a][b], adj[b][c], adj[c][a]]
        f = [adj[a][c], adj[c][b], adj[b][a]]
        if all(e) or all(f):
            tri += 1
    sccs = [[signals[i] for i in comp] for comp in tournament_sccs(adj)]
    print(f"vertices={signals}")
    print(f"axes={axes}")
    print(f"score_hist={dict(sorted(Counter(out_scores.values()).items()))}")
    print(f"scores={out_scores}")
    print(f"directed_3cycles={tri}")
    print(f"SCCs={sccs}")
    print(f"Hamiltonian_path_count={hamiltonian_paths_tournament(adj)}")
    path = sorted(signals, key=lambda name: -out_scores[name])
    print("one priority path=" + " > ".join(path))
    print("\nAssumption challenge:")
    print("  Vertices considered: runners, gaps, sector boundaries, wall crossings,")
    print("  residues, cover arcs, Fourier modes, PGF roots, quartic potential wells,")
    print("  root-collision events, ear graph components, and proof obligations.")
    print("  Chosen vertices: signals/proof obligations.  This preserves the question")
    print("  'which coordinate must a proof retain next?' and destroys row identity.")


def main() -> None:
    print_named_rows()
    print_consecutive_ladder()
    metrics, by_row = analyze_bank()
    print_ear_map(metrics, by_row)
    print_signal_tournament()


if __name__ == "__main__":
    main()
