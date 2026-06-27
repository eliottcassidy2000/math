#!/usr/bin/env python3
"""HYP-3115 scout: Minkowski, circuit complexity, Ising, and de Moivre cues.

This extends the HYP-3108/HYP-3109 Lee-Yang root bank with four extra
controlled-forgetting sidecars:

* Minkowski: short affine-relation pressure as a finite proxy for relation
  lattice/successive-minima data.
* Circuit complexity: small threshold circuits over proof signals, measuring
  which sidecars are needed to isolate the observed extremal row without using
  p0 as an input.
* Ising model: the #real-root stratum becomes a spin field on the one-swap
  graph; root-collision edges are domain walls.
* de Moivre quintic: the derivative G'(z) is a quintic.  Its distance from the
  translated de Moivre normal form is used as an algebraic-solvability signal,
  not as a proof.

Tournament Analysis uses proof sidecars/signals as vertices, not runners.
"""

from __future__ import annotations

import itertools
import math
from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from statistics import median


SECTORS = 7
DOMAIN = tuple(range(1, 14))
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


def pair_relation_count(row: tuple[int, ...]) -> int:
    """Count low-height affine walls a+b=c+d among distinct unordered pairs."""
    counts: Counter[int] = Counter()
    for a, b in itertools.combinations(row, 2):
        counts[a + b] += 1
    return sum(n * (n - 1) // 2 for n in counts.values())


def relation_covolume_proxy(row: tuple[int, ...]) -> float:
    # For a primitive one-form v, the kernel lattice covolume is ||v||.  This
    # is only a proxy because the LRC relation sidecar has extra packet data.
    return math.sqrt(sum(v * v for v in row if v != 0))


def axis_ising(q: tuple[Fraction, ...]) -> tuple[int, int, int]:
    vals = [float(x) for x in q]
    med = median(vals)
    spins = [1 if v >= med else -1 for v in vals]
    walls = sum(1 for a, b in zip(spins, spins[1:]) if a != b)
    energy = -sum(a * b for a, b in zip(spins, spins[1:]))
    magnetization = sum(spins)
    return walls, energy, magnetization


def shift_poly_desc(coeffs: list[float], shift: float) -> list[float]:
    """Return coefficients of P(y+shift), descending order."""
    deg = len(coeffs) - 1
    low = list(reversed(coeffs))
    out = [0.0 for _ in range(deg + 1)]
    for power, coeff in enumerate(low):
        for k in range(power + 1):
            out[k] += coeff * math.comb(power, k) * (shift ** (power - k))
    return list(reversed(out))


def de_moivre_quintic_residual(q: tuple[Fraction, ...]) -> tuple[float, float, float]:
    """Distance of G'(z) from translated x^5 + a*x^3 + a^2/5*x + b form."""
    deriv_low = [float(i * q[i]) for i in range(1, SECTORS)]
    if abs(deriv_low[-1]) < 1e-15:
        return 999.0, 999.0, 999.0
    desc = [c / deriv_low[-1] for c in reversed(deriv_low)]
    c4 = desc[1]
    depressed = shift_poly_desc(desc, -c4 / 5.0)
    a = depressed[2]
    forbidden_y2 = depressed[3]
    linear_gap = depressed[4] - (a * a / 5.0)
    scale = math.sqrt(1.0 + sum(c * c for c in depressed[2:]))
    residual = math.sqrt(forbidden_y2 * forbidden_y2 + linear_gap * linear_gap) / scale
    return residual, forbidden_y2 / scale, linear_gap / scale


@dataclass(frozen=True)
class Metrics:
    row: tuple[int, ...]
    q: tuple[Fraction, ...]
    p0: float
    lyk8: float
    n_real: int
    nearest_modulus: float
    segment_clearance: float
    apex7_error: float
    pair_relations: int
    relation_covolume: float
    minkowski_pressure: float
    axis_walls: int
    axis_energy: int
    axis_magnetization: int
    dm_residual: float
    dm_y2: float
    dm_linear_gap: float


def root_sidecars(row: tuple[int, ...], q: tuple[Fraction, ...]) -> tuple[int, float, float, float]:
    roots = pgf_roots(q)
    n_real = sum(1 for z in roots if abs(z.imag) < 1e-7)
    nearest = min(abs(z) for z in roots)
    clearance = min(distance_to_segment(z) for z in roots)
    upper_angles = sorted(
        abs(math.degrees(math.atan2(z.imag, z.real)))
        for z in roots
        if z.imag >= -1e-7
    )
    upper_angles = [a for a in upper_angles if 1e-7 < a < 180.0 - 1e-7]
    targets = [360.0 / 7.0, 2.0 * 360.0 / 7.0, 3.0 * 360.0 / 7.0]
    if len(upper_angles) >= 3:
        apex = sum(abs(a - b) for a, b in zip(upper_angles[:3], targets)) / 3.0
    else:
        apex = 999.0
    return n_real, nearest, clearance, apex


def metrics_for(row: tuple[int, ...]) -> Metrics:
    q = missdist(row)
    n_real, nearest, clearance, apex = root_sidecars(row, q)
    pair_rel = pair_relation_count(row)
    covol = relation_covolume_proxy(row)
    walls, axis_energy, axis_mag = axis_ising(q)
    dm_res, dm_y2, dm_linear = de_moivre_quintic_residual(q)
    return Metrics(
        row=row,
        q=q,
        p0=float(q[0]),
        lyk8=sum(float(q[t]) * (3.0 ** t) for t in range(SECTORS)),
        n_real=n_real,
        nearest_modulus=nearest,
        segment_clearance=clearance,
        apex7_error=apex,
        pair_relations=pair_rel,
        relation_covolume=covol,
        minkowski_pressure=pair_rel / max(covol, 1e-12),
        axis_walls=walls,
        axis_energy=axis_energy,
        axis_magnetization=axis_mag,
        dm_residual=dm_res,
        dm_y2=dm_y2,
        dm_linear_gap=dm_linear,
    )


def bank_rows() -> list[tuple[int, ...]]:
    return [tuple([0] + list(rest)) for rest in itertools.combinations(DOMAIN, 7)]


def one_swap_edges(rows: list[tuple[int, ...]]) -> list[tuple[tuple[int, ...], tuple[int, ...]]]:
    row_set = set(rows)
    edges = []
    for row in rows:
        s = set(row)
        for rem in row:
            if rem == 0:
                continue
            for add in DOMAIN:
                if add in s:
                    continue
                nb = tuple(sorted((s - {rem}) | {add}))
                if nb in row_set and row < nb:
                    edges.append((row, nb))
    return edges


@dataclass(frozen=True)
class Predicate:
    name: str
    values: frozenset[tuple[int, ...]]


def threshold_predicates(metrics: list[Metrics]) -> list[Predicate]:
    rows = [m.row for m in metrics]

    def high(name: str, vals: dict[tuple[int, ...], float], thresholds: list[float]) -> list[Predicate]:
        out = []
        for th in thresholds:
            out.append(Predicate(f"{name}>={th:.5g}", frozenset(r for r in rows if vals[r] >= th)))
        return out

    def low(name: str, vals: dict[tuple[int, ...], float], thresholds: list[float]) -> list[Predicate]:
        out = []
        for th in thresholds:
            out.append(Predicate(f"{name}<={th:.5g}", frozenset(r for r in rows if vals[r] <= th)))
        return out

    by = {m.row: m for m in metrics}
    preds = [
        Predicate("real_roots==0", frozenset(m.row for m in metrics if m.n_real == 0)),
        Predicate("axis_walls<=1", frozenset(m.row for m in metrics if m.axis_walls <= 1)),
        Predicate("axis_walls<=2", frozenset(m.row for m in metrics if m.axis_walls <= 2)),
    ]
    clearance = {r: by[r].segment_clearance for r in rows}
    nearest = {r: by[r].nearest_modulus for r in rows}
    apex = {r: by[r].apex7_error for r in rows}
    pair_rel = {r: float(by[r].pair_relations) for r in rows}
    mink = {r: by[r].minkowski_pressure for r in rows}
    dm = {r: by[r].dm_residual for r in rows}

    preds += high("segment_clearance", clearance, [0.25, 0.5, 0.75, 0.9])
    preds += high("nearest_root", nearest, [0.75, 1.0, 1.25, 1.45])
    preds += low("apex7_error", apex, [5.0, 8.0, 12.0, 20.0])
    preds += high("pair_relations", pair_rel, [8.0, 12.0, 16.0, 20.0])
    preds += high("minkowski_pressure", mink, [0.25, 0.5, 0.75, 1.0])
    preds += low("de_moivre_residual", dm, [0.05, 0.1, 0.2, 0.4])
    return [p for p in preds if p.values]


def find_small_circuit(metrics: list[Metrics]) -> tuple[list[str], int, int, int]:
    rows = frozenset(m.row for m in metrics)
    max_p0 = max(m.p0 for m in metrics)
    target = frozenset(m.row for m in metrics if abs(m.p0 - max_p0) < 1e-12)
    preds = threshold_predicates(metrics)
    best_combo: tuple[Predicate, ...] | None = None
    best_false_pos = len(rows)
    best_false_neg = len(target)
    for size in range(1, 5):
        for combo in itertools.combinations(preds, size):
            selected = rows
            for pred in combo:
                selected = selected & pred.values
            fp = len(selected - target)
            fn = len(target - selected)
            if fp == 0 and fn == 0:
                return [p.name for p in combo], size, 0, 0
            if (fn, fp, size) < (best_false_neg, best_false_pos, len(best_combo or ()) or 99):
                best_combo = combo
                best_false_pos = fp
                best_false_neg = fn
    return [p.name for p in (best_combo or ())], len(best_combo or ()), best_false_pos, best_false_neg


def tournament_sccs(adj: list[list[bool]]) -> list[list[int]]:
    n = len(adj)
    index = 0
    stack: list[int] = []
    on_stack: set[int] = set()
    idx: dict[int, int] = {}
    low: dict[int, int] = {}
    comps: list[list[int]] = []

    def visit(v: int) -> None:
        nonlocal index
        idx[v] = index
        low[v] = index
        index += 1
        stack.append(v)
        on_stack.add(v)
        for w in range(n):
            if not adj[v][w]:
                continue
            if w not in idx:
                visit(w)
                low[v] = min(low[v], low[w])
            elif w in on_stack:
                low[v] = min(low[v], idx[w])
        if low[v] == idx[v]:
            comp = []
            while True:
                w = stack.pop()
                on_stack.remove(w)
                comp.append(w)
                if w == v:
                    break
            comps.append(comp)

    for v in range(n):
        if v not in idx:
            visit(v)
    return comps


def hamiltonian_paths_tournament(adj: list[list[bool]]) -> int:
    n = len(adj)
    dp: dict[tuple[int, int], int] = {(1 << i, i): 1 for i in range(n)}
    for size in range(1, n):
        new: dict[tuple[int, int], int] = {}
        for (mask, last), count in dp.items():
            if mask.bit_count() != size:
                new[(mask, last)] = new.get((mask, last), 0) + count
                continue
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                if adj[last][nxt]:
                    key = (mask | (1 << nxt), nxt)
                    new[key] = new.get(key, 0) + count
            new[(mask, last)] = new.get((mask, last), 0) + count
        dp = new
    full = (1 << n) - 1
    return sum(count for (mask, _), count in dp.items() if mask == full)


def signal_tournament() -> None:
    signals = [
        "root_zero_locus",
        "ising_domain_wall_tension",
        "minkowski_relation_lattice",
        "circuit_depth_sidecar",
        "de_moivre_quintic_residual",
        "bravais_q_lattice_rank",
        "observer_gluing_certificate",
        "single_scalar_p0",
    ]
    axes = [
        "preserves_lrc_predicate",
        "retains_finite_witness",
        "detects_phase_wall",
        "certifiable",
        "controlled_forgetting",
        "computable",
        "algebraic_structure",
        "lower_bound_signal",
    ]
    score = {
        "root_zero_locus":              [2, 1, 2, 1, 2, 2, 1, 1],
        "ising_domain_wall_tension":    [1, 1, 2, 1, 1, 2, 0, 1],
        "minkowski_relation_lattice":   [2, 2, 1, 2, 2, 1, 1, 1],
        "circuit_depth_sidecar":        [1, 2, 0, 2, 2, 2, 0, 2],
        "de_moivre_quintic_residual":   [1, 1, 1, 1, 1, 2, 2, 1],
        "bravais_q_lattice_rank":       [1, 2, 0, 1, 2, 2, 1, 1],
        "observer_gluing_certificate":  [2, 2, 1, 2, 2, 1, 0, 1],
        "single_scalar_p0":             [1, 0, 0, 0, 0, 2, 0, 0],
    }
    tie = {name: i for i, name in enumerate(signals)}
    n = len(signals)
    adj = [[False] * n for _ in range(n)]
    for i, a in enumerate(signals):
        for j, b in enumerate(signals):
            if i >= j:
                continue
            av = sum(1 for k in range(len(axes)) if score[a][k] > score[b][k])
            bv = sum(1 for k in range(len(axes)) if score[b][k] > score[a][k])
            if av > bv or (av == bv and tie[a] < tie[b]):
                adj[i][j] = True
            else:
                adj[j][i] = True
    out_scores = {signals[i]: sum(adj[i]) for i in range(n)}
    tri = 0
    for a, b, c in itertools.combinations(range(n), 3):
        if (adj[a][b] and adj[b][c] and adj[c][a]) or (adj[a][c] and adj[c][b] and adj[b][a]):
            tri += 1
    sccs = [[signals[i] for i in comp] for comp in tournament_sccs(adj)]
    print("Tournament Analysis:")
    print(f"  vertices={signals}")
    print(f"  axes={axes}")
    print(f"  score_hist={dict(sorted(Counter(out_scores.values()).items()))}")
    print(f"  scores={out_scores}")
    print(f"  directed_3cycles={tri}")
    print(f"  SCCs={sccs}")
    print(f"  Hamiltonian_path_count={hamiltonian_paths_tournament(adj)}")
    print("  priority_path=" + " > ".join(sorted(signals, key=lambda x: -out_scores[x])))
    print("  assumption_challenge=vertices considered included runners, gaps,")
    print("    relation-lattice vectors, Boolean gates, spin domains, stationary")
    print("    quintic coefficients, PGF roots, Bravais cells, and proof obligations;")
    print("    chosen vertices are proof sidecars, preserving predicate-retention")
    print("    order and destroying raw row identity.")


def main() -> None:
    rows = bank_rows()
    print("HYP-3115 Minkowski / circuit / Ising / de Moivre bridge scout")
    print(f"anchored_bank_rows={len(rows)}")
    metrics = [metrics_for(row) for row in rows]
    by_row = {m.row: m for m in metrics}
    top = max(metrics, key=lambda m: m.p0)
    consec = by_row[tuple(range(8))]
    competitor = by_row[(0, 2, 4, 6, 7, 8, 10, 12)]
    print("\nNamed rows:")
    for label, m in [("consec_8", consec), ("nearest_root_competitor", competitor), ("top_p0", top)]:
        print(
            f"  {label}: row={m.row}, p0={m.p0:.5f}, n_real={m.n_real}, "
            f"clearance={m.segment_clearance:.5f}, nearest={m.nearest_modulus:.5f}, "
            f"apex7={m.apex7_error:.3f}, pair_rel={m.pair_relations}, "
            f"minkowski_pressure={m.minkowski_pressure:.5f}, axis_walls={m.axis_walls}, "
            f"dm_residual={m.dm_residual:.5f}"
        )
    print("\nMinkowski-style relation sidecar:")
    xs = [m.p0 for m in metrics]
    print(f"  corr(p0, pair_relation_count)={pearson(xs, [m.pair_relations for m in metrics]):+.3f}")
    print(f"  corr(p0, minkowski_pressure)={pearson(xs, [m.minkowski_pressure for m in metrics]):+.3f}")
    rel_leaders = sorted(metrics, key=lambda m: (-m.pair_relations, -m.p0))[:5]
    for m in rel_leaders:
        print(f"  relation_leader row={m.row} pair_rel={m.pair_relations} p0={m.p0:.5f}")
    print("\nIsing root-stratum graph:")
    edges = one_swap_edges(rows)
    spins = {m.row: (1 if m.n_real == 0 else -1) for m in metrics}
    transitions = Counter()
    energy = 0
    for a, b in edges:
        sa, sb = spins[a], spins[b]
        transitions[(sa, sb)] += 1
        energy += -(sa * sb)
    walls = sum(n for (a, b), n in transitions.items() if a != b)
    same = len(edges) - walls
    print(f"  one_swap_edges={len(edges)}")
    print(f"  spin_counts={dict(Counter(spins.values()))}  (+1 means #real=0)")
    print(f"  same_edges={same}, domain_wall_edges={walls}, wall_fraction={walls/len(edges):.5f}")
    print(f"  Ising_energy_J1={energy}, magnetization={sum(spins.values())}")
    print(f"  transition_edges={dict(sorted(transitions.items()))}")
    print("\nCircuit-complexity sidecar:")
    combo, size, fp, fn = find_small_circuit(metrics)
    binary_depth = math.ceil(math.log2(max(size, 1))) if size else 0
    print(f"  target=max_p0_rows, max_p0={top.p0:.5f}, target_row={top.row}")
    print(f"  minimal_observed_non_p0_conjunction={combo}")
    print(f"  input_literals={size}, false_positives={fp}, false_negatives={fn}")
    print(f"  binary_AND_depth={binary_depth}, unbounded_fanin_AND_depth={1 if size else 0}")
    print("  note=This is an observed finite-bank circuit, not a lower bound.")
    print("\nde Moivre derivative-quintic sidecar:")
    print(f"  corr(p0, -de_moivre_residual)={pearson(xs, [-m.dm_residual for m in metrics]):+.3f}")
    dm_leaders = sorted(metrics, key=lambda m: (m.dm_residual, -m.p0))[:5]
    for m in dm_leaders:
        print(f"  dm_leader row={m.row} residual={m.dm_residual:.5f} p0={m.p0:.5f} n_real={m.n_real}")
    print("  interpretation=G'(z) is quintic; residual measures distance to")
    print("    translated de Moivre normal form x^5+a*x^3+(a^2/5)*x+b.")
    print("\nCross-signal correlations:")
    print(f"  corr(segment_clearance, -Ising_axis_walls)={pearson([m.segment_clearance for m in metrics], [-m.axis_walls for m in metrics]):+.3f}")
    print(f"  corr(pair_relations, -de_moivre_residual)={pearson([m.pair_relations for m in metrics], [-m.dm_residual for m in metrics]):+.3f}")
    print(f"  corr(n_real, de_moivre_residual)={pearson([m.n_real for m in metrics], [m.dm_residual for m in metrics]):+.3f}")
    print()
    signal_tournament()


if __name__ == "__main__":
    main()
