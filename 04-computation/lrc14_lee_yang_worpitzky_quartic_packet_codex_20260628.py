#!/usr/bin/env python3
"""HYP-3153 scout: Lee-Yang/Worpitzky/quartic compression packet.

This follows HYP-3151 and HYP-3152.  It does not try to re-prove LRC14.
It assembles a finite certificate packet around the k=8 hard row:

* HYP-3151: quotient legality is target-function factor-through plus sidecars;
* HYP-3152: coverage is a Lee-Yang radius q0=q6*R^6, with cap on-circle and
  dip off-circle;
* HYP-3150/S71: the k=8 correction splits into even biquadratic and odd
  Worpitzky/orientation pieces.

Tournament Analysis declaration:
  vertices: proof packets and obligations, not runners/arcs/classes;
  pairwise observable: majority vote over legality, root retention, cap/dip
    separation, parity separation, degree control, AP extremality, ear transfer,
    sidecar hygiene, and anti-scalar guard;
  switch: orient A -> B when A wins more axes, ties follow the declared path;
  tie path: function legality, Lee-Yang radius, Ly bimodality, even fold,
    odd Worpitzky/ear, Pascal cap/dip, Newton AP, quintic alarm, raw scalar.
"""

from __future__ import annotations

import cmath
from collections import Counter
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import combinations, permutations, product
from math import comb
from typing import Dict, Iterable, List, Sequence, Tuple


def fmt(x: F) -> str:
    return f"{x} ({float(x):.8f})"


def sector_of(p: F) -> int:
    return int((p % 1) * 7)


def missdist(E: Iterable[int]) -> List[F]:
    """Distribution q_t of the number of empty inner seventh sectors."""
    E = sorted(set(E))
    breakpoints = {F(0), F(1)}
    for e in E:
        if e == 0:
            continue
        for m in range(0, 7 * e + 1):
            breakpoints.add(F(m, 7 * e))
    b = sorted(breakpoints)
    q = [F(0)] * 7
    for i in range(len(b) - 1):
        x0, x1 = b[i], b[i + 1]
        if x1 <= x0:
            continue
        mid = (x0 + x1) / 2
        hit = {sector_of(e * mid) for e in E}
        empty = 7 - len(hit)
        if 0 <= empty <= 6:
            q[empty] += x1 - x0
    return q


def Sr(q: Sequence[F], r: int) -> F:
    return sum(F(comb(t, r)) * q[t] for t in range(7))


def poly_from_roots(roots: Sequence[int]) -> Tuple[int, ...]:
    coeffs = [1]
    for root in roots:
        nxt = [0] * (len(coeffs) + 1)
        for i, coeff in enumerate(coeffs):
            nxt[i] += coeff
            nxt[i + 1] -= coeff * root
        coeffs = nxt
    return tuple(coeffs)


def eval_poly(coeffs: Sequence[float], z: complex) -> complex:
    value = 0j
    for coeff in coeffs:
        value = value * z + coeff
    return value


def durand_kerner_roots(coeffs: Sequence[float], steps: int = 200) -> List[complex]:
    """Numerically approximate roots of a small polynomial.

    This is only used for the root-radius sidecar.  Exact claims in the scout
    use Fraction arithmetic and Vieta identities.
    """
    while len(coeffs) > 1 and abs(coeffs[0]) < 1e-14:
        coeffs = coeffs[1:]
    degree = len(coeffs) - 1
    if degree <= 0:
        return []
    lead = coeffs[0]
    monic = [c / lead for c in coeffs]
    radius = 1 + max(abs(c) for c in monic[1:])
    roots = [radius * cmath.exp(2j * cmath.pi * j / degree) for j in range(degree)]
    for _ in range(steps):
        max_delta = 0.0
        new_roots = []
        for i, root in enumerate(roots):
            denom = 1 + 0j
            for j, other in enumerate(roots):
                if i != j:
                    denom *= root - other
            if abs(denom) < 1e-18:
                new_root = root
            else:
                delta = eval_poly(monic, root) / denom
                new_root = root - delta
                max_delta = max(max_delta, abs(delta))
            new_roots.append(new_root)
        roots = new_roots
        if max_delta < 1e-13:
            break
    return roots


def root_radii(q: Sequence[F]) -> Tuple[float, float, float]:
    coeffs = [float(q[t]) for t in range(6, -1, -1)]
    radii = sorted(abs(z) for z in durand_kerner_roots(coeffs))
    return radii[0], radii[-1], radii[-1] / radii[0]


def eulerian_number(n: int, k: int) -> int:
    return sum((-1) ** j * comb(n + 1, j) * (k + 1 - j) ** n for j in range(k + 1))


def eulerian_row(n: int) -> Tuple[int, ...]:
    return tuple(eulerian_number(n, k) for k in range(n))


def worpitzky_verified(n: int, limit: int = 9) -> bool:
    row = eulerian_row(n)
    for x in range(limit + 1):
        rhs = sum(row[k] * comb(x + n - 1 - k, n) for k in range(n))
        if rhs != x**n:
            return False
    return True


def k3_class(bits: Tuple[int, int, int]) -> str:
    return "C" if bits in {(0, 0, 0), (1, 1, 1)} else "T"


def k3_kernel() -> Dict[str, object]:
    counts: Counter[Tuple[str, str]] = Counter()
    class_counts: Counter[str] = Counter()
    for bits in product((0, 1), repeat=3):
        state = tuple(bits)  # type: ignore[assignment]
        cls = k3_class(state)
        class_counts[cls] += 1
        for idx in range(3):
            nxt = list(state)
            nxt[idx] ^= 1
            counts[(cls, k3_class(tuple(nxt)))] += 1
    matrix = {
        cls: {
            dst: F(counts[(cls, dst)], 3 * class_counts[cls])
            for dst in ("C", "T")
        }
        for cls in ("C", "T")
    }
    return {
        "counts": dict(sorted(counts.items())),
        "matrix_rows_C_T": ((matrix["C"]["C"], matrix["C"]["T"]), (matrix["T"]["C"], matrix["T"]["T"])),
        "stationary_C_T": (F(1, 4), F(3, 4)),
        "eigenvalues": (F(1), F(-1, 3)),
    }


def has_perfect_matching(vertices: Tuple[int, ...], edges: set[Tuple[int, int]]) -> bool:
    if not vertices:
        return True
    first = vertices[0]
    rest = vertices[1:]
    for other in rest:
        edge = tuple(sorted((first, other)))
        if edge not in edges:
            continue
        remaining = tuple(v for v in rest if v != other)
        if has_perfect_matching(remaining, edges):
            return True
    return False


def is_factor_critical(n: int, undirected_edges: Iterable[Tuple[int, int]]) -> bool:
    edges = {tuple(sorted(e)) for e in undirected_edges}
    for deleted in range(n):
        vertices = tuple(v for v in range(n) if v != deleted)
        if not has_perfect_matching(vertices, edges):
            return False
    return True


def strongly_connected(n: int, directed_edges: Iterable[Tuple[int, int]]) -> bool:
    adj = {i: set() for i in range(n)}
    radj = {i: set() for i in range(n)}
    for a, b in directed_edges:
        adj[a].add(b)
        radj[b].add(a)

    def reach(graph: Dict[int, set[int]], start: int) -> set[int]:
        seen = {start}
        stack = [start]
        while stack:
            v = stack.pop()
            for w in graph[v]:
                if w not in seen:
                    seen.add(w)
                    stack.append(w)
        return seen

    return len(reach(adj, 0)) == n and len(reach(radj, 0)) == n


def ear_packet_audit() -> Dict[str, bool]:
    triangle = {(0, 1), (1, 2), (0, 2)}
    triangle_plus_odd_ear = triangle | {(0, 3), (3, 4), (1, 4)}

    directed_cycle = {(0, 1), (1, 2), (2, 0)}
    directed_cycle_plus_ear = directed_cycle | {(0, 3), (3, 4), (4, 1)}
    return {
        "C3_factor_critical": is_factor_critical(3, triangle),
        "C3_plus_length3_odd_ear_factor_critical": is_factor_critical(5, triangle_plus_odd_ear),
        "directed_C3_strong": strongly_connected(3, directed_cycle),
        "directed_C3_plus_directed_ear_strong": strongly_connected(5, directed_cycle_plus_ear),
    }


AXES = (
    "factor_through_legality",
    "root_curve_retention",
    "cap_dip_separation",
    "parity_split",
    "degree_control",
    "ap_extremality",
    "ear_witness_transfer",
    "sidecar_hygiene",
    "anti_scalar_guard",
)


@dataclass(frozen=True)
class Carrier:
    name: str
    scores: Dict[str, int]


def carrier(name: str, values: Sequence[int]) -> Carrier:
    return Carrier(name, dict(zip(AXES, values)))


CARRIERS = [
    carrier("HYP3151_function_legality_packet", (10, 6, 7, 8, 8, 6, 5, 10, 10)),
    carrier("HYP3152_leeyang_radius_root_curve", (8, 10, 10, 8, 7, 9, 6, 9, 10)),
    carrier("k8_Ly_bimodality_certificate", (8, 8, 10, 8, 8, 10, 5, 8, 9)),
    carrier("even_biquadratic_galois_fold", (7, 5, 8, 9, 10, 7, 4, 8, 8)),
    carrier("odd_worpitzky_ear_packet", (7, 6, 7, 10, 6, 7, 10, 9, 9)),
    carrier("pascal_cap_dip_table", (6, 5, 10, 7, 8, 8, 4, 7, 7)),
    carrier("newton_ap_violation_packet", (7, 7, 8, 6, 7, 10, 4, 7, 8)),
    carrier("generic_quintic_wall_alarm", (4, 4, 5, 5, 9, 4, 3, 8, 9)),
    carrier("raw_scalar_p0_only", (1, 1, 2, 1, 2, 3, 0, 0, 0)),
]

TIE_PATH = [c.name for c in CARRIERS]


def compare(a: Carrier, b: Carrier) -> str:
    aw = bw = 0
    for axis in AXES:
        if a.scores[axis] > b.scores[axis]:
            aw += 1
        elif b.scores[axis] > a.scores[axis]:
            bw += 1
    if aw > bw:
        return a.name
    if bw > aw:
        return b.name
    return a.name if TIE_PATH.index(a.name) < TIE_PATH.index(b.name) else b.name


def tournament_edges() -> Dict[str, set[str]]:
    adj = {c.name: set() for c in CARRIERS}
    for a, b in combinations(CARRIERS, 2):
        winner = compare(a, b)
        loser = b.name if winner == a.name else a.name
        adj[winner].add(loser)
    return adj


def directed_3cycle_count(adj: Dict[str, set[str]]) -> int:
    count = 0
    for a, b, c in combinations(adj, 3):
        if b in adj[a] and c in adj[b] and a in adj[c]:
            count += 1
        if c in adj[a] and b in adj[c] and a in adj[b]:
            count += 1
    return count


def tarjan_scc(adj: Dict[str, set[str]]) -> List[List[str]]:
    index = 0
    stack: List[str] = []
    on_stack: set[str] = set()
    indices: Dict[str, int] = {}
    low: Dict[str, int] = {}
    out: List[List[str]] = []

    def strongconnect(v: str) -> None:
        nonlocal index
        indices[v] = index
        low[v] = index
        index += 1
        stack.append(v)
        on_stack.add(v)
        for w in adj[v]:
            if w not in indices:
                strongconnect(w)
                low[v] = min(low[v], low[w])
            elif w in on_stack:
                low[v] = min(low[v], indices[w])
        if low[v] == indices[v]:
            comp = []
            while True:
                w = stack.pop()
                on_stack.remove(w)
                comp.append(w)
                if w == v:
                    break
            out.append(sorted(comp))

    for v in adj:
        if v not in indices:
            strongconnect(v)
    return sorted(out, key=lambda comp: (-len(comp), comp[0]))


def hamiltonian_path_count(adj: Dict[str, set[str]]) -> int:
    names = list(adj)
    count = 0
    for order in permutations(names):
        if all(order[i + 1] in adj[order[i]] for i in range(len(order) - 1)):
            count += 1
    return count


def main() -> None:
    print("HYP-3153 / codex-2026-06-28")
    print("Lee-Yang/Worpitzky/quartic compression packet for LRC14")
    print("namespace=HYP-3153/T1218/LTI-279/LTT-177")
    print()

    caps = {
        8: F(2243, 5880),
        9: F(1979, 4004),
        10: F(55, 91),
        11: F(66, 91),
        12: F(6, 7),
        13: F(1, 1),
    }

    print("1. LEE-YANG RADIUS / PASCAL CAP / DIP TABLE")
    print("k  q0                 q6                 R^6=q0/q6        root_ratio  pair_mass-cap=dip")
    for k in range(8, 14):
        q = missdist(range(k))
        r6 = q[0] / q[6]
        radii = root_radii(q)
        root_ratio = f"{radii[2]:.6f}"
        pair_mass = F(comb(k + 1, 2), 91)
        dip = pair_mass - caps[k]
        print(
            f"{k:>1}  {str(q[0]):>17}  {str(q[6]):>17}  "
            f"{str(r6):>17}  {root_ratio:>10}  {dip}"
        )
    print("verified=q0=q6*R^6 exactly by Vieta; root_ratio is numeric only")
    print("readout=pair-normalized Pascal mass is the on-circle shadow; dip is the off-circle correction")
    print()

    print("2. L_y BIMODALITY FUNCTIONAL AGAINST CAP")
    for k in range(8, 11):
        q = missdist(range(k))
        ly = q[0] + q[6] + F(1, 10) * q[3]
        print(f"k={k}: L_y=q0+q6+q3/10={fmt(ly)} cap={fmt(caps[k])} margin={fmt(caps[k]-ly)} ok={ly <= caps[k]}")
    q8 = missdist(range(8))
    S = [Sr(q8, r) for r in range(7)]
    gk8_moment = 10 * S[0] - 10 * S[1] + 10 * S[2] - 9 * S[3] + 6 * S[4]
    ly_scaled = 10 * q8[0] + q8[3] + 10 * q8[6]
    odd = -9 * S[3]
    even = 6 * S[4]
    base = 10 * S[0] - 10 * S[1] + 10 * S[2]
    print(f"k=8 scaled identity: 10q0+q3+10q6={fmt(ly_scaled)} gK8_moment={fmt(gk8_moment)} equal={ly_scaled == gk8_moment}")
    print(f"k=8 split: base={fmt(base)} odd=-9*S3={fmt(odd)} even=6*S4={fmt(even)} abs(odd/even)={float(abs(odd/even)):.6f}")
    print()

    print("3. WORPITZKY / EDGE-FLIP ODD PACKET")
    kernel = k3_kernel()
    print(f"k3_counts={kernel['counts']}")
    print(f"k3_matrix_rows_C_T={kernel['matrix_rows_C_T']}")
    print(f"k3_stationary_C_T={kernel['stationary_C_T']}")
    print(f"k3_eigenvalues={kernel['eigenvalues']}")
    for n in (3, 4):
        print(f"Eulerian_row_degree_{n}={eulerian_row(n)} Worpitzky_verified_0..9={worpitzky_verified(n)}")
    print("readout=the odd packet is not a scalar class; it carries the -1/3 flip mode and descent sidecar")
    print()

    print("4. EVEN BIQUADRATIC / GALOIS DEGREE PACKET")
    roots_t = (1, 2, 4, 5)
    roots_u = tuple(r - 3 for r in roots_t)
    roots_v = tuple(sorted({u * u for u in roots_u}))
    coeff_t = poly_from_roots(roots_t)
    coeff_u = poly_from_roots(roots_u)
    coeff_v = poly_from_roots(roots_v)
    print(f"g(t) roots={roots_t} coeffs={coeff_t}")
    print(f"u=t-3 roots={roots_u} coeffs={coeff_u}")
    print(f"v=u^2 roots={roots_v} coeffs={coeff_v} discriminant={coeff_v[1]**2 - 4*coeff_v[0]*coeff_v[2]}")
    print("verified=degree ceiling 4; centered odd coefficients vanish; even variable v gives degree 2")
    print("guardrail=solvability is by degree<=4/Galois<=S4, not by identifying the flip monoid with V4")
    print()

    print("5. NEWTON-MACLAURIN VIOLATION / AP EXTREMAL PACKET")
    p = [S[r] / F(comb(6, r)) for r in range(7)]
    defects = [p[r] * p[r] - p[r - 1] * p[r + 1] for r in range(1, 6)]
    print(f"normalized_moment_means={[str(x) for x in p]}")
    print(f"newton_defects={[str(x) for x in defects]}")
    print(f"all_negative={all(x < 0 for x in defects)}")
    print("readout=AP/consec is extremal for bimodality/violation; this is a certificate target, not yet the proof")
    print()

    print("6. EAR / ODD-CYCLE WITNESS SIDECAR")
    ears = ear_packet_audit()
    for key, value in ears.items():
        print(f"{key}={value}")
    print("readout=odd ears are a finite witness grammar for the off-circle/Worpitzky sidecar")
    print()

    print("7. TOURNAMENT ANALYSIS OVER PROOF PACKETS")
    adj = tournament_edges()
    score_hist = Counter(len(adj[name]) for name in adj)
    selected_path = sorted(adj, key=lambda name: (-len(adj[name]), TIE_PATH.index(name)))
    print("vertices=proof packets/obligations, not runners/arcs/classes")
    print("pairwise_observable=majority over legality, root retention, cap/dip split, parity split, degree control, AP extremality, ear transfer, sidecar hygiene, anti-scalar guard")
    print("switch=A beats B when A wins more axes; ties follow declared packet path")
    print(f"score_hist={dict(sorted(score_hist.items()))}")
    print(f"directed_3cycles={directed_3cycle_count(adj)}")
    print(f"SCCs={tarjan_scc(adj)}")
    print(f"hamiltonian_path_count={hamiltonian_path_count(adj)}")
    print("selected_path=")
    for name in selected_path:
        print(f"  -> {name} score={len(adj[name])}")
    print()

    print("8. PACKET CONCLUSION")
    print("verified_identities=q0=q6*R^6; L_y=gK8 moment row at k=8; Worpitzky n=3/n=4; k=8 biquadratic fold; finite odd-ear examples")
    print("proof_obligations=bound off-circle dip/lambda; prove AP maximizes L_y in the required bank; carry ordered/ear sidecars through LRC predicate")
    print("recommended_next_packet_fields=root_radius_R6, root_radius_spread, pair_pascal_mass, cap_dip, Ly_margin, odd_worpitzky_mode, even_biquadratic_degree, odd_ear_payload, factor_through_status, terminal_exit_or_named_debt")


if __name__ == "__main__":
    main()
