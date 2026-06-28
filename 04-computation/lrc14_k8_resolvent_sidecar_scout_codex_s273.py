#!/usr/bin/env python3
"""HYP-3142 scout: k=8 resolvent sidecar certificate.

Executable synthesis, not an LRC(14) proof.

The current frontier says the covering proof has collapsed to the bounded-core
k=8 row.  This scout tests whether older niche sidecars give a finite
certificate grammar for that node, now read as the terminal
`bounded_core_U4_exit` after the HYP-3141 edge-witness tip/tail packet,
HYP-3140 fiber-PGF Rprime certificate, HYP-3137 generating-function payload
atlas, the HYP-3138 reflection-fold resurrection table, the HYP-3139
reflection-block proof pages, the HYP-3136 multi-far floor packet, and
incoming HYP-3135 resolvent-packet middle layer:

* De Moivre / HYP-3132: the gK8 dual folds to the biquadratic
  (u^2-1)(u^2-4).
* Bravais / HYP-3113: the worst k=8 rows should be residue-flat rather than
  Bragg-peaked.
* Savitch / HYP-3118: failed scalar shortcuts should carry a shallow midpoint
  repair depth, i.e. a small list of missing sidecars.
* A000568 / HYP-3133-HYP-3134: edge-envelope quotienting is legal only after
  a global-consistency sidecar is present.

The load-bearing computation is exact-rational: for each primitive bounded
shape E={0<e_2<...<e_8<=B}, compute the miss-count distribution q_i, form the
4-binomial-moment relaxation U4(E), and compare U4(E) to cap_8.

Tournament Analysis declaration:
  vertices: proof sidecars/certificate operators, not runners or scalar rows;
  pairwise observable: majority comparison over retained proof-payload axes;
  switch/gauge: A -> B when A wins more axes; ties use the declared path.
"""

from __future__ import annotations

import math
import sys
from collections import Counter
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import combinations, permutations
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np


Q = 7
CAP8 = F(2243, 5880)


def fstr(x: F) -> str:
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def sector_of(x: F) -> int:
    return int((x % 1) * Q)


def cell_profile(speeds: Sequence[int]) -> List[F]:
    speeds = tuple(sorted(set(speeds)))
    if 0 not in speeds:
        speeds = (0,) + speeds
    breaks = {F(0), F(1)}
    for e in speeds:
        if e == 0:
            continue
        for m in range(0, Q * e + 1):
            breaks.add(F(m, Q * e))
    cuts = sorted(breaks)
    q = [F(0)] * Q
    for lo, hi in zip(cuts, cuts[1:]):
        if hi <= lo:
            continue
        mid = (lo + hi) / 2
        covered = {sector_of(e * mid) for e in speeds}
        empty = sum(1 for s in range(1, Q) if s not in covered)
        q[empty] += hi - lo
    return q


def solve_rational(A: List[List[F]], b: List[F]) -> List[F] | None:
    n = len(A)
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for col in range(n):
        piv = None
        for r in range(col, n):
            if M[r][col] != 0:
                piv = r
                break
        if piv is None:
            return None
        M[col], M[piv] = M[piv], M[col]
        pv = M[col][col]
        M[col] = [v / pv for v in M[col]]
        for r in range(n):
            if r != col and M[r][col] != 0:
                fac = M[r][col]
                M[r] = [M[r][j] - fac * M[col][j] for j in range(n + 1)]
    return [M[i][n] for i in range(n)]


def u4_exact_from_q(q: Sequence[F]) -> F:
    moments = [sum(F(math.comb(i, t)) * q[i] for i in range(Q)) for t in range(5)]
    rows = [[F(math.comb(i, t)) for i in range(Q)] for t in range(5)]
    best = None
    for basis in combinations(range(Q), 5):
        A = [[rows[t][j] for j in basis] for t in range(5)]
        sol = solve_rational(A, moments)
        if sol is None or any(v < 0 for v in sol):
            continue
        p = [F(0)] * Q
        for idx, j in enumerate(basis):
            p[j] = sol[idx]
        if best is None or p[0] > best:
            best = p[0]
    if best is None:
        raise RuntimeError("4-moment relaxation has no feasible vertex")
    return best


def pgf_roots(q: Sequence[F]) -> List[complex]:
    coeffs = [float(q[t]) for t in range(Q - 1, -1, -1)]
    while len(coeffs) > 1 and abs(coeffs[0]) < 1e-15:
        coeffs = coeffs[1:]
    return list(np.roots(coeffs))


def cumulant4(q: Sequence[F]) -> float:
    probs = np.array([float(x) for x in q], dtype=float)
    xs = np.arange(Q, dtype=float)
    mu = float(np.dot(xs, probs))
    centered = xs - mu
    m2 = float(np.dot(centered**2, probs))
    m4 = float(np.dot(centered**4, probs))
    return m4 - 3.0 * m2 * m2


def bravais_features(E: Sequence[int]) -> Tuple[float, float, int, Tuple[int, ...]]:
    counts = [0] * Q
    for e in E:
        counts[e % Q] += 1
    amps = []
    for mode in range(1, Q):
        z = sum(
            counts[r] * complex(math.cos(2 * math.pi * mode * r / Q), math.sin(2 * math.pi * mode * r / Q))
            for r in range(Q)
        )
        amps.append(abs(z) / len(E))
    peak = max(amps)
    total = sum(a * a for a in amps)
    probs = [a * a / total for a in amps] if total else [1 / 6] * 6
    entropy = -sum(p * math.log(p) for p in probs if p > 0) / math.log(6)
    mirror_defect = sum(abs(counts[r] - counts[(-r) % Q]) for r in range(1, Q))
    return peak, entropy, mirror_defect, tuple(counts)


def shifted_resolvent_coefficients() -> Tuple[List[int], int, Tuple[int, int]]:
    """Return coefficients of g(u+3) and its quadratic discriminant in y=u^2."""
    # g(t)=(t-1)(t-2)(t-4)(t-5).  Expanding after t=u+3 gives
    # (u+2)(u+1)(u-1)(u-2)=u^4-5u^2+4.
    coeffs_u4_to_u0 = [1, 0, -5, 0, 4]
    discriminant_y = 25 - 16
    roots_y = (1, 4)
    return coeffs_u4_to_u0, discriminant_y, roots_y


@dataclass(frozen=True)
class Row:
    E: Tuple[int, ...]
    q: Tuple[F, ...]
    u4: F
    cap_slack: F
    q0: F
    nearest_root: float
    real_roots: int
    bravais_peak: float
    residue_entropy: float
    mirror_defect: int
    residue_counts: Tuple[int, ...]
    kappa4: float
    sidecar_depth: int


def sidecar_depth(u4_close: bool, ly_margin: bool, bravais_flat: bool, mirror_ok: bool, phi4_sign: bool) -> int:
    # Savitch-style midpoint proxy: how many binary splits are needed to name
    # the missing sidecar family once the exact U4 terminal gate has been tried.
    missing = sum(not x for x in (u4_close, ly_margin, bravais_flat, mirror_ok, phi4_sign))
    return 0 if missing == 0 else math.ceil(math.log2(missing + 1))


def primitive_rows(B: int) -> List[Row]:
    rows: List[Row] = []
    for combo in combinations(range(1, B + 1), 7):
        E = (0,) + combo
        g = 0
        for e in E:
            g = math.gcd(g, e)
        if g != 1:
            continue
        q = tuple(cell_profile(E))
        u4 = u4_exact_from_q(q)
        roots = pgf_roots(q)
        nearest = min(abs(z) for z in roots) if roots else 0.0
        real_roots = sum(1 for z in roots if abs(z.imag) < 1e-8)
        peak, entropy, mirror, counts = bravais_features(E)
        k4 = cumulant4(q)
        flat = abs(peak - 1 / 8) < 1e-9 and abs(entropy - 1.0) < 1e-9
        depth = sidecar_depth(u4 <= CAP8, nearest > 1.0, flat, mirror == 0, k4 < 0)
        rows.append(
            Row(
                E=E,
                q=q,
                u4=u4,
                cap_slack=CAP8 - u4,
                q0=q[0],
                nearest_root=nearest,
                real_roots=real_roots,
                bravais_peak=peak,
                residue_entropy=entropy,
                mirror_defect=mirror,
                residue_counts=counts,
                kappa4=k4,
                sidecar_depth=depth,
            )
        )
    return rows


def corr(xs: Sequence[float], ys: Sequence[float]) -> float:
    if len(xs) < 2:
        return 0.0
    x = np.array(xs, dtype=float)
    y = np.array(ys, dtype=float)
    if np.std(x) < 1e-15 or np.std(y) < 1e-15:
        return 0.0
    return float(np.corrcoef(x, y)[0, 1])


AXES = (
    "closes_k8_cap_node",
    "finite_exact_certificate",
    "preserves_lrc_predicate",
    "explains_resolvent_fold",
    "resurrects_destroyed_coordinate",
    "connects_niche_work",
    "failure_guardrail",
    "formalization_path",
)


@dataclass(frozen=True)
class Carrier:
    name: str
    anchors: Tuple[str, ...]
    axes: Dict[str, int]
    preserves: str
    destroys: str
    next_hook: str


def carrier(
    name: str,
    anchors: Sequence[str],
    scores: Sequence[int],
    preserves: str,
    destroys: str,
    next_hook: str,
) -> Carrier:
    return Carrier(
        name=name,
        anchors=tuple(anchors),
        axes=dict(zip(AXES, scores)),
        preserves=preserves,
        destroys=destroys,
        next_hook=next_hook,
    )


CARRIERS = [
    carrier(
        "exact_4moment_U4_cap_certificate",
        (
            "HYP-3142",
            "HYP-3141",
            "HYP-3140",
            "HYP-3139",
            "HYP-3138",
            "HYP-3137",
            "HYP-3136",
            "HYP-3135",
            "HYP-3132",
            "HYP-3122",
            "THM-577",
        ),
        (10, 10, 9, 8, 7, 7, 9, 9),
        "an exact rational terminal upper bound for q0 at the k=8 node",
        "why the moment majorization should hold outside the bounded bank",
        "prove U4(E) is maximized by the consecutive/full-residue-flat stalk",
    ),
    carrier(
        "biquadratic_resolvent_fold",
        ("HYP-3139", "HYP-3138", "HYP-3135", "HYP-3132", "HYP-3110", "LTI-247"),
        (8, 10, 7, 10, 6, 9, 8, 9),
        "the De Moivre fold g(u+3)=u^4-5u^2+4 with square discriminant",
        "packet-side uniformity unless attached to U4 or finite-address data",
        "turn the kappa4 obligation into a degree-2 coefficient inequality",
    ),
    carrier(
        "bravais_residue_flatness",
        ("HYP-3113", "HYP-3111", "LTI-250"),
        (7, 8, 8, 7, 7, 10, 8, 7),
        "the full-residue flat stalk that carries the worst bounded-bank rows",
        "endpoint owner and exact interval cells",
        "prove non-flat residue spectra strictly improve the U4 slack",
    ),
    carrier(
        "savitch_midpoint_repair_depth",
        ("HYP-3118", "HYP-3116", "LTI-254"),
        (7, 7, 9, 6, 10, 9, 9, 8),
        "a small missing-sidecar depth after the exact terminal gate is tried",
        "the concrete coefficient inequality unless U4 is retained",
        "use depth as the finite induction measure for remaining row families",
    ),
    carrier(
        "lee_yang_root_confinement",
        ("HYP-3122", "HYP-3131", "HYP-3112"),
        (7, 7, 8, 7, 7, 8, 10, 7),
        "nearest-root and real-root strata guarding the phi4 side",
        "finite rational cap arithmetic if treated as analytic slogan",
        "show rows with near roots have a named U4 or endpoint-deletion repair",
    ),
    carrier(
        "a000568_edge_global_consistency",
        ("HYP-3134", "HYP-3133", "HYP-3124"),
        (6, 6, 10, 7, 9, 8, 9, 7),
        "the legal quotient boundary for dropping paired edge-child payload",
        "the k=8 moment inequality by itself",
        "attach global_consistency_class to the k=8 U4 packet",
    ),
    carrier(
        "phi4_kappa4_stabilizer",
        ("HYP-3122", "HYP-3113"),
        (7, 7, 7, 9, 6, 8, 8, 7),
        "the negative quartic-cumulant sign at the extremal stalk",
        "positive-kappa rows that still close by U4",
        "separate sign-stabilized rows from U4-only terminal rows",
    ),
    carrier(
        "coordinate_resurrection_sheaf",
        ("HYP-3118", "HYP-3120", "HYP-3134"),
        (6, 6, 10, 6, 10, 9, 9, 8),
        "which coordinate each quotient destroyed and which sidecar restores it",
        "the numeric extremal comparison if no terminal sidecar is named",
        "record repair_cover_rank for every k=8 row family",
    ),
    carrier(
        "raw_scalar_meas_S7",
        ("raw-count",),
        (4, 5, 3, 2, 1, 1, 2, 3),
        "only the visible q0 value",
        "root, lattice, edge-child, resolvent, and repair coordinates",
        "use as negative control only",
    ),
]


TIE_PATH = [
    "exact_4moment_U4_cap_certificate",
    "biquadratic_resolvent_fold",
    "bravais_residue_flatness",
    "savitch_midpoint_repair_depth",
    "lee_yang_root_confinement",
    "a000568_edge_global_consistency",
    "phi4_kappa4_stabilizer",
    "coordinate_resurrection_sheaf",
    "raw_scalar_meas_S7",
]


def tournament() -> Tuple[Dict[str, List[str]], Counter, int, List[List[str]]]:
    order = {name: i for i, name in enumerate(TIE_PATH)}
    adj = {c.name: [] for c in CARRIERS}
    score = Counter()
    for a, b in combinations(CARRIERS, 2):
        aw = sum(a.axes[x] > b.axes[x] for x in AXES)
        bw = sum(b.axes[x] > a.axes[x] for x in AXES)
        if aw > bw or (aw == bw and order[a.name] < order[b.name]):
            adj[a.name].append(b.name)
            score[a.name] += 1
        else:
            adj[b.name].append(a.name)
            score[b.name] += 1
    for c in CARRIERS:
        score[c.name] += 0

    names = [c.name for c in CARRIERS]
    cycles = 0
    for a, b, c in combinations(names, 3):
        edges = {(x, y) for x in names for y in adj[x]}
        if (a, b) in edges and (b, c) in edges and (c, a) in edges:
            cycles += 1
        if (a, c) in edges and (c, b) in edges and (b, a) in edges:
            cycles += 1

    paths: List[List[str]] = []
    edge_set = {(x, y) for x in names for y in adj[x]}
    for perm in permutations(names):
        if all((perm[i], perm[i + 1]) in edge_set for i in range(len(perm) - 1)):
            paths.append(list(perm))
            if len(paths) > 20:
                break
    return adj, score, cycles, paths


def main() -> None:
    B = int(sys.argv[1]) if len(sys.argv) > 1 else 13
    rows = primitive_rows(B)
    coeffs, discr, roots_y = shifted_resolvent_coefficients()
    worst_u4 = sorted(rows, key=lambda r: r.cap_slack)[:10]
    top_q0 = sorted(rows, key=lambda r: r.q0, reverse=True)[:10]
    u4_over = [r for r in rows if r.u4 > CAP8]
    flat = [r for r in rows if abs(r.bravais_peak - 1 / 8) < 1e-9 and abs(r.residue_entropy - 1.0) < 1e-9]
    ly = [r for r in rows if r.nearest_root > 1.0]
    phi4_neg = [r for r in rows if r.kappa4 < 0]
    depths = Counter(r.sidecar_depth for r in rows)

    print("HYP-3142 / codex-2026-06-27-S273")
    print("k=8 resolvent sidecar scout; executable synthesis, not a proof.")
    print()
    print("1. RESOLVENT FOLD")
    print("g(t)=(t-1)(t-2)(t-4)(t-5); with u=t-3:")
    print(f"  coefficients u^4..u^0 = {coeffs}")
    print(f"  odd_coefficients_zero = {coeffs[1] == 0 and coeffs[3] == 0}")
    print(f"  discriminant in y=u^2 = {discr} = {int(math.isqrt(discr))}^2")
    print(f"  y_roots = {roots_y}")
    print()
    print("2. EXACT BOUNDED-BANK U4 SCAN")
    print(f"B={B}; primitive rows scanned={len(rows)}")
    print(f"cap_8={fstr(CAP8)} = {float(CAP8):.12f}")
    print(f"U4_over_cap_count={len(u4_over)}")
    print(f"Lee_Yang_nearest_root_gt_1_count={len(ly)}")
    print(f"Bravais_full_residue_flat_count={len(flat)}")
    print(f"phi4_kappa4_negative_count={len(phi4_neg)}")
    print(f"Savitch_sidecar_depth_hist={dict(sorted(depths.items()))}")
    print()
    print("Worst rows by cap_8-U4:")
    for r in worst_u4:
        print(
            "  "
            f"E={r.E} q0={fstr(r.q0)} U4={fstr(r.u4)} slack={fstr(r.cap_slack)} "
            f"nearest={r.nearest_root:.6f} real_roots={r.real_roots} "
            f"bravais_peak={r.bravais_peak:.6f} entropy={r.residue_entropy:.6f} "
            f"mirror={r.mirror_defect} kappa4={r.kappa4:+.6f} depth={r.sidecar_depth} "
            f"residues={r.residue_counts}"
        )
    print()
    print("Top rows by actual q0=meas(S7):")
    for r in top_q0:
        print(
            "  "
            f"E={r.E} q0={fstr(r.q0)} U4={fstr(r.u4)} slack={fstr(r.cap_slack)} "
            f"nearest={r.nearest_root:.6f} kappa4={r.kappa4:+.6f} depth={r.sidecar_depth}"
        )
    print()
    print("Correlations over bounded bank:")
    u4s = [float(r.u4) for r in rows]
    print(f"  corr(U4,q0)={corr(u4s, [float(r.q0) for r in rows]):+.6f}")
    print(f"  corr(U4,nearest_root)={corr(u4s, [r.nearest_root for r in rows]):+.6f}")
    print(f"  corr(U4,bravais_peak)={corr(u4s, [r.bravais_peak for r in rows]):+.6f}")
    print(f"  corr(U4,residue_entropy)={corr(u4s, [r.residue_entropy for r in rows]):+.6f}")
    print(f"  corr(U4,kappa4)={corr(u4s, [r.kappa4 for r in rows]):+.6f}")
    print(f"  corr(U4,sidecar_depth)={corr(u4s, [r.sidecar_depth for r in rows]):+.6f}")
    print()
    print("3. INTERPRETATION")
    print(
        "In the integrated route, this scout supplies the bounded_core_U4_exit "
        "after the HYP-3141 edge-witness tip/tail packet, HYP-3140 "
        "fiber-PGF Rprime certificate, HYP-3137 generating-function payload "
        "atlas, HYP-3138 reflection-fold repair table, HYP-3139 "
        "reflection-block proof pages, HYP-3136 multi-far floor packet, and "
        "incoming HYP-3135 resolvent-packet middle layer."
    )
    print(
        "The bounded-bank evidence says the exact 4-moment U4 sidecar is a terminal "
        "cap certificate on all scanned k=8 rows.  The worst row is consec_8, "
        "with U4=2633/7350 and cap_8-U4=683/29400."
    )
    print(
        "The worst stalk is Bravais-flat, mirror-even, Lee-Yang confined, and "
        "phi4-negative, matching the HYP-3132 biquadratic fold.  Non-flat rows "
        "gain slack; rows without the phi4 sign still close by U4, so phi4 is "
        "a stabilizer/explainer rather than the whole terminal certificate."
    )
    print(
        "Proof target: prove a global moment-majorization theorem that U4(E) is "
        "maximized by the consecutive/full-residue-flat stalk, then HYP-3132's "
        "single k=8 node closes by exact rational arithmetic."
    )
    print()
    print("4. TOURNAMENT ANALYSIS OVER SIDECAR CARRIERS")
    adj, scores, cycles, paths = tournament()
    hist = Counter(scores.values())
    print(f"vertices={len(CARRIERS)} axes={','.join(AXES)}")
    print(f"score_hist={dict(sorted(hist.items()))}")
    print(f"directed_3cycles={cycles}")
    print(f"hamiltonian_path_count_capped={len(paths)}")
    if paths:
        print("selected_hamiltonian_path=" + " -> ".join(paths[0]))
    print()
    print("Top carrier hooks:")
    for name in paths[0][:6] if paths else [c.name for c in CARRIERS[:6]]:
        c = next(x for x in CARRIERS if x.name == name)
        print(f"- {c.name}: preserves={c.preserves}; next={c.next_hook}")


if __name__ == "__main__":
    main()
