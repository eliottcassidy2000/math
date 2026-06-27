#!/usr/bin/env python3
"""HYP-3125 scout: multi-far Rprime floor through edge witnesses.

This is an executable synthesis scout, not an LRC14 proof.

Goal: connect the open covering-case floor

    Rprime = mu(R-safe and Q-safe) / (mu(R-safe) * mu(Q-safe)) >= c

for the r=2..6 multi-far branch to the completed HYP-3124 edge-witness
recursion.  The experiment treats Rprime as a diagonal edge-sector mass:
the tail is an R-safe packet, the tip is a Q-safe packet, and the directed
edge tail -> tip carries a four-sector deck plus recursive deletion children.

Tournament Analysis declaration:
  vertices: proof carriers / repair operators, not runners, arcs, residues,
            primes, Gaussian kernels, or raw scalar ratios;
  pairwise observable: majority comparison over retained proof-payload axes;
  switch/gauge: A beats B if it wins more axes, with a declared Hamiltonian
                tie path;
  tie Hamiltonian path: Rprime floor, two-ended edge witness, Gaussian
                        smoothing, EH-style level of distribution, Asano
                        contraction, phi4 stabilizer, Cech finite ruler,
                        chiral guard, H7 lift, raw Bonferroni.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from itertools import combinations
from math import sqrt
from typing import Dict, Iterable, List, Sequence, Tuple


THRESHOLD = 1.0 / 14.0
GRID_STEPS = 14 * 12000

AXES = (
    "lrc_predicate_retention",
    "attacks_multifar_floor",
    "tail_tip_recursion",
    "decorrelation_power",
    "zero_free_contraction",
    "finite_packet_fields",
    "formalization_readiness",
    "failure_guard",
)


@dataclass(frozen=True)
class Carrier:
    name: str
    anchors: Tuple[str, ...]
    axes: Dict[str, int]
    preserves: str
    destroys: str
    next_hook: str


@dataclass(frozen=True)
class RowAudit:
    name: str
    speeds: Tuple[int, ...]
    r: int
    R_size: int
    R_safe: float
    Q_safe: float
    joint: float
    bonferroni_margin: float
    Rprime: float
    tail_child_min: float
    tail_child_max: float
    tip_child_min: float
    tip_child_max: float
    edge_ratio_min: float
    edge_ratio_max: float
    edge_ratio_worst: Tuple[int, int]
    edge_ratio_best: Tuple[int, int]
    eh_proxy_max: float


def carrier(
    name: str,
    anchors: Sequence[str],
    scores: Sequence[int],
    preserves: str,
    destroys: str,
    next_hook: str,
) -> Carrier:
    assert len(scores) == len(AXES)
    return Carrier(
        name=name,
        anchors=tuple(anchors),
        axes=dict(zip(AXES, scores)),
        preserves=preserves,
        destroys=destroys,
        next_hook=next_hook,
    )


CARRIERS: List[Carrier] = [
    carrier(
        "uniform_multifar_Rprime_floor",
        ("HYP-3121", "HYP-2968", "HYP-2861", "HYP-2856", "HYP-2871"),
        (9, 10, 7, 10, 4, 8, 7, 8),
        "the exact covering predicate R-safe intersect Q-safe is nonempty",
        "the local reason a joint sector is positive unless sidecars are kept",
        "turn every r=2..6 covering packet into an Rprime floor row with named exception set",
    ),
    carrier(
        "edge_witness_two_ended_RQ_packet",
        ("HYP-3124", "HYP-3050", "HYP-3054", "HYP-3106"),
        (9, 9, 10, 8, 5, 10, 8, 9),
        "tail packet, tip packet, four-sector deck, and both deletion children",
        "global distribution estimates if the packet is kept purely local",
        "attach HYP-3124 edge_witness_certificate to R-safe -> Q-safe packet edges",
    ),
    carrier(
        "gaussian_heat_minorant_smoothing",
        ("HYP-3032", "HYP-3110", "OPEN-Q-104"),
        (7, 8, 5, 8, 6, 7, 7, 8),
        "a smoothed positive lower bound before finite-ruler desmoothing",
        "sharp endpoint ownership and exact forbidden-center collisions",
        "build a heat-kernel/Selberg minorant whose loss is paid by finite-ruler thresholds",
    ),
    carrier(
        "EH_level_distribution_proxy",
        ("HYP-3024", "HYP-3032", "HYP-3084"),
        (7, 8, 4, 9, 4, 7, 6, 8),
        "average residue/frequency cancellation outside named major arcs",
        "individual packet witnesses; Elliott-Halberstam is only an analogy here",
        "measure an L2/BV level-of-distribution proxy over edge-sector residue bins",
    ),
    carrier(
        "asano_lee_yang_edge_contraction",
        ("HYP-3122", "HYP-3112", "HYP-3108", "THM-082"),
        (8, 7, 8, 6, 10, 8, 6, 7),
        "zero-free diagonal contraction of a two-ended edge partition packet",
        "measure-level positivity if zero-free status is not converted back to mass",
        "test whether tail/tip sector polynomials admit Lee-Yang/Asano contraction sidecars",
    ),
    carrier(
        "phi4_quartic_cumulant_stabilizer",
        ("HYP-3122", "HYP-3103", "HYP-3113"),
        (7, 6, 5, 6, 8, 7, 7, 8),
        "quartic cumulant sign controlling the hard binding-row dip",
        "the exact R/Q witness coordinate unless coupled to edge sectors",
        "add kappa4/kappa2^2 and Lee-Yang root-margin columns to edge-floor packets",
    ),
    carrier(
        "normal_fan_cech_finite_ruler",
        ("HYP-3123", "HYP-3101", "THM-573", "THM-565"),
        (8, 7, 4, 6, 4, 9, 8, 8),
        "component/barcode data in normalized finite-ruler coordinates",
        "tail/tip orientation unless the chiral/edge sidecar is joined",
        "use Cech component payload as the desmoothing receiver for the Gaussian route",
    ),
    carrier(
        "chiral_cross_sector_guard",
        ("HYP-3123", "HYP-3124", "HYP-3106"),
        (8, 6, 8, 5, 5, 8, 7, 9),
        "mirror/converse/rootless orientation around edge sectors",
        "quantitative lower bounds if not coupled to Rprime",
        "require cross_sector_orientation_word on every R-safe -> Q-safe witness edge",
    ),
    carrier(
        "H7_state_lift_zero_kernel",
        ("THM-572", "THM-029", "THM-343"),
        (9, 4, 5, 3, 4, 7, 8, 10),
        "terminal contradiction route for zero-open-mass residual kernels",
        "positive open decorrelation floor information",
        "route packets with no positive smoothed mass to the H=7 state-lift constructor",
    ),
    carrier(
        "raw_Bonferroni_or_scalar_p0",
        ("HYP-3121", "HYP-3116"),
        (2, 2, 1, 2, 1, 1, 2, 1),
        "a quick negative-control scalar check",
        "the quasi-independence, endpoint owner, and tail/tip witness payload",
        "keep only as a failure alarm; it is known to fail on few-apex rows",
    ),
]


TIE_PATH = [carrier.name for carrier in CARRIERS]


def mod_distance(x: float) -> float:
    y = x % 1.0
    return y if y <= 0.5 else 1.0 - y


def safe_R_speed(speed: int, u: float) -> bool:
    return mod_distance(speed * u / 14.0) >= THRESHOLD


def safe_Q_speed(speed: int, u: float) -> bool:
    return mod_distance(speed * u) >= THRESHOLD


def split_row(speeds: Sequence[int]) -> Tuple[List[int], List[int]]:
    R = [s for s in speeds if s % 14 != 0]
    Q = [s // 14 for s in speeds if s % 14 == 0]
    return R, Q


def ratio_for_parts(R: Sequence[int], Q: Sequence[int]) -> Tuple[float, float, float, float, float]:
    r_count = q_count = joint_count = 0
    for j in range(GRID_STEPS):
        u = 14.0 * (j + 0.5) / GRID_STEPS
        r_ok = all(safe_R_speed(s, u) for s in R)
        q_ok = all(safe_Q_speed(q, u) for q in Q)
        r_count += int(r_ok)
        q_count += int(q_ok)
        joint_count += int(r_ok and q_ok)
    R_safe = r_count / GRID_STEPS
    Q_safe = q_count / GRID_STEPS
    joint = joint_count / GRID_STEPS
    Rprime = joint / (R_safe * Q_safe) if R_safe and Q_safe else 0.0
    return R_safe, Q_safe, joint, R_safe + Q_safe - 1.0, Rprime


def edge_ratio(R_speed: int, Q_speed: int) -> Tuple[float, Tuple[int, int, int, int]]:
    tail_count = tip_count = both_count = 0
    deck = [0, 0, 0, 0]
    for j in range(GRID_STEPS):
        u = 14.0 * (j + 0.5) / GRID_STEPS
        tail_ok = safe_R_speed(R_speed, u)
        tip_ok = safe_Q_speed(Q_speed, u)
        tail_count += int(tail_ok)
        tip_count += int(tip_ok)
        both_count += int(tail_ok and tip_ok)
        if tail_ok and tip_ok:
            deck[0] += 1
        elif tail_ok and not tip_ok:
            deck[1] += 1
        elif (not tail_ok) and tip_ok:
            deck[2] += 1
        else:
            deck[3] += 1
    tail = tail_count / GRID_STEPS
    tip = tip_count / GRID_STEPS
    both = both_count / GRID_STEPS
    ratio = both / (tail * tip) if tail and tip else 0.0
    return ratio, tuple(deck)


def eh_distribution_proxy(R: Sequence[int], Q: Sequence[int], bins: int) -> float:
    totals = [0] * bins
    hits = [0] * bins
    for j in range(GRID_STEPS):
        u = 14.0 * (j + 0.5) / GRID_STEPS
        b = j % bins
        ok = all(safe_R_speed(s, u) for s in R) and all(safe_Q_speed(q, u) for q in Q)
        totals[b] += 1
        hits[b] += int(ok)
    global_density = sum(hits) / sum(totals)
    mean_square = 0.0
    for h, total in zip(hits, totals):
        local = h / total if total else 0.0
        mean_square += (local - global_density) ** 2
    return sqrt(mean_square / bins)


def audit_row(name: str, speeds: Sequence[int]) -> RowAudit:
    R, Q = split_row(speeds)
    R_safe, Q_safe, joint, bonferroni_margin, Rprime = ratio_for_parts(R, Q)

    tail_children = []
    for s in R:
        child_R = [x for x in R if x != s]
        tail_children.append(ratio_for_parts(child_R, Q)[4])

    tip_children = []
    for q in Q:
        child_Q = list(Q)
        child_Q.remove(q)
        tip_children.append(ratio_for_parts(R, child_Q)[4])
    if not tip_children:
        tip_children = [1.0]

    ratios: List[Tuple[float, int, int]] = []
    for s in R:
        for q in Q:
            ratio, _deck = edge_ratio(s, q)
            ratios.append((ratio, s, q))
    if ratios:
        worst = min(ratios)
        best = max(ratios)
        edge_min = worst[0]
        edge_max = best[0]
        worst_pair = (worst[1], worst[2])
        best_pair = (best[1], best[2])
    else:
        edge_min = edge_max = 1.0
        worst_pair = best_pair = (0, 0)

    eh_proxy = max(eh_distribution_proxy(R, Q, bins) for bins in (14, 28, 98))

    return RowAudit(
        name=name,
        speeds=tuple(speeds),
        r=len(Q),
        R_size=len(R),
        R_safe=R_safe,
        Q_safe=Q_safe,
        joint=joint,
        bonferroni_margin=bonferroni_margin,
        Rprime=Rprime,
        tail_child_min=min(tail_children),
        tail_child_max=max(tail_children),
        tip_child_min=min(tip_children),
        tip_child_max=max(tip_children),
        edge_ratio_min=edge_min,
        edge_ratio_max=edge_max,
        edge_ratio_worst=worst_pair,
        edge_ratio_best=best_pair,
        eh_proxy_max=eh_proxy,
    )


def compare(a: Carrier, b: Carrier) -> int:
    wins_a = wins_b = 0
    for axis in AXES:
        if a.axes[axis] > b.axes[axis]:
            wins_a += 1
        elif b.axes[axis] > a.axes[axis]:
            wins_b += 1
    if wins_a > wins_b:
        return 1
    if wins_b > wins_a:
        return -1
    return 1 if TIE_PATH.index(a.name) < TIE_PATH.index(b.name) else -1


def adjacency(carriers: Sequence[Carrier]) -> Dict[str, List[str]]:
    out = {carrier.name: [] for carrier in carriers}
    for a, b in combinations(carriers, 2):
        if compare(a, b) > 0:
            out[a.name].append(b.name)
        else:
            out[b.name].append(a.name)
    return out


def directed_three_cycles(out: Dict[str, List[str]]) -> List[Tuple[str, str, str]]:
    names = list(out)
    cycles = []
    edge = {(a, b) for a, bs in out.items() for b in bs}
    for a, b, c in combinations(names, 3):
        if (a, b) in edge and (b, c) in edge and (c, a) in edge:
            cycles.append((a, b, c))
        elif (a, c) in edge and (c, b) in edge and (b, a) in edge:
            cycles.append((a, c, b))
    return cycles


def scc_sizes(out: Dict[str, List[str]]) -> List[int]:
    names = list(out)
    reverse = {name: [] for name in names}
    for a, bs in out.items():
        for b in bs:
            reverse[b].append(a)

    seen = set()
    order: List[str] = []

    def dfs(v: str) -> None:
        seen.add(v)
        for w in out[v]:
            if w not in seen:
                dfs(w)
        order.append(v)

    for name in names:
        if name not in seen:
            dfs(name)

    seen.clear()
    sizes = []

    def rdfs(v: str, acc: List[str]) -> None:
        seen.add(v)
        acc.append(v)
        for w in reverse[v]:
            if w not in seen:
                rdfs(w, acc)

    for name in reversed(order):
        if name not in seen:
            acc: List[str] = []
            rdfs(name, acc)
            sizes.append(len(acc))
    return sorted(sizes, reverse=True)


def hamiltonian_path_count(out: Dict[str, List[str]]) -> int:
    names = list(out)
    n = len(names)
    index = {name: i for i, name in enumerate(names)}
    edge = [[False] * n for _ in range(n)]
    for a, bs in out.items():
        for b in bs:
            edge[index[a]][index[b]] = True

    dp = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for mask in range(1 << n):
        for last in range(n):
            count = dp[mask][last]
            if not count:
                continue
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                if edge[last][nxt]:
                    dp[mask | (1 << nxt)][nxt] += count
    return sum(dp[(1 << n) - 1])


def selected_path(out: Dict[str, List[str]]) -> List[str]:
    remaining = set(out)
    path = []
    while remaining:
        best = max(
            remaining,
            key=lambda name: (
                sum(1 for b in out[name] if b in remaining),
                -TIE_PATH.index(name),
            ),
        )
        path.append(best)
        remaining.remove(best)
    return path


def tournament_fingerprint(carriers: Sequence[Carrier]) -> Dict[str, object]:
    out = adjacency(carriers)
    scores = {name: len(bs) for name, bs in out.items()}
    return {
        "score_hist": dict(sorted(Counter(scores.values()).items())),
        "directed_3cycles": len(directed_three_cycles(out)),
        "scc_sizes": scc_sizes(out),
        "hamiltonian_path_count": hamiltonian_path_count(out),
        "selected_hamiltonian_path": selected_path(out),
        "scores": scores,
    }


def fmt(x: float) -> str:
    return f"{x:.6f}"


def main() -> None:
    rows = [
        ("single_far_control_r1", tuple(list(range(1, 12)) + [13, 84])),
        ("two_far_multifar_r2", tuple(list(range(1, 11)) + [13, 84, 154])),
        ("three_far_probe_r3", tuple(list(range(1, 10)) + [13, 84, 154, 224])),
        ("four_far_probe_r4", tuple(list(range(1, 9)) + [13, 84, 154, 224, 294])),
    ]
    audits = [audit_row(name, speeds) for name, speeds in rows]
    fp = tournament_fingerprint(CARRIERS)

    print("HYP-3125 / codex-2026-06-27-S271")
    print("Multi-far edge-witness Rprime floor scout; executable synthesis, not a proof.")
    print()
    print("1. ASSUMPTION CHALLENGE")
    print("candidate vertex sets considered=runners,gaps,sections,boundaries,wall-crossing events,residues,cover arcs,Fourier modes,Gaussian kernels,Asano variables,edge witnesses,tail packets,tip packets,proof obligations")
    print("chosen vertices=proof carriers / repair operators; edge diagnostics use R-safe -> Q-safe witness edges as measured packets")
    print("preserved predicate=existence of a time with all original speeds at distance >= 1/14 after the u=14t lift")
    print("destroyed by scalarization=tail/tip endpoint owner, four-sector deck, deletion child, residue exception set, zero-free status, finite-ruler desmoothing debt")
    print()

    print("2. GRID AUDIT OF COVERING ROWS")
    print(f"grid_steps={GRID_STEPS}; u in [0,14); threshold=1/14; midpoint rule")
    print("row r |R| R_safe Q_safe joint Bonferroni Rprime tail_child_range tip_child_range edge_ratio_range worst_edge best_edge EH_proxy")
    for a in audits:
        print(
            f"{a.name} {a.r} {a.R_size} "
            f"{fmt(a.R_safe)} {fmt(a.Q_safe)} {fmt(a.joint)} "
            f"{fmt(a.bonferroni_margin)} {fmt(a.Rprime)} "
            f"[{fmt(a.tail_child_min)},{fmt(a.tail_child_max)}] "
            f"[{fmt(a.tip_child_min)},{fmt(a.tip_child_max)}] "
            f"[{fmt(a.edge_ratio_min)},{fmt(a.edge_ratio_max)}] "
            f"{a.edge_ratio_worst} {a.edge_ratio_best} {fmt(a.eh_proxy_max)}"
        )
    print()
    print("Readout: Bonferroni is negative on the named r=1/r=2 rows, but the diagonal ratio Rprime remains positive.  The tail-deletion child usually improves the floor, while the tip-deletion child exposes the multiple-of-14 side as the sharper recursion.  This supports treating Rprime as an edge witness with both endpoint children, not as a one-sided safe-mass scalar.")
    print()

    print("3. PROPOSED EDGE-FLOOR PACKET")
    print("edge_floor_packet = (")
    for field in [
        "edge_witness_id",
        "tail_R_safe_packet",
        "tip_Q_safe_packet",
        "outside_four_sector_deck",
        "edge_tail_tip_sector_word",
        "tail_deletion_child_Rprime",
        "tip_deletion_child_Rprime",
        "Rprime_lower_bound_candidate",
        "multifar_r_bucket",
        "EH_level_distribution_proxy",
        "major_arc_residue_exception_set",
        "gaussian_heat_kernel_width",
        "finite_ruler_desmoothing_threshold",
        "asano_contraction_status",
        "lee_yang_zero_free_after_contraction",
        "phi4_kappa4_sign",
        "normal_fan_chamber_id",
        "chiral_guard_word",
        "terminal_exit_or_named_debt",
    ]:
        print(f"  {field},")
    print(")")
    print()

    print("4. TOURNAMENT ANALYSIS OVER REPAIR OPERATORS")
    print("vertices=proof carriers / repair operators, not runners, primes, Gaussian samples, or scalar ratios")
    print("pairwise_observable=majority over axes " + ",".join(AXES))
    print("switch=A->B when A wins more axes; ties use " + " -> ".join(TIE_PATH))
    print(f"score_hist={fp['score_hist']}")
    print(f"directed_3cycles={fp['directed_3cycles']}")
    print(f"scc_sizes={fp['scc_sizes']}")
    print(f"hamiltonian_path_count={fp['hamiltonian_path_count']}")
    print("selected_hamiltonian_path=" + " -> ".join(fp["selected_hamiltonian_path"]))
    print()
    print("Top carrier hooks:")
    for name in fp["selected_hamiltonian_path"][:6]:
        carrier_obj = next(c for c in CARRIERS if c.name == name)
        print(f"- {name}: preserves={carrier_obj.preserves}; next={carrier_obj.next_hook}")
    print()

    print("5. SYNTHESIS")
    print("EH is used only as a level-of-distribution metaphor: the target is an LRC-specific BV/L2 average over edge-sector residue bins, not a prime-distribution theorem.  Gaussian smoothing is the candidate positivity engine before finite-ruler desmoothing.  Asano contraction is the candidate zero-free legality check for contracting the tail and tip variables of an edge-sector partition packet.  The working proof target is a dichotomy: either the Gaussian/EH edge packet gives a uniform positive Rprime floor after named major-arc exceptions, or the exception carries an Asano/Lee-Yang/phi4/normal-fan sidecar that routes to finite-address discharge or H=7 state-lift debt.")


if __name__ == "__main__":
    main()
