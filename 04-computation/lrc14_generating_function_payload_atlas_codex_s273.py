#!/usr/bin/env python3
"""HYP-3137 scout: generating-function payload atlas for the LRC14 endgame.

This is a synthesis and hypothesis generator, not an LRC proof.  It compares
generating-function carriers by the proof coordinates they preserve before a
scalar evaluation destroys information.

Tournament Analysis declaration:
  vertices: generating-function proof carriers, not runners, roots, arcs, or
            scalar values;
  pairwise observable: majority comparison over retained proof-payload axes;
  switch/gauge: orient A -> B when A wins more axes, with a fixed tie path;
  tie path: signed SPEC, edge-recursion, A000568 quotient, PGF roots,
            resolvent middle layer, Walsh/OCF EGFs, repair Hilbert series,
            hard-core partition, IO determinant, modular tails, raw scalar.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from itertools import combinations
from typing import Dict, Iterable, List, Sequence, Tuple


AXES: Tuple[str, ...] = (
    "lrc_predicate_retention",
    "coefficient_layer_payload",
    "root_zero_locus",
    "signed_tail_certificate",
    "quotient_legality",
    "edge_recursion_gluing",
    "finite_constant_chase",
    "closed_form_telescoping",
    "scalar_collapse_guard",
    "proof_exit_readiness",
)

OBLIGATIONS: Tuple[str, ...] = (
    "Q_apex_floor",
    "R_safe_positive",
    "Rprime_signed_low",
    "SPEC_parseval_tail",
    "bounded_core_k8_phi4",
    "far_push_zero_motion",
    "edge_tail_tip_children",
    "global_consistency_quotient",
    "finite_address_observer_gluing",
    "closed_form_uniform_constants",
    "coefficient_middle_layer",
    "destroyed_coordinate_guard",
)


@dataclass(frozen=True)
class Carrier:
    name: str
    gf_type: str
    anchors: Tuple[str, ...]
    axes: Dict[str, int]
    covers: Tuple[str, ...]
    preserves: str
    destroys: str
    proposed_signal: str
    next_hook: str
    wild_guess: str


def carrier(
    name: str,
    gf_type: str,
    anchors: Sequence[str],
    scores: Sequence[int],
    covers: Sequence[str],
    preserves: str,
    destroys: str,
    proposed_signal: str,
    next_hook: str,
    wild_guess: str,
) -> Carrier:
    if len(scores) != len(AXES):
        raise ValueError(f"{name}: expected {len(AXES)} scores")
    unknown = set(covers) - set(OBLIGATIONS)
    if unknown:
        raise ValueError(f"{name}: unknown obligations {sorted(unknown)}")
    return Carrier(
        name=name,
        gf_type=gf_type,
        anchors=tuple(anchors),
        axes=dict(zip(AXES, scores)),
        covers=tuple(covers),
        preserves=preserves,
        destroys=destroys,
        proposed_signal=proposed_signal,
        next_hook=next_hook,
        wild_guess=wild_guess,
    )


CARRIERS: List[Carrier] = [
    carrier(
        "signed_SPEC_resonance_series",
        "Fourier / singular-series GF",
        ("HYP-3129", "HYP-3136", "HYP-3135", "LTI-260"),
        (10, 8, 4, 10, 7, 7, 10, 8, 10, 10),
        (
            "Rprime_signed_low",
            "SPEC_parseval_tail",
            "closed_form_uniform_constants",
            "coefficient_middle_layer",
            "destroyed_coordinate_guard",
        ),
        "exact low resonance on 14Z plus a Parseval-controlled signed tail",
        "edge-local child ownership if it is reported as one Fourier scalar",
        "SPEC_support_sieve = (resonance_lattice, low_modes, L2_tail, sign_budget)",
        "turn HYP-3129 representative constants into a closed-form finite packet theorem",
        "The LRC analogue of the Walsh support theorem is not PGF roots, but the support law that only the right signed modes survive.",
    ),
    carrier(
        "edge_witness_recursion_OGF",
        "ordinary GF of edge witness decks by recursion depth",
        ("HYP-3124", "HYP-3125", "HYP-3136", "ce3529219"),
        (9, 7, 2, 5, 8, 10, 8, 6, 10, 8),
        (
            "edge_tail_tip_children",
            "finite_address_observer_gluing",
            "destroyed_coordinate_guard",
            "R_safe_positive",
        ),
        "tail/tip deletion children, recursive child decks, and orientation words",
        "global economy of the quotient unless joined to A000568 or observer gluing",
        "edge_recursion_depth_PGF = sum_d (#safe edge decks at depth d) u^d",
        "attach depth-1/depth-2 edge packets to every finite SPEC row before quotienting",
        "The ear-decomposition theorems hint that the missing witness is a typed ear word for proof packets, not a new scalar invariant.",
    ),
    carrier(
        "A000568_cycle_index_quotient",
        "unlabeled tournament OGF / cycle-index quotient",
        ("HYP-3133", "HYP-3134", "HYP-3135", "LTI-262"),
        (8, 5, 2, 4, 10, 10, 8, 6, 9, 8),
        (
            "global_consistency_quotient",
            "edge_tail_tip_children",
            "finite_address_observer_gluing",
            "destroyed_coordinate_guard",
        ),
        "the middle quotient between raw sector decks and paired child packets",
        "endpoint roles unless the quotient map is marked",
        "global_consistency_class = orbit(edge_packet)/legal_gluing",
        "stratify HYP-3129 finite rows by sector -> A000568 shadow -> paired children",
        "A000568 is useful because it is neither local nor fully resolved: it is exactly a controlled forgetting tax.",
    ),
    carrier(
        "miss_count_PGF_root_locus",
        "danger-count PGF and Lee-Yang root locus",
        ("HYP-3128", "HYP-3131", "HYP-3132", "HYP-3136"),
        (8, 6, 10, 3, 5, 5, 7, 5, 8, 7),
        (
            "Q_apex_floor",
            "bounded_core_k8_phi4",
            "far_push_zero_motion",
            "R_safe_positive",
        ),
        "the whole miss-count distribution and zero motion under far additions",
        "signed coupling and tail constants if reduced to P(M=0) or nearest root",
        "PGF_root_trajectory = sorted roots of E[z^M] under far-push and tail deletion",
        "prove the bounded-core implication rho>1 => Rprime>=c or name the gap",
        "The whole PGF curve matters, but the load-bearing coordinate may be its zero-motion derivative under endpoint deletion.",
    ),
    carrier(
        "resolvent_elementary_symmetric_GF",
        "De Moivre / resolvent elementary-symmetric GF",
        ("HYP-3138", "HYP-3132", "HYP-3135", "THM-577", "HYP-3122"),
        (8, 10, 8, 4, 6, 5, 9, 9, 9, 8),
        (
            "bounded_core_k8_phi4",
            "coefficient_middle_layer",
            "closed_form_uniform_constants",
            "far_push_zero_motion",
        ),
        "pair/triple middle layers such as -e2=120 and e3=320",
        "edge quotient legality unless it is joined to the packet ledger",
        "middle_layer_vector = (e2, e3, kappa2, kappa4, phi4_sign, k8_even_fold)",
        "use HYP-3138's fold-adjoint lookup to turn the k=8 biquadratic fold into a coefficient inequality",
        "The proof may close by bounding one reflected u^2 coefficient, with HYP-3138 supplying the finite odd-coordinate resurrection table.",
    ),
    carrier(
        "walsh_ocf_support_EGF",
        "Walsh support EGF / OCF factorization",
        ("THM-077", "THM-076", "THM-059", "05-knowledge/variables/W-polynomial.md"),
        (6, 9, 3, 5, 5, 6, 6, 10, 10, 7),
        (
            "coefficient_middle_layer",
            "closed_form_uniform_constants",
            "destroyed_coordinate_guard",
        ),
        "which coefficient supports survive: even path components and telescoping blocks",
        "direct LRC predicate unless translated into resonance supports",
        "support_survival_word = (allowed supports, cancellation involution, telescoper)",
        "copy the THM-076 constant-term/hockey-stick pattern into the HYP-3129 mode sum",
        "Savitch-style reachability may appear as repeated halving of proof obligations over surviving coefficient supports.",
    ),
    carrier(
        "trivariate_master_GF",
        "W/F/master polynomial and Eulerian deformation GF",
        ("THM-059", "THM-114", "THM-293", "05-knowledge/variables/F-polynomial.md"),
        (5, 9, 5, 4, 4, 5, 6, 9, 9, 6),
        (
            "coefficient_middle_layer",
            "closed_form_uniform_constants",
            "destroyed_coordinate_guard",
        ),
        "a two-axis deformation where evaluation is visibly lossy",
        "LRC signs and edge gluing unless a sidecar restores them",
        "evaluation_loss_vector = (axis, derivative, killed coefficients)",
        "use as the warning template for P(M=0), Rprime, and raw A000568 evaluation",
        "The right LRC GF may be a two-variable deformation: one axis for danger count, one for edge/observer ownership.",
    ),
    carrier(
        "repair_cover_Hilbert_series",
        "Hilbert/Poincare series of coordinate-resurrection covers",
        ("HYP-3118", "HYP-3116", "HYP-3120", "LTI-254"),
        (7, 7, 3, 3, 8, 9, 7, 6, 10, 7),
        (
            "finite_address_observer_gluing",
            "destroyed_coordinate_guard",
            "global_consistency_quotient",
        ),
        "minimal sidecar covers, proof obligations, and quotient sections",
        "analytic floor constants unless joined to SPEC/PGF carriers",
        "repair_cover_H(q) = sum_c (#minimal repair covers of size c) q^c",
        "make every scalar quotient report the minimal sidecar that resurrects its lost coordinate",
        "Bravais-lattice language belongs here: quotient cells need a lattice address, not just a shape name.",
    ),
    carrier(
        "hard_core_independence_partition",
        "hard-core / polymer partition function",
        ("THM-209", "THM-114", "THM-341", "HYP-3128"),
        (5, 8, 6, 5, 4, 5, 5, 8, 8, 5),
        (
            "Q_apex_floor",
            "coefficient_middle_layer",
            "far_push_zero_motion",
        ),
        "power-of-two independent-cycle hierarchy and polymer cluster expansions",
        "edge endpoints and LRC resonance lattice if used generically",
        "polymer_cluster_signature = (activity, incompatible pairs, Mayer blocks)",
        "test whether edge-floor rows have a small polymer cluster expansion after quotienting",
        "The phi4 density exp(-lambda*S^4-b*S^2)dS suggests the hard node is a cumulant-stabilized polymer, not a bare PGF root.",
    ),
    carrier(
        "IO_walk_determinant_GF",
        "Irving-Omar determinant/permanent walk GF",
        ("Irving-Omar arXiv:2412.10572", "T118", "INV-046", "THM-467"),
        (4, 7, 7, 4, 3, 5, 5, 7, 5, 4),
        (
            "coefficient_middle_layer",
            "closed_form_uniform_constants",
        ),
        "det/per and odd-trace log-derivative data from walk/cycle covers",
        "endpoint conditioning and LRC predicate without extra variables",
        "odd_trace_logderivative = coefficients of tr(A^(2k+1))/(2k+1)",
        "use only as a matrix-algebra analogy unless endpoint-conditioned variables are added",
        "If it becomes useful, it will be through a noncommutative endpoint-refined version, not the commutative determinant alone.",
    ),
    carrier(
        "q_pochhammer_modular_tail_GF",
        "q-series/product tail guardrail",
        ("LTT-125", "q-series principal-part ledgers", "modular cusp guardrails"),
        (3, 5, 4, 7, 2, 2, 4, 6, 6, 3),
        (
            "SPEC_parseval_tail",
            "closed_form_uniform_constants",
        ),
        "finite principal parts and named polar exits for infinite products",
        "local LRC packet data if imported as modular decoration",
        "principal_part_tail_certificate = (finite part, pole ledger, tail bound)",
        "keep as a guardrail for any product-form SPEC or theta tail expression",
        "Useful only if the finite constant chase turns into a product tail with a named principal part.",
    ),
    carrier(
        "raw_scalar_evaluation",
        "negative-control scalar",
        ("negative-control",),
        (2, 1, 1, 0, 0, 0, 2, 0, 0, 1),
        tuple(),
        "the visible number: P(M=0), Rprime, A000568(n), root radius, or a residual",
        "almost every proof coordinate needed by HYP-3136",
        "scalar_shadow_alarm = scalar_value_without_payload",
        "never use as a proof carrier unless a sidecar restores the lost coordinates",
        "The single number is the smoke, not the engine.",
    ),
]

TIE_PATH: Tuple[str, ...] = tuple(car.name for car in CARRIERS)


def total_score(car: Carrier) -> int:
    return sum(car.axes.values())


def compare(a: Carrier, b: Carrier) -> str:
    a_votes = 0
    b_votes = 0
    for axis in AXES:
        if a.axes[axis] > b.axes[axis]:
            a_votes += 1
        elif b.axes[axis] > a.axes[axis]:
            b_votes += 1
    if a_votes > b_votes:
        return a.name
    if b_votes > a_votes:
        return b.name
    return a.name if TIE_PATH.index(a.name) < TIE_PATH.index(b.name) else b.name


def tournament_edges(vertices: Sequence[Carrier]) -> Dict[str, set[str]]:
    adj: Dict[str, set[str]] = {v.name: set() for v in vertices}
    for a, b in combinations(vertices, 2):
        winner = compare(a, b)
        loser = b.name if winner == a.name else a.name
        adj[winner].add(loser)
    return adj


def tarjan_scc(adj: Dict[str, set[str]]) -> List[List[str]]:
    index = 0
    stack: List[str] = []
    on_stack: set[str] = set()
    indices: Dict[str, int] = {}
    lowlink: Dict[str, int] = {}
    out: List[List[str]] = []

    def strongconnect(v: str) -> None:
        nonlocal index
        indices[v] = index
        lowlink[v] = index
        index += 1
        stack.append(v)
        on_stack.add(v)

        for w in adj[v]:
            if w not in indices:
                strongconnect(w)
                lowlink[v] = min(lowlink[v], lowlink[w])
            elif w in on_stack:
                lowlink[v] = min(lowlink[v], indices[w])

        if lowlink[v] == indices[v]:
            comp: List[str] = []
            while True:
                w = stack.pop()
                on_stack.remove(w)
                comp.append(w)
                if w == v:
                    break
            out.append(sorted(comp, key=TIE_PATH.index))

    for vertex in adj:
        if vertex not in indices:
            strongconnect(vertex)
    return out


def directed_3cycles(adj: Dict[str, set[str]]) -> List[Tuple[str, str, str]]:
    cycles: List[Tuple[str, str, str]] = []
    for a, b, c in combinations(adj, 3):
        # Canonical orientation cycle a->b->c->a or reverse.
        if b in adj[a] and c in adj[b] and a in adj[c]:
            cycles.append((a, b, c))
        elif c in adj[a] and b in adj[c] and a in adj[b]:
            cycles.append((a, c, b))
    return cycles


def hamiltonian_path_count(adj: Dict[str, set[str]]) -> int:
    vertices = list(adj)
    n = len(vertices)
    idx = {v: i for i, v in enumerate(vertices)}
    dp: Dict[Tuple[int, int], int] = {}
    for v in vertices:
        dp[(1 << idx[v], idx[v])] = 1
    for mask in range(1 << n):
        for last in range(n):
            count = dp.get((mask, last), 0)
            if not count:
                continue
            last_name = vertices[last]
            for nxt_name in adj[last_name]:
                nxt = idx[nxt_name]
                if not (mask & (1 << nxt)):
                    dp[(mask | (1 << nxt), nxt)] = dp.get((mask | (1 << nxt), nxt), 0) + count
    full = (1 << n) - 1
    return sum(dp.get((full, last), 0) for last in range(n))


def hamiltonian_path(adj: Dict[str, set[str]], vertices: Sequence[Carrier]) -> List[str]:
    """Construct a genuine directed Hamiltonian path in any tournament."""
    order = sorted(vertices, key=lambda c: (-total_score(c), c.name))
    path: List[str] = []
    for car in order:
        inserted = False
        for pos, existing in enumerate(path):
            if existing in adj[car.name]:
                path.insert(pos, car.name)
                inserted = True
                break
        if not inserted:
            path.append(car.name)
    return path


def cover_score(names: Iterable[str]) -> Tuple[int, int]:
    selected = [car for car in CARRIERS if car.name in set(names)]
    covered = set().union(*(set(car.covers) for car in selected)) if selected else set()
    return len(covered), sum(total_score(car) for car in selected)


def minimal_covers(limit: int = 8) -> List[Tuple[Tuple[str, ...], int]]:
    target = set(OBLIGATIONS)
    covers: List[Tuple[Tuple[str, ...], int]] = []
    names = [car.name for car in CARRIERS if car.covers]
    for size in range(1, min(6, len(names)) + 1):
        for combo in combinations(names, size):
            got = set()
            for name in combo:
                got.update(next(car.covers for car in CARRIERS if car.name == name))
            if got >= target:
                covers.append((combo, sum(total_score(next(car for car in CARRIERS if car.name == name)) for name in combo)))
        if covers:
            break
    covers.sort(key=lambda item: (-item[1], item[0]))
    return covers[:limit]


def print_map_1() -> None:
    print("MAP 1: GENERATING-FUNCTION CARRIERS BY RETAINED PAYLOAD")
    print("carrier | total | covers | preserves | destroys")
    for car in sorted(CARRIERS, key=lambda c: (-total_score(c), c.name)):
        covers = ", ".join(car.covers) if car.covers else "NONE"
        print(f"- {car.name} [{car.gf_type}] total={total_score(car)}")
        print(f"  anchors: {', '.join(car.anchors)}")
        print(f"  covers: {covers}")
        print(f"  preserves: {car.preserves}")
        print(f"  destroys if scalarized: {car.destroys}")
        print(f"  signal: {car.proposed_signal}")
        print(f"  next: {car.next_hook}")


def print_map_2() -> None:
    print()
    print("MAP 2: LRC14 REMAINING OBLIGATIONS AND BEST CARRIER SUPPORT")
    by_obligation: Dict[str, List[Carrier]] = {ob: [] for ob in OBLIGATIONS}
    for car in CARRIERS:
        for ob in car.covers:
            by_obligation[ob].append(car)

    for ob in OBLIGATIONS:
        carriers = sorted(by_obligation[ob], key=lambda c: (-total_score(c), c.name))
        names = ", ".join(car.name for car in carriers) if carriers else "NONE"
        print(f"- {ob}: {names}")

    print()
    print("SMALLEST FULL-COVER PACKETS")
    covers = minimal_covers()
    for combo, score in covers:
        got = set()
        for name in combo:
            got.update(next(car.covers for car in CARRIERS if car.name == name))
        print(f"- size={len(combo)} score={score}: {' + '.join(combo)}")
        print(f"  covers {len(got)}/{len(OBLIGATIONS)} obligations")


def print_tournament() -> None:
    print()
    print("TOURNAMENT ANALYSIS OVER GF CARRIERS")
    adj = tournament_edges(CARRIERS)
    score_hist = Counter(len(out) for out in adj.values())
    print(f"vertices={len(CARRIERS)} edges={sum(len(v) for v in adj.values())}")
    print(f"score histogram={dict(sorted(score_hist.items()))}")
    print(f"directed 3-cycles={len(directed_3cycles(adj))}")
    print(f"SCCs={tarjan_scc(adj)}")
    print(f"Hamiltonian path count={hamiltonian_path_count(adj)}")
    path = hamiltonian_path(adj, CARRIERS)
    print("selected Hamiltonian path:")
    print("  " + " -> ".join(path))

    print()
    print("PAIRWISE EDGE FLIPS AGAINST RAW SCALAR")
    raw = next(car for car in CARRIERS if car.name == "raw_scalar_evaluation")
    for car in sorted(CARRIERS, key=lambda c: (-total_score(c), c.name)):
        if car.name == raw.name:
            continue
        winner = compare(car, raw)
        print(f"- {car.name} {'beats' if winner == car.name else 'loses to'} raw_scalar_evaluation")


def print_hypotheses() -> None:
    print()
    print("NEW SIGNALS TO MEASURE")
    rows = [
        (
            "SPEC_support_sieve",
            "Walsh support theorem -> HYP-3129 resonance support: record which signed modes survive, not only their sum.",
        ),
        (
            "edge_recursion_depth_PGF",
            "Ear-decomposition analogy: count proof-edge decks by recursive depth and mark tail/tip ownership.",
        ),
        (
            "global_consistency_class",
            "A000568 wedge as quotient legality: local sectors < global quotient < paired children.",
        ),
        (
            "PGF_root_trajectory_derivative",
            "Measure how zeros move under far addition and endpoint deletion, not just nearest-root radius.",
        ),
        (
            "middle_layer_vector",
            "Track pair/triple elementary symmetric layers, cumulants, and phi4 sign as the bounded-core payload.",
        ),
        (
            "k8_reflection_fold_adjoint",
            "Use HYP-3138's even fold as a finite lookup that preserves gK8 while resurrecting odd leakage or naming debt.",
        ),
        (
            "repair_cover_H(q)",
            "Hilbert series for minimal sidecars needed to resurrect coordinates destroyed by each quotient.",
        ),
        (
            "Bravais_resonance_cell",
            "Record the lattice cell/chamber of low SPEC modes; quotient shapes without a cell address are unsafe.",
        ),
        (
            "Savitch_packet_depth",
            "Try repeated-squaring style proof reachability: can edge packets be joined in logarithmic depth?",
        ),
    ]
    for name, reading in rows:
        print(f"- {name}: {reading}")

    print()
    print("BOLD CURRENT HYPOTHESES")
    for car in sorted(CARRIERS, key=lambda c: (-total_score(c), c.name))[:8]:
        print(f"- {car.name}: {car.wild_guess}")

    print()
    print("WORKING REFRAIN")
    print(
        "The single evaluation is never the whole GF curve.  The proof packet wants "
        "a small, typed payload: support law + signed tail + edge quotient + "
        "middle coefficient layer + finite observer exit."
    )


def main() -> None:
    print("HYP-3137 / codex-2026-06-27-S273")
    print("LRC14 GENERATING-FUNCTION PAYLOAD ATLAS")
    print("status: executable synthesis scout; not a proof")
    print_map_1()
    print_map_2()
    print_tournament()
    print_hypotheses()


if __name__ == "__main__":
    main()
