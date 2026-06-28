#!/usr/bin/env python3
"""LRC14 proof-angle scout: sheaf exactness plus cusp transfer.

This is a synthesis script, not a proof.  It deliberately avoids making the
next pass another scalar extremality, observability-matrix, or Morse-descent
run.  The two selected proof programs are:

1. a first-obstruction sheaf exactness / holonomy route, and
2. a Farey-cusp renormalization transfer route for the covering residual.

Tournament Analysis is over proof carriers and chart transitions, not runners.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from itertools import combinations, permutations
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parents[1]


WEIGHTS = {
    "boundary_open_status": 9,
    "chart_gluing_exactness": 9,
    "covering_residual_exit": 9,
    "finite_checkability": 7,
    "scale_normality": 7,
    "endpoint_owner": 7,
    "topological_hole": 7,
    "arithmetic_period": 7,
    "state_lift_exit": 7,
    "bulk_floor_margin": 8,
    "residue_magnitude_split": 7,
    "formal_packet_schema": 6,
}


KEY_FILES = [
    "00-navigation/LRC-LENS-MAP.md",
    "00-navigation/LRC-TECHNIQUE-INDEX.md",
    "00-navigation/LRC-TOURNAMENT-TECHNIQUE-INDEX.md",
    "00-navigation/OPEN-QUESTIONS.md",
    "05-knowledge/hypotheses/INDEX.md",
    "05-knowledge/hypotheses/HYP-3253-lrc14-contact-holonomy-curvature-sheaf.md",
    "05-knowledge/hypotheses/HYP-3265-lrc14-equioscillation-contact-graph-case-split.md",
    "05-knowledge/hypotheses/HYP-3300-lrc14-observability-morse-proof-angles.md",
    "07-reflections/lrc14-first-obstruction-cocycle-generation-codex-s259.md",
    "07-reflections/lrc14-boundary-moment-adjunction-missing-picture-codex.md",
    "07-reflections/lrc14-universal-scale-invariance-recursion-ledger-codex-20260628.md",
    "07-reflections/lrc-iterated-logs-are-the-inverse-hyperoperation-tower-s597.md",
    "07-reflections/minkowski-circuit-ising-demoivre-carrier-atlas-codex-s264.md",
    "comms/POKE-COORDINATION.md",
    "comms/POKE-FORUM.md",
]


@dataclass(frozen=True)
class Angle:
    key: str
    name: str
    family: str
    preserves: tuple[str, ...]
    destroys: tuple[str, ...]
    proof_test: str
    cross_links: tuple[str, ...]
    risk: str
    priority_bonus: int = 0


ANGLES: tuple[Angle, ...] = (
    Angle(
        "B01",
        "first_obstruction_sheaf_exactness",
        "Cech/local-system exactness route",
        (
            "chart_gluing_exactness",
            "boundary_open_status",
            "topological_hole",
            "endpoint_owner",
            "finite_checkability",
            "state_lift_exit",
        ),
        ("raw_chart_scalar",),
        (
            "Build the quotient/observer overlap cochain complex; prove every "
            "first obstruction is exact, killed by contact holonomy, lifted to "
            "an endpoint chamber, descended, stopped at AP/GW equality, or "
            "named as state-lift debt."
        ),
        (
            "HYP-3253",
            "HYP-3242",
            "HYP-3234",
            "HYP-3102",
            "Cech cohomology",
            "local systems",
            "Hodge decomposition",
            "LDPC syndrome",
        ),
        "a sheaf label is empty unless each chart overlap gets an explicit payload",
        12,
    ),
    Angle(
        "B02",
        "farey_cusp_renormalization_transfer",
        "modular/Farey scale-transfer route",
        (
            "covering_residual_exit",
            "scale_normality",
            "arithmetic_period",
            "bulk_floor_margin",
            "formal_packet_schema",
            "residue_magnitude_split",
        ),
        ("single_denominator_view",),
        (
            "Turn qdiv>14 covering rows into a scale-normal Farey-cusp "
            "transition: exact-period boundary -> boundary-moment image -> "
            "positive floor, impossible AP/GW kernel, K33/H7 debt, or a named "
            "new zero-open kernel."
        ),
        (
            "HYP-2963",
            "HYP-2969",
            "HYP-3230",
            "HYP-3231",
            "THM-573",
            "Farey tree",
            "modular cusp",
            "renormalization group",
        ),
        "must not turn the cusp tree into a metaphor without exact packet labels",
        11,
    ),
    Angle(
        "B03",
        "boundary_moment_chain_map",
        "exact-period to moment-dual route",
        (
            "covering_residual_exit",
            "arithmetic_period",
            "bulk_floor_margin",
            "finite_checkability",
        ),
        ("endpoint_owner", "topological_hole"),
        (
            "Prove the signed exact-period packet boundary maps to the HYP-2704 "
            "missed-depth/gK8 currency with only AP/GW, K33, H7, or positive "
            "image kernels."
        ),
        (
            "HYP-2954",
            "HYP-2969",
            "HYP-2704",
            "gK8",
            "Ramanujan exact-period packets",
        ),
        "moment positivity is unsafe if the exact-period owner labels are dropped",
        6,
    ),
    Angle(
        "B04",
        "contact_graph_kernel_classifier",
        "unit-contact equality-core route",
        (
            "boundary_open_status",
            "topological_hole",
            "residue_magnitude_split",
            "finite_checkability",
        ),
        ("covering_residual_exit",),
        (
            "Classify surviving six-unit contact graphs into AP/GW equality, "
            "strict off-unit margin, killed-unit covering route, or debt."
        ),
        (
            "HYP-3265",
            "HYP-3257",
            "HYP-3255",
            "THM-523",
            "Chebyshev equioscillation",
        ),
        "unit contacts do not see same-residue height moves or covering floors",
        4,
    ),
    Angle(
        "B05",
        "qsqrt7_residue_to_magnitude_lift",
        "Galois residue/magnitude split route",
        (
            "residue_magnitude_split",
            "boundary_open_status",
            "arithmetic_period",
            "finite_checkability",
        ),
        ("bulk_floor_margin",),
        (
            "Use the C3 orbit over Q(cos 2pi/7) as residue skeleton, then "
            "force a magnitude-level certificate for tight/loose separation."
        ),
        (
            "HYP-3255",
            "HYP-3257",
            "Q(sqrt(-7))",
            "Galois orbit",
            "Hurwitz-7",
        ),
        "residue symmetry organizes the equality core but does not prove the floor",
        4,
    ),
    Angle(
        "B06",
        "schwarz_christoffel_contact_polygon",
        "polygon-angle / accessory-parameter route",
        (
            "boundary_open_status",
            "endpoint_owner",
            "topological_hole",
        ),
        ("arithmetic_period", "bulk_floor_margin"),
        (
            "Read contact times as polygon vertices and destroyed off-unit "
            "geometry as accessory-parameter debt; useful only if the polygon "
            "word is tied to endpoint owners and exact period packets."
        ),
        (
            "Schwarz-Christoffel mapping",
            "HYP-3265",
            "normal-fan Cech bars",
            "contact graph",
        ),
        "beautiful polygon grammar can hide the arithmetic route labels",
        2,
    ),
    Angle(
        "B07",
        "loglog_mertens_scale_entropy",
        "scale-hierarchy entropy route",
        (
            "scale_normality",
            "formal_packet_schema",
            "arithmetic_period",
        ),
        ("endpoint_owner", "topological_hole", "state_lift_exit"),
        (
            "Use Mertens/Hardy-Ramanujan loglog behavior as a complexity bound "
            "on obstruction-prime ledgers, not as a witness inequality."
        ),
        (
            "Mertens",
            "Hardy-Ramanujan",
            "loglog",
            "renormalization depth",
            "HYP-2145",
        ),
        "complexity thinness is not an LRC certificate without packet exits",
        1,
    ),
    Angle(
        "B08",
        "hyperbolic_reciprocal_defect",
        "Fermat-Catalan / orbifold deficit route",
        (
            "formal_packet_schema",
            "scale_normality",
            "residue_magnitude_split",
        ),
        ("boundary_open_status", "endpoint_owner"),
        (
            "Treat sum(1/a_i)<1 as a hyperbolic-defect warning for triple "
            "cusp widths: useful to classify route curvature, not to close "
            "the LRC inequality directly."
        ),
        (
            "Fermat-Catalan",
            "Hurwitz theorem",
            "Markov numbers",
            "hyperbolic orbifolds",
        ),
        "the reciprocal bound becomes numerology unless tied to a packet triple",
        1,
    ),
    Angle(
        "B09",
        "proof_circuit_uniform_gate_test",
        "complexity missing-input route",
        (
            "formal_packet_schema",
            "finite_checkability",
            "chart_gluing_exactness",
        ),
        ("bulk_floor_margin", "arithmetic_period"),
        (
            "Use circuit inputs as a regression test: any proof shortcut that "
            "drops finite address, observer gluing, or exact-period packet "
            "fields cannot compute the residual predicate."
        ),
        (
            "HYP-3116",
            "circuit complexity",
            "proof DAG",
            "missing-input vector",
        ),
        "organizes proof obligations but rarely supplies the inequality itself",
        0,
    ),
    Angle(
        "B10",
        "raw_cross_discipline_analogy",
        "discarded baseline",
        ("formal_packet_schema",),
        (
            "boundary_open_status",
            "endpoint_owner",
            "covering_residual_exit",
            "topological_hole",
            "bulk_floor_margin",
        ),
        (
            "Name external theorems without declaring preserved predicate, "
            "destroyed coordinate, and terminal handoff."
        ),
        ("warning baseline",),
        "this is the failure mode the scout is trying to avoid",
        -20,
    ),
)


EXACTNESS_COLUMNS = [
    "unit_contact_matching",
    "danger_nerve_hole",
    "endpoint_arrangement_cell",
    "zeta7_contact_holonomy",
    "shell_lag_commutator",
    "index_degree_packet",
    "floor_margin_certificate",
    "signed_address_word",
    "first_obstruction_cocycle",
    "repair_cover_rank",
    "Qsqrt7_residue_orbit",
    "magnitude_jump_anchor",
    "unit_nullspace",
    "exact_period_boundary",
    "boundary_moment_image",
    "Farey_parent_interval",
    "cusp_principal_part",
    "renormalization_depth",
    "hensel_unit",
    "K33_kernel_label",
    "state_lift_H7",
    "schwarz_christoffel_turning_word",
    "accessory_parameter_debt",
]


EXACTNESS_ROWS = {
    "unit_contact_to_danger_nerve": {
        "unit_contact_matching",
        "danger_nerve_hole",
        "endpoint_arrangement_cell",
        "zeta7_contact_holonomy",
    },
    "lag_residue_to_shell_magic": {
        "zeta7_contact_holonomy",
        "shell_lag_commutator",
        "endpoint_arrangement_cell",
        "first_obstruction_cocycle",
    },
    "index_degree_to_floor_gap": {
        "index_degree_packet",
        "floor_margin_certificate",
        "danger_nerve_hole",
        "zeta7_contact_holonomy",
    },
    "residue_to_magnitude_split": {
        "Qsqrt7_residue_orbit",
        "magnitude_jump_anchor",
        "unit_nullspace",
        "floor_margin_certificate",
    },
    "chart_change_to_first_obstruction": {
        "signed_address_word",
        "first_obstruction_cocycle",
        "repair_cover_rank",
    },
    "exact_period_to_boundary_moment": {
        "exact_period_boundary",
        "boundary_moment_image",
        "Farey_parent_interval",
        "hensel_unit",
    },
    "covering_to_state_lift_kernel": {
        "exact_period_boundary",
        "boundary_moment_image",
        "K33_kernel_label",
        "state_lift_H7",
    },
    "farey_cusp_to_scale_normal": {
        "Farey_parent_interval",
        "cusp_principal_part",
        "renormalization_depth",
        "hensel_unit",
    },
    "polygon_contact_to_endpoint_owner": {
        "schwarz_christoffel_turning_word",
        "accessory_parameter_debt",
        "unit_contact_matching",
        "endpoint_arrangement_cell",
    },
    "renormalized_kernel_to_named_debt": {
        "renormalization_depth",
        "boundary_moment_image",
        "K33_kernel_label",
        "state_lift_H7",
        "repair_cover_rank",
    },
}


RENORMALIZATION_STATES = {
    "F0_qdiv_witness": {
        "potential": (0, 0, 0, 0),
        "exit": "terminal loose witness with M >= 1/qdiv > 1/14",
    },
    "F1_q14_unit_contact_core": {
        "potential": (1, 0, 0, 0),
        "exit": "AP/GW equality or strict off-unit margin by contact classifier",
    },
    "F2_q14_killed_contact_covering": {
        "potential": (2, 1, 0, 1),
        "exit": "promote to Phi_14d or covering floor sidecar",
    },
    "F3_qgt14_exact_period_cover": {
        "potential": (3, 2, 1, 1),
        "exit": "exact-period boundary feeds boundary-moment map",
    },
    "F4_positive_boundary_moment_image": {
        "potential": (1, 1, 0, 0),
        "exit": "strict covering floor witness",
    },
    "F5_AP_GW_kernel_in_covering": {
        "potential": (2, 1, 1, 0),
        "exit": "forbidden by qdiv>14 unless scale labels were lost",
    },
    "F6_K33_H7_named_kernel": {
        "potential": (2, 1, 1, 1),
        "exit": "route to K33/state-lift theorem debt",
    },
    "F7_unknown_zero_open_kernel": {
        "potential": (4, 3, 2, 2),
        "exit": "the actual remaining target: prove empty or name new sidecar",
    },
}


RENORMALIZATION_EDGES = [
    ("F3_qgt14_exact_period_cover", "F4_positive_boundary_moment_image", "B_D positive image"),
    ("F3_qgt14_exact_period_cover", "F5_AP_GW_kernel_in_covering", "kernel test"),
    ("F3_qgt14_exact_period_cover", "F6_K33_H7_named_kernel", "state-lift exception"),
    ("F3_qgt14_exact_period_cover", "F7_unknown_zero_open_kernel", "unspanned obstruction"),
    ("F2_q14_killed_contact_covering", "F3_qgt14_exact_period_cover", "covering promotion"),
    ("F1_q14_unit_contact_core", "F0_qdiv_witness", "strict off-unit peak"),
    ("F1_q14_unit_contact_core", "F2_q14_killed_contact_covering", "unit contacts killed"),
]


def corpus_hits(angle: Angle) -> int:
    words = {angle.name.replace("_", " ")}
    words.update(angle.cross_links)
    words.update(angle.preserves)
    total = 0
    for rel in KEY_FILES:
        path = ROOT / rel
        if not path.exists():
            continue
        text = path.read_text(errors="ignore").lower()
        for word in words:
            token = word.lower()
            if token:
                total += min(text.count(token), 8)
    return total


def angle_score(angle: Angle) -> int:
    preserved = sum(WEIGHTS.get(x, 0) for x in angle.preserves)
    destroyed = 5 * len(angle.destroys)
    hits = min(corpus_hits(angle), 18)
    return 2 * preserved - destroyed + hits + angle.priority_bonus


def gf2_rank(rows: list[list[int]]) -> int:
    matrix = [row[:] for row in rows]
    rank = 0
    col = 0
    while matrix and rank < len(matrix) and col < len(matrix[0]):
        pivot = None
        for r in range(rank, len(matrix)):
            if matrix[r][col]:
                pivot = r
                break
        if pivot is None:
            col += 1
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        for r in range(len(matrix)):
            if r != rank and matrix[r][col]:
                matrix[r] = [a ^ b for a, b in zip(matrix[r], matrix[rank])]
        rank += 1
        col += 1
    return rank


def row_vector(row_name: str, cols: Iterable[str] = EXACTNESS_COLUMNS) -> list[int]:
    support = EXACTNESS_ROWS[row_name]
    return [1 if col in support else 0 for col in cols]


def hitting_sets(max_size: int = 6) -> list[tuple[str, ...]]:
    row_sets = list(EXACTNESS_ROWS.values())
    hits: list[tuple[str, ...]] = []
    for r in range(1, max_size + 1):
        for combo in combinations(EXACTNESS_COLUMNS, r):
            cset = set(combo)
            if all(cset & row for row in row_sets):
                hits.append(combo)
        if hits:
            return hits
    return hits


def tournament_edges(scored: list[tuple[Angle, int]]) -> dict[tuple[str, str], str]:
    score = {angle.key: value for angle, value in scored}
    destroy_count = {angle.key: len(angle.destroys) for angle, _ in scored}
    edges: dict[tuple[str, str], str] = {}
    for a, _ in scored:
        for b, _ in scored:
            if a == b:
                continue
            if score[a.key] > score[b.key]:
                winner = a.key
            elif score[a.key] < score[b.key]:
                winner = b.key
            elif destroy_count[a.key] < destroy_count[b.key]:
                winner = a.key
            elif destroy_count[a.key] > destroy_count[b.key]:
                winner = b.key
            else:
                winner = min(a.key, b.key)
            edges[(a.key, b.key)] = winner
    return edges


def directed_3cycles(scored: list[tuple[Angle, int]]) -> int:
    edges = tournament_edges(scored)
    keys = [angle.key for angle, _ in scored]
    count = 0
    for a, b, c in combinations(keys, 3):
        ab = edges[(a, b)]
        bc = edges[(b, c)]
        ca = edges[(c, a)]
        if ab == a and bc == b and ca == c:
            count += 1
        if ab == b and bc == c and ca == a:
            count += 1
    return count


def scc_sizes(scored: list[tuple[Angle, int]]) -> list[int]:
    keys = [angle.key for angle, _ in scored]
    edges = tournament_edges(scored)
    graph: dict[str, list[str]] = {k: [] for k in keys}
    rev: dict[str, list[str]] = {k: [] for k in keys}
    for a, b in combinations(keys, 2):
        winner = edges[(a, b)]
        loser = b if winner == a else a
        graph[winner].append(loser)
        rev[loser].append(winner)

    seen: set[str] = set()
    order: list[str] = []

    def dfs(v: str) -> None:
        seen.add(v)
        for w in graph[v]:
            if w not in seen:
                dfs(w)
        order.append(v)

    for key in keys:
        if key not in seen:
            dfs(key)

    seen.clear()
    sizes: list[int] = []

    def rdfs(v: str, acc: list[str]) -> None:
        seen.add(v)
        acc.append(v)
        for w in rev[v]:
            if w not in seen:
                rdfs(w, acc)

    for key in reversed(order):
        if key not in seen:
            acc: list[str] = []
            rdfs(key, acc)
            sizes.append(len(acc))
    return sorted(sizes, reverse=True)


def hamiltonian_path_count(scored: list[tuple[Angle, int]]) -> int:
    keys = [angle.key for angle, _ in scored]
    edges = tournament_edges(scored)
    count = 0
    for perm in permutations(keys):
        ok = True
        for a, b in zip(perm, perm[1:]):
            if edges[(a, b)] != a:
                ok = False
                break
        if ok:
            count += 1
    return count


def exactness_rows() -> list[list[int]]:
    return [row_vector(row) for row in EXACTNESS_ROWS]


def renormalization_potential_deltas() -> list[tuple[str, str, tuple[int, ...], str]]:
    deltas = []
    for src, dst, label in RENORMALIZATION_EDGES:
        sp = RENORMALIZATION_STATES[src]["potential"]
        dp = RENORMALIZATION_STATES[dst]["potential"]
        delta = tuple(a - b for a, b in zip(sp, dp))
        deltas.append((src, dst, delta, label))
    return deltas


def print_angle_table(scored: list[tuple[Angle, int]]) -> None:
    print("ANGLE PRIORITY TABLE")
    for angle, score in scored:
        print(f"  {angle.key} score={score:3d} {angle.name}")
        print(f"      family={angle.family}")
        print(f"      preserves={','.join(angle.preserves)}")
        print(f"      destroys={','.join(angle.destroys)}")
        print(f"      test={angle.proof_test}")
        print(f"      cross_links={', '.join(angle.cross_links)}")
        print(f"      risk={angle.risk}")


def print_exactness_packet() -> None:
    matrix = exactness_rows()
    rank = gf2_rank(matrix)
    hits = hitting_sets()
    print("\nFIRST-OBSTRUCTION / SHEAF EXACTNESS PACKET")
    print(f"  rows={len(EXACTNESS_ROWS)} cols={len(EXACTNESS_COLUMNS)} gf2_rank={rank}")
    print("  columns=" + ", ".join(EXACTNESS_COLUMNS))
    for row_name, vector in zip(EXACTNESS_ROWS, matrix):
        support = [c for c, bit in zip(EXACTNESS_COLUMNS, vector) if bit]
        print(f"  {row_name}: " + ", ".join(support))
    print("\n  minimal chart-sidecar hitting sets:")
    if not hits:
        print("    none up to requested size")
    else:
        for combo in hits[:12]:
            print("    " + ", ".join(combo))
        if len(hits) > 12:
            print(f"    ... {len(hits) - 12} more")
    print(
        "  readout=the sheaf route asks for exactness of overlap payloads. "
        "A quotient is legal only if its first obstruction is exact, repaired "
        "by holonomy, lifted to an endpoint chamber, descended, stopped at "
        "AP/GW, or named as debt."
    )


def print_renormalization_packet() -> None:
    print("\nFAREY-CUSP RENORMALIZATION TRANSFER PACKET")
    print("  potential=(covering_depth,cusp_height,kernel_unknown,sidecar_debt)")
    for name, data in RENORMALIZATION_STATES.items():
        print(f"  {name}: potential={data['potential']} exit={data['exit']}")
    print("\n  transition deltas:")
    for src, dst, delta, label in renormalization_potential_deltas():
        print(f"    {src} -> {dst}: delta={delta} via {label}")
    print(
        "  readout=the only bad transition is the unknown zero-open kernel. "
        "Every useful cusp move either lowers the packet to a witness, emits "
        "positive boundary-moment slack, proves the AP/GW kernel impossible in "
        "qdiv>14, or routes to K33/H7 debt."
    )


def print_cross_discipline() -> None:
    print("\nCROSS-DISCIPLINARY SYNTHESIS")
    bullets = [
        (
            "sheaf theory",
            "Cech overlaps become the actual proof surface; first forgotten "
            "payloads are cocycles, and a proof is an exactness statement with "
            "named harmonic/state-lift residues.",
        ),
        (
            "local systems",
            "zeta7 contact holonomy is a connection coordinate; it repairs "
            "lag/residue curvature locally without pretending to be a global "
            "terminal quotient.",
        ),
        (
            "modular/Farey geometry",
            "qdiv and exact-period packets should move on a cusp tree; the "
            "principal part is the packet handoff, not a decorative modular "
            "analogy.",
        ),
        (
            "Schwarz-Christoffel",
            "unit contacts look like polygon vertices; the accessory parameter "
            "is precisely the off-unit endpoint-owner debt that a polygon word "
            "would otherwise hide.",
        ),
        (
            "Mertens/loglog",
            "scale hierarchies bound proof complexity through obstruction-prime "
            "ledgers, but loglog thinness is not a witness.",
        ),
        (
            "Fermat-Catalan/Hurwitz/Markov",
            "sum reciprocal < 1 is best read as a hyperbolic triple-cusp "
            "defect; in LRC it can label route curvature only after a concrete "
            "triple of packet widths is declared.",
        ),
    ]
    for name, text in bullets:
        print(f"  {name}: {text}")


def main() -> None:
    scored = sorted(
        ((angle, angle_score(angle)) for angle in ANGLES),
        key=lambda item: (-item[1], item[0].key),
    )
    print("HYP-3301 LRC14 SHEAF EXACTNESS + CUSP TRANSFER PROOF-ANGLE SCOUT")
    print("=" * 78)
    print(
        "goal=pick two remaining proof angles different from HYP-3300's "
        "observability/Morse pair: first-obstruction sheaf exactness and "
        "Farey-cusp renormalization transfer."
    )
    print_angle_table(scored)
    print_exactness_packet()
    print_renormalization_packet()
    print_cross_discipline()

    score_hist = Counter(score for _, score in scored)
    path = " -> ".join(angle.name for angle, _ in scored)
    print("\nTOURNAMENT ANALYSIS")
    print("  vertices=proof carriers and chart transitions, not runners or arcs")
    print(
        "  pairwise_observable=which carrier preserves boundary/open status, "
        "chart gluing, covering exits, scale normality, owner/topology, "
        "residue/magnitude split, and formal packet schema while destroying "
        "fewer coordinates"
    )
    print(
        "  switch/gauge=A->B iff A has larger weighted proof-carrier score; "
        "ties prefer fewer destroyed coordinates and then stable key order"
    )
    print(f"  score_hist={dict(sorted(score_hist.items()))}")
    print(f"  directed_3cycles={directed_3cycles(scored)}")
    print(f"  scc_sizes={scc_sizes(scored)}")
    print(f"  hamiltonian_path_count={hamiltonian_path_count(scored)}")
    print(f"  priority_path={path}")

    print("\nASSUMPTION-CHALLENGE")
    print(
        "  alternate vertices considered: runners, gaps, unit contacts, "
        "endpoint sections, wall crossings, residues, exact-period packets, "
        "cover arcs, Fourier modes, matroid cocircuits, Cech overlaps, "
        "cusp states, boundary-moment kernels, and proof obligations."
    )
    print(
        "  chosen vertices: proof carriers plus chart transitions.  This "
        "preserves LRC status only through explicit boundary/open, exactness, "
        "covering-exit, scale-normal, owner/topology, residue/magnitude, and "
        "state-lift fields."
    )
    print(
        "  destroyed by runner/arc or raw analogy quotients: endpoint owner, "
        "accessory parameter, exact-period owner, cusp principal part, first "
        "obstruction cocycle, AP/GW kernel status, and H7 debt."
    )
    print(
        "  challenged assumption: the remaining proof is just a stronger "
        "extremal inequality.  This scout suggests an exactness theorem plus "
        "a cusp-transfer theorem."
    )


if __name__ == "__main__":
    main()
