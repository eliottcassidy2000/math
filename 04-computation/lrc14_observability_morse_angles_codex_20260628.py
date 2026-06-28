#!/usr/bin/env python3
"""LRC14 remaining proof-angle scout: observability plus Morse descent.

This is a synthesis script, not a proof.  It turns the current LRC14 frontier
into two testable proof programs:

1. a sidecar observability matrix whose columns are retained coordinates, and
2. a discrete-Morse / Lyapunov descent packet on finite proof chambers.

Tournament Analysis is over proof angles and sidecar carriers, not runners.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parents[1]


OBLIGATION_WEIGHTS = {
    "strict_open_witness": 8,
    "boundary_equality_core": 8,
    "bulk_equidistribution": 7,
    "tight_locus_finiteness": 8,
    "compression_legality": 7,
    "sidecar_reconstructability": 7,
    "finite_discharge": 8,
    "state_lift_exit": 6,
    "formalizable_certificate": 6,
    "odd_sign_payload": 7,
    "scale_lift_stability": 6,
}


KEY_FILES = [
    "00-navigation/LRC-LENS-MAP.md",
    "00-navigation/CONCEPT-MAP.md",
    "00-navigation/OPEN-QUESTIONS.md",
    "00-navigation/LRC-TECHNIQUE-INDEX.md",
    "00-navigation/LRC-TOURNAMENT-TECHNIQUE-INDEX.md",
    "05-knowledge/hypotheses/INDEX.md",
    "05-knowledge/hypotheses/HYP-3048-expanded-tournament-matrix-atlas.md",
    "05-knowledge/hypotheses/HYP-3069-medianized-route-center-gate-lrc14.md",
    "05-knowledge/hypotheses/HYP-3070-route-triple-center-control-lrc14.md",
    "05-knowledge/hypotheses/HYP-3108-lee-yang-savitch-ear-lattice-extremality.md",
    "05-knowledge/hypotheses/HYP-3124-lrc14-tournament-edge-witness-recursion.md",
    "05-knowledge/hypotheses/HYP-3225-lrc14-green-lorentzian-trap-fingerprints.md",
    "05-knowledge/hypotheses/HYP-3236-lrc14-green-conductance-connectivity.md",
    "07-reflections/lonely-runner-as-chebyshev-equioscillation.md",
    "05-knowledge/hypotheses/HYP-3243-lrc14-topology-geometry-graph-routes.md",
    "05-knowledge/hypotheses/HYP-3244-tiling-half-tiling-interlocking-recursions.md",
    "05-knowledge/hypotheses/HYP-3245-lrc14-equioscillation-autocorrelation-atlas.md",
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
        "A01",
        "sidecar_observability_matrix",
        "linear-control / coding-syndrome route",
        (
            "sidecar_reconstructability",
            "compression_legality",
            "formalizable_certificate",
            "odd_sign_payload",
            "state_lift_exit",
        ),
        ("raw_scalar_shadow",),
        (
            "Build a finite residual-pair by sidecar-column matrix; every "
            "route/status-changing pair must be separated, reconstructed, "
            "dual-annihilated, descended, or routed to named debt."
        ),
        (
            "HYP-3048",
            "HYP-3071",
            "coding syndromes",
            "control observability",
            "Farkas/KKT certificates",
            "MacWilliams/Delsarte",
        ),
        "needs real residual-pair rows, not just proposed columns",
        9,
    ),
    Angle(
        "A02",
        "finite_chamber_discrete_morse_descent",
        "topological-combinatorial descent route",
        (
            "finite_discharge",
            "tight_locus_finiteness",
            "strict_open_witness",
            "boundary_equality_core",
            "state_lift_exit",
        ),
        ("continuous_measure_shadow",),
        (
            "Define a vector-valued Morse energy on chamber packets; prove "
            "every nonterminal chamber has a strictly descending legal wall."
        ),
        (
            "HYP-3246",
            "HYP-3247",
            "HYP-3248",
            "HYP-3249",
            "HYP-3243",
            "Forman discrete Morse",
            "oriented matroids",
            "Cech barcodes",
            "Gauss-sum index",
            "Borsuk-Ulam forcing gap",
            "chip-firing",
            "Green/Rayleigh energy",
        ),
        "must prove acyclicity of the matching, not just local descent",
        8,
    ),
    Angle(
        "A03",
        "green_toeplitz_thomson_energy",
        "electrical-network / moment-cone route",
        (
            "finite_discharge",
            "strict_open_witness",
            "formalizable_certificate",
        ),
        ("negative_covariance_payload", "chart_change_class"),
        (
            "Use Toeplitz normal-fan slack plus Green resistance excess as "
            "a two-coordinate certificate for the finite trap boundary."
        ),
        (
            "HYP-3225",
            "HYP-3236",
            "Thomson principle",
            "Rayleigh monotonicity",
            "M-matrices",
        ),
        "positive-part conductance loses odd/negative leakage unless retained",
        5,
    ),
    Angle(
        "A04",
        "lag_autocorrelation_transport",
        "harmonic-analysis / signal-processing route",
        (
            "strict_open_witness",
            "bulk_equidistribution",
            "formalizable_certificate",
        ),
        ("endpoint_owner", "state_lift_exit"),
        (
            "Turn HYP-3245 low-lag deficit plus outward-lag surplus into a "
            "transport inequality feeding AP support or Green/Toeplitz slack."
        ),
        (
            "HYP-3245",
            "Fejer kernel",
            "Delsarte LP",
            "Remez equiripple",
            "Golay/Barker sidecars",
        ),
        "single lag curve cannot see endpoint owners or H=7 state-lifts",
        4,
    ),
    Angle(
        "A05",
        "odd_negative_duality_resurrection",
        "parity/sign-isotypic route",
        (
            "odd_sign_payload",
            "boundary_equality_core",
            "sidecar_reconstructability",
        ),
        ("positive_even_compression",),
        (
            "Promote the odd/negative coordinate to a resurrection sidecar "
            "so false positive-even terminal rows become visible debts."
        ),
        (
            "HYP-3238",
            "HYP-3239",
            "Borsuk-Ulam degree",
            "Hermite-Biehler",
            "sign representation",
        ),
        "can become a label-only guard unless attached to residual rows",
        4,
    ),
    Angle(
        "A06",
        "tiling_half_tiling_descent_span",
        "finite-cover / quotient-descent route",
        (
            "compression_legality",
            "sidecar_reconstructability",
            "strict_open_witness",
        ),
        ("raw_A000568_count", "fixed_path_cube_count"),
        (
            "Lift witnesses in the fixed-path tiling cover and compress only "
            "after fiber, parent-aut, rectangle/hourglass, and deletion "
            "sidecars certify descent."
        ),
        (
            "HYP-3244",
            "A000568",
            "Burnside fibers",
            "rectangle coboundaries",
            "partial cubes",
        ),
        "finite quotient evidence must be tied to LRC predicate preservation",
        4,
    ),
    Angle(
        "A07",
        "roth_halasz_hensel_lift_packet",
        "analytic discrepancy plus local p-adic route",
        (
            "bulk_equidistribution",
            "scale_lift_stability",
            "formalizable_certificate",
        ),
        ("finite_chamber_address", "odd_sign_payload"),
        (
            "Treat Roth-Halasz discrepancy as a bulk descent bound and "
            "Hensel-Krasner as the local unit that prevents scale-lift drift."
        ),
        (
            "S289",
            "Roth-Halasz discrepancy",
            "Hensel-Krasner",
            "Ramanujan-Soldner boundary",
            "rigid analytic lifting",
        ),
        "incoming coordination entry is schematic and still needs row-level data",
        2,
    ),
    Angle(
        "A08",
        "lee_yang_ear_root_motion",
        "PGF roots / ear-decomposition route",
        (
            "state_lift_exit",
            "sidecar_reconstructability",
            "finite_discharge",
        ),
        ("bulk_equidistribution",),
        (
            "Use root-motion words and one-runner ear payloads as observability "
            "columns for trap exits and forbidden H=7 state-lifts."
        ),
        (
            "HYP-3108",
            "HYP-3109",
            "HYP-3112",
            "Lee-Yang zeros",
            "ear decompositions",
            "Savitch midpoints",
        ),
        "root signatures are diagnostics unless tied to legal packet exits",
        3,
    ),
    Angle(
        "A09",
        "median_route_center_closure",
        "median algebra / scheduler route",
        (
            "compression_legality",
            "sidecar_reconstructability",
            "formalizable_certificate",
        ),
        ("numeric_extremality_signal",),
        (
            "Close proof-state triples under legal sidecar closure; a nonunique "
            "or empty center emits the first missing sidecar or named debt."
        ),
        (
            "HYP-3069",
            "HYP-3070",
            "median graphs",
            "Horn closure",
            "partial cubes",
        ),
        "scheduler route organizes proof obligations but may not produce a new inequality",
        1,
    ),
    Angle(
        "A10",
        "raw_single_scalar_extremality",
        "discarded baseline",
        ("strict_open_witness",),
        (
            "endpoint_owner",
            "odd_sign_payload",
            "finite_chamber_address",
            "state_lift_exit",
            "compression_legality",
        ),
        (
            "Try to prove AP by a single scalar such as p0, lambda2, lag cost, "
            "or root radius."
        ),
        ("prior false-terminal audits",),
        "kept only as a warning; too many false terminals in prior scouts",
        -20,
    ),
)


SIDECAR_COLUMNS = [
    "endpoint_owner",
    "boundary_cocircuit",
    "phi_witness_address",
    "dilation_grid",
    "toeplitz_normal_fan_slack",
    "green_resistance_slack",
    "odd_negative_payload",
    "lee_yang_root_word",
    "tiling_descent_packet",
    "lag_transport_signature",
    "unit_equioscillation_index",
    "binding_complement_pair_word",
    "analytic_topological_index_equalizer",
    "gauss_sum_index_word",
    "borsuk_ulam_forcing_gap",
    "roth_halasz_discrepancy",
    "hensel_krasner_unit",
    "state_lift_H7",
    "cech_betti_hole",
    "ear_payload",
]


OBSERVABILITY_ROWS = {
    "open_safe_vs_closed_boundary": {
        "endpoint_owner",
        "boundary_cocircuit",
        "cech_betti_hole",
    },
    "AP_vs_GW_phi14_equality": {
        "phi_witness_address",
        "boundary_cocircuit",
        "odd_negative_payload",
        "unit_equioscillation_index",
        "binding_complement_pair_word",
    },
    "chebyshev_units_vs_covering_floor": {
        "unit_equioscillation_index",
        "binding_complement_pair_word",
        "roth_halasz_discrepancy",
        "cech_betti_hole",
    },
    "index_prediction_vs_naive_bu_gap": {
        "unit_equioscillation_index",
        "analytic_topological_index_equalizer",
        "gauss_sum_index_word",
        "borsuk_ulam_forcing_gap",
        "cech_betti_hole",
    },
    "base_phi14_vs_dilation_phi14d": {
        "phi_witness_address",
        "dilation_grid",
        "hensel_krasner_unit",
    },
    "toeplitz_pair_plucker_trap": {
        "toeplitz_normal_fan_slack",
        "green_resistance_slack",
        "lag_transport_signature",
    },
    "green_low_connectivity_trap": {
        "green_resistance_slack",
        "endpoint_owner",
        "tiling_descent_packet",
    },
    "odd_negative_false_terminal": {
        "odd_negative_payload",
        "lee_yang_root_word",
        "ear_payload",
    },
    "lee_yang_root_collision_vs_confinement": {
        "lee_yang_root_word",
        "ear_payload",
        "state_lift_H7",
    },
    "tiling_lift_vs_half_tiling_descent": {
        "tiling_descent_packet",
        "endpoint_owner",
        "boundary_cocircuit",
    },
    "bulk_discrepancy_vs_core_witness": {
        "roth_halasz_discrepancy",
        "phi_witness_address",
        "unit_equioscillation_index",
        "cech_betti_hole",
        "lag_transport_signature",
    },
    "state_lift_H7_residual": {
        "state_lift_H7",
        "hensel_krasner_unit",
        "ear_payload",
    },
    "scale_lift_unit_vs_chart_drift": {
        "hensel_krasner_unit",
        "dilation_grid",
        "tiling_descent_packet",
    },
}


MORSE_COORDINATES = [
    "tope_boundary_depth",
    "toeplitz_deficit_order",
    "green_resistance_excess",
    "lag_transport_cost",
    "unit_equioscillation_defect",
    "index_equality_defect",
    "gauss_sum_word_depth",
    "bu_forcing_gap_status",
    "root_collision_depth",
    "odd_negative_debt",
    "chart_change_depth",
    "hensel_unit_failure",
    "state_lift_distance",
]


DESCENT_CASES = [
    (
        "strict_open_tope",
        "terminal safe witness",
        ("tope_boundary_depth", "cech_betti_hole"),
    ),
    (
        "chebyshev_unit_equioscillation",
        "unit-group boundary packet or covering-floor handoff",
        ("unit_equioscillation_index", "binding_complement_pair_word"),
    ),
    (
        "index_theorem_prediction_boundary",
        "analytic/topological index equality or named Borsuk-Ulam forcing gap",
        (
            "analytic_topological_index_equalizer",
            "gauss_sum_index_word",
            "borsuk_ulam_forcing_gap",
        ),
    ),
    (
        "AP_or_GW_phi14_boundary",
        "closed equality atom, not a trap",
        ("phi_witness_address", "unit_equioscillation_index", "odd_negative_payload"),
    ),
    (
        "phi14d_dilation_witness",
        "constructed equality on finer witness grid",
        ("dilation_grid", "hensel_krasner_unit"),
    ),
    (
        "finite_toeplitz_green_trap",
        "descend by moment-cone or electrical slack",
        ("toeplitz_normal_fan_slack", "green_resistance_slack"),
    ),
    (
        "lag_transport_outward_mass",
        "push outward surplus into AP support or trap sidecar",
        ("lag_transport_signature", "green_resistance_slack"),
    ),
    (
        "odd_negative_false_terminal",
        "retain parity/sign debt before positive compression",
        ("odd_negative_payload", "lee_yang_root_word"),
    ),
    (
        "bulk_discrepancy_chamber",
        "Roth-Halasz bulk descent unless core witness takes over",
        ("roth_halasz_discrepancy", "phi_witness_address"),
    ),
    (
        "local_scale_lift_chamber",
        "Hensel-Krasner unit prevents scale-lift drift",
        ("hensel_krasner_unit", "dilation_grid"),
    ),
    (
        "state_lift_H7_exit",
        "forbidden tournament-state contradiction or named debt",
        ("state_lift_H7", "ear_payload"),
    ),
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
                total += min(text.count(token), 7)
    return total


def angle_score(angle: Angle) -> int:
    preserved = sum(OBLIGATION_WEIGHTS.get(x, 0) for x in angle.preserves)
    destroyed = len(angle.destroys) * 4
    hits = min(corpus_hits(angle), 16)
    return 2 * preserved - destroyed + hits + angle.priority_bonus


def gf2_rank(rows: list[list[int]]) -> int:
    matrix = [row[:] for row in rows]
    rank = 0
    col = 0
    while rank < len(matrix) and col < len(matrix[0]):
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


def row_vector(row_name: str, cols: Iterable[str] = SIDECAR_COLUMNS) -> list[int]:
    active = OBSERVABILITY_ROWS[row_name]
    return [1 if col in active else 0 for col in cols]


def hitting_sets(max_size: int = 4) -> list[tuple[str, ...]]:
    row_sets = list(OBSERVABILITY_ROWS.values())
    cols = SIDECAR_COLUMNS
    hits: list[tuple[str, ...]] = []
    for r in range(1, max_size + 1):
        for combo in combinations(cols, r):
            cset = set(combo)
            if all(cset & row for row in row_sets):
                hits.append(combo)
        if hits:
            return hits
    return hits


def tournament_edges(scored: list[tuple[Angle, int]]) -> dict[tuple[str, str], str]:
    scores = {angle.key: score for angle, score in scored}
    destroy_counts = {angle.key: len(angle.destroys) for angle, _ in scored}
    edges: dict[tuple[str, str], str] = {}
    for a, _ in scored:
        for b, _ in scored:
            if a == b:
                continue
            if scores[a.key] > scores[b.key]:
                winner = a.key
            elif scores[a.key] < scores[b.key]:
                winner = b.key
            elif destroy_counts[a.key] < destroy_counts[b.key]:
                winner = a.key
            elif destroy_counts[a.key] > destroy_counts[b.key]:
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


def hamiltonian_path_count(scored: list[tuple[Angle, int]]) -> int:
    # For this score-oriented tournament, ties are deterministically resolved,
    # so the sorted score order is the canonical path.  Count whether any equal
    # score blocks create additional compatible permutations.
    counts = Counter(score for _, score in scored)
    paths = 1
    for block in counts.values():
        for x in range(2, block + 1):
            paths *= x
    return paths


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


def print_observability_matrix() -> None:
    rows = list(OBSERVABILITY_ROWS)
    matrix = [row_vector(row) for row in rows]
    rank = gf2_rank(matrix)
    hits = hitting_sets(max_size=5)
    print("\nSIDECAR OBSERVABILITY MATRIX")
    print(f"  rows={len(rows)} cols={len(SIDECAR_COLUMNS)} gf2_rank={rank}")
    print("  columns=" + ", ".join(SIDECAR_COLUMNS))
    for row_name, vector in zip(rows, matrix):
        support = [c for c, bit in zip(SIDECAR_COLUMNS, vector) if bit]
        print(f"  {row_name}: " + ", ".join(support))
    print("\n  minimal hitting sets among proposed columns:")
    if not hits:
        print("    none up to requested size")
    else:
        for combo in hits[:12]:
            print("    " + ", ".join(combo))
        if len(hits) > 12:
            print(f"    ... {len(hits) - 12} more")
    print(
        "  readout=the proof matrix is not asking for every scalar at once; "
        "it asks for the smallest column basis that separates every row where "
        "a quotient could change the LRC status."
    )


def print_morse_packet() -> None:
    print("\nDISCRETE MORSE / LYAPUNOV DESCENT PACKET")
    print("  coordinates=" + ", ".join(MORSE_COORDINATES))
    for case, exit_text, columns in DESCENT_CASES:
        print(f"  {case}: exit={exit_text}; retained={', '.join(columns)}")
    print(
        "  proposed energy order="
        "(strict-open deficit, equality-core distance, finite trap slack, "
        "odd/negative debt, scale-lift debt, state-lift distance)"
    )
    print(
        "  proof obligation=construct an acyclic matching on chamber packets. "
        "Critical cells must be AP/GW Phi14 equality, Phi14d dilation equality, "
        "strict open safe cells, state-lift H=7 contradictions, or named debt."
    )


def print_cross_discipline_synthesis() -> None:
    print("\nCROSS-DISCIPLINARY SYNTHESIS")
    bullets = [
        (
            "control/coding",
            "observability columns are syndrome bits; residual pairs are error "
            "classes; MacWilliams/Delsarte language explains why autocorrelation "
            "and weight enumerators are useful but not terminal alone.",
        ),
        (
            "optimization",
            "Farkas/KKT/Schur certificates make each sidecar column a missing "
            "dual variable instead of a metaphor.",
        ),
        (
            "oriented matroids/topology",
            "topes, cocircuits, Cech bars, and danger-cover holes provide the "
            "finite chamber complex on which the Morse descent should live.",
        ),
        (
            "electrical networks",
            "Toeplitz curvature and Green resistance define a Thomson/Rayleigh "
            "energy, while negative covariance is leakage payload.",
        ),
        (
            "arithmetic dynamics",
            "Roth-Halasz controls bulk discrepancy; Hensel-Krasner is the local "
            "unit rule preventing scale-lift drift at the 7 x 2 apex.",
        ),
        (
            "circuit complexity",
            "the bounded proof circuit should expose missing-input vectors: "
            "small columns first, named sidecar debt when a gate would forget a "
            "coordinate.",
        ),
    ]
    for name, text in bullets:
        print(f"  {name}: {text}")


def main() -> None:
    scored = sorted(((angle, angle_score(angle)) for angle in ANGLES), key=lambda x: (-x[1], x[0].key))
    print("HYP-3300 LRC14 OBSERVABILITY + MORSE PROOF-ANGLE SCOUT")
    print("=" * 78)
    print(
        "goal=pick two remaining LRC14 proof angles that differ from the "
        "recent chamber-cover atlas, while importing the useful columns from "
        "HYP-3243/HYP-3244/HYP-3245 and S289."
    )
    print_angle_table(scored)
    print_observability_matrix()
    print_morse_packet()
    print_cross_discipline_synthesis()

    score_hist = Counter(score for _, score in scored)
    path = " -> ".join(angle.name for angle, _ in scored)
    print("\nTOURNAMENT ANALYSIS")
    print("  vertices=proof angles and sidecar carriers, not runners, arcs, or sectors")
    print(
        "  pairwise_observable=which angle preserves more LRC proof obligations "
        "while destroying fewer scarce coordinates"
    )
    print(
        "  switch/gauge=A->B iff A has larger weighted obligation score; ties "
        "prefer fewer destroyed coordinates and then stable key order"
    )
    print(f"  score_hist={dict(sorted(score_hist.items()))}")
    print(f"  directed_3cycles={directed_3cycles(scored)}")
    print("  scc_sizes=[1,1,1,1,1,1,1,1,1,1]")
    print(f"  hamiltonian_path_count={hamiltonian_path_count(scored)}")
    print(f"  priority_path={path}")

    print("\nASSUMPTION-CHALLENGE")
    print(
        "  alternate vertices considered: runners, speed gaps, fixed circle "
        "sections, endpoint boundaries, wall crossings, residues, cover arcs, "
        "Fourier modes, matroid circuits, proof obligations, sidecar columns, "
        "finite chambers, and state-lift exits."
    )
    print(
        "  chosen vertices: proof angles plus sidecar carriers.  This preserves "
        "the LRC predicate only through explicit strict-open, boundary-equality, "
        "finite-discharge, compression-legality, bulk/core, and state-lift exits."
    )
    print(
        "  destroyed by the tempting runner/arc quotient: endpoint owners, "
        "odd/negative sign, finite chamber address, tiling descent packet, "
        "root-motion word, Hensel unit, and H=7 state-lift obligation."
    )
    print(
        "  challenged assumption: the next proof step should be another scalar "
        "extremality.  This scout suggests a column-basis theorem plus an "
        "acyclic chamber descent theorem instead."
    )


if __name__ == "__main__":
    main()
