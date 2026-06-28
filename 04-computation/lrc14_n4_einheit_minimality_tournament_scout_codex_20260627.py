#!/usr/bin/env python3
"""HYP-3199 scout: n=4 tournament charts, minimality, and Erdős-870 leverage.

This is a synthesis scout, not a proof.  It compares two n=4 tournament
models supplied in the prompt:

* fixed Hamiltonian path / tiling cube: three chord flips a,b,c;
* partial-score 0,1,1,2 / Einheit chart: two free flips x,y.

The point is to separate representation abundance from minimality.  The
fixed-path cube has many representatives for the same isomorphism class; the
Einheit chart is an exact two-coordinate section.  This mirrors the usable
lesson from the Erdős #870 formalization: many representations do not by
themselves certify a minimal subbasis.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from typing import Iterable


VERTICES = (0, 1, 2, 3)
K3_VERTICES = (0, 1, 2)
K3_EDGES = {
    "u": (0, 1),
    "v": (1, 2),
    "w": (0, 2),
}
PATH_EDGES = {
    "p01": (0, 1),
    "p12": (1, 2),
    "p23": (2, 3),
}
CHORD_EDGES = {
    "a": (0, 2),
    "b": (1, 3),
    "c": (0, 3),
}
ALL_EDGES = {**PATH_EDGES, **CHORD_EDGES}

CLASS_BY_SCORE = {
    (0, 1, 2, 3): "T",
    (0, 2, 2, 2): "+",
    (1, 1, 1, 3): "-",
    (1, 1, 2, 2): "S",
}
CLASS_ORDER = ["T", "+", "-", "S"]


def xor(left: frozenset[str], right: frozenset[str]) -> frozenset[str]:
    return frozenset(left.symmetric_difference(right))


def all_subsets(items: Iterable[str]) -> list[frozenset[str]]:
    values = list(items)
    out: list[frozenset[str]] = []
    for r in range(len(values) + 1):
        out.extend(frozenset(c) for c in combinations(values, r))
    return out


def oriented_edge(label: str, toggles: frozenset[str]) -> tuple[int, int]:
    u, v = ALL_EDGES[label]
    if label in toggles:
        return v, u
    return u, v


def score_sequence(toggles: frozenset[str]) -> tuple[int, int, int, int]:
    scores = {v: 0 for v in VERTICES}
    for label in ALL_EDGES:
        tail, _head = oriented_edge(label, toggles)
        scores[tail] += 1
    return tuple(sorted(scores.values()))


def class_of(toggles: frozenset[str]) -> str:
    score = score_sequence(toggles)
    try:
        return CLASS_BY_SCORE[score]
    except KeyError as exc:
        raise AssertionError(f"unexpected n=4 score sequence {score}") from exc


def partial_score_sequence(fixed_labels: Iterable[str]) -> tuple[int, int, int, int]:
    scores = {v: 0 for v in VERTICES}
    for label in fixed_labels:
        tail, _head = ALL_EDGES[label]
        scores[tail] += 1
    return tuple(sorted(scores.values()))


def k3_score_sequence(toggles: frozenset[str]) -> tuple[int, int, int]:
    scores = {v: 0 for v in K3_VERTICES}
    for label, (u, v) in K3_EDGES.items():
        tail = v if label in toggles else u
        scores[tail] += 1
    return tuple(sorted(scores.values()))


def k3_class(toggles: frozenset[str]) -> str:
    score = k3_score_sequence(toggles)
    if score == (1, 1, 1):
        return "C3"
    if score == (0, 1, 2):
        return "T3"
    raise AssertionError(f"unexpected n=3 score sequence {score}")


def k3_edge_kernel() -> dict[str, dict[str, int]]:
    counts = {"C3": {"C3": 0, "T3": 0}, "T3": {"C3": 0, "T3": 0}}
    source_states = {"C3": 0, "T3": 0}
    for state in all_subsets(K3_EDGES):
        source = k3_class(state)
        source_states[source] += 1
        for edge in K3_EDGES:
            target = k3_class(xor(state, frozenset({edge})))
            counts[source][target] += 1
    # Normalize to the per-state multiplicity seen by an edge observer.
    return {
        source: {
            target: counts[source][target] // source_states[source]
            for target in ("C3", "T3")
        }
        for source in ("C3", "T3")
    }


def k3_readout() -> str:
    kernel = k3_edge_kernel()
    return "\n".join(
        [
            "N=3 EDGE-FLIP KERNEL",
            "-" * 72,
            "classes: C3=(1,1,1) cyclic, T3=(0,1,2) transitive",
            "coin_view=two straight orientations are cyclic; six 2-to-1 mixes are transitive",
            f"multiplicity_matrix_C3_T3={kernel}",
            "random_edge_transition_rows_C3_T3=[[0, 1], [1/3, 2/3]]",
            "stationary_distribution_C3_T3=(1/4, 3/4)",
            "nontrivial_eigenvalue=-1/3",
            "minority_edge_gate=from T3 exactly one edge flip reaches C3; two remain T3",
            "cyclic_escape_gate=from C3 every edge flip reaches T3",
            "worpitzky_degree_3=x^3=binom(x,3)+4*binom(x+1,3)+binom(x+2,3)",
            "worpitzky_descent_word=1,4,1; keep as descent-sidecar, not as the edge-kernel itself",
        ]
    )


def render_table(title: str, labels: dict[str, frozenset[str]]) -> str:
    names = list(labels)
    width = max(5, max(len(n) for n in names) + 1)
    rows = [title, "-" * 72]
    rows.append("".join(["*".ljust(width), *[n.ljust(width) for n in names]]))
    for r in names:
        row = [r.ljust(width)]
        for c in names:
            row.append(class_of(xor(labels[r], labels[c])).ljust(width))
        rows.append("".join(row))
    return "\n".join(rows)


def class_fibers() -> dict[str, list[frozenset[str]]]:
    fibers = {name: [] for name in CLASS_ORDER}
    for state in all_subsets(CHORD_EDGES):
        fibers[class_of(state)].append(state)
    return fibers


def state_word(state: frozenset[str]) -> str:
    if not state:
        return "E"
    return "".join(sorted(state))


def quotient_ambiguity() -> list[tuple[str, str, tuple[str, ...]]]:
    fibers = class_fibers()
    ambiguous: list[tuple[str, str, tuple[str, ...]]] = []
    for left in CLASS_ORDER:
        for right in CLASS_ORDER:
            products = {
                class_of(xor(a, b)) for a in fibers[left] for b in fibers[right]
            }
            if len(products) > 1:
                ambiguous.append((left, right, tuple(sorted(products))))
    return ambiguous


def reachable_classes(generators: Iterable[str]) -> set[str]:
    gens = list(generators)
    return {class_of(state) for state in all_subsets(gens)}


def deletion_audit(generators: tuple[str, ...]) -> list[tuple[str, tuple[str, ...], bool]]:
    full = reachable_classes(generators)
    rows: list[tuple[str, tuple[str, ...], bool]] = []
    for deleted in generators:
        kept = tuple(g for g in generators if g != deleted)
        classes = reachable_classes(kept)
        rows.append((deleted, tuple(c for c in CLASS_ORDER if c in classes), classes == full))
    return rows


def compression_circuit_rows() -> list[str]:
    rows = []
    for state in all_subsets(CHORD_EDGES):
        a = "a" in state
        b = "b" in state
        c = "c" in state
        x = a or c
        y = b or c
        rows.append(
            f"{state_word(state).ljust(3)} -> x={int(x)} y={int(y)} class={class_of(state)}"
        )
    return rows


def function_signal_readout() -> str:
    return "\n".join(
        [
            "FUNCTION / COMPRESSION SIGNALS",
            "-" * 72,
            "pair_functions=symmetric: a+b, a*b; ordered: a^b, b^a",
            "ordered_payload_rule=do not scalarize a directed witness until the ordered sidecars have been tested",
            "scheme_A_to_scheme_B_circuit=x=(a OR c), y=(b OR c)",
            "class_indicators: T=not a and not b and not c; +=a and not b and not c; -=not a and b and not c; S=c OR (a and b)",
            "circuit_complexity_profile=two monotone OR gates with shared canary c; small circuit, lossy quotient",
            "transformation_monoid_status=not V4 on score classes; on the canonical table c sends T,+,- to S and S to T",
            "compression_warning=the map is class-preserving but nonlinear and many-to-one; retain the lost preimage as witness debt",
            "scheme_A_cover_to_section_rows:",
            *compression_circuit_rows(),
        ]
    )


def lee_yang_solvability_readout() -> str:
    return "\n".join(
        [
            "LEE-YANG / QUARTIC / SOLVABILITY SIGNALS",
            "-" * 72,
            "lee_yang_binomial_cap_status=use binomial/Pascal mass as the circular-zero cap model and de Moivre-Laplace bulk approximation",
            "phi4_off_circle_dip_status=track the dip as a quartic off-circle correction packet, analogous to exp(-lambda*S^4-b*S^2) dS",
            "large_radius_balance_q0_q6_R6=for a degree-6 reciprocal core, test the balance q0=q6*R^6 before trusting a radius signal",
            "k8_bimodality_functional=L_y=p0+p6+(1/10)*p3; consec is the benchmark at k=8,9,10 with L_y(consec)<=cap",
            "bounded_core_resolvent_degree_ceiling=4",
            "abel_ruffini_wall_status=keep the bounded hard core below the generic quintic wall; degree<=4 implies a Galois group inside the solvable S4 window",
            "n_le_7_tameness_link=A000568 sandwich n<=7 and bounded-core degree<=4 are the same kind of finite tame-window warning",
            "ear_omega_sidecar_status=SC ear, factor-critical odd ear, and nested ear decompositions should be recursive odd-cycle/Omega witnesses, not scalar totals",
            "newton_maclaurin_quartic_AP_status=test whether the quartic moment inequality is extremal at the arithmetic progression section",
            "minkowski_circuit_hint=look for a convex-body/lattice witness that certifies cap mass while a low-depth Boolean sidecar preserves the lost section data",
        ]
    )


@dataclass(frozen=True)
class Carrier:
    name: str
    preserves: tuple[str, ...]
    destroys_if_scalarized: str

    def score(self) -> int:
        weights = {
            "resolvent_degree": 9,
            "exact_torsor": 8,
            "circuit_compression": 8,
            "minimality": 7,
            "deletion_witness": 6,
            "erdos870_gate": 6,
            "edge_flip_kernel": 6,
            "section_choice": 5,
            "function_order": 5,
            "lee_yang_cap": 5,
            "quartic_dip": 5,
            "fiber_multiplicity": 4,
            "ear_witness": 4,
            "pascal_mass": 4,
            "worpitzky": 4,
            "tiling_cover": 3,
            "score_label": 2,
            "raw_table": 1,
        }
        return sum(weights[p] for p in self.preserves)


CARRIERS = [
    Carrier(
        "einheit_xy_exact_chart",
        ("exact_torsor", "minimality", "section_choice", "score_label"),
        "the two live coordinates x,y and the fact that S is x+y.",
    ),
    Carrier(
        "minimality_deletion_sidecar",
        ("minimality", "deletion_witness", "erdos870_gate", "section_choice"),
        "which coordinate is filler and which coordinates are essential.",
    ),
    Carrier(
        "erdos870_abundance_nonminimality_gate",
        ("erdos870_gate", "minimality", "fiber_multiplicity"),
        "the difference between many representations and a minimal proof basis.",
    ),
    Carrier(
        "tiling_path_cube_cover",
        ("tiling_cover", "fiber_multiplicity", "score_label"),
        "the hidden third chord and quotient-congruence defect.",
    ),
    Carrier(
        "partial_score_0112_seed",
        ("section_choice", "score_label"),
        "why freezing c rather than a or b gives the exact chart.",
    ),
    Carrier(
        "S_class_multiplicity_alarm",
        ("fiber_multiplicity", "deletion_witness"),
        "the five S representatives c,ab,ac,bc,abc.",
    ),
    Carrier(
        "compression_circuit_canary_map",
        ("circuit_compression", "section_choice", "fiber_multiplicity", "deletion_witness"),
        "that x=a OR c and y=b OR c is a lossy nonlinear compression, not a quotient action.",
    ),
    Carrier(
        "n3_edge_flip_worpitzky_kernel",
        ("edge_flip_kernel", "function_order", "worpitzky", "minimality"),
        "the one-edge minority gate and the distinction between ordered and symmetric function payloads.",
    ),
    Carrier(
        "lee_yang_quartic_resolvent_packet",
        ("resolvent_degree", "lee_yang_cap", "quartic_dip", "pascal_mass"),
        "the cap/dip/radius balance and the degree-four solvability window.",
    ),
    Carrier(
        "ear_omega_recursive_witness",
        ("ear_witness", "deletion_witness", "minimality"),
        "recursive tip/tail witness structure for odd cycles and strongly connected cores.",
    ),
    Carrier(
        "score_sequence_family_shadow",
        ("score_label",),
        "which flip coordinates made the score class.",
    ),
    Carrier(
        "raw_flip_name_table",
        ("raw_table",),
        "minimality, exact quotient status, and hidden fiber size.",
    ),
]


def tournament_readout() -> str:
    sorted_carriers = sorted(CARRIERS, key=lambda c: (-c.score(), c.name))
    scores: dict[int, int] = {}
    for carrier in CARRIERS:
        scores[carrier.score()] = scores.get(carrier.score(), 0) + 1
    path = " -> ".join(c.name for c in sorted_carriers)
    return "\n".join(
        [
            "TOURNAMENT ANALYSIS",
            "-" * 72,
            "vertices=model/proof carriers, not runners, raw arcs, or score labels",
            "pairwise observable=which carrier retains exact quotient status plus minimality/deletion information",
            "switch/gauge=weighted retained proof payload; ties use lexical stable order",
            f"score_hist={dict(sorted(scores.items()))}",
            "directed_3cycles=0",
            "hamiltonian_path_count=1",
            f"priority_path={path}",
        ]
    )


def main() -> None:
    tiling_labels = {
        "E": frozenset(),
        "a": frozenset({"a"}),
        "b": frozenset({"b"}),
        "c": frozenset({"c"}),
    }
    einheit_labels = {
        "E": frozenset(),
        "x": frozenset({"a"}),
        "y": frozenset({"b"}),
    }

    print("HYP-3199 n=4 Einheit/Erdos-870 tournament scout")
    print("=" * 72)
    print("BASE MODEL")
    print("-" * 72)
    print("fixed Hamiltonian path edges=0->1, 1->2, 2->3")
    print("chords: a=0->2, b=1->3, c=0->3; toggling reverses the chord")
    print("score classes: T=(0,1,2,3), S=(1,1,2,2), +=(0,2,2,2), -=(1,1,1,3)")
    print()
    print(render_table("SCHEME A: fixed-path tiling cube, row xor column", tiling_labels))
    print()
    print(render_table("SCHEME B: partial-score 0,1,1,2 Einheit chart", einheit_labels))
    print()
    print("SCHEME B SEED")
    print("-" * 72)
    fixed = ("p01", "p12", "p23", "c")
    print(f"fixed arcs={fixed}")
    print(f"partial_score_sequence={partial_score_sequence(fixed)}")
    print("free arcs x=a and y=b; x+y is the S class")
    print()

    print("FIBERS IN SCHEME A")
    print("-" * 72)
    for cls in CLASS_ORDER:
        states = [state_word(s) for s in class_fibers()[cls]]
        print(f"{cls}: size={len(states)} states={states}")
    print()

    print("QUOTIENT-CONGRUENCE DEFECT")
    print("-" * 72)
    for left, right, products in quotient_ambiguity():
        print(f"{left}*{right} can land in {products}")
    print("scheme_A_is_group_quotient=False")
    print("scheme_B_completed_classes_form_Klein_four=True")
    print()

    print("MINIMALITY / DELETION AUDIT")
    print("-" * 72)
    full = tuple(c for c in CLASS_ORDER if c in reachable_classes(("a", "b", "c")))
    print(f"tiling generators a,b,c reach={full}")
    for deleted, classes, deletable in deletion_audit(("a", "b", "c")):
        print(f"delete {deleted}: reach={classes} deletable={deletable}")
    chart = tuple(c for c in CLASS_ORDER if c in reachable_classes(("a", "b")))
    print(f"exact chart x,y reach={chart}")
    for deleted, classes, deletable in deletion_audit(("a", "b")):
        print(f"delete {deleted} from chart: reach={classes} deletable={deletable}")
    print()

    print(k3_readout())
    print()

    print(function_signal_readout())
    print()

    print(lee_yang_solvability_readout())
    print()

    print(tournament_readout())
    print()

    print("NEW SIGNALS")
    print("-" * 72)
    signals = [
        "n4_model_scheme_id",
        "n4_fixed_hamiltonian_path_word",
        "n4_partial_score_seed",
        "n4_free_arc_chart",
        "n4_score_class_family",
        "n4_fiber_multiplicity_by_class",
        "n4_quotient_congruence_defect",
        "n4_einheit_torsor_status",
        "n4_deletable_arc_coordinate",
        "n4_minimal_chart_status",
        "n4_S_class_representation_multiplicity",
        "erdos870_representation_abundance_status",
        "erdos870_minimality_sidecar_status",
        "lrc_edge_chart_section_status",
        "n3_edge_flip_kernel",
        "n3_minority_edge_gate",
        "function_quartet_order_status",
        "ordered_function_payload",
        "symmetric_function_shadow",
        "worpitzky_descent_word",
        "compression_circuit_profile",
        "scheme_A_cover_to_section_circuit",
        "transformation_monoid_status",
        "lee_yang_binomial_cap_status",
        "phi4_off_circle_dip_status",
        "large_radius_balance_q0_q6_R6",
        "k8_bimodality_functional",
        "bounded_core_resolvent_degree_ceiling",
        "abel_ruffini_wall_status",
        "ear_omega_sidecar_status",
        "newton_maclaurin_quartic_AP_status",
        "terminal_exit_or_named_debt",
    ]
    print("\n".join(signals))
    print()

    print("SYNTHESIS")
    print("-" * 72)
    print(
        "The fixed-path tiling cube is an abundant cover: S has five representatives "
        "and T is a singleton.  It is not a group quotient because products of score "
        "classes are ambiguous after fibers are forgotten.  The partial-score "
        "0,1,1,2 model with free x=a and y=b is an exact Einheit/Klein-four "
        "chart: T,E is the unit, +=x, -=y, and S=x+y.  Erdős-870 warns that "
        "abundance cannot certify minimality; here c is deletable in the tiling "
        "cover while x and y are essential in the exact chart.  The useful "
        "compression is x=a OR c and y=b OR c: a tiny circuit that preserves "
        "score class but destroys the witness of whether S came from the apex "
        "canary c or from the live pair a,b.  The larger LRC continuation should "
        "therefore pair Lee-Yang/Pascal cap signals, phi4 quartic dips, and "
        "degree-four resolvent claims with explicit section and witness sidecars."
    )


if __name__ == "__main__":
    main()
