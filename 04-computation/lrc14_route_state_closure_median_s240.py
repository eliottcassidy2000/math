#!/usr/bin/env python3
"""Route-state closure median scout for LRC14 proof carriers.

This is not a proof of LRC14.  It is a small exact interface test for the
current proof-stack idea:

    serious route triples should have a unique legal center after sidecars
    are attached.

The state space is intentionally finite and declarative.  A route/certificate
state is a set of named packet, route, certificate, sidecar, and discharge
coordinates.  Legal sidecars close a state under local rules.  The median of a
triple is the coordinate-wise majority state.  A triple passes only when the
closed median is legal; otherwise the scout reports the missing sidecars or
explicit debt exits.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from typing import Dict, FrozenSet, Iterable, List, Sequence, Set, Tuple


BASE = frozenset(
    {
        "lrc_predicate",
        "packet_scale_exact",
        "endpoint_owner",
        "topology_safe",
    }
)

FIELDS = [
    "lrc_predicate",
    "packet_scale_exact",
    "endpoint_owner",
    "topology_safe",
    "q_witness",
    "ap_gw_boundary_atom",
    "endpoint_owner_cycle",
    "route_word_automatic",
    "native_transition_word",
    "magnitude_cocycle",
    "fibbinary_partial_cube",
    "partial_cube_theta",
    "theta_class_word",
    "gated_subcube_status",
    "median_interval_status",
    "moser_two_lane_split",
    "doubled_triangular_bridge",
    "simplex_face_rank",
    "bridge_rank_line_id",
    "toeplitz_square_scale_gate",
    "ordered_quad_collapse_mode",
    "positive_scale_noncollapse",
    "fejer_interval_certificate",
    "haar_zeta_square",
    "observer_cut_payload",
    "cross_sector_orientation",
    "deletion_fiber_profile",
    "rectangle_hourglass_residue",
    "certificate_cycle_generators",
    "cycle_generator_checked",
    "cycle_class_image_status",
    "hodge_positive_shadow",
    "cochain_closedness_status",
    "hodge_cycle_image",
    "relation_lattice_covolume",
    "low_height_wall_ledger",
    "exceptional_approximant",
    "residual_tail_signed",
    "residual_capacitor_cut",
    "residual_tooth",
    "state_lift_obligation",
    "f7_state_lift_target",
    "residual_debt_exit",
    "discharge_resolved",
    "discharge_named_debt",
]


@dataclass(frozen=True)
class State:
    name: str
    fields: FrozenSet[str]


def fs(*items: str) -> FrozenSet[str]:
    return frozenset(items)


STATES: Dict[str, State] = {
    "q_witness_exact_scale": State(
        "q_witness_exact_scale",
        BASE
        | fs(
            "q_witness",
            "magnitude_cocycle",
            "discharge_resolved",
        ),
    ),
    "ap_gw_boundary_atom": State(
        "ap_gw_boundary_atom",
        BASE
        | fs(
            "ap_gw_boundary_atom",
            "endpoint_owner_cycle",
            "discharge_resolved",
        ),
    ),
    "automatic_word_shadow": State(
        "automatic_word_shadow",
        BASE
        | fs(
            "route_word_automatic",
            "fibbinary_partial_cube",
            "moser_two_lane_split",
        ),
    ),
    "partial_cube_bridge_rank": State(
        "partial_cube_bridge_rank",
        BASE
        | fs(
            "fibbinary_partial_cube",
            "doubled_triangular_bridge",
            "partial_cube_theta",
        ),
    ),
    "fejer_interval_certificate": State(
        "fejer_interval_certificate",
        BASE
        | fs(
            "fejer_interval_certificate",
            "hodge_positive_shadow",
            "certificate_cycle_generators",
        ),
    ),
    "toeplitz_square_scale": State(
        "toeplitz_square_scale",
        BASE
        | fs(
            "toeplitz_square_scale_gate",
            "ordered_quad_collapse_mode",
            "positive_scale_noncollapse",
            "hodge_positive_shadow",
            "certificate_cycle_generators",
        ),
    ),
    "hodge_positive_shadow": State(
        "hodge_positive_shadow",
        BASE
        | fs(
            "hodge_positive_shadow",
            "cochain_closedness_status",
        ),
    ),
    "fejer_interval_cycle": State(
        "fejer_interval_cycle",
        BASE
        | fs(
            "fejer_interval_certificate",
            "hodge_positive_shadow",
            "certificate_cycle_generators",
            "cycle_generator_checked",
        ),
    ),
    "toeplitz_square_cycle": State(
        "toeplitz_square_cycle",
        BASE
        | fs(
            "toeplitz_square_scale_gate",
            "ordered_quad_collapse_mode",
            "positive_scale_noncollapse",
            "hodge_positive_shadow",
            "certificate_cycle_generators",
            "cycle_generator_checked",
        ),
    ),
    "hodge_generated_cycle": State(
        "hodge_generated_cycle",
        BASE
        | fs(
            "hodge_positive_shadow",
            "cochain_closedness_status",
            "certificate_cycle_generators",
            "cycle_generator_checked",
            "hodge_cycle_image",
        ),
    ),
    "haar_zipper_square": State(
        "haar_zipper_square",
        BASE
        | fs(
            "haar_zeta_square",
            "observer_cut_payload",
            "certificate_cycle_generators",
            "cycle_generator_checked",
        ),
    ),
    "observer_cut_sector": State(
        "observer_cut_sector",
        BASE
        | fs(
            "observer_cut_payload",
            "cross_sector_orientation",
            "deletion_fiber_profile",
            "rectangle_hourglass_residue",
        ),
    ),
    "residual_capacitor_cut": State(
        "residual_capacitor_cut",
        BASE
        | fs(
            "observer_cut_payload",
            "residual_capacitor_cut",
            "residual_tooth",
            "residual_debt_exit",
        ),
    ),
    "k33_state_lift": State(
        "k33_state_lift",
        BASE
        | fs(
            "state_lift_obligation",
            "bridge_rank_line_id",
            "residual_debt_exit",
            "f7_state_lift_target",
        ),
    ),
    "f7_residual_debt": State(
        "f7_residual_debt",
        BASE
        | fs(
            "residual_debt_exit",
            "f7_state_lift_target",
            "discharge_named_debt",
        ),
    ),
    "hodge_phantom_debt_state": State(
        "hodge_phantom_debt_state",
        BASE
        | fs(
            "hodge_positive_shadow",
            "cochain_closedness_status",
            "residual_debt_exit",
            "discharge_named_debt",
        ),
    ),
    "roth_minkowski_wall": State(
        "roth_minkowski_wall",
        BASE
        | fs(
            "relation_lattice_covolume",
            "low_height_wall_ledger",
            "exceptional_approximant",
            "residual_tail_signed",
        ),
    ),
}


RULES: Sequence[Tuple[FrozenSet[str], FrozenSet[str], str]] = [
    (
        fs("route_word_automatic"),
        fs("native_transition_word", "magnitude_cocycle"),
        "automatic words must retain native transition and magnitude data",
    ),
    (
        fs("fibbinary_partial_cube"),
        fs("theta_class_word", "gated_subcube_status", "median_interval_status"),
        "fibbinary/Moser data must be a gated partial-cube sidecar",
    ),
    (
        fs("moser_two_lane_split"),
        fs("native_transition_word", "magnitude_cocycle"),
        "Moser split is not usable without the sequence transition",
    ),
    (
        fs("doubled_triangular_bridge"),
        fs("bridge_rank_line_id", "simplex_face_rank"),
        "doubled triangular bridge must name its K_{k,k+1} line",
    ),
    (
        fs("toeplitz_square_scale_gate"),
        fs("ordered_quad_collapse_mode", "positive_scale_noncollapse"),
        "Toeplitz square-peg carrier must keep the noncollapse scale",
    ),
    (
        fs("fejer_interval_certificate"),
        fs("certificate_cycle_generators", "packet_scale_exact"),
        "Fejer interval certificate must name its cycle generators",
    ),
    (
        fs("hodge_positive_shadow"),
        fs("cochain_closedness_status"),
        "positive Hodge shadow must first be a closed packet cochain",
    ),
    (
        fs("cycle_generator_checked", "certificate_cycle_generators", "hodge_positive_shadow"),
        fs("cycle_class_image_status", "hodge_cycle_image"),
        "checked cycle generators discharge the positive shadow",
    ),
    (
        fs("haar_zeta_square", "cycle_generator_checked"),
        fs("cycle_class_image_status"),
        "Haar zipper squares must report their cycle-class image",
    ),
    (
        fs("observer_cut_payload"),
        fs("packet_scale_exact", "endpoint_owner"),
        "observer-cut payloads must retain packet scale and endpoint owner",
    ),
    (
        fs("state_lift_obligation"),
        fs("f7_state_lift_target", "residual_debt_exit"),
        "state-lift obligations must name the forbidden/debt target",
    ),
]


def close_fields(seed: FrozenSet[str]) -> FrozenSet[str]:
    """Attach all legal forced sidecars."""
    fields = set(seed)
    changed = True
    while changed:
        changed = False
        for antecedent, consequent, _reason in RULES:
            if antecedent <= fields and not consequent <= fields:
                fields.update(consequent)
                changed = True
    return frozenset(fields)


def missing_legalities(fields: FrozenSet[str]) -> List[str]:
    missing: List[str] = []
    for base in sorted(BASE):
        if base not in fields:
            missing.append(f"missing base coordinate `{base}`")
    for antecedent, consequent, reason in RULES:
        if antecedent <= fields and not consequent <= fields:
            lack = ", ".join(sorted(consequent - fields))
            missing.append(f"{reason}: missing {lack}")
    if (
        "hodge_positive_shadow" in fields
        and "hodge_cycle_image" not in fields
        and "residual_debt_exit" not in fields
    ):
        missing.append(
            "positive closed shadow has no cycle-class image and no named residual debt exit"
        )
    if (
        "observer_cut_payload" in fields
        and "cross_sector_orientation" not in fields
        and "cycle_class_image_status" not in fields
        and "residual_debt_exit" not in fields
    ):
        missing.append(
            "observer-cut payload has no cross-sector, cycle-image, or debt sidecar"
        )
    return missing


def is_legal(fields: FrozenSet[str]) -> bool:
    return not missing_legalities(fields)


def median(a: FrozenSet[str], b: FrozenSet[str], c: FrozenSet[str]) -> FrozenSet[str]:
    center = {x for x in FIELDS if (x in a) + (x in b) + (x in c) >= 2}
    return frozenset(center)


def format_fields(fields: Iterable[str]) -> str:
    return ", ".join(sorted(fields)) or "(none)"


TRIPLES = [
    (
        "automatic_partial_cube_router",
        (
            "automatic_word_shadow",
            "partial_cube_bridge_rank",
            "q_witness_exact_scale",
        ),
        "Moser/fibbinary plus bridge rank must become a gated partial cube, not a raw automaton quotient.",
    ),
    (
        "hodge_toeplitz_fejer_phantom",
        (
            "hodge_positive_shadow",
            "fejer_interval_certificate",
            "toeplitz_square_scale",
        ),
        "Positivity plus two analytic certificates is still phantom unless cycle image or debt is named.",
    ),
    (
        "hodge_toeplitz_fejer_repaired",
        (
            "hodge_generated_cycle",
            "fejer_interval_cycle",
            "toeplitz_square_cycle",
        ),
        "Checked cycle generators should medianize the Hodge/Toeplitz/Fejer route.",
    ),
    (
        "named_residual_debt_exit",
        (
            "hodge_phantom_debt_state",
            "k33_state_lift",
            "f7_residual_debt",
        ),
        "A phantom positive shadow can pass only by becoming an explicit F7/THM-572 debt exit.",
    ),
    (
        "observer_cut_incompatible_repairs",
        (
            "observer_cut_sector",
            "haar_zipper_square",
            "residual_capacitor_cut",
        ),
        "Observer-cut payload alone is not enough when its three repair lanes do not agree.",
    ),
]


def evaluate_triple(names: Sequence[str]) -> Dict[str, object]:
    raw = [STATES[name].fields for name in names]
    closed = [close_fields(fields) for fields in raw]
    raw_center = median(*raw)
    closed_center = median(*closed)
    raw_missing = missing_legalities(raw_center)
    closed_missing = missing_legalities(closed_center)
    return {
        "raw_center": raw_center,
        "closed_center": closed_center,
        "raw_missing": raw_missing,
        "closed_missing": closed_missing,
        "raw_pass": not raw_missing,
        "closed_pass": not closed_missing,
        "added_sidecars": [
            (name, closed_fields - STATES[name].fields)
            for name, closed_fields in zip(names, closed)
            if closed_fields - STATES[name].fields
        ],
    }


TOURNAMENT_TIE_PATH = [
    "hodge_generated_cycle",
    "toeplitz_square_cycle",
    "fejer_interval_cycle",
    "partial_cube_bridge_rank",
    "automatic_word_shadow",
    "roth_minkowski_wall",
    "observer_cut_sector",
    "haar_zipper_square",
    "residual_capacitor_cut",
    "k33_state_lift",
    "hodge_phantom_debt_state",
    "f7_residual_debt",
    "ap_gw_boundary_atom",
    "q_witness_exact_scale",
    "toeplitz_square_scale",
    "fejer_interval_certificate",
    "hodge_positive_shadow",
]

WEIGHTS = {field: 1 for field in FIELDS}
WEIGHTS.update(
    {
        "hodge_cycle_image": 4,
        "cycle_class_image_status": 4,
        "residual_debt_exit": 4,
        "bridge_rank_line_id": 3,
        "gated_subcube_status": 3,
        "median_interval_status": 3,
        "observer_cut_payload": 2,
        "cross_sector_orientation": 2,
        "toeplitz_square_scale_gate": 2,
        "fejer_interval_certificate": 2,
        "discharge_resolved": 2,
        "discharge_named_debt": 2,
    }
)


def retention_score(fields: FrozenSet[str]) -> int:
    legal_bonus = 10 if is_legal(fields) else 0
    return legal_bonus + sum(WEIGHTS.get(field, 1) for field in fields)


def orient(a: str, b: str, closed: bool) -> Tuple[str, str]:
    fa = close_fields(STATES[a].fields) if closed else STATES[a].fields
    fb = close_fields(STATES[b].fields) if closed else STATES[b].fields
    sa = retention_score(fa)
    sb = retention_score(fb)
    if sa > sb:
        return a, b
    if sb > sa:
        return b, a
    ia = TOURNAMENT_TIE_PATH.index(a)
    ib = TOURNAMENT_TIE_PATH.index(b)
    return (a, b) if ia < ib else (b, a)


def tournament(closed: bool) -> Dict[str, object]:
    names = list(STATES)
    edges: Dict[Tuple[str, str], str] = {}
    outdegree = {name: 0 for name in names}
    for a, b in combinations(names, 2):
        winner, loser = orient(a, b, closed)
        edges[(winner, loser)] = winner
        outdegree[winner] += 1
    cycles = []
    for a, b, c in combinations(names, 3):
        pair_winners = {
            frozenset((a, b)): orient(a, b, closed)[0],
            frozenset((a, c)): orient(a, c, closed)[0],
            frozenset((b, c)): orient(b, c, closed)[0],
        }
        if (
            pair_winners[frozenset((a, b))] == a
            and pair_winners[frozenset((b, c))] == b
            and pair_winners[frozenset((a, c))] == c
        ) or (
            pair_winners[frozenset((a, b))] == b
            and pair_winners[frozenset((b, c))] == c
            and pair_winners[frozenset((a, c))] == a
        ):
            cycles.append((a, b, c))
    score_hist: Dict[int, int] = {}
    for value in outdegree.values():
        score_hist[value] = score_hist.get(value, 0) + 1
    order = sorted(names, key=lambda name: (-outdegree[name], TOURNAMENT_TIE_PATH.index(name)))
    return {
        "edges": edges,
        "outdegree": outdegree,
        "score_hist": dict(sorted(score_hist.items())),
        "directed_3cycles": cycles,
        "hamiltonian_path_count": count_hamiltonian_paths(names, closed),
        "tie_hamiltonian_path": order,
    }


def count_hamiltonian_paths(names: Sequence[str], closed: bool) -> int:
    n = len(names)
    idx = {name: i for i, name in enumerate(names)}
    beats = [[False] * n for _ in range(n)]
    for a, b in combinations(names, 2):
        winner, loser = orient(a, b, closed)
        beats[idx[winner]][idx[loser]] = True
    dp: Dict[Tuple[int, int], int] = {}
    for i in range(n):
        dp[(1 << i, i)] = 1
    for mask in range(1 << n):
        for last in range(n):
            count = dp.get((mask, last), 0)
            if count == 0:
                continue
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                if beats[last][nxt]:
                    dp[(mask | (1 << nxt), nxt)] = dp.get((mask | (1 << nxt), nxt), 0) + count
    full = (1 << n) - 1
    return sum(dp.get((full, last), 0) for last in range(n))


def edge_flips() -> List[Tuple[str, str, str, str]]:
    flips = []
    for a, b in combinations(STATES, 2):
        raw = orient(a, b, closed=False)
        closed = orient(a, b, closed=True)
        if raw != closed:
            flips.append((a, b, raw[0], closed[0]))
    return flips


def main() -> None:
    print("LRC14 route-state closure median scout (codex S240)")
    print("=" * 58)
    print()
    print("State coordinates:", len(FIELDS))
    print("States:", len(STATES))
    print("Closure rules:", len(RULES))
    print()
    print("Assumption challenge")
    print(
        "  Vertices are proof states, not runners: packet states, route obligations, "
        "certificate generators, sidecar coordinates, and discharge exits."
    )
    print(
        "  Preserved predicate: whether a serious route triple has a unique legal "
        "majority center after sidecars are attached."
    )
    print(
        "  Destroyed data: raw runner names, raw time geometry, and the distinction "
        "between two witnesses that agree on all retained proof coordinates."
    )
    print(
        "  Challenged assumption: discharge mode is not a label on a route; it must "
        "be a coordinate, otherwise median centers can look legal while hiding F7/THM-572 debt."
    )
    print()

    print("Route triples")
    print("-" * 58)
    for label, names, note in TRIPLES:
        result = evaluate_triple(names)
        print(label)
        print(f"  states: {', '.join(names)}")
        print(f"  note: {note}")
        print(
            f"  raw median: {'PASS' if result['raw_pass'] else 'FAIL'}; "
            f"{format_fields(result['raw_center'])}"
        )
        if result["raw_missing"]:
            for missing in result["raw_missing"]:
                print(f"    raw missing: {missing}")
        print(
            f"  closed median: {'PASS' if result['closed_pass'] else 'FAIL'}; "
            f"{format_fields(result['closed_center'])}"
        )
        if result["closed_missing"]:
            for missing in result["closed_missing"]:
                print(f"    closed missing: {missing}")
        if result["added_sidecars"]:
            for name, added in result["added_sidecars"]:
                print(f"    sidecars attached to {name}: {format_fields(added)}")
        print()

    raw_t = tournament(closed=False)
    closed_t = tournament(closed=True)
    flips = edge_flips()

    print("Tournament Analysis")
    print("-" * 58)
    print(
        "  Pairwise observable: weighted retention score of proof-state "
        "coordinates, with a legal-center bonus."
    )
    print(
        "  Binary relation: A -> B when A retains the higher score; ties follow "
        "the declared Hamiltonian path."
    )
    print("  Raw score histogram:", raw_t["score_hist"])
    print("  Closed score histogram:", closed_t["score_hist"])
    print("  Raw directed 3-cycles:", len(raw_t["directed_3cycles"]))
    print("  Closed directed 3-cycles:", len(closed_t["directed_3cycles"]))
    print("  Raw Hamiltonian path count:", raw_t["hamiltonian_path_count"])
    print("  Closed Hamiltonian path count:", closed_t["hamiltonian_path_count"])
    print("  Edge flips after legal sidecars:", len(flips))
    for a, b, raw_winner, closed_winner in flips[:12]:
        print(f"    flip {a} vs {b}: raw {raw_winner}, closed {closed_winner}")
    if len(flips) > 12:
        print(f"    ... {len(flips) - 12} more flips")
    print("  Closed tie Hamiltonian path:")
    print("    " + " > ".join(closed_t["tie_hamiltonian_path"]))
    print()

    print("Proof-interface readout")
    print("-" * 58)
    print(
        "  Angle 1, automatic/partial-cube bridge: the raw automaton/Moser/fibbinary "
        "triple fails because the median remembers only `fibbinary_partial_cube`. "
        "After sidecars, the center becomes the gated partial-cube state with "
        "`theta_class_word`, `gated_subcube_status`, `median_interval_status`, "
        "and `magnitude_cocycle`; the non-majority witnesses still carry the "
        "native-transition and bridge-rank sidecars needed to audit the closure."
    )
    print(
        "  Angle 2, Hodge/Toeplitz/Fejer certificates: positivity and analytic "
        "certificates do not medianize unless the center has either "
        "`hodge_cycle_image` or `residual_debt_exit`.  The repaired triple passes "
        "only because checked cycle generators force `cycle_class_image_status`."
    )
    print(
        "  Interface claim to test next: every serious LRC14 route triple should be "
        "representable in this coordinate cube, and every failed median should name "
        "one of: missing gated sidecar, missing cycle image, missing observer-cut "
        "repair, or explicit F7/THM-572 debt."
    )


if __name__ == "__main__":
    main()
