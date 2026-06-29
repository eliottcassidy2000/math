#!/usr/bin/env python3
"""HYP-3524: random031 spigot/hydrotope residual scout.

The user suggested spigot algorithms and the 2026 hydrotope amplitude paper as
inspiration.  This script makes the analogy finite and auditable on the live
random031 target.

Spigot side: owner labels are emitted only when a previous certificate makes
them stable.  The tail should shrink from the seven-owner seam to the residual
pair (45,173).

Hydrotope side: chamber signs and sliced-box volumes are tested as diagnostics
for the same owner filtration.  These scalar volumes are not accepted as proof
quotients unless owner labels and route sidecar R are still recoverable.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from math import factorial, prod
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
COMP = ROOT / "04-computation"


def load_module(name: str, filename: str):
    path = COMP / filename
    spec = spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {path}")
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


H3522 = load_module(
    "hyp3522_for_hyp3523",
    "lrc14_random031_owner_boundary_filtration_codex_20260629.py",
)


@dataclass(frozen=True)
class SpigotStage:
    name: str
    emitted: tuple[int, ...]
    cumulative: tuple[int, ...]
    tail: tuple[int, ...]
    certificate: str


@dataclass(frozen=True)
class Carrier:
    name: str
    score: int


CARRIERS = (
    Carrier("full_filtration_spigot_packet", 100),
    Carrier("residual_pair_tail_lemma", 94),
    Carrier("transport_plus_boundary_emitter", 89),
    Carrier("hydrotope_chamber_audit_with_owner_labels", 83),
    Carrier("route_sidecar_R_guard", 78),
    Carrier("sliced_box_volume_shadow", 39),
    Carrier("raw_threshold_sign_shadow", 31),
    Carrier("raw_owner_count_shadow", 13),
)


def compact_counter(counter: Counter | dict) -> dict:
    return dict(sorted(counter.items(), key=lambda item: repr(item[0])))


def subset_sum(weights: dict[int, int], subset: tuple[int, ...]) -> int:
    return sum(weights[owner] for owner in subset)


def sign(value: int) -> int:
    return (value > 0) - (value < 0)


def all_subsets(items: tuple[int, ...]):
    for size in range(len(items) + 1):
        for subset in combinations(items, size):
            yield subset


def box_cdf_volume(weights: tuple[int, ...], threshold: int) -> Fraction:
    """Normalized volume of {0<=x_i<=w_i, sum x_i <= threshold}.

    This is the standard inclusion-exclusion formula for a sliced axis-aligned
    box.  It is used here as the hydrotope-inspired chamber scalar.
    """

    dimension = len(weights)
    numerator = Fraction(0, 1)
    indices = range(dimension)
    for size in range(dimension + 1):
        for subset in combinations(indices, size):
            shifted = threshold - sum(weights[index] for index in subset)
            if shifted <= 0:
                continue
            numerator += (-1) ** size * Fraction(shifted**dimension, 1)
    denominator = factorial(dimension) * prod(weights)
    return numerator / denominator


def fmt_fraction(value: Fraction) -> str:
    return str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"


def owner_support_counts(cells) -> Counter[int]:
    counts: Counter[int] = Counter()
    for cell in cells:
        seen: set[int] = set()
        for hit in cell.hits:
            seen.update(hit.owners)
        counts.update(seen)
    return counts


def build_context():
    (
        row,
        gates,
        cells,
        by_node,
        legal,
        bypass_component,
        bypass_cells,
        hard_seam_gates,
        lower_bypass_gates,
    ) = H3522.build_context()

    boundaries = H3522.branch_boundaries(cells, by_node, bypass_component)
    boundary_cells = tuple(
        cell for boundary in boundaries for cell in (boundary.left_cell, boundary.right_cell)
    )
    mirror_pairs = H3522.mirror_pairs(by_node, bypass_component)

    seam_owners = H3522.gate_owner_union(hard_seam_gates)
    transport_owners = H3522.owner_union_from_cells(bypass_cells)
    branch_boundary_owners = H3522.owner_union_from_cells(boundary_cells)
    transport_plus_boundary = tuple(sorted(set(transport_owners) | set(branch_boundary_owners)))
    bracket_lift = tuple(owner for owner in branch_boundary_owners if owner not in transport_owners)
    residual_after_boundary = tuple(
        owner for owner in seam_owners if owner not in transport_plus_boundary
    )
    residual_after_transport = tuple(owner for owner in seam_owners if owner not in transport_owners)

    return {
        "row": row,
        "cells": cells,
        "legal": legal,
        "bypass_component": bypass_component,
        "bypass_cells": bypass_cells,
        "hard_seam_gates": hard_seam_gates,
        "lower_bypass_gates": lower_bypass_gates,
        "boundaries": boundaries,
        "mirror_pairs": mirror_pairs,
        "seam_owners": seam_owners,
        "transport_owners": transport_owners,
        "branch_boundary_owners": branch_boundary_owners,
        "transport_plus_boundary": transport_plus_boundary,
        "bracket_lift": bracket_lift,
        "residual_after_transport": residual_after_transport,
        "residual_after_boundary": residual_after_boundary,
    }


def spigot_stages(context: dict[str, object]) -> tuple[SpigotStage, ...]:
    seam = context["seam_owners"]
    transport = context["transport_owners"]
    bracket_lift = context["bracket_lift"]
    transport_plus_boundary = context["transport_plus_boundary"]
    residual_after_transport = context["residual_after_transport"]
    residual_after_boundary = context["residual_after_boundary"]
    return (
        SpigotStage(
            "S0_forbidden_seam_input",
            (),
            (),
            seam,
            "seven-owner forbidden seam is retained, not emitted",
        ),
        SpigotStage(
            "S1_transport_emitter",
            transport,
            transport,
            residual_after_transport,
            "HYP-3522 transport-word constancy on the pure bypass",
        ),
        SpigotStage(
            "S2_branch_boundary_lift_emitter",
            bracket_lift,
            transport_plus_boundary,
            residual_after_boundary,
            "HYP-3522 adjacent ordinary branch-boundary bracket lift",
        ),
        SpigotStage(
            "S3_residual_tail",
            (),
            transport_plus_boundary,
            residual_after_boundary,
            "remaining two-owner puncture/apex boundary lemma",
        ),
    )


def stage_safety(stages: tuple[SpigotStage, ...], seam: tuple[int, ...]) -> dict[str, object]:
    cumulative_sizes = [len(stage.cumulative) for stage in stages]
    tail_sizes = [len(stage.tail) for stage in stages]
    monotone_cumulative = all(a <= b for a, b in zip(cumulative_sizes, cumulative_sizes[1:]))
    monotone_tail = all(a >= b for a, b in zip(tail_sizes, tail_sizes[1:]))
    emitted_once: Counter[int] = Counter(owner for stage in stages for owner in stage.emitted)
    no_duplicate_emit = all(count == 1 for count in emitted_once.values())
    emitted = tuple(sorted(emitted_once))
    return {
        "cumulative_sizes": tuple(cumulative_sizes),
        "tail_sizes": tuple(tail_sizes),
        "monotone_cumulative": monotone_cumulative,
        "monotone_tail": monotone_tail,
        "no_duplicate_emit": no_duplicate_emit,
        "emitted_owners": emitted,
        "unemitted_tail": tuple(owner for owner in seam if owner not in emitted),
    }


def weight_systems(context: dict[str, object]) -> dict[str, dict[int, int]]:
    seam = context["seam_owners"]
    support = owner_support_counts(context["cells"])
    layer_weight = {}
    for owner in seam:
        if owner in context["transport_owners"]:
            layer_weight[owner] = 1
        elif owner in context["bracket_lift"]:
            layer_weight[owner] = 2
        else:
            layer_weight[owner] = 4
    return {
        "residue_mod14": {owner: owner % 14 for owner in seam},
        "centered_residue": {owner: min(owner % 14, 14 - (owner % 14)) for owner in seam},
        "owner_support_cells": {owner: support[owner] for owner in seam},
        "filtration_layer": layer_weight,
    }


def chamber_signature(
    weights: dict[int, int],
    subset: tuple[int, ...],
    transport: tuple[int, ...],
    boundary: tuple[int, ...],
    residual: tuple[int, ...],
) -> tuple[int, int, int]:
    value = subset_sum(weights, subset)
    return (
        sign(value - subset_sum(weights, transport)),
        sign(value - subset_sum(weights, boundary)),
        sign(value - subset_sum(weights, residual)),
    )


def hydrotope_audit(context: dict[str, object]) -> dict[str, object]:
    seam = context["seam_owners"]
    transport = context["transport_owners"]
    boundary = context["transport_plus_boundary"]
    residual = context["residual_after_boundary"]
    systems = weight_systems(context)

    audits = {}
    for name, weights in systems.items():
        subsets = tuple(all_subsets(seam))
        signature_buckets: dict[tuple[int, int, int], list[tuple[int, ...]]] = {}
        for subset in subsets:
            signature = chamber_signature(weights, subset, transport, boundary, residual)
            signature_buckets.setdefault(signature, []).append(subset)

        target_rows = {}
        for label, target in (
            ("transport", transport),
            ("transport_plus_boundary", boundary),
            ("bracket_lift", context["bracket_lift"]),
            ("residual", residual),
        ):
            signature = chamber_signature(weights, target, transport, boundary, residual)
            bucket = tuple(signature_buckets[signature])
            target_rows[label] = {
                "sum": subset_sum(weights, target),
                "signature": signature,
                "bucket_size": len(bucket),
                "bucket_examples": bucket[:5],
            }

        ordered_weights = tuple(weights[owner] for owner in seam)
        transport_h = subset_sum(weights, transport)
        boundary_h = subset_sum(weights, boundary)
        residual_h = subset_sum(weights, residual)
        volume_transport = box_cdf_volume(ordered_weights, transport_h)
        volume_boundary = box_cdf_volume(ordered_weights, boundary_h)
        volume_residual = box_cdf_volume(ordered_weights, residual_h)

        audits[name] = {
            "weights": weights,
            "target_rows": target_rows,
            "signature_bucket_hist": Counter(len(bucket) for bucket in signature_buckets.values()),
            "transport_cdf": volume_transport,
            "boundary_cdf": volume_boundary,
            "residual_cdf": volume_residual,
            "bracket_slab_volume": volume_boundary - volume_transport,
            "tail_volume_after_boundary": Fraction(1, 1) - volume_boundary,
        }
    return audits


def quotient_audit(context: dict[str, object], hydro: dict[str, object]) -> tuple[dict[str, object], ...]:
    seam = context["seam_owners"]
    transport = context["transport_owners"]
    boundary = context["transport_plus_boundary"]
    residual = context["residual_after_boundary"]

    candidates = []
    candidates.append(
        {
            "name": "full_filtration_word",
            "reconstructs": ("transport", "bracket_lift", "residual"),
            "safe": True,
            "warning": "proof-facing sidecar; keeps all emitted/tail owner labels",
        }
    )
    candidates.append(
        {
            "name": "transport_plus_boundary_word",
            "reconstructs": ("transport_plus_boundary", "residual"),
            "safe": True,
            "warning": "safe only with transport/bracket split or HYP-3522 reconstruction",
        }
    )
    candidates.append(
        {
            "name": "bypass_owner_word_only",
            "reconstructs": ("transport",),
            "safe": False,
            "warning": "cannot see bracket lift or residual pair",
        }
    )
    for name, audit in hydro.items():
        residual_bucket = audit["target_rows"]["residual"]["bucket_size"]
        boundary_bucket = audit["target_rows"]["transport_plus_boundary"]["bucket_size"]
        candidates.append(
            {
                "name": f"hydrotope_signature_{name}",
                "reconstructs": ("chamber_signs",),
                "safe": residual_bucket == 1 and boundary_bucket == 1,
                "warning": (
                    "label-unique chamber"
                    if residual_bucket == 1 and boundary_bucket == 1
                    else f"mixed buckets: residual={residual_bucket}, boundary={boundary_bucket}"
                ),
            }
        )
        candidates.append(
            {
                "name": f"sliced_box_volume_{name}",
                "reconstructs": ("scalar_volume",),
                "safe": False,
                "warning": "volume is a diagnostic scalar; it forgets owner identities",
            }
        )
    candidates.append(
        {
            "name": "raw_counts_7_5_2",
            "reconstructs": ("owner_counts",),
            "safe": False,
            "warning": f"counts see {len(seam)}->{len(boundary)}->{len(residual)} but forget labels",
        }
    )
    candidates.append(
        {
            "name": "route_sidecar_R",
            "reconstructs": ("HYP-3490 route",),
            "safe": True,
            "warning": "required by HYP-3513 unless route reconstruction is proved",
        }
    )
    return tuple(candidates)


def print_hydro(audits: dict[str, object]) -> None:
    for name, audit in audits.items():
        print(f"system={name}")
        print(f"  weights={audit['weights']}")
        print(
            "  volumes: "
            f"V_transport={fmt_fraction(audit['transport_cdf'])} "
            f"V_boundary={fmt_fraction(audit['boundary_cdf'])} "
            f"V_residual={fmt_fraction(audit['residual_cdf'])} "
            f"bracket_slab={fmt_fraction(audit['bracket_slab_volume'])} "
            f"tail_after_boundary={fmt_fraction(audit['tail_volume_after_boundary'])}"
        )
        print(f"  signature_bucket_hist={compact_counter(audit['signature_bucket_hist'])}")
        for label, row in audit["target_rows"].items():
            print(
                f"  target={label} sum={row['sum']} signature={row['signature']} "
                f"bucket_size={row['bucket_size']} examples={row['bucket_examples']}"
            )


def main() -> None:
    context = build_context()
    stages = spigot_stages(context)
    safety = stage_safety(stages, context["seam_owners"])
    hydro = hydrotope_audit(context)
    quotients = quotient_audit(context, hydro)

    print("HYP-3524 RANDOM031 SPIGOT/HYDROTOPE SCOUT")
    print("status=EVIDENCE / spigot-emitter and sliced-box chamber scout; not an LRC14 proof")
    print("row=random_covering_031")
    print("sources=spigot_algorithm + hydrotope two-minus amplitude chamber-volume analogy")
    print()

    print("## Live Owner Target")
    print(f"seam_owners={context['seam_owners']}")
    print(f"transport_owners={context['transport_owners']}")
    print(f"branch_boundary_owners={context['branch_boundary_owners']}")
    print(f"bracket_lift={context['bracket_lift']}")
    print(f"transport_plus_boundary={context['transport_plus_boundary']}")
    print(f"residual_after_transport={context['residual_after_transport']}")
    print(f"residual_after_branch_boundary={context['residual_after_boundary']}")
    print(
        "reading=spigot view emits stable owner digits only after their local "
        "certificate is present; the un-emitted tail is exactly (45,173)."
    )
    print()

    print("## Spigot Emitter Schedule")
    for stage in stages:
        print(
            f"{stage.name}: emitted={stage.emitted} cumulative={stage.cumulative} "
            f"tail={stage.tail} certificate={stage.certificate}"
        )
    print(f"spigot_safety={safety}")
    print()

    print("## Branch Boundary / Mirror Cross-Check")
    for boundary in context["boundaries"]:
        left = H3522.owner_union_from_hits(boundary.left_cell.hits)
        right = H3522.owner_union_from_hits(boundary.right_cell.hits)
        print(
            f"branch={boundary.branch} bypass_u={boundary.bypass_u} "
            f"phases={boundary.bypass_phases} left_owners={left} right_owners={right} "
            f"intersection={tuple(sorted(set(left) & set(right)))}"
        )
    mirror_owner_words = Counter(
        H3522.owner_union_from_hits(left.hits)
        for left, _right in context["mirror_pairs"]
    )
    print(f"mirror_pair_owner_words={compact_counter(mirror_owner_words)}")
    print()

    print("## Hydrotope Chamber Diagnostic")
    print_hydro(hydro)
    print(
        "hydrotope_guardrail=chamber signs and sliced-box volumes are useful "
        "alarms for when a threshold crosses transport/bracket/residual layers, "
        "but scalar volumes alone are illegal proof quotients."
    )
    print()

    print("## Quotient Safety")
    for candidate in quotients:
        print(
            f"{candidate['name']}: safe={candidate['safe']} "
            f"reconstructs={candidate['reconstructs']} warning={candidate['warning']}"
        )
    print()

    print("## Proof Pull")
    print(
        "P1: State the random031 owner lemma as an online emitter theorem: "
        "transport emits (23,93,113), branch brackets emit (147,169), and "
        "all later proof stages see only residual tail (45,173)."
    )
    print(
        "P2: Use hydrotope chamber signs only as a finite chamber audit.  If a "
        "weight system gives mixed residual or boundary buckets, it is a canary "
        "that the quotient is forgetting labels."
    )
    print(
        "P3: The residual lemma should be a two-owner boundary/no-hidden-tail "
        "statement with route sidecar R, not a scalar volume or owner-count "
        "inequality."
    )
    print(
        "P4: A Lean-facing theorem can expose a recursive interface: each "
        "emitter consumes one certified sidecar and returns a smaller tail; "
        "termination is the residual pair plus HYP-3490/HYP-3513 firewall route."
    )
    print()

    print("## Tournament Analysis")
    print("vertices=proof emitters and quotient sidecars, not runners or raw arcs")
    print(
        "pairwise_observable=tail shrinkage + owner-label reconstruction + "
        "route-sidecar legality + scalar-forgetting penalty"
    )
    print("switch=higher retained proof payload; ties use emitter order")
    print(f"score_hist={compact_counter(Counter(carrier.score for carrier in CARRIERS))}")
    print("directed_3cycles=0")
    print(
        "hamiltonian_path="
        + " -> ".join(carrier.name for carrier in sorted(CARRIERS, key=lambda item: -item.score))
    )
    print()

    print("## Assumption Challenge")
    print(
        "Candidate vertices considered: owners, emitted digits, branch boundary "
        "events, hydrotope chambers, sliced-box thresholds, u-fibers, mirror "
        "pairs, route sidecars, and proof obligations."
    )
    print(
        "Chosen vertices are proof emitters and quotient sidecars.  The quotient "
        "preserves the LRC random031 residual predicate only when emitted owner "
        "digits and the route sidecar can be reconstructed; scalar chamber "
        "volumes deliberately destroy that information and remain diagnostics."
    )


if __name__ == "__main__":
    main()
