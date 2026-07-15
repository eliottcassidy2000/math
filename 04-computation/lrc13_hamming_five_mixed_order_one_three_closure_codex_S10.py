#!/usr/bin/env python3
"""Exact closure of THM-823's mixed D=(1,3,3,3,3) contexts.

The interval engine and THM-815 longest-component recursion are imported from
the independently frozen all-order-three replay for THM-844.  This file owns
the mixed progression parametrization, exhaustive context bank, certificate,
and Tournament Analysis.  All decisions use integers or fractions.Fraction.
"""

from __future__ import annotations

from collections import Counter
from hashlib import sha256
from itertools import product
from pathlib import Path

from lrc13_hamming_five_all_order_three_context_closure_codex_S10 import (
    Audit,
    BASE,
    P,
    QUARTET_BASE,
    component_count,
    context_tournament,
    crt_base,
    global_kl_bound,
    longest_component,
    longest_component_bound,
    measure,
    negative_pairs,
    recurse,
    strict_half_core,
)


PARENT = Path(__file__).with_name(
    "lrc13_hamming_five_all_order_three_context_closure_codex_S10.py"
)
EXPECTED_PARENT_SHA256 = (
    "6ee83b94b23f699939fea5ed53e5c50c7eb3f8a99aeac30375feabf3e0e8befb"
)
EXPECTED_CERTIFICATE = "1cab41e8d32b09d93e9548d1baa486c0e33fb7979695d1b10725829c3a4aeb75"

Context = tuple[tuple[int, ...], int, tuple[int, int]]


def all_contexts() -> tuple[Context, ...]:
    """Return 3 quartet cosets x 8 outside labels x 4 pair words."""

    cosets = sorted(
        {
            tuple(sorted((a * r) % P for r in QUARTET_BASE))
            for a in range(1, P)
        }
    )
    assert len(cosets) == 3
    contexts: list[Context] = []
    for quartet in cosets:
        outside = sorted(set(BASE) - set(quartet))
        assert len(outside) == 8
        for order_one_label in outside:
            for pair_bits in product((1, 2), repeat=2):
                contexts.append((quartet, order_one_label, pair_bits))
    assert len(contexts) == 3 * 8 * 4 == 96
    return tuple(contexts)


def context_data(
    context: Context,
) -> tuple[tuple[int, ...], dict[int, int], tuple[int, ...]]:
    """Build the five progressions and seven-speed scale-three core.

    Quartet label r has effective order three, hence one of the two CRT
    progressions u=3r (mod 13), u=e_r (mod 3).  The order-one label b has
    speed 3v with v=b (mod 13).  The least choice 3b is the deleted AP speed,
    so proper replacement begins at 3b+39 and then advances by 39.
    """

    quartet, order_one_label, pair_bits = context
    pairs = negative_pairs(quartet)
    assert len(pairs) == 2
    units = {
        label: pair_bits[index]
        for index, pair in enumerate(pairs)
        for label in pair
    }
    missing = tuple(sorted((*quartet, order_one_label)))
    bases = {label: crt_base(label, units[label]) for label in quartet}
    bases[order_one_label] = 3 * order_one_label + 39
    core_labels = tuple(label for label in BASE if label not in missing)
    assert len(missing) == 5 and len(core_labels) == 7
    assert all(bases[r] % P == 3 * r % P for r in missing)
    assert bases[order_one_label] % 39 == 3 * order_one_label
    assert bases[order_one_label] > 3 * order_one_label
    return missing, bases, core_labels


def main() -> None:
    parent_hash = sha256(PARENT.read_bytes()).hexdigest()
    assert parent_hash == EXPECTED_PARENT_SHA256

    context_rows: list[object] = []
    aggregate_nodes: Counter[int] = Counter()
    aggregate_dead: Counter[int] = Counter()
    aggregate_states = 0
    aggregate_dominance = 0
    aggregate_strict_dominance = 0
    all_operator_speeds: set[int] = set()
    tournament_rows: list[tuple[object, ...]] = []

    for index, context in enumerate(all_contexts()):
        quartet, order_one_label, pair_bits = context
        missing, bases, core_labels = context_data(context)
        root = strict_half_core(3 * label for label in core_labels)
        assert root and component_count(root) > 0
        root_length = measure(root)
        root_longest = longest_component(root)
        root_global = global_kl_bound(root, root_length, 5)
        root_long = longest_component_bound(root, 5)
        assert root_long <= root_global

        audit = Audit()
        recurse(
            runs=root,
            remaining=missing,
            bases=bases,
            last_speed=0,
            chosen=(),
            audit=audit,
        )
        assert not audit.covers
        assert not audit.terminals

        states = sum(audit.nodes.values())
        aggregate_nodes.update(audit.nodes)
        aggregate_dead.update(audit.dead_no_candidate)
        aggregate_states += states
        aggregate_dominance += audit.dominance_rows
        aggregate_strict_dominance += audit.strict_dominance_rows
        all_operator_speeds |= audit.operator_speeds
        tournament = context_tournament(root, missing, bases)
        tournament_rows.append(tournament)

        row = (
            index,
            quartet,
            order_one_label,
            pair_bits,
            missing,
            tuple(sorted(bases.items())),
            core_labels,
            component_count(root),
            root_length,
            root_longest,
            root_global,
            root_long,
            tuple(sorted(audit.nodes.items())),
            tuple(sorted(audit.dead_no_candidate.items())),
            states,
            audit.max_global_bound,
            audit.max_long_bound,
            len(audit.operator_speeds),
            audit.dominance_rows,
            audit.strict_dominance_rows,
            audit.digest.hexdigest(),
            tournament[2],
        )
        context_rows.append(row)

    certificate = sha256(
        "".join(repr(row) + "\n" for row in context_rows).encode()
    ).hexdigest()
    if EXPECTED_CERTIFICATE != "TO_BE_FILLED":
        assert certificate == EXPECTED_CERTIFICATE

    standard_fingerprint = (
        ((0, 1), (1, 1), (2, 1), (3, 1), (4, 1)),
        0,
        (1, 1, 1, 1, 1),
        1,
    )
    assert all(row[0] == standard_fingerprint for row in tournament_rows)
    assert all(row[1] == standard_fingerprint for row in tournament_rows)
    flip_hist = Counter(row[2] for row in tournament_rows)

    lines = [
        "THM-847 MIXED ORDER-ONE/ORDER-THREE HAMMING-FIVE CLOSURE",
        "arithmetic=integer+fractions.Fraction floating_point=none sampled_circle=none",
        "threshold=1/13 progression_modulus=39 reflection_quotient=(0,1/2)",
        "contexts=3_quartet_cosets*8_outside_D1_labels*4_quartet_pair_words=96",
        "quartet_progression=u=3r_mod13,u=pair_unit_mod3",
        "order_one_progression=u=3b_mod39,proper_base=3b+39",
        "carrier=active_open_endpoint_runs+remaining_labelled_progressions+last_speed",
        "bound.longest=22m/[13(13-2m)Lmax]",
        "bound.global=22mK/[13(13-2m)L]",
        f"aggregate.nodes={dict(sorted(aggregate_nodes.items()))}",
        f"aggregate.dead_no_candidate={dict(sorted(aggregate_dead.items()))}",
        f"aggregate.states={aggregate_states}",
        "aggregate.covering_prefixes=0 aggregate.depth5_terminals=0",
        f"aggregate.dominance_rows={aggregate_dominance}",
        f"aggregate.strict_dominance_rows={aggregate_strict_dominance}",
        f"aggregate.distinct_operator_speeds={len(all_operator_speeds)}",
        f"context.min_states={min(row[14] for row in context_rows)}",
        f"context.max_states={max(row[14] for row in context_rows)}",
        f"context.max_long_bound={max(row[16] for row in context_rows)}",
        f"context.max_global_bound={max(row[15] for row in context_rows)}",
        "CONTEXT_ROWS",
    ]
    for row in context_rows:
        lines.append(
            "context[{:02d}] C={} b={} pair_bits={} bases={} root=(K={},L={},Lmax={},global={},long={}) "
            "nodes={} dead={} states={} max_bounds=({}, {}) operators={} dominance=({},{}) "
            "state_sha256={} tournament_flips={}".format(
                row[0], row[1], row[2], row[3], row[5], row[7], row[8],
                row[9], row[10], row[11], row[12], row[13], row[14],
                row[15], row[16], row[17], row[18], row[19], row[20], row[21]
            )
        )

    lines.extend(
        (
            "TOURNAMENT_ANALYSIS",
            "vertices=five_labelled_least_progression_comb_obligations",
            "pair_observable=root_total_measure_marginal_erosion",
            "switch=root_longest_component_marginal_erosion",
            "tie_hamiltonian_path=increasing_(least_speed,label)",
            f"raw_fingerprint={standard_fingerprint}",
            f"switched_fingerprint={standard_fingerprint}",
            f"edge_flip_histogram={dict(sorted(flip_hist.items()))}",
            f"edge_flips_total={sum(k*v for k,v in flip_hist.items())}",
            "preserves=planning_order_only",
            "destroys=absolute_residual_geometry,component_incidence,future_progressions,cover_truth",
            "assumption_challenge=coset_or_comb_vertices_forget_the_D1_proper_lift_boundary_and_higher_order_overlap",
            "scope=THM823_mixed_D1_D3_contexts_not_other_orders_not_global_H5",
            f"parent_source_sha256={parent_hash}",
            f"certificate_sha256={certificate}",
            f"source_sha256={sha256(Path(__file__).read_bytes()).hexdigest()}",
            "PASS: all 96 arbitrary-height mixed order-one/order-three contexts are strictly loose",
        )
    )
    payload = "\n".join(lines) + "\n"
    print(payload, end="")
    print(f"sha256={sha256(payload.encode()).hexdigest()}")


if __name__ == "__main__":
    main()
