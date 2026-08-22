#!/usr/bin/env python3
"""Exact q=15 boundary probe for the first effective three-phase mode.

This is outside the literal LRC(14) body/divisor atlas: no speed in 1..14 is
q-divisible when q=15.  It tests the first local mode size not present at
q<=14 and separates local availability from global effectivity through rank
six.  Higher-rank physical edges are deliberately outside this probe and are
not inferred from its truncated edge-generated profile.

Runtime gates survive python -O.
"""

from __future__ import annotations

import ast
from collections import Counter
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from math import gcd
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MODE_PATH = ROOT / "04-computation/lrc14_q8_q14_finite_mode_clutter_probe_20260814.py"
PINNED = (
    (
        "q8-q14-script",
        MODE_PATH,
        "e4afeeef7ea6dc2b9b64d3085390b0b6cc53a461dd9bde817c473688dca136f3",
    ),
    (
        "q8-q14-output",
        ROOT / "05-knowledge/results/lrc14_q8_q14_finite_mode_clutter_probe_20260814.out",
        "2fdede2a04ad2b9ccdfbae8b597ecbf82ad71274fe6ba0c221abd7635aed140f",
    ),
)

Q = 15
VERTICES = tuple(range(1, 15))
SEARCH_MAX_RANK = 6
CANONICAL_EDGE = (1, 2, 3, 4, 5, 7)
DOMINO_SUFFICIENT_EDGES = (
    (1, 6, 7, 10, 11, 13),
    (1, 7, 10, 11, 12, 13),
)
EXPECTED_SEMANTIC_DIGEST = "bae0a83a2ca7f19906309b1bab951bd9e4dac241e3a3103726da3bb8099e1eb8"


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def lf_hash(path):
    return sha256(path.read_bytes().replace(bytes((13, 10)), bytes((10,)))).hexdigest()


class ExactDigest:
    def __init__(self):
        self._hash = sha256()

    def add(self, value):
        self._hash.update(repr(value).encode("ascii"))
        self._hash.update(bytes((10,)))

    def hexdigest(self):
        return self._hash.hexdigest()


def load_modes():
    spec = spec_from_file_location("finite_modes", MODE_PATH)
    require(spec is not None and spec.loader is not None, "mode module spec")
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def restricted_mode_search(module, subset, max_phase_size):
    def bank(speed):
        return tuple(
            mode
            for mode in module.owner_modes(Q, speed)
            if mode[3] <= max_phase_size
        )

    ordered = tuple(sorted(subset, key=lambda speed: (len(bank(speed)), speed)))
    banks = tuple(bank(speed) for speed in ordered)
    max_sheet_sizes = tuple(max(len(mode[0]) for mode in modes) for modes in banks)

    def visit(index, assigned, covered):
        if len(covered) + sum(max_sheet_sizes[index:]) < Q:
            return None
        if index == len(ordered):
            if len(covered) != Q:
                return None
            cochain = module.potential_witness(Q, ordered, tuple(assigned))
            if cochain is None:
                return None
            return ordered, tuple(assigned), cochain

        speed = ordered[index]
        for mode in banks[index]:
            if any(
                not module.gap_values(
                    Q, ordered[prior], speed, assigned[prior], mode
                )
                for prior in range(index)
            ):
                continue
            result = visit(
                index + 1,
                assigned + [mode],
                covered | set(mode[0]),
            )
            if result is not None:
                return result
        return None

    return visit(0, [], set())


def sign_representative(residue):
    residue %= Q
    return min(residue, (-residue) % Q)


def main():
    for name, path, expected in PINNED:
        require(lf_hash(path) == expected, ("dependency changed", name, path))
    module = load_modes()

    source = Path(__file__)
    tree = ast.parse(source.read_text(encoding="utf-8"))
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)), "assert node")
    require(
        not any(
            isinstance(node, ast.Constant) and isinstance(node.value, float)
            for node in ast.walk(tree)
        ),
        "floating literal",
    )

    minimal_edges = []
    lower_cover_counts = []
    for size in range(1, SEARCH_MAX_RANK + 1):
        cover_count = 0
        for subset in combinations(VERTICES, size):
            if any(set(edge) <= set(subset) for edge in minimal_edges):
                cover_count += 1
                continue
            if module.has_cover(Q, subset):
                cover_count += 1
                minimal_edges.append(subset)
        lower_cover_counts.append(cover_count)
    minimal_edges = tuple(minimal_edges)
    require(tuple(lower_cover_counts[:5]) == (0, 0, 0, 0, 0), lower_cover_counts)
    require(len(minimal_edges) == 157, ("edge count", len(minimal_edges)))
    require(all(len(edge) == 6 for edge in minimal_edges), "edge rank")
    require(CANONICAL_EDGE == minimal_edges[0], ("first edge", minimal_edges[0]))

    witnesses = []
    domino_sufficient = []
    triphase_required = []
    for edge in minimal_edges:
        witness = module.mode_search(Q, edge)
        require(witness is not None, ("missing mode witness", edge))
        require(
            set().union(*(set(mode[0]) for mode in witness[1])) == set(range(Q)),
            ("mode cover", edge),
        )
        witnesses.append((edge, witness))
        if restricted_mode_search(module, edge, 2) is None:
            triphase_required.append(edge)
        else:
            domino_sufficient.append(edge)
    witnesses = tuple(witnesses)
    domino_sufficient = tuple(domino_sufficient)
    triphase_required = tuple(triphase_required)
    require(domino_sufficient == DOMINO_SUFFICIENT_EDGES, domino_sufficient)
    require(len(triphase_required) == 155, len(triphase_required))

    canonical_witness = dict(witnesses)[CANONICAL_EDGE]
    expected_order = (5, 3, 1, 2, 4, 7)
    expected_blocks = (
        (0, 3, 6, 9, 12),
        (0, 5, 10),
        (0, 1, 14),
        (0, 7, 8),
        (0, 4, 11),
        (0, 2, 13),
    )
    require(canonical_witness[0] == expected_order, canonical_witness[0])
    require(tuple(mode[0] for mode in canonical_witness[1]) == expected_blocks, canonical_witness[1])
    require(tuple(mode[1] for mode in canonical_witness[1]) == (0,) * 6, "centres")
    require(tuple(mode[2] for mode in canonical_witness[1]) == (15, 15, 1, 1, 1, 1), "widths")
    require(tuple(mode[3] for mode in canonical_witness[1]) == (1, 1, 3, 3, 3, 3), "phase sizes")
    require(all(gap[2] == 0 for gap in canonical_witness[2]), "zero cochain")

    zero_divisors = {residue for residue in range(Q) if gcd(residue, Q) != 1}
    kernel_union = set(expected_blocks[0]) | set(expected_blocks[1])
    require(kernel_union == zero_divisors, (kernel_union, zero_divisors))
    units = set(range(Q)) - zero_divisors
    unit_blocks = expected_blocks[2:]
    require(set().union(*(set(block) - {0} for block in unit_blocks)) == units, unit_blocks)
    require(
        tuple(
            set(block) - {0}
            for block in unit_blocks
        )
        == ({1, 14}, {7, 8}, {4, 11}, {2, 13}),
        unit_blocks,
    )
    unit_representatives = (1, 2, 4, 7)
    require(
        tuple(sign_representative(2 * value) for value in unit_representatives)
        == (2, 4, 7, 1),
        "unit sign quotient cycle",
    )

    for index, block in enumerate(unit_blocks):
        private = set(block) - {0}
        other_union = kernel_union | set().union(
            *(set(other) for other_index, other in enumerate(unit_blocks) if other_index != index)
        )
        require(private.isdisjoint(other_union), ("private unit sheets", block, other_union))
    require(restricted_mode_search(module, CANONICAL_EDGE, 2) is None, "canonical domino hostile")

    edge_sets = tuple(frozenset(edge) for edge in minimal_edges)
    independence_profile = tuple(
        sum(
            not any(edge <= set(subset) for edge in edge_sets)
            for subset in combinations(VERTICES, size)
        )
        for size in range(len(VERTICES) + 1)
    )
    gcd_profile = tuple(
        sorted(
            Counter(
                tuple(sorted((gcd(speed, Q) for speed in edge), reverse=True))
                for edge in minimal_edges
            ).items()
        )
    )
    require(gcd_profile == (((3, 1, 1, 1, 1, 1), 1), ((5, 3, 1, 1, 1, 1), 156)), gcd_profile)
    require((Q + 6) // 7 == 3 and (14 + 6) // 7 == 2, "first local trimode boundary")
    require(not any(speed % Q == 0 for speed in VERTICES), "outside literal core atlas")

    semantic = ExactDigest()
    semantic.add(("edges", minimal_edges))
    semantic.add(("witnesses", witnesses))
    semantic.add(("domino_sufficient", domino_sufficient))
    semantic.add(("triphase_required", triphase_required))
    semantic.add(("independence", independence_profile))
    semantic.add(("gcd_profile", gcd_profile))
    digest = semantic.hexdigest()
    require(digest == EXPECTED_SEMANTIC_DIGEST, ("semantic digest", digest))

    print("Q15 FIRST EFFECTIVE THREE-PHASE MODE EXACT PROBE")
    print(f"source_sha256_lf={lf_hash(source)}")
    print(f"dependency_sha256_lf={tuple((name, expected) for name, _, expected in PINNED)}")
    print("status=FINITE-EXACT q15 rank<=6 boundary model plus analytic private-sheet hostile;unnumbered_and_not_canon")
    print(f"vertices={VERTICES};cover_counts_rank1_to_6={tuple(lower_cover_counts)}")
    print(f"minimal_edges_through_rank6={len(minimal_edges)};rank_le6_profile=((6,157),);edge_sha256={sha256(repr(minimal_edges).encode('ascii')).hexdigest()}")
    print(f"gcd_profile={gcd_profile}")
    print(f"profile_generated_by_rank_le6_edges={independence_profile}")
    print(f"triphase_required={len(triphase_required)};domino_sufficient={domino_sufficient}")
    print(f"canonical_edge={CANONICAL_EDGE};ordered={canonical_witness[0]};blocks={expected_blocks}")
    print("canonical_mode_profile=(1,1,3,3,3,3);all_centres_h=0;all_15_cochain_gaps=0")
    print("crt_partition=ker(5)Uker(3)=zero_divisors;four_triphases_minus_zero=unit_sign_pairs")
    print("unit_sign_quotient=(Z/15Z)^x/{+-1}=(1,2,4,7);times2_cycle=(1,2,4,7,1)")
    print("hostile=each_coprime_triphase_has_two_private_unit_sheets;every_domino_subblock_loses_one;no_size<=2_mode_cover")
    print("boundary=q<=14_has_only_singleton_or_domino_witness_modes;q15_first_local_and_rank6_effective_triphase")
    print("scope=q15_rank_le6_boundary;q15_has_no_rank<=5_cover_and_no_q_divisible_literal_speed;higher_rank_physical_clutter_not_enumerated;not_a_T3387_body_row;no_LRC14_consequence")
    print(f"semantic_sha256={digest}")
    print("verdict=PASS")


if __name__ == "__main__":
    main()
