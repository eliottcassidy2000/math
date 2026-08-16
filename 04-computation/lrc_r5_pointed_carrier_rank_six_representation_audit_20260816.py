#!/usr/bin/env python3
"""Identify the common rank-six endpoint quotient with the pointed-six carrier.

This derived audit replays the pinned current-branch x root-difference joint
construction and then performs a basis-free representation comparison not used
by the parent probe.  It checks that difference marginalization is injective on
the six pointed response rows, that their image is the previously pinned common
six-dimensional endpoint space, and that pointed-path reversal induces the same
abstract involution as current branch/path reversal.

The joint integrand is imported from its candidate implementation, so this is
not a clean-room audit of the four-way physical construction.  It is an exact
audit of the stated representation-theoretic consequence.
"""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from hashlib import sha256
import importlib
import json
from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

JOINT = importlib.import_module(
    "lrc_r5_owner_node_inverse_branch_root_difference_four_way_probe_20260816"
)

P = 13
V = 4
PRIME = JOINT.PRIME
EXPECTED_COMMON_W = "6e9083f15408f6d2d85fb3f2747ba0bd1f987e83ce4b836cb7298aaccc84e0c4"
EXPECTED_POINTED_REVERSAL = (5, 4, 3, 2, 1, 0)
EXPECTED_ROWSPACE_DIGESTS = (
    "6014f99bee469da11163238341e815f6901d1460ec2f3d1ac62515b65470fda5",
    EXPECTED_COMMON_W,
    EXPECTED_COMMON_W,
    "67850c2586ac9a4f242ada7173eba2fe1b3c5e6d1ae20a5e2b6d7ad0dc7d3f88",
    "67850c2586ac9a4f242ada7173eba2fe1b3c5e6d1ae20a5e2b6d7ad0dc7d3f88",
)
EXPECTED_SEMANTIC_SHA256 = (
    "5e491136809fc164bbbfc7aeb9c272b7aff05020992985a39f291bf31297903e"
)


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def digest_json(value: object) -> str:
    return sha256(
        json.dumps(value, separators=(",", ":"), sort_keys=True).encode("ascii")
    ).hexdigest()


def rank(rows) -> int:
    return JOINT.C.rank_mod(tuple(tuple(value % PRIME for value in row) for row in rows))


def canonical_row_basis(rows):
    matrix = [list(value % PRIME for value in row) for row in rows]
    if not matrix:
        return ()
    columns = len(matrix[0])
    require(all(len(row) == columns for row in matrix), "ragged matrix")
    pivot_row = 0
    for column in range(columns):
        pivot = next(
            (row for row in range(pivot_row, len(matrix)) if matrix[row][column]),
            None,
        )
        if pivot is None:
            continue
        matrix[pivot_row], matrix[pivot] = matrix[pivot], matrix[pivot_row]
        inverse = pow(matrix[pivot_row][column], -1, PRIME)
        matrix[pivot_row] = [value * inverse % PRIME for value in matrix[pivot_row]]
        for row in range(len(matrix)):
            if row == pivot_row or not matrix[row][column]:
                continue
            factor = matrix[row][column]
            matrix[row] = [
                (left - factor * right) % PRIME
                for left, right in zip(matrix[row], matrix[pivot_row])
            ]
        pivot_row += 1
        if pivot_row == len(matrix):
            break
    return tuple(tuple(row) for row in matrix[:pivot_row])


def rowspace_digest(rows) -> str:
    return digest_json(canonical_row_basis(rows))


def flatten_bank(chunks, position: int):
    return tuple(
        row
        for chunk in chunks
        for beta_rows in chunk[position]
        for row in beta_rows
    )


def owner_rows(owner_tensor):
    return tuple(
        owner_tensor[state][branch]
        for state in range(V)
        for branch in range(P)
    )


def pointed_expanded_rows(pointed_tensor):
    return tuple(
        tuple(
            pointed_tensor[point][difference][relation]
            for difference in range(P)
            for relation in range(P)
        )
        for point in range(len(JOINT.POINTED_STATES))
    )


def pointed_relation_rows(pointed_tensor):
    return tuple(
        tuple(
            sum(pointed_tensor[point][difference][relation] for difference in range(P))
            % PRIME
            for relation in range(P)
        )
        for point in range(len(JOINT.POINTED_STATES))
    )


def pointed_reversal():
    lookup = {pair: index for index, pair in enumerate(JOINT.POINTED_STATES)}
    return tuple(
        lookup[(state ^ 2, P - 1 - root)]
        for state, root in JOINT.POINTED_STATES
    )


def owner_reversal(rows):
    return tuple(
        rows[(state ^ 2) * P + (P - 1 - branch)]
        for state in range(V)
        for branch in range(P)
    )


def graph(rows, reflected):
    return tuple(tuple(left) + tuple(right) for left, right in zip(rows, reflected))


def eigenspace_ranks(rows, reflected):
    plus = tuple(
        tuple((left + right) % PRIME for left, right in zip(row, image))
        for row, image in zip(rows, reflected)
    )
    minus = tuple(
        tuple((left - right) % PRIME for left, right in zip(row, image))
        for row, image in zip(rows, reflected)
    )
    return rank(plus), rank(minus)


def main() -> None:
    with ProcessPoolExecutor(max_workers=4) as pool:
        chunks = tuple(pool.map(JOINT.worker, range(P)))
    require(tuple(chunk[0] for chunk in chunks) == tuple(range(P)), "worker order")

    weighted_gamma = flatten_bank(chunks, 1)
    pointed_gamma = flatten_bank(chunks, 3)
    require(digest_json(pointed_gamma) == JOINT.POINTED_GAMMA_SHA256,
            "pointed gamma drift")

    zeta = JOINT.C.context()["zeta"]
    actual = JOINT.inverse_tensor(weighted_gamma, zeta)
    pointed = JOINT.inverse_pointed_tensor(pointed_gamma, zeta)
    require(digest_json(pointed) == JOINT.POINTED_TENSOR_SHA256,
            "pointed tensor drift")
    require(JOINT.pointed_carrier_record(actual, pointed)
            == ((1, 1, 1, 1), (2, 2, 2, 2),
                (1, 1, 1, 1), (2, 2, 2, 2)),
            "statewise pointed factorization")

    owner_tensor = JOINT.tensor_difference_marginal(actual)
    require(digest_json(owner_tensor) == JOINT.R.EXPECTED_DIGESTS[1],
            "owner marginal drift")
    owner = owner_rows(owner_tensor)
    expanded = pointed_expanded_rows(pointed)
    pointed_relation = pointed_relation_rows(pointed)

    relation_record = (
        rank(expanded),
        rank(pointed_relation),
        rank(owner),
        rank(tuple(pointed_relation) + tuple(owner)),
    )
    require(relation_record == (6, 6, 6, 6),
            ("pointed marginal identification", relation_record))

    reversal = pointed_reversal()
    require(reversal == EXPECTED_POINTED_REVERSAL, ("pointed reversal", reversal))
    pointed_reflected = tuple(pointed_relation[index] for index in reversal)
    owner_reflected = owner_reversal(owner)
    pointed_graph = graph(pointed_relation, pointed_reflected)
    owner_graph = graph(owner, owner_reflected)
    representation_record = (
        rank(pointed_graph),
        rank(owner_graph),
        rank(tuple(pointed_graph) + tuple(owner_graph)),
        eigenspace_ranks(pointed_relation, pointed_reflected),
        eigenspace_ranks(owner, owner_reflected),
    )
    require(representation_record == (6, 6, 6, (3, 3), (3, 3)),
            ("pointed representation", representation_record))

    rowspace_digests = (
        rowspace_digest(expanded),
        rowspace_digest(pointed_relation),
        rowspace_digest(owner),
        rowspace_digest(pointed_graph),
        rowspace_digest(owner_graph),
    )
    require(rowspace_digests[1] == rowspace_digests[2] == EXPECTED_COMMON_W,
            ("common endpoint quotient", rowspace_digests))
    require(rowspace_digests[3] == rowspace_digests[4],
            ("common involution graph", rowspace_digests))
    require(rowspace_digests == EXPECTED_ROWSPACE_DIGESTS,
            ("row-space semantic drift", rowspace_digests))

    record = (
        JOINT.POINTED_STATES,
        reversal,
        relation_record,
        representation_record,
        rowspace_digests,
        digest_json(actual),
        digest_json(pointed),
    )
    semantic = digest_json(record)
    require(semantic == EXPECTED_SEMANTIC_SHA256,
            ("semantic drift", semantic, EXPECTED_SEMANTIC_SHA256))

    print("== r=5 pointed-carrier representation audit ==")
    print(f"pointed_states={JOINT.POINTED_STATES}")
    print(f"pointed_path_reversal={reversal}")
    print("relation_record_order=(expanded_pointed_rank,marginal_pointed_rank,owner_rank,union_rank)")
    print(f"difference_marginal_identification={relation_record}")
    print("representation_order=(point_graph_rank,owner_graph_rank,joint_graph_rank,point_+-,owner_+-)")
    print(f"pointed_permutation_module={representation_record}")
    print(f"rowspace_digests_(expanded_pointed,marginal_pointed,owner,point_graph,owner_graph)={rowspace_digests}")
    print(f"tensor_digests_(joint,pointed)=({digest_json(actual)},{digest_json(pointed)})")
    print(f"semantic_sha256={semantic}")
    print("verdict=common rank-six endpoint quotient is the injective difference marginal of the six-pointed carrier")
    print("typing=pointed cut-state incidences with path reversal;not tournament vertices or radix branches")
    print("scope=derived exact representation audit of candidate joint integrand;physical four-way independent audit pending")
    print("no_exact_address_chronology_current_row_exclusion_or_LRC14")


if __name__ == "__main__":
    main()
