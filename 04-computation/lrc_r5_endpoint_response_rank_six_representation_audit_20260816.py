#!/usr/bin/env python3
"""Exact representation audit for the two rank-six branch response tables.

This probe compares two already typed refinements of the same owner-visible
Boolean-square endpoint table:

* ``b_source=floor(n/13)``, the high digit of the left source-time ``P_169``
  inverse branch; and
* ``r_owner=a mod 13``, the last digit of the current-leg ``P_(13^5)`` sheet.

The source table is rebuilt through the candidate implementation whose exact
table digest was independently reproduced by the clean-room audit.  The
current-leg table is rebuilt through its disjoint radix-fold implementation.
Both tables are then treated only as row spaces in the common, typed endpoint
relation carrier ``Fun(F_13,F_p)``.  The decisive hostile is the rank of their
stack: equality at rank six would exhibit one common response sector, while a
larger rank falsifies that identification.

No joint ancestry table, temporal intertwiner, exact THM-2334 address, physical
current, row exclusion, or LRC(14) conclusion is constructed here.
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

SOURCE = importlib.import_module(
    "lrc_r5_ufull_owner_node_boolean_square_branch_sheet_sidecar_probe_20260816"
)
OWNER = importlib.import_module(
    "lrc_r5_ufull_owner_node_boolean_square_inverse_owner_branch_probe_20260816"
)

P = 13
V = 4
PRIME = OWNER.PRIME

EXPECTED_RAW_PAIR = (6, 6, 4, 6, 6, 2)
EXPECTED_PURE_PAIR = (6, 6, 0, 6, 6, 6)
EXPECTED_WITHIN = (7, 7)
EXPECTED_DESCENT = (
    (6, 3, 3), (6, 3, 3), 6,
    (6, 3, 3), (6, 3, 3), 6,
)
EXPECTED_CENTERING = (
    (6, 6, 7, 7, 7),
    (6, 6, 7, 7, 7),
)
EXPECTED_DIGESTS = (
    "504e91431cf4d55f6e3cb5aa7c6e6f2c34db538cb692c67857a5c7627b252261",
    "2a195fac7fbd60a4bbd2597bf34baf87302ead16c7c0e8c67c0667b0d320dfba",
    "5f4b9609faaa5f148d112a7cde5cfba0ab2c1385b4c53ea9c4bcfc6e93d106fc",
    "2a289e3f4de379985083c48841557da6f6777d21b6f6e7354b83b1b45ecb8e29",
    "fc362ff7a4fe98ba4226fa1624235e134b43d88a76af8c1cfa09948094159c38",
)
EXPECTED_CANONICAL_SPACES = (
    "6e9083f15408f6d2d85fb3f2747ba0bd1f987e83ce4b836cb7298aaccc84e0c4",
    "6e9083f15408f6d2d85fb3f2747ba0bd1f987e83ce4b836cb7298aaccc84e0c4",
    "ae86ad2a3fd03bea95c823b2454b78f244581aa048cb2da63a03a75f484cc596",
    "ae86ad2a3fd03bea95c823b2454b78f244581aa048cb2da63a03a75f484cc596",
    "1b4fef00a23dcb79dc52ace662bae2f858ce3da27b6ef19b902ae40f5a79e755",
    "67850c2586ac9a4f242ada7173eba2fe1b3c5e6d1ae20a5e2b6d7ad0dc7d3f88",
    "67850c2586ac9a4f242ada7173eba2fe1b3c5e6d1ae20a5e2b6d7ad0dc7d3f88",
    "579a14c27833b79eb5cdc045c00bd93d05e9b7bee191ced128436e15406800e8",
)
EXPECTED_SEMANTIC_SHA256 = (
    "7baeb128a4c4d5998342611fdcf821d002ffb4622692dd72c03e6f11c8d9825a"
)


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def digest_json(value: object) -> str:
    return sha256(
        json.dumps(value, separators=(",", ":"), sort_keys=True).encode("ascii")
    ).hexdigest()


def rank(rows) -> int:
    return OWNER.C.rank_mod(tuple(tuple(value % PRIME for value in row) for row in rows))


def canonical_row_basis(rows):
    """Return the unique reduced-row-echelon basis over the audit field."""
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


def build_source_table():
    """Rebuild the source-time high-digit table with its pinned exact digest."""
    with ProcessPoolExecutor(max_workers=4) as pool:
        chunks = tuple(pool.map(SOURCE.worker, range(P)))
    require(tuple(chunk[0] for chunk in chunks) == tuple(range(P)), "source worker order")
    gamma = tuple(row for chunk in chunks for beta_rows in chunk[1] for row in beta_rows)
    table = SOURCE.inverse_table(gamma, SOURCE.B.context()["zeta"])
    require(digest_json(table) == SOURCE.EXPECTED_DIGESTS[2], "source table drift")
    return table


def build_owner_tensor():
    """Rebuild the current-leg last-digit tensor with its pinned exact digest."""
    with ProcessPoolExecutor(max_workers=4) as pool:
        chunks = tuple(pool.map(OWNER.worker, range(P)))
    require(tuple(chunk[0] for chunk in chunks) == tuple(range(P)), "owner worker order")
    gamma = tuple(
        row for chunk in chunks for beta_rows in chunk[1] for row in beta_rows
    )
    tensor = OWNER.inverse_tensor(gamma, OWNER.C.context()["zeta"])
    require(digest_json(tensor) == OWNER.EXPECTED_DIGESTS[1], "owner tensor drift")
    return tensor


def owner_rows(tensor):
    return tuple(
        tensor[state][branch]
        for state in range(V)
        for branch in range(P)
    )


def source_parent(table):
    return tuple(
        tuple(
            sum(table[branch * V + state][relation] for branch in range(P)) % PRIME
            for relation in range(P)
        )
        for state in range(V)
    )


def owner_parent(tensor):
    return tuple(
        tuple(
            sum(tensor[state][branch][relation] for branch in range(P)) % PRIME
            for relation in range(P)
        )
        for state in range(V)
    )


def permute_relation(rows, multiplier: int, shift: int):
    """Pull back rows along t |-> multiplier*t+shift on F_13."""
    return tuple(
        tuple(row[(multiplier * relation + shift) % P] for relation in range(P))
        for row in rows
    )


def affine_stabilizer(rows):
    base_rank = rank(rows)
    return tuple(
        (multiplier, shift)
        for multiplier in range(1, P)
        for shift in range(P)
        if rank(tuple(rows) + permute_relation(rows, multiplier, shift)) == base_rank
    )


def reflection_splits(rows):
    base_rank = rank(rows)
    answer = []
    for shift in range(P):
        reflected = permute_relation(rows, -1, shift)
        if rank(tuple(rows) + reflected) != base_rank:
            continue
        plus = tuple(
            tuple((left + right) % PRIME for left, right in zip(row, image))
            for row, image in zip(rows, reflected)
        )
        minus = tuple(
            tuple((left - right) % PRIME for left, right in zip(row, image))
            for row, image in zip(rows, reflected)
        )
        answer.append((shift, rank(plus), rank(minus)))
    return tuple(answer)


def simultaneous_reflections(tensor, source_order: bool):
    """Test the typed branch/path reversal against every affine target map.

    The path involution is ``state -> state XOR 2``.  Both radix labels use
    ``digit -> 12-digit``.  This is deliberately not a tournament action.
    """
    def value(state: int, branch: int, relation: int) -> int:
        if source_order:
            return tensor[branch * V + state][relation]
        return tensor[state][branch][relation]

    solutions = []
    for multiplier in range(1, P):
        for shift in range(P):
            for sign in (1, -1):
                if all(
                    value(state, branch, relation)
                    == sign
                    * value(
                        state ^ 2,
                        P - 1 - branch,
                        (multiplier * relation + shift) % P,
                    )
                    % PRIME
                    for state in range(V)
                    for branch in range(P)
                    for relation in range(P)
                ):
                    solutions.append((multiplier, shift, sign))
    return tuple(solutions)


def domain_reflection(rows, source_order: bool):
    """Apply the typed state/digit reversal to the domain basis labels."""
    if source_order:
        return tuple(
            rows[(P - 1 - branch) * V + (state ^ 2)]
            for branch in range(P)
            for state in range(V)
        )
    return tuple(
        rows[(state ^ 2) * P + (P - 1 - branch)]
        for state in range(V)
        for branch in range(P)
    )


def involution_descent(rows, source_order: bool):
    """Test whether domain reversal descends to the six-dimensional image.

    For the response map ``phi(e_i)=row_i``, descent is equivalent to
    ``ker(phi)`` being reversal-stable.  In matrix language the columns of the
    reflected response matrix must lie in the original column space, tested by
    one exact horizontal block rank.
    """
    reflected = domain_reflection(rows, source_order)
    horizontal = tuple(tuple(left) + tuple(right) for left, right in zip(rows, reflected))
    plus = tuple(
        tuple((left + right) % PRIME for left, right in zip(row, image))
        for row, image in zip(rows, reflected)
    )
    minus = tuple(
        tuple((left - right) % PRIME for left, right in zip(row, image))
        for row, image in zip(rows, reflected)
    )
    return rank(horizontal), rank(plus), rank(minus), horizontal


def parent_involution_descent(parent):
    reflected = tuple(parent[state ^ 2] for state in range(V))
    horizontal = tuple(tuple(left) + tuple(right) for left, right in zip(parent, reflected))
    plus = tuple(
        tuple((left + right) % PRIME for left, right in zip(row, image))
        for row, image in zip(parent, reflected)
    )
    minus = tuple(
        tuple((left - right) % PRIME for left, right in zip(row, image))
        for row, image in zip(parent, reflected)
    )
    return rank(horizontal), rank(plus), rank(minus), horizontal


def centre_relation(rows):
    inverse = pow(P, -1, PRIME)
    return tuple(
        tuple((value - sum(row) * inverse) % PRIME for value in row)
        for row in rows
    )


def centering_record(raw_rows, pure_rows):
    centred = centre_relation(raw_rows)
    ones = ((1,) * P,)
    return (
        rank(centred),
        rank(tuple(centred) + tuple(pure_rows)),
        rank(tuple(raw_rows) + tuple(centred)),
        rank(tuple(raw_rows) + ones),
        rank(tuple(raw_rows) + tuple(pure_rows)),
    )


def interaction_source(table):
    return SOURCE.flatten_tensor(SOURCE.three_way_interaction(table))


def interaction_owner(tensor):
    return owner_rows(OWNER.three_way_interaction(tensor))


def pair_record(left, right, parent):
    left_rank = rank(left)
    right_rank = rank(right)
    parent_rank = rank(parent)
    stack_rank = rank(tuple(left) + tuple(right))
    intersection = left_rank + right_rank - stack_rank
    quotient_union = stack_rank - parent_rank
    return (
        left_rank,
        right_rank,
        parent_rank,
        stack_rank,
        intersection,
        quotient_union,
    )


def main() -> None:
    require(PRIME == SOURCE.B.PRIME, ("field mismatch", PRIME, SOURCE.B.PRIME))
    source = build_source_table()
    owner_tensor = build_owner_tensor()
    owner = owner_rows(owner_tensor)

    source_base = source_parent(source)
    owner_base = owner_parent(owner_tensor)
    require(source_base == owner_base, "typed Boolean-square parents differ")
    require(rank(source_base) == 4, "Boolean-square parent rank")

    source_three = interaction_source(source)
    owner_three = interaction_owner(owner_tensor)

    raw = pair_record(source, owner, source_base)
    pure = pair_record(source_three, owner_three, ())
    within = (
        rank(tuple(source) + tuple(source_three)),
        rank(tuple(owner) + tuple(owner_three)),
    )

    source_stabilizer = affine_stabilizer(source)
    owner_stabilizer = affine_stabilizer(owner)
    common_stabilizer = tuple(sorted(set(source_stabilizer) & set(owner_stabilizer)))
    source_reflections = reflection_splits(source)
    owner_reflections = reflection_splits(owner)
    source_three_reflections = reflection_splits(source_three)
    owner_three_reflections = reflection_splits(owner_three)

    source_simultaneous = simultaneous_reflections(source, source_order=True)
    owner_simultaneous = simultaneous_reflections(owner_tensor, source_order=False)

    source_descent = involution_descent(source, source_order=True)
    owner_descent = involution_descent(owner, source_order=False)
    source_three_descent = involution_descent(source_three, source_order=True)
    owner_three_descent = involution_descent(owner_three, source_order=False)
    raw_graph_stack_rank = rank(
        tuple(source_descent[3]) + tuple(owner_descent[3])
    )
    pure_graph_stack_rank = rank(
        tuple(source_three_descent[3]) + tuple(owner_three_descent[3])
    )
    descent_record = (
        source_descent[:3],
        owner_descent[:3],
        raw_graph_stack_rank,
        source_three_descent[:3],
        owner_three_descent[:3],
        pure_graph_stack_rank,
    )
    parent_descent_full = parent_involution_descent(source_base)
    parent_graph_in_raw_rank = rank(
        tuple(source_descent[3]) + tuple(parent_descent_full[3])
    )
    parent_descent = (
        parent_descent_full[0],
        parent_descent_full[1],
        parent_descent_full[2],
        parent_graph_in_raw_rank,
        source_descent[1] - parent_descent_full[1],
        source_descent[2] - parent_descent_full[2],
    )
    centering = (
        centering_record(source, source_three),
        centering_record(owner, owner_three),
    )

    digests = (
        digest_json(source),
        digest_json(owner_tensor),
        digest_json(source_base),
        digest_json(source_three),
        digest_json(owner_three),
    )
    canonical_spaces = (
        rowspace_digest(source),
        rowspace_digest(owner),
        rowspace_digest(source_three),
        rowspace_digest(owner_three),
        rowspace_digest(source_base),
        rowspace_digest(source_descent[3]),
        rowspace_digest(owner_descent[3]),
        rowspace_digest(parent_descent_full[3]),
    )
    record = (
        PRIME,
        raw,
        pure,
        within,
        source_stabilizer,
        owner_stabilizer,
        common_stabilizer,
        source_reflections,
        owner_reflections,
        source_three_reflections,
        owner_three_reflections,
        source_simultaneous,
        owner_simultaneous,
        descent_record,
        parent_descent,
        centering,
        digests,
        canonical_spaces,
    )
    semantic = digest_json(record)

    require(raw == EXPECTED_RAW_PAIR, ("raw pair", raw))
    require(pure == EXPECTED_PURE_PAIR, ("pure pair", pure))
    require(within == EXPECTED_WITHIN, ("within", within))
    require(source_stabilizer == owner_stabilizer == ((1, 0),),
            ("affine stabilizers", source_stabilizer, owner_stabilizer))
    require(not source_reflections and not owner_reflections,
            ("raw affine reflections", source_reflections, owner_reflections))
    require(not source_three_reflections and not owner_three_reflections,
            ("pure affine reflections", source_three_reflections, owner_three_reflections))
    require(not source_simultaneous and not owner_simultaneous,
            ("entrywise affine reflection", source_simultaneous, owner_simultaneous))
    require(descent_record == EXPECTED_DESCENT, ("abstract descent", descent_record))
    require(parent_descent == (4, 2, 2, 6, 1, 1),
            ("parent subrepresentation", parent_descent))
    require(centering == EXPECTED_CENTERING, ("centering", centering))
    require(canonical_spaces[0] == canonical_spaces[1], "raw row-space hash mismatch")
    require(canonical_spaces[2] == canonical_spaces[3], "pure row-space hash mismatch")
    require(canonical_spaces[5] == canonical_spaces[6], "involution graph hash mismatch")
    require(digests == EXPECTED_DIGESTS, ("table digests", digests))
    require(canonical_spaces == EXPECTED_CANONICAL_SPACES,
            ("canonical row spaces", canonical_spaces))
    require(semantic == EXPECTED_SEMANTIC_SHA256,
            ("semantic drift", semantic, EXPECTED_SEMANTIC_SHA256))

    print("== r=5 shared rank-six endpoint-response representation audit ==")
    print("coordinates=(source b_source,state,relation);(owner r_owner,state,relation)")
    print(f"field=(prime={PRIME},relation_carrier=Fun(F13,Fp),dimension=13)")
    print("typed_parent=(same Boolean-square table): PASS")
    print("pair_record_order=(source_rank,owner_rank,parent_rank,stack_rank,intersection_dim,union_mod_parent_dim)")
    print(f"raw_pair={raw}")
    print("pure_pair_order=(source_rank,owner_rank,empty_parent_rank,stack_rank,intersection_dim,stack_rank)")
    print(f"pure_three_way_pair={pure}")
    print(f"raw_with_own_pure_stack_ranks_(source,owner)={within}")
    print(f"affine_stabilizers_(source,owner,common)=({source_stabilizer},{owner_stabilizer},{common_stabilizer})")
    print(f"raw_reflection_splits_(source,owner)=({source_reflections},{owner_reflections})")
    print(f"pure_reflection_splits_(source,owner)=({source_three_reflections},{owner_three_reflections})")
    print(f"typed_simultaneous_path_branch_target_affine_(source,owner)=({source_simultaneous},{owner_simultaneous})")
    print("descent_order=(source(horizontal,+,-),owner(horizontal,+,-),joint_graph_rank;then_pure)")
    print(f"abstract_involution_descent={descent_record}")
    print("parent_descent_order=(horizontal,+,-,graph_in_raw,quotient_+,quotient_-)")
    print(f"Boolean_parent_subrepresentation={parent_descent}")
    print("centering_order=(centred_rank,centred_plus_pure_rank,raw_plus_centred_rank,raw_plus_constants_rank,raw_plus_pure_rank)")
    print(f"augmentation_projection_(source,owner)={centering}")
    print(f"digests_(source,owner,parent,source_pure,owner_pure)={digests}")
    print(f"canonical_rowspace_digests_(source,owner,source_pure,owner_pure,parent,source_graph,owner_graph,parent_graph)={canonical_spaces}")
    print(f"semantic_sha256={semantic}")
    print("typing=path involution on states;radix-tree reversal on digits;endpoint affine action on relation")
    print("not_typing=cut arcs are marginalized;tournament relation absent;no joint ancestry/current map")
    print("scope=finite exact target-row-space comparison only;no physical current,row exclusion,or LRC14")


if __name__ == "__main__":
    main()
