#!/usr/bin/env python3
"""Expose the pointed-six carrier as a bidirected P4 arc cocycle.

This is a representation sidecar to the finite two-current-digit response
tensor.  It rebuilds that tensor from its pre-integration workers, extracts
the 169 diagonal address-conditioned carrier maps, and asks exactly what
structure survives after the six pointed lines are retyped as directed arcs.

The four owner states are the endpoint-touching intervals

    0:{0}, 1:{0,6}, 3:{6,12}, 2:{12}.

Every retained root belongs to exactly two adjacent states.  Sending a
pointed incidence (state,root) to the arc from that state to the other state
containing the root identifies the six carrier lines with the six arcs of

    0 <-> 1 <-> 3 <-> 2.

The alternating state/root incidence graph is P7.  Its H1 is zero.  Adding
the missing Boolean-square closure edge 2--0 produces C4 and a one-dimensional
H1; the script records this as a typing gate, not as a physical current.

For each pointed arc e, the diagonal scalars k_e(r0,r1) form a 13 by 13
row-sum-one matrix over the exact split field.  These are static address
weights, not probabilities or a chronological Markov chain.  Arc reversal,
chamber reflection, support subsets, and harmonic-sum loss are tested
explicitly.
"""

from __future__ import annotations

from collections import defaultdict, deque
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
COMPUTATION = ROOT / "04-computation"
sys.path.insert(0, str(COMPUTATION))

import lrc_r5_two_current_digits_pointed_root_difference_carrier_transition_probe_20260816 as M


PARENT_PATH = (
    COMPUTATION
    / "lrc_r5_two_current_digits_pointed_root_difference_carrier_transition_probe_20260816.py"
)
PARENT_SHA256 = "9d1671e0f823fdbaa9ab79915ba05dbb4dda4c6eabb97fe4484baf4c2e3205f2"
EXPECTED_PROFILE_SHA256 = (
    "d1c7e561538ac6abb7c631c642d547aea7ca8dcfa33296e849a22a5061d2f595"
)
EXPECTED_SEMANTIC_SHA256 = (
    "b2ba313f88fbab0d36e95a63ade832492743ed6007fa0480434083d7dab0ecd3"
)


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


def lf_sha256(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def digest_json(value: object) -> str:
    return sha256(
        json.dumps(value, separators=(",", ":"), sort_keys=True).encode("ascii")
    ).hexdigest()


def histogram(values):
    result = defaultdict(int)
    for value in values:
        result[value] += 1
    return tuple(sorted(result.items()))


def build_actual_tensor():
    """Rebuild the parent tensor without importing any post-integration table."""

    _profiles, boundaries, cells, source_record = M.T2.two_digit_context()
    source_digest = M.digest_json((source_record, boundaries))
    require(source_digest == M.T2.EXPECTED_SOURCE_SHA256,
            ("source digest", source_digest, M.T2.EXPECTED_SOURCE_SHA256))

    tau_actual = [[[[[0 for _s in range(M.P)] for _r1 in range(M.P)]
                    for _r0 in range(M.P)] for _point in M.POINTS]
                  for _tau in range(M.P)]
    alpha_merkle = []
    with ProcessPoolExecutor(max_workers=4) as pool:
        for expected_alpha, chunk in enumerate(pool.map(M.worker, range(M.P))):
            require(chunk[0] == expected_alpha,
                    ("worker order", expected_alpha, chunk[0]))
            M.add_core(tau_actual, chunk[1])
            alpha_merkle.append(chunk[5])

    require(M.digest_json(tuple(alpha_merkle))
            == M.EXPECTED_JOINT_ALPHA_MERKLE_SHA256,
            "alpha Merkle drift")
    actual = M.inverse_core(tau_actual, M.C.context()["zeta"])
    tensor_digest = M.digest_json(actual)
    require(tensor_digest == M.EXPECTED_TENSOR_DIGESTS[0][0],
            ("tensor drift", tensor_digest, M.EXPECTED_TENSOR_DIGESTS[0][0]))
    return cells, actual, source_digest, tensor_digest


def extract_kernels(tensor):
    """Return k[point][r0][r1] from the unique diagonal 6-line maps."""

    parent = M.tensor_one_digit_pointed(tensor)
    profile = []
    for r0 in range(M.P):
        parent_rows = tuple(
            tuple(parent[point][r0][difference][relation]
                  for difference in range(M.P) for relation in range(M.P))
            for point in range(len(M.POINTS))
        )
        require(M.rank_rows(parent_rows) == 6, ("parent rank", r0))
        r0_profile = []
        for r1 in range(M.P):
            child_rows = tuple(
                tuple(tensor[point][r0][r1][difference][relation]
                      for difference in range(M.P) for relation in range(M.P))
                for point in range(len(M.POINTS))
            )
            matrix, witness = M.coordinate_matrix(parent_rows, child_rows)
            require(witness is None, ("coordinate failure", r0, r1, witness))
            require(all(matrix[i][j] == 0
                        for i in range(6) for j in range(6) if i != j),
                    ("nondiagonal map", r0, r1))
            r0_profile.append(tuple(matrix[point][point] for point in range(6)))
        profile.append(tuple(r0_profile))

    profile = tuple(profile)
    profile_digest = M.digest_json(profile)
    require(profile_digest == EXPECTED_PROFILE_SHA256,
            ("diagonal profile drift", profile_digest, EXPECTED_PROFILE_SHA256))
    kernels = tuple(
        tuple(tuple(profile[r0][r1][point] for r1 in range(M.P))
              for r0 in range(M.P))
        for point in range(6)
    )
    return profile, kernels, profile_digest


def transpose(matrix):
    return tuple(tuple(matrix[row][column] for row in range(len(matrix)))
                 for column in range(len(matrix[0])))


def subtract_identity(matrix):
    return tuple(
        tuple((value - int(row == column)) % M.PRIME
              for column, value in enumerate(line))
        for row, line in enumerate(matrix)
    )


def scc_count(matrix):
    """Count SCCs in the nonzero directed support graph."""

    n = len(matrix)
    graph = tuple(tuple(j for j, value in enumerate(row) if value)
                  for row in matrix)
    reverse = [[] for _ in range(n)]
    for i, row in enumerate(graph):
        for j in row:
            reverse[j].append(i)

    seen = set()
    order = []

    def dfs(start):
        stack = [(start, 0)]
        seen.add(start)
        while stack:
            vertex, index = stack[-1]
            if index < len(graph[vertex]):
                nxt = graph[vertex][index]
                stack[-1] = (vertex, index + 1)
                if nxt not in seen:
                    seen.add(nxt)
                    stack.append((nxt, 0))
            else:
                order.append(vertex)
                stack.pop()

    for vertex in range(n):
        if vertex not in seen:
            dfs(vertex)

    seen.clear()
    count = 0
    for start in reversed(order):
        if start in seen:
            continue
        count += 1
        queue = deque((start,))
        seen.add(start)
        while queue:
            vertex = queue.popleft()
            for nxt in reverse[vertex]:
                if nxt not in seen:
                    seen.add(nxt)
                    queue.append(nxt)
    return count


def kernel_record(kernels):
    records = []
    for point, matrix in enumerate(kernels):
        row_sums = tuple(sum(row) % M.PRIME for row in matrix)
        column_sums = tuple(sum(matrix[row][column] for row in range(M.P)) % M.PRIME
                            for column in range(M.P))
        require(row_sums == (1,) * M.P, ("row partition", point, row_sums))
        rank = M.rank_rows(matrix)
        right_one = M.P - M.rank_rows(subtract_identity(matrix))
        left_one = M.P - M.rank_rows(subtract_identity(transpose(matrix)))
        records.append((
            point,
            M.POINTS[point],
            rank,
            right_one,
            left_one,
            scc_count(matrix),
            histogram(sum(value != 0 for value in row) for row in matrix),
            len(set(matrix)),
            sum(value == 1 for value in column_sums),
            digest_json(matrix),
        ))
    return tuple(records)


def chain_record():
    """Verify P7/P4/C4 ranks and the explicit closure-cycle generator."""

    # State order is the Gray path 0--1--3--2.  Root vertices are inserted
    # between their two incident states, giving P7.
    state_order = (0, 1, 3, 2)
    subsets = ((0,), (0, 6), (6, 12), (12,))
    arcs = ((0, 1, 0), (1, 0, 0), (1, 3, 6),
            (3, 1, 6), (3, 2, 12), (2, 3, 12))
    require(tuple((arc[0], arc[2]) for arc in arcs) == M.POINTS,
            ("point/arc bijection", arcs, M.POINTS))

    p7_edges = ((0, 1), (2, 1), (2, 3), (4, 3), (4, 5), (6, 5))
    boundary_p7 = tuple(
        tuple((1 if vertex == head else -1 if vertex == tail else 0) % M.PRIME
              for tail, head in p7_edges)
        for vertex in range(7)
    )
    p7_rank = M.rank_rows(boundary_p7)

    p4_edges = ((0, 1), (1, 2), (2, 3))
    c4_edges = p4_edges + ((3, 0),)

    def boundary(vertex_count, edges):
        return tuple(
            tuple((1 if vertex == head else -1 if vertex == tail else 0) % M.PRIME
                  for tail, head in edges)
            for vertex in range(vertex_count)
        )

    p4_rank = M.rank_rows(boundary(4, p4_edges))
    c4_matrix = boundary(4, c4_edges)
    c4_rank = M.rank_rows(c4_matrix)
    cycle = (1, 1, 1, 1)
    require(all(sum(c4_matrix[row][column] * cycle[column]
                    for column in range(4)) % M.PRIME == 0
                for row in range(4)), "C4 cycle")
    return (
        state_order,
        subsets,
        arcs,
        ("both_way", 3, "one_way", 0, "missing", 3),
        ("P7", 7, 6, p7_rank, 6 - p7_rank),
        ("P4", 4, 3, p4_rank, 3 - p4_rank),
        ("C4_after_closure", 4, 4, c4_rank, 4 - c4_rank, cycle),
    )


def symmetry_record(kernels):
    arc_reverse = (1, 0, 3, 2, 5, 4)
    chamber_reflect = (5, 4, 3, 2, 1, 0)

    reverse_equal_entries = 0
    reverse_descending_addresses = 0
    reverse_failures = []
    for r0 in range(M.P):
        for r1 in range(M.P):
            equal_here = True
            failed_pairs = []
            for point in range(6):
                if kernels[point][r0][r1] == kernels[arc_reverse[point]][r0][r1]:
                    reverse_equal_entries += 1
                else:
                    equal_here = False
                    if point < arc_reverse[point]:
                        failed_pairs.append((point, arc_reverse[point]))
            reverse_descending_addresses += int(equal_here)
            if not equal_here:
                reverse_failures.append((r0, r1, tuple(failed_pairs)))

    chamber_equal_entries = 0
    chamber_full = True
    first_chamber_failure = None
    for point in range(6):
        for r0 in range(M.P):
            for r1 in range(M.P):
                left = kernels[point][r0][r1]
                right = kernels[chamber_reflect[point]][12 - r0][12 - r1]
                if left == right:
                    chamber_equal_entries += 1
                else:
                    chamber_full = False
                    if first_chamber_failure is None:
                        first_chamber_failure = (point, r0, r1, left, right)

    interior = tuple(index for index in range(M.P) if index not in (0, 6, 12))
    source_middle = tuple((r0, r1) for r0 in (3, 9) for r1 in interior)
    endpoint_residual = ((0, 11), (0, 12), (6, 5),
                         (6, 7), (12, 0), (12, 1))
    observed_failure_pairs = tuple((r0, r1) for r0, r1, _pairs in reverse_failures)
    require(set(observed_failure_pairs) == set(source_middle + endpoint_residual),
            ("arc-reversal failure decomposition", observed_failure_pairs))

    return (
        ("arc_reversal", arc_reverse, reverse_equal_entries, 6 * M.P * M.P,
         reverse_descending_addresses, M.P * M.P, tuple(reverse_failures),
         ("source_middle", source_middle),
         ("endpoint_residual", endpoint_residual)),
        ("chamber_reflection", chamber_reflect, chamber_equal_entries,
         6 * M.P * M.P, chamber_full, first_chamber_failure),
    )


def support_and_harmonic_record(cells):
    actual_supports = set()
    for state, u_values, _v_values, _mask, addresses in cells:
        if state is None:
            continue
        for root, active_root in enumerate(u_values):
            if not active_root:
                continue
            for r0 in range(M.P):
                support = tuple(r1 for r1, value in enumerate(addresses[root][r0])
                                if value)
                if support:
                    actual_supports.add(support)

    def harmonic_sum(support):
        return sum((Fraction(1, index + 1) for index in support), Fraction())

    actual_by_sum = defaultdict(list)
    for support in sorted(actual_supports):
        actual_by_sum[harmonic_sum(support)].append(support)
    actual_collisions = tuple(
        tuple(group) for _value, group in sorted(actual_by_sum.items())
        if len(group) > 1
    )

    all_by_sum = defaultdict(list)
    for mask in range(1 << M.P):
        support = tuple(index for index in range(M.P) if mask & (1 << index))
        all_by_sum[harmonic_sum(support)].append(support)
    all_collisions = tuple(group for group in all_by_sum.values() if len(group) > 1)
    first_collision = min(
        (tuple(group[:2]) for group in all_collisions),
        key=lambda pair: (sum(map(len, pair)), pair),
    )
    require(first_collision == ((1,), (2, 5)),
            ("canonical harmonic collision", first_collision))

    return (
        ("actual_supports", len(actual_supports), len(actual_by_sum),
         len(actual_collisions), max((len(group) for group in actual_collisions),
                                     default=1),
         digest_json(tuple(sorted(actual_supports)))),
        ("all_subsets_of_1_to_13", 1 << M.P, len(all_by_sum),
         len(all_collisions), max(map(len, all_collisions)), first_collision),
    )


def exception_record(cells):
    _kernels, records = M.source_ratio_kernels(cells)
    actual_groups = records[0][2]
    require(len(actual_groups) == 6, ("exception group count", actual_groups))
    r0_by_point = [None] * 6
    r1_exception_set = tuple(index for index in range(M.P) if index not in (0, 6, 12))
    for (point, pair, r0), exceptions in actual_groups:
        require(pair == M.POINTS[point], ("exception point", point, pair))
        require(tuple(r1 for r1, _size in exceptions) == r1_exception_set,
                ("exception r1 support", point, exceptions))
        r0_by_point[point] = r0
    r0_by_point = tuple(r0_by_point)
    drops = tuple(r0_by_point[index] - r0_by_point[index + 1]
                  for index in range(5))
    require(r0_by_point == (12, 12, 9, 3, 0, 0), r0_by_point)
    require(drops == (0, 3, 6, 3, 0), drops)
    require(all(r0_by_point[index] + r0_by_point[5 - index] == 12
                for index in range(6)), "exception reflection")
    return (
        r0_by_point,
        drops,
        ("drop_generating_polynomial", "3*x*(1+x)^2"),
        ("exception_r1", r1_exception_set, "boundary_survivors", (0, 6, 12)),
    )


def main() -> None:
    require(lf_sha256(PARENT_PATH) == PARENT_SHA256,
            ("parent source drift", lf_sha256(PARENT_PATH), PARENT_SHA256))
    cells, actual, source_digest, tensor_digest = build_actual_tensor()
    profile, kernels, profile_digest = extract_kernels(actual)

    chains = chain_record()
    kernels_record = kernel_record(kernels)
    symmetries = symmetry_record(kernels)
    supports = support_and_harmonic_record(cells)
    exceptions = exception_record(cells)

    block_rank = sum(record[2] for record in kernels_record)
    record = (
        PARENT_SHA256,
        M.EXPECTED_SEMANTIC_SHA256,
        source_digest,
        tensor_digest,
        profile_digest,
        chains,
        kernels_record,
        block_rank,
        symmetries,
        supports,
        exceptions,
    )
    semantic = digest_json(record)
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_PINNED":
        require(semantic == EXPECTED_SEMANTIC_SHA256,
                ("semantic drift", semantic, EXPECTED_SEMANTIC_SHA256))

    print("== pointed-six bidirected-P4 address-cocycle sidecar ==")
    print(f"parent=(source_sha256={PARENT_SHA256},semantic_sha256={M.EXPECTED_SEMANTIC_SHA256})")
    print(f"rebuild=(source_digest={source_digest},tensor_digest={tensor_digest},diagonal_profile_digest={profile_digest})")
    print(f"chain_dictionary={chains}")
    print(f"kernel_record=(point,(state,root),rank,right_eigen1_nullity,left_eigen1_nullity,support_sccs,row_support_histogram,distinct_rows,unit_column_sums,digest)={kernels_record}")
    print(f"six_kernel_block_rank={block_rank}/78;all_78_rows_sum_to_one=PASS")
    print(f"symmetry_and_descent={symmetries}")
    print(f"support_harmonic_scalar_loss={supports}")
    print(f"exception_boundary_profile={exceptions}")
    print(f"semantic_sha256={semantic}")
    print("interpretation=six pointed lines are the arcs of bidirected P4 and edges of alternating P7;the static chain has H1=0")
    print("closure_gate=adding the missing Boolean-square edge gives C4 with H1 dimension 1;no physical closure/current is asserted")
    print("address_gate=the six diagonal 13x13 kernels are row-sum-one finite-field 2-block weights,not positive probabilities or chronology")
    print("harmonic_gate=retaining a subset retains its indicator word;retaining only the reciprocal scalar sum is noninjective")
    print("scope=FINITE-EXACT representation sidecar on one owner base;no U_clock,no arrival ancestry,no physical current,no D5 bridge,no row exclusion,no LRC(14)")
    print("all exact checks passed")


if __name__ == "__main__":
    main()
