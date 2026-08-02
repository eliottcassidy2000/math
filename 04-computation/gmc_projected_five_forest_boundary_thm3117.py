#!/usr/bin/env python3
"""Exact projected forest-boundary surjectivity for THM-3117.

For N=6,7,8,9, form the integer matrix

    T_N = P_4 partial_5,

whose columns are five-edge forests and whose rows are rank-four component
partitions.  Reduction modulo two turns every boundary sign into one.  A
deterministic bitset Gaussian elimination proves full row rank.  Full rank
over F_2 implies full row rank over Q for the original signed integer matrix.
"""

from collections import defaultdict
from functools import lru_cache
from itertools import combinations


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


@lru_cache(maxsize=None)
def set_partitions(items):
    answer = []

    def recurse(position, blocks):
        if position == len(items):
            answer.append(tuple(tuple(block) for block in blocks))
            return
        item = items[position]
        for index in range(len(blocks)):
            blocks[index].append(item)
            recurse(position + 1, blocks)
            blocks[index].pop()
        blocks.append([item])
        recurse(position + 1, blocks)
        blocks.pop()

    recurse(0, [])
    return tuple(answer)


def canonical_partition(blocks):
    return tuple(
        sorted(
            (tuple(sorted(block)) for block in blocks),
            key=lambda block: (block[0], len(block), block),
        )
    )


def component_partition(vertex_count, edges):
    """Return the component partition, or None if the edge set is cyclic."""

    parent = list(range(vertex_count))

    def find(vertex):
        while parent[vertex] != vertex:
            parent[vertex] = parent[parent[vertex]]
            vertex = parent[vertex]
        return vertex

    for left, right in edges:
        root_left = find(left)
        root_right = find(right)
        if root_left == root_right:
            return None
        parent[root_right] = root_left

    blocks = defaultdict(list)
    for vertex in range(vertex_count):
        blocks[find(vertex)].append(vertex)
    return canonical_partition(blocks.values())


def projected_boundary_rank_mod_two(vertex_count):
    """Return target size, exact F_2 rank, and columns read to full rank."""

    target_partitions = tuple(
        partition
        for partition in set_partitions(tuple(range(vertex_count)))
        if vertex_count - len(partition) == 4
    )
    row_index = {
        partition: index for index, partition in enumerate(target_partitions)
    }
    target_size = len(target_partitions)
    pivots = {}
    columns_read = 0
    complete_edges = tuple(combinations(range(vertex_count), 2))

    for forest in combinations(complete_edges, 5):
        if component_partition(vertex_count, forest) is None:
            continue

        # Modulo two, all five alternating boundary signs are one.  The five
        # deleted forests have distinct edge sets; repeated component flats,
        # if any, cancel with their exact parity in this bit vector.
        column = 0
        for position in range(5):
            face = forest[:position] + forest[position + 1 :]
            partition = component_partition(vertex_count, face)
            require(partition in row_index, "boundary face escaped rank four")
            column ^= 1 << row_index[partition]

        while column:
            pivot = column.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = column
                break
            column ^= pivots[pivot]

        columns_read += 1
        if len(pivots) == target_size:
            break

    return target_size, len(pivots), columns_read


EXPECTED = {
    6: (31, 31, 433),
    7: (301, 301, 3187),
    8: (1701, 1701, 13683),
    9: (6951, 6951, 44524),
}


records = {}
for vertex_count in EXPECTED:
    record = projected_boundary_rank_mod_two(vertex_count)
    require(record == EXPECTED[vertex_count], "projected boundary rank drift")
    records[vertex_count] = record


print("projected_five_forest_boundary=exact_mod2_rank_certificate")
for vertex_count, (rows, rank, columns) in records.items():
    print(
        f"K{vertex_count}=rank4_partitions:{rows}:rank_F2:{rank}:"
        f"lex_five_forests_to_full_rank:{columns}"
    )
print("actual_product_gamma_spaces=K8,K9:full_row_rank_over_F2")
print("integer_signed_matrix=full_row_rank_over_Q:odd_maximal_minor")
print("every_rational_rank4_current=projection_of_rational_forest_cycle")
print("all_exact_checks=PASS")
