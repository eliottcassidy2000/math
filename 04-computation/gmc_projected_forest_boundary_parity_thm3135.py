#!/usr/bin/env python3
"""Exact controls for THM-3135 projected forest-boundary parity.

For K_n and 0 <= r <= n-2, the integer map under study is

    T_(n,r) = P_r partial_(r+1),

from canonically oriented (r+1)-edge forests to rank-r set partitions.
The proof in THM-3135 is general.  This companion independently checks every
rank for 3 <= n <= 7 over F_2 and exhausts the endpoint-swap identity on every
elementary one-point transfer in the same universe.
"""

from collections import defaultdict, deque
from functools import lru_cache
from itertools import combinations


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def canonical_partition(blocks):
    return tuple(sorted(
        (tuple(sorted(block)) for block in blocks),
        key=lambda block: (block[0], len(block), block),
    ))


@lru_cache(maxsize=None)
def set_partitions(vertex_count):
    answer = []

    def recurse(vertex, blocks):
        if vertex == vertex_count:
            answer.append(canonical_partition(blocks))
            return
        for block in blocks:
            block.append(vertex)
            recurse(vertex + 1, blocks)
            block.pop()
        blocks.append([vertex])
        recurse(vertex + 1, blocks)
        blocks.pop()

    recurse(0, [])
    return tuple(answer)


def component_partition(vertex_count, edges):
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


def canonical_forest(partition):
    edges = []
    for block in partition:
        root = block[0]
        edges.extend(tuple(sorted((root, vertex))) for vertex in block[1:])
    return tuple(sorted(edges))


def boundary_column_mod_two(vertex_count, forest, row_index):
    column = 0
    for position in range(len(forest)):
        face = forest[:position] + forest[position + 1:]
        partition = component_partition(vertex_count, face)
        require(partition in row_index, "boundary face escaped target rank")
        column ^= 1 << row_index[partition]
    return column


def matrix_rank_mod_two(vertex_count, rank):
    targets = tuple(
        partition for partition in set_partitions(vertex_count)
        if vertex_count - len(partition) == rank
    )
    row_index = {partition: index for index, partition in enumerate(targets)}
    pivots = {}
    columns = 0
    complete_edges = tuple(combinations(range(vertex_count), 2))

    for forest in combinations(complete_edges, rank + 1):
        if component_partition(vertex_count, forest) is None:
            continue
        column = boundary_column_mod_two(vertex_count, forest, row_index)
        while column:
            pivot = column.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = column
                break
            column ^= pivots[pivot]
        columns += 1

    return targets, row_index, len(pivots), columns


def elementary_transfers(partition):
    """Return distinct neighbours obtained by moving one point."""

    neighbours = set()
    blocks = [set(block) for block in partition]
    for donor_index, donor in enumerate(blocks):
        if len(donor) < 2:
            continue
        for point in donor:
            for target_index in range(len(blocks)):
                if target_index == donor_index:
                    continue
                moved = [set(block) for block in blocks]
                moved[donor_index].remove(point)
                moved[target_index].add(point)
                neighbours.add(canonical_partition(moved))
    return tuple(sorted(neighbours))


def swap_forests(partition, neighbour):
    """Construct the two columns isolating one transfer pair mod two."""

    old_blocks = [set(block) for block in partition]
    new_blocks = [set(block) for block in neighbour]
    moved = None
    donor = target = None
    for point in range(sum(map(len, old_blocks))):
        old_index = next(index for index, block in enumerate(old_blocks)
                         if point in block)
        new_index = next(index for index, block in enumerate(new_blocks)
                         if point in block)
        if old_blocks[old_index] != new_blocks[new_index]:
            candidate_donor = old_blocks[old_index]
            candidate_target = new_blocks[new_index]
            if (len(candidate_donor) >= 2
                    and point in candidate_target
                    and candidate_donor - {point} in new_blocks
                    and candidate_target - {point} in old_blocks):
                moved = point
                donor = candidate_donor
                target = candidate_target - {point}
                break
    require(moved is not None, "failed to recover elementary transfer")

    donor_remainder = donor - {moved}
    refinement_blocks = []
    consumed = (donor, target)
    for block in old_blocks:
        if block not in consumed:
            refinement_blocks.append(block)
    refinement_blocks.extend((donor_remainder, target, {moved}))
    refinement = canonical_partition(refinement_blocks)
    base = list(canonical_forest(refinement))
    left = min(donor_remainder)
    right = min(target)
    bridge = tuple(sorted((left, right)))
    base.append(bridge)
    base = tuple(sorted(base))

    attach_left = tuple(sorted((left, moved)))
    attach_right = tuple(sorted((right, moved)))
    first = tuple(sorted(base + (attach_left,)))
    second = tuple(sorted(base + (attach_right,)))
    return first, second


EXPECTED = {
    3: ((1, 1), (3, 2)),
    4: ((1, 1), (6, 5), (7, 7)),
    5: ((1, 1), (10, 9), (25, 25), (15, 14)),
    6: ((1, 1), (15, 14), (65, 65), (90, 89), (31, 31)),
    7: ((1, 1), (21, 20), (140, 140), (350, 349),
        (301, 301), (63, 62)),
}


records = {}
swap_checks = transfer_edges = 0
for vertex_count, expected_rows in EXPECTED.items():
    rank_records = []
    for rank, expected in enumerate(expected_rows):
        targets, row_index, matrix_rank, column_count = matrix_rank_mod_two(
            vertex_count, rank)
        require((len(targets), matrix_rank) == expected,
                "projected-boundary rank census drift")
        expected_rank = len(targets) if rank % 2 == 0 else len(targets) - 1
        require(matrix_rank == expected_rank, "parity rank law failed")

        if rank:
            adjacency = {partition: set() for partition in targets}
            seen_pairs = set()
            for partition in targets:
                for neighbour in elementary_transfers(partition):
                    require(neighbour in row_index,
                            "transfer escaped fixed-rank partition layer")
                    pair = tuple(sorted((partition, neighbour)))
                    if pair in seen_pairs:
                        continue
                    seen_pairs.add(pair)
                    adjacency[partition].add(neighbour)
                    adjacency[neighbour].add(partition)
                    first, second = swap_forests(partition, neighbour)
                    first_column = boundary_column_mod_two(
                        vertex_count, first, row_index)
                    second_column = boundary_column_mod_two(
                        vertex_count, second, row_index)
                    target_pair = ((1 << row_index[partition])
                                   ^ (1 << row_index[neighbour]))
                    require(first_column ^ second_column == target_pair,
                            "endpoint-swap identity failed")
                    swap_checks += 1

            start = targets[0]
            reached = {start}
            queue = deque([start])
            while queue:
                current = queue.popleft()
                for neighbour in adjacency[current]:
                    if neighbour not in reached:
                        reached.add(neighbour)
                        queue.append(neighbour)
            require(len(reached) == len(targets),
                    "elementary-transfer graph disconnected")
            transfer_edges += len(seen_pairs)

        rank_records.append((len(targets), matrix_rank, column_count))
    records[vertex_count] = tuple(rank_records)


print("projected_forest_boundary=exact_all_rank_parity_controls")
for vertex_count, rank_records in records.items():
    payload = ";".join(
        f"r{rank}:rows{rows}:rank{matrix_rank}:forests{columns}"
        for rank, (rows, matrix_rank, columns) in enumerate(rank_records)
    )
    print(f"K{vertex_count}={payload}")
print(f"elementary_transfer_edges={transfer_edges}")
print(f"endpoint_swap_identities={swap_checks}")
print("even_rank=image_full;odd_rank=image_augmentation_zero")
print("general_proof=swap_pair_span+fixed_rank_transfer_connectivity+parity")
print("all_exact_checks=PASS")
