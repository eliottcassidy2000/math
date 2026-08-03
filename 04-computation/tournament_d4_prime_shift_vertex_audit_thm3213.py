#!/usr/bin/env python3
"""Independent vertex-level audit of the U_2 prime-shift constant."""

from itertools import permutations
from math import factorial


def require(condition, payload):
    if not condition:
        raise RuntimeError(payload)


def build_u2():
    # Blocks 0,1,2; each internal T_2 is 0 -> 1; quotient is 0 -> 1 -> 2 -> 0.
    vertices = tuple(range(6))

    def block(v):
        return v // 2

    def local(v):
        return v % 2

    def arc(v, w):
        if block(v) == block(w):
            return local(v) < local(w)
        return (block(v), block(w)) in ((0, 1), (1, 2), (2, 0))

    return vertices, arc


def set_partitions(items):
    if not items:
        yield ()
        return
    first, rest = items[0], items[1:]
    for partition in set_partitions(rest):
        yield ((first,),) + partition
        for index in range(len(partition)):
            block = tuple(sorted((first,) + partition[index]))
            yield partition[:index] + (block,) + partition[index + 1 :]


def canonical_partition(partition):
    return tuple(sorted((tuple(sorted(block)) for block in partition), key=lambda b: b[0]))


def directed_path_count(block, arc):
    answer = 0
    for order in permutations(block):
        if all(arc(order[i], order[i + 1]) for i in range(len(order) - 1)):
            answer += 1
    return answer


vertices, arc = build_u2()
partitions = {canonical_partition(p) for p in set_partitions(vertices)}
unordered_profile = [0] * 7
partition_counts = [0] * 7
for partition in partitions:
    d = len(partition)
    partition_counts[d] += 1
    contribution = 1
    for block in partition:
        contribution *= directed_path_count(block, arc)
    unordered_profile[d] += contribution

ordered_profile = tuple(factorial(d) * unordered_profile[d] for d in range(1, 7))
require(tuple(partition_counts[1:]) == (1, 31, 90, 65, 15, 1), partition_counts)
require(tuple(unordered_profile[1:]) == (45, 171, 186, 81, 15, 1), unordered_profile)
require(ordered_profile == (45, 342, 1116, 1944, 1800, 720), ordered_profile)

print("U2 DIRECT VERTEX PATH-COVER AUDIT")
print(f"set_partition_counts={tuple(partition_counts[1:])}")
print(f"unordered_path_cover_profile={tuple(unordered_profile[1:])}")
print(f"ordered_factorial_profile={ordered_profile}")
print("d4_constant=1944=2^3*3^5")
print("implementation=canonical set partitions plus direct induced-path permutations")
