#!/usr/bin/env python3
"""THM-3202 independent vertex audit with no transfer-kernel imports."""

from array import array
from functools import lru_cache


def require(condition, payload):
    if not condition:
        raise AssertionError(payload)


def transitive(n):
    return tuple(tuple(i < j for j in range(n)) for i in range(n))


def directed_triangle():
    return (
        (False, True, False),
        (False, False, True),
        (True, False, False),
    )


def order_join(left, right):
    a, b = len(left), len(right)
    out = [[False] * (a + b) for _ in range(a + b)]
    for i in range(a):
        for j in range(a):
            out[i][j] = left[i][j]
    for i in range(b):
        for j in range(b):
            out[a + i][a + j] = right[i][j]
    for i in range(a):
        for j in range(b):
            out[i][a + j] = True
    return tuple(tuple(row) for row in out)


def repeated_join(base, r):
    require(r >= 1, r)
    result = base
    for _ in range(r - 1):
        result = order_join(result, base)
    return result


def substitute(quotient, blocks):
    offsets = []
    total = 0
    for block in blocks:
        offsets.append(total)
        total += len(block)
    out = [[False] * total for _ in range(total)]
    for block_index, block in enumerate(blocks):
        offset = offsets[block_index]
        for i in range(len(block)):
            for j in range(len(block)):
                out[offset + i][offset + j] = block[i][j]
    for i in range(len(blocks)):
        for j in range(len(blocks)):
            if i == j or not quotient[i][j]:
                continue
            for u in range(len(blocks[i])):
                for v in range(len(blocks[j])):
                    out[offsets[i] + u][offsets[j] + v] = True
    return tuple(tuple(row) for row in out)


def c3_equal_lift(block):
    c3 = directed_triangle()
    return substitute(c3, (block, block, block))


def held_karp_subset_counts(tournament):
    """Return H(T[S]) for every subset using a flat unsigned-64 endpoint DP."""
    n = len(tournament)
    size = 1 << n
    incoming = [0] * n
    for v in range(n):
        incoming[v] = sum(1 << u for u in range(n) if tournament[u][v])
    ending = array("Q", [0]) * (size * n)
    for v in range(n):
        ending[(1 << v) * n + v] = 1
    hamilton = [0] * size
    for mask in range(1, size):
        bits = mask
        total = 0
        while bits:
            bit = bits & -bits
            v = bit.bit_length() - 1
            previous = mask ^ bit
            if previous:
                predecessors = previous & incoming[v]
                value = 0
                while predecessors:
                    predecessor_bit = predecessors & -predecessors
                    u = predecessor_bit.bit_length() - 1
                    value += ending[previous * n + u]
                    predecessors ^= predecessor_bit
                ending[mask * n + v] = value
            else:
                value = 1
            total += value
            bits ^= bit
        hamilton[mask] = total
    return tuple(hamilton)


def hamilton_count(tournament):
    return held_karp_subset_counts(tournament)[-1]


def path_cover_profile(tournament):
    n = len(tournament)
    hamilton = held_karp_subset_counts(tournament)

    @lru_cache(maxsize=None)
    def covers(mask):
        if mask == 0:
            return (1,)
        least = mask & -mask
        remainder = mask ^ least
        result = [0] * (mask.bit_count() + 1)
        sub = remainder
        while True:
            component = sub | least
            rest = covers(mask ^ component)
            weight = hamilton[component]
            for d, count in enumerate(rest):
                result[d + 1] += weight * count
            if sub == 0:
                break
            sub = (sub - 1) & remainder
        return tuple(result)

    return covers((1 << n) - 1)


def main():
    c3 = directed_triangle()

    # Independent full profiles at total lift order at most nine.
    profile_controls = (
        ("K1", 1, c3_equal_lift(transitive(1)), (3, 3, 1)),
        ("K1", 2, c3_equal_lift(transitive(2)), (45, 171, 186)),
        ("K1", 3, c3_equal_lift(transitive(3)), (2721, 18345, 37135)),
        ("C3", 1, c3_equal_lift(c3), (3159, 21303, 42201)),
    )
    profile_checks = 0
    for seed, r, tournament, expected in profile_controls:
        profile = path_cover_profile(tournament)
        require(
            tuple(profile[1:4]) == expected,
            (seed, r, tuple(profile[1:4]), expected),
        )
        profile_checks += 3

    # Hamiltonian controls extend to 12 vertices for K1 and 18 for C3.
    k1_expected = (3, 45, 2721, 421425)
    k1_values = []
    for r, expected in enumerate(k1_expected, start=1):
        value = hamilton_count(c3_equal_lift(transitive(r)))
        require(value == expected, ("K1", r, value, expected))
        k1_values.append(value)

    c3_expected = (3159, 82069875945)
    c3_values = []
    for r, expected in enumerate(c3_expected, start=1):
        block = repeated_join(c3, r)
        value = hamilton_count(c3_equal_lift(block))
        require(value == expected, ("C3", r, value, expected))
        c3_values.append(value)

    print("C3 REPEATED-JOIN MOVING JET -- THM-3202 INDEPENDENT AUDIT")
    print(
        "status=VERIFIED-EXACT;"
        "no_transfer_kernel_or_finite_difference_imports=true"
    )
    print("engine=flat_Held_Karp_endpoints_plus_canonical_set_partitions")
    print(
        f"profile_controls={len(profile_controls)};"
        f"d=1..3;checks={profile_checks};PASS"
    )
    print(f"K1_H_r1_to_r4={tuple(k1_values)};max_vertices=12;PASS")
    print(f"C3_H_r1_to_r2={tuple(c3_values)};max_vertices=18;PASS")
    print("cyclic_cube_control=3159;vertex_level=PASS")


if __name__ == "__main__":
    main()
