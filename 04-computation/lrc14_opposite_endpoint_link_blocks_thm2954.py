#!/usr/bin/env python3
"""Exact controls for the proposed THM-2954 endpoint-link theorem."""

from fractions import Fraction as F
from math import gcd


N = 7
LIMIT = 200


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def valuation(n: int, prime: int) -> int:
    answer = 0
    while n % prime == 0:
        answer += 1
        n //= prime
    return answer


def unit14(n: int) -> int:
    return (
        n // (2 ** valuation(n, 2) * 7 ** valuation(n, 7))
    ) % 14


def opposition_pair(residue: int) -> tuple[int, int]:
    return tuple(sorted((residue, (-residue) % 14)))


def opposite_link(u: int, v: int) -> bool:
    return (u + v) % (14 * gcd(u, v)) == 0


def block_link(u: int, v: int) -> bool:
    return (
        valuation(u, 2) == valuation(v, 2)
        and valuation(u, 7) == valuation(v, 7)
        and (unit14(u) + unit14(v)) % 14 == 0
    )


def bipartite(adjacency: list[list[int]]) -> bool:
    colours: list[int | None] = [None] * len(adjacency)
    for start in range(len(adjacency)):
        if colours[start] is not None:
            continue
        colours[start] = 0
        stack = [start]
        while stack:
            i = stack.pop()
            for j, edge in enumerate(adjacency[i]):
                if not edge:
                    continue
                if colours[j] is None:
                    colours[j] = 1 - colours[i]  # type: ignore[operator]
                    stack.append(j)
                elif colours[j] == colours[i]:
                    return False
    return True


def circle_norm(x: F) -> F:
    residue = x % 1
    return min(residue, 1 - residue)


def danger(speed: int, x: F) -> bool:
    return circle_norm(speed * x) < F(1, 14)


link_count = 0
for u in range(1, LIMIT + 1):
    for v in range(1, LIMIT + 1):
        observed = opposite_link(u, v)
        require(
            observed == block_link(u, v),
            f"block classification failed at {(u, v)}",
        )
        link_count += int(observed)

# Verify the complete-bipartite block description on the whole finite bank.
block_count = 0
signatures = {
    (
        valuation(u, 2),
        valuation(u, 7),
        opposition_pair(unit14(u)),
    )
    for u in range(1, LIMIT + 1)
}
for a, b, pair in sorted(signatures):
    left = [
        u
        for u in range(1, LIMIT + 1)
        if valuation(u, 2) == a
        and valuation(u, 7) == b
        and unit14(u) == pair[0]
    ]
    right = [
        u
        for u in range(1, LIMIT + 1)
        if valuation(u, 2) == a
        and valuation(u, 7) == b
        and unit14(u) == pair[1]
    ]
    if not left or not right:
        continue
    block_count += 1
    for u in left:
        for v in right:
            require(opposite_link(u, v), f"missing cross edge at {(u, v)}")
    for side in (left, right):
        for i, u in enumerate(side):
            for v in side[i:]:
                require(not opposite_link(u, v), f"illegal same-side edge {(u, v)}")

# A directed C7-circulant zero-diagonal support is a subset of the six
# nonzero differences.  Every nonempty support contains a translated
# seven-cycle after orientation is forgotten, and therefore is not bipartite.
for mask in range(1, 1 << (N - 1)):
    adjacency = [[0] * N for _ in range(N)]
    for d in range(1, N):
        if mask & (1 << (d - 1)):
            for i in range(N):
                j = (i + d) % N
                adjacency[i][j] = 1
                adjacency[j][i] = 1
    require(not bipartite(adjacency), f"circulant support {mask} is bipartite")

# Sharp seam: D_1 hands to D_13 almost everywhere at x=1/14, while both
# strict-open sets omit the seam itself.
seam = F(1, 14)
epsilon = F(1, 10000)
require(danger(1, seam - epsilon), "left owner absent before seam")
require(not danger(13, seam - epsilon), "right owner entered before seam")
require(not danger(1, seam), "D_1 contains its strict endpoint")
require(not danger(13, seam), "D_13 contains its strict endpoint")
require(not danger(1, seam + epsilon), "left owner persists past seam")
require(danger(13, seam + epsilon), "right owner absent after seam")

print("THM-2954 OPPOSITE ENDPOINT-LINK CONTROL")
print(f"ordered_pairs={LIMIT * LIMIT};links={link_count};classification=PASS")
print(f"nonempty_complete_bipartite_blocks={block_count};block_census=PASS")
print("nonempty_directed_c7_supports=63;all_nonbipartite=PASS")
print("strict_seam=(1,13,1/14);a.e._switch=yes;pointwise_cover=no")
print("all_exact_checks=PASS")
