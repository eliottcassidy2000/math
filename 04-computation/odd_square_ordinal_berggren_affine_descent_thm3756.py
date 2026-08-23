#!/usr/bin/env python3
"""Exact audit for THM-3756.

Universe and controls
---------------------
* triangular fold: -4000 <= z <= 4000 and 1 <= h <= 9;
* every odd-square ordinal pair (r,s) with 2 <= r <= 600;
* every Berggren node through depth 10 (88,573 nodes in levels 0 through 10);
* ambient and compressed natural-number addresses through shell r=600.

The global statements in THM-3756 are proved by the displayed algebra and
strict descent.  This companion checks constants, branch conventions,
consequence objects, equality boundaries, and hostile quotients.  It uses a
runtime ``require`` rather than ``assert`` so ``python`` and ``python -O``
exercise the same gates.
"""

from __future__ import annotations

from bisect import bisect_right
from collections import Counter, defaultdict
from functools import lru_cache
from hashlib import sha256
from math import gcd, isqrt
import sys


MAX_Z = 4_000
MAX_H = 9
MAX_R = 600
TREE_DEPTH = 10


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def triangular(z: int) -> int:
    return z * (z + 1) // 2


def centered_triangular_difference(z: int, h: int) -> int:
    return triangular(z + h) - triangular(z - h)


def odd_root(rank: int) -> int:
    return 2 * rank - 1


def admissible(rank: int, inner_rank: int) -> bool:
    return (
        rank > inner_rank >= 1
        and gcd(odd_root(rank), odd_root(inner_rank)) == 1
    )


def cone_pair(rank: int, inner_rank: int) -> bool:
    return rank > inner_rank >= 1


def ambient_triple(rank: int, inner_rank: int) -> tuple[int, int, int]:
    require(cone_pair(rank, inner_rank), "pair outside the triangular cone")
    q = odd_root(rank)
    d = odd_root(inner_rank)
    return q * d, (q * q - d * d) // 2, (q * q + d * d) // 2


def ordinal_triple(rank: int, inner_rank: int) -> tuple[int, int, int]:
    """Return (odd leg, even leg, hypotenuse)."""
    require(admissible(rank, inner_rank), "inadmissible ordinal pair")
    return ambient_triple(rank, inner_rank)


def ordinal_from_triple(a: int, b: int, c: int) -> tuple[int, int]:
    require(a > 0 and b > 0 and c > 0, "triple must be positive")
    require(b % 2 == 0, "B must be the ordered even leg")
    require(a * a + b * b == c * c, "not Pythagorean")
    q = isqrt(b + c)
    d = isqrt(c - b)
    require(q * q == b + c, "B+C is not a square")
    require(d * d == c - b, "C-B is not a square")
    require(q % 2 == 1 and d % 2 == 1, "square roots must be odd")
    pair = ((q + 1) // 2, (d + 1) // 2)
    require(ordinal_triple(*pair) == (a, b, c), "inverse mismatch")
    return pair


def parameters(rank: int, inner_rank: int) -> tuple[int, int]:
    """THM-3357 convention: 0 < m < n and Psi=(n^2-m^2,2mn,n^2+m^2)."""
    return rank - inner_rank, rank + inner_rank - 1


def ordinal_from_parameters(m: int, n: int) -> tuple[int, int]:
    require(0 < m < n, "parameters outside the positive chamber")
    require((m + n) % 2 == 1, "parameters must have opposite parity")
    return (m + n + 1) // 2, (n - m + 1) // 2


def parameter_triple(m: int, n: int) -> tuple[int, int, int]:
    return n * n - m * m, 2 * m * n, n * n + m * m


def children(pair: tuple[int, int]) -> dict[str, tuple[int, int]]:
    r, s = pair
    require(cone_pair(r, s), "children need a positive-cone parent")
    return {
        "L": (r + 2 * s - 1, s),
        "M": (2 * r + s - 1, r),
        "R": (2 * r - s, r),
    }


def parameter_children(m: int, n: int) -> dict[str, tuple[int, int]]:
    return {
        "L": (n, 2 * n - m),
        "M": (n, 2 * n + m),
        "R": (m, 2 * m + n),
    }


def parent(pair: tuple[int, int]) -> tuple[str, tuple[int, int]] | None:
    """Return the intrinsic branch label and strict predecessor."""
    r, s = pair
    require(cone_pair(r, s), "parent needs a positive-cone node")
    if r == 3 * s - 1:
        return None
    if r >= 3 * s:
        result = "L", (r - 2 * s + 1, s)
    elif 2 * s <= r <= 3 * s - 2:
        result = "M", (s, r - 2 * s + 1)
    elif s < r <= 2 * s - 1:
        result = "R", (s, 2 * s - r)
    else:
        raise AssertionError(f"uncovered primitive cone boundary: {pair}")
    require(cone_pair(*result[1]), "descent left the positive ordinal cone")
    require(result[1][0] < r, "descent did not lower the outer rank")
    require(children(result[1])[result[0]] == pair, "inverse branch mismatch")
    return result


@lru_cache(maxsize=None)
def totient(n: int) -> int:
    require(n >= 1, "totient domain")
    result = n
    p = 2
    remaining = n
    while p * p <= remaining:
        if remaining % p == 0:
            while remaining % p == 0:
                remaining //= p
            result -= result // p
        p += 1 if p == 2 else 2
    if remaining > 1:
        result -= result // remaining
    return result


@lru_cache(maxsize=None)
def admissible_fibre(rank: int) -> tuple[int, ...]:
    return tuple(s for s in range(1, rank) if admissible(rank, s))


@lru_cache(maxsize=None)
def compressed_offset(rank: int) -> int:
    return sum(totient(2 * k - 1) // 2 for k in range(2, rank))


def ambient_address(rank: int, inner_rank: int) -> int:
    require(rank > inner_rank >= 1, "ambient address domain")
    return triangular(rank - 2) + inner_rank


def selected_index(rank: int, inner_rank: int) -> int:
    require(admissible(rank, inner_rank), "selected index domain")
    return bisect_right(admissible_fibre(rank), inner_rank)


def compressed_address(rank: int, inner_rank: int) -> int:
    require(admissible(rank, inner_rank), "compressed address domain")
    return compressed_offset(rank) + selected_index(rank, inner_rank)


def component_content(pair: tuple[int, int]) -> int:
    return gcd(odd_root(pair[0]), odd_root(pair[1]))


def component_root(pair: tuple[int, int]) -> tuple[int, int]:
    g = component_content(pair)
    return (3 * g + 1) // 2, (g + 1) // 2


def replay_word(upward_path: list[str], root: tuple[int, int] = (2, 1)) -> tuple[int, int]:
    node = root
    for branch in reversed(upward_path):
        node = children(node)[branch]
    return node


def main() -> None:
    # Make the stored transcript a literal byte match on every platform.
    sys.stdout.reconfigure(newline="\n")

    # The integer triangular fold and its orientation-sensitive sidecar.
    triangular_checks = 0
    for z in range(-MAX_Z, MAX_Z + 1):
        require(triangular(z) >= 0, "integer triangular value is negative")
        require(triangular(z) == triangular(-z - 1), "fold reflection failed")
        for h in range(1, MAX_H + 1):
            require(
                centered_triangular_difference(z, h) == h * (2 * z + 1),
                "centered triangular difference failed",
            )
            triangular_checks += 1

    fold_iff_checks = 0
    for n in range(0, MAX_R + 1):
        for z in range(-MAX_Z, MAX_Z + 1):
            require(
                2 * (triangular(z) - triangular(n))
                == (z - n) * (z + n + 1),
                "triangular collision factorization failed",
            )
            require(
                (triangular(z) == triangular(n))
                == (z == n or z == -n - 1),
                "triangular fold fibre failed",
            )
            fold_iff_checks += 1

    odd_square_checks = 0
    for r in range(1, MAX_R + 1):
        extracted_root = centered_triangular_difference(r - 1, 2) // 2
        require(
            extracted_root == odd_root(r),
            "odd-root extraction failed",
        )
        require(
            extracted_root * extracted_root == odd_root(r) ** 2,
            "odd-square ordinal extraction failed",
        )
        odd_square_checks += 1

    # The two-ordinal chart, exact fibres, natural-number addresses, and
    # parameter/Berggren branch convention.
    pairs: list[tuple[int, int]] = []
    fibre_sizes: dict[int, int] = {}
    primitive_ambient: set[int] = set()
    compressed: set[int] = set()
    max_descent = 0
    branch_checks = 0
    for r in range(2, MAX_R + 1):
        fibre = list(admissible_fibre(r))
        fibre_sizes[r] = len(fibre)
        require(len(fibre) == totient(odd_root(r)) // 2, "totient fibre failed")
        require(fibre[0] == 1, "canonical U-spine section missing")
        for s in fibre:
            pair = (r, s)
            pairs.append(pair)
            a, b, c = ordinal_triple(r, s)
            require(a * a + b * b == c * c, "ordinal triple identity failed")
            require(gcd(gcd(a, b), c) == 1, "ordinal triple is not primitive")
            require(b + c == odd_root(r) ** 2, "outer odd square failed")
            require(c - b == odd_root(s) ** 2, "inner odd square failed")
            require(ordinal_from_triple(a, b, c) == pair, "ordinal inverse failed")
            m, n = parameters(r, s)
            require(gcd(m, n) == 1 and (m + n) % 2 == 1, "parameter boundary failed")
            require(parameter_triple(m, n) == (a, b, c), "parameter chart mismatch")
            require(ordinal_from_parameters(m, n) == pair, "parameter inverse failed")

            primitive_ambient.add(ambient_address(r, s))
            compressed.add(compressed_address(r, s))

            ordinal_children = children(pair)
            raw_children = parameter_children(m, n)
            for label in "LMR":
                child = ordinal_children[label]
                require(admissible(*child), "branch left the ordinal domain")
                require(parameters(*child) == raw_children[label], "branch convention mismatch")
                require(parent(child) == (label, pair), "branch/parent round trip failed")
                branch_checks += 1

            ranks = {label: child[0] for label, child in ordinal_children.items()}
            require(ranks["M"] > ranks["L"], "middle child must outrank left")
            require(ranks["M"] > ranks["R"], "middle child must outrank right")
            require(
                ranks["L"] - ranks["R"] == 3 * s - r - 1,
                "outer-child comparison failed",
            )

            upward: list[str] = []
            cursor = pair
            while cursor != (2, 1):
                step = parent(cursor)
                require(step is not None, "nonroot has no parent")
                label, cursor = step
                upward.append(label)
            require(replay_word(upward) == pair, "descent word failed to replay")
            max_descent = max(max_descent, len(upward))

    require(len(pairs) == len(set(pairs)), "ordinal pairs duplicated")
    total_compressed = sum(fibre_sizes.values())
    require(
        compressed == set(range(1, total_compressed + 1)),
        "compressed selected addresses are not contiguous",
    )
    all_ambient = set(range(1, triangular(MAX_R - 1) + 1))
    first_ambient_hole = min(all_ambient - primitive_ambient)
    require(first_ambient_hole == 8, "first nonprimitive ambient hole moved")

    # The entire triangular cone is a forest of odd-square-scaled primitive
    # Berggren trees.  Branches preserve the odd root content, and each
    # component terminates on its unique r=3s-1 root.
    forest_branch_checks = 0
    forest_max_descent = 0
    component_roots: set[tuple[int, int]] = set()
    for r in range(2, MAX_R + 1):
        for s in range(1, r):
            pair = (r, s)
            g = component_content(pair)
            a, b, c = ambient_triple(r, s)
            require(a * a + b * b == c * c, "ambient triple identity failed")
            require(gcd(gcd(a, b), c) == g * g, "ambient content law failed")
            for label, child in children(pair).items():
                require(component_content(child) == g, "branch changed component content")
                require(parent(child) == (label, pair), "forest branch inverse failed")
                forest_branch_checks += 1
            upward: list[str] = []
            cursor = pair
            while True:
                step = parent(cursor)
                if step is None:
                    break
                label, cursor = step
                upward.append(label)
            root = component_root(pair)
            require(cursor == root, "forest descent found the wrong component root")
            require(replay_word(upward, root) == pair, "forest word failed to replay")
            component_roots.add(root)
            forest_max_descent = max(forest_max_descent, len(upward))

    # The boundary r=3s-1 is primitive only at the content-one root.
    for s in range(1, MAX_R // 3 + 1):
        r = 3 * s - 1
        require(
            admissible(r, s) == ((r, s) == (2, 1)),
            "root boundary classification failed",
        )

    # Exact tree breadth, uniqueness, and rank collisions through depth ten.
    levels: list[set[tuple[int, int]]] = [{(2, 1)}]
    seen = {(2, 1)}
    rank_histograms: list[Counter[int]] = [Counter({2: 1})]
    for depth in range(1, TREE_DEPTH + 1):
        level: set[tuple[int, int]] = set()
        for node in levels[-1]:
            for label, child in children(node).items():
                require(parent(child) == (label, node), "tree inverse failed")
                require(child not in seen, "Berggren child duplicated an earlier node")
                require(child not in level, "Berggren siblings/parents collided")
                level.add(child)
        require(len(level) == 3**depth, "ternary level size failed")
        seen.update(level)
        levels.append(level)
        rank_histograms.append(Counter(r for r, _ in level))

    depth_equality_nodes: list[list[tuple[int, int]]] = []
    for depth, level in enumerate(levels):
        require(all(depth <= r - 2 for r, _ in level), "depth/rank bound failed")
        equality = sorted(pair for pair in level if pair[0] == depth + 2)
        expected = (
            [(2, 1)]
            if depth == 0
            else [(depth + 2, 1), (depth + 2, depth + 1)]
        )
        require(equality == expected, "depth/rank equality boundary failed")
        depth_equality_nodes.append(equality)

    # Hostiles: the coarse outer rank collides immediately and is neither
    # depth nor child-transition sufficient.
    rank_fibres: dict[int, list[tuple[int, int]]] = defaultdict(list)
    for pair in pairs:
        rank_fibres[pair[0]].append(pair)
    first_collision_rank = min(r for r, fibre in rank_fibres.items() if len(fibre) > 1)
    require(first_collision_rank == 3, "minimal outer-rank collision moved")
    rank_three = [ordinal_triple(*pair) for pair in rank_fibres[3]]
    require(rank_three == [(5, 12, 13), (15, 8, 17)], "rank-three hostile moved")
    signatures = {
        pair: tuple(children(pair)[label][0] for label in "LMR")
        for pair in rank_fibres[3]
    }
    require(len(set(signatures.values())) == 2, "outer rank unexpectedly determines children")
    require(
        any(4 in histogram for histogram in rank_histograms[1:2])
        and any(4 in histogram for histogram in rank_histograms[2:3]),
        "outer rank unexpectedly determines depth",
    )
    hostile = ambient_triple(5, 2)
    require(hostile == (27, 36, 45), "nonprimitive hostile construction moved")
    require(
        hostile[0] ** 2 + hostile[1] ** 2 == hostile[2] ** 2
        and gcd(gcd(*hostile[:2]), hostile[2]) == 9,
        "nonprimitive hostile moved",
    )

    # Prime outer roots have the maximal r-1 overlap fibre.
    prime_shells = [
        r
        for r in range(2, MAX_R + 1)
        if totient(odd_root(r)) == odd_root(r) - 1
    ]
    for r in prime_shells:
        require(fibre_sizes[r] == r - 1, "prime shell is not full")

    semantic = {
        "triangular_checks": triangular_checks,
        "fold_iff_checks": fold_iff_checks,
        "odd_square_checks": odd_square_checks,
        "pair_count": len(pairs),
        "branch_checks": branch_checks,
        "forest_branch_checks": forest_branch_checks,
        "forest_max_descent": forest_max_descent,
        "component_roots": sorted(component_roots)[:20],
        "max_descent": max_descent,
        "tree_level_sizes": [len(level) for level in levels],
        "rank_histograms": [sorted(hist.items()) for hist in rank_histograms[:6]],
        "depth_equality_nodes": depth_equality_nodes,
        "rank_three": rank_three,
        "rank_three_child_signatures": sorted(signatures.items()),
        "fibre_sizes_2_20": [fibre_sizes[r] for r in range(2, 21)],
        "first_ambient_hole": first_ambient_hole,
        "compressed_total": total_compressed,
        "prime_shells_2_50": [r for r in prime_shells if r <= 50],
    }
    digest = sha256(repr(semantic).encode("ascii")).hexdigest()

    print("THM-3756 exact audit")
    print(f"triangular_checks={triangular_checks}")
    print(f"fold_iff_checks={fold_iff_checks}")
    print(f"odd_square_checks={odd_square_checks}")
    print(f"ordinal_pairs_r_le_{MAX_R}={len(pairs)}")
    print(f"branch_round_trips={branch_checks}")
    print(f"maximum_descent_length={max_descent}")
    print(f"forest_branch_round_trips={forest_branch_checks}")
    print(f"forest_maximum_descent_length={forest_max_descent}")
    print(f"forest_component_roots_r_le_{MAX_R}={len(component_roots)}")
    print(f"tree_level_{TREE_DEPTH}_nodes={len(levels[-1])}")
    print(f"depth_rank_equality_two_boundary_rays_through={TREE_DEPTH}")
    print(f"first_coarse_rank_collision={first_collision_rank}: {rank_three}")
    print(f"rank_three_child_signatures={signatures}")
    print(f"first_ambient_primitivity_hole={first_ambient_hole}")
    print(f"compressed_addresses=1..{total_compressed}")
    print(f"prime_shells_through_50={[r for r in prime_shells if r <= 50]}")
    print(f"semantic_sha256={digest}")
    print("PASS")


if __name__ == "__main__":
    main()
