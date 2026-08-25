"""Independent Calkin--Wilf branch/fibre audit for THM-4057."""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations
from math import gcd


def require(flag: bool, payload: object) -> None:
    if not flag:
        raise AssertionError(payload)


def left(pair: tuple[int, int]) -> tuple[int, int]:
    p, q = pair
    return p, p + q


def right(pair: tuple[int, int]) -> tuple[int, int]:
    p, q = pair
    return p + q, q


def reverse(pair: tuple[int, int]) -> tuple[int, int]:
    return pair[1], pair[0]


def leg_swap(x: Fraction) -> Fraction:
    return (1 - x) / (1 + x)


def normalized_unordered_triple(p: int, q: int) -> tuple[int, int, int]:
    require(0 < p < q and gcd(p, q) == 1, (p, q))
    content = 2 if p % 2 and q % 2 else 1
    legs = sorted(((q * q - p * p) // content, 2 * p * q // content))
    return legs[0], legs[1], (q * q + p * p) // content


def children(m: int, n: int) -> tuple[tuple[str, int, int], ...]:
    return (("A", 2 * m - n, m), ("B", 2 * m + n, m), ("C", m + 2 * n, n))


def cw_pair(word: str) -> tuple[int, int]:
    p, q = 1, 1
    for letter in word:
        if letter == "L":
            q += p
        else:
            p += q
    return p, q


def word_index(word: str) -> int:
    return int("1" + word.translate(str.maketrans("LR", "01")), 2)


def complement(word: str) -> str:
    return word.translate(str.maketrans("LR", "RL"))


def reflected_index(k: int) -> int:
    return 3 * (1 << (k.bit_length() - 1)) - 1 - k


def cycle_after_one_flip(vertices: tuple[int, int, int], edge: tuple[int, int]) -> bool:
    x, y, z = sorted(vertices)
    arcs = {(x, y), (x, z), (y, z)}
    lo, hi = sorted(edge)
    arcs.remove((lo, hi))
    arcs.add((hi, lo))
    return (
        {(x, y), (y, z), (z, x)} <= arcs
        or {(y, x), (z, y), (x, z)} <= arcs
    )


def main() -> None:
    monoid_checks = 0
    for p in range(1, 101):
        for q in range(1, 101):
            if gcd(p, q) != 1:
                continue
            require(reverse(left((p, q))) == right(reverse((p, q))), (p, q))
            require(reverse(right((p, q))) == left(reverse((p, q))), (p, q))
            monoid_checks += 2

    word_checks = 0
    for depth in range(16):
        for bits in range(1 << depth):
            word = format(bits, f"0{depth}b").translate(str.maketrans("01", "LR")) if depth else ""
            k = word_index(word)
            require(word_index(complement(word)) == reflected_index(k), word)
            p, q = cw_pair(word)
            require(cw_pair(complement(word)) == (q, p), word)
            word_checks += 1

    seen: dict[tuple[int, int, int], list[Fraction]] = {}
    collision_pairs = 0
    for q in range(2, 241):
        for p in range(1, q):
            if gcd(p, q) != 1:
                continue
            x = Fraction(p, q)
            triple = normalized_unordered_triple(p, q)
            for y in seen.get(triple, []):
                require(x == leg_swap(y), (x, y, triple))
                collision_pairs += 1
            seen.setdefault(triple, []).append(x)
    require(max(map(len, seen.values())) <= 2, "triple fibre")

    # Standard m>n branch triangles.  The closure gcd is controlled by the
    # endpoint parity; the cycle-producing extreme edge is closure for A/B
    # and child for C.
    levels = [[(2, 1)]]
    branch_collisions = []
    primitive_counts = []
    scale2_counts = []
    for depth in range(9):
        primitive = 0
        scale2 = 0
        next_level = []
        for m, n in levels[depth]:
            require(m > n > 0 and gcd(m, n) == 1 and (m - n) % 2 == 1, (m, n))
            parent = Fraction(n, m)
            for branch, a, b in children(m, n):
                require(a > b > 0 and gcd(a, b) == 1, (branch, a, b))
                if branch in ("A", "B"):
                    closure = (n, a)
                    values = (parent, Fraction(b, a), Fraction(n, a))
                    require(gcd(*closure) == gcd(n, 2), (m, n, branch))
                    extreme_edge = closure
                else:
                    closure = (m, a)
                    values = (parent, Fraction(b, a), Fraction(m, a))
                    require(gcd(*closure) == gcd(m, 2), (m, n, branch))
                    extreme_edge = (n, a)
                vertices = (n, m, a)
                require(cycle_after_one_flip(vertices, extreme_edge), (m, n, branch, extreme_edge))
                all_edges = {
                    tuple(sorted((n, m))),
                    tuple(sorted((m, a))),
                    tuple(sorted((n, a))),
                }
                for short_edge in all_edges - {tuple(sorted(extreme_edge))}:
                    require(not cycle_after_one_flip(vertices, short_edge), (m, n, branch, short_edge))
                closure_gcd = gcd(*closure)
                primitive += closure_gcd == 1
                scale2 += closure_gcd == 2
                for i, j in combinations(range(3), 2):
                    if values[j] == leg_swap(values[i]):
                        branch_collisions.append((depth, m, n, branch, i, j, values))
                next_level.append((a, b))
        require(primitive == (3 ** (depth + 1) + (-1) ** depth) // 2, (depth, primitive))
        require(scale2 == (3 ** (depth + 1) - (-1) ** depth) // 2, (depth, scale2))
        primitive_counts.append(primitive)
        scale2_counts.append(scale2)
        levels.append(next_level)

    require(
        branch_collisions
        == [(0, 2, 1, "A", 0, 2, (Fraction(1, 2), Fraction(2, 3), Fraction(1, 3)))],
        branch_collisions,
    )

    # Every orientation on [5] is recovered by its (primitive pair,scale) bits.
    scale_selector_checks = 0
    edges = [(a, b) for a in range(1, 6) for b in range(a + 1, 6)]
    for mask in range(1 << len(edges)):
        selector = {}
        for i, (a, b) in enumerate(edges):
            g = gcd(a, b)
            selector[a // g, b // g, g] = (mask >> i) & 1
        rebuilt = sum(selector[a // gcd(a, b), b // gcd(a, b), gcd(a, b)] << i for i, (a, b) in enumerate(edges))
        require(rebuilt == mask, mask)
        scale_selector_checks += 1

    print("status=INDEPENDENT FINITE-EXACT Calkin--Wilf/Berggren audit")
    print(f"pair_monoid_checks={monoid_checks};cw_word_checks={word_checks}")
    print(f"pythagorean_collision_pairs={collision_pairs};unique_branch_collision=1")
    print(f"primitive_K3_counts={primitive_counts}")
    print(f"scale2_closure_counts={scale2_counts}")
    print(f"scale_selector_tournaments_N5={scale_selector_checks}")
    print("PASS")


if __name__ == "__main__":
    main()
