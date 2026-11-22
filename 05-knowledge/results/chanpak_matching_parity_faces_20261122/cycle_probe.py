"""Exact marked-cycle tests for complete-graph and Latin matching rows."""
from collections import Counter
from fractions import Fraction
from itertools import combinations, permutations, product
from math import factorial, gcd

from probe import factorization, row_bits, subsets


def cyclic_orders(n, bits):
    for tail in permutations(range(1, n)):
        cycle = (0,) + tail + (0,)
        forced = backward = 0
        for a, b in zip(cycle, cycle[1:]):
            edge = bits[tuple(sorted((a, b)))]
            forced |= edge
            if a > b:
                backward |= edge
        yield forced, backward


def row_union_sign(bits, union):
    """Return zero for a disconnected union, traversal sign otherwise."""
    neighbors = {v: [] for e in bits for v in e}
    for (a, b), mask in bits.items():
        if union & mask:
            neighbors[a].append(b)
            neighbors[b].append(a)
    seen, edge_sign = {0}, 1
    a, b = 0, neighbors[0][0]
    while True:
        if a > b:
            edge_sign *= -1
        if b == 0:
            break
        seen.add(b)
        u, v = neighbors[b]
        a, b = b, (v if u == a else u)
    return edge_sign if len(seen) == len(neighbors) else 0


def complete_graph_check(n, bits):
    rows = factorization(n, bits)
    hist = Counter()
    for forced, backward in cyclic_orders(n, bits):
        free = (2 ** len(bits) - 1) ^ forced
        for other in subsets(free):
            hist[row_bits(backward | other, rows)] += 1
    total = sum(hist.values())
    signs = {(i, j): row_union_sign(bits, rows[i] | rows[j])
             for i, j in combinations(range(n - 1), 2)}
    sectors = list(product((0, 1), repeat=len(rows)))
    for u in sectors:
        correction = sum(e * (-1) ** (u[i] + u[j]) for (i, j), e in signs.items())
        predicted = (1 + Fraction(2 * correction, factorial(n - 1))) / (2 ** len(rows))
        assert Fraction(hist[u], total) == predicted
    print(f"K{n}: {len(sectors)} sectors, {total} marked pairs, "
          f"{sum(hist[u] == 0 for u in sectors)} empty sectors")


def latin_check(m):
    n = 2 * m
    bits = {e: 1 << i for i, e in enumerate(combinations(range(n), 2))}
    rows = [sum(bits[c, m + (c + i) % m] for c in range(m)) for i in range(m)]
    moments = Counter()
    for forced, backward in cyclic_orders(n, bits):
        for i, j in combinations(range(m), 2):
            union = rows[i] | rows[j]
            if union & forced == union:
                moments[i, j] += (-1) ** ((backward & union).bit_count())
    for i, j in combinations(range(m), 2):
        assert moments[i, j] == (2 * (-1) ** m if gcd(j - i, m) == 1 else 0)
    print(f"cyclic Latin order {m}: {sum(v != 0 for v in moments.values())} detected "
          f"row pairs, signed correlation {Fraction(2 * (-1) ** m, factorial(n - 1))}")


if __name__ == "__main__":
    for n in (4, 6):
        bits = {e: 1 << i for i, e in enumerate(combinations(range(n), 2))}
        complete_graph_check(n, bits)
    for m in (2, 3, 4):
        latin_check(m)
