"""Exact small-order tests for mixed flips and prescribed matching parities."""
from collections import Counter
from fractions import Fraction
from itertools import combinations, permutations, product
from math import factorial


def subsets(mask):
    child = mask
    while child:
        yield child
        child = (child - 1) & mask
    yield 0


def paths(n, bit):
    for order in permutations(range(n)):
        forced, backward = 0, 0
        for a, b in zip(order, order[1:]):
            edge = bit[tuple(sorted((a, b)))]
            forced |= edge
            if a > b:
                backward |= edge
        yield forced, backward


def counts(n, bit):
    full = (1 << len(bit)) - 1
    h = [0] * (full + 1)
    for forced, backward in paths(n, bit):
        for other in subsets(full ^ forced):
            h[backward | other] += 1
    return h


def nonzero_mixed_faces(h):
    w = h.copy()
    step = 1
    while step < len(w):
        for start in range(0, len(w), 2 * step):
            for j in range(start, start + step):
                a, b = w[j], w[j + step]
                w[j], w[j + step] = a + b, a - b
        step *= 2
    active = [x != 0 for x in w]
    edge = 1
    while edge < len(active):
        for mask in range(len(active)):
            if not mask & edge:
                active[mask] |= active[mask | edge]
        edge *= 2
    return active


def admissible(n, bit, mask):
    neighbors = [[] for _ in range(n)]
    for (a, b), e in bit.items():
        if e & mask:
            neighbors[a].append(b)
            neighbors[b].append(a)
    if any(len(adj) > 2 for adj in neighbors):
        return False
    unseen = set(range(n))
    odd_component = False
    while unseen:
        todo = [unseen.pop()]
        vertices, incidents = 0, 0
        while todo:
            a = todo.pop()
            vertices += 1
            incidents += len(neighbors[a])
            for b in neighbors[a]:
                if b in unseen:
                    unseen.remove(b)
                    todo.append(b)
        if incidents != 2 * (vertices - 1):
            return False
        odd_component |= bool(vertices % 2)
    return odd_component


def factorization(n, bit):
    # Circle construction on the moving vertices modulo n-1.
    rows = []
    for k in range(n - 1):
        pairs = [(k, n - 1)] + [
            tuple(sorted(((k + i) % (n - 1), (k - i) % (n - 1))))
            for i in range(1, n // 2)
        ]
        rows.append(sum(bit[e] for e in pairs))
    return rows


def row_bits(mask, rows):
    return tuple((mask & row).bit_count() % 2 for row in rows)


def verify_sampler(n, bit, rows, h):
    # Exact history probabilities, scaled by n! 2^(d-n+1-k).
    # At order four this scale can be a fraction.
    edges = {e: pair for pair, e in bit.items()}
    full = (1 << len(bit)) - 1
    all_paths = list(paths(n, bit))
    for target in product((0, 1), repeat=len(rows)):
        histogram = Counter()
        for forced, backward in all_paths:
            saturated = [i for i, row in enumerate(rows) if row & forced == row]
            assert len(saturated) <= 1
            if saturated:
                i = saturated[0]
                if row_bits(backward, [rows[i]])[0] != target[i]:
                    a, b = edges[rows[i] & -rows[i]]
                    swap = {a: b, b: a}
                    new_forced, new_backward = 0, 0
                    for e, (u, v) in edges.items():
                        if forced & e:
                            x, y = swap.get(u, u), swap.get(v, v)
                            new_e = bit[tuple(sorted((x, y)))]
                            new_forced |= new_e
                            if bool(backward & e) != (x > y):
                                new_backward |= new_e
                    forced, backward = new_forced, new_backward
            for other in subsets(full ^ forced):
                mask = backward | other
                if row_bits(mask, rows) == target:
                    histogram[mask] += Fraction(1, 2 ** len(saturated))
        assert all(histogram[mask] == h[mask] for mask in range(len(h))
                   if row_bits(mask, rows) == target)
    print("n=4 sampler gives exactly H(T) mass in all eight sectors")


def main():
    for n in range(2, 7):
        bit = {pair: 1 << i for i, pair in enumerate(combinations(range(n), 2))}
        h = counts(n, bit)
        active = nonzero_mixed_faces(h)
        assert all(x == admissible(n, bit, mask) for mask, x in enumerate(active))
        print(f"n={n}: {sum(active)} of {len(h)} flip sets detected")
        if n % 2 == 0:
            rows = factorization(n, bit)
            hist = Counter()
            for mask, count in enumerate(h):
                hist[row_bits(mask, rows)] += count
            assert set(hist.values()) == {sum(h) // (2 ** len(rows))}
            print(f"  {len(hist)} row sectors, marked mass {min(hist.values())} each")
            if n == 4:
                for sector in [(0, 0, 0), (0, 1, 0)]:
                    data = Counter(x for mask, x in enumerate(h) if row_bits(mask, rows) == sector)
                    print(f"  sector {sector}: count multiplicities {sorted(data.items())}")
                verify_sampler(n, bit, rows, h)
        else:
            m = (n - 1) // 2
            rows = [sum(bit[(i, i + 1)] for i in range(j, n - 1, 2)) for j in (0, 1)]
            M = rows[0]
            def even_sum(mask):
                return (-1) ** (mask.bit_count() % 2)
            squares = [sum(even_sum(b) * h[x | b] for b in subsets(M)) ** 2
                       for x in subsets((len(h) - 1) ^ M)]
            mean_square = Fraction(sum(squares), len(squares))
            correlation = Fraction(sum(v * even_sum(x & (rows[0] | rows[1]))
                                        for x, v in enumerate(h)), sum(h))
            assert mean_square == Fraction(factorial(m + 1), 2 ** (m - 1))
            assert correlation == Fraction(2, factorial(n))
            print(f"  mixed-square mean {mean_square}, row correlation {correlation}")


if __name__ == "__main__":
    main()
