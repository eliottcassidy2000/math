"""Exact divisor-layer / simplex-pair audit.  Standard library only.

Run: python -X utf8 -B 04-computation/experiments/divisor_layers_pair_geometry_20261004.py
All checks are explicit and survive python -O.  No imported local experiment.
"""

from collections import Counter
from itertools import combinations_with_replacement, permutations, product
from math import prod


def need(condition, message):
    if not condition:
        raise ArithmeticError(message)


def factor(n):
    out = []
    p = 2
    while p * p <= n:
        a = 0
        while n % p == 0:
            n //= p
            a += 1
        if a:
            out.append((p, a))
        p += 1
    if n > 1:
        out.append((n, 1))
    return out


def divisor_stats(n):
    ds = [d for d in range(2, n) if n % d == 0]
    return len(ds), sum(all(a == 1 for _, a in factor(d)) for d in ds), sum(
        factor(d) == [(d, 1)] for d in ds
    )


def exponent_points(k, m=2):
    return tuple(
        v for v in product(range(k + 1), *([range(2)] * m))
        if v != (0,) * (m + 1) and v != (k,) + (1,) * m
    )


def character(s, bits):
    return (-1) ** sum(x * y for x, y in zip(s, bits))


def convolution(left, right):
    result = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            result[i + j] += a * b
    return result


def walsh_rank_direct(k, s):
    result = [0] * (k + len(s) + 1)
    for v in exponent_points(k, len(s)):
        result[sum(v)] += character(s, v[1:])
    return tuple(result)


def walsh_rank_product(k, s):
    result = [1] * (k + 1)
    for bit in s:
        result = convolution(result, [1, (-1) ** bit])
    result[0] -= 1
    result[-1] -= (-1) ** sum(s)
    return tuple(result)


def leq(v, w):
    return all(x <= y for x, y in zip(v, w))


def rank_chart(v):
    return (sum(v), *v[1:])


def inverse_rank_chart(v):
    return (v[0] - sum(v[1:]), *v[1:])


def paired_alphabet_map(v):
    """Inherited labelled p^2qr -> Sym^2({0,p=1,q=2,r=3}) map."""
    table = {
        (1, 0, 0): (0, 1), (0, 1, 0): (0, 2), (0, 0, 1): (0, 3),
        (1, 1, 0): (1, 2), (1, 0, 1): (1, 3), (0, 1, 1): (2, 3),
        (1, 1, 1): (0, 0), (2, 0, 0): (1, 1),
        (2, 1, 0): (2, 2), (2, 0, 1): (3, 3),
    }
    return table[v]


def pair_histogram(pair, size=4):
    return tuple(pair.count(i) for i in range(size))


def involutions(n):
    for pi in permutations(range(n)):
        if all(pi[pi[i]] == i for i in range(n)):
            yield pi


def poset_automorphisms(points):
    """Independent finite order-matrix backtracking, not a coordinate ansatz."""
    rel = tuple(tuple(leq(v, w) for w in points) for v in points)
    signatures = [(sum(rel[j][i] for j in range(len(points))), sum(rel[i]))
                  for i in range(len(points))]
    candidates = [[j for j, t in enumerate(signatures) if t == s]
                  for s in signatures]
    order = sorted(range(len(points)), key=lambda i: len(candidates[i]))
    solutions = []

    def visit(position, assignment, used):
        if position == len(order):
            solutions.append(tuple(assignment[i] for i in range(len(points))))
            return
        i = order[position]
        for j in candidates[i]:
            if j in used:
                continue
            if all(rel[i][a] == rel[j][b] and rel[a][i] == rel[b][j]
                   for a, b in assignment.items()):
                assignment[i] = j
                visit(position + 1, assignment, used | {j})
                del assignment[i]

    visit(0, {}, set())
    return solutions


def main():
    print("DIVISOR LAYERS / PAIR GEOMETRY: EXACT AUDIT")
    print("Universe: N=2..500; k=2..20 for p^kqr; k=1..8,m=1..5 for general Walsh;")
    print("all simplex vertex involutions v=1..8; independent proper-poset aut k=2..6;")
    print("all words length 1..5 for ten-letter charges, using exact convolution.")
    print()

    balanced = []
    for n in range(2, 501):
        fsu = divisor_stats(n)
        aa = [a for _, a in factor(n)]
        r = len(aa)
        predicted = (prod(a + 1 for a in aa) - 2,
                     2 ** r - 1 - int(all(a == 1 for a in aa)),
                     r - int(aa == [1]))
        need(fsu == predicted, ("stats", n, fsu, predicted))
        classified = aa == [1] or aa == [3] or sorted(aa) == [1, 1, 2]
        need((fsu[0] == fsu[1] + fsu[2]) == classified, ("classification", n))
        if classified:
            balanced.append(n)
    print(f"Inherited F=S+U direct audit: 499 integers, {len(balanced)} balanced.")
    for n in (2, 8, 60, 63, 126):
        print(f"  N={n}: (F,S,U)={divisor_stats(n)}")

    general_count = 0
    for k in range(1, 9):
        for m in range(1, 6):
            for s in product(range(2), repeat=m):
                direct = walsh_rank_direct(k, s)
                need(direct == walsh_rank_product(k, s), ("product", k, m, s))
                sign = (-1) ** sum(s)
                need(tuple(reversed(direct)) == tuple(sign * x for x in direct),
                     ("complement", k, m, s))
                expected = (k + 1) * 2 ** m - 2 if not any(s) else -1 - sign
                need(sum(direct) == expected, ("integrated spectrum", k, m, s))
                general_count += 1
    print(f"General rank-Walsh product/duality/integrated spectrum: {general_count} cases.")

    order_checks = 0
    for k in range(2, 21):
        points = exponent_points(k)
        need(len(points) == 4 * k + 2, ("count", k))
        ranks = Counter(map(sum, points))
        need([ranks[j] for j in range(1, k + 2)] == [3] + [4] * (k - 1) + [3],
             ("rank sizes", k))
        fibres = Counter((v[1], v[2]) for v in points)
        need([fibres[x] for x in product(range(2), repeat=2)] == [k, k + 1, k + 1, k],
             ("charge fibres", k))
        for v in points:
            image = rank_chart(v)
            need(inverse_rank_chart(image) == v, ("rank inverse", k, v))
            comp = (k - v[0], 1 - v[1], 1 - v[2])
            need(rank_chart(comp) == (k + 2 - image[0], 1 - image[1], 1 - image[2]),
                 ("rank complement", k, v))
            for w in points:
                iw = rank_chart(w)
                expected = image[1] <= iw[1] and image[2] <= iw[2] and (
                    image[0] - image[1] - image[2] <= iw[0] - iw[1] - iw[2])
                need(leq(v, w) == expected, ("rank order", k, v, w))
                order_checks += 1
        for s in ((0, 1), (1, 0), (1, 1)):
            expected = [0] * (k + 3)
            expected[1] = -1 if s == (1, 1) else 1
            expected[k + 1] = -1
            need(walsh_rank_direct(k, s) == tuple(expected), ("rank simplification", k, s))
    print(f"p^kqr rank chart, inverse and order: k=2..20, {order_checks} pair comparisons.")
    for k in range(2, 7):
        points = exponent_points(k)
        found = poset_automorphisms(points)
        expected = {tuple(range(len(points))),
                    tuple(points.index((a, c, b)) for a, b, c in points)}
        need(set(found) == expected, ("automorphisms", k))
    print("Proper-poset automorphisms: exactly identity and q/r exchange for k=2..6.")

    pairs = tuple(combinations_with_replacement(range(4), 2))
    points = exponent_points(2)
    need(set(map(paired_alphabet_map, points)) == set(pairs), "pair bijection")
    need(len(set(map(pair_histogram, pairs))) == 10, "histogram bijection")
    height_layers = Counter(sum(i >= 2 for i in pair) for pair in pairs)
    need([height_layers[j] for j in range(3)] == [3, 4, 3], "tetrahedral 2+2 layers")
    layer_map = {}
    for j in range(1, 4):
        source = sorted(v for v in points if sum(v) == j)
        target = sorted(pair for pair in pairs if sum(i >= 2 for i in pair) == j - 1)
        layer_map.update(zip(source, target))
    need(len(set(layer_map.values())) == 10, "graded set bijection")
    need(any(layer_map[v] != paired_alphabet_map(v) for v in points), "different rank map")
    pair_charge = Counter(a ^ b for a, b in pairs)
    natural_charge = Counter(2 * b + c for a, b, c in points)
    need([pair_charge[i] for i in range(4)] == [4, 2, 2, 2], "pair charge")
    need([natural_charge[i] for i in range(4)] == [2, 3, 3, 2], "natural charge")
    need(sorted(pair_charge.values()) != sorted(natural_charge.values()), "no charge relabel")
    print("Ten-point bridge: inherited pair map is bijective; histograms are 4 corners + 6 midpoints.")
    print("  A 2+2 vertex partition gives parallel height layers=(3,4,3); rankwise sorted matching")
    print("  is a graded-set bijection, different from the inherited squarefree/pair map.")
    print("  Pair-XOR multiplicities=(4,2,2,2), spectrum=(10,2,2,2).")
    print("  q/r exponent multiplicities=(2,3,3,2), spectrum=(10,0,0,-2).")
    print("  k=2 rank-Walsh rows:")
    for s in product(range(2), repeat=2):
        print(f"    {s}: {walsh_rank_direct(2, s)}")

    def pair_xor(v):
        x, y = paired_alphabet_map(v)
        return x ^ y

    u, v, uv = (2, 0, 0), (0, 1, 0), (2, 1, 0)
    need((pair_xor(u) ^ pair_xor(v)) != pair_xor(uv), "multiplication hostile")
    print("  Arithmetic hostile: pair-XOR(4) XOR pair-XOR(3)=0 XOR 2=2, but pair-XOR(12)=0 for N=60.")

    counts = Counter({0: 1})
    for length in range(1, 6):
        next_counts = Counter()
        for old, amount in counts.items():
            for charge in (2 * b + c for a, b, c in points):
                next_counts[old ^ charge] += amount
        counts = next_counts
        expected = [(10 ** length + (-2) ** length * (-1) ** i.bit_count()) // 4
                    for i in range(4)]
        need([counts[i] for i in range(4)] == expected, ("word convolution", length))
        print(f"  Natural word charges length {length}: {expected}")

    inv_count = 0
    for n in range(1, 9):
        letters = tuple(combinations_with_replacement(range(n), 2))
        histogram = Counter()
        for pi in involutions(n):
            fixed_vertices = sum(pi[i] == i for i in range(n))
            fixed_pairs = sum(tuple(sorted((pi[a], pi[b]))) == (a, b) for a, b in letters)
            need(2 * fixed_pairs == n + fixed_vertices ** 2, ("fixed count", n, pi))
            need(fixed_pairs >= (n + 1) // 2, ("lower bound", n, pi))
            histogram[fixed_pairs] += 1
            inv_count += 1
        print(f"Simplex v={n}: involution fixed-letter distribution={dict(sorted(histogram.items()))}")
    print(f"All involutions through 8 vertices: {inv_count} exact checks.")
    need(Counter((2 - a, 1 - b, 1 - c) for a, b, c in points) == Counter(points),
         "divisor complement permutation")
    need(not any((a, b, c) == (2 - a, 1 - b, 1 - c) for a, b, c in points),
         "complement fixed-free")
    print("Complement obstruction: divisor ten-set has 0 fixed; tetrahedral involutions have 2,4,10.")
    print("Sharp boundary: p^4 has 3 proper divisors; exponents 1,2,3 map to 0,1/2,1,")
    print("and complement becomes segment reflection (one fixed midpoint).")
    print("PASS: inherited classifications, new all-height formulas, exact maps and hostile controls.")


if __name__ == "__main__":
    main()
