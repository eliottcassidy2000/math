"""Exploratory exact rank probe for one cubic vertex plus paired ports."""
from itertools import combinations, product
from collections import Counter
import sympy as sp


def matchings(items):
    if not items:
        yield ()
        return
    a = items[0]
    for i in range(1, len(items)):
        for rest in matchings(items[1:i] + items[i + 1:]):
            yield ((a, items[i]),) + rest


def canonical_word(word):
    mapping = {}
    return tuple(mapping.setdefault(c, len(mapping)) for c in word)


def rank_mod(rows, p=101):
    a = [list(row) for row in rows]
    r = 0
    for col in range(len(a[0])):
        pivot = next((i for i in range(r, len(a)) if a[i][col] % p), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        inv = pow(a[r][col] % p, -1, p)
        a[r] = [(x * inv) % p for x in a[r]]
        for i in range(r + 1, len(a)):
            scale = a[i][col] % p
            if scale:
                a[i] = [(x - scale * y) % p for x, y in zip(a[i], a[r])]
        r += 1
    return r


for n in (3, 5, 7, 9):
    words = sorted({canonical_word(word) for word in product(range(3), repeat=n)
                    if all(word.count(c) % 2 for c in range(3))})
    kernels = []
    for triple in combinations(range(n), 3):
        rest = tuple(v for v in range(n) if v not in triple)
        for pairs in matchings(rest):
            kernels.append([int(len({word[v] for v in triple}) == 3
                                and all(word[a] == word[b] for a, b in pairs))
                            for word in words])
    rank = rank_mod(kernels)
    print("n", n, "kernels", len(kernels), "parity_orbits", len(words),
          "rank_mod101", rank, flush=True)
    if n <= 7:
        print("rank_Q", sp.Matrix(kernels).rank(),
              "column_sums", dict(Counter(map(sum, zip(*kernels)))), flush=True)
