"""Exact ordered-defect / automorphism checks. No external dependencies.

Run from any directory with Python 3.10+; -O must give identical output.
The defect statistic is relative to the displayed total order, NOT invariant
under relabelling. Its cycle-block formula is not a weighted Burnside census.
"""
from collections import Counter
from itertools import permutations
from math import gcd, lcm


def need(value, message):
    if not value:
        raise RuntimeError(message)


def multiply(p, q):
    r = Counter()
    for i, a in p.items():
        for j, b in q.items():
            r[i + j] += a * b
    return +r


def odd_partitions(n, least=1):
    if n == 0:
        yield ()
    for v in range(least, n + 1, 2):
        for tail in odd_partitions(n - v, v):
            yield (v,) + tail


def cycle_permutation(parts):
    g = []
    for size in parts:
        start = len(g)
        g.extend(list(range(start + 1, start + size)) + [start])
    return tuple(g)


def cycles(g):
    seen, result = set(), []
    for v in range(len(g)):
        if v not in seen:
            c, w = [], v
            while w not in seen:
                seen.add(w)
                c.append(w)
                w = g[w]
            result.append(tuple(c))
    return result


def block_polynomial(parts):
    p = Counter({0: 1})
    for a in parts:
        if a % 2 == 0:
            return Counter()
        for d in range(1, (a + 1) // 2):
            p = multiply(p, Counter({d: 1, a - d: 1}))
    for i, a in enumerate(parts):
        for b in parts[i + 1:]:
            for _ in range(gcd(a, b)):
                p = multiply(p, Counter({0: 1, lcm(a, b): 1}))
    return p


def orbit_polynomial(g):
    """Independent signed pair-orbit walk, arbitrary vertex permutation."""
    remaining = {(i, j) for i in range(len(g)) for j in range(i + 1, len(g))}
    p = Counter({0: 1})
    while remaining:
        first = min(remaining)
        directed, a, b = [], *first
        while tuple(sorted((a, b))) not in [tuple(sorted(e)) for e in directed]:
            directed.append((a, b))
            remaining.remove(tuple(sorted((a, b))))
            a, b = g[a], g[b]
        if (a, b) == first[::-1]:
            return Counter()
        need((a, b) == first, "pair orbit did not close at its start")
        k = sum(a > b for a, b in directed)
        choices = Counter()
        choices[k] += 1
        choices[len(directed) - k] += 1
        p = multiply(p, choices)
    return p


def brute_polynomial(g):
    """Enumerate all labelled tournaments and check the actual arc predicate."""
    n = len(g)
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    index = {edge: k for k, edge in enumerate(edges)}
    image = [(index[tuple(sorted((g[i], g[j])))], g[i] > g[j])
             for i, j in edges]
    p = Counter()
    for mask in range(1 << len(edges)):
        if all(((mask >> k) & 1) == (((mask >> z) & 1) ^ flip)
               for k, (z, flip) in enumerate(image)):
            p[mask.bit_count()] += 1
    return p


def main():
    brute_cases = 0
    for n in range(7):
        for part in odd_partitions(n):
            g = cycle_permutation(part)
            expected = block_polynomial(part)
            need(expected == brute_polynomial(g), f"brute polynomial {part}")
            need(expected == orbit_polynomial(g), f"orbit polynomial {part}")
            brute_cases += 1
    print("Labelled-tournament brute check: all odd cycle types through n=6:",
          brute_cases)
    arbitrary = 0
    hostile = None
    for n in range(1, 8):
        for g in permutations(range(n)):
            cs = cycles(g)
            part = tuple(sorted(map(len, cs)))
            p = orbit_polynomial(g)
            if any(a % 2 == 0 for a in part):
                need(not p, "even cycle fixes a tournament")
            else:
                lower = sum((a * a - 1) // 8 for a in part)
                need(min(p) >= lower, "regular-cycle defect bound")
                moved = sum(a for a in part if a > 1)
                need(moved <= 3 * min(p), "moved-support bound")
                if hostile is None and p != block_polynomial(part):
                    hostile = (n, g, part, dict(sorted(p.items())),
                               dict(sorted(block_polynomial(part).items())))
            arbitrary += 1
    print("Arbitrary permutations through n=7:", arbitrary,
          "; all sparse support / even-cycle controls pass")
    print("First conjugacy-invariance hostile:", hostile)
    checked = 0
    for n in range(31):
        for part in odd_partitions(n):
            p = block_polynomial(part)
            degree = n * (n - 1) // 2
            e = sum((a - 1) // 2 for a in part)
            e += sum(gcd(a, b) for i, a in enumerate(part) for b in part[i + 1:])
            need(sum(p.values()) == 2 ** e, "Davis specialization")
            need(min(p) == sum((a * a - 1) // 8 for a in part), "sharp bound")
            need(p == Counter({degree - k: v for k, v in p.items()}),
                 "converse reciprocity")
            checked += 1
    print("All odd cycle types through n=30:", checked,
          "; Davis specialization, sharp minimum, reciprocity pass")
    for part in [(3,), (3, 3), (5,), (3, 5)]:
        p = block_polynomial(part)
        print("Block type", part, "minimum", min(p), "count", sum(p.values()),
              "polynomial", dict(sorted(p.items())))
    # The sharp 3k support control is an ordered join of k directed triangles.
    for k in range(1, 31):
        p = block_polynomial((3,) * k)
        need(min(p) == k, "sharp triangle family")
    print("Ordered joins of k triangles: support=3k, minimum=k, k=1..30")
    print("PASS")


if __name__ == "__main__":
    main()
