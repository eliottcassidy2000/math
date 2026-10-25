#!/usr/bin/env python3
"""Independent bounded monoid-generation audit, without importing check.py."""
from collections import deque
from itertools import combinations, permutations, product
from pathlib import Path
import hashlib
import json


def code(labels):
    return frozenset(tuple((a & x).bit_count() % 2 for x in labels[1:])
                     for a in range(8))


def main():
    codes = {code((1,)+p) for p in permutations(range(2, 8))}
    zero = (0,)*6
    gates = 0
    pairs = 0
    index_valuations = set()
    for c, d in combinations(codes, 2):
        if c & d != {zero}:
            continue
        pairs += 1
        generators = c | d | {tuple(2*int(i == j) for i in range(6)) for j in range(6)}
        reached = {zero}
        queue = deque([zero])
        while queue:
            v = queue.popleft()
            for g in generators:
                w = tuple(a+b for a, b in zip(v, g))
                if max(w) <= 2 and w not in reached:
                    reached.add(w)
                    queue.append(w)
        # Decode only after generation, then compare all 729 possible monomials.
        minima = {tuple((a+b) % 2 for a, b in zip(p, q)):
                  tuple(a+b for a, b in zip(p, q)) for p in c for q in d}
        for v in product(range(3), repeat=6):
            expected = all(a >= b for a, b in zip(v, minima[tuple(a % 2 for a in v)]))
            if (v in reached) != expected:
                raise RuntimeError('monoid-generation mismatch')
            gates += 1
        indices = tuple(sum((v[i] - (v[i] % 2))//2 for v in minima.values()) for i in range(6))
        index_valuations.add(indices)
        if indices != (16,)*6:
            raise RuntimeError('inclusion determinant valuation')
        # All six root factors are needed: each absent factor admits a hole
        # after multiplication by a normalized monomial in the unit cube.
        for i in range(6):
            a = tuple(int(j != i) for j in range(6))
            if not any(tuple(x+y for x,y in zip(a,w)) not in reached
                       for w in product(range(2), repeat=6)):
                raise RuntimeError('missing factor did not expose a monoid hole')
    print(json.dumps({'status': 'FINITE-EXACT', 'pairs': pairs, 'membership_gates': gates,
                     'box': '[0,2]^6', 'inclusion_determinant_branch_valuations': list(next(iter(index_valuations))),
                     'scope': 'Independent additive-monoid BFS for all 120 pairs in chart 0; no algebraic theorem inferred from a bounded box.',
                     'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
