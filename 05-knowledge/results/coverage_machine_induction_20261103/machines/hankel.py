#!/usr/bin/env python3
"""A small exact lower bound on native-code linear observer dimension."""
import json
from collections import Counter
from functools import lru_cache
from fractions import Fraction
from itertools import product
from check import build, complement, need, path_count, strong


def odd_cycle_count(rows):
    """Independent cycle enumeration and vertex-disjoint odd-cycle expansion."""
    n = len(rows)
    cycles = Counter()
    for start in range(n):
        def visit(last, mask, length):
            if length >= 3 and length % 2 and rows[last] >> start & 1:
                cycles[mask] += 1
            candidates = rows[last] & ~mask & ~((1 << (start+1))-1)
            while candidates:
                bit = candidates & -candidates
                visit(bit.bit_length()-1, mask | bit, length+1)
                candidates ^= bit
        visit(start, 1 << start, 1)
    at = [[] for _ in range(n)]
    for mask, multiplicity in cycles.items():
        for v in range(n):
            if mask >> v & 1:
                at[v].append((mask, multiplicity))
    @lru_cache(None)
    def pack(mask):
        if not mask:
            return 1
        bit = mask & -mask
        v = bit.bit_length()-1
        return pack(mask ^ bit) + sum(2*c*pack(mask ^ cyc)
                    for cyc, c in at[v] if cyc & mask == cyc)
    return pack((1 << n)-1), sum(cycles.values())


def determinant(matrix):
    a = [list(map(Fraction, row)) for row in matrix]
    answer = Fraction(1)
    for j in range(len(a)):
        pivot = next((i for i in range(j, len(a)) if a[i][j]), None)
        if pivot is None:
            return 0
        if pivot != j:
            a[pivot], a[j] = a[j], a[pivot]
            answer *= -1
        p = a[j][j]
        answer *= p
        for k in range(j, len(a)):
            a[j][k] /= p
        for i in range(j+1, len(a)):
            c = a[i][j]
            for k in range(j, len(a)):
                a[i][k] -= c*a[j][k]
    need(answer.denominator == 1, "integral determinant")
    return answer.numerator


def main():
    words = [()] + [(0,), (1,)] + list(product((0, 1), repeat=2))
    cache = {}
    records = []
    for length in range(5):
        for suffix in product((0, 1), repeat=length):
            code = (1,) + suffix
            source = build(code)[0]
            native = complement(source)
            need(strong(source) and strong(native), "strong common-source domain")
            h = path_count(native)
            h2, cycles = odd_cycle_count(native)
            need(h == h2, "independent odd-cycle count")
            cache[suffix] = h
            records.append(dict(code=list(code), order=len(native), rows=native,
                                H=h, directed_odd_cycles=cycles))
    matrix = [[cache[u+v] for v in words] for u in words]
    d = determinant(matrix)
    need(d == -66127356928, "Hankel determinant")
    print(json.dumps(dict(
        status="FINITE-EXACT certificate; implies linear dimension at least7",
        scope="31 anchored strong native sources, orders5..13; all suffixes length0..4",
        prefix_and_suffix_words=[list(w) for w in words],
        matrix=matrix, determinant=d, records=records), sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
