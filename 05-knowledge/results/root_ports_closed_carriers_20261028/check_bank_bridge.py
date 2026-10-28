#!/usr/bin/env python3
"""Explicit six-field-coordinate -> K4 readout -> nine-prime K3,3 map."""
from itertools import combinations
from math import isqrt
import json


def need(ok, why):
    if not ok:
        raise RuntimeError(why)


def xor(xs):
    z = 0
    for x in xs: z ^= x
    return z


def product(xs):
    z = 1
    for x in xs: z *= x
    return z


EDGES = tuple(combinations(range(4), 2))
BASIS = (58, 20, 1, 26, 6, 11)  # F64 polynomial masks assigned to K4 edges
PROBES = (48, 29, 25, 56)
PORTS = (1, 2, 4, 5, 8, 16, 24, 32, 34)
L = (1, 2, 12, 7)  # columns of the invertible four-bit output change
PRIME = {1: 7, 2: 2, 4: 19, 5: 89, 9: 13, 6: 59, 10: 83, 12: 139, 11: 41}


def linear(x, cols):
    return xor(c for i, c in enumerate(cols) if x >> i & 1)


def observation(x):
    return sum((sum((x >> i) & 1 for i, e in enumerate(EDGES) if v in e) % 2) << v
               for v in range(3)) | ((x.bit_count() % 2) << 3)


def mul(x, y):
    z = 0
    while y:
        if y & 1: z ^= x
        y >>= 1; x <<= 1
        if x & 64: x ^= 67  # t^6+t+1
    return z


def power(x, k):
    z = 1
    while k:
        if k & 1: z = mul(z, x)
        x = mul(x, x); k >>= 1
    return z


def trace(x):
    z = 0
    for _ in range(6):
        z ^= x; x = mul(x, x)
    need(z in (0, 1), 'absolute trace lands in F2')
    return z


def squareclass(x):
    result = 1; p = 2
    while p*p <= x:
        parity = 0
        while x % p == 0:
            x //= p; parity ^= 1
        if parity: result *= p
        p += 1
    return result*x


def main():
    need(len({linear(x, BASIS) for x in range(64)}) == 64, 'six field coordinates form basis')
    need(len({linear(x, L) for x in range(16)}) == 16, 'four-bit output change invertible')
    for x in range(64):
        field = linear(x, BASIS)
        probe = sum(trace(mul(field, power(p, 8))) << i for i, p in enumerate(PROBES))
        need(probe == observation(x), 'literal F64 trace equals K4 readout')
    need(all(power(linear(x, BASIS), 9) == 1 for x in PORTS), 'all nine ports norm one')
    observed = {x: linear(observation(x), L) for x in PORTS}
    need(set(observed.values()) == set(PRIME), 'nine prime square-class bank image')
    old_lines = [t for t in combinations(PORTS, 3) if xor(t) == 0]
    all_lines = [t for t in combinations(PORTS, 3) if xor(observed[x] for x in t) == 0]
    new_lines = [t for t in all_lines if t not in old_lines]
    need(len(old_lines) == len(new_lines) == 3, 'two rulings of three lines')
    need(all(len(set(t) & set(s)) == 1 for t in old_lines for s in new_lines), 'literal K3,3 incidence')
    need({xor(t) for t in new_lines} == {30, 45, 51}, 'new dependencies lift to K4 four-cycles')
    masks = {x: i for i, x in enumerate(PORTS)}
    graph_edges = [(next(i for i, t in enumerate(old_lines) if x in t),
                    3+next(j for j, t in enumerate(new_lines) if x in t)) for x in PORTS]
    cuts = {sum(((chosen >> u & 1) ^ (chosen >> v & 1)) << i for i, (u, v) in enumerate(graph_edges))
            for chosen in range(64)}
    relations = set(); field_relations = set(); histogram = {}
    for selection in range(512):
        xs = [x for i, x in enumerate(PORTS) if selection >> i & 1]
        syndrome = xor(observed[x] for x in xs)
        sigma = product(PRIME[observed[x]]+1 for x in xs)
        square = isqrt(sigma)**2 == sigma
        need(square == (syndrome == 0), 'exact divisor sum square test on all512 subsets')
        if not syndrome:
            relations.add(selection)
            histogram[len(xs)] = histogram.get(len(xs), 0)+1
        if not xor(xs): field_relations.add(selection)
    need(relations == cuts and len(relations) == 32, 'square-total subsets equal all K3,3 cuts')
    need(len(field_relations) == 8 and field_relations < relations, 'two new relation bits, not a field isomorphism')
    rows = []
    for x in PORTS:
        vector = observed[x]; p = PRIME[vector]
        sf = product(q for i, q in enumerate((2, 3, 5, 7)) if vector >> i & 1)
        need(squareclass(p+1) == sf, 'prime p+1 square class')
        rows.append(dict(edge_mask=x, field_mask=linear(x, BASIS), four_probe_mask=observation(x),
                         square_class_mask=vector, square_class=sf, prime=p))
    return dict(status='FINITE-EXACT explicit linear bridge', edge_order=EDGES,
                output_change_columns=L, table=rows, original_field_lines=old_lines,
                added_prime_bank_lines=new_lines, added_kernel_lifts=[xor(t) for t in new_lines],
                prime_bank=sorted(PRIME.values()), square_total_subset_count=32,
                weight_enumerator=histogram, field_relation_count=8,
                lost_relation_dimensions=2, complete_subset_universe=512)


if __name__ == '__main__':
    print(json.dumps(main(), indent=2, sort_keys=True))
