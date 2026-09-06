#!/usr/bin/env python3
"""Independent exact symbolic and small-universe signed-cycle proof controls.

The separate packed-bitset C++ companion replays the entire K8 finite base.
This file imports no repository mathematics and allocates no K9 character array.
"""
from fractions import Fraction as Q
from itertools import combinations, permutations
from math import factorial, comb
import sys

CHECKS = 0


def need(ok, label):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(label)


def trim(p):
    p = list(map(Q, p))
    while len(p) > 1 and not p[-1]:
        p.pop()
    return tuple(p)


def add(p, q):
    return trim(tuple((p[j] if j < len(p) else 0) +
                      (q[j] if j < len(q) else 0)
                      for j in range(max(len(p), len(q)))))


def scale(p, q):
    return trim(tuple(q*x for x in p))


def mul(p, q):
    out = [Q(0)]*(len(p)+len(q)-1)
    for i, x in enumerate(p):
        for j, y in enumerate(q):
            out[i+j] += x*y
    return trim(out)


def shift(p, c):
    return add(p, (Q(c),))


def falling(p, length):
    result = (Q(1),)
    for j in range(length):
        result = mul(result, shift(p, -j))
    return result


def symbolic_controls():
    n = (Q(9), Q(1))
    nm2, nm3, nm4 = (shift(n, -j) for j in (2, 3, 4))
    middle = mul(nm2, nm3)
    threshold = mul(middle, nm4)
    qlo, qhi = scale(nm4, 2), scale(mul(nm2, nm2), Q(1, 4))
    for q, claimed in (
        (qlo, mul(nm4, (Q(2), Q(5), Q(1)))),
        (qhi, scale(mul(mul(nm2, nm4), (Q(1), Q(6), Q(1))), Q(1, 8))),
    ):
        excess = add(mul(q, add(middle, scale(q, -2))), scale(threshold, -1))
        need(excess == claimed, "strict transposition endpoint polynomial identity")
        need(all(x > 0 for x in excess), "strict positive-coefficient tail")
    n = (Q(8), Q(1))
    excess = (Q(0),)
    for k in (3, 5, 7):
        excess = add(excess, scale(falling(n, k), Q(1, 2*k)))
    for k in range(3, 9):
        excess = add(excess, scale(falling(shift(n, -2), k-2), -1))
    claimed = (Q(1652), Q(244667, 105), Q(2723, 2), Q(1328, 3),
               Q(189, 2), Q(73, 5), Q(3, 2), Q(1, 14))
    need(excess == claimed, "complete antibalanced S7 polynomial")
    for k in range(4, 42, 2):
        n = (Q(k), Q(1))
        excess = add(mul(n, shift(n, -1)), scale(shift(n, -k+2), -2*(k-1)))
        need(excess == (Q((k-1)*(k-4)), Q(1), Q(1)), "odd/even cumulative pair identity")
    return tuple(map(str, claimed))


def edge_index(n):
    pairs = tuple(combinations(range(n), 2))
    return pairs, {pair:j for j, pair in enumerate(pairs)}


def bit(index, u, v):
    return 1 << index[tuple(sorted((u, v)))]


def cycle_masks(n, k):
    _, index = edge_index(n)
    out = []
    for vertices in combinations(range(n), k):
        root = vertices[0]
        for tail in permutations(vertices[1:]):
            if tail[0] > tail[-1]:
                continue
            word = (root,) + tail
            out.append(sum(bit(index, word[j], word[(j+1) % k]) for j in range(k)))
    need(len(out) == factorial(n)//(2*k*factorial(n-k)), "unoriented simple cycle universe")
    need(len(out) == len(set(out)), "simple cycle masks injective")
    return out


def negative_count(h, cycles):
    return sum((h & cycle).bit_count() % 2 for cycle in cycles)


def two_star(n, r, index, u=0, v=1):
    chosen = [x for x in range(n) if x not in (u, v)][:r]
    return sum(bit(index, u, x) + bit(index, v, x) for x in chosen)


def star_count(n, r):
    s = n-2-r
    return 2*r*s*((n-2)*(n-3)-2*r*s)*factorial(n-5)


def direct_star_controls():
    rows = []
    for n in range(6, 11):
        _, index = edge_index(n)
        cycles = cycle_masks(n, n)
        counts = []
        for r in range(n-1):
            actual = negative_count(two_star(n, r, index), cycles)
            need(actual == star_count(n, r), ("two-star direct cycle count", n, r))
            counts.append(actual)
            if n >= 9 and 2 <= r <= n-4:
                need(actual > 2*factorial(n-2), ("strict interior threshold", n, r))
        rows.append((n, tuple(counts)))
    need(star_count(8, 3) == 1296 < 2*factorial(6), "K8 analytic-threshold hostile")
    return rows


def crossing_gap_controls():
    cases = 0
    for m in range(4, 13):
        first, second, sizes = [0]*(m+1), [0]*(m+1), [0]*(m+1)
        # Fix the labelled cycle and vary its R-set; this is independent of
        # enumerating Hamilton cycles with a fixed R-set.
        for mask in range(1 << m):
            r = mask.bit_count()
            x = sum(((mask >> j) ^ (mask >> ((j+1) % m))) & 1 for j in range(m))
            first[r] += x
            second[r] += x*x
            sizes[r] += 1
        for r in range(m+1):
            s = m-r
            need(sizes[r] == comb(m, r), "complete fixed-cycle color universe")
            need(Q(first[r], sizes[r]) == Q(2*r*s, m-1), "crossing-gap first moment")
            need(Q(second[r], sizes[r]) == Q(4*r*s*(r*s-1), (m-1)*(m-2)),
                 "crossing-gap second moment")
            inserted = (m-1)*factorial(m-2)*(Q((m+1)*first[r]-second[r], sizes[r]))
            need(inserted == star_count(m+2, r), "deleted-cycle insertion identity")
            cases += 1
    return cases


def graph_rigidity_controls():
    rows = []
    for m in (5, 6):
        pairs = tuple(combinations(range(m), 2))
        end = 1 << len(pairs)
        expected = {0, end-1}
        for j in range(len(pairs)):
            expected.update((1 << j, (end-1) ^ (1 << j)))
        for v in range(m):
            star = sum(1 << j for j, pair in enumerate(pairs) if v in pair)
            expected.update((star, (end-1) ^ star))
        found = set()
        allowed = {0, 1, m-2, m-1}
        for h in range(end):
            rows_of_g = [0]*m
            for j, (u, v) in enumerate(pairs):
                if (h >> j) & 1:
                    rows_of_g[u] |= 1 << v
                    rows_of_g[v] |= 1 << u
            if any(row.bit_count() not in allowed for row in rows_of_g):
                continue
            if all(((rows_of_g[u] ^ rows_of_g[v]) & ~((1 << u) | (1 << v))).bit_count()
                   in allowed for u, v in pairs):
                found.add(h)
        need(found == expected, ("complete six-family disagreement classification", m))
        rows.append((m, end, len(found)))
    return rows


def deletion_and_transposition_controls():
    n = 6
    pairs, index = edge_index(n)
    nongauge = tuple(pair for pair in pairs if 0 not in pair)
    end = 1 << len(nongauge)
    physical = [bit(index, *pair) for pair in nongauge]
    single = {1 << j for j in range(len(nongauge))}
    for v in range(1, n):
        single.add(sum(1 << j for j, pair in enumerate(nongauge) if v in pair))
    negsingle = {(end-1) ^ h for h in single}
    deletions = []
    for v in range(n):
        masks = []
        for tri in combinations([x for x in range(n) if x != v], 3):
            masks.append(sum(bit(index, *pair) for pair in combinations(tri, 2)))
        deletions.append(masks)
    hamilton = cycle_masks(n, n)
    for gauge in range(end):
        h = sum(physical[j] for j in range(len(nongauge)) if (gauge >> j) & 1)
        types = [negative_count(h, masks) for masks in deletions]
        balanced = types.count(0)
        antibalanced = types.count(comb(n-1, 3))
        need(not (balanced and antibalanced), "opposite exceptional deletion types cannot overlap")
        if gauge != 0 and balanced >= 2:
            need(gauge in single, "two balanced deletions force a single-edge switching class")
        if gauge != end-1 and antibalanced >= 2:
            need(gauge in negsingle, "two antibalanced deletions force globally negative single edge")
        w = negative_count(h, hamilton)
        for u, v in pairs:
            swap = list(range(n))
            swap[u], swap[v] = swap[v], swap[u]
            transposed = sum(bit(index, swap[a], swap[b]) for a, b in pairs if h & bit(index, a, b))
            differing = [x for x in range(n) if x not in (u, v)
                         and bool(h & bit(index, u, x)) != bool(h & bit(index, v, x))]
            difference = sum(bit(index, u, x)+bit(index, v, x) for x in differing)
            need(h ^ transposed == difference, "literal transposition gives the two-star signing")
            wd = negative_count(difference, hamilton)
            need(wd == star_count(n, len(differing)) and wd <= 2*w,
                 "two-star count and Hamming triangle inequality")
    c5_apex = sum(bit(index, j, (j+1) % 5) for j in range(5))
    need(negative_count(c5_apex, hamilton) == 20 < factorial(4), "negative C5 plus positive apex")
    return end, end*len(pairs)


def main():
    sys.stdout.reconfigure(newline="\n")
    polynomial = symbolic_controls()
    stars = direct_star_controls()
    moments = crossing_gap_controls()
    graphs = graph_rigidity_controls()
    deletions, transpositions = deletion_and_transposition_controls()
    print("status=INDEPENDENT_EXACT_CONTROLS; all-order proof audited in companion report")
    print("two_star_direct_rows="+str(stars))
    print("crossing_moment_cases="+str(moments))
    print("complete_disagreement_graph_rows=(m,all_graphs,six_family_graphs):"+str(graphs))
    print("K6_complete_switching_classes="+str(deletions)+"; transpositions="+str(transpositions))
    print("antibalanced_S7_coefficients_at_n8_plus_t="+str(polynomial))
    print("hostiles=K6_negative_C5_apex_20_below_24; K8_two_star_1296_below_1440")
    print("K9_character_array=NOT_ALLOCATED; required_K8_base=separate_bitset_audit")
    print("checks="+str(CHECKS))
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
