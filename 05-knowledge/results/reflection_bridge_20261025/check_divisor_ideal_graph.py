"""Exact finite controls for divisor_ideal_graph.md; standard library only."""
from collections import Counter
from itertools import product
from math import comb
import json


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def encode(k, d):
    a, b, c = d
    return ((1 << a)-1) | (b << k) | (c << (k+1))


def decode(k, mask):
    low = mask & ((1 << k)-1)
    a = low.bit_count()
    require(low == (1 << a)-1, "not a prefix ideal")
    return a, (mask >> k) & 1, (mask >> (k+1)) & 1


def permute_bits(mask, permutation):
    return sum(((mask >> i) & 1) << j for i, j in enumerate(permutation))


def edge_image(edges, permutation):
    return frozenset(tuple(sorted((permutation[i], permutation[j])))
                     for i, j in edges)


def main():
    records = []
    pair_checks = quotient_checks = action_checks = 0
    for k in range(2, 7):
        full = list(product(range(k+1), range(2), range(2)))
        proper = [d for d in full if d not in [(0, 0, 0), (k, 1, 1)]]
        universe = (1 << (k+2))-1
        J = tuple(range(k-1, -1, -1))+(k, k+1)
        R = tuple(range(k))+(k+1, k)
        complement = lambda d: (k-d[0], 1-d[1], 1-d[2])
        exchange = lambda d: (d[0], d[2], d[1])
        chart = lambda d: (2*d[0]-k*(d[1]+d[2]), d[1]-d[2])
        for d in full:
            x = encode(k, d)
            require(decode(k, x) == d, "inverse")
            require(x.bit_count() == sum(d), "rank")
            require(encode(k, complement(d)) == universe ^ permute_bits(x, J),
                    "twisted complement")
            require(encode(k, exchange(d)) == permute_bits(x, R), "prime swap")
            require(permute_bits(permute_bits(x, J), R) ==
                    permute_bits(permute_bits(x, R), J), "commuting operations")
            action_checks += 1
            for e in full:
                y = encode(k, e)
                require(all(a <= b for a, b in zip(d, e)) == ((x & y) == x),
                        "order equivalence")
                require(encode(k, tuple(min(a, b) for a, b in zip(d, e))) == x & y,
                        "gcd/intersection")
                require(encode(k, tuple(max(a, b) for a, b in zip(d, e))) == x | y,
                        "lcm/union")
                pair_checks += 1
        require(len({chart(d) for d in proper}) == len(proper), "chart injectivity")
        fixes = dict(C=sum(complement(d) == d for d in proper),
                     Q=sum(exchange(d) == d for d in proper),
                     CQ=sum(complement(exchange(d)) == d for d in proper))
        require(fixes == dict(C=0, Q=2*k, CQ=2 if k % 2 == 0 else 0), "fixed sets")
        occupancy = lambda x: ((x & ((1 << k)-1)).bit_count(),
                               (x >> k) & 1, (x >> (k+1)) & 1)
        fibres = Counter()
        for x in range(universe+1):
            d = occupancy(x)
            fibres[d] += 1
            require(occupancy(universe ^ x) == complement(d), "quotient complement")
            require(occupancy(permute_bits(x, R)) == exchange(d), "quotient swap")
            quotient_checks += 1
        require(all(fibres[d] == comb(k, d[0]) for d in full), "binomial fibres")
        require(occupancy(1) == occupancy(2) and occupancy(1 & 2) != occupancy(1),
                "intersection hostile")
        require(universe ^ encode(k, (1, 0, 0)) not in {encode(k, d) for d in full},
                "bare complement hostile")
        records.append(dict(k=k, proper_sources=len(proper), edge_slots=k+2,
                            fixed_points=fixes,
                            rank_histogram=dict(sorted(Counter(map(sum, proper)).items())),
                            g_histogram=dict(sorted(Counter(d[1]-d[2] for d in proper).items()))))
    host = frozenset([(1, 2), (0, 3), (0, 1), (2, 3)])
    slots = [(1, 2), (0, 3), (0, 1), (2, 3)]
    Jv, Rv = (1, 0, 3, 2), (3, 2, 1, 0)
    def square(d):
        mask = encode(2, d)
        return frozenset(edge for i, edge in enumerate(slots) if mask & (1 << i))
    square_checks = 0
    for d in product(range(3), range(2), range(2)):
        cd, qd = (2-d[0], 1-d[1], 1-d[2]), (d[0], d[2], d[1])
        require(square(cd) == host - edge_image(square(d), Jv), "square actual C")
        require(square(qd) == edge_image(square(d), Rv), "square actual Q")
        square_checks += 1
    halfturn = tuple(Jv[Rv[i]] for i in range(4))
    diagonals = [(0, 2), (1, 3)]
    require(all(edge_image([edge], halfturn) == frozenset([edge]) for edge in diagonals),
            "fixed excluded diagonals")
    require(all(edge_image([edge], halfturn) != frozenset([edge]) for edge in host),
            "no fixed host edge")
    d2 = [d for d in product(range(3), range(2), range(2))
          if d not in [(0, 0, 0), (2, 1, 1)]]
    rank_partition = {frozenset(d for d in d2 if sum(d) == rank)
                      for rank in range(1, 4)}
    g_partition = {frozenset(d for d in d2 if d[1]-d[2] == g)
                   for g in range(-1, 2)}
    require(rank_partition != g_partition, "rank/g partitions must differ")
    print(json.dumps(dict(status="FINITE-EXACT", source_pair_checks=pair_checks,
                          source_action_checks=action_checks, quotient_checks=quotient_checks,
                          square_actual_vertex_checks=square_checks, rows=records),
                     indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
