"""
Compute edge orbits via Burnside for tournament flip graph.
Compare: edge orbits vs quotient edges (pairs of vertex orbits connected by an edge).
"""
import itertools
from collections import defaultdict
import math

def canonical_form(bits, n, pairs, pair_index):
    best = bits
    for perm in itertools.permutations(range(n)):
        new_bits = 0
        for k, (i, j) in enumerate(pairs):
            pi, pj = perm[i], perm[j]
            if pi < pj:
                orig_k = pair_index[(pi, pj)]
                if bits & (1 << orig_k):
                    new_bits |= (1 << k)
            else:
                orig_k = pair_index[(pj, pi)]
                if not (bits & (1 << orig_k)):
                    new_bits |= (1 << k)
        if new_bits < best:
            best = new_bits
    return best

for n in [3, 4, 5]:
    print(f"\n=== n = {n} ===")
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    m = len(pairs)
    pair_index = {p: k for k, p in enumerate(pairs)}
    total = 2 ** m

    canon_map = {}
    orbit_sizes = defaultdict(int)
    for bits in range(total):
        c = canonical_form(bits, n, pairs, pair_index)
        canon_map[bits] = c
        orbit_sizes[c] += 1

    iso_classes = sorted(orbit_sizes.keys())

    # Quotient edges
    quotient_edges = set()
    for bits in range(total):
        c1 = canon_map[bits]
        for k in range(m):
            flipped = bits ^ (1 << k)
            c2 = canon_map[flipped]
            if c1 != c2:
                quotient_edges.add((min(c1, c2), max(c1, c2)))

    # Burnside edge orbit count
    edge_fix_sum = 0
    self_fix_sum = 0
    cross_fix_sum = 0

    for perm in itertools.permutations(range(n)):
        perm_action = {}
        for bits in range(total):
            new_bits = 0
            for kk, (i, j) in enumerate(pairs):
                pi, pj = perm[i], perm[j]
                if pi < pj:
                    orig_k = pair_index[(pi, pj)]
                    if bits & (1 << orig_k):
                        new_bits |= (1 << kk)
                else:
                    orig_k = pair_index[(pj, pi)]
                    if not (bits & (1 << orig_k)):
                        new_bits |= (1 << kk)
            perm_action[bits] = new_bits

        for bits in range(total):
            for k in range(m):
                flipped = bits ^ (1 << k)
                if flipped > bits:  # each edge once
                    pb = perm_action[bits]
                    pf = perm_action[flipped]
                    if (pb == bits and pf == flipped) or (pb == flipped and pf == bits):
                        edge_fix_sum += 1
                        c1 = canon_map[bits]
                        c2 = canon_map[flipped]
                        if c1 == c2:
                            self_fix_sum += 1
                        else:
                            cross_fix_sum += 1

    nfact = math.factorial(n)
    edge_orbits = edge_fix_sum // nfact
    self_orbits = self_fix_sum // nfact
    cross_orbits = cross_fix_sum // nfact

    print(f"  Vertex orbits: {len(iso_classes)}")
    print(f"  Quotient edges: {len(quotient_edges)}")
    print(f"  Total edge orbits (Burnside): {edge_orbits}")
    print(f"  Self edge orbits: {self_orbits}")
    print(f"  Cross edge orbits: {cross_orbits}")
    if len(quotient_edges) > 0:
        print(f"  Ratio cross_orbits/quotient_edges: {cross_orbits / len(quotient_edges):.4f}")
