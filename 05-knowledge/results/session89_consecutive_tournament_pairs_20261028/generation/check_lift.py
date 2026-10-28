#!/usr/bin/env python3
"""Finite full-bank audit of the consecutive-pair/self-converse lift.

Through n=3 (lift order 7): compare image bijectively with the entire ordinary
S88 bank filtered by full-order converse canonicalization, not a count lookup.
"""
from collections import Counter
from hashlib import sha256
import argparse
import json
from time import perf_counter
import check_pairs as pairs


def lift(n, tcode, scode, mask):
    N = 2*n+1
    rows = [0]*N
    p = lambda i: i+1
    q = lambda i: n+i+1
    def orient(i, j, bit):
        rows[i if bit else j] |= 1 << (j if bit else i)
    for i in range(n):
        orient(p(i), q(i), 1)
        orient(0, p(i), (mask >> i) & 1)
        orient(0, q(i), 1-((mask >> i) & 1))
    k = 0
    for i in range(n):
        for j in range(i+1, n):
            t, s = (tcode >> k) & 1, (scode >> k) & 1
            orient(p(i), p(j), t)
            orient(q(i), q(j), 1-t)
            orient(p(i), q(j), t ^ s)
            orient(p(j), q(i), t ^ s)
            k += 1
    sigma = (0,) + tuple(q(i) for i in range(n)) + tuple(p(i) for i in range(n))
    return tuple(rows), sigma


def converse(rows):
    n = len(rows)
    return tuple(((1 << n)-1) ^ (1 << i) ^ rows[i] for i in range(n))


def audit(old, block, max_n=3):
    ordinary = {0: ((),)}
    bank = [()]
    for size in range(1, 2*max_n+2):
        result = block.build(bank, size-1, 'blockfill')
        bank = list(result['bank'])
        ordinary[size] = tuple(bank)
    records = []
    for n in range(max_n+1):
        start = perf_counter()
        N = 2*n+1
        prep = pairs.prepare(n, old, block)
        paired, _ = pairs.core_first(prep, old)
        engine = old.ParentCanonicalizer()
        target = {engine.maps(W)[0] for W in ordinary[N]
                  if engine.maps(W)[0] == engine.maps(converse(W))[0]}
        images = {}
        profile = Counter()
        for key, (rooted_key, group) in paired.items():
            W, sigma = lift(n, *key)
            pairs.need(all(((W[sigma[i]] >> sigma[j]) & 1) == (1-((W[i] >> j) & 1))
                           for i in range(N) for j in range(N) if i != j), 'lift not anti-involutive')
            code, maps = engine.maps(W)
            pairs.need(code not in images, 'lift identified two aligned pair classes')
            images[code] = key
            base_inverse = pairs.inverse(maps[0])
            aut = tuple(sorted(tuple(p[base_inverse[i]] for i in range(N)) for p in maps))
            centralizer = tuple(g for g in aut if all(g[sigma[i]] == sigma[g[i]] for i in range(N)))
            lifted_pair_group = tuple(sorted((0,) + tuple(1+p[i] for i in range(n))
                                             + tuple(1+n+p[i] for i in range(n)) for p in group))
            pairs.need(centralizer == lifted_pair_group, 'full anti-involution centralizer mismatch')
            anti_involutions = sum(all(g[sigma[g[sigma[i]]]] == i for i in range(N)) for g in aut)
            pairs.need(anti_involutions * len(centralizer) == len(aut),
                       'anti-involution conjugacy count mismatch')
            profile[(len(group), len(aut), anti_involutions)] += 1
        pairs.need(set(images) == target, 'lift image differs from full self-converse bank')
        records.append({'n': n, 'lift_order': N, 'aligned_pairs': len(paired),
                        'ordinary_tournaments': len(ordinary[N]),
                        'self_converse_tournaments': len(target),
                        'full_image_equals_full_self_converse_bank': True,
                        'every_lift_anti_involutive': True,
                        'every_full_pair_group_equals_lift_aut_centralizer': True,
                        'image_key_dictionary_sha256': pairs.digest(tuple(sorted(images.items()))),
                        'group_profiles': [{'pair_aut': a, 'lift_aut': b, 'anti_involutions': c, 'classes': z}
                                           for (a,b,c),z in sorted(profile.items())],
                        'seconds': perf_counter()-start})
    return {'status': 'PASS', 'scope': 'n<=3, entire ordinary lift bank through order7; no order9 bank',
            'results': records, 'semantic_sha256': pairs.semantic_hash(records)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--baseline-root')
    args = parser.parse_args()
    old, block = pairs.load_baseline(args.baseline_root)
    return audit(old, block)

if __name__ == '__main__':
    print(json.dumps(main(), indent=2, sort_keys=True))
