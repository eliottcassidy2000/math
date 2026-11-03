#!/usr/bin/env python3
"""Finite ascending-prime live boundaries; not a minimum-pathwidth claim."""
from math import isqrt
import json


def need(ok, why):
    if not ok:
        raise RuntimeError(why)


groups = [[2,3,5,7,11,13], [17,19,23,31,47], [29,41,59,71],
          [37,43,53,61,67], [73,79,83,89,97]]
primes = [p for p in range(2, 101) if all(p % d for d in range(2, isqrt(p)+1))]
need(sorted(sum(groups, [])) == primes and len(primes) == 25, 'actual prime partition')
edges = {(p,q) for i,p in enumerate(primes) for q in primes[i+1:]
         if isqrt(p+q)**2 == p+q}
by_roots = {(p,s*s-p) for s in range(2, 15) for p in primes
            if s*s-p in primes and p < s*s-p}
need(edges == by_roots and len(edges) == 24, 'independent edge enumeration')


def ledger(vertices, relation):
    rows = []
    for k in range(1, len(vertices)):
        later = set(vertices[k:])
        live = [p for p in vertices[:k]
                if any((min(p,q), max(p,q)) in relation for q in later)]
        witnesses = [[p, min(q for q in later if (min(p,q),max(p,q)) in relation)] for p in live]
        rows.append({'after': vertices[k-1], 'before': vertices[k],
                     'live': live, 'future_neighbor_witnesses': witnesses})
    width = max(map(lambda r: len(r['live']), rows))
    return {'width': width, 'maximizing_cuts': [r for r in rows if len(r['live']) == width], 'cuts': rows}


full = ledger(primes, edges)
odd = [p for p in primes if p != 2]
odd_edges = {(p,q) for p,q in edges if p != 2}
odd_ledger = ledger(odd, odd_edges)
scaled_vertices = [49*p for p in odd]
scaled_edges = {(49*p,49*q) for p,q in odd_edges}
scaled = ledger(scaled_vertices, scaled_edges)
need(full['width'] == 8 and odd_ledger['width'] == scaled['width'] == 7, 'exact ordered widths')
for original, image in zip(odd_ledger['cuts'], scaled['cuts']):
    need(image['live'] == [49*p for p in original['live']], 'all-cut square-carry transport')
need(all(isqrt(a+b)**2 == a+b for a,b in scaled_edges), 'square roots scale by7')
print(json.dumps({'scope': 'given ascending order on all primes through100; no optimization over orders',
                  'groups': groups, 'primes': primes, 'edges': sorted(edges),
                  'full': full, 'odd': odd_ledger,
                  'square_carry': {'factor': 49, 'vertices': scaled_vertices,
                                   'ordered_width': scaled['width'],
                                   'all_cut_transport_checked': len(scaled['cuts'])}}, indent=2, sort_keys=True))
