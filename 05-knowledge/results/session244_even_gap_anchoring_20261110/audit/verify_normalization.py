#!/usr/bin/env python3
"""Independent adjacent-parent normalization tests on every history through order7."""
import argparse
import itertools
import json
import math


def need(test, message):
    if not test:
        raise RuntimeError(message)


def presentation(a, order):
    levels = [0]*len(a)
    edges = 0
    for i, v in enumerate(order):
        parents = [u for u in order[:i] if a[v] >> u & 1]
        need(len(parents) <= 1, 'not a forest presentation')
        if parents:
            levels[v] = 1+levels[parents[0]]
            edges += 1
    return levels, edges


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--fail', action='store_true')
    args = ap.parse_args()
    need(not args.fail, 'intentional failure control')
    rows = []
    for n in range(1,8):
        cases, swaps_total, maximum_swaps = 0,0,0
        for choice in itertools.product(*(range(v+1) for v in range(1,n))):
            a = [0]*n
            for v,q in enumerate(choice,1):
                for u in range(v):
                    if q == u+1:
                        a[v] |= 1 << u
                    else:
                        a[u] |= 1 << v
            order = list(range(n))
            depth, edges = presentation(a,order)
            original_edges = edges
            swaps = 0
            while True:
                where = next((i for i in range(1,n) if a[order[i]] >> order[i-1] & 1), None)
                if where is None:
                    break
                old_depth, old_edges = depth, edges
                order[where-1],order[where] = order[where],order[where-1]
                depth, edges = presentation(a,order)
                need(edges == old_edges-1, 'exactly one reversal disappears')
                need(all(x<=y for x,y in zip(depth,old_depth)), 'depth never increases at any vertex')
                swaps += 1
            need(swaps == original_edges-edges and swaps <= n-1, 'termination bound')
            need(all(not(a[order[i]] >> order[i-1] & 1) for i in range(1,n)), 'normalized boundary')
            cases += 1
            swaps_total += swaps
            maximum_swaps = max(maximum_swaps,swaps)
        need(cases == math.factorial(n), 'exhaustive original histories')
        rows.append({'n':n,'all_histories':cases,'total_swaps':swaps_total,'maximum_swaps':maximum_swaps,
                     'normalized_history_universe_size':math.factorial(n-1)})
    print(json.dumps({'status':'PASS','method':'direct graph order swaps, every vertex depth and every reversal recounted after each swap; no producer imports',
                      'orders':rows,'scope':'Finite controls of the separately proved all-order graph-preserving normalization.'},indent=2,sort_keys=True))


if __name__ == '__main__':
    main()
