#!/usr/bin/env python3
"""Independent literal count and forced-depth proof for the order-12 exception."""
import argparse
from functools import lru_cache
import json


def need(test, message):
    if not test:
        raise RuntimeError(message)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--fail', action='store_true')
    args = parser.parse_args()
    need(not args.fail, 'intentional failure control')
    n, code = 12, 25112126
    choices, a, depth = [0]*n, [0]*n, [0]*n
    c = code
    for v in range(n, 1, -1):
        choices[v-1], c = c % v, c // v
    need(c == 0, 'factorial code')
    for v in range(1, n):
        q = choices[v]
        depth[v] = 1+depth[q-1] if q else 0
        for u in range(v):
            if q == u+1:
                a[v] |= 1 << u
            else:
                a[u] |= 1 << v
    need(a == [3910,2044,3576,3057,4001,4033,3728,3841,3648,3076,2056,2], 'literal masks')
    full = (1 << n)-1

    @lru_cache(None)
    def ends(mask, last):
        if mask == 1 << last:
            return 1
        rest = mask^(1 << last)
        return sum(ends(rest, u) for u in range(n) if rest >> u & 1 and a[u] >> last & 1)

    end_count = sum(ends(full, v) for v in range(n))
    starts = [[0]*n for _ in range(full+1)]
    for mask in range(1, full+1):
        for v in range(n):
            if mask >> v & 1:
                rest = mask^(1 << v)
                starts[mask][v] = 1 if not rest else sum(starts[rest][u] for u in range(n) if rest >> u & 1 and a[v] >> u & 1)
    start_count = sum(starts[full])
    need(end_count == start_count == 97309, 'two independent endpoint directions')
    expected = [(11,1),(10,3),(9,2),(8,6),(7,0),(6,4),(5,0),(4,0)]
    mask = full
    for v, parent in expected:
        candidates = [u for u in range(n) if mask >> u & 1 and (a[u]&mask).bit_count() <= 1]
        need(candidates == [v] and a[v]&mask == 1 << parent, 'forced reverse deletion')
        mask ^= 1 << v
    need(max(depth) == 3, 'attaining depth')
    # All reverse orders contain 8->6->4->0; no nonnegative depth can make
    # that chain have maximum depth less than three.
    table = {(1,0): 1}
    for size in range(1,n):
        nxt = {}
        for (roots, ones), value in table.items():
            for state, multiplicity in [((roots+1,ones),1),((roots,ones+1),roots),((roots,ones),ones)]:
                nxt[state] = nxt.get(state,0)+value*multiplicity
        table = nxt
    need(sum(table.values()) == 123489888, 'independent shallow-universe size')
    print(json.dumps({'status':'PASS','order':n,'code':code,'H':end_count,'masks':a,
        'choices':choices[1:],'depths':depth,'minimum_presentation_depth_of_this_graph':3,
        'forced_reverse_deletions_vertex_parent':expected,'forced_depth_three_chain':[8,6,4,0],
        'depth_at_most_two_order12_histories':sum(table.values()),
        'scope':'Exact positive witness and graph-depth proof. Absence of every other shallow graph with H97309 requires the separate exhaustive scalar census.'},indent=2,sort_keys=True))


if __name__ == '__main__':
    main()
