#!/usr/bin/env python3
"""Permutation-only audit of seven-coordinate deletion interfaces; no producer imports."""
import argparse
import itertools
import json


def need(test, message):
    if not test:
        raise RuntimeError(message)


def counts(a, cut):
    n = len(a)
    size = 1 << n
    H, starts, ends = [0]*size, [0]*size, [0]*size
    H[0] = 1
    for k in range(1, n+1):
        for path in itertools.permutations(range(n), k):
            if all(a[u] >> v & 1 for u, v in zip(path, path[1:])):
                mask = sum(1 << v for v in path)
                H[mask] += 1
                starts[mask] += bool(cut >> path[0] & 1)
                ends[mask] += bool(cut >> path[-1] & 1)
    return H, starts, ends


def interface(a, cut):
    H, starts, ends = counts(a, cut)
    full = (1 << len(a))-1
    return (H[full], sum(H[s]*H[full^s] for s in range(full+1)),
            sum(H[s]*starts[full^s] for s in range(full+1)),
            sum(ends[s]*H[full^s] for s in range(full+1)),
            sum(ends[s]*starts[full^s] for s in range(full+1)),
            starts[full], ends[full])


def build(a, r, word):
    a = list(a)
    for c in word:
        n = len(a)
        a.append(0 if c == 'O' else 1 << r)
        for i in range(n):
            if c == 'O' or i != r:
                a[i] |= 1 << n
        if c == 'J':
            r = n
    keep = [v for v in range(len(a)) if v != r]
    b = [sum(1 << j for j, w in enumerate(keep) if a[v] >> w & 1) for v in keep]
    cut = sum(1 << j for j, w in enumerate(keep) if a[r] >> w & 1)
    return b, cut


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--fail', action='store_true')
    args = parser.parse_args()
    need(not args.fail, 'intentional failure control')
    records = []
    for a, r in [([0],0), ([2,0],0), ([2,0],1), ([2,4,1],0)]:
        H, starts, ends = counts(a, 1 << r)
        full = (1 << len(a))-1
        D = H[full]
        U = sum(starts[s]*H[full^s] for s in range(full+1))
        V = sum(ends[s]*H[full^s] for s in range(full+1))
        Q = sum(H[s]*H[full^s] for s in range(full+1))
        W = 0
        for x in range(full+1):
            rest = full^x
            y = rest
            while True:
                W += starts[x]*H[y]*H[rest^y]
                if not y:
                    break
                y = (y-1)&rest
        lower, upper = [build(a,r,w) for w in ('OJJJJOZ','OJOZJOZ')]
        lo, hi = [interface(*module) for module in (lower, upper)]
        need(tuple(y-x for x,y in zip(lo,hi)) == (0,0,0,-12*W,-2*U,0,0), 'seven-coordinate difference')
        need(lo[0] == 3*(D+2*U) and lo[1] == 4*(24*Q+42*W-6*V), 'h/q identities')
        need(lo[2] == 18*D+44*U and lo[6] == 0, 'A/E identities')
        need(lo[4] == 4*D+10*U and hi[4] == 4*D+8*U, 'individual T')
        n = len(a)
        perm = list(range(len(lower[0])))
        perm[n+2],perm[n+3] = perm[n+3],perm[n+2]
        need(all(bool(lower[0][u] >> v & 1) == bool(upper[0][perm[u]] >> perm[v] & 1)
                 for u in range(len(perm)) for v in range(len(perm))), 'literal deletion-graph isomorphism')
        E = ends[full]
        raw = [interface(*build(a,r,w)) for w in ('JJJJ','JOZJ')]
        need(tuple(y-x for x,y in zip(*raw)) == (0,0,E,-4*W,0,0,-U), 'raw singleton-cut interface')
        need(raw[0][4] == raw[1][4] == 0, 'singleton-cut T is zero')
        escaped = [interface(*build(a,r,w)) for w in ('JJJJOZ','JOZJOZ')]
        need(tuple(y-x for x,y in zip(*escaped)) == (0,0,4*E,-4*W,E-U,0,0), 'bounded-error seven-coordinate difference')
        need(escaped[0][0] == 3*D-E+3*U and escaped[0][1] == 4*(12*Q+14*W-6*V), 'bounded-error h/q')
        records.append({'source':a,'root':r,'D':D,'E':E,'U':U,'Q':Q,'V':V,'W':W,'lower':lo,'upper':hi,
                        'raw_lower':raw[0], 'raw_upper':raw[1],
                        'bounded_error_lower':escaped[0], 'bounded_error_upper':escaped[1]})
    print(json.dumps({'status':'PASS','method':'permutations of all subsets, explicit newborn graphs and deletion isomorphism; no producer imports',
                      'source_cases':records},indent=2,sort_keys=True))


if __name__ == '__main__':
    main()
