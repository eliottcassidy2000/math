#!/usr/bin/env python3
"""Exact literal controls and finite-word checks for the paired-language proof."""
from itertools import product
import json

PAIRS = [
    (677, (3, 3), 2, (2,12,17,101,75,23,167,63), (2,12,17,165,43,71,31,119)),
    (759, (4, 3), 1, (2,44,49,133,11,89,31,119), (2,52,25,35,105,197,15,95)),
    (1547, (4, 3), 0, (2,148,97,7,45,75,27,125), (2,148,41,67,45,75,23,125)),
    (821, (5, 5), 2, (4,9,50,85,163,75,151,47), (4,9,50,85,99,139,167,31)),
]
EXPECTED_DELTAS = [(-10,-1,0,9,-12,0,-88,12,-57),
                   (0,2,0,0,0,0,-168,0,-84),
                   (0,2,0,0,0,0,-24,0,-12),
                   (2,-1,0,-3,0,0,8,0,3)]

def need(ok, message):
    if not ok:
        raise RuntimeError(message)

def remark(rows, root):
    order = [root] + [i for i in range(len(rows)) if i != root]
    return tuple(sum(((rows[u] >> v) & 1) << j for j,v in enumerate(order)) for u in order)

def profile(rows):
    n = len(rows)
    full = (1 << n)-1
    for i in range(n):
        need(rows[i] & ~full == 0 and not (rows[i] >> i & 1), 'invalid row')
        for j in range(i):
            need(((rows[i] >> j) & 1) + ((rows[j] >> i) & 1) == 1, 'invalid edge')
    ends = [[0]*n for _ in range(full+1)]
    starts = [[0]*n for _ in range(full+1)]
    counts = [0]*(full+1)
    counts[0] = 1
    for mask in range(1,full+1):
        for v in range(n):
            if not (mask >> v & 1):
                continue
            old = mask ^ (1 << v)
            if old == 0:
                ends[mask][v] = starts[mask][v] = 1
            else:
                ends[mask][v] = sum(ends[old][u] for u in range(n)
                                   if old >> u & 1 and rows[u] >> v & 1)
                starts[mask][v] = sum(starts[old][u] for u in range(n)
                                     if old >> u & 1 and rows[v] >> u & 1)
        counts[mask] = sum(ends[mask])
        need(counts[mask] == sum(starts[mask]), 'start/end totals')
    covers = [0]*(full+1)
    for mask in range(full+1):
        sub = mask
        while True:
            covers[mask] += counts[sub]*counts[mask ^ sub]
            if sub == 0:
                break
            sub = (sub-1) & mask
    L = sum(ends[s][0]*counts[full ^ s] for s in range(full+1))
    R = sum(starts[s][0]*counts[full ^ s] for s in range(full+1))
    T = sum(starts[s][0]*covers[full ^ s] for s in range(full+1))
    return (counts[full],starts[full][0],ends[full][0],L,R,
            counts[full ^ 1],covers[full],covers[full ^ 1],T)

def move(p, c):
    h,s,e,L,R,d,C,b,T = p
    if c == 'O':
        return h,s,0,L,2*R,d,2*C,2*b,3*T
    if c == 'Z':
        return h+R-e,s,d,L+b,2*R-d,d,2*C+2*T-2*L,2*b,3*T-b
    if c == 'J':
        return h+R-e,s,h-e,C-L,h+R,h,2*C+2*T-2*L,C,T+C
    raise ValueError(c)

def graph_move(rows, c):
    n = len(rows)
    result = [row | (1 << n) for row in rows] + [0]
    if c != 'O':
        result[0] ^= 1 << n
        result[n] = 1
    return remark(result,n) if c == 'J' else tuple(result)

def cover(p):
    h,s,e,L,R,d,C,b,T = p
    return L,C-2*L,T-b,b

def ordinary(p):
    h,s,e,L,R,d,C,b,T = p
    return h-e-d,R-h+e-d,e,d-e

def cover_move(p,c):
    L,A,W,b = p
    return {'O':(L,2*L+2*A,3*W+b,2*b),
            'Z':(L+b,2*A+2*W,3*W,2*b),
            'J':(L+A,2*W+2*b,W+b,2*L+A)}[c]

def ordinary_move(p,c):
    a,f,e,p = p
    return {'O':(a+e,a+2*f+e+2*p,0,e+p),
            'Z':(2*a+f+e+p,f,e+p,0),
            'J':(f+p,0,a+e+p,e)}[c]

def base(p,j):
    return p[0]+p[3]-p[1]+2*j

def main():
    symbolic = 0
    for i in range(9):
        p = tuple(int(k == i) for k in range(9))
        for c in 'OZJ':
            need(cover(move(p,c)) == cover_move(cover(p),c), 'cover identity')
            need(ordinary(move(p,c)) == ordinary_move(ordinary(p),c), 'ordinary identity')
            symbolic += 2
    for fn in (cover_move,ordinary_move):
        for i in range(4):
            for c in 'OZJ':
                need(min(fn(tuple(int(k == i) for k in range(4)),c)) >= 0, 'positive matrix')
    report = {'status':'PASS','scope':'All-order claim rests on PROOF.md; finite controls only.',
              'symbolic_basis_checks':symbolic,'pair_controls':[]}
    for idx,(X,js,root,a,b) in enumerate(PAIRS):
        rows = [remark(a,root),remark(b,root)]
        profiles = [profile(r) for r in rows]
        need(tuple(x-y for x,y in zip(*profiles)) == EXPECTED_DELTAS[idx], 'seed delta')
        for r,(p,j) in enumerate(zip(profiles,js),1):
            need(p[1] == r+3*j and 0 <= j <= 8 and base(p,j) == X, 'seed certificate')
        literal = 0
        for side in range(2):
            queue = [(rows[side],profiles[side],0)]
            while queue:
                rr,pp,depth = queue.pop()
                need(profile(rr) == pp, 'literal profile mismatch')
                literal += 1
                if depth < 2:
                    queue.extend((graph_move(rr,c),move(pp,c),depth+1) for c in 'OZJ')
        total = sync = separated = 0
        queue = [(profiles[0],profiles[1],False,0)]
        while queue:
            p,q,has_j,depth = queue.pop()
            D = base(p,js[0])-base(q,js[1])
            need((D == 0) == (not has_j), 'language classification')
            if has_j:
                need(D > 0 if idx == 3 else D < 0, 'separation sign')
                separated += 1
            else:
                sync += 1
            total += 1
            if depth < 10:
                queue.extend((move(p,c),move(q,c),has_j or c == 'J',depth+1) for c in 'OZJ')
        report['pair_controls'].append({'base':X,'profiles':profiles,
            'difference':EXPECTED_DELTAS[idx],'literal_graph_profiles':literal,
            'word_length_cap':10,'all_words':total,'synchronized':sync,'separated':separated})
    print(json.dumps(report,sort_keys=True,indent=2))

if __name__ == '__main__':
    main()
