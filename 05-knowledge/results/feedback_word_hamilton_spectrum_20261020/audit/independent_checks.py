#!/usr/bin/env python3
"""Clean-room, bounded proof controls; imports no producer/repository code."""
from itertools import product, permutations, combinations
from collections import Counter
from math import comb, factorial
import json


def check(test, context):
    if not test:
        raise RuntimeError(context)


def graph(word, k, code):
    # Feedback vertices first, chain vertices last (different producer labels).
    n = k + len(word)
    adj = [[False] * n for _ in range(n)]
    for bit, (i, j) in enumerate(combinations(range(k), 2)):
        adj[i][j] = bool(code & (1 << bit))
        adj[j][i] = not adj[i][j]
    for r, a in enumerate(word):
        for i in range(k):
            adj[i][k+r] = bool(a & (1 << i))
            adj[k+r][i] = not adj[i][k+r]
        for s in range(r+1, len(word)):
            adj[k+r][k+s] = True
    return adj


def literal(adj):
    return sum(all(adj[p[i]][p[i+1]] for i in range(len(p)-1))
               for p in permutations(range(len(adj))))


def gap_assignments(word, k, code):
    adj = graph(word, k, code)
    answer = 0
    for sigma in permutations(range(k)):
        for f in product(range(k+1), repeat=len(word)):
            p = []
            for g in range(k+1):
                p += [k+i for i in range(len(word)) if f[i] == g]
                if g < k:
                    p.append(sigma[g])
            answer += all(adj[p[i]][p[i+1]] for i in range(len(p)-1))
    return answer


def inversion_count(word):
    return 1 + sum((2 if j-i == 1 else 2**(j-i-1))
                   for i in range(len(word)) for j in range(i+1, len(word))
                   if word[i] and not word[j])


def projector(s, a, j, k):
    # Here feedback permutation is identity: all others are bit relabelings.
    if s is None:
        return None
    if s[j] == 0 and j and not (a & (1 << (j-1))):
        return None
    z = list(s)
    z[j] = 1 if j == k or not (a & (1 << j)) else 2
    return tuple(z)


def apply_letter(row, a, k):
    out = Counter()
    for s, c in row.items():
        for j in range(k+1):
            t = projector(s, a, j, k)
            if t is not None:
                out[t] += c
    return +out


def surj(r, t):
    return sum((-1)**j * comb(t, j) * (t-j)**r for j in range(t+1))


def run_expansion(s, a, k, r):
    row = Counter()
    for mask in range(1 << (k+1)):
        t = s
        for j in range(k+1):
            if mask & (1 << j):
                t = projector(t, a, j, k)
        if t is not None:
            row[t] += surj(r, mask.bit_count())
    return +row


def main():
    result = {}
    literal_controls = 0
    for k, max_m in [(0, 5), (1, 5), (2, 4), (3, 2)]:
        for code in range(1 << (k*(k-1)//2)):
            for m in range(max_m+1):
                for w in product(range(1 << k), repeat=m):
                    h = literal(graph(w, k, code))
                    check(h == gap_assignments(w, k, code), ('literal gap bijection', k, code, w))
                    if k == 1:
                        check(h == inversion_count(w), ('one-apex inversion', w))
                    literal_controls += 1
    result['literal_permutation_vs_literal_gap_assignments'] = literal_controls

    extrema = []
    for m in range(2, 11):
        for r in range(1, m):
            words = [(1,)+w+(0,) for w in product(range(2), repeat=m-2)
                     if sum(w) == r-1]
            vals = {w: inversion_count(w) for w in words}
            mn = (1,)+(0,)*(m-r-1)+(1,)*(r-1)+(0,)
            mx = (1,)*r+(0,)*(m-r)
            check([w for w in words if vals[w] == min(vals.values())] == [mn], ('unique min', m, r))
            check([w for w in words if vals[w] == max(vals.values())] == [mx], ('unique max', m, r))
        vals = {w: inversion_count(w) for w in ((1,)+a+(0,) for a in product(range(2), repeat=m-2))}
        extrema.append({'n': m+1, 'minimum': min(vals.values()), 'maximum': max(vals.values()),
                        'minimizing_words': [''.join(map(str,w)) for w,v in vals.items() if v == min(vals.values())]})
    check(extrema[2]['minimizing_words'] == ['1000','1010','1110'], 'order-five boundary tie')
    result['fixed_order_extrema'] = extrema

    h13_rows = (22,12,24,17,2)
    h13_adj = [[bool(h13_rows[i] & (1 << j)) for j in range(5)] for i in range(5)]
    check(literal(h13_adj) == 13, 'unrestricted H13 literal witness')
    for i in range(5):
        seen = {i}
        for _ in range(5):
            seen |= {j for u in seen for j in range(5) if h13_adj[u][j]}
        check(len(seen) == 5, 'H13 witness strong')
    result['unrestricted_H13_witness'] = dict(rows=h13_rows, literal_H=13, strong=True)

    idem = comm = expansion = annihilators = 0
    for k in range(5):
        states = list(product(range(3), repeat=k+1))
        for a in range(1 << k):
            for s in states:
                for i in range(k+1):
                    p = projector(s,a,i,k)
                    check(projector(p,a,i,k) == p, ('idempotent',k,a,s,i))
                    idem += 1
                    for j in range(k+1):
                        check(projector(p,a,j,k) == projector(projector(s,a,j,k),a,i,k), ('commutation',k,a,s,i,j))
                        comm += 1
                # Exact integer polynomial product (M-jI), no eigenvalue numerics.
                row = {s: 1}
                for j in range(k+2):
                    out = Counter()
                    for t,c in row.items():
                        for g in range(k+1):
                            dest = projector(t,a,g,k)
                            if dest is not None:
                                out[dest] += c
                        out[t] -= j*c
                    row = {t:c for t,c in out.items() if c}
                check(not row, ('annihilator',k,a,s))
                annihilators += 1
            probes = set([states[0], states[-1], (1,)*(k+1), tuple(j%3 for j in range(k+1))])
            for s in probes:
                row = Counter({s:1})
                for r in range(8):
                    check(row == run_expansion(s,a,k,r), ('positive expansion',k,a,s,r))
                    expansion += 1
                    row = apply_letter(row,a,k)
    result['projector_checks'] = dict(k_range=[0,4], idempotence=idem, commutation=comm,
                                      annihilator_columns=annihilators, stirling_rows=expansion)

    rows = (502,508,504,497,160,192,144,256,112)
    adj = [[bool(rows[i] & (1 << j)) for j in range(9)] for i in range(9)]
    check(all(not adj[i][i] for i in range(9)) and all(adj[i][j] != adj[j][i] for i,j in combinations(range(9),2)), 'hostile is tournament')
    ts = [frozenset(t) for t in combinations(range(9),3) if all(sum(adj[x][y] for y in t) == 1 for x in t)]
    check(ts == list(map(frozenset,[(0,1,3),(0,2,3),(4,5,6),(4,7,8),(5,7,8),(6,7,8)])), 'hostile exact triangles')
    check(any(x.isdisjoint(y) for x,y in combinations(ts,2)) and
          not any(x.isdisjoint(y) and x.isdisjoint(z) and y.isdisjoint(z)
                  for x,y,z in combinations(ts,3)), 'hostile triangle matching two')
    remainder = {2,7,8}
    poison = [(u,v) for u in sorted(remainder) for v in sorted(remainder-{u})
              if all(v in t for t in ts if u in t)]
    check(poison == [(7,8),(8,7)], 'hostile poisoning arcs')
    safe = [v for v in range(9) if all(any(w in t and v not in t for t in ts) for w in range(9) if w != v)]
    check(2 in safe, 'safe deletion remains true')
    check(literal(adj) == 75, 'hostile literal H')
    result['poisoning_hostile'] = dict(triangles=[sorted(t) for t in ts], safe_deletions=safe, literal_H=75,
                                      maximum_triangle_matching=2, poisoned_pair=[7,8])

    # The prime-counting argument's prime-power valuation bound, in small controls.
    prime_checks = 0
    for n in range(1,101):
        primes = [p for p in range(2,2*n+1) if all(p%d for d in range(2,int(p**0.5)+1))]
        c = comb(2*n,n)
        check(c*(2*n+1) >= 4**n, ('central binomial lower',n))
        for p in primes:
            v = 0; rem=c
            while rem%p == 0:
                rem //= p; v += 1
            check(p**v <= 2*n, ('prime power bound',n,p,v))
            prime_checks += 1
        check(c <= (2*n)**len(primes), ('prime count lower bound',n))
    result['elementary_prime_lower_bound_controls'] = prime_checks
    result['status'] = 'PASS'
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
