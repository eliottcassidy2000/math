#!/usr/bin/env python3
"""Independent small trinomial controls and cheap all-shape hostile census."""
from hashlib import sha256
from math import gcd
import json

checks = cases = edges = 0
digest = sha256()

def require(ok):
    global checks
    checks += 1
    if not ok:
        raise RuntimeError('independent one-swap check failed')

def carry(w):
    ones = [i for i, bit in enumerate(w) if bit]
    return sum(2**i * 3**(len(ones)-j-1) for j, i in enumerate(ones))

def rotate(w, r):
    return w[r:]+w[:r]

for q in range(4, 15):
    for p in range(1, q):
        if gcd(q, p) != 1:
            continue
        c = tuple((k+1)*p//q-k*p//q for k in range(q))
        D = 2**q-3**p
        inverse = pow(p, -1, q)
        for i in range(q):
            ii = (i+1) % q
            if c[i] == c[ii]:
                continue
            w = list(c)
            w[i], w[ii] = w[ii], w[i]
            w = tuple(w)
            epsilon = 1 if c[i] else -1
            exceptional = ii
            t_exceptional = ii*p % q
            g = gcd(D, carry(w))
            E = {}
            for r in range(q):
                j = (i-r) % q
                require((j == q-1) == (r == exceptional))
                if r == exceptional:
                    continue
                cr = rotate(c, r)
                s = sum(cr[j+2:])
                E[r] = epsilon*2**j*3**s
                require(carry(rotate(w, r))-carry(cr) == E[r])
            retained = 0
            bound = abs(D)
            for t in range(q-1):
                r, rr = t*inverse % q, (t+1)*inverse % q
                k = ((q-1-t)*inverse % q)-1
                h = p-1-(k+1)*p//q
                M = 2**k*3**h
                require(carry(rotate(c, r))-carry(rotate(c, rr)) == M)
                require((r == exceptional or rr == exceptional)
                        == (t == t_exceptional or t+1 == t_exceptional))
                if r == exceptional or rr == exceptional:
                    continue
                T = M+E[r]-E[rr]
                require(T == carry(rotate(w, r))-carry(rotate(w, rr)))
                require(T != 0 and T % g == 0)
                bound = gcd(bound, T)
                retained += 1
            require(retained == q-1-(1 if t_exceptional in (0,q-1) else 2))
            require(bound % g == 0)
            digest.update(repr((q,p,i,g,bound,retained)).encode())
            cases += 1
            edges += retained

probe = nonunit = 0
first = maximum = None
positive_integer_candidates = []
for q in range(3,81):
    for p in range(1,q):
        if gcd(q,p) != 1:
            continue
        c = [(k+1)*p//q-k*p//q for k in range(q)]
        D = 2**q-3**p
        for i in range(q):
            ii = (i+1) % q
            if c[i] == c[ii]:
                continue
            w = c.copy()
            w[i],w[ii] = w[ii],w[i]
            B = carry(w)
            g = gcd(D,B)
            probe += 1
            if g > 1:
                nonunit += 1
                row = dict(q=q,p=p,swap_start=i,D=D,B=B,g=g,
                           word=''.join(map(str,w)))
                if first is None:
                    first = row
                if maximum is None or g > maximum['g']:
                    maximum = row
            if D > 1 and g == D:
                positive_integer_candidates.append([q,p,i])
require((probe, nonunit) == (52792,3821))
require(first == dict(q=7,p=2,swap_start=3,D=119,B=112,g=7,word='0000101'))
require(maximum['g'] == 11243 and (maximum['q'],maximum['p']) == (66,35))
require(not positive_integer_candidates)
print(json.dumps(dict(status='PASS',gate_lengths='4 through 14',
    canonical_word='zero-intercept lower mechanical word',gate_cases=cases,
    retained_edges=edges,checks=checks,certificate_digest_sha256=digest.hexdigest(),
    probe_lengths='3 through 80',probe_cases=probe,nonunit_occurrences=nonunit,
    first_nonunit=first,maximum_content=maximum,
    positive_integer_candidates_D_gt_1=positive_integer_candidates,
    boundary='Finite probe does not establish the all-q proper-divisor target.'),
    indent=2,sort_keys=True))
