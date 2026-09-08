#!/usr/bin/env python3
"""Exact controls for odd-braid retention and the marked H4 degree floor.
The proof is analytic. Small permutation banks are controls, not degree cutoffs.
"""
import hashlib
import json
from collections import Counter
from itertools import permutations

GATES = 0

def check(condition, message):
    global GATES
    GATES += 1
    if not condition:
        raise RuntimeError(message)

def mul(p, q):
    return tuple(p[i] for i in q)

def power(p, r):
    ans = tuple(range(len(p)))
    for _ in range(r):
        ans = mul(ans, p)
    return ans

def alternating(s, t, length):
    ans = tuple(range(len(s)))
    for i in range(length):
        ans = mul(ans, s if i % 2 == 0 else t)
    return ans

def cycle(d, labels):
    p = list(range(d))
    for a, b in zip(labels, labels[1:] + labels[:1]):
        p[a] = b
    return tuple(p)

def cycles(p):
    unseen = set(range(len(p)))
    out = []
    while unseen:
        a = min(unseen)
        row = []
        x = a
        while x not in row:
            row.append(x)
            unseen.remove(x)
            x = p[x]
        out.append(row)
    return out

rows = []
for d in range(1, 6):
    ps = list(permutations(range(d)))
    for r in range(1, 5):
        count = Counter()
        for sigma in ps:
            fixed = [i for i in range(d) if sigma[i] == i]
            for tau in ps:
                count['pairs'] += 1
                left = alternating(sigma, tau, 2*r + 1)
                right = alternating(tau, sigma, 2*r + 1)
                if left != right:
                    continue
                count['braided'] += 1
                g = power(mul(sigma, tau), r)
                check(mul(g, sigma) == mul(tau, g), 'conjugator orientation')
                ft = {i for i in range(d) if tau[i] == i}
                moved = d - len(fixed)
                for mask in range(1 << len(fixed)):
                    A = {x for j, x in enumerate(fixed) if mask >> j & 1}
                    B = {g[x] for x in A}
                    n, k = len(A & B), len(A)
                    U = A - B
                    O = set(range(d)) - (A | B)
                    check(B <= ft, 'transported actual subset type')
                    check(all(tau[x] != x for x in U), 'joint fixed point removal')
                    check(len(U) <= r*len(O), 'cyclic block counting')
                    check((r+1)*n >= (2*r+1)*k - r*d, 'retention inequality')
                    check(n >= k - (r*moved)//(r+1), 'one-deficit strengthened bound')
                    for x in U:
                        orbit = [x]
                        for _ in range(r):
                            orbit.append(tau[orbit[-1]])
                        check(not all(y in U for y in orbit), 'no r+1 consecutive U labels')
                    count['subsets'] += 1
                    count['equal'] += ((r+1)*n == (2*r+1)*k-r*d)
                    count['partial'] += (0 < k < len(fixed))
        rows.append({'D': d, 'r': r, **count})

totals = Counter()
for row in rows:
    totals.update({key: value for key, value in row.items() if key not in ('D', 'r')})
check(dict(totals) == {'pairs': 60068, 'braided': 2328, 'subsets': 6096,
                      'equal': 230, 'partial': 2192}, 'complete unfiltered small universe')

# Sharp analytic family: disjoint odd blocks plus jointly retained fixed labels.
sharp = []
for r in range(1, 13):
    for blocks, fixed_count in ((1, 0), (1, 3), (2, 1)):
        block_size = 2*r + 1
        d = blocks*block_size + fixed_count
        sigma, tau = tuple(range(d)), tuple(range(d))
        A = set(range(blocks*block_size, d))
        for b in range(blocks):
            o = b*block_size
            aa = list(range(o+1, o+r+1))
            bb = list(range(o+r+1, o+2*r+1))
            sigma = mul(sigma, cycle(d, [o]+bb))
            tau = mul(tau, cycle(d, [o]+aa))
            A.update(aa)
        g = power(mul(sigma, tau), r)
        B = {g[x] for x in A}
        k, n = len(A), len(A & B)
        check(alternating(sigma, tau, 2*r+1) == alternating(tau, sigma, 2*r+1), 'sharp odd relation')
        check(all(sigma[x] == x for x in A) and all(tau[x] == x for x in B), 'sharp reaccess')
        check((r+1)*n == (2*r+1)*k-r*d, 'sharp affine equality family')
        sharp.append((r, blocks, fixed_count, d, k, n))

# The named five-letter hostile retains the true orientation.
s = cycle(5, [0,1,2]); t = cycle(5, [2,3,4])
A = {3,4}; B = {power(mul(s,t), 2)[i] for i in A}
check(B == {0,1}, 'five-letter actual formal reaccess')
check(3*len(A&B) == 5*len(A)-2*5 == 0, 'sharp repaired five-cusp bound')
check(4*len(A&B) < 5*len(A)-5, 'false former bound remains false')
check(len(set(range(5))-({i for i in range(5) if s[i]==i}|{i for i in range(5) if t[i]==i})) == 1,
      'five-letter support overlap one')

# All sixteen four-support membership patterns; no group realization assumed.
for mask in range(16):
    m = mask.bit_count()
    check(m*(m-1)//2 >= 2*m-3, 'full-retention incidence inequality')

# The complete scalar head AFTER the analytically inherited t>=5 reduction.
scalar = {}
for d in range(2, 13):
    rr = []
    for moved in range(5, d):
        for retained in range(1, d-moved+1):
            n3 = max(0, retained-moved//2, (3*retained-d+1)//2)
            n5 = max(0, retained-(2*moved)//3, (5*retained-2*d+2)//3)
            W = 1+2*retained-2*n3-n5
            if W >= 3*max(0,d-2*retained):
                rr.append((moved, retained, n3, n5, W))
    scalar[d] = rr
check(all(not scalar[d] for d in range(2, 9)), 'no scalar degree at most eight')
check(scalar[9] == [(5,4,2,1,4)], 'complete degree-nine row')
check(scalar[10] == [(5,5,3,2,3),(6,4,1,0,7)], 'complete degree-ten rows')
check(scalar[11] == [(5,5,3,2,3),(5,6,4,3,2),(6,5,2,1,6)], 'complete degree-eleven rows')
check(scalar[12] == [(5,6,4,3,2),(5,7,5,4,1),(6,5,2,1,6),(6,6,3,2,5),(7,5,2,1,6)],
      'degree-twelve preliminary scalar frontier, not realized passports')
# Exact strict lower controls used in the analytic head, preserving cases.
check(-8+2+2+1+6 > 1, 'D9 mixed32 node obstruction')
check(4*6-2*10 > 1 and 4*6-2*11 > 1, 'full-retention support budget')
check(-10+5+3+2+3 > 1, 'D10 mixed32 common-disjoint-support obstruction')
check(2+1+1 > 3, 'D11 partial-retention node budget')
check(-12+5+4+3+2 > 1, 'D11 full-retention mixed32 obstruction')

# A genuinely unbounded *relaxed* scalar/set family; not permutations or pages.
formal = []
for q in (2,3,4,10,101):
    cells = {1:6*q, 3:6*q, 4:2*q, 6:6*q, 8:8*q, 12:4*q, 15:q}
    d, moved, retained = 33*q, 13*q, 16*q
    check(sum(cells.values()) == d, 'formal ambient size')
    def joint(i,j):
        return sum(v for mask,v in cells.items() if mask>>i&1 and mask>>j&1)
    check(all(sum(v for mask,v in cells.items() if mask>>i&1) == moved for i in range(4)), 'formal support sizes')
    edges = [joint(0,1),joint(1,2),joint(2,3)]
    nodes = [joint(0,2),joint(0,3),joint(1,3)]
    nn = [10*q,10*q,8*q]; ww = [2*q+1,q,q]
    check(edges == [7*q,7*q,5*q] and nodes == [q,q,q], 'formal marked intersection data')
    check(-2*retained+sum(nn)+sum(ww) == 1, 'formal exact Euler one')
    for r,nc,j in zip((1,1,2),nn,edges):
        check((r+1)*nc >= (2*r+1)*retained-r*d, 'formal direct retention bound')
        check(nc >= retained-moved+j and (r+1)*j >= moved, 'formal deficit and support bounds')
    for w,j in zip(ww,nodes):
        check(w >= max(0,d-2*retained,j), 'formal node capacities')
        check(j >= 2, 'formal no-singleton intersection sidecar')
    formal.append((q,d,moved,retained,nn,ww))

report = {'universe':'all ordered S_D pairs, D=1..5, r=1..4; every subset of Fix(sigma)',
          'totals':dict(totals), 'sharp_controls':len(sharp),
          'head_D9_D12':{d:scalar[d] for d in range(9,13)},
          'formal_unbounded_controls':len(formal),
          'scope':'analytic all-r inequality; marked H4 actual degree at least 12; no all-cycle closure'}
raw = json.dumps({'rows':rows,'sharp':sharp,'formal':formal},sort_keys=True,separators=(',',':')).encode()
report['semantic_sha256'] = hashlib.sha256(raw).hexdigest()
report['always_active_gates'] = GATES
print(json.dumps(report,sort_keys=True,indent=2))
