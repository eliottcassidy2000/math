#!/usr/bin/env python3
"""Exact, deterministic controls for three repo-derived mathematical bridges.
Run: python3 -B verify_unprompted_bridges.py
No assertion statements, random choices, or floating-point arithmetic.
"""
from collections import Counter
from fractions import Fraction as Q
from itertools import combinations
from math import comb
import json
import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ


def require(ok, text):
    if not ok:
        raise RuntimeError(text)


def rank(a):
    if not a:
        return 0
    a = [[Q(x) for x in row] for row in a]
    i = 0
    for j in range(len(a[0])):
        p = next((p for p in range(i, len(a)) if a[p][j]), None)
        if p is None:
            continue
        a[i], a[p] = a[p], a[i]
        c = a[i][j]
        a[i] = [x/c for x in a[i]]
        for h in range(i+1, len(a)):
            c = a[h][j]
            if c:
                a[h] = [x-c*y for x,y in zip(a[h],a[i])]
        i += 1
        if i == len(a):
            break
    return i


def unique_solve(a,b,n):
    a = [[Q(x) for x in row] + [Q(y)] for row,y in zip(a,b)]
    i = 0
    pivots = []
    for j in range(n):
        p = next((p for p in range(i,len(a)) if a[p][j]),None)
        if p is None:
            continue
        a[i],a[p] = a[p],a[i]
        c = a[i][j]
        a[i] = [x/c for x in a[i]]
        for h in range(len(a)):
            if h != i and a[h][j]:
                c = a[h][j]
                a[h] = [x-c*y for x,y in zip(a[h],a[i])]
        pivots.append(j)
        i += 1
    if any(all(x == 0 for x in row[:n]) and row[-1] for row in a):
        return None
    if len(pivots) != n:
        return None
    out = [Q(0)]*n
    for h,j in enumerate(pivots):
        out[j] = a[h][-1]
    return out


def product_linear(points):
    p = [Q(1)]
    for x in points:
        a = [Q(0)]*(len(p)+1)
        for j,c in enumerate(p):
            a[j] += c
            a[j+1] -= x*c
        p = a
    return p


def runs(s):
    s = sorted(s)
    return [(x, next((y for y in range(x,max(s)+2) if y not in s),max(s)+1))
            for x in s if x-1 not in s] if s else []


def decode_moments(moments,k):
    # R(z)=prod(1-start*z)/prod(1-exit*z), log R=sum D_j*z^j/j.
    d = [Q(0)] + [sum(Q(comb(j,m))*moments[m] for m in range(j))
                         for j in range(1,2*k+1)]
    c = [Q(1)]
    for n in range(1,2*k+1):
        c.append(sum(d[j]*c[n-j] for j in range(1,n+1))/n)
    for degree in range(k+1):
        a = [[c[n-j] for j in range(1,degree+1)]
             for n in range(degree+1,2*k+1)]
        b = [-c[n] for n in range(degree+1,2*k+1)]
        sol = unique_solve(a,b,degree)
        if sol is not None:
            den = [Q(1)]+sol
            num = [sum(den[j]*c[n-j] for j in range(min(n,degree)+1))
                   for n in range(degree+1)]
            return num,den
    raise RuntimeError('No bounded-run rational decoder')


def interval_checks():
    total = 0
    counts = {}
    for k in (1,2,3):
        n = 12
        signatures = {}
        checked = 0
        for mask in range(1 << n):
            s = [i+1 for i in range(n) if mask >> i & 1]
            rr = runs(s)
            if len(rr)>k:
                continue
            mom = tuple(sum(x**j for x in s) for j in range(2*k))
            require(mom not in signatures or signatures[mom] == mask,
                    'moment collision at sufficient budget')
            signatures[mom] = mask
            num,den = decode_moments(mom,k)
            require(num == product_linear(a for a,b in rr), 'wrong start polynomial')
            require(den == product_linear(b for a,b in rr), 'wrong exit polynomial')
            checked += 1
        counts[str(k)] = checked
        total += checked
    shifted_checks = 0
    for mask in range(1 << 10):
        s = [i-5 for i in range(10) if mask >> i & 1]
        rr = runs(s)
        if len(rr)>3:
            continue
        moments = [sum(x**j for x in s) for j in range(6)]
        num,den = decode_moments(moments,3)
        require(num == product_linear(a for a,b in rr), 'shifted start including zero')
        require(den == product_linear(b for a,b in rr), 'shifted exit including zero')
        # Homogenize to common degree d. Trailing zero coefficient records
        # a genuine zero endpoint; it is not discarded as a lost root.
        d = len(num)-1
        require(d == len(rr), 'zero-endpoint degree loss')
        shifted_checks += 1
    a,b = [1,5,6],[2,3,7]
    ma = [sum(x**j for x in a) for j in range(4)]
    mb = [sum(x**j for x in b) for j in range(4)]
    require(ma[:3] == mb[:3] and ma[3] != mb[3], 'hostile failed')
    require(len(runs(a)) == len(runs(b)) == 2, 'hostile run count')
    return {'universe':'all subsets of {1,...,12} with at most k runs',
            'decode_counts_by_k':counts, 'total_positive_decodes':total,
            'shifted_universe':'all subsets of {-5,...,4} with at most 3 runs',
            'shifted_and_zero_endpoint_decodes':shifted_checks,
            'three_moment_hostile':{'A':a,'B':b,'moments_A':ma,'moments_B':mb}}


def evaluations(vectors,q):
    exponents = [(a,b,q-a-b) for a in range(q+1) for b in range(q-a+1)]
    return [[x**a*y**b*z**c for a,b,c in exponents] for x,y,z in vectors]


def conductor_checks():
    cases = []
    for r in range(6,19):
        a = [(1,t,t*t) for t in range(r)]
        b = a[:-1]+[(1,r-1,(r-1)**2+1)]
        for tag,v in [('conic',a),('off_conic',b)]:
            require(all(sp.Matrix([v[i] for i in ids]).det() != 0
                        for ids in combinations(range(r),3)), 'not uniform matroid')
            hs = []
            separators = [None]*r
            for q in range((r+1)//2+1):
                e = evaluations(v,q)
                h = rank(e)
                require(h == sp.Matrix(e).rank(), 'independent rank disagreement')
                expected = min(r,2*q+1) if tag == 'conic' or q<2 else min(r,2*q+2)
                require(h == expected,'Hilbert formula failed')
                hs.append(h)
                for i in range(r):
                    ei = [int(j == i) for j in range(r)]
                    if separators[i] is None and rank([row+[x] for row,x in zip(e,ei)]) == h:
                        separators[i] = q
            defect = sum(r-h for h in hs[1:])
            expected_defect = (r-2)**2//4
            expected_sep = [(r-1+1)//2]*r
            if tag == 'off_conic':
                expected_defect -= (r-4)//2
                expected_sep = [(r-2+1)//2]*(r-1)+[2]
            require(defect == expected_defect,'defect formula')
            require(separators == expected_sep,'branch conductor formula')
            cases.append({'r':r,'configuration':tag,'hilbert_ranks':hs,
                          'seminormal_defect':defect,'branch_conductor_exponents':separators})
    a = [(1,t,t*t) for t in range(6)]
    lam = [-1,5,-10,10,-5,1]
    e = evaluations(a,2)
    require(all(sum(lam[i]*e[i][j] for i in range(6)) == 0 for j in range(len(e[0]))),
            'quadratic fifth-difference relation')
    return {'r_range':[6,18], 'all_triples_independent':True,
            'quadratic_relation_r6':lam, 'cases':cases}


EDGES = sorted([(2,13),(2,17),(2,19),(2,21),(3,11),(3,16),
                (7,22),(10,11),(16,17),(18,19),(21,22),
                (2,10),(3,22),(7,18),(10,16),(10,22),(11,13),
                (13,18),(13,22),(16,21),(17,22),(19,22)])
SMALL = {2,11,16,18,22}


def laplacian(extra):
    vertices = sorted(set(v for e in EDGES for v in e))
    index = {v:i for i,v in enumerate(vertices)}
    m = sp.zeros(len(vertices))
    for u,v in EDGES:
        i,j = index[u],index[v]
        w = 1+int((u,v) in extra)
        m[i,i] += w; m[j,j] += w; m[i,j] -= w; m[j,i] -= w
    return vertices,m


def critical_checks():
    degree = Counter(v for e in EDGES for v in e)
    odd = {v for v,d in degree.items() if d%2}
    require(len(EDGES) == 22 and len(degree) == 12,'core sizes')
    require(all((u in SMALL) != (v in SMALL) for u,v in EDGES),'bipartition')
    swap = lambda e:tuple(sorted(21 if v==17 else 17 if v==21 else v for v in e))
    require({swap(e) for e in EDGES} == set(EDGES),'involution')
    repairs = {}; symmetric = {}
    for k in range(7):
        rr = []
        ss = []
        for es in combinations(EDGES,k):
            dd = Counter(v for e in es for v in e)
            if {v for v,d in dd.items() if d%2} != odd:
                continue
            rr.append(es)
            if {swap(e) for e in es} == set(es):
                ss.append(es)
        repairs[k] = rr; symmetric[k] = ss
    require(not any(repairs[k] for k in range(4)) and len(repairs[4]) == 2,'four-repair census')
    require(not any(symmetric[k] for k in range(6)) and len(symmetric[6]) == 27,'symmetric census')
    chosen = [(2,17),(3,11),(18,19),(21,22)]
    chosen_sym = [(2,10),(2,17),(2,21),(3,11),(10,22),(18,19)]
    smith = {}
    for tag,es in [('original',[]),('optimal',chosen),('symmetric',chosen_sym)]:
        vv,L = laplacian(es)
        c = sp.Matrix([int(v in SMALL) for v in vv])
        harmonic = all(int(x)%2 == 0 for x in L*c)
        require(harmonic == (tag != 'original'),'character criterion')
        reduced = L[:-1,:-1]
        snf = smith_normal_form(reduced,domain=ZZ)
        inv = [abs(int(snf[i,i])) for i in range(len(vv)-1)]
        require(sp.prod(inv) == abs(int(reduced.det())),'Smith/determinant')
        smith[tag] = {'nonunit_invariants':[x for x in inv if x != 1],
                      'group_order':abs(int(reduced.det())), 'harmonic_cut':harmonic}
    require(smith['original']['nonunit_invariants'] == [74748], 'inherited core checksum')
    return {'degrees':dict(sorted(degree.items())), 'odd_vertices':sorted(odd),
            'optimal_cost':4,'optimal_repairs':repairs[4],
            'involution_preserving_cost':6,'involution_preserving_repair_count':27,
            'symmetric_witness':chosen_sym,'smith_controls':smith}


if __name__ == '__main__':
    result = {'interval_moments':interval_checks(),
              'conductor_veronese':conductor_checks(),
              'critical_character_repair':critical_checks(),
              'status':'PASS: all exact checks completed'}
    print(json.dumps(result,indent=2,sort_keys=True))
