#!/usr/bin/env python3
"""One complete actual 5+8 entry hostile to forcing the cross-divisor tests."""
from fractions import Fraction as F
from itertools import combinations
from math import gcd, lcm, prod
from pathlib import Path
import hashlib
import json
import subprocess

Q = 91**6
T = 4912751
GATES = 0


def check(truth, label):
    global GATES
    GATES += 1
    if not truth:
        raise RuntimeError(label)


def atlas(A, B):
    d = gcd(A, B)
    s = (A+B)//d
    if s > 356:
        return False
    p = 2
    while p*p <= s:
        exponent = 0
        while s % p == 0:
            s //= p
            exponent += 1
        if exponent and (p % 3 != 2 or exponent > 2):
            return False
        p += 1
    return s == 1 or s % 3 == 2


def graph_components(row):
    todo = set(range(len(row)))
    parts = []
    while todo:
        seed = min(todo)
        todo.remove(seed)
        seen, pending = {seed}, [seed]
        while pending:
            i = pending.pop()
            for j in list(todo):
                if atlas(row[i], row[j]):
                    todo.remove(j)
                    seen.add(j)
                    pending.append(j)
        parts.append(sorted(seen))
    return parts


def ceildiv(a, b):
    return -((-a)//b)


def support(A, B, Y):
    A, B = sorted((A, B))
    D = gcd(A, B)
    a, b = A//D, B//D
    delta = gcd(D, Y)
    c, x = D//delta, Y//delta
    check(b <= Q, 'internal primitive height')
    r = (pow(a, -1, b)*x) % b
    s = (x-a*r)//b
    lower = max(ceildiv(-Q-r,b),ceildiv(s-Q,a))
    upper = min((Q-r)//b,(s+Q)//a)
    positive = c <= Q and lower <= upper
    return dict(pair=[A,B], distinguished=Y, pair_gcd=D,
                normalized=[a,b], delta=delta, coefficient=c,
                target=x, lower=lower, upper=upper, positive=positive)


def main():
    qs = (337,343,349,355)
    C = prod(qs)
    V = tuple(sorted((C,)+tuple((356-q)*(C//q) for q in qs)))
    primes = (2,3,11,17,23,29)
    P = prod(primes)
    U = tuple(sorted((P,3*P)+tuple(P//p for p in primes)))
    row = tuple(T*v for v in V)+U
    check(V == (40341259,287243635,542783995,807423715,14321146945), 'V literal')
    check(U == (25806,32538,44022,68034,249458,374187,748374,2245122), 'U literal')
    check(gcd(*V) == gcd(*U) == gcd(*row) == 1, 'primitive normalizations')
    check(len(set(row)) == 13 and min(row) > 0, 'physical distinctness')
    check(1 not in V and 1 not in U, 'both normalized components unitless')
    check(gcd(T,P) == 1 and T % 2 == 1, 'scale coprime and odd')
    check(all(gcd(v,u) == 1 for v in V for u in U), 'disjoint prime supports')
    check(sum(row) == 78598806272076840 < Q*Q, 'full physical box')
    check(graph_components(row) == [list(range(5)),list(range(5,13))], 'actual 5+8 graph')
    for q in qs:
        check(atlas(C,(356-q)*(C//q)), 'actual sum356 star edge')
    for p,q in ((2,3),(3,17),(17,29),(2,23),(23,11)):
        check(atlas(P//p,P//q), 'actual cofactor tree edge')
    check(atlas(P//3,P) and atlas(P,3*P), 'actual multiplier edges')

    minimum_D = min(gcd(v,w) for v,w in combinations(V,2))
    maximum_sum = max((u+w)//gcd(u,w) for u,w in combinations(U,2))
    check(minimum_D == 115591, 'minimum smaller pair gcd')
    check(maximum_sum == 88, 'maximum larger primitive pair sum')
    check(T*minimum_D == 567869800841 > Q, 'all opposite coefficient gates fail')
    check(T*min(V) > Q*maximum_sum, 'all larger pair supports too short')
    records = []
    counts = []
    for same, other, kind in ((tuple(T*v for v in V),U,'coefficient'),
                              (U,tuple(T*v for v in V),'support')):
        count = 0
        for A,B in combinations(same,2):
            for Y in other:
                record = support(A,B,Y)
                check(not record['positive'], 'complete mixed support rejection')
                if kind == 'coefficient':
                    check(record['coefficient'] > Q, 'literal coefficient obstruction')
                else:
                    check(record['target'] > Q*sum(record['normalized']),
                          'literal full support obstruction')
                records.append(record)
                count += 1
        counts.append(count)
    check(counts == [80,140], 'complete support orientations')

    native = []
    endpoint = []
    for u,w in combinations(U,2):
        D = gcd(u,w)
        a,b = u//D,w//D
        R = Q*(a+b)-(a-1)*(b-1)
        for v in V:
            delta = gcd(D,T*v)
            check(delta == 1, 'no cross-divisor cancellation')
            score = F(5*delta*R,63*max(U)*v)
            check(score < 1, 'every exact native width test fails')
            native.append((score,u,w,v,D))
            if w == max(U):
                d = gcd(D,v)
                check(30*(D//d) <= Q, 'endpoint coefficient comparison does pass')
                check(63*lcm(D,v) > 5*Q, 'every endpoint lcm test fails')
                endpoint.append(lcm(D,v))
    best = max(native)
    check(best[0] == F(4730272820,108022718367), 'exact best native score')
    check(min(endpoint) == 1041046529754, 'exact endpoint minimum')
    check(5*T < 63*max(U), 'universal U-arc gluing test fails')
    check(8 < 42*max(V), 'universal V-arc gluing test fails')

    # Complete inherited profile comparison, independent of the scalar caps.
    rel = '05-knowledge/results/lrc14_joint_shadow_empty_core_next_sep06.json'
    base = Path(__file__).resolve().parent.parent
    file = base/rel
    raw = file.read_bytes() if file.exists() else subprocess.check_output(
        ['git','show','HEAD:'+rel],cwd=base)
    check(hashlib.sha256(raw).hexdigest() ==
          '935f3f687b6d7c89cc099e536f536238fd753bcc4c1747906d213cef387ca93f', 'profile pin')
    bank = {int(d): {(c,tuple(gs)) for c,gs in level['profiles']}
            for d,level in json.loads(raw)['levels'].items()}
    maxima, profiles = [], set()
    for tails in range(1,7):
        top = 1
        for body in combinations(range(13),13-tails):
            occupied = set(body)
            c = gcd(*(row[i] for i in occupied))
            gs = tuple(sorted(gcd(c,row[i]) for i in range(13) if i not in occupied))
            check((c,gs) in bank[tails], 'every inherited profile retained')
            top = max(top,c)
            if c > 1:
                profiles.add((tails,c,gs))
        maxima.append(top)
    check(maxima == [1,1,1,1,1,29], 'subset gcd maxima')
    check(profiles == {(6,p,(1,)*6) for p in primes}, 'only nontrivial profiles')

    phase = F(1403643,9825502)
    distances = [min((n*phase)%1,1-(n*phase)%1) for n in row]
    check((T*phase)%1 == F(1,2), 'exact smaller half phase')
    check(all(v % 2 for v in V), 'every smaller coordinate odd')
    check(abs(phase-F(1,7)) == F(1,14*T), 'favorable actual grid alignment')
    check(all(gcd(u,7) == 1 for u in U), 'larger clock7 phase')
    check(min(distances) == F(696962,4912751) > F(1,14), 'exact full safe phase')

    manifest = dict(Q=Q,V=V,U=U,t=T,g=1,physical_sum=sum(row),
                    minimum_V_pair_gcd=minimum_D,maximum_U_pair_sum=maximum_sum,
                    mixed_support_counts=counts,subset_gcd_maxima=maxima,
                    nontrivial_profiles=sorted(profiles),
                    best_native_score=str(best[0]),best_native_labels=best[1:],
                    minimum_endpoint_lcm=min(endpoint),phase=str(phase),
                    clearance=str(min(distances)),
                    support_digest=hashlib.sha256(json.dumps(records,sort_keys=True,
                                      separators=(',',':')).encode()).hexdigest())
    print('FINITE-EXACT actual equality entry; cross-divisor necessity REFUTED; row SAFE')
    print(json.dumps(manifest,sort_keys=True,separators=(',',':')))
    print('EXPLICIT_GATES',GATES)
    print('SEMANTIC_SHA256',hashlib.sha256(json.dumps(manifest,sort_keys=True,
                                separators=(',',':')).encode()).hexdigest())


if __name__ == '__main__':
    main()
