#!/usr/bin/env python3
"""Independent hostile checks. Does not import or alter the audited packet.

Coordinates are exponent residues, with -1 (never 0) the absorbing zero.
Support-wise canonical generation, full lift-pair checks, explicit cycle
polynomials, and cached functional-graph decomposition are independent paths.
"""
from collections import Counter
from itertools import combinations, permutations, product
from math import gcd, lcm
from pathlib import Path
import json
import random

OUT = Path(__file__).resolve().parent
STATS = Counter()

def check(condition, detail):
    if not condition:
        raise RuntimeError(detail)

def pts(m, k):
    result = []
    for size in range(1, k + 1):
        for support in combinations(range(k), size):
            for phase in product(range(m), repeat=size - 1):
                p = [-1] * k
                for i, a in zip(support, (0,) + phase):
                    p[i] = a
                result.append(tuple(p))
    return result

def project(x, m):
    first = next(a for a in x if a >= 0)
    return tuple(-1 if a < 0 else (a - first) % m for a in x)

def raw(x, m, b, perm, r):
    inverse = [perm.index(i) for i in range(len(perm))]
    return tuple(-1 if x[inverse[i]] < 0 else (b[i] + r*x[inverse[i]]) % m
                 for i in range(len(x)))

def cycles(perm):
    seen = set()
    out = []
    for start in range(len(perm)):
        if start in seen:
            continue
        orbit, i = [], start
        while i not in seen:
            seen.add(i)
            orbit.append(i)
            i = perm[i]
        out.append(orbit)
    return out

def polynomial_count(m, k, b, perm, r, rotations=0):
    cs = cycles(perm)
    data = []
    for cyc in cs:
        h = rotations % len(cyc)
        cyc = cyc[h:] + cyc[:h]
        ell = len(cyc)
        # r^0 is the constant polynomial coefficient, including at r=0.
        weights = [1]
        for _ in range(ell-1):
            weights.append(weights[-1]*r)
        A = sum(weights)
        B = sum(weights[ell-j] * b[cyc[j % ell]] for j in range(1,ell+1))
        data.append((gcd(m, r**ell-1), A, B))
    numerator = 0
    for scalar in range(m):
        term = 1
        for g,A,B in data:
            term *= 1 + g * ((A*scalar-B) % g == 0)
        numerator += term-1
    check(numerator % m == 0, ('nonintegral',m,k,b,perm,r))
    return numerator//m

def lifted_pair_audit(m,k,b,perm,r,points):
    for x in points:
        by_scalar = Counter()
        for shift in range(m):
            lift = tuple(-1 if a < 0 else (a+shift)%m for a in x)
            y = raw(lift,m,b,perm,r)
            for scalar in range(m):
                if y == tuple(-1 if a < 0 else (a+scalar)%m for a in lift):
                    by_scalar[scalar] += 1
        fixed = project(raw(x,m,b,perm,r),m) == x
        check(sum(by_scalar.values()) == m*fixed, ('lift cardinality',m,k,b,perm,r,x,by_scalar))
        if fixed:
            # Stronger hostile check: t need not be free. Each attained scalar
            # has gcd(m,r-1) lifts, and exactly m/gcd(m,r-1) scalars occur.
            g = gcd(m,r-1)
            check(set(by_scalar.values()) == {g} and len(by_scalar) == m//g,
                  ('scalar distribution',m,k,b,perm,r,x,by_scalar))
        STATS['projective_lift_classes'] += 1

def mu(n):
    primes=[]
    for p in range(2,n+1):
        if n%p == 0 and all(p%q for q in range(2,p)):
            if n%(p*p) == 0:
                return 0
            primes.append(p)
    return (-1)**len(primes)

def countN(m,k):
    return sum(len(list(combinations(range(k),s)))*m**(s-1) for s in range(1,k+1))

def exact_formula(m,k,r,t,n):
    total=0
    for divisor in range(1,n+1):
        if n%divisor:
            continue
        a=countN(gcd(m,r**t*(r**divisor-1)),k)
        if t:
            a-=countN(gcd(m,r**(t-1)*(r**divisor-1)),k)
        total+=mu(n//divisor)*a
    return total

def graph_profile(points, f):
    labels={}
    for x in points:
        if x in labels:
            continue
        path=[]
        at={}
        y=x
        while y not in at and y not in labels:
            at[y]=len(path)
            path.append(y)
            y=f(y)
        if y in labels:
            t,n=labels[y]
            for z in reversed(path):
                t+=1
                labels[z]=(t,n)
        else:
            start=at[y]
            n=len(path)-start
            for z in path[start:]:
                labels[z]=(0,n)
            for i in range(start-1,-1,-1):
                labels[path[i]]=(start-i,n)
    return Counter(labels.values())

def main():
    rng=random.Random(59261025)
    # All nonunits and zero powers, with repeated and unequal cycle lengths.
    selected=[]
    for m in range(1,17):
        for k in range(1,5):
            for repeat in range(6):
                perm=list(range(k))
                rng.shuffle(perm)
                b=tuple(rng.randrange(m) for _ in range(k))
                r=[0,m,2,4,6,rng.randrange(3*m+1)][repeat]
                selected.append((m,k,b,tuple(perm),r))
    for m,k,b,perm,r in selected:
        points=pts(m,k)
        actual=sum(project(raw(x,m,b,perm,r),m)==x for x in points)
        predicted=polynomial_count(m,k,b,perm,r)
        check(actual==predicted,('semilinear',m,k,b,perm,r,actual,predicted))
        for rotation in range(1,k):
            check(polynomial_count(m,k,b,perm,r,rotation)==actual,
                  ('rotation',m,k,b,perm,r,rotation))
            STATS['cycle_start_rotations']+=1
        # Diagonal shift gauge changes no projective map.
        shifted=tuple((v+1)%m for v in b)
        check(polynomial_count(m,k,shifted,perm,r)==actual,
              ('diagonal gauge',m,k,b,perm,r))
        if m<=8 and k<=3:
            lifted_pair_audit(m,k,b,perm,r,points)
        STATS['semilinear_cases']+=1
    for m in range(1,41):
        for k in range(0,4):
            points=pts(m,k)
            for r in sorted({1,2,3,4,6,12,m,m+1,2*m}):
                profile=graph_profile(points,lambda x:tuple(-1 if a<0 else r*a%m for a in x))
                maxt=max((t for t,n in profile),default=0)
                maxn=max((n for t,n in profile),default=1)
                for t in range(maxt+2):
                    for n in range(1,maxn+2):
                        check(exact_formula(m,k,r,t,n)==profile[t,n],
                              ('tail-period',m,k,r,t,n,profile[t,n]))
                        STATS['tail_period_cells']+=1
                STATS['functional_graph_cases']+=1
    # Descent, including composite q=6 and d strictly larger than ord(sigma).
    for q in range(2,9):
        for d in range(1,5):
            m=q**d-1
            for k in range(1,4):
                if countN(m,k)>25000:
                    continue
                points=pts(m,k)
                for perm in permutations(range(k)):
                    if any(d%len(c) for c in cycles(perm)):
                        continue
                    expected=(q**k-1)//(q-1)
                    fixed=[x for x in points if project(raw(x,m,(0,)*k,perm,q),m)==x]
                    check(len(fixed)==expected,('descent',q,d,k,perm,len(fixed),expected))
                    for x in fixed:
                        # Exactly q-1 genuinely fixed lifts per projective point.
                        fixed_lifts=0
                        for shift in range(m):
                            lift=tuple(-1 if a<0 else (a+shift)%m for a in x)
                            fixed_lifts+=raw(lift,m,(0,)*k,perm,q)==lift
                        check(fixed_lifts==q-1,('descent lift',q,d,k,perm,x,fixed_lifts))
                        STATS['descent_fixed_lift_classes']+=1
                    STATS['descent_cases']+=1
    # Build coefficient morphisms as actual cyclic-group image assignments.
    for n in range(1,17):
        for m in range(1,17):
            maps=[tuple(j*a%m for j in range(n)) for a in range(m) if n*a%m==0]
            check(len(set(maps))==gcd(m,n),('coefficient map count',n,m))
            for mapping in maps:
                check(all(mapping[(i+j)%n]==(mapping[i]+mapping[j])%m
                          for i in range(n) for j in range(n)),('coefficient group law',n,m))
            STATS['coefficient_hom_cases']+=1
    # Exact correction discovered in the frozen prose. Its substantive
    # nonadditivity conclusion survives replacing 12Z by 24Z.
    gaps=(4056,3000,3696,2016)
    subgroup=0
    for value in gaps:
        subgroup=gcd(subgroup,value)
    check(subgroup==24,('D65 subgroup generator',subgroup))
    check(not [(x,y,z) for x in gaps for y in gaps for z in gaps if x+y==z],
          'unexpected D65 Schur triple')
    STATS['arithmetic_boundary_controls']+=2
    # If sigma^d=id is dropped, the descent total need not hold.
    hostile=pts(2,2)
    wrong_hypothesis_fixed=sum(project(raw(x,2,(0,0),(1,0),3),2)==x for x in hostile)
    check(wrong_hypothesis_fixed==2 and wrong_hypothesis_fixed!=(3**2-1)//(3-1),
          'descent permutation-hypothesis negative control')
    STATS['descent_negative_controls']+=1
    report={'status':'PASS','seed':59261025,'statistics':dict(sorted(STATS.items())),
            'scope':'Independent finite-exact controls; universal claims rely on audited proofs.',
            'source_packet_untouched':True,
            'frozen_packet_corrections':{'D65_generated_subgroup':'24Z, not the frozen prose 12Z',
                                         'automorphism_degeneracy':'k=0 as well as k=1'}}
    payload=json.dumps(report,indent=2,sort_keys=True)+'\n'
    (OUT/'independent_verification.json').write_text(payload)
    print(payload,end='')

if __name__=='__main__':
    main()
