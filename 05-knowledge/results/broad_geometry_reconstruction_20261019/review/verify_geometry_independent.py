#!/usr/bin/env python3
"""Independent integer controls. No source verifier is imported.
Theorems are proved in the review; finite checks do not verify analytic realizability.
"""
from collections import Counter
from itertools import combinations, permutations
from fractions import Fraction
from math import gcd, prod
from functools import reduce
import json


def check(condition, message):
    if not condition:
        raise RuntimeError(message)


def factor(n):
    result = {}
    q = 2
    while q*q <= n:
        while n % q == 0:
            result[q] = result.get(q, 0) + 1
            n //= q
        q += 1
    if n > 1:
        result[n] = 1
    return result


def squarefree_divisors(n):
    result = [(1, 1)]
    for q in factor(n):
        result += [(q*d, -mu) for d, mu in list(result)]
    return result


def valuation(p, q):
    if not p:
        return float('inf')
    p = abs(p)
    e = 0
    while p % q == 0:
        e += 1
        p //= q
    return e


def local_count(a, b, p):
    answer = 1
    for q in factor(a*b):
        aa, bb, pp = valuation(a,q), valuation(b,q), valuation(p,q)
        e = min(aa,bb)
        if aa != bb:
            if pp != e:
                return 0
            answer *= 1 if e == 0 else (q-1)*q**(e-1)
        elif pp < e:
            return 0
        else:
            answer *= (q-1 if pp > e else q-2)*q**(e-1)
    return answer


def mobius_terms(a,b):
    return [(d*e*gcd(a//d,b//e), mu*nu*gcd(a//d,b//e))
            for d,mu in squarefree_divisors(a)
            for e,nu in squarefree_divisors(b)]


def normalize(v,a,b):
    z,x,y=v
    q,x=divmod(x,a)
    r,y=divmod(y,b)
    return (z-q-r,x,y)


def gauge_checks():
    residues = signed = torsion_points = 0
    for a in range(2,25):
        for b in range(2,25):
            g = gcd(a,b)
            all_pairs = Counter((-b*x-a*y)%(a*b) for x in range(a) for y in range(b))
            units = Counter((-b*x-a*y)%(a*b) for x in range(a) for y in range(b)
                            if gcd(x,a) == gcd(y,b) == 1)
            terms = mobius_terms(a,b)
            for p in range(-a*b,a*b):
                expected = sum(w for modulus,w in terms if p % modulus == 0)
                check(expected == units[p%(a*b)], f'Mobius mismatch {(a,b,p)}')
                check(expected == local_count(a,b,p), f'local mismatch {(a,b,p)}')
                check(all_pairs[p%(a*b)] == (g if p%g == 0 else 0), 'full fibre count')
                signed += 1
            residues += a*b
            orbit = {normalize((0,r*(a//g),-r*(b//g)),a,b) for r in range(g)}
            check(len(orbit) == g, 'torsion orbit is not exact')
            check(normalize((0,g*(a//g),-g*(b//g)),a,b) == (0,0,0), 'torsion period')
            for z,x,y in orbit:
                check(a*b*z-b*x-a*y == 0, 'torsion not in clutch kernel')
            torsion_points += len(orbit)
    hostile = [(1,1,2),(1,2,1),(2,4,4)]
    check({normalize(v,5,5) for v in hostile} == set(hostile), 'hostile gauge collision')
    for z,x,y in hostile:
        check(25*z-5*x-5*y == 10,'hostile clutch')
        check(gcd(gcd(5,x),gcd(5,y-5*z)) == 1,'clutch cyclicity')
    return {'ordered_a_b_range':[2,24], 'ordered_pairs':23**2,
            'residue_classes':residues, 'signed_integer_tests':signed,
            'torsion_orbit_points':torsion_points,
            'independent_methods':['normalized residue enumeration','Mobius inclusion-exclusion','primewise formula'],
            'hostile_a_b_p':[5,5,10], 'hostile_representatives':hostile}


def determinant(a):
    n=len(a)
    return sum((-1)**sum(p[i]>p[j] for i in range(n) for j in range(i+1,n))
               *prod(a[i][p[i]] for i in range(n)) for p in permutations(range(n)))


def eye(n):
    return [[int(i==j) for j in range(n)] for i in range(n)]


def multiply(a,b):
    return [[sum(x*y for x,y in zip(row,col)) for col in zip(*b)] for row in a]


def inverse(a):
    n=len(a)
    aug=[[Fraction(x) for x in row]+[Fraction(i==j) for j in range(n)]
         for i,row in enumerate(a)]
    for j in range(n):
        p=next(i for i in range(j,n) if aug[i][j])
        aug[p],aug[j]=aug[j],aug[p]
        pivot=aug[j][j]
        aug[j]=[x/pivot for x in aug[j]]
        for i in range(n):
            if i!=j:
                c=aug[i][j]
                aug[i]=[x-c*y for x,y in zip(aug[i],aug[j])]
    check(all(x.denominator==1 for row in aug for x in row[n:]),'nonintegral inverse')
    return [[int(x) for x in row[n:]] for row in aug]


def rank(a):
    a=[[Fraction(x) for x in row] for row in a]
    rr=0
    for j in range(len(a[0])):
        p=next((i for i in range(rr,len(a)) if a[i][j]),None)
        if p is None:
            continue
        a[p],a[rr]=a[rr],a[p]
        pv=a[rr][j]
        a[rr]=[x/pv for x in a[rr]]
        for i in range(rr+1,len(a)):
            c=a[i][j]
            a[i]=[x-c*y for x,y in zip(a[i],a[rr])]
        rr+=1
        if rr==len(a):
            break
    return rr


def exterior(a,q):
    ids=list(combinations(range(len(a)),q))
    return [[determinant([[a[i][j] for j in jj] for i in ii]) for jj in ids] for ii in ids]


def determinantal_profile(a):
    r=rank(a)
    ds=[1]
    for k in range(1,r+1):
        d=0
        for ii in combinations(range(len(a)),k):
            for jj in combinations(range(len(a[0])),k):
                d=gcd(d,determinant([[a[i][j] for j in jj] for i in ii]))
        ds.append(d)
    return [ds[i+1]//ds[i] for i in range(r)]


def exterior_checks():
    matrices=0
    conjugate_cases=0
    examples=[]
    for s in range(1,8):
        for m in range(1,9):
            t=s*m
            shear=eye(4); shear[0][2]=s; shear[1][3]=t
            # Full lattice mixing, including the kernel and complementary coordinates.
            change=eye(4)
            for k in range(1,10):
                i=(s+k)%4; j=(i+1+(m+k)%3)%4
                change[i]=[x+((-1)**k)*(k%3+1)*y for x,y in zip(change[i],change[j])]
            check(abs(determinant(change))==1,'change not unimodular')
            monodromy=multiply(multiply(change,shear),inverse(change))
            conjugate_cases+=1
            profiles=[]
            for q in range(5):
                w=exterior(monodromy,q)
                difference=[[x-int(i==j) for j,x in enumerate(row)] for i,row in enumerate(w)]
                profile=determinantal_profile(difference)
                expected=[[],[s,t],[s,s],[s,t],[]][q]
                check(profile==expected,f'exterior profile {(s,t,q)}: {profile}')
                profiles.append(profile)
                matrices+=1
            if s==1 and m in (1,2,7):
                examples.append({'s':s,'t':t,'profiles':profiles})
    free=[1,3,6,6,3,1]
    check(sum((-1)**k*b for k,b in enumerate(free))==0,'mapping torus Euler')
    return {'smith_method':'gcds of minors after rational rank; no Smith-library call',
            'full_unimodular_conjugate_cases':conjugate_cases,
            'exterior_matrices_checked':matrices,'range_s':[1,7],'range_t_over_s':[1,8],
            'mapping_torus_betti_numbers':free,'examples':examples}


def toric_combinatorial_checks():
    rows=[]
    for B in ([[1,0],[0,d]] for d in (1,2,3,5,12,31)):
        d=abs(determinant(B))
        adj=[[B[1][1],-B[0][1]],[-B[1][0],B[0][0]]]
        classes={tuple(sum(adj[i][j]*v[j] for j in range(2))%d for i in range(2))
                 for v in ((x,y) for x in range(d) for y in range(d))}
        check(len(classes)==d,'lattice quotient count')
        # The A2 triangulation has one vertex, three edge and two triangle
        # types per lattice translate. This checks the finite orbit census only.
        counts=[len(classes),3*len(classes),2*len(classes)]
        check(counts[0]-counts[1]+counts[2]==0,'real torus Euler')
        rows.append({'B':B,'lattice_index':d,'A2_cell_orbits':counts,
                     'conditional_toric_Euler':counts[2]})
    return {'scope':'A2 cell orbits; no analytic quotient or specialization is constructed',
            'controls':rows}


if __name__=='__main__':
    result={'gauge':gauge_checks(),'exterior':exterior_checks(),
            'toric_combinatorics':toric_combinatorial_checks(),
            'result':'PASS'}
    print(json.dumps(result,indent=2,sort_keys=True))
