#!/usr/bin/env python3
"""Exact exhaustive classifier for magic arrangements of f(k), k=-4,...,4.
The coefficient domain is REAL, not a bounded integer search.
f(k)=A*k^3+B*k^2+C*k+D. Rational Gaussian elimination enumerates
all real kernels of 945 exact opposite-pair systems. No external modules.
Run: python3 full_cubic_audit.py > full_cubic_audit.json
The checks use explicit exceptions and survive python -O.
"""
from collections import Counter, defaultdict
from fractions import Fraction as Q
from itertools import combinations, permutations, product
from math import gcd, lcm
import json

K=tuple(range(-4,5))

def need(condition, message):
    if not condition:
        raise RuntimeError(message)

def pairings(items):
    if not items:
        yield ()
        return
    a=items[0]
    for i in range(1,len(items)):
        b=items[i]
        for tail in pairings(items[1:i]+items[i+1:]):
            yield ((a,b),)+tail

def rref_nullspace(matrix):
    rows=[[Q(x) for x in row] for row in matrix]
    width=3
    pivots=[]
    r=0
    for c in range(width):
        choices=[i for i in range(r,len(rows)) if rows[i][c]]
        if not choices:
            continue
        i=choices[0]
        rows[r],rows[i]=rows[i],rows[r]
        divisor=rows[r][c]
        rows[r]=[x/divisor for x in rows[r]]
        for i in range(len(rows)):
            if i!=r and rows[i][c]:
                scalar=rows[i][c]
                rows[i]=[a-scalar*b for a,b in zip(rows[i],rows[r])]
        pivots.append(c)
        r+=1
    basis=[]
    for free in range(width):
        if free in pivots:
            continue
        v=[Q(0)]*width
        v[free]=Q(1)
        for i,pivot in enumerate(pivots):
            v[pivot]=-rows[i][free]
        basis.append(tuple(v))
    for v in basis:
        need(all(sum(Q(x)*y for x,y in zip(row,v))==0 for row in matrix),
             'Gaussian basis does not lie in kernel')
    return tuple(tuple(row) for row in rows[:r]),tuple(basis)

def determinant3(rows):
    (a,b,c),(d,e,f),(g,h,i)=rows
    return a*(e*i-f*h)-b*(d*i-f*g)+c*(d*h-e*g)

def rank_by_minors(rows):
    if any(determinant3(triple)!=0 for triple in combinations(rows,3)):
        return 3
    for x,y in combinations(rows,2):
        for i,j in combinations(range(3),2):
            if x[i]*y[j]-x[j]*y[i]:
                return 2
    return int(any(x for row in rows for x in row))

def primitive(v):
    denom=lcm(*(x.denominator for x in v))
    z=[int(x*denom) for x in v]
    divisor=0
    for x in z:
        divisor=gcd(divisor,abs(x))
    need(divisor!=0,'Zero vector is not a ray')
    z=[x//divisor for x in z]
    if next(x for x in z if x)!=abs(next(x for x in z if x)):
        z=[-x for x in z]
    return tuple(z)

def fraction_vectors(vs):
    return [[str(x) for x in v] for v in vs]

def centered_quartet(ray):
    A,B,C=ray
    mu=Q(20*B,3)
    values=[A*k**3+B*k*k+C*k-mu for k in K]
    counts=Counter(values)
    need(counts[0]%2==1,'The center-zero multiplicity is not odd')
    need(all(counts[x]==counts[-x] for x in counts),'Not centrally symmetric')
    quartet=sorted(x for x in values if x>0)
    quartet += [Q(0)]*((counts[0]-1)//2)
    quartet.sort()
    need(len(quartet)==4,'Incorrect quartet size')
    return quartet

def check_quartet(p):
    witnesses=[]
    for order in permutations(range(4)):
        for signs in product((-1,1),repeat=4):
            u,v,z,w=[signs[i]*p[order[i]] for i in range(4)]
            if u+v==z and u-v==w:
                witnesses.append((u,v,z,w))
    return witnesses

# Gaussian and minor-rank positive controls, including every possible rank.
for rows,expected_rank in [([(0,0,0)],0),([(1,2,3)],1),
                          ([(1,0,1),(0,1,2)],2),
                          ([(1,0,0),(0,1,0),(0,0,1)],3)]:
    rr,ns=rref_nullspace(rows)
    need(len(rr)==expected_rank==rank_by_minors(rows),'Rank control failed')
    need(len(ns)==3-expected_rank,'Nullity control failed')

systems=0
ranks=Counter()
centers=[]
nonodd_rays=set()
odd_plane_count=0
for center in K:
    matching_list=list(pairings(tuple(k for k in K if k!=center)))
    need(len(matching_list)==105 and len(set(matching_list))==105,
         'The matching universe is not exactly 105')
    spaces=Counter()
    examples={}
    local_ranks=Counter()
    for pairs in matching_list:
        need(sorted(x for pair in pairs for x in pair)==[k for k in K if k!=center],
             'Matching does not cover the eight noncentral indices')
        rows=[(3*center**3,3*center*center-20,3*center)]
        rows += [(u**3+v**3-2*center**3,
                  u*u+v*v-2*center*center,u+v-2*center) for u,v in pairs]
        rr,ns=rref_nullspace(rows)
        rank=rank_by_minors(rows)
        need(rank==len(rr) and len(ns)==3-rank,'Independent minor rank mismatch')
        systems+=1;ranks[rank]+=1;local_ranks[rank]+=1
        if not ns:
            continue
        if len(ns)==2:
            need(ns==((Q(1),Q(0),Q(0)),(Q(0),Q(0),Q(1))),
                 'Unexpected two-dimensional coefficient space')
            key=('plane_B_zero',)
            odd_plane_count+=1
        elif len(ns)==1:
            ray=primitive(ns[0]);key=('ray',)+ray
            if ray[1]:
                nonodd_rays.add(ray)
        else:
            raise RuntimeError('Unexpected full three-dimensional kernel')
        spaces[key]+=1
        examples.setdefault(key,pairs)
    descriptions=[]
    for key,count in sorted(spaces.items()):
        row={'type':key[0],'matching_count':count,
             'example_matching':examples[key]}
        if key[0]=='ray':
            row['primitive_coefficients_A_B_C']=list(key[1:])
        else:
            row['basis']=[[1,0,0],[0,0,1]]
        descriptions.append(row)
    centers.append({'center_index':center,'matrix_rank_counts':dict(sorted(local_ranks.items())),
                    'kernel_spaces':descriptions})
need(systems==945,'Incorrect total system count')
need(odd_plane_count==1,'Unexpected count for full odd-cubic plane')
expected_nonodd={(1,-3,-18),(1,3,-18),(2,-3,-39),(2,3,-39)}
need(nonodd_rays==expected_nonodd,'Nonodd symmetric rays changed')
nonodd=[]
for ray in sorted(nonodd_rays):
    p=centered_quartet(ray)
    witnesses=check_quartet(p)
    need(not witnesses,'A nonodd ray gives a magic square')
    # Independent obstruction: largest must be sum of two other members.
    need(p[-1]<p[0]+p[1],'Expected strict small-pair obstruction failed')
    nonodd.append({'primitive_coefficients_A_B_C':list(ray),
                   'positive_centered_quartet':[str(x) for x in p],
                   'magic_witness_count':len(witnesses),
                   'obstruction':'largest < sum of two smallest'})

# Classify the full B=0 plane. A=0 is linear and is a positive control.
# For A!=0 divide entries by A: f(k)=k^3+t*k, t=C/A.
odd_survivors=[];odd_unsatisfiable=0;parameter_free=0
for p in permutations(range(1,5)):
    for s in product((-1,1),repeat=4):
        a,b,c,d=p;sa,sb,sc,sd=s
        eqs=[(sa*a**3+sb*b**3-sc*c**3,sa*a+sb*b-sc*c),
             (sa*a**3-sb*b**3-sd*d**3,sa*a-sb*b-sd*d)]
        values=set();bad=False
        for const,coef in eqs:
            if coef:
                values.add(Q(-const,coef))
            elif const:
                bad=True
        if bad or len(values)>1:
            odd_unsatisfiable+=1
        elif not values:
            parameter_free+=1
        else:
            t=next(iter(values))
            odd_survivors.append({'t':str(t),'permutation':p,'signs':s})
need(odd_unsatisfiable+parameter_free+len(odd_survivors)==384,'Odd branch count failed')
need(parameter_free==0 and len(odd_survivors)==8 and
     {r['t'] for r in odd_survivors}=={'-15'},'Unexpected odd-cubic classification')

# Actual consequence objects: all eight line sums, not only a proxy.
def line_sums(matrix):
    return [sum(row) for row in matrix]+[sum(matrix[i][j] for i in range(3)) for j in range(3)]+[
        sum(matrix[i][i] for i in range(3)),sum(matrix[i][2-i] for i in range(3))]
linear=[[3,-4,1],[-2,0,2],[-1,4,-3]]
cubic=[[4,-22,18],[14,0,-14],[-18,22,-4]]
need(sorted(x for row in linear for x in row)==list(K),'Linear sample mismatch')
need(sorted(x for row in cubic for x in row)==sorted(k**3-15*k for k in K),'Cubic sample mismatch')
need(line_sums(linear)==[0]*8==line_sums(cubic),'Positive control not magic')
print(json.dumps({'status':'PASS','universe':{'coefficient_domain':'real',
    'polynomial':'A*k^3+B*k^2+C*k+D','arguments':list(K),
    'repeated_entries_allowed':True,'center_pair_systems':systems,
    'matchings_per_center':105,'odd_signed_quartet_branches':384},
    'classification':'Magic iff B=0 and (A=0 or C=-15*A); D arbitrary',
    'matrix_rank_counts':dict(sorted(ranks.items())),
    'centers':centers,'nonodd_ray_obstructions':nonodd,
    'odd_cubic_survivors':odd_survivors,
    'odd_parameter_free_branches':parameter_free,
    'positive_controls':{'linear_grid':linear,'linear_line_sums':line_sums(linear),
        'cubic_grid':cubic,'cubic_line_sums':line_sums(cubic)},
    'scope_boundary':'These are numerical label assignments, not square-valued assignments; general magic squares of squares remain open.'},indent=2))
