#!/usr/bin/env python3
"""Independent audit by enumerating every rank-three subgroup of C2^6."""
from itertools import combinations
from collections import Counter
from math import gcd
import json
from pathlib import Path

def need(c,m):
    if not c: raise RuntimeError(m)
def span(vs):
    s={0}
    for v in vs:s|={x^v for x in s}
    return frozenset(s)
def cross(a,b):return (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])
def dot(a,b):return sum(x*y for x,y in zip(a,b))
def proj(v):
    g=0
    for x in v:g=gcd(g,abs(x))
    v=tuple(x//g for x in v)
    return v if next(x for x in v if x)!=-abs(next(x for x in v if x)) else tuple(-x for x in v)
def matchings(items):
    if not items:yield ();return
    a=items[0]
    for b in items[1:]:
        for tail in matchings([x for x in items if x not in (a,b)]):yield tuple(sorted(((a,b),)+tail))
def even(k):return k if k.bit_count()%2==0 else k^127

rows=[(1,0,0),(1,1,0),(1,-1,0),(1,0,1),(1,0,-1),(1,1,1),(1,-1,-1)]
points={proj(cross(rows[i],rows[j])) for i,j in combinations(range(7),2)}
zero_sets={tuple(i for i,r in enumerate(rows) if dot(r,p)==0):p for p in points}
triples=sorted(s for s in zero_sets if len(s)==3)
need(triples==[(0,1,2),(0,3,4),(0,5,6)],'exact triple centers')
need(Counter(map(len,zero_sets))=={2:12,3:3},'intersection census')
meridians=[1<<i for i in range(6)]+[63]
# This route does not enumerate line label permutations or weight-four generators.
all_subgroups={span(b) for b in combinations(range(1,64),3) if len(span(b))==8}
need(len(all_subgroups)==1395,'all rank-three subgroups')
forbidden=set(meridians)|{a^b for a,b in combinations(meridians,2)}
valid=sorted((tuple(sorted(k)) for k in all_subgroups if not(k&forbidden)))
need(len(valid)==30,'distinct nonzero inertia labels')
need(len(forbidden)==28,'seven inertia and twenty-one collision exclusions')
records=[]
for kt in valid:
    k=set(kt)
    even_code={even(x) for x in k}
    need(Counter(x.bit_count() for x in even_code)=={0:1,4:7},'simplex even representatives')
    fano={tuple(i for i in range(7) if not(c>>i&1)) for c in even_code if c}
    need(len(fano)==7 and all(sum(i in t and j in t for t in fano)==1 for i,j in combinations(range(7),2)), 'Fano pair axiom')
    fixed=[]
    for g in k-{0}:
        pos={i for i in range(7) if g>>i&1}
        neg=set(range(7))-pos
        places=[]
        for zs,p in zero_sets.items():
            if pos<=set(zs) or neg<=set(zs):places.append(zs)
        if places:
            need(len(places)==1 and len(places[0])==3,'kernel fixed locus is only one triple fiber')
            fixed.append(places[0])
    dependent=[t for t in triples if (meridians[t[0]]^meridians[t[1]]^meridians[t[2]]) in k]
    need(sorted(fixed)==dependent==sorted(set(triples)&fano),'three independent fixed/dependency tests')
    d=len(dependent)
    # Universal triple fiber is a D/local-inertia coset set of size 8.
    node_orbits=[]
    for t in triples:
        local=span([meridians[i] for i in t])
        need(len(local)==8,'independent universal inertia')
        stabilizer=k&local
        orbit_size=8//len(stabilizer)
        number_of_orbits=8//orbit_size
        need((number_of_orbits,len(stabilizer))==((2,2) if t in dependent else (1,1)), 'node quotient orbit count')
        node_orbits.append({'triple':t,'quotient_points':number_of_orbits,'stabilizer_order':len(stabilizer)})
    records.append({'kernel6':kt,'even_code':sorted(even_code),'d':d,'fixed_triples':fixed,'node_orbits':node_orbits})
need(Counter(r['d'] for r in records)=={0:16,1:12,3:2},'kernel distribution')
ms=list(matchings(list(range(1,7))))
actual={(1,2),(3,4),(5,6)}
matching_hist=Counter(len(set(m)&actual) for m in ms)
need(len(ms)==15 and matching_hist=={0:8,1:6,3:1},'independent matching distribution')
# Direct intersections on Bl_3(P2): (h,e1,e2,e3), Ei^2=-1.
def meet(a,b):return a[0]*b[0]-sum(x*y for x,y in zip(a[1:],b[1:]))
KY=(-3,1,1,1)
eigen_checks=0
for subset in range(8):
    negative_L=(-2,)+tuple((subset>>i)&1 for i in range(3))
    negative_L_minus_K=tuple(a-b for a,b in zip(negative_L,KY))
    need(1+meet(negative_L,negative_L_minus_K)//2==0,'all possible eigensheaf Euler characteristics vanish')
    need(negative_L[0]<0 and KY[0]-negative_L[0]<0,'H0 and dual H0 excluded by nef H')
    eigen_checks+=1
# Independent action of the realized twelve-element arrangement group.
rotation=(0,5,6,2,1,3,4)
reflection=(0,3,4,1,2,5,6)
def compose(p,q):return tuple(p[q[i]] for i in range(7))
perms={tuple(range(7))}
while True:
    enlarged=perms|{compose(g,h) for g in perms for h in (rotation,reflection)}
    if enlarged==perms:break
    perms=enlarged
need(len(perms)==12,'realized arrangement group order')
def permute_mask(v,p):return sum(((v>>i)&1)<<p[i] for i in range(7))
def act(code,p):return tuple(sorted(permute_mask(v,p) for v in code))
bycode={tuple(r['even_code']):r for r in records}
remaining=set(bycode)
orbits=[]
while remaining:
    code=min(remaining)
    orbit={act(code,p) for p in perms}
    need(orbit<=set(bycode),'closed action on thirty kernels')
    remaining-=orbit
    stab=sorted(p for p in perms if act(code,p)==code)
    d=bycode[code]['d']
    need({bycode[k]['d'] for k in orbit}=={d},'dependency count invariant')
    fixed_quotient_counts=[]
    if d==0 and len(stab)==3:
        # Even signs are D; cosets modulo the code give the rank-three deck quotient.
        evens=[v for v in range(128) if v.bit_count()%2==0]
        cosets={min(v^c for c in code) for v in evens}
        need(len(cosets)==8,'quotient deck order')
        for p in stab:
            if p==tuple(range(7)):continue
            number=sum(min(permute_mask(v,p)^c for c in code)==v for v in cosets)
            need(number==2,'order-three deck action fixes one binary dimension')
            fixed_quotient_counts.append(number)
    orbits.append({'d':d,'size':len(orbit),'stabilizer_order':len(stab),'stabilizer':stab,
                   'nonidentity_C3_fixed_deck_elements':fixed_quotient_counts})
need(sorted((o['d'],o['size'],o['stabilizer_order']) for o in orbits)==[(0,4,3),(0,12,1),(1,6,2),(1,6,2),(3,2,6)],'D12 orbit census')
out={'status':'PASS','route':'all 1395 rank-three subgroups of C2^6, not 5040 line-label permutations','intersection_histogram':dict(Counter(map(len,zero_sets))),
     'triple_centers':[{'indices':t,'point':zero_sets[t]} for t in triples], 'all_rank3_subgroups':1395,'valid_Fano_kernels':30,'dependency_histogram':dict(Counter(r['d'] for r in records)),
     'matching_histogram':dict(matching_hist),'eigendisor_RR_cases':eigen_checks,'records':records,'D12_orbits':orbits}
Path(__file__).with_name('independent_certificate.json').write_text(json.dumps(out,sort_keys=True,indent=2)+'\n')
print('PASS: all 1395 rank-three subgroups; 30 Fano kernels; d=0/1/3 has 16/12/2.')
print('PASS: exact line intersections, all kernel fixed loci, triple-fiber quotient orbits.')
print('PASS: independent 15-matchings census; 8 eigensheaf Riemann-Roch possibilities.')
print('PASS: independent D12 orbits; nontrivial C3 action on rank-three quotient deck.')
