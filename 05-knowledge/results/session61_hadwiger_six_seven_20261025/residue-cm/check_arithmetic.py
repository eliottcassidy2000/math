#!/usr/bin/env python3
"""Exact finite controls for the residue/CM bridge. No external dependencies.
Universe: every unit class mod55 and165; every prime 5..199; all lattice
unit directions; all residues in F4; all branch triples at2,3,7,13.
Assertions are explicit runtime checks and survive python -O.
"""
import json, math
from itertools import combinations, product

checks=0
def check(condition, message):
    global checks
    checks+=1
    if not condition: raise RuntimeError(message)
def prime(n):
    return n>=2 and all(n%d for d in range(2,math.isqrt(n)+1))
def cyc(a,m):
    out=[]; x=1
    while True:
        x=x*a%m; out.append(x)
        if x==1: return out
        check(len(out)<=m, 'cyclic subgroup terminates')
def fixed(m):
    return [(a,len(cyc(a,m)),cyc(a,m).index(m-1)+1)
            for a in range(m) if math.gcd(a,m)==1 and m-1 in cyc(a,m)]
def chi(n,p):
    r=pow(n%p,(p-1)//2,p)
    return -1 if r==p-1 else r

def norm(a,b):return a*a+a*b+b*b
units=sorted((a,b) for a in range(-3,4) for b in range(-3,4) if norm(a,b)==1)
check(units==[(-1,0),(-1,1),(0,-1),(0,1),(1,-1),(1,0)],'six lattice directions')
check(all((a+2*b)%3 for a,b in units),'three-coloring')
check(sorted((a+5*b)%7 for a,b in units)==list(range(1,7)),'norm7 quotient K7')
# F4 = F2[w]/(w²+w+1), encoded a+2b. Six units cover all three nonzero values.
check(set((a%2)+2*(b%2) for a,b in units)=={1,2,3},'mod2 quotient K4')
check(norm(2,1)==norm(3,-1)==7,'conjugate norms')
check((2+5)%7==0 and (3-5)%7!=0,'pole numerator at first prime')
check((2+3)%7!=0 and (3-3)%7==0,'pole denominator at conjugate prime')
check((1+2*0)%3!=0 and (0+2*1)%3!=0 and ((1-0)+2*(0-1))%3!=0,'source triangle three colors')

fixed55=fixed(55);fixed165=fixed(165)
check(fixed55==[(19,10,5),(24,10,5),(29,10,5),(39,10,5),(54,2,1)],'corrected mod55 classes')
check(fixed165==[(29,10,5),(74,10,5),(134,10,5),(149,10,5),(164,2,1)],'mod165 classes')
# Independent CRT screen: mod55 fixed, mod3=-1. All half-orders are odd.
crt=[a for a in range(165) if a%55 in {v[0] for v in fixed55} and a%3==2]
check(crt==[v[0] for v in fixed165],'independent CRT check')
ramified={p:cyc(p,165//p) for p in (3,5,11)}
check(all(165//p-1 not in ramified[p] for p in ramified),'no ramified fixed prime')
check(all(not prime(n) for n in range(164,659,165)) and prime(659),'least degree2 fixed prime')
check(29**5>659,'degree10 torus parameter larger')
check(660**2>4*659*12**2,'Weil real lower bound strictly greater than13')
check(660**2<4*659*13**2,'Weil real lower bound below14')
# p=109 is a hostile control: fixed mod55, not fixed mod165.
check(54 in {a for a,f,h in fixed55} and 109 not in {a for a,f,h in fixed165},'109 drops after adjoining zeta3')
# Order20 is NOT enough: both corrected near-miss cases lack conjugation.
check(len(cyc(2,55))==20 and 54 not in cyc(2,55),'order20 hostile at2')
check(len(cyc(3,55))==20 and 54 not in cyc(3,55),'ramified3 hostile')

norm_graph_triangle_counts={}
for p in (3,5,7,11,13,17,19):
    d=next(t for t in range(2,p) if chi(t,p)==-1)
    n=sum(1 for a,b in product(range(p),repeat=2) if (2*a-1)%p==0 and (a*a-d*b*b-1)%p==0)
    check(n==(1 if p==3 else 2 if p%3==2 else 0),'norm graph triangle screen')
    norm_graph_triangle_counts[p]=n

elliptic=[]
for p in range(5,200):
    if not prime(p):continue
    roots=[r for r in range(p) if (r*r-r+1)%p==0]
    check(len(roots)==(2 if p%3==1 else 0),'splitting polynomial')
    # E_-: y²=x³-1. Point counts by square enumeration and Euler character.
    count_roots=1+sum(1 for x,y in product(range(p),repeat=2) if (y*y-x*x*x+1)%p==0)
    count_char=p+1+sum(chi(x*x*x-1,p) for x in range(p))
    check(count_roots==count_char,'independent elliptic counts')
    a=p+1-count_char;m=(p-1)//2
    h=(math.comb(m,(p-1)//3)*(-1)**(m-(p-1)//3))%p if (p-1)%3==0 else 0
    check(a%p==h,'Hasse coefficient trace congruence')
    check((a==0)==(p%3==2),'supersingular trace criterion')
    # E_+ twist V²=1-U³, plus the genus-two C6 count, check degree2 quotient trace sum.
    count_plus=p+1+sum(chi(1-x*x*x,p) for x in range(p))
    count_c6=p+2+sum(chi(pow(x,6,p)-1,p) for x in range(p)) # two rational points at infinity
    check(p+1-count_c6==(p+1-count_char)+(p+1-count_plus),'genus2 quotient trace decomposition')
    elliptic.append({'p':p,'split_in_Eisenstein':bool(roots),'a_p_Eminus':a,'hasse':h,'Eminus_points':count_char,'C6_points':count_c6})

normals=[(1,0,0),(1,1,0),(1,-1,0),(1,0,1),(1,0,-1),(1,1,1),(1,-1,-1)]
def det3(m):
    a,b,c=m;return a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0])
dets={abs(det3([normals[i] for i in ids])) for ids in combinations(range(7),3)}
check(dets=={0,1,2,3},'branch incidence bad-prime support')
incidence={}
for p in (2,3,7,13):
    collapsed=len({tuple(x%p for x in n) for n in normals})
    triples=[list(ids) for ids in combinations(range(1,7),3) if det3([normals[i] for i in ids])%p==0]
    incidence[p]={'distinct_normals':collapsed,'central_dependent_triples':triples}
check(incidence[2]['distinct_normals']==4,'char2 branch collisions')
check(incidence[3]['central_dependent_triples']==[[1,3,6],[2,4,5]],'char3 new triples')
check(all(not incidence[p]['central_dependent_triples'] for p in (7,13)),'good central configuration at7,13')

out={'status':'FINITE-EXACT controls, not exhaustive all-prime proof','checks':checks,'units_omega6_convention':units,'fixed_classes55':fixed55,'fixed_classes165':fixed165,'ramified_subgroups165':ramified,'least_residue_parameter165':659,'Weil_integer_lower_bound165':14,'norm_graph_triangle_counts':norm_graph_triangle_counts,'elliptic_controls':elliptic,'branch_incidence_controls':incidence,'hostiles':['order20 does not imply conjugation fixed','109 fails after adjoining zeta3','K7 quotient source has chi3','split prime has norm-one vector pole','bad branch primes2,3']}
print(json.dumps(out,sort_keys=True,indent=2))
