#!/usr/bin/env python3
"""Independent proof-review controls; no producer imports or removable asserts.

Direct projective root tuples at q=3,5,7 are a third count path. Cremona
polynomial arithmetic is an independent integer-dictionary implementation.
The normal/optimized outputs are deterministic and intentionally concise.
"""
from collections import defaultdict
from fractions import Fraction
from itertools import combinations, product
from pathlib import Path
import hashlib
import json

HERE=Path(__file__).resolve().parent

def need(test,message):
    if not test:
        raise RuntimeError(message)

def pp(q,n):
    for pivot in range(n):
        for tail in product(range(q),repeat=n-pivot-1):
            yield (0,)*pivot+(1,)+tail

def roots_valid(v,q):
    sq=[x*x%q for x in v]
    e=sq[0]
    return ((sq[1]+sq[2]-2*e)%q==0
            and (sq[3]+sq[4]-2*e)%q==0
            and (sq[5]+sq[6]-2*e)%q==0
            and (sq[5]-sq[1]-sq[3]+e)%q==0)

def char(a,q):
    a%=q
    if not a:
        return 0
    return 1 if pow(a,(q-1)//2,q)==1 else -1

root_rows=[]
for q in (3,5,7):
    direct=sum(roots_valid(v,q) for v in pp(q,7))
    expanded=0
    fibre_counts=[]
    for e,x,z in pp(q,3):
        forms=(e,e+x,e-x,e+z,e-z,e+x+z,e-x-z)
        terms=0
        for mask in range(128):
            if mask.bit_count()%2:
                continue
            term=1
            for j in range(7):
                if mask>>j&1:
                    term*=forms[j]
            terms+=char(term,q)
        # Explicit scalar/root multiplicity average, no character expansion.
        raw=0
        for scalar in range(1,q):
            ways=1
            for form in forms:
                ways*=sum(u*u%q==scalar*form%q for u in range(q))
            raw+=ways
        need(raw==(q-1)*terms,'local scalar/root fibre disagreement')
        expanded+=terms
        fibre_counts.append(terms)
    need(direct==expanded,'direct P6/character total disagreement')
    root_rows.append(dict(q=q,projective_ambient_points=sum(1 for _ in pp(q,7)),
                          root_points=direct,character_points=expanded,
                          fibre_values=sorted(set(fibre_counts))))

# Integer polynomial dictionary, disjoint from the producer's SymPy path.
lines=((-0,-2,-1),(2,0,-1),(0,1,-2),(2,-3,0),(-1,0,-2),(3,-2,0))
def multiply(p,q):
    out=defaultdict(int)
    for a,c in p.items():
        for b,d in q.items():
            out[tuple(x+y for x,y in zip(a,b))]+=c*d
    return {m:c for m,c in out.items() if c}

def linear(v):
    return {tuple(int(i==j) for i in range(3)):c for j,c in enumerate(v) if c}

F={(0,0,0):1}
for line in lines:
    F=multiply(F,linear(line))
pull=defaultdict(int)
for (a,b,c),value in F.items():
    pull[(b+c,a+c,a+b)]+=(-1)**c*value
pull={m:c for m,c in pull.items() if c}
expected={tuple(a+2 for a in m):c for m,c in F.items()}
need(pull==expected,'independent sextic identity failed')

def cross(a,b):
    return (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])
def norm(v):
    lead=next(x for x in v if x)
    return tuple(Fraction(x,lead) for x in v)
def dot(a,b):
    return sum(x*y for x,y in zip(a,b))

nodes={norm(cross(lines[i],lines[j])) for i,j in combinations(range(6),2)}
need(len(nodes)==15,'wrong line crossing count')
need(all(sum(dot(line,p)==0 for line in lines)==2 for p in nodes),'unexpected triple')
centers={(Fraction(1),Fraction(0),Fraction(0)),
         (Fraction(0),Fraction(1),Fraction(0)),
         (Fraction(0),Fraction(0),Fraction(1))}
remaining=nodes-centers
need(len(remaining)==12,'wrong unresolved node count')
need(all(all(coordinate for coordinate in p) for p in remaining),'node on boundary')
orbits=set()
for p in remaining:
    a,b,c=p
    image=norm((b*c,a*c,-a*b))
    need(image in remaining and image!=p,'bad node involution')
    orbits.add(tuple(sorted((p,image))))
need(len(orbits)==6,'wrong node orbit count')

# At F17 every geometric fixed fiber splits, giving a finite-field control
# of the 8/0 assertion without pretending it proves the complex result.
q=17
fixed=[]
for u,v in product(range(1,q),repeat=2):
    if (u*u-1)%q or (v*v+1)%q:
        continue
    value=1
    for row in lines:
        value=value*dot(row,(1,u,v))%q
    value=value*pow((u*v)**2,-1,q)%q
    sheets=[a for a in range(q) if a*a%q==value]
    need(value!=0 and len(sheets)==2,'fixed fiber not split/etale at F17')
    fixed.append(dict(u=u,v=v,branch_value=value,sheets=sheets))
need(len(fixed)==4,'wrong base fixed locus')
need(sum(len(row['sheets']) for row in fixed)==8,'wrong plus fixed count')
need(not any(a==(-a)%q for row in fixed for a in row['sheets']),
     'minus lift unexpectedly fixes a sheet')

result=dict(status='PASS / FINITE-EXACT independent controls',
            direct_projective_root_enumeration=root_rows,
            cremona=dict(polynomial_terms=len(F),integer_identity=True,
                         branch_crossings=len(nodes),remaining_node_pairs=len(orbits),
                         fixed_fibres_F17=fixed),
            source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
encoded=json.dumps(result,indent=2,sort_keys=True)+'\n'
(HERE/'character_support_independent.json').write_text(encoded,encoding='utf-8')
print(encoded,end='')
