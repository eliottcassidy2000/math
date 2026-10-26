#!/usr/bin/env python3
"""Independent whole-distribution involution audit plus weighted divisor checks."""
from collections import Counter
from itertools import combinations, product
from math import comb, isqrt, prod
from pathlib import Path
import hashlib
import json
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
checks = 0
def require(ok, detail):
    global checks
    checks += 1
    if not ok:
        raise RuntimeError(detail)

def ds(n):
    return [d for d in range(1,n+1) if n%d == 0]

def factor(n):
    result=[]
    for p in range(2,n+1):
        if p*p>n:
            break
        e=0
        while n%p==0:
            n//=p
            e+=1
        if e:
            result.append((p,e))
    if n>1:
        result.append((n,1))
    return result

distributions=[]
orientation_count=0
for f in range(7):
    for b in range((6-f)//2+1):
        n=f+2*b
        perm=list(range(n))
        for i in range(f,n,2):
            perm[i],perm[i+1]=perm[i+1],perm[i]
        edges=list(combinations(range(n),2))
        where={edge:i for i,edge in enumerate(edges)}
        image=[]
        for u,v in edges:
            x,y=perm[u],perm[v]
            image.append((where[tuple(sorted((x,y)))], int(x>y)))
        hist=Counter()
        for mask in range(1<<len(edges)):
            defect=sum(((mask>>i)&1)^((mask>>j)&1)^flip for i,(j,flip) in enumerate(image))
            hist[defect]+=1
            orientation_count+=1
        M=b*(f+b-1)
        expected={b+2*j:2**(comb(f,2)+b+M)*comb(M,j) for j in range(M+1)}
        require(dict(hist)==expected, ('whole polynomial',f,b,dict(hist),expected))
        require(sum(hist.values())==1<<comb(n,2), ('total orientations',f,b))
        distributions.append({'f':f,'b':b,'coefficients':dict(sorted(hist.items()))})

for N in range(1,601):
    divisors=ds(N)
    factors=factor(N)
    tuples=list(product(*(range(e+1) for p,e in factors)))
    weights=[prod(p**j for (p,e),j in zip(factors,exps)) for exps in tuples]
    require(sorted(weights)==divisors, ('polynomial evaluation',N))
    require(sum(weights)-N==sum(divisors)-N, ('aliquot',N))
    require(len(weights)-1==len(divisors)-1, ('augmentation',N))
    if N>1:
        interior=divisors[1:-1]
        require(set(N//d for d in interior)==set(interior), ('complement carrier',N))
        f=int(isqrt(N)**2==N)
        require(sum(N//d==d for d in interior)==f, ('fixedpoints',N))
        require(all(N//a>N//b for a,b in combinations(interior,2)), ('numeric reversal',N))

for n in range(1,101):
    s=sum(ds(n))-n
    for k in range(2,31):
        actual=sum(ds(k*n))-k*n
        require(actual>=s+n, ('divisor injection',n,k))
        require(actual>=k*s+1, ('scaled injection',n,k))

# Execute frozen producer in both modes without touching its source directory.
path=ROOT/'source-snapshot/divisor_carrier.py'
normal=subprocess.check_output([sys.executable,'-I','-B',str(path)])
optimized=subprocess.check_output([sys.executable,'-I','-B','-O',str(path)])
require(normal==optimized, 'root normal/optimized mismatch')
producer=json.loads(normal)
arcs={tuple(e) for e in producer['N210']['minimum_orientation_arcs']}
interior=ds(210)[1:-1]
require(len(arcs)==comb(14,2), '210 orientation cardinality')
require(all(((u,v) in arcs)+((v,u) in arcs)==1 for u,v in combinations(interior,2)), '210 exactly one orientation')
require(sum((210//u,210//v) not in arcs for u,v in arcs)==7, '210 actual minimizing defect')
require(sum(ds(210))-210==366, 's210')
require({d for d in ds(210)[:-1]}!={210//d for d in ds(210)[:-1]}, 'proper carrier hostile')
require(sum(ds(4))-4==3 and sum(ds(9))-9==4, 'prime label hostile')
require(2**49*2**42==2**91, '210 polynomial normalization')

print(json.dumps({
    'status':'PASS', 'explicit_checks':checks,
    'exhaustive_orientations':orientation_count,
    'involution_distributions':distributions,
    'carrier_range':[1,600], 'strict_divisibility_n':[1,100], 'strict_divisibility_k':[2,30],
    'N210':{'proper_count':15,'interior_count':14,'minimum_defect':7,'numeric_defect':91,'aliquot':366},
    'producer_normal_optimized_equal':True,
    'frozen_producer_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
    'root_replay_output_sha256':hashlib.sha256(normal).hexdigest()
},sort_keys=True,indent=2))
