#!/usr/bin/env python3
"""Complete 12-row/8-row Smith tables and independent finite-group probes."""
from itertools import product,combinations
from functools import reduce
from collections import Counter
from math import gcd,prod,lcm
from pathlib import Path
import json,hashlib
HERE=Path(__file__).parent
D12=((1,0),(0,1),(1,1),(1,-1),(1,2),(1,-2),(2,1),(2,-1),(2,0),(0,2),(2,2),(2,-2))
D8=((2,0),(0,2),(2,2),(2,-2),(2,1),(2,-1),(1,2),(1,-2))
def req(c,m):
    if not c:raise RuntimeError(m)
def invariants(rows):
    if not rows:return ()
    a=reduce(gcd,(abs(x) for row in rows for x in row),0)
    b=reduce(gcd,(abs(r[0]*s[1]-r[1]*s[0]) for r,s in combinations(rows,2)),0)
    return (a,b//a) if b else(a,)
def census(rows):
    table=Counter();subsets=[]
    for mask in range(1<<len(rows)):
        ids=[i for i in range(len(rows)) if mask>>i&1];ds=invariants([rows[i] for i in ids])
        table[ds]+=(-1)**len(ids)
        subsets.append({'mask':mask,'invariants':ds})
    return {d:c for d,c in table.items() if c},subsets

def I_formula(N,h):return N*N-4*(1+h[2])*N-2-4*h[2]+4*h[4]+2*h[5]+4*h[6]+2*h[2]**2+h[2]*h[4]
def test_group(moduli,full_centers=False):
    group=tuple(product(*(range(m) for m in moduli)));N=len(group);zero=(0,)*len(moduli)
    add=lambda a,b:tuple((x+y)%m for x,y,m in zip(a,b,moduli))
    mul=lambda c,a:tuple((c*x)%m for x,m in zip(a,moduli))
    h={d:sum(mul(d,x)==zero for x in group) for d in [1,2,3,4,5,6]}
    req(h[6]==h[2]*h[3],('coprime torsion product',moduli))
    I=0;center_count=0;witness=None;corners=0;line_counts=[0]*4;pair_counts=Counter();triple_count=0
    for x,z in product(group,repeat=2):
        offsets={add(mul(a,x),mul(b,z)) for a,b in product([-1,0,1],repeat=2)}
        distinct=len(offsets)==9
        test12=all(add(mul(a,x),mul(b,z))!=zero for a,b in D12)
        test8=all(add(mul(a,x),mul(b,z))!=zero for a,b in D8)
        req(distinct==test12==test8,('collision reduction',moduli,x,z))
        u,v=mul(2,x),mul(2,z)
        corner=all(a!=zero for a in [u,v,add(u,v),add(u,mul(-1,v))])
        corners+=corner
        if corner:
            lines=[z==mul(2,x),z==mul(-2,x),x==mul(2,z),x==mul(-2,z)]
            for i,ok in enumerate(lines):line_counts[i]+=ok
            for i,j in combinations(range(4),2):pair_counts[(i,j)]+=lines[i] and lines[j]
            triple_count+=sum(lines)>=3
        I+=distinct
        if distinct and witness is None:witness=[x,z]
        if distinct and full_centers:
            center_count+=sum(all(add(e,o)!=zero for o in offsets) for e in group)
    predicted=I_formula(N,h);req(I==predicted,('exact I',moduli,I,predicted))
    C=N*N-4*h[2]*N+2*h[2]*h[2]+h[2]*h[4]
    req(corners==C,('corner count',moduli))
    req(line_counts==[N-h[4]-h[6]+h[2]]*4,('four line sizes',moduli))
    req(triple_count==0,('no triple overlap',moduli))
    for pair,count in pair_counts.items():
        req(count==(h[5]-1 if pair in [(0,3),(1,2)] else 0),('five torsion overlap',moduli,pair))
    req(I==C-4*(N-h[4]-h[6]+h[2])+2*(h[5]-1),('human proof formula',moduli))
    if full_centers:req(center_count==(N-9)*I,('center exclusion',moduli))
    twoA=len({mul(2,x) for x in group})
    req((I>0)==(N>=9 and twoA>=4),('existence criterion',moduli))
    return {'cyclic_factors':moduli,'order':N,'torsion_sizes':h,'doubled_image_size':twoA,'injective_gap_pairs':I,
            'injective_magic_arrays':N*I,'nonzero_injective_arrays':(N-9)*I,'first_gap_witness':witness}

def all_abelian_models(bound):
    # Invariant-factor decompositions d1|...|dr with product<=bound, each di>=2.
    ans=[()]
    def rec(prefix,last,size):
        for n in range(last,bound//size+1,last):
            p=prefix+(n,);ans.append(p);rec(p,n,size*n)
    for n in range(2,bound+1):
        ans.append((n,));rec((n,),n,n)
    return sorted(set(ans),key=lambda a:(prod(a),a))

def classification_formula_models(bound):
    models=all_abelian_models(bound);bad=[]
    for factors in models:
        N=prod(factors);h={d:prod(gcd(m,d) for m in factors) for d in [2,3,4,5,6]}
        I=I_formula(N,h);twoA=prod(m//gcd(2,m) for m in factors)
        req(I>=0 and (I>0)==(N>=9 and twoA>=4),('formula classification',factors))
    return {'maximum_order':bound,'invariant_factor_models':len(models),'status':'PASS'}

def main():
    t12,s12=census(D12);t8,s8=census(D8)
    req(t12==t8,'independent 4096/256 signed tables')
    expected={():1,(1,):-4,(2,):-4,(1,1):-2,(1,2):-4,(1,4):4,(1,5):2,(1,6):4,(2,2):2,(2,4):1}
    req(t12==expected,'exact aggregate')
    # Every difference in the full 3x3 grid occurs in +/-D12.
    grid=list(product([-1,0,1],repeat=2))
    differences={(a-c,b-d) for a,b in grid for c,d in grid if (a,b)!=(c,d)}
    req(differences==set(D12)|{(-a,-b) for a,b in D12},'complete collision directions')
    cases=[(n,) for n in range(1,61)]
    cases +=[(2,2),(2,4),(4,4),(2,8),(2,2,2),(2,2,4),(3,3),(3,6),(5,5),(7,7),(2,2,2,2),(2,2,2,4),(2,2,3),(3,3,3),(2,2,2,3)]
    groups=[test_group(g,prod(g)<=25) for g in cases]
    result={'scope':'Complete finite matrix certification proves all finite abelian-group formulas; finite group probes corroborate.',
            'rows12':D12,'rows8':D8,'signed_snf':[{'invariants':d,'coefficient':c} for d,c in sorted(t12.items())],
            'subsets12':s12,'subsets8':s8,'direct_group_cases':groups,
            'formula_classification_check':classification_formula_models(512),
            'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (HERE/'distinct-certificate.json').write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    print(json.dumps({'status':'PASS','full_collision_subsets':4096,'reduced_collision_subsets':256,
                      'direct_group_cases':len(groups),'classification_models':result['formula_classification_check']['invariant_factor_models']},sort_keys=True))
if __name__=='__main__':main()
