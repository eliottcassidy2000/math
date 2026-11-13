"""Exact controls for the all-order half-repair growth lemma.

Run from the repository root with Python normally and with -O.
No assert statements: failed controls remain active in optimized mode.
"""
from pathlib import Path
from itertools import permutations
import gzip,json,sys
import numpy as np

HERE=Path(__file__).resolve().parent
PRIOR=HERE.parent/'endpoint_context_recursion_20261113'
sys.path.insert(0,str(PRIOR))
from port_compiler import compiler,mul
from typed_glue import seven,glue
from feedback_probe import all_profiles

D=(0,1,2,1)
SEED=((2,4,9,3),(2,8,3,5))

def need(ok,msg):
    if not ok:raise RuntimeError(msg)

def join(*blocks):
    out=[];offset=0;N=sum(map(len,blocks))
    for B in blocks:
        after=((1<<N)-1)^((1<<(offset+len(B)))-1)
        out.extend((row<<offset)|after for row in B)
        offset+=len(B)
    return tuple(out)

def tt(n):return tuple(((1<<n)-1)^((1<<(i+1))-1) for i in range(n))

def is_return(d):return d[0]==0 and d[1]!=0 and d[2]==2*d[1] and d[3]==d[1]

def brute(T,r):
    n=len(T);h=s=q=F=0
    for order in permutations(range(n)):
        good=all(T[a]>>b&1 for a,b in zip(order,order[1:]))
        h+=good;s+=good and order[0]==r
        cut=order.index(r)+1
        F+=all(T[a]>>b&1 for j,(a,b) in enumerate(zip(order,order[1:])) if j+1!=cut)
    other=[v for v in range(n) if v!=r]
    for order in permutations(other):q+=all(T[a]>>b&1 for a,b in zip(order,order[1:]))
    return h,s,q,F

def run():
    records=[];physical=compiled=literal=0
    for n in range(2,11):
        contexts=[('transitive',tt(n))]
        if n>=3:
            contexts.extend((f'C3:{a}:{n-3-a}',join(tt(a),(2,4,1),tt(n-3-a))) for a in range(n-2))
        for label,A in contexts:
            profiles=all_profiles(np.array(A,dtype=np.int64))
            for r,raw in enumerate(profiles):
                h,s,e,q,F,G,P=map(int,raw)
                if G+e!=2*h:continue
                Ts=[glue(A,r,B,3) for B in SEED]
                Z=[all_profiles(np.array(T,dtype=np.int64)) for T in Ts]
                images=[]
                for u in range(n):
                    # Seed[0]-Seed[1] = D.
                    d=tuple(int(Z[0][u,j]-Z[1][u,j]) for j in (0,1,3,4))
                    need(not is_return(d),'low-H context retained a nonzero half repair')
                    images.append(d);physical+=1
                    if n<=5:
                        M=compiler(A,r,u)
                        need(mul(M,D)==d,'abstract compiler versus actual pair')
                        compiled+=1
                    if n<=3:
                        zs=[brute(T,u) for T in Ts]
                        need(tuple(x-y for x,y in zip(*zs))==d,'permutation versus subset DP')
                        literal+=1
                records.append(dict(label=label,rows=A,root=r,h=h,e=e,images=images))
    # Every catalogued half-return must have the claimed universal target row.
    catalog=json.loads(gzip.decompress((PRIOR/'renewal_catalog.json.gz').read_bytes()))
    catalog_checks=[]
    for row in catalog['records']:
        h,s,e,q,F,G,P=seven(tuple(row['rows']),row['entry'])
        need(G+e==2*h and h>=5,'catalog neutrality and lower bound')
        need(tuple(row['target_matrix'][0])==(e,h-e),'universal target row')
        catalog_checks.append(dict(order=row['n'],h=h,e=e,target_row=(e,h-e)))
    special=(10,12,9,0)
    images=[mul(compiler(special,0,u),D) for u in range(4)]
    U=compiler((2,4,9,3),3,2)
    need(mul(U,D)==(0,-1,-2,-1),'positive U return control')
    result=dict(status='FINITE-EXACT controls; all-order proof in GROWTH.md',universe='All H=1 and H=3 SCC templates of orders2..10, every entry satisfying G+e=2h, every context output vertex; compiler crosscheck through5, literal permutation crosscheck through3',physical=physical,compiled=compiled,literal=literal,records=records,catalog_checks=catalog_checks,special=dict(rows=special,entry=0,images=images),positive_U=dict(rows=(2,4,9,3),entry=3,exit=2,matrix=U))
    suffix='' if __debug__ else '.optimized'
    (HERE/('growth'+suffix+'.json')).write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k not in ('records','positive_U')}),flush=True)

if __name__=='__main__':run()
