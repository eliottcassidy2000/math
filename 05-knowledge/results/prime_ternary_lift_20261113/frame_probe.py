"""All two-run transitive-module frames through eight vertices.

Two-run means the contracted vertex's out-to-in bipartite graph has no
matching of size two. This is an exact structural filter, not interpolation.
"""
from pathlib import Path
from collections import Counter
import gzip,json,sys
import numpy as np

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'endpoint_context_recursion_20261113'))
from context_head import ordinary_bank
from feedback_probe import all_profiles
from inflation_probe import end_table

def need(ok,msg):
    if not ok:raise RuntimeError(msg)

def clone(T,r,m):
    old=[v for v in range(len(T)) if v!=r];labels=old+[r]*m
    return tuple(sum(int(i!=j and ((a==b and i<j) or (a!=b and T[a]>>b&1)))<<j for j,b in enumerate(labels)) for i,a in enumerate(labels))

def hp(T):return int(end_table(np.array(T,dtype=np.int64))[-1].sum())

def valid(T):
    n=len(T)
    return all(not(T[i]>>i&1) and T[i]<(1<<n) for i in range(n)) and all(((T[i]>>j)&1)+((T[j]>>i)&1)==1 for i in range(n) for j in range(i+1,n))

def two_run(T,r):
    I=[v for v in range(len(T)) if T[v]>>r&1];O=[v for v in range(len(T)) if T[r]>>v&1]
    arcs=[(o,i) for o in O for i in I if T[o]>>i&1]
    return not any(a!=c and b!=d for a,b in arcs for c,d in arcs)

def strong(T):
    n=len(T);full=(1<<n)-1
    for r in range(n):
        seen=1<<r
        while True:
            nxt=seen
            for v in range(n):
                if seen>>v&1:nxt|=T[v]
            if nxt==seen:break
            seen=nxt
        if seen!=full:return False
    return True

def valuation2(a):
    k=0
    while a%2==0:a//=2;k+=1
    return k

def run():
    catalog={};occurrences=Counter();controls=0
    for (n,code),T in sorted(ordinary_bank().items()):
        if n<3 or not strong(T):continue
        h=hp(T)
        for r in range(n):
            if not two_run(T,r):continue
            h2=hp(clone(T,r,2));alpha=(h2-h)//2;beta=h-2*alpha
            need(h2-h==2*alpha,'even twin difference')
            if alpha==0:continue
            key=(alpha,beta);occurrences[key]+=1
            if key not in catalog:
                need(all(valid(clone(T,r,m)) for m in (1,2,3,4)),'valid tournament after cloning, including internal clone arcs')
                vals=[hp(clone(T,r,m)) for m in (1,2,3,4)]
                need(vals==[alpha*2**m+beta for m in (1,2,3,4)],'structurally guarded affine response')
                controls+=4
                catalog[key]=dict(n=n,rows=T,root=r,alpha=alpha,beta=beta,checks=vals)
    rays={}
    for (a,b),rec in catalog.items():
        v=valuation2(a);key=(a>>v,b)
        old=rays.get(key)
        cand=dict(**rec,shift=v,min_exponent=v+1)
        if old is None or (cand['min_exponent'],cand['n'])<(old['min_exponent'],old['n']):rays[key]=cand
    transitions=[]
    for (a,b),source in sorted(rays.items()):
        for digit in (2,4):
            target=rays.get((3*a,3*b+digit))
            if target:
                transitions.append(dict(odd_amplitude=a,offset=b,digit=digit,min_exponent=max(source['min_exponent'],target['min_exponent']),source=source,target=target))
    result=dict(status='FINITE-EXACT catalog; all-exponent identities use the two-run theorem',scope='All strong ordinary frames3..8, every module vertex with matching number<=1; retain one least-order witness per affine response, then normalize powers of2 in the amplitude. Only exact frame replacements with the same normalized exponent are listed.',profile_count=len(catalog),ray_count=len(rays),occurrences=sum(occurrences.values()),controls=controls,profiles=[dict(**r,occurrences=occurrences[k]) for k,r in sorted(catalog.items())],rays=list(rays.values()),transitions=transitions)
    suffix='' if __debug__ else '.optimized'
    (HERE/('frames'+suffix+'.json.gz')).write_bytes(gzip.compress(json.dumps(result,separators=(',',':')).encode(),mtime=0))
    print('PROFILES',len(catalog),'RAYS',len(rays),'TRANSITIONS',Counter(r['digit'] for r in transitions),'CONTROLS',controls,flush=True)
    print('NEGATIVE_ONE_RAYS',[(a,rec['min_exponent'],rec['rows'],rec['root']) for (a,b),rec in sorted(rays.items()) if b==-1 and a in (3,5,9,15,27,45,81)],flush=True)
    print('SMALL_TRANSITIONS',[(r['odd_amplitude'],r['offset'],r['digit'],r['min_exponent'],r['source']['n'],r['target']['n']) for r in transitions[:35]],flush=True)

if __name__=='__main__':run()
