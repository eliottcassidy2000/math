"""Complete recursive colouring carriers; native complement and joint ports."""
from collections import Counter
from functools import lru_cache
from itertools import product
from pathlib import Path
import importlib.util
import json

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('weighted',HERE/'check_weighted.py')
w=importlib.util.module_from_spec(spec); spec.loader.exec_module(w)
CHECKS=0
def need(ok,msg):
    global CHECKS
    CHECKS+=1
    if not ok: raise RuntimeError(msg)

@lru_cache(None)
def ca(values,k):
    if not k: return values[0]
    d=3**(k-1); children=tuple(ca(values[i*d:(i+1)*d],k-1) for i in range(3))
    return min(children,children[1:]+children[:1],children[2:]+children[:2])

@lru_cache(None)
def cc(values,k):
    if not k: return values[0]
    d=3**(k-1); t=w.negation(k-1)
    paired=tuple((values[d+i],values[2*d+t[i]]) for i in range(d))
    return (cc(values[:d],k-1),ca(paired,k-1))

def moved(values,p):
    out=[None]*len(p)
    for i,x in enumerate(p): out[x]=values[i]
    return tuple(out)

def complement(values,k):
    t=w.negation(k)
    return tuple(1-values[t[i]] for i in range(len(values)))

def parent_rows(k):
    n=3**k; rows=[]
    for x in range(n):
        mask=0
        for y in range(n):
            if x==y: continue
            for d in reversed(range(k)):
                a,b=(x//3**d)%3,(y//3**d)%3
                if a!=b:
                    if (b-a)%3==1: mask|=1<<y
                    break
        rows.append(mask)
    return tuple(rows)

def strong(rows):
    n=len(rows)
    for source in range(n):
        reached=1<<source
        while True:
            after=reached
            for v in range(n):
                if reached>>v&1: after|=rows[v]
            if after==reached: break
            reached=after
        if reached!=(1<<n)-1: return False
    return True

def main():
    rows=[]
    for k in range(3):
        n=3**k; ga=w.permutations_a(k); t=w.negation(k)
        gc=[p for p in ga if all(p[t[i]]==t[p[i]] for i in range(n))]
        keys_a={}; keys_c={}; complement_keys={}; seen_c=set(); poly=Counter(); parent=parent_rows(k)
        for values in product([0,1],repeat=n):
            ka,kc=ca(values,k),cc(values,k)
            actual_a=min(moved(values,p) for p in ga)
            actual_c=min(moved(values,p) for p in gc)
            need(keys_a.setdefault(ka,actual_a)==actual_a,'A canonical completeness')
            need(keys_c.setdefault(kc,actual_c)==actual_c,'C canonical completeness')
            flipped=cc(complement(values,k),k)
            need(flipped!=kc,'complement free on colouring orbits')
            need(complement_keys.setdefault(kc,flipped)==flipped,'complement descends')
            if k:
                child=tuple(row|((1-values[i])<<n) for i,row in enumerate(parent))+(sum(v<<i for i,v in enumerate(values)),)
                need(strong(child)==(0<sum(values)<n),'actual strong-atom attachment')
            if kc not in seen_c: seen_c.add(kc); poly[sum(values)]+=1
        expected=w.color_poly(w.actual_index(gc))
        need([poly[i] for i in range(n+1)]==expected,'canonical and weighted counts')
        root_zero={cc(values,k) for values in product([0,1],repeat=n) if values[0]==0}
        need(len(root_zero)*2==len(keys_c),'one root bit selects complement orientation')
        ports={}
        for r in range(1,4):
            signatures=set()
            for locations in product(range(n),repeat=r):
                values=tuple(sum(1<<i for i,x in enumerate(locations) if x==v) for v in range(n))
                signatures.add(cc(values,k))
            ports[r]=len(signatures)
        need([ports[i] for i in range(1,4)]==[2*k+1,8*k*k+1,48*k**3-40*k*k+18*k+1], 'joint named ports')
        rows.append({'depth':k,'vertices':n,'A_coloring_orbits':len(keys_a),'C_coloring_orbits':len(keys_c),
                     'complement_orbits':len(root_zero),'ordered_port_orbits':ports})
    # Same one-port types, but the paired diagonal action sees relative offset.
    n=9
    def portkey(xs): return cc(tuple(sum(1<<i for i,x in enumerate(xs) if x==v) for v in range(n)),2)
    need(portkey((3,))==portkey((4,)) and portkey((6,))==portkey((7,)),'equal marginal port types')
    need(portkey((3,6))!=portkey((3,7)),'missing joint phase')
    result={'status':'FINITE-EXACT canonical checks for proved all-depth recurrence',
            'rows':rows,'joint_hostile':{'pairs':[[3,6],[3,7]],'equal_marginals':True,'distinct_joint':True},
            'checks':CHECKS}
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__': main()
