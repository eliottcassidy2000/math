#!/usr/bin/env python3
"""Literal finite controls for fixed-block shared-root entropy theorem."""
import json, random, math
from functools import lru_cache

def need(ok, message):
    if not ok: raise RuntimeError(message)

def induced(rows, vertices):
    return tuple(sum(1<<j for j,v in enumerate(vertices) if rows[u]>>v&1) for u in vertices)

def hamilton(rows):
    n=len(rows)
    if not n: return 1
    dp=[None]*(1<<n)
    for i in range(n):
        a=[0]*n; a[i]=1; dp[1<<i]=a
    for mask in range(1,1<<n):
        a=dp[mask]
        if a is None: continue
        for i,v in enumerate(a):
            if not v: continue
            avail=rows[i]&~mask
            while avail:
                bit=avail&-avail; avail-=bit; j=bit.bit_length()-1
                nxt=mask|bit
                if dp[nxt] is None: dp[nxt]=[0]*n
                dp[nxt][j]+=v
    return sum(dp[-1])

def bins(rows):
    n=len(rows); allmask=(1<<n)-1
    val={mask:hamilton(induced(rows,[i for i in range(n) if mask>>i&1])) for mask in range(1<<n)}
    return sum(val[mask]*val[allmask^mask] for mask in val)

def glue(core,m,rootout):
    s=len(core); n=1+s*m; rows=[0]*n
    for i in range(s*m):
        v=i+1
        if rootout>>i&1: rows[0]|=1<<v
        else: rows[v]|=1
        block,j=divmod(i,s)
        rows[v]|=sum(1<<(1+block*s+k) for k in range(s) if core[j]>>k&1)
        rows[v]|=sum(1<<(1+k) for k in range((block+1)*s,s*m))
    return tuple(rows)

def root_component(rows):
    n=len(rows); reach=[rows[i]|1<<i for i in range(n)]
    for k in range(n):
        for i in range(n):
            if reach[i]>>k&1: reach[i]|=reach[k]
    return [i for i in range(n) if reach[0]>>i&1 and reach[i]&1]

def test(core,m,masks):
    s=len(core); h=hamilton(core); beta=bins(core)
    strong=0; minimum=None; maximum=0
    for mask in masks:
        rows=glue(core,m,mask); val=hamilton(rows); comp=root_component(rows)
        nonroot=[v-1 for v in comp if v]
        blocks=sorted({v//s for v in nonroot}); ell=len(blocks)
        if blocks:
            need(blocks==list(range(blocks[0],blocks[-1]+1)), 'noncontiguous root component')
            need(len(nonroot)==s*ell,'partial strong block in component')
        cv=hamilton(induced(rows,comp))
        need(val==h**(m-ell)*cv, 'SCC count factorization')
        if ell>=2: need(cv>=beta**(ell-2),'strong two-bin floor')
        if len(comp)==len(rows):
            strong+=1; minimum=val if minimum is None else min(minimum,val)
        maximum=max(maximum,val)
    return {'blocks':m,'presentations':len(masks),'strong_presentations':strong,'minimum_strong_h':minimum,'maximum_h':maximum}

def main():
    cores=[('C3',(2,4,1),3,12,3),('D4',(8,1,3,6),5,32,2)]
    results={}; checks=0
    for name,rows,h,beta,limit in cores:
        need(hamilton(rows)==h,'core H');need(bins(rows)==beta,'core two-bin count')
        need(len(root_component(rows))==len(rows),'strong core')
        need(beta>2**len(rows),'strict entropy gap')
        tests=[]
        for m in range(limit+1):
            masks=list(range(1<<(m*len(rows))))
            tests.append(test(rows,m,masks));checks+=len(masks)
        results[name]={'rows':rows,'order':len(rows),'H':h,'P2':beta,'rho_approximation_only':len(rows)*math.log(2)/math.log(beta),'complete_universes':tests}
    rng=random.Random(20261106)
    masks=sorted({0,(1<<12)-1,1,(1<<11),0x555,0xaaa}|set(rng.randrange(1<<12) for _ in range(32)))
    results['D4']['additional_declared_sample']=test((8,1,3,6),3,masks);checks+=len(masks)
    print(json.dumps({'status':'PASS finite literal controls, not target coverage','graph_presentations_checked':checks,'results':results},indent=2,sort_keys=True))
if __name__=='__main__': main()
