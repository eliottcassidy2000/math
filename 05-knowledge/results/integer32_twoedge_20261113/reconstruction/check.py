#!/usr/bin/env python3
"""Exact quadratic-jet image, repair, and independent inherited controls.
No producer imports. Prior immutable class metadata is the only census input.
"""
from functools import lru_cache
from itertools import combinations, permutations, product
from collections import Counter
from pathlib import Path
import hashlib
import json

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
DATA=ROOT/'05-knowledge/results/incidence_closure_20261113/closure/data.json'
DATA_SHA='c0484a99712edff88e8158b5f6000097f32f88874150d059da5aefa20f3de8f7'
CHECKS=0

def need(ok, label):
    global CHECKS
    CHECKS+=1
    if not ok: raise RuntimeError(label)

@lru_cache(None)
def edges(n): return tuple(combinations(range(n),2))

@lru_cache(None)
def wedges(n):
    return tuple((v,u,w) for v in range(n) for u,w in combinations([a for a in range(n) if a!=v],2))

@lru_cache(None)
def structure(n):
    E=edges(n); ids={e:i for i,e in enumerate(E)}
    links=tuple((ids[tuple(sorted((v,u)))],ids[tuple(sorted((v,w)))]) for v,u,w in wedges(n))
    base=sum(int(u<v<w)<<k for k,(v,u,w) in enumerate(wedges(n)))
    return links,base

@lru_cache(None)
def constraints(n):
    wi={w:i for i,w in enumerate(wedges(n))}
    def at(v,u,w): return 1<<wi[v,min(u,w),max(u,w)]
    out=[]
    for v in range(n):
        for a,b,c in combinations([a for a in range(n) if a!=v],3):
            out.append((at(v,a,b)^at(v,a,c)^at(v,b,c),0))
    for a,b,c in combinations(range(n),3):
        out.append((at(a,b,c)^at(b,a,c)^at(c,a,b),1))
    return tuple(out)

def encode(n,x):
    """x_e=1 reverses reference i->j (i<j); output full local wedge bits."""
    links,base=structure(n)
    return base ^ sum((((x>>i)^(x>>j))&1)<<k for k,(i,j) in enumerate(links))

def decode(n,beta):
    """Return canonical reversal mask x with bit0=0, or None if unrealizable."""
    if not 0<=beta<(1<<len(wedges(n))): return None
    if any((beta&row).bit_count()%2!=rhs for row,rhs in constraints(n)): return None
    M=len(edges(n))
    if not M:return 0
    links,base=structure(n);z=beta^base
    adj=[[] for _ in range(M)]
    for k,(i,j) in enumerate(links):
        adj[i].append((j,(z>>k)&1));adj[j].append((i,(z>>k)&1))
    values={0:0};todo=[0]
    while todo:
        i=todo.pop()
        for j,d in adj[i]:
            val=values[i]^d
            if j in values:
                if values[j]!=val:return None
            else:values[j]=val;todo.append(j)
    if len(values)!=M:return None
    x=sum(v<<i for i,v in values.items())
    return x if encode(n,x)==beta else None

@lru_cache(None)
def repair_paths(n):
    """2n-4 edge-disjoint length<=3 line-graph paths for every wedge."""
    E=edges(n);ei={e:i for i,e in enumerate(E)}
    links,_=structure(n);wi={tuple(sorted(pair)):k for k,pair in enumerate(links)}
    def e(a,b):return ei[min(a,b),max(a,b)]
    out=[]
    for v,u,w in wedges(n):
        a,b=e(v,u),e(v,w)
        paths=[(a,b),(a,e(u,w),b)]
        for z in range(n):
            if z not in (v,u,w):
                paths.extend(((a,e(v,z),b),(a,e(u,z),e(w,z),b)))
        out.append(tuple(sum(1<<wi[tuple(sorted((a,b)))] for a,b in zip(path,path[1:])) for path in paths))
    return tuple(out)

def repair(n,noisy):
    """Certified bounded-distance decoder. Returns x or None, radius n-3."""
    if n<3:return decode(n,noisy)
    if not 0<=noisy<(1<<len(wedges(n))):return None
    _,base=structure(n);z=noisy^base;corrected=0
    for k,paths in enumerate(repair_paths(n)):
        votes=sum((z&mask).bit_count()%2 for mask in paths)
        if votes==n-2:return None
        if votes>n-2:corrected|=1<<k
    x=decode(n,corrected^base)
    return x if x is not None and (encode(n,x)^noisy).bit_count()<=n-3 else None

def binary_rank(rows):
    pivots={}
    for row in rows:
        while row:
            k=row.bit_length()-1
            if k in pivots:row^=pivots[k]
            else:pivots[k]=row;break
    return len(pivots)

def image_controls():
    rows=[]
    for n in range(3,7):
        M=len(edges(n));L=len(wedges(n));base=encode(n,0)
        counts=Counter();images=set()
        for x in range(0,1<<M,2):
            beta=encode(n,x);images.add(beta)
            need(decode(n,beta)==x,'exact image inverse')
            counts[(beta^base).bit_count()]+=1
        need(len(images)==1<<(M-1),'all labelled converse pairs')
        distance=min(k for k in counts if k)
        need(distance==2*n-4,'exact minimum code distance')
        if n<=4:
            recognized={b for b in range(1<<L) if decode(n,b) is not None}
            need(recognized==images,'exhaustive arbitrary-field realizability')
        need(binary_rank([r for r,_ in constraints(n)])==L-M+1,'local-check codimension')
        rows.append({'n':n,'length':L,'dimension':M-1,'minimum_distance':distance,
                     'image_size':len(images),'weight_distribution':dict(sorted(counts.items()))})
    # Separate the star and vertex-triangle guards.
    need(not any((0&r).bit_count()%2 for r,b in constraints(3) if b==0),'all-zero passes stars')
    need(decode(3,0) is None,'all-zero triangle hostile')
    n=4;bad=encode(n,0);wi={w:i for i,w in enumerate(wedges(n))}
    bad^=(1<<wi[0,1,2])^(1<<wi[1,0,2])
    need(all((bad&r).bit_count()%2==b for r,b in constraints(n) if b==1),'triangle-only hostile passes')
    need(any((bad&r).bit_count()%2!=b for r,b in constraints(n) if b==0),'triangle-only hostile fails star')
    return rows

def error_controls():
    rows=[]
    for n in range(3,7):
        L=len(wedges(n));base=encode(n,0);paths=repair_paths(n)
        for family in paths:
            need(len(family)==2*n-4,'path family cardinality')
            union=0
            for path in family:
                need(not union&path and 1<=path.bit_count()<=3,'edge-disjoint short paths')
                union|=path
        total=0
        for t in range(n-2):
            for subset in combinations(range(L),t):
                error=sum(1<<i for i in subset)
                need(repair(n,base^error)==0,'all errors in guaranteed radius')
                total+=1
        # Translation/covariance positive controls around nonzero codewords.
        M=len(edges(n))
        probes=[2,((1<<M)-1)^1,sum(1<<i for i in range(1,M,2))]
        sample_errors=[0]+[sum(1<<i for i in subset) for subset in list(combinations(range(L),n-3))[:12]]
        for x in probes:
            for error in sample_errors:need(repair(n,encode(n,x)^error)==x,'nonzero-codeword covariance')
        diff=encode(n,1)^base;support=[i for i in range(L) if diff>>i&1]
        halfway=base^sum(1<<i for i in support[:n-2])
        need((halfway^base).bit_count()==(halfway^encode(n,1)).bit_count()==n-2,'sharp ambiguity radius')
        need(repair(n,halfway) is None,'ambiguity is not falsely certified')
        rows.append({'n':n,'all_error_patterns_at_zero':total,'correctable_radius':n-3,
                     'first_ambiguous_radius':n-2})
    return rows

def path_tables(n,code):
    E=edges(n);idx={e:i for i,e in enumerate(E)}
    def arc(u,v):return bool(code>>idx[min(u,v),max(u,v)]&1)==(u<v)
    table={0:[()]}
    for mask in range(1,1<<n):
        vertices=[i for i in range(n) if mask>>i&1]
        table[mask]=[p for p in permutations(vertices) if all(arc(u,v) for u,v in zip(p,p[1:]))]
    return table

def closure(adj,vertices):
    vertices=set(vertices);reach={u:({u}|(set(adj[u])&vertices)) for u in vertices}
    for k in vertices:
        for i in vertices:
            if k in reach[i]:reach[i]|=reach[k]
    return reach

def decode_core(n,positive,observed):
    E=edges(n);links,base=structure(n);used=[i for i in range(len(E)) if positive>>i&1]
    if not used:return 0
    adj={i:[] for i in used}
    for k,(i,j) in enumerate(links):
        if observed>>k&1:
            d=1^((base>>k)&1)
            adj[i].append((j,d));adj[j].append((i,d))
    vals={used[0]:0};todo=[used[0]]
    while todo:
        i=todo.pop()
        for j,d in adj[i]:
            v=vals[i]^d
            if j in vals:need(vals[j]==v,'one-slot orientation consistency')
            else:vals[j]=v;todo.append(j)
    need(len(vals)==len(used),'one-slot positive wedge graph connected')
    A=[set() for _ in range(n)]
    for i in used:
        u,v=E[i]
        if vals[i]:u,v=v,u
        A[u].add(v)
    reach=closure(A,range(n));remaining=set(range(n));blocks=[]
    while remaining:
        u=min(remaining);c={v for v in remaining if v in reach[u] and u in reach[v]};blocks.append(c);remaining-=c
    unsorted_blocks=tuple(blocks)
    blocks.sort(key=lambda c:sum(next(iter(c)) in reach[next(iter(d))] for d in unsorted_blocks if d is not c))
    out=[set(s) for s in A]
    for i,B in enumerate(blocks):
        for C in blocks[i+1:]:
            for u in B:out[u]|=C
        if len(B)>1:
            deleted={z:closure(A,B-{z}) for z in B}
            for u,v in combinations(sorted(B),2):
                if v in A[u] or u in A[v]:continue
                choices=set()
                for z in B-{u,v}:
                    uv=v in deleted[z][u];vu=u in deleted[z][v]
                    if uv!=vu:choices.add((u,v) if uv else (v,u))
                need(len(choices)==1,'dark arc deleted-reachability guard')
                a,b=choices.pop();out[a].add(b)
    code=sum(int(v in out[u])<<i for i,(u,v) in enumerate(E))
    return code

def actual_controls():
    raw=DATA.read_bytes();need(hashlib.sha256(raw).hexdigest()==DATA_SHA,'immutable census data pin')
    universe=json.loads(raw);rows=[]
    for r in universe['rows']:
        n=r['n']
        if n>6:continue
        E=edges(n);ei={e:i for i,e in enumerate(E)};wi={w:i for i,w in enumerate(wedges(n))};mask=(1<<len(E))-1
        htotal=0;dark=0
        for node in r['ordinary_nodes']:
            c=node['code'];tables=path_tables(n,c);hp=tables[(1<<n)-1];htotal+=len(hp)
            need(len(hp)==node['H'],'literal Hamilton paths vs immutable metadata')
            positive=observed=0
            for p in hp:
                for u,v in zip(p,p[1:]):positive|=1<<ei[min(u,v),max(u,v)]
                for u,v,w in zip(p,p[1:],p[2:]):observed|=1<<wi[v,min(u,w),max(u,w)]
            got=decode_core(n,positive,observed)
            need(got in (c,c^mask),'one-slot recovery including dark arcs')
            # Independent two-slot enumeration by disjoint subset and path products.
            beta=0
            full=(1<<n)-1
            for S in range(1<<n):
                for p,q in product(tables[S],tables[full^S]):
                    for word in (p,q):
                        for u,v,w in zip(word,word[1:],word[2:]):beta|=1<<wi[v,min(u,w),max(u,w)]
            need(beta==encode(n,c^mask),'literal two-slot support equals local mixed incidences')
            dark+=len(E)-positive.bit_count()
        rows.append({'n':n,'all_ordinary_classes':len(r['ordinary_nodes']),'Hamilton_paths':htotal,'omitted_positive_arcs_total':dark})
    # One-slot positivity is NOT the full two-slot code: transitive order four.
    n=4;c=(1<<6)-1;hp=path_tables(n,c)[15];wi={w:i for i,w in enumerate(wedges(n))}
    observed=sum(1<<wi[v,u,w] for u,v,w in zip(hp[0],hp[0][1:],hp[0][2:]))
    need(decode(4,observed) is None,'one-slot field lies outside full two-slot code')
    return rows

def module_controls():
    n=6;E=edges(n)
    A=[{1,3,4},{2,3,4},{0,3,4},{4,5},{5},{0,1,2}]
    B=[set(a) for a in A]
    for i,j in combinations(range(3),2):
        if j in B[i]:B[i].remove(j);B[j].add(i)
        else:B[j].remove(i);B[i].add(j)
    tables=[]
    for T in (A,B):
        code=sum(int(j in T[i])<<k for k,(i,j) in enumerate(E))
        tables.append(path_tables(n,code))
    need([len(tables[0][S]) for S in range(1<<n)]==[len(tables[1][S]) for S in range(1<<n)],'all induced H under module reversal')
    wi={w:i for i,w in enumerate(wedges(n))}
    def moments(table,q):
        count=0;edge=Counter();start=Counter();end=Counter();wedge=Counter()
        for assignment in product(range(q),repeat=n):
            masks=[sum(1<<v for v in range(n) if assignment[v]==i) for i in range(q)]
            for words in product(*(table[S] for S in masks)):
                count+=1
                for i,word in enumerate(words):
                    if word:start[i,word[0]]+=1;end[i,word[-1]]+=1
                    for u,v in zip(word,word[1:]):edge[i,min(u,v),max(u,v)]+=1
                    for u,v,w in zip(word,word[1:],word[2:]):wedge[i,wi[v,min(u,w),max(u,w)]]+=1
        return count,edge,start,end,wedge
    results=[]
    for q in (1,2,3):
        left,right=(moments(t,q) for t in tables)
        need(left[:4]==right[:4],'independent per-slot first-moment control')
        need(left[4]!=right[4],'wedge marks see reversed module')
        results.append({'q':q,'slot_covers':left[0],'first_moments_equal':True,'wedge_moments_different':True})
    return {'scope':'six-vertex labelled toy with a self-converse C3 module; not the nonisomorphic order-24 witness','results':results}

def main():
    result={'status':'FINITE-EXACT controls; new statements RESERVED pending independent review',
            'image':image_controls(),'error_correction':error_controls(),'actual_tournaments':actual_controls(),'module_invariance':module_controls(),
            'input_sha256':DATA_SHA}
    result['checks']=CHECKS;result['source_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    print(json.dumps(result,sort_keys=True,indent=2))
if __name__=='__main__':main()
