#!/usr/bin/env python3
import json, math, itertools, hashlib
from pathlib import Path
from collections import Counter
from sympy import Matrix, ZZ
from sympy.matrices.normalforms import smith_normal_form
ROOT=next(p for p in [*Path(__file__).resolve().parents, Path(__file__).resolve().parent.parent/'math'] if (p/'AGENTS.md').is_file() and (p/'01-canon').is_dir())
SRC=ROOT/'05-knowledge/results/merged_metagraph_black_quotient_bridges_20261020'
OUT=Path(__file__).parent

def rank2(rows):
    basis={}
    for z in rows:
        while z:
            p=z.bit_length()-1
            if p in basis: z^=basis[p]
            else: basis[p]=z; break
    return len(basis)

def lap(n,es):
    a=[[0]*n for _ in range(n)]
    for (u,v),w in es.items():
        if u==v:continue
        a[u][u]+=w;a[v][v]+=w;a[u][v]-=w;a[v][u]-=w
    return a

def comps(n,es):
    adj=[set() for _ in range(n)]
    for (u,v),w in es.items():
        if w:adj[u].add(v);adj[v].add(u)
    seen=set();cs=[]
    for u in range(n):
        if u in seen:continue
        c=set();q=[u]
        while q:
            v=q.pop()
            if v in c:continue
            c.add(v);seen.add(v);q.extend(adj[v]-c)
        cs.append(sorted(c))
    return cs

def poly_trees(n,B,Q):
    # Independent explicit tree subsets in the graph whose black edges are R-orbits.
    edges=[(u,v,w,0) for (u,v),w in B.items() if u!=v]+[(u,v,w,1) for (u,v),w in Q.items() if u!=v]
    p=Counter()
    for subset in itertools.combinations(edges,n-1):
        uf=list(range(n));weight=1;k=0
        def get(x):
            while x!=uf[x]:x=uf[x]
            return x
        for u,v,w,b in subset:
            ru,rv=get(u),get(v)
            if ru==rv:break
            uf[ru]=rv;weight*=w;k+=b
        else:p[k]+=weight
    return dict(sorted(p.items()))

def direct_fibres(n,tiles):
    size=1<<len(tiles); canon=[];perms=list(itertools.permutations(range(n)))
    for t in range(size):
        arcs={(i,i+1) for i in range(n-1)}
        for k,(x,y) in enumerate(tiles):
            i,j=n-x,n-y
            arcs.add((j,i) if (t>>k)&1 else (i,j))
        fibre=set()
        for p in perms:
            if all((p[i],p[i+1]) in arcs for i in range(n-1)):
                fibre.add(sum(int((p[n-y],p[n-x]) in arcs)<<k for k,(x,y) in enumerate(tiles)))
        canon.append(min(fibre))
    return canon

records=[]
for n in range(3,8):
    path=SRC/f'census_n{n}.json';data=json.loads(path.read_text())
    nodes=data['nodes'];N=len(nodes);pi=data['pi'];tiles=[tuple(t) for t in data['tiles']]; ix={p:i for i,p in enumerate(tiles)}
    ref=[ix[(n+1-y,n+1-x)] for x,y in tiles];full=len(pi)-1
    B=Counter();K=Counter()
    for t in range(len(pi)//2):
        rt=sum(((t>>i)&1)<<j for i,j in enumerate(ref));pair=tuple(sorted((pi[t],pi[full^t])))
        (B if rt==t else K)[pair]+=1
        if pi[rt]!=pi[t]:raise RuntimeError('Owner reflection mismatch')
    if B!=Counter({(u,v):w for u,v,w in data['blue_edges']}) or K!=Counter({(u,v):w for u,v,w in data['black_edges']}):raise RuntimeError('Edge mismatch')
    if any(w%2 for w in K.values()):raise RuntimeError('Unpaired black')
    Q=Counter({e:w//2 for e,w in K.items()});G=B+K
    if len(comps(N,G))!=1:raise RuntimeError('Expected finite connected control')
    if n<=5:
        cc=direct_fibres(n,tiles)
        known=[data['classes'][ci]['representative'] for ci in data['class_by_mask']]
        if cc!=known:raise RuntimeError('Direct permutation classifier failed')
    LB=lap(N,B);LG=lap(N,G);cb=len(comps(N,B));NS=sum(not r['SC'] for r in nodes)
    rank=rank2([sum((v%2)<<i for i,v in enumerate(row)) for row in LB]);bicycle=N-cb-rank;rho=N-1-rank
    rank_reduced=rank2([sum((v%2)<<i for i,v in enumerate(row[1:])) for row in LG[1:]])
    if rho != N-1-rank_reduced:raise RuntimeError('Reduced/full parity rank mismatch')
    r={'n':n,'vertices':N,'NS':NS,'blue_components':cb,'blue_mod2_rank':rank,'blue_bicycle_dimension':bicycle,'critical_group_2_rank':rho,'source_sha256':hashlib.sha256(path.read_bytes()).hexdigest()}
    if n<=6:
        cof=Matrix([row[1:] for row in LG[1:]])
        det=int(cof.det(method='domain-ge'));r['spanning_trees']=det;r['spanning_tree_v2']=(det&-det).bit_length()-1
        snf=smith_normal_form(cof,domain=ZZ);inv=[abs(int(snf[i,i])) for i in range(N-1)]
        if sum(x%2==0 for x in inv)!=rho:raise RuntimeError('Smith mismatch')
        r['smith_invariants']=inv
    if n<=5:
        p=poly_trees(N,B,Q);r['tree_polynomial_B_plus_xQ']=p
        if sum(v*2**k for k,v in p.items())!=r['spanning_trees']:raise RuntimeError('Tree polynomial mismatch')
        bc=comps(N,B);owner={v:i for i,c in enumerate(bc) for v in c};cq=Counter()
        for (u,v),w in Q.items():
            a,b=sorted((owner[u],owner[v]));cq[a,b]+=w
        TC=int(Matrix([row[1:] for row in lap(cb,cq)[1:]]).det()) if cb>1 else 1
        prod=TC
        for c in bc:
            arr=[[LB[u][v] for v in c[1:]] for u in c[1:]]
            prod*=int(Matrix(arr).det()) if arr else 1
        if p.get(cb-1,0)!=prod:raise RuntimeError('Leading coefficient mismatch')
        r['minimum_black_degree']=min(p);r['leading_tree_coefficient']=prod
    records.append(r)
    print(json.dumps(r),flush=True)
(OUT/'probe.json').write_text(json.dumps({'records':records},indent=2)+'\n')
