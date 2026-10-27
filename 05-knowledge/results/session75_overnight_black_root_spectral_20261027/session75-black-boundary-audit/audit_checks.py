#!/usr/bin/env python3
"""Independent Phase2 matrix primitives. No producer imports; explicit errors."""
from itertools import combinations

def need(ok, message):
    if not ok: raise RuntimeError(message)

def from_rows(rows):
    n=len(rows)
    need(all(isinstance(x,int) and 0<=x<1<<n for x in rows),'invalid row encoding')
    a=tuple(tuple(bool(rows[i]>>j&1) for j in range(n)) for i in range(n))
    need(all(not a[i][i] for i in range(n)) and all(a[i][j]!=a[j][i] for i in range(n) for j in range(i)),'invalid tournament')
    return a

def rows(a): return [sum(1<<j for j,v in enumerate(row) if v) for row in a]
def opposite(a): return tuple(zip(*a))
def degrees(a): return [sum(row) for row in a]
def strong(a, vertices=None):
    vertices=tuple(range(len(a))) if vertices is None else tuple(vertices)
    if not vertices:return False
    for converse in (False,True):
        seen={vertices[0]};stack=[vertices[0]]
        while stack:
            u=stack.pop()
            for v in vertices:
                if v not in seen and (a[v][u] if converse else a[u][v]):seen.add(v);stack.append(v)
        if len(seen)!=len(vertices):return False
    return True

def core(a):
    n=len(a)
    return [v for v,d in enumerate(degrees(a)) if d in (0,n-1) and strong(a,[u for u in range(n) if u!=v])]

def anti_failure(a,p):
    mirror=dict(zip(p,reversed(p)))
    for u in range(len(a)):
        for v in range(u):
            if a[u][v] != a[mirror[v]][mirror[u]]:return [u,v,mirror[v],mirror[u]]
    return None

def legal(a,p):
    need(sorted(p)==list(range(len(a))),'path not a permutation')
    need(all(a[u][v] for u,v in zip(p,p[1:])),'path not directed')

def reverse(a,p):
    legal(a,p);b=[list(row) for row in a]
    for u,v in zip(p,p[1:]):b[u][v]=False;b[v][u]=True
    return tuple(tuple(row) for row in b)

def native_k(a,p):
    legal(a,p);edges={frozenset((u,v)) for u,v in zip(p,p[1:])}
    n=len(a)
    return tuple(tuple(False if u==v else (a[u][v] if frozenset((u,v)) in edges else not a[u][v]) for v in range(n)) for u in range(n))

def paths(a, end=None):
    n=len(a)
    def extend(p,left):
        if not left:
            if end is None or p[-1]==end:yield tuple(p)
            return
        for v in sorted(left):
            if end==v and len(left)!=1:continue
            if not p or a[p[-1]][v]:yield from extend(p+[v],left-{v})
    yield from extend([],set(range(n)))

def cut(a,S):return sum(a[u][v] for u in S for v in range(len(a)) if v not in S)

def cut_certificate(a):
    """Exact min of each cardinality and every minimizing shore via cutoff ties."""
    n=len(a);ds=degrees(a);sd=sorted(ds);prefix=0;by_size=[]
    for k in range(1,n):
        prefix+=sd[k-1];by_size.append(prefix-k*(k-1)//2)
    lam=min(by_size);shores=[]
    for k,value in enumerate(by_size,1):
        if value!=lam:continue
        threshold=sd[k-1];fixed=[v for v in range(n) if ds[v]<threshold];ties=[v for v in range(n) if ds[v]==threshold]
        for C in combinations(ties,k-len(fixed)):shores.append(sum(1<<v for v in fixed+list(C)))
    return lam,min(ds+[n-1-d for d in ds]),sorted(shores),by_size

def replay(a,words, require_core=False, detail=False):
    """Native named start is a; no quotient or unrecorded source converse."""
    native=a;original=a;out=[]
    for j,p in enumerate(words):
        p=tuple(p);legal(a,p);bad=anti_failure(a,p);need(bad is not None,'blue R move')
        need(strong(a),'nonstrong intermediate')
        b=reverse(a,p);need(anti_failure(b,p[::-1]) is not None,'R inverse blue')
        q=p if j%2==0 else p[::-1]
        need(native==(a if j%2==0 else opposite(a)),'native gauge start mismatch')
        kn=native_k(native,q);need(kn==(opposite(b) if j%2==0 else b),'native K versus R mismatch')
        need(anti_failure(native,q) is not None and anti_failure(kn,q) is not None,'K endpoint blue')
        ds,de=degrees(a),degrees(b)
        need(de==[ds[v]+int(v==p[-1])-int(v==p[0]) for v in range(len(a))],'endpoint degree current mismatch')
        item={'R_path':p,'K_path':q,'black_failed_anti':bad,'target_strong':strong(b)}
        if detail:item.update(R_source_rows=rows(a),R_target_rows=rows(b),K_source_rows=rows(native),K_target_rows=rows(kn))
        out.append(item);a=b;native=kn
    need(not strong(a),'terminal remains strong');need(not strong(native),'native terminal remains strong')
    if require_core:need(core(a) and core(native),'terminal is not singleton plus strong core')
    return a,out
