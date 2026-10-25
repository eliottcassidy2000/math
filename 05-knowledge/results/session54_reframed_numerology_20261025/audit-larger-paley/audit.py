#!/usr/bin/env python3
"""Independent deterministic certificate audit. No producer code imported."""
import json,pathlib,math
ROOT=pathlib.Path(__file__).resolve().parent

def require(value,message):
    if not value: raise ValueError(message)

def prime(n):
    return n>=2 and all(n%d for d in range(2,math.isqrt(n)+1))

def column_determinant(matrix,p):
    # Column rather than row elimination: combine whole columns to clear a row.
    columns=[list(x) for x in zip(*matrix)]; n=len(columns); answer=1
    for i in range(n):
        k=next((k for k in range(i,n) if columns[k][i]%p),None)
        if k is None: return 0
        if k!=i: columns[i],columns[k]=columns[k],columns[i]; answer=-answer
        v=columns[i][i]%p; answer=answer*v%p; inv=pow(v,p-2,p)
        for j in range(i+1,n):
            factor=columns[j][i]*inv%p
            if factor:
                columns[j][i:]=[(x-factor*y)%p for x,y in zip(columns[j][i:],columns[i][i:])]
    return answer%p

def verify(cert):
    n=cert['n']; p=cert['modulus']; require(n in (19,23,31),'order scope'); require(p==101 and prime(p) and prime(n),'prime fields')
    r=(n-1)//2; lam=(n-3)//4; V=set(range(n)); Q={a*a%n for a in range(1,n)}
    E=[(u,v) for u in range(n) for v in range(n) if (v-u)%n in Q]; edges=set(E)
    require(len(Q)==r and Q.isdisjoint({-a%n for a in Q}) and len(E)==n*r,'tournament')
    out=[{v for u,v in E if u==x} for x in range(n)]; ins=[{u for u,v in E if v==x} for x in range(n)]
    require(all(len(x)==r for x in out+ins),'regularity')
    require(all(len(ins[a]&ins[b])==lam and len(out[a]&out[b])==lam for a in V for b in V if a!=b),'double regularity')
    group=[tuple((a*x+b)%n for x in range(n)) for a in sorted(Q) for b in range(n)]
    require(len(set(group))==n*r and tuple(range(n)) in group,'group order identity')
    require(all(a*b%n in Q and pow(a,-1,n) in Q for a in Q for b in Q),'group composition and inverses')
    require(all({(g[u],g[v]) for u,v in E}==edges for g in group),'all affine maps preserve arcs')
    require({(g[0],g[1]) for g in group}==edges,'arc transitivity')
    adj=[set() for _ in range(2*n)]
    for u,v in E: adj[u].add(n+v); adj[n+v].add(u)
    seen={0}; queue=[0]; tree=[]
    for u in queue:
        for v in sorted(adj[u]-seen): seen.add(v);queue.append(v);tree.append([u,v])
    require(len(seen)==2*n and len(tree)==2*n-1,'connected support')
    def vector(cycle):
        require(len(cycle)==n and set(cycle)==V,'Hamilton cycle vertex set')
        used=set(zip(cycle,cycle[1:]+cycle[:1])); require(used<=edges,'Hamilton cycle arcs')
        return [int(e in used) for e in E]
    base=cert['base_cycle']; b=vector(base); require(cert['seed_cycles']==[base],'one seed only')
    orbit={tuple(g[x] for x in base) for g in group}
    require(all(tuple(c) in orbit for c in cert['independent_cycles']),'selected cycles in first seed affine orbit')
    rows=[[x-y for x,y in zip(vector(c),b)] for c in cert['independent_cycles']]
    target=n*r-2*n+1; require(cert['attained_rank']==cert['target_rank']==len(rows)==target,'rank target')
    # Explicit annihilation gives the matching rational upper bound via connectivity.
    for row in rows:
        require(all(sum(row[i] for i,(u,v) in enumerate(E) if u==x)==0 and sum(row[i] for i,(u,v) in enumerate(E) if v==x)==0 for x in V),'row column annihilation')
    pivots=cert['pivot_columns']; require(len(pivots)==len(set(pivots))==target and all(type(j)==int and 0<=j<len(E) for j in pivots),'minor indices')
    determinant=column_determinant([[row[j] for j in pivots] for row in rows],p)
    require(determinant!=0,'nonzero minor')
    path=cert['endpoint_path']; require(len(path)==n and set(path)==V and path[0]==0 and path[-1]==1 and set(zip(path,path[1:]))<=edges,'forward endpoint witness')
    require(all(r*(r-2)+lam*(t-1)*(2*r-1-t)>0 for t in range(1,2*r)),'deficiency one interior excluded')
    require(all(r*(r-3)+lam*(t-1)*(2*r-2-t)>0 for t in range(1,2*r-1)),'deficiency two interior excluded')
    # Family distinctness is a separate elementary finite check, NOT blocker testing.
    family=set()
    for a in range(n):
        for b0 in range(a):
            family.add(frozenset((u,v) for u,v in E if u in (a,b0)))
            family.add(frozenset((u,v) for u,v in E if v in (a,b0)))
        family.add(frozenset((u,v) for u,v in E if u==a or v==a))
    require(len(family)==n*n and all(len(f)==n-1 for f in family),'Hall star family sizes and distinctness')
    return dict(n=n,arcs=len(E),affine_group_order=len(group),listed_cycles=1+len(rows),rank_Q=target,minor_mod101=determinant,support_tree=tree,beta_HP=n-1,theorem_minimum_blockers=n*n,two_source=n*(n-1)//2,two_sink=n*(n-1)//2,stars=n,discovery='heuristic seeded randomized DFS; not replayed',verification='deterministic polynomial finite checks; no Hamilton path counting or blocker enumeration')

def main():
    results=[verify(json.loads((ROOT/'inputs/session54-larger-paley-certificates'/f'certificate-{n}.json').read_text())) for n in (19,23,31)]
    print(json.dumps({'status':'PASS','results':results},indent=2,sort_keys=True))
if __name__=='__main__': main()
