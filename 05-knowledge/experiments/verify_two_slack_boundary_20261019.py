#!/usr/bin/env python3
"""Three prescribed fixtures; no census or numeric probability inference."""
from itertools import combinations
import json
checks=0

def require(ok,msg):
    global checks
    checks+=1
    if not ok:raise RuntimeError(msg)

def close(n,edges):
    p=[0]*n
    for a,b in edges:p[b]|=1<<a
    for k in range(n):
        for b in range(n):
            if p[b]>>k&1:p[b]|=p[k]
    require(all(not(p[i]>>i&1) for i in range(n)),'acyclic')
    return p

def data(p):
    n=len(p);u=[sum(1<<j for j in range(n) if p[j]>>i&1) for i in range(n)]
    inc=[((1<<n)-1)&~(p[i]|u[i]|1<<i) for i in range(n)]
    return u,inc

def chain(bits,p):
    return all((p[a]>>b|p[b]>>a)&1 for a,b in combinations([i for i in range(len(p)) if bits>>i&1],2))

def vg(a,b,p):
    u,inc=data(p)
    return bool(inc[a]>>b&1) and ((p[a]==p[b] and chain(u[a]&~u[b],p) and chain(u[b]&~u[a],p)) or (u[a]==u[b] and chain(p[a]&~p[b],p) and chain(p[b]&~p[a],p)))

def antichains(p):
    u,inc=data(p);found=[0]
    for i in range(len(p)):
        found += [s|1<<i for s in found if s&~inc[i]==0]
    w=max(s.bit_count() for s in found)
    return w,[s for s in found if s.bit_count()==w]

def fixture(nD,nU,lo,hi,lower_edges=()):
    A=list(range(6));D=list(range(6,6+nD));U=list(range(6+nD,6+nD+nU))
    edges=[(d,u) for d in D for u in U]
    edges += [(d,a) for a in A for j,d in enumerate(D) if j not in lo[a]]
    edges += [(a,u) for a in A for j,u in enumerate(U) if j not in hi[a]]
    edges += [(D[i],D[j]) for i,j in lower_edges]
    return close(6+nD+nU,edges),A,D,U

def check(name,p,A,D,U,unique=False,expected_lower_comparable=0):
    ups,inc=data(p)
    require(max(v.bit_count() for v in inc)<=7,'seven-thin')
    require(not any(vg(a,b,p) for a,b in combinations(A,2)),'no very good pair in A')
    require(any(vg(a,b,p) for a,b in combinations(range(len(p)),2)),'ambient very good pair exists')
    supports={v:set(a for a in A if inc[v]>>a&1) for v in D+U}
    require(all(len(x)<=2 for x in supports.values()),'boundary exterior degree')
    require(all(sum(inc[a]>>v&1 for v in D+U)<=2 for a in A),'A defect degree')
    require(all(p[u]>>d&1 for d in D for u in U),'D below U')
    dc=[(a,b) for a,b in combinations(D,2) if (p[a]>>b|p[b]>>a)&1]
    uc=[(a,b) for a,b in combinations(U,2) if (p[a]>>b|p[b]>>a)&1]
    require(len(dc)==expected_lower_comparable and len(uc)==0,'almost antichain boundary')
    width,maxa=antichains(p)
    require(width==6,'maximum antichain six')
    if unique:require(maxa==[sum(1<<a for a in A)],'unique maximum antichain')
    out={'name':name,'vertices':len(p),'predecessor_masks':p,'degrees':[i.bit_count() for i in inc],'width':width,'maximum_antichains':len(maxa),'very_good_pairs_in_A':0,'ambient_very_good_pairs':sum(vg(a,b,p) for a,b in combinations(range(len(p)),2)),'lower_comparable_pairs':len(dc)}
    return out

def main():
    # The unique maximum antichain has five mixed rows and one pure upper row.
    p,A,D,U=fixture(5,5,[{i} for i in range(5)]+[set()],[{i} for i in range(5)]+[{0,1}])
    out=[check('unique_maximum_antichain_local_obstruction',p,A,D,U,True)]
    require(all(vg(a,b,p) for a,b in combinations(D,2)),'lower boundary clique')
    # All-mixed six-by-six crown interfaces exhibit the strict downward move.
    p,A,D,U=fixture(6,6,[{i} for i in range(6)],[{i} for i in range(6)])
    out.append(check('all_mixed_descent',p,A,D,U))
    idealA=set(A)|set().union(*[{j for j in range(len(p)) if p[a]>>j&1} for a in A])
    idealD=set(D)|set().union(*[{j for j in range(len(p)) if p[d]>>j&1} for d in D])
    require(idealD<idealA,'strict downideal descent')
    # A chain doubleton really is allowed by local no-very-good-pair tests.
    p,A,D,U=fixture(7,5,[{0,1}]+[{i+1} for i in range(1,6)],[set()]+[{i-1} for i in range(1,6)],[(0,1)])
    out.append(check('permitted_single_chain_doubleton',p,A,D,U,expected_lower_comparable=1))
    print(json.dumps({'status':'ALL_SELECTED_CONTROLS_PASSED','guards':checks,'fixtures':out},indent=2))

if __name__=='__main__':main()
