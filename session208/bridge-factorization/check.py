"""Exact finite controls for the bridge factorization theorem."""
from itertools import combinations
from fractions import Fraction
import json

def components(n, edges, keep):
    adj=[[] for _ in range(n)]
    for e in keep:
        u,v=edges[e]; adj[u].append(v);adj[v].append(u)
    labels=[None]*n; result=[]
    for v in range(n):
        if labels[v] is not None: continue
        i=len(result); stack=[v]; labels[v]=i; group=[]
        while stack:
            u=stack.pop();group.append(u)
            for w in adj[u]:
                if labels[w] is None: labels[w]=i;stack.append(w)
        result.append(group)
    return labels,result

def paths(n,edges,tree):
    adj=[[] for _ in range(n)]
    for e in tree:
        u,v=edges[e];adj[u].append((v,e));adj[v].append((u,e))
    out=[]
    for u,v in edges:
        stack=[(u,-1,())]
        while stack:
            x,parent,p=stack.pop()
            if x==v:out.append(p);break
            for y,e in adj[x]:
                if y!=parent:stack.append((y,x,p+(e,)))
    return out

def test(n,edges):
    m=len(edges); all_e=set(range(m))
    bridges={e for e in all_e if len(components(n,edges,all_e-{e})[1])>1}
    labels, groups=components(n,edges,all_e-bridges)
    lengths=[2*e+1 for e in range(m)];caps=[3*e+2 for e in range(m)]
    count=0
    for tup in combinations(range(m),n-1):
        tree=set(tup)
        if len(components(n,edges,tree)[1])!=1:continue
        count+=1;assert bridges<=tree
        ps=paths(n,edges,tree)
        loads=[sum(caps[e] for e in range(m) if f in ps[e]) for f in range(m)]
        for b in bridges:assert ps[b]==(b,) and loads[b]==caps[b]
        for i,vertices in enumerate(groups):
            internal={e for e,(u,v) in enumerate(edges) if labels[u]==i==labels[v]}
            assert len(tree&internal)==len(vertices)-1
            for e in internal: assert set(ps[e])<=internal
            for f in internal:
                assert loads[f]==sum(caps[e] for e in internal if f in ps[e])
    assert count
    return dict(vertices=n,edges=m,bridges=len(bridges),trees=count)

def coupling():
    distributions=[[Fraction(1,3),Fraction(2,3)], [Fraction(1,4),Fraction(1,2),Fraction(1,4)], [Fraction(2,5),Fraction(3,5)]]
    ends=[]
    for ds in distributions:
        cum=[Fraction(0)]
        for p in ds:cum.append(cum[-1]+p)
        ends.append(cum)
    bounds=sorted(set(x for c in ends for x in c));marginals=[[Fraction(0) for _ in d] for d in distributions]
    for a,b in zip(bounds,bounds[1:]):
        for i,c in enumerate(ends):
            j=next(j for j in range(len(c)-1) if c[j]<=a<c[j+1]);marginals[i][j]+=b-a
    assert marginals==distributions
    assert len(bounds)-1<=1+sum(len(d)-1 for d in distributions)
    return len(bounds)-1

cases=[(2,[(0,1)]),(4,[(0,1),(1,2),(2,3)]),(5,[(0,1),(1,2),(2,0),(2,3),(3,4)]),(7,[(0,1),(1,2),(2,0),(2,3),(3,4),(4,5),(5,3),(5,6)]),(4,[(0,1),(0,1),(1,2),(2,3),(2,3)])]
print(json.dumps({'status':'PASS','cases':[test(n,e) for n,e in cases],'rational_coupling_atoms':coupling()},indent=2))
