#!/usr/bin/env python3
"""Exact controls; explicit checks remain enabled under python -O."""
import itertools,json,collections,hashlib
from pathlib import Path

def require(x,msg):
    if not x: raise RuntimeError(msg)

def controls(n):
    edges=list(itertools.combinations(range(n),2)); index={e:i for i,e in enumerate(edges)}; full=(1<<len(edges))-1
    perms=list(itertools.permutations(range(n)))
    def inv(p):
        rank={v:i for i,v in enumerate(p)}
        return sum(1<<k for k,(a,b) in enumerate(edges) if rank[a]>rank[b])
    maps=[]
    for p in perms:
        maps.append(([index[tuple(sorted((p[a],p[b])))] for a,b in edges],inv(p)))
    def graph_action(x,mp): return sum(1<<mp[k] for k in range(len(edges)) if (x>>k)&1)
    def cyc(x):
        deg=[0]*n
        for k,(a,b) in enumerate(edges): deg[b if (x>>k)&1 else a]+=1
        return n*(n-1)*(n-2)//6-sum(d*(d-1)//2 for d in deg)
    gc={};tc={}
    for x in range(full+1):
        gv=[graph_action(x,mp) for mp,_ in maps]
        tv=[a^s for a,(_,s) in zip(gv,maps)]
        gc[x]=min(gv);tc[x]=min(tv)
    gp=lambda x:min(gc[x],gc[full^x])
    tp=lambda x:min(tc[x],tc[full^x])
    overlap=collections.Counter((gp(x),tp(x)) for x in range(full+1))
    neighbors=collections.defaultdict(set)
    for (a,b),w in overlap.items(): neighbors[a].add(b)
    graphpairs=set(gp(x) for x in range(full+1));tourpairs=set(tp(x) for x in range(full+1))
    require(all(len(v)>=2 for k,v in neighbors.items() if k!=0),'nontrivial graph pair unexpectedly unique')
    # Connectedness of the bipartite overlap graph, including converse quotient.
    adj=collections.defaultdict(set)
    for a,b in overlap: adj['g',a].add(('t',b));adj['t',b].add(('g',a))
    seen={('g',0)};todo=list(seen)
    while todo:
        for y in adj[todo.pop()]-seen: seen.add(y);todo.append(y)
    require(len(seen)==len(graphpairs)+len(tourpairs),'overlap graph disconnected')
    four_count=0
    for x in range(1,full):
        def bit(a,b):return (x>>index[tuple(sorted((a,b)))])&1
        witness=next(((a,b,c) for a,b,c in itertools.permutations(range(n),3) if bit(a,c)!=bit(b,c)),None)
        require(witness is not None,'missing mixed triple')
        a,b,c=witness;rest=tuple(v for v in range(n) if v not in witness)
        orders=[(a,b,c)+rest,(b,a,c)+rest,(c,a,b)+rest,(c,b,a)+rest]
        vals=[cyc(x^inv(p)) for p in orders]
        require(abs((vals[1]-vals[0])-(vals[3]-vals[2]))==2,'four-order difference identity')
        require(len(set(vals))>1,'cyclic count did not change')
        four_count+=1
    # Pure and affine adjacent-transposition actions differ by one edge flip.
    if n>=2:
        p=(1,0)+tuple(range(2,n));mp,s=maps[perms.index(p)]
        require(s==1,'adjacent inversion must be single first edge')
        for x in range(full+1):
            y=graph_action(graph_action(x,mp),mp)^s
            require(y==x^1,'composition is not edge flip')
    return {'n':n,'graphs':len(set(gc.values())),'complement_graph_pairs':len(graphpairs),'tournaments':len(set(tc.values())),'converse_tournament_pairs':len(tourpairs),'overlap_edges':len(overlap),'graph_pair_degree_histogram':dict(sorted(collections.Counter(map(len,neighbors.values())).items())),'four_order_witnesses':four_count,'overlap_weights':[[a,b,w] for (a,b),w in sorted(overlap.items())]}

if __name__=='__main__':
    out={'status':'PASS','scope':'all labeled graphs n=2..5; all permutations; complete overlap incidence; all-size claims proved in note','orders':[controls(n) for n in range(2,6)]}
    print(json.dumps(out,indent=2,sort_keys=True))
