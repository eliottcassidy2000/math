#!/usr/bin/env python3
"""Independent finite corroboration. No producer imports or assert statements."""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations, product
from math import isqrt, comb
import json

def require(c,label):
    if not c: raise RuntimeError(label)

def divisors(n):
    ds=[]
    for d in range(1,isqrt(n)+1):
        if n%d==0:
            ds.append(d)
            if d*d!=n: ds.append(n//d)
    return sorted(ds)

def aliquot():
    hits=[]; classes=Counter()
    for n in range(1,10001):
        ds=divisors(n); s=sum(ds)-n
        pairs=[(a,n//a) for a in ds if a<=n//a and 1+a+n//a==s]
        require(bool(pairs)==(len(ds)==4 or s==n+2),f'factor compression iff {n}')
        require(len(pairs)<=1,f'compression witness unique {n}')
        if len(ds)==6 and pairs: hits.append(n)
        if s==n+2: classes[len(ds)-1]+=1
    require(hits==[20],'five-divisor overlap')
    selected={n:divisors(n)[:-1] for n in (7,15,20,24,26,694,704,130304)}
    require(sum(selected[20])==22 and sum(selected[20][1:])==21,'20 corrections')
    require(704-694==10,'subtraction correction')
    require(len(selected[130304])==17 and sum(selected[130304])==130306,'arity boundary')
    fibers={v:[n for n in range(1,(v-1)**2+1) if sum(divisors(n))-n==v] for v in (7,15,16,20,21,22,24,26,28,36)}
    require(fibers=={7:[8],15:[16,33],16:[12,26],20:[34,361],21:[18,51,91],22:[20,38],24:[529],26:[46],28:[28],36:[24]},'complete small s fibers')
    sigma_fibers={v:[n for n in range(1,v) if sum(divisors(n))==v] for v in (24,42,60)}
    require(sigma_fibers=={24:[14,15,23],42:[20,26,41],60:[24,38,59]},'complete sigma fibers')
    return {'finite_sources':'1..10000; square-root divisor enumeration, no shape filter', 'five_proper_divisor_compression_sources':hits,'selected_proper_divisors':selected,'complete_s_fibers':fibers,'complete_sigma_fibers':sigma_fibers,'excess_two_arity_counts':dict(classes)}

def component_data(n,edges):
    seen=set(); out=[]
    adj=[set() for _ in range(n)]
    for u,v in edges: adj[u].add(v);adj[v].add(u)
    for root in range(n):
        if root in seen: continue
        todo=[root];seen.add(root);vs=[]
        while todo:
            u=todo.pop();vs.append(u)
            for v in adj[u]:
                if v not in seen: seen.add(v);todo.append(v)
        out.append((len(vs),sum(len(adj[u]) for u in vs)//2))
    return out

def graph_tests():
    counts={}
    for n in range(1,6):
        bank=list(combinations(range(n),2)); M=len(bank)
        pseudo={}; edge_counts={}; exact_layers={}
        for mask in range(1<<M):
            edges=[e for i,e in enumerate(bank) if mask>>i&1]; m=len(edges)
            data=component_data(n,edges)
            pseudo[mask]=all(e<=v for v,e in data);edge_counts[mask]=m
            t=sum(e==v-1 for v,e in data)
            best_delete=m
            sub=mask
            while True:
                if pseudo.get(sub,False): best_delete=min(best_delete,m-edge_counts[sub])
                if sub==0:break
                sub=(sub-1)&mask
            require(best_delete==m-n+t,f'edge defect n{n} mask{mask}')
            orientation_min=n
            for bits in range(1<<m):
                deg=[0]*n
                for j,(u,v) in enumerate(edges): deg[v if bits>>j&1 else u]+=1
                orientation_min=min(orientation_min,max(deg,default=0))
            density_max=0
            for W in range(1,1<<n):
                v=W.bit_count();e=sum(bool(W>>u&1 and W>>w&1) for u,w in edges)
                density_max=max(density_max,(e+v-1)//v)
            require(orientation_min==density_max,f'layer criterion n{n} mask{mask}')
            if edges:
                u,v=edges[0];subdiv=edges[1:]+[(u,n),(v,n)]
                d2=component_data(n+1,subdiv)
                require(len(subdiv)-(n+1)+sum(e==v-1 for v,e in d2)==best_delete,'subdivision defect')
        counts[n]=1<<M
    return {'universe':'every labelled simple graph on 1..5 vertices; every orientation and every induced subset', 'graph_counts':counts,'checks':['minimum deletion by all spanning edge subsets','minimum layer count by all orientations','induced density criterion','one-edge subdivision invariant']}

def main():
    print(json.dumps({'status':'PASS','assert_statements_used':False,'aliquot':aliquot(),'functional_graphs':graph_tests()},indent=2,sort_keys=True))
if __name__=='__main__':main()
