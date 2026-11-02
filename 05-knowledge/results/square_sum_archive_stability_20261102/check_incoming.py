#!/usr/bin/env python3
"""Pinned incoming witness, independent Q25 incidence census, next-order hostile."""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
from math import isqrt
import subprocess,json,hashlib

BASE='d7bbbf3f05'
PATH='05-knowledge/results/machines_automata_fractran_tournaments_dictionary_20261101/data/prime_path_N617059.txt'
def need(b,m):
    if not b:raise RuntimeError(m)
def prime_witness():
    raw=subprocess.check_output(['git','show',BASE+':'+PATH])
    xs=list(map(int,raw.split()));N=617059
    sieve=bytearray(b'\1')*N;sieve[0]=sieve[1]=0
    for p in range(2,isqrt(N-1)+1):
        if sieve[p]:
            for m in range(p*p,N,p):sieve[m]=0
    ps=[p for p in range(N) if sieve[p]]
    need(len(xs)==len(set(xs))==50396 and sorted(xs)==ps,'prime universe')
    sums=[a+b for a,b in zip(xs,xs[1:])]
    need(all(isqrt(s)**2==s for s in sums),'prime square edges')
    return {'status':'verified construction','source_revision':BASE,'path':PATH,'sha256':hashlib.sha256(raw).hexdigest(),'vertices':len(xs),'max_prime':max(xs),'edges_checked':len(sums),'minimality':'inherited finite-exact nontrivial classification; not rerun here','trivial_exception':[2]}

def square25():
    n=25;edges=[(i,j) for i in range(1,n+1) for j in range(i+1,n+1) if isqrt(i+j)**2==i+j]
    adjacency={i:[] for i in range(1,n+1)}
    for i,j in edges:adjacency[i].append(j);adjacency[j].append(i)
    leaves=[v for v in adjacency if len(adjacency[v])==1]
    need(leaves==[18],'required starting leaf')
    paths=[]
    def dfs(path,used):
        if len(path)==n:paths.append(tuple(path));return
        for v in adjacency[path[-1]]:
            if not used>>v&1:dfs(path+[v],used|1<<v)
    dfs([18],1<<18)
    paths+= [p[::-1] for p in paths[:]]
    need(len(paths)==20 and len(set(paths))==20,'complete oriented paths')
    arcs=[set(zip(p,p[1:])) for p in paths]
    compatible=[]
    for i,j in combinations(range(20),2):
        if not any((b,a) in arcs[j] for a,b in arcs[i]):compatible.append((i,j))
    need(len(compatible)==4 and len(set(sum(([i,j] for i,j in compatible),[])))==8,'matching compatibility graph')
    sizes=[len(arcs[i]|arcs[j]) for i,j in compatible]
    need(sizes==[26]*4 and len(edges)==32,'union incidence')
    first=len(paths)*2**(32-24);second=sum(2**(32-s) for s in sizes)
    histogram={0:2**32-first+second,1:first-2*second,2:second}
    need(histogram=={0:4294962432,1:4608,2:256},'guarded spectrum')
    return {'status':'finite-exact','n':n,'guard_edges':len(edges),'directed_paths':len(paths),'compatible_pairs':compatible,'histogram':histogram,'outside_free_bits':268,'scope':'guarded H_G, not raw tournament H'}

def eleven_hostile():
    # Exact matching polynomial of P3 plus four independent edges.
    poly=[1,-2]
    for _ in range(4):
        nxt=[0]*(len(poly)+1)
        for j,x in enumerate(poly):nxt[j]+=x;nxt[j+1]-=x
        poly=nxt
    P=1-sum(F(a,4**j) for j,a in enumerate(poly))
    V=F(6,4);P21=F(4,64);W=F(1,16);k=F(781,1280)
    residual=P+P21-k*V
    need(residual==-F(7,640) and -residual/W==F(7,40),'eleven hostile arithmetic')
    return {'status':'refuted direct extension only','graph':'P3 disjoint-union 4 K2','P':str(P),'V':str(V),'P21':str(P21),'residual':str(residual),'minimum_W_coefficient':str(-residual/W)}

def main():
    print(json.dumps({'status':'passed','prime_witness':prime_witness(),'guarded_square25':square25(),'eleven_hostile':eleven_hostile()},sort_keys=True,indent=2))
if __name__=='__main__':main()
