#!/usr/bin/env python3
"""Independent literal-path test of the all-subset/port rigidity lemma."""
from itertools import combinations,permutations
from collections import defaultdict
from pathlib import Path
import json

def check(ok,msg):
    if not ok:raise RuntimeError(msg)

def make(n,code):
    A=[[0]*n for _ in range(n)]
    for j,(u,v) in enumerate(combinations(range(n),2)):
        A[u][v]=(code>>j)&1;A[v][u]=1-A[u][v]
    K=[[[0]*n for _ in range(n)] for _ in range(1<<n)]
    for k in range(1,n+1):
        for p in permutations(range(n),k):
            if all(A[u][v] for u,v in zip(p,p[1:])):
                S=sum(1<<v for v in p);K[S][p[0]][p[-1]]+=1
    H=tuple(sum(map(sum,M)) for M in K)
    starts=[[sum(row) for row in M] for M in K]
    ends=[[sum(M[u][v] for u in range(n)) for v in range(n)] for M in K]
    return (A,K,H,starts,ends)

def main():
    result=[]
    for n in (3,4,5):
        states=[make(n,c) for c in range(1<<(n*(n-1)//2))]
        groups=defaultdict(list)
        for i,t in enumerate(states):groups[tuple(t[2][S] for S in range(1<<n) if S.bit_count()==3)].append(i)
        considered=0;tested=0;ports=0;bilinears=0;module_checks=0
        for group in groups.values():
            for i,j in combinations(group,2):
                considered+=1
                # Exhaust all through n4; deterministic sample of all H-equal n5 pairs.
                if n==5 and (i*1031+j)%17:continue
                tested+=1;A,K,H,S,E=states[i];B,L,G,U,V=states[j]
                check(H==G,'C3 profile implies all subset H in finite control')
                D=[[A[u][v]!=B[u][v] for v in range(n)] for u in range(n)]
                pending=set(range(n));comps=[]
                while pending:
                    seed=min(pending);C={seed};pending.remove(seed);stack=[seed]
                    while stack:
                        u=stack.pop()
                        new={v for v in pending if D[u][v]}
                        C|=new;pending-=new;stack.extend(new)
                    comps.append(C)
                for C in comps:
                    for u in set(range(n))-C:
                        check(len({A[u][v] for v in C})==1,'common module T')
                        check(len({B[u][v] for v in C})==1,'common module Tprime')
                        module_checks+=1
                cuts=[]
                for r in range(1<<n):
                    if all(sum((S[m][v]-U[m][v])*((r>>v)&1) for v in range(n))==0 for m in range(1<<n) if m.bit_count()==2):
                        cuts.append(r);ports+=1
                        check(all(sum((S[m][v]-U[m][v])*((r>>v)&1) for v in range(n))==0 for m in range(1<<n)),'all start entries from degree2')
                        for C in comps:check(len({(r>>v)&1 for v in C})==1,'port module')
                        check(all(sum((E[m][v]-V[m][v])*(1-((r>>v)&1)) for v in range(n))==0 for m in range(1<<n)),'incoming endpoint')
                for r in cuts:
                    for q in cuts:
                        for m in range(1<<n):
                            delta=sum((K[m][u][v]-L[m][u][v])*((r>>u)&1)*(1-((q>>v)&1)) for u in range(n) for v in range(n))
                            check(delta==0,'paired endpoint inventory')
                        bilinears+=1
        result.append(dict(n=n,tournaments=len(states),C3_equal_pairs=considered,pairs_tested=tested,
                           accepted_port_cuts=ports,ordered_port_pairs=bilinears,module_checks=module_checks,
                           scope='exhaustive pairs' if n<5 else 'deterministic 1/17 pair sample'))
    print(json.dumps({'status':'PASS','universes':result,'method':'literal permutation paths, no path DP or algebraic producer'},indent=2))
if __name__=='__main__':main()
