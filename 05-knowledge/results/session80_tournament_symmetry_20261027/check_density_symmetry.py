#!/usr/bin/env python3
"""C3-only maximal-prime-density family: full Aut-order distribution.
Independent cyclic-tree generator through b9; actual tournament perms through n7.
"""
from collections import Counter
from fractions import Fraction as F
from itertools import permutations,product
import json
M=15
C=[Counter({0:F(1)})]
for b in range(1,M+1):
    out=Counter()
    for i in range(b):
        for j in range(b-i):
            k=b-1-i-j
            for a,v in C[i].items():
                for d,w in C[j].items():
                    for e,q in C[k].items(): out[a+d+e]+=v*w*q/3
    if (b-1)%3==0:
        for a,v in C[(b-1)//3].items(): out[3*a]-=v/3;out[3*a+1]+=v
    out=Counter({a:v for a,v in out.items() if v})
    if any(v.denominator!=1 or v<0 for v in out.values()): raise RuntimeError('Polynomial integrality')
    if 0 in out: raise RuntimeError('Rigid nontrivial C3 tree')
    if out.get(1)!=1 or max(out)>b: raise RuntimeError('Chain or order bound')
    power=2*b+1
    while power%3==0: power//=3
    if out.get(b,0)!=(1 if power==1 else 0): raise RuntimeError('Complete-tree equality boundary')
    C.append(out)
# Independent cyclic-tree canonicalization (not weighted Burnside).
trees=[{()}]
def autexp(T): return 0 if not T else sum(map(autexp,T))+(T[0]==T[1]==T[2])
def canon(T):return min(T,T[1:]+T[:1],T[2:]+T[:2])
def tournament(T):
    if not T:return [[False]]
    children=list(map(tournament,T));sizes=list(map(len,children));offset=[0,sizes[0],sizes[0]+sizes[1]];n=sum(sizes)
    A=[[False]*n for _ in range(n)]
    for i,child in enumerate(children):
        for x in range(sizes[i]):
            for y in range(sizes[i]):A[offset[i]+x][offset[i]+y]=child[x][y]
        j=(i+1)%3
        for x in range(sizes[i]):
            for y in range(sizes[j]): A[offset[i]+x][offset[j]+y]=True
    return A
actual=[]
for b in range(1,10):
    bank=set()
    for i in range(b):
        for j in range(b-i):
            k=b-1-i-j
            for T in product(trees[i],trees[j],trees[k]):bank.add(canon(T))
    trees.append(bank)
    histogram=Counter(autexp(T) for T in bank)
    if histogram!=C[b]:raise RuntimeError(('Tree oracle',b,histogram,C[b]))
    if b<=3:
        for T in sorted(bank):
            A=tournament(T);n=len(A)
            count=sum(all(A[i][j]==A[p[i]][p[j]] for i in range(n) for j in range(i)) for p in permutations(range(n)))
            if count!=3**autexp(T):raise RuntimeError(('Actual group order',b,count,autexp(T)))
            actual.append({'prime_nodes':b,'vertices':n,'actual_aut_order':count})
# One binary linear node: independent scalar E=C^2/(1-z C^2).
counts=[sum(poly.values()) for poly in C]
square=[sum(counts[i]*counts[b-i] for i in range(b+1)) for b in range(M+1)]
E=[F(0)]*(M+1)
for b in range(M+1): E[b]=square[b]+sum(square[i]*E[b-1-i] for i in range(b))
if E[:8] != [1,3,8,23,69,208,636,1963]:
    raise RuntimeError('Independent one-linear-node tree oracle mismatch')
if any(v.denominator!=1 or v<0 for v in E):
    raise RuntimeError('One-linear-node integrality')
print(json.dumps({'checks':'PASS: exact weighted recurrence; independent cyclic-tree canonicalization b<=9; actual tournament automorphisms n<=7; no rigid nontrivial members','counts_by_b_then_log3_aut':{b:{a:int(v) for a,v in sorted(poly.items())} for b,poly in enumerate(C)},'actual_tournament_checks':actual,'one_linear_node_counts':list(map(int,E))},indent=2))
