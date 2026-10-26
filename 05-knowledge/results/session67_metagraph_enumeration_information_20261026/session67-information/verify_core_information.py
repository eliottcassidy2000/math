#!/usr/bin/env python3
"""Independent n<=5 information audit of maximal-transitive-module encoding."""
from itertools import combinations,permutations
from collections import defaultdict
import json,math
from pathlib import Path
from verify_information import graph,require
ROOT=Path(__file__).resolve().parent

def core_code(A):
    n=len(A)
    if not n:return (0,(),())
    modules=[]
    for mask in range(1,1<<n):
        ins=[i for i in range(n) if (mask>>i)&1];out=[i for i in range(n) if not (mask>>i)&1]
        if any(len({A[v][u] for u in ins})>1 for v in out):continue
        if any(A[i][j]*A[j][k]*A[k][i]+A[i][k]*A[k][j]*A[j][i] for i,j,k in combinations(ins,3)):continue
        modules.append(mask)
    maximal=[m for m in modules if not any(m!=r and m&r==m for r in modules)]
    require(sum(m.bit_count() for m in maximal)==n,'maximal transitive modules overlap')
    require(sum(maximal)==(1<<n)-1,'maximal modules not a partition')
    blocks=[[i for i in range(n) if m>>i&1] for m in maximal];q=len(blocks)
    Q=[[0 if i==j else A[blocks[i][0]][blocks[j][0]] for j in range(q)] for i in range(q)]
    for i,j in combinations(range(q),2):require(any(Q[v][i]!=Q[v][j] for v in range(q) if v not in [i,j]),'quotient has transitive pair module')
    possibilities=[]
    for p in permutations(range(q)):
        bits=tuple(Q[p[i]][p[j]] for i,j in combinations(range(q),2))
        weights=tuple(len(blocks[p[i]]) for i in range(q))
        possibilities.append((q,bits,weights))
    return min(possibilities)

def canonical(A):
    n=len(A)
    return min(tuple(A[p[i]][p[j]] for i,j in combinations(range(n),2)) for p in permutations(range(n)))

out=[]
for n in range(6):
    d=defaultdict(set);e=defaultdict(set);sizes=defaultdict(int)
    for b in range(1<<math.comb(n,2)):
        A=graph(n,b);c=canonical(A);k=core_code(A);d[k].add(c);e[c].add(k)
    require(all(len(x)==1 for x in d.values()),'weighted core loses class')
    require(all(len(x)==1 for x in e.values()),'weighted core not intrinsic')
    for k in d:sizes[k[0]]+=1
    out.append({'n':n,'labelled_tournaments':1<<math.comb(n,2),'weighted_core_orbits':len(d),'conditional_class_entropy_given_weighted_core_exact':0,'classes_by_core_order':dict(sorted(sizes.items()))})
(ROOT/'core_information_verification.json').write_text(json.dumps({'status':'PASS','scope':'FINITE-EXACT n=0..5; does not replace all-order maximal-module theorem','rows':out},indent=2)+'\n')
print(json.dumps({'status':'PASS','rows':out},indent=2))
