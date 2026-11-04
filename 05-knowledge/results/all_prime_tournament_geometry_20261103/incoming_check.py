#!/usr/bin/env python3
"""Literal controls linking the incoming norm wrapper to AP3's H11 seed."""
import itertools as it
import json
from check import need, cyclic, profile, slot, coupled, path_counts, val, ordinal


def mirror(rows):
    n=len(rows)
    rev=tuple(sum(1<<j for j in range(n) if rows[j]>>i&1) for i in range(n))
    graph=coupled([rows,rev])
    return tuple(row & ((1<<(2*n+1))-1) for row in graph[:-1])


records=[]
for rows in [(0,),(2,0),cyclic(3),(8,1,3,6),(24,1,3,6,14),cyclic(5)]:
    P=profile(rows);a,b=P[1],slot(P,2)
    graph=mirror(rows)
    H=path_counts(graph)[-1]
    need(H==(b-a)**2+2*a*a and H%8==3,'literal incoming norm')
    records.append({'rows':rows,'order_wrapper':len(graph),'H':H,'norm_coordinates':[b-a,a]})
B=(24,1,3,6,14);target=mirror((2,0))
perms=[p for p in it.permutations(range(5)) if all((B[i]>>j&1)==(target[p[i]]>>p[j]&1) for i in range(5) for j in range(5))]
need(perms and profile(B)==profile(target),'H11 exact mirrored source')
controls=0
for p in (5,7,13,23,29,31):
    need(pow(p-2,(p-1)//2,p)==p-1,'inert prime control')
    for x,y in it.product(range(-30,31),repeat=2):
        if not (x or y):continue
        need(val(x*x+2*y*y,p)==2*min(val(x,p),val(y,p)),'norm valuation')
        controls+=1
T=ordinal((8,1,3,6),B);P=profile(T);a,b=P[1],slot(P,2)
norm=(b-a)**2+2*a*a
need((a,b,norm,val(norm,5))==(55,3200,9897075,2),'shared slot norm bridge')
print(json.dumps({'status':'FINITE-EXACT inherited-norm integration',
  'literal_wrappers':records,'H11_isomorphism_to_mirrored_TT2':perms[0],
  'H55_mirrored_norm':{'a':a,'b':b,'H':norm,'v5':val(norm,5),'count_method':'inherited universal norm formula'},
  'inert_norm_controls':controls},indent=2,sort_keys=True))
