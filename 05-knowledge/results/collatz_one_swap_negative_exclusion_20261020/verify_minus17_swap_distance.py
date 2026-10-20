#!/usr/bin/env python3
"""Exact cyclic adjacent-swap distance between two rotation orbits."""
from collections import deque,Counter
from math import comb
import json

def need(ok,message):
    if not ok:raise ArithmeticError(message)

def carry(word):
    B=0
    for i,bit in enumerate(word):B=(3 if bit else 1)*B+bit*2**i
    return B

q,p=11,7
source=tuple(((i+1)*p)//q-(i*p)//q for i in range(q))
target=tuple(map(int,'11110111000'))
D={source:0};queue=deque([source])
while queue:
    u=queue.popleft()
    for i in range(q):
        j=(i+1)%q
        if u[i]==u[j]:continue
        v=list(u);v[i],v[j]=v[j],v[i];v=tuple(v)
        if v not in D:
            D[v]=D[u]+1;queue.append(v)
need(len(D)==comb(q,p)==330,'incomplete binary-word universe')
rows=[]
for r in range(q):
    v=target[r:]+target[:r]
    s=0;prefix=[]
    for a,b in zip(source,v):s+=a-b;prefix.append(s)
    need(s==0,'unequal weights')
    median=sorted(prefix)[q//2]
    lower=sum(abs(a-median) for a in prefix)
    need(lower==D[v],'flow lower bound differs from BFS')
    rows.append({'rotation':r,'word':''.join(map(str,v)),'bfs_distance':D[v],
                 'prefix_difference':prefix,'median':median,'flow_lower_bound':lower})
path=['01011011011','10011011011','10010111011','10001111011']
edges=[(0,1),(4,5),(3,4)]
for u,v,(i,j) in zip(path,path[1:],edges):
    w=list(u);need(w[i]!=w[j],'equal-bit path edge')
    w[i],w[j]=w[j],w[i]
    need(''.join(w)==v,'bad path edge')
need(tuple(map(int,path[0]))==source,'wrong path source')
need(tuple(map(int,path[-1])) in {target[r:]+target[:r] for r in range(q)},'wrong path target')
clock=2**q-3**p
B=carry(tuple(map(int,path[-1])))
need(B==12649 and clock==-139 and B==clock*(-91),'wrong signed fixed point')
need(min(row['bfs_distance'] for row in rows)==3,'wrong minimum')
print(json.dumps({'status':'PASS','distance_definition':'minimum cyclic unequal-adjacent-swap path length between the rotation orbits of source and target',
    'source':''.join(map(str,source)),'target_reference':''.join(map(str,target)),
    'binary_word_universe':len(D),'bfs_layer_counts':dict(sorted(Counter(D.values()).items())),
    'rotation_rows':rows,'minimum_distance':3,'explicit_path':path,'explicit_edges_zero_based':edges,
    'final_word_carry':B,'clock':clock,'final_word_fixed_point':-91},indent=2))
