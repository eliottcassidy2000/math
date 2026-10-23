#!/usr/bin/env python3
from itertools import product,permutations
from collections import Counter
from pathlib import Path
import hashlib,json
from compile_attachments import compile_pairs,encode,decode,feasible,block_frame_relations
c=Counter()
def need(v,label):
    if not v:raise RuntimeError(label)
    c[label]+=1
pairs=list(product((0,1),repeat=2));edges=((0,1),(1,2),(0,2))
rels=[tuple(p for k,p in enumerate(pairs) if m>>k&1) for m in range(16)]
for R in product(rels,repeat=3):
    model=compile_pairs(((0,1),)*3,dict(zip(edges,R)))
    for x in product((0,1),repeat=3):
        ok=all((x[i],x[j]) in rr for (i,j),rr in zip(edges,R))
        z=encode(model,x)
        need(feasible(model,z)==ok,'compiler_integer_rows_match_assignments')
        if ok:need(decode(model,z)==x,'compiler_encode_decode_inverse')
# Actual tournament block comparisons: two directed C3 cells, all6 matchings.
cells=((0,1,2),(3,4,5))
frames=tuple(tuple(tuple(C[(k+a)%3] for k in range(3)) for a in range(3)) for C in cells)
def tournament(q):
    A=[[0]*6 for _ in range(6)]
    for C in cells:
        for k in range(3):A[C[k]][C[(k+1)%3]]=1
    for i in range(3):
        for j in range(3):A[i][3+j]=int(j==q[i]);A[3+j][i]=int(j!=q[i])
    return A
for q,r in product(permutations(range(3)),repeat=2):
    A=tournament(q);B=tournament(r)
    D,R=block_frame_relations(A,B,cells,frames);model=compile_pairs(D,R)
    for x in product(range(3),repeat=2):
        g=frames[0][x[0]]+frames[1][x[1]]
        direct=all(A[u][v]==B[g[u]][g[v]] for u in range(6) for v in range(6))
        need(feasible(model,encode(model,x))==direct,'actual_six_vertex_block_compiler')
try:
    block_frame_relations(tournament((0,1,2)),tournament((0,1,2)),cells,(frames[0]+frames[0][:1],frames[1]))
except ValueError:
    need(True,'duplicate_whole_frame_rejected')
else:
    raise RuntimeError('duplicate whole frame accepted')
print(json.dumps({'status':'PASS','counts':dict(c),'hashes':{p:hashlib.sha256(Path(__file__).with_name(p).read_bytes()).hexdigest() for p in ('compile_attachments.py','check_compiler.py')}},sort_keys=True,indent=2))
