#!/usr/bin/env python3
"""Verify a common base path and one full geometric-band cycle."""
import json
from pathlib import Path
from math import isqrt
ROOT=Path(__file__).resolve().parent
release=ROOT.parent/'seeded_eventual_cycle_theorem_20261019'
source=(release/'construct_seeded_cycles.py').read_text();env={'__file__':str(release/'construct_seeded_cycles.py')};exec(source[:source.index('p=argparse.ArgumentParser()')],env)
Q=49;k=2;b=54;h=(Q**k-1)//2
B=env['base_pair'](b);P=[3,46,35,29,52,48,33,16,20,44,37,12]
if any(path[:13]!=[1]+P for path in B):raise RuntimeError('common prefix')
for j,v in enumerate(P):
 if B[0].index(v)%2!=B[1].index(v)%2:raise RuntimeError('common parity')
A_index=49*b+24;N=49*A_index+24;A=env['transition'](B,'Aplus',A_index)
orders=tuple(1 if path.index(2)<path.index(7) else -1 for path in A)
table=json.loads((release.parent/'seeded_final_split_table_20261019/final_split_glue_certificates_20261019.json').read_text())
f=next(x for x in table if(x['parity'],x['residue'],x['order0'],x['order1'])==(A_index%2,24,*orders))
cuts=[min(path.index(2),path.index(7))+1 for path in A]
required=set();path_vertices=set()
for t in range(-h,h+1):
 path=[Q**k*v+(-1)**j*t for j,v in enumerate(P)]
 if any(v in path_vertices for v in path):raise RuntimeError('copy disjointness')
 path_vertices.update(path)
 required.update(tuple(sorted((x,y))) for x,y in zip(path,path[1:]))
initial_edges=len(required);seen=bytearray(N+1);first=last=None;count=seed=checks=0
for atom in f['word']:
 if atom[0]=='single':part=[atom[1]]
 else:
  _,c,i,seg,o=atom;start,end=(0,cuts[i]) if seg==0 else(cuts[i],len(A[i]))
  indices=range(start,end) if o==1 else range(end-1,start-1,-1)
  part=(49*A[i][j]+(-1)**(j+1)*c for j in indices)
 for v in part:
  checks+=1
  if not(1<=v<=N) or seen[v]:raise RuntimeError('coverage')
  seen[v]=1;count+=1
  if first is None:first=v
  if last is not None:
   checks+=1
   if isqrt(last+v)**2!=last+v:raise RuntimeError('square edge')
   required.discard(tuple(sorted((last,v))))
   if(last==3362 and v in(482,6242))or(v==3362 and last in(482,6242)):seed+=1
  last=v
required.discard(tuple(sorted((last,first))))
checks+=4
if count!=N or sum(seen)!=N or isqrt(last+first)**2!=last+first or seed!=2 or required:raise RuntimeError('final cycle/path check')
print(json.dumps({'status':'passed','checks':checks,'base_index':b,'common_prefix':[1]+P,'retained_path':P,'k':k,'N':N,'band':[Q**k*b+h,Q**k*(b+1)+h-1],'prescribed_paths':Q**k,'prescribed_vertices':len(path_vertices),'prescribed_edges':initial_edges,'seed_edges':seed,'used_transitions':env['used'],'cycle_stored':False},indent=2,sort_keys=True))
