#!/usr/bin/env python3
"""Verify a common base path and one full geometric-band cycle."""
import argparse,json
from pathlib import Path
from math import isqrt
ROOT=Path(__file__).resolve().parent
release=ROOT.parent/'seeded_eventual_cycle_theorem_20261019'
source=(release/'construct_seeded_cycles.py').read_text();env={'__file__':str(release/'construct_seeded_cycles.py')};exec(source[:source.index('p=argparse.ArgumentParser()')],env)
parser=argparse.ArgumentParser();parser.add_argument('N',nargs='?',type=int,default=6294221);N=parser.parse_args().N
if N<6294221:raise RuntimeError('N below proved threshold')
Q=49;k=2;h=1200;P=[1,3,6];A_index=(N-24)//49;d2=N-49*A_index;b=(A_index-24)//49;weak_index=(b-24)//49
weak=env['build_B'](weak_index)
t=json.loads((ROOT/'witnesses'/f'B_p{weak_index%2}_d{b-49*weak_index}.json').read_text());B=[]
for word in t['words']:
 out=[]
 for atom in word:
  if atom[0]=='single':out.append(atom[1])
  elif atom[0]=='head':out.append(49-atom[1])
  else:
   kind,c,i,o=atom;part=[49*v+(-1)**(j+1)*c for j,v in enumerate(weak[i])]
   if kind=='tail':part=part[1:]
   out.extend(part if o==1 else part[::-1])
 if out[:3]!=P:raise RuntimeError('strong prefix')
 B.append(out)
env['check_pair'](B,'B')
family='Aplus' if (A_index+int(d2>=57))%2==0 else 'Aminus';A=env['transition'](B,family,A_index)
orders=tuple(1 if path.index(2)<path.index(7) else -1 for path in A)
table=json.loads((release.parent/'seeded_final_split_table_20261019/final_split_glue_certificates_20261019.json').read_text())
f=next(x for x in table if(x['parity'],x['residue'],x['order0'],x['order1'])==(A_index%2,d2,*orders))
cuts=[min(path.index(2),path.index(7))+1 for path in A]
required=set();path_vertices=set()
for t in range(-h,h+1):
 if -24<=t<=24 or 74<=t<=122:continue
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
print(json.dumps({'status':'passed','checks':checks,'base_index':b,'retained_path':P,'k':k,'N':N,'threshold':6294221,'prescribed_paths':2303,'prescribed_vertices':len(path_vertices),'prescribed_edges':initial_edges,'seed_edges':seed,'used_transitions':env['used'],'cycle_stored':False},indent=2,sort_keys=True))
