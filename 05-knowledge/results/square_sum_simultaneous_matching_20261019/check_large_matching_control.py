#!/usr/bin/env python3
"""Stream one6.47m-vertex cycle, checking the level3 matching without storing it."""
import json
from pathlib import Path
from math import isqrt
ROOT=Path(__file__).resolve().parent
release=ROOT.parent/'seeded_eventual_cycle_theorem_20261019'
source=(release/'construct_seeded_cycles.py').read_text();env={'__file__':str(release/'construct_seeded_cycles.py')};exec(source[:source.index('p=argparse.ArgumentParser()')],env)
Q=49;k=3;power=Q**k;threshold=(107*power-1)//2
B_index=2694;A_index=49*B_index+24;N=49*A_index+24
B=env['build_B'](B_index);A=env['transition'](B,'Aplus',A_index)
orders=tuple(1 if path.index(2)<path.index(7) else -1 for path in A)
table=json.loads((release.parent/'seeded_final_split_table_20261019/final_split_glue_certificates_20261019.json').read_text())
f=next(x for x in table if(x['parity'],x['residue'],x['order0'],x['order1'])==(A_index%2,24,*orders))
cuts=[min(path.index(2),path.index(7))+1 for path in A]
lo=(power+1)//2;hi=(3*power-1)//2;gap=(Q**(k-1)-1)//2
seen=bytearray(N+1);found=bytearray(hi-lo+1);first=None;last=None;count=0;seed=0;checks=0
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
   a=min(last,v)
   if last+v==4*power and lo<=a<=hi and abs(a-power)>gap:found[a-lo]=1
   if(last==3362 and v in(482,6242))or(v==3362 and last in(482,6242)):seed+=1
  last=v
checks+=4
if count!=N or sum(seen)!=N or isqrt(last+first)**2!=last+first or seed!=2:raise RuntimeError('final cyclic check')
required=power-Q**(k-1)
if sum(found)!=required:raise RuntimeError(('missing matching',sum(found),required))
print(json.dumps({'status':'passed','k':k,'threshold':threshold,'N':N,'checks':checks,'matching_edges_checked':required,'matching_lower_vertex_range':[lo,hi],'excluded_central_interval':[power-gap,power+gap],'square_sum':4*power,'square_root':2*7**k,'seed_edges_preserved':seed,'B_index':B_index,'A_index':A_index,'used_transitions':env['used'],'large_cycle_stored':False},indent=2,sort_keys=True))
