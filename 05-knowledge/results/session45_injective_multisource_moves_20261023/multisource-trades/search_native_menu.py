#!/usr/bin/env python3
"""Exhaust the declared nonselfconverse, trivial-aut source menu for targets1965/2234."""
import argparse,csv,gzip,json,sys
from pathlib import Path
from array import array
from hashlib import sha256
from collections import Counter
from itertools import combinations
HERE=Path(__file__).resolve().parent
p=argparse.ArgumentParser();p.add_argument('--repo',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args();base=a.repo/'05-knowledge/results/metagraph_order8_20261023'
for f in json.loads((HERE/'SOURCE_PINS.json').read_text())['files']:
 if sha256((a.repo/f['path']).read_bytes()).hexdigest()!=f['sha256']:raise RuntimeError('source hash mismatch '+f['path'])
nodes=list(csv.DictReader((base/'merged_nodes.tsv').open(),delimiter='\t'));D=json.loads((a.repo/'05-knowledge/results/metagraph_order8_electrical_20261023/receipt.json').read_text())['order8'];X=D['higher_layers'][-1]['lifts'][0];Y=D['radical_lifts_mod4'];Z=D['eighth_radical_lifts_mod8'][0]
owner=array('H');owner.frombytes(gzip.decompress((base/'owner_map.u16.gz').read_bytes()))
if sys.byteorder!='little':owner.byteswap()
merged={int(row[k]):int(row['merged_id']) for row in nodes for k in ('class_a','class_b')};chords=[(i,j) for i in range(8) for j in range(i+2,8)];pairs=list(combinations(range(8),2));full=(1<<21)-1
colored={(int(r['source']),int(r['target'])):(int(r['blue']),int(r['black'])) for r in csv.DictReader(gzip.decompress((base/'merged_edges.tsv.gz').read_bytes()).decode().splitlines(),delimiter='\t')}
counts=Counter();records=[];U,V=1965,2234
for s,row in enumerate(nodes):
 if row['class_a']==row['class_b']:continue
 counts['nonselfconverse_sources']+=1;canon=int(row['canon_a']);A=[[0]*8 for _ in range(8)]
 for bit,(i,j) in enumerate(pairs):A[i][j]=canon>>(27-bit)&1;A[j][i]=1-A[i][j]
 paths=[]
 def dfs(word,left):
  if not left:paths.append(word);return
  for v in sorted(left):
   if not word or A[word[-1]][v]:dfs(word+(v,),left-{v})
 dfs((),set(range(8)))
 # The original fibre counts path words, namely actual paths divided by aut,
 # summed over the two distinct ordinary converse classes.
 if 2*len(paths)!=int(row['fibre']):continue
 counts['trivial_automorphism_sources']+=1
 if s in (U,V):continue
 counts['disjoint_support_sources']+=1;fibres={U:[],V:[]}
 for path in paths:
  q=sum(A[path[i]][path[j]]<<bit for bit,(i,j) in enumerate(chords));t=merged[owner[full^q]]
  if owner[q]!=int(row['class_a']):raise RuntimeError('actual ordinary source owner mismatch')
  if t in fibres:fibres[t].append(dict(path=path,mask=q,target=t))
 best=None;options=0
 for r in fibres[U]:
  for t in fibres[V]:
   pp,qq=r['path'],t['path']
   if pp[0]!=qq[0] or pp[-1]!=qq[-1]:continue
   for k in range(1,6):
    if pp[k:k+2]!=qq[k:k+2]:continue
    options+=1;cost=sum(x!=y for x,y in zip(pp,qq));pins=(pp[0],pp[-1],k,pp[k],pp[k+1]);key=(cost,pins,pp,qq)
    if best is None or key<best[0]:best=(key,r,t)
 if best is None:continue
 for target in (U,V):
  if colored[tuple(sorted((s,target)))]!=(0,2*len(fibres[target])):raise RuntimeError('actual paired-black convention mismatch')
 key,r,t=best;records.append(dict(source=s,targets=[U,V],source_label=[*[y[s]%2 for y in Y[:4]],X[s]%128],cost=key[0],pins=key[1],before=r,after=t,all_indexed_arc_option_records=options,complete_target_capacities={str(k):len(v) for k,v in fibres.items()}))
counts['eligible_sources']=len(records);counts['distinct_source_labels']=len({tuple(r['source_label']) for r in records})
if counts['eligible_sources']!=28 or counts['distinct_source_labels']!=28:raise RuntimeError('unexpected frozen finite menu census')
a.out.write_text(json.dumps(dict(status='EXHAUSTIVE DECLARED FINITE MENU',scope='All original nonselfconverse canonical source owners with trivial automorphism group, disjoint from targets1965/2234; endpoints plus at least one indexed internal arc; same original phase automatically because target X values equal.',counts=counts,records=records),indent=2,sort_keys=True)+'\n');print(dict(counts))
