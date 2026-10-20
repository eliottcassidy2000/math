#!/usr/bin/env python3
"""Independently construct generator actions on directed pair orbits, then BFS
all cube vertices while propagating the graph sign character. No normalizer
listing, coordinate-action deduplication, or census pair-action code is imported.
"""
import json,subprocess
from pathlib import Path
from itertools import combinations
import independent_engine as q
ROOT=Path(__file__).resolve().parent

def require(x,s):
 if not x:raise RuntimeError(s)

def run(row):
 n=row['n'];H=q.generated_group(row['generators'],n);left=set(range(n));blocks=[];fixed=[]
 while left:
  u=min(left);B=sorted({h[u] for h in H});left-=set(B)
  if len(B)>1:blocks.append(B)
  else:fixed+=B
 ngens=[]
 for B in blocks:
  p=len(B);mult=3 if p==7 else 2;g=list(range(n))
  for i,u in enumerate(B):g[u]=B[mult*i%p]
  ngens.append(g)
 if len(fixed)>1:
  g=list(range(n));g[fixed[0]],g[fixed[1]]=fixed[1],fixed[0];ngens.append(g)
  g=list(range(n))
  for i,u in enumerate(fixed):g[u]=fixed[(i+1)%len(fixed)]
  ngens.append(g)
 unseen=set(combinations(range(n),2));orbits=[];directed=[]
 while unseen:
  u,v=min(unseen);D=frozenset((h[u],h[v]) for h in H);O=frozenset(tuple(sorted(e)) for e in D)
  require(len(O)==len(D),'odd orbit orientation');orbits.append(O);directed.append(D);unseen-=O
 actions=[]
 for g in ngens:
  im=[];flip=sign=0
  for i,(O,D) in enumerate(zip(orbits,directed)):
   out=frozenset(tuple(sorted((g[u],g[v]))) for u,v in O);j=orbits.index(out);im.append(j)
   outd=frozenset((g[u],g[v]) for u,v in D)
   if outd!=directed[j]:
    require(outd==frozenset((v,u) for u,v in directed[j]),'whole orbit orientation');flip|=1<<j
   if sum(g[u]>g[v] for u,v in O)%2:sign|=1<<i
  actions.append((flip,sign,im))
 tag=f'{n}_{len(H)}_{len(orbits)}';inp=ROOT/(tag+'.generator-actions')
 with inp.open('w') as f:
  f.write(f'{len(orbits)} {len(actions)}\n')
  for flip,sign,im in actions:f.write(' '.join(map(str,[flip,sign]+im))+'\n')
 results=[]
 for T in (False,True):
  mode='T' if T else 'G';out=ROOT/(tag+'_'+mode+'.bfs-orbits')
  subprocess.run([str(ROOT/'generator_orbits'),str(inp),str(out),mode],check=True)
  got=[tuple(map(int,line.split())) for line in out.read_text().splitlines()]
  expected=[(int(s),int(size)) for s,stab,size in (line.split() for line in (ROOT/(tag+'_'+mode+'.orbits')).read_text().splitlines())]
  require(got==expected,'independent BFS exact record match')
  results.append({'kind':mode,'accepted_orbits':len(got),'accepted_labelled_objects':sum(size for s,size in got)})
 return {'class':row['class'],'pair_orbits':len(orbits),'normalizer_generators':ngens,'records':results}

def main():
 rows=[run(r) for r in json.loads((ROOT/'large13.json').read_text())['rows']]
 (ROOT/'generator_action_replay.json').write_text(json.dumps({'status':'PASS','rows':rows},indent=2)+'\n')
if __name__=='__main__':main()
