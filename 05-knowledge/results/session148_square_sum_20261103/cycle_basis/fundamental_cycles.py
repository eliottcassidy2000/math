#!/usr/bin/env python3
"""Check fundamental-cycle bases on every enumerated Hamilton path."""
from collections import deque
from orientation import paths

def require(ok,msg):
 if not ok:raise ValueError(msg)

def check(n):
 E,P=paths(n); idx={e:i for i,e in enumerate(E)}; beta=len(E)-n+1
 require(P,'no Hamilton paths')
 for path in P:
  path_edges={tuple(sorted((u,v))) for u,v in zip(path,path[1:])}
  chords=set(E)-path_edges
  require(len(path_edges)==n-1 and len(chords)==beta,'chord count')
  positions={v:i for i,v in enumerate(path)}; cycles=[]
  for a,b in sorted(chords):
   lo,hi=sorted((positions[a],positions[b])); c=1<<idx[(a,b)]
   for u,v in zip(path[lo:hi],path[lo+1:hi+1]):c^=1<<idx[tuple(sorted((u,v)))]
   require(sum(bool(c>>idx[e]&1) for e in chords)==1,'chord projection')
   parity=0
   for i,(u,v) in enumerate(E):
    if c>>i&1:parity^=(1<<(u-1))|(1<<(v-1))
   require(parity==0,'cycle not Eulerian');cycles.append(c)
  generated={0}
  for c in cycles:generated|={x^c for x in generated}
  require(len(generated)==2**beta,'basis dependence')
  if n in (15,16,17):
   require(chords=={(1,3)},'unique path chord')
   cycle_edges={e for i,e in enumerate(E) if cycles[0]>>i&1}
   require(cycle_edges=={(1,3),(3,6),(6,10),(10,15),(1,15)},'5-cycle')
 print('n=%d m=%d beta=%d U=%d verified bases=%d vectors=%d path-incidences=%d'%(n,len(E),beta,len(P),len(P),2**beta,len(P)*2**(beta+1)))

for n in (15,16,17,23,25,26,27):check(n)
print('All fundamental-cycle checks passed.')
