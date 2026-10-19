"""Standalone exact replay of the stopped relaxed-cone witness."""
from pathlib import Path
from itertools import permutations,combinations
from fractions import Fraction as Q
import json
R=Path(__file__).resolve().parent
s=json.loads((R/'simple_obstacle.json').read_text());p=s['pred'];n=len(p);o=s['order'];rank={v:i for i,v in enumerate(o)}
def ck(x,m):
 if not x:raise RuntimeError(m)
words=[]
for w in permutations(range(n)):
 mask=0
 for v in w:
  if p[v]&~mask:break
  mask|=1<<v
 else:words.append(w)
ck(len(words)==s['extensions'],'extension count')
I=[{v for v in range(n) if u!=v and not(p[u]>>v&1 or p[v]>>u&1)} for u in range(n)]
inc=[(a,b) if rank[a]<rank[b] else (b,a) for a,b in combinations(range(n),2) if b in I[a]]
ck(inc==list(map(tuple,s['forward_pairs'])),'orientations');d=rank[1]-rank[0];ck(d==4 and all(abs(rank[a]-rank[b])<=d for a,b in inc),'maximum span')
ck(all(len(I[a]&I[b])<=len(I[0]&I[1]) for a,b in inc if abs(rank[a]-rank[b])==d),'canonical common screen')
margins=[Q(0)]*len(inc);mass=Q(0)
def active(w,g):
 a,b=tuple(g['prefix']),tuple(g['suffix']);return w[:len(a)]==a and (not b or w[-len(b):]==b)
for g in s['mixture']:
 ww=[w for w in words if active(w,g)];ck(len(ww)==g['E'],'conditional count');q=Q(g['weight']);ck(q>0,'positive mass');mass+=q
 for i,(a,b) in enumerate(inc):margins[i]+=q*Q(3*sum(w.index(a)<w.index(b) for w in ww)-2*len(ww),len(ww))
ck(mass==1 and min(margins)>0 and list(map(str,margins))==s['margins'],'strict exact obstacle')
cells={}
for w in words:
 a,b=sorted((w.index(0),w.index(1)));key=(w[:a],w[b+1:])
 if key in [((),()),((2,),()),((),(8,)),((2,),(8,))]:
  val=sum(Q(g['weight'])/g['E'] for g in s['mixture'] if active(w,g));ck(key not in cells or cells[key][1]==val,'cell weight');cells[key]=(w,val)
ck(len(cells)==4,'rectangle realized')
det=cells[((),())][1]*cells[((2,),(8,))][1]-cells[((2,),())][1]*cells[((),(8,))][1];ck(det<0,'rank-one violation')
patterns=[]
for m in range(3,n+1):
 for S in combinations(range(n),m):
  if all(b in I[a] for a,b in combinations(S,2)):
   outside=set.intersection(*(I[v] for v in S))
   if len(outside)>=7-m:patterns.append((S,sorted(outside)))
ck(not patterns,'unexpected theorem141 pattern')
out={'status':'EXACT_RELAXED_OBSTACLE_REPLAY_PASSED','extensions':len(words),'minimum_margin':str(min(margins)),'rank_one_rectangle_determinant':str(det),'rectangle':[{'prefix':k[0],'suffix':k[1],'word':v[0],'extension_weight':str(v[1])} for k,v in cells.items()],'core_degrees':list(map(len,I)),'theorem141_core_patterns':patterns}
(R/'simple_obstacle_replay.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
