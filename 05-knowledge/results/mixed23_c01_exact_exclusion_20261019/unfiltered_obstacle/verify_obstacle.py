#!/usr/bin/env python3
from pathlib import Path
from fractions import Fraction
from itertools import permutations,combinations
import json
ROOT=Path(__file__).resolve().parent

def check(x,m):
 if not x:raise RuntimeError(m)
def main():
 r=json.loads((ROOT/'first_obstacle.json').read_text());p=r['pred'];n=len(p);ob=r['obstruction'];words=[]
 for o in permutations(range(n)):
  rank={v:i for i,v in enumerate(o)}
  if all(rank[u]<rank[v] for v in range(n) for u in range(n) if p[v]>>u&1):words.append(o)
 pairs=[tuple(x) for x in ob['forward_pairs']];M=ob['order'];rank={v:i for i,v in enumerate(M)};inc={tuple(x) for x in combinations(range(n),2) if not(p[x[1]]>>x[0]&1 or p[x[0]]>>x[1]&1)}
 check({tuple(sorted(x)) for x in pairs}==inc,'pair set');check(all(rank[a]<rank[b] for a,b in pairs),'orientation');span=rank[1]-rank[0];check(span>0 and all(abs(rank[a]-rank[b])<=span for a,b in inc),'span')
 margins=[Fraction(0)]*len(pairs);mass=Fraction(0)
 for g in ob['mixture']:
  u=tuple(g['prefix']);v=tuple(g['suffix']);w=Fraction(g['weight']);check(w>0,'weight');mass+=w
  selected=[o for o in words if o[:len(u)]==u and (not v or o[-len(v):]==v) and all(o.index(z)>=len(u) and o.index(z)<n-len(v) for z in (0,1))]
  E=len(selected);check(E==g['E'],'generator count')
  for j,(a,b) in enumerate(pairs):margins[j]+=w*(3*Fraction(sum(o.index(a)<o.index(b) for o in selected),E)-2)
 check(mass==1,'mass');check(margins==list(map(Fraction,ob['margins'])),'margin equality');check(min(margins)>0,'not strict')
 out={'status':'EXACT_CONE_OBSTACLE_REPLAYED','core_extensions':len(words),'projected_root_span':span,'minimum_3p_minus_2':str(min(margins)),'scope':'Core-only relaxed cone obstruction, not an actual counterexample. Root span3 is excluded globally using the known six-thin theorem.'}
 (ROOT/'obstacle_replay.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))
if __name__=='__main__':main()
