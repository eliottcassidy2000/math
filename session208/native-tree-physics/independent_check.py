import itertools,json
from collections import defaultdict
from pathlib import Path
src=Path(__file__).with_name('RESULT.json')
w=json.load(open(src))['4']['energy_hostile'];pairs=list(itertools.combinations(range(4),2));idx={p:i for i,p in enumerate(pairs)};E=set()
for p in itertools.permutations(range(4)):
 forced={idx[tuple(sorted((a,b)))]:int(a<b) for a,b in zip(p,p[1:])};free=[k for k in range(6) if k not in forced]
 for bits in itertools.product((0,1),repeat=3):
  t=sum(v<<k for k,v in forced.items())+sum(v<<k for k,v in zip(free,bits));u=t^sum(1<<k for k in free);E.add(tuple(sorted((t,u))))
assert len(E)==96
T=defaultdict(list)
for a,b,c in w['tree']:
 assert (a,b) in E;T[a].append(b);T[b].append(a)
assert len(w['tree'])==63
for a,b,c in w['tree']:
 S={a};todo=[a]
 for v in todo:
  for u in T[v]:
   if {v,u}=={a,b}:continue
   if u not in S:S.add(u);todo.append(u)
 assert b not in S
 assert c==sum((v in S)!=(u in S) for v,u in E)
x=w['potential'];eg=sum((x[a]-x[b])**2 for a,b in E);et=sum(c*(x[a]-x[b])**2 for a,b,c in w['tree']);vg=sum(abs(x[a]-x[b]) for a,b in E);vt=sum(c*abs(x[a]-x[b]) for a,b,c in w['tree']);assert (eg,et,vg,vt)==(31,23,23,23)
print(json.dumps({'edges':96,'tree_edges':63,'graph_energy':eg,'tree_energy':et,'graph_TV':vg,'tree_TV':vt,'status':'independent exact check passed'},indent=2))
