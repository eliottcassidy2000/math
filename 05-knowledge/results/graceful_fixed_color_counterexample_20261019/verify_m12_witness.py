import json,itertools
from pathlib import Path
w=json.loads((Path(__file__).parent/'witness.json').read_text());m=w['m'];F={int(d):tuple(e) for d,e in w['fixed_forest'].items()};free=w['free'];checks=0

def check(ok,msg):
 global checks
 checks+=1
 if not ok:raise RuntimeError(msg)
def parts(edges):
 g=[set() for _ in range(m+1)]
 for a,b in edges:g[a].add(b);g[b].add(a)
 unseen=set(range(m+1));out=[]
 while unseen:
  seed=min(unseen);seen={seed};todo=[seed];unseen.remove(seed)
  while todo:
   for b in g[todo.pop()]:
    if b in unseen:unseen.remove(b);seen.add(b);todo.append(b)
  out.append(sorted(seen))
 return out
check(set(F)==set(range(1,m+1))-set(free),'fixed colors exact')
check(all(b-a==d for d,(a,b) in F.items()),'fixed differences')
comps=parts(F.values());check(len(comps)==4,'fixed forest')
color_of={v:i for i,c in enumerate(comps) for v in c};contacts={};choices={}
for d in free:
 choices[d]=[(a,a+d) for a in range(m-d+1)]
 contacts[d]=sorted({tuple(sorted((color_of[a],color_of[b]))) for a,b in choices[d] if color_of[a]!=color_of[b]})
 check(len(contacts[d])==2,'two contact types')
states=[]
for es in itertools.product(*(choices[d] for d in free)):
 if len(parts(list(F.values())+list(es)))==1:states.append(es)
check(len(states)>1,'multiple completions')
adj=[{j for j,s in enumerate(states) if sum(a!=b for a,b in zip(t,s))<=2} for t in states]
unseen=set(range(len(states)));cc=[]
while unseen:
 i=min(unseen);seen={i};todo=[i];unseen.remove(i)
 while todo:
  for j in adj[todo.pop()]:
   if j in unseen:unseen.remove(j);seen.add(j);todo.append(j)
 cc.append(sorted(seen))
check(len(cc)==2,'exactly two width-two components')
print(json.dumps({'checks':checks,'candidate_free_tuples':len(list(itertools.product(*(choices[d] for d in free)))),'m':m,'free_colors':free,'fixed_forest':F,'fixed_components':comps,'contacts':contacts,'fiber_states':states,'width_two_components':cc,'component_sizes':list(map(len,cc))},indent=2,sort_keys=True))
