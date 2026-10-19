import itertools,sys,json
from exact_utils import partition
E=list(itertools.combinations(range(4),2));support=[]
for mask in range(1,64):
 es=[E[i] for i in range(6) if mask>>i&1]
 if 4-len(partition(es,3))<=2:support.append(tuple(es))
counts={};bad=[]
for sets in itertools.product(support,repeat=3):
 states=[x for x in itertools.product(*sets) if len(partition(x,3))==1]
 if not states:continue
 seen={0};stack=[0]
 while stack:
  i=stack.pop()
  for j,t in enumerate(states):
   if j not in seen and any(a==b for a,b in zip(states[i],t)):seen.add(j);stack.append(j)
 if len(seen)<len(states):
  union=set().union(*map(set,sets));expected=len(union)==3 and len(partition(union,3))==1 and all(len(s)==2 for s in sets) and len(set(sets))==3
  if not expected:bad.append((sets,states))
  key=(len(union),tuple(sorted(map(len,sets))),len(states));counts[str(key)]=counts.get(str(key),0)+1
if bad: raise RuntimeError('unexpected disconnected support triple')
if len(support)!=25 or sum(counts.values())!=96: raise RuntimeError('classification count mismatch')
print(json.dumps({'supports':len(support),'ordered_support_triples':len(support)**3,'disconnected_by_type':counts,'unexpected':len(bad)},sort_keys=True,indent=2))
