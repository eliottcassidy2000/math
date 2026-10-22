from itertools import permutations,combinations
import json

def cycle_count(out,l):
 c=0
 for sub in combinations(range(len(out)),l):
  u=min(sub)
  for tail in permutations([v for v in sub if v!=u]):
   p=(u,)+tail
   c+=all(out[p[i]]>>p[(i+1)%l]&1 for i in range(l))
 return c
results=[]
for core in [[4,1,2],[2,4,1]]:
 out=core[:]+[0,0]
 for i,mask in enumerate([1,2]):
  for v in range(3):
   if mask>>v&1:out[v]|=1<<(3+i)
   else:out[3+i]|=1<<v
 out[3]|=1<<4
 paths=[p for p in permutations(range(5)) if all(out[p[i]]>>p[i+1]&1 for i in range(4))]
 results.append({'out_masks':out,'scores':[x.bit_count() for x in out],'H':len(paths),'c3':cycle_count(out,3),'c5':cycle_count(out,5)})
print(json.dumps(results,indent=2))
