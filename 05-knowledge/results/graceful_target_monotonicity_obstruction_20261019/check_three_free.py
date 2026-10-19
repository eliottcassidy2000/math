import sys,json,itertools
from collections import deque
from exact_utils import enumerate_fiber
for m in range(4,9):
 trees=enumerate_fiber(m,m,{})
 ts=[tuple(t[d] for d in range(1,m+1)) for t in trees]
 found=None;groups_checked=0
 for free in itertools.combinations(range(m),3):
  fixed=[d for d in range(m) if d not in free];buckets={}
  for t in ts:buckets.setdefault(tuple(t[d] for d in fixed),[]).append(t)
  for frozen,group in buckets.items():
   if len(group)<2:continue
   groups_checked+=1
   # With 3 free colors, width2 adjacency means at least 1 equal free coordinate.
   p=list(range(len(group)))
   def root(i):
    while p[i]!=i:i=p[i]
    return i
   by={}
   for i,t in enumerate(group):
    for d in free:
     key=(d,t[d])
     if key in by:p[root(i)]=root(by[key])
     else:by[key]=i
   classes={}
   for i in range(len(group)):classes.setdefault(root(i),[]).append(i)
   if len(classes)>1:
    reps=[group[x[0]] for x in classes.values()]
    found={'m':m,'free_colors':[d+1 for d in free],'fixed':{d+1:e for d,e in zip(fixed,frozen)},'states':len(group),'components':[len(x) for x in classes.values()],'source':reps[0],'target':reps[1]};break
  if found:break
 print(json.dumps({'m':m,'groups_checked':groups_checked,'obstruction':found}),flush=True)
 if found:break
