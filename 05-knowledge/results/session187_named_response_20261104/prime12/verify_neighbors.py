import json
from verify import endpoints,op,mods,need
bases=[[1284,1297,1322,1411,1325,1355,1311,1399,1536,1279,2048,1023],[3706,680,3603,1508,3914,1300,294,1397,519,1256,2370,1002]]
out=[]
for k,base in enumerate(bases):
 counts=[0,0,0];survivors=[]
 for ed in [None]+[(u,v) for u in range(12) for v in range(u+1,12)]:
  a=base.copy()
  if ed:
   u,v=ed;a[u]^=1<<v;a[v]^=1<<u
  counts[0]+=1;E=endpoints(a);S=endpoints(op(a));f=4095
  if any((x-y)%2 for x,y in zip(E[f],S[f])):continue
  counts[1]+=1;counts[2]+=not mods(a)
  # Independent exact two-cover obstruction, rather than path parity.
  bad=[]
  for u in range(12):
   for v in range(u+1,12):
    c=sum(E[m][u]*S[f^m][v]-E[m][v]*S[f^m][u] for m in range(1,f))
    if c%4:bad.append([u,v,c])
  need(bad,'unexpected joint survivor');survivors.append({'edge':ed,'first_cover_skew':bad[0]})
 need(counts==([67,7,6] if k==0 else [67,1,1]),'independent radius-one counts')
 out.append({'base':k,'counts':counts,'endpoint_balanced':survivors})
print(json.dumps(out,indent=2,sort_keys=True))
