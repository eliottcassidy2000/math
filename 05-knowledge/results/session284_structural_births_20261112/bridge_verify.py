#!/usr/bin/env python3
import itertools,math,json
rows=[]
for m in range(15):
 counts={}
 for b in itertools.product((0,1),repeat=m):
  r=sum(b[i:i+2]==(1,0) for i in range(max(0,m-1)))
  counts[r]=counts.get(r,0)+1
 for r,c in counts.items():
  if c!=math.comb(m+1,2*r+1):raise RuntimeError((m,r,c))
 rows.append({'spine_length':m,'intersection_marked':sum(c for r,c in counts.items() if r<=1)})
print(json.dumps({'status':'PASS','scope':'descent-word count control, not independent D membership census','rows':rows},indent=2))
