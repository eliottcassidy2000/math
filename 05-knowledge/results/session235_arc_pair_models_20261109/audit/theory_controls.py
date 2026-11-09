#!/usr/bin/env python3
import itertools,json
from pathlib import Path

def count(a):
 n=len(a)
 def walk(v,mask):
  if mask==(1<<n)-1:return 1
  return sum(walk(w,mask|1<<w) for w in range(n) if a[v][w] and not(mask>>w&1))
 return sum(walk(v,1<<v) for v in range(n))
def strong(a):
 n=len(a)
 for s in range(n):
  seen={s}; todo=[s]
  while todo:
   v=todo.pop()
   for w in range(n):
    if a[v][w] and w not in seen:seen.add(w);todo.append(w)
  if len(seen)!=n:return False
 return True
checks=[]
for h in range(3,9):
 for m in range(3,8):
  for mask in range(1,(1<<h)-1):
   # h core-cycle vertices, hub h, m-1 pendant vertices.
   n=h+m;a=[[0]*n for _ in range(n)]
   for i in range(h):a[i][(i+1)%h]=1
   for i in range(h):
    if mask>>i&1:a[i][h]=1
    else:a[h][i]=1
   for i in range(h,n-1):a[i][i+1]=1
   a[n-1][h]=1
   assert strong(a) and count(a)==h,(h,m,mask,count(a))
   checks.append((h,m,mask))
bouquets=0
for r in range(3,9):
 for s in range(3,9):
  n=r+s-1;a=[[0]*n for _ in range(n)]
  for c in ([0]+list(range(1,r)),[0]+list(range(r,n))):
   for u,v in zip(c,c[1:]+c[:1]):a[u][v]=1
  assert strong(a) and count(a)==2
  bouquets+=1
parity=0
for n in range(1,5):
 arcs=[(i,j) for i in range(n) for j in range(n) if i!=j]
 for mask in range(1<<len(arcs)):
  a=[[0]*n for _ in range(n)]
  for k,(i,j) in enumerate(arcs):a[i][j]=mask>>k&1
  b=[[int(i!=j and not a[j][i]) for j in range(n)] for i in range(n)]
  assert (count(a)-count(b))%2==0
  parity+=1
report={'cycle_hub_pendant_cases':len(checks),'h_range':[3,8],'pendant_cycle_length_range':[3,7],'all_nontrivial_partitions':True,'two_cycle_bouquets':bouquets,'parity_all_digraphs_through_order4':parity,'status':'PASS'}
p=Path(__file__).with_suffix('.json');p.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
