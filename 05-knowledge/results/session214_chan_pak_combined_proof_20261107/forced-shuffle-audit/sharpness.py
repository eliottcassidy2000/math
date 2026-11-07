#!/usr/bin/env python3
import json

def count(rows):
 n=len(rows);d=[[0]*n for _ in range(1<<n)]
 for v in range(n):d[1<<v][v]=1
 for mask in range(1,1<<n):
  for last in range(n):
   z=d[mask][last]
   if not z:continue
   todo=rows[last]&~mask
   while todo:
    bit=todo&-todo;todo-=bit
    d[mask|bit][bit.bit_length()-1]+=z
 return sum(d[-1]),d[-1]

cases=[]
for n in range(3,15):
 a=[sum(1<<j for j in range(i+1,n)) for i in range(n)]
 a[0]^=1<<(n-1);a[-1]|=1
 h,e=count(a);a0=n-2
 child=[r|((1<<n) if i!=a0 else 0) for i,r in enumerate(a)]+[1<<a0]
 hy,_=count(child)
 if h!=2**(n-2)+1 or e[a0]!=2**(n-3) or hy!=(5 if n==3 else 3*2**(n-3)+3):raise RuntimeError((n,h,e,hy))
 cases.append({'n':n,'parent':h,'penultimate_end_count':e[a0],'child':hy})
print(json.dumps({'status':'PASS','sharpness_family':cases},indent=2,sort_keys=True))
