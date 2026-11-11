#!/usr/bin/env python3
"""Exact finite controls for the all-order forest-mean and birth-sector proofs."""
from itertools import combinations
from collections import Counter
import json

def edges(n): return [(i,j) for i in range(n) for j in range(i+2,n)]
def cuts(n):
 e=edges(n)
 return sorted({sum((((s>>i)^(s>>j))&1)<<k for k,(i,j) in enumerate(e)) for s in range(1<<(n-1))})
def rows(n,t):
 r=[0]*n
 for i in range(n-1):r[i]|=1<<(i+1)
 for k,(i,j) in enumerate(edges(n)):
  u,v=(j,i) if (t>>k)&1 else (i,j)
  r[u]|=1<<v
 return r

def profile(n,t):
 r=rows(n,t); full=1<<n
 dp=[[0]*n for _ in range(full)]
 for i in range(n):dp[1<<i][i]=1
 for s in range(1,full):
  for v in range(n):
   if not dp[s][v]:continue
   out=r[v]&~s
   while out:
    z=out&-out;out-=z;w=z.bit_length()-1
    dp[s|z][w]+=dp[s][v]
 h=[sum(a) for a in dp]
 c=[[0]*(n+1) for _ in range(full)];c[0][0]=1
 for s in range(1,full):
  z=s&-s; sub=s
  while sub:
   if sub&z:
    for k in range(1,n+1):c[s][k]+=h[sub]*c[s^sub][k-1]
   sub=(sub-1)&s
 return tuple(c[-1][1:])

def birth(n,t,a):
 old=rows(n,t);N=n+1; r=[0]*N
 pos=lambda i:i+(i>=a)
 for i in range(n):
  for j in range(n):
   if old[i]>>j&1:r[pos(i)]|=1<<pos(j)
  if i==a:r[a]|=1<<pos(i)
  else:r[pos(i)]|=1<<a
 return sum(((r[j]>>i)&1)<<k for k,(i,j) in enumerate(edges(N)))

out={}
for n in range(4,7):
 cs=cuts(n);m=len(edges(n));seen=set(); sectors=[]
 for t in range(1<<m):
  if t in seen:continue
  orb={t^c for c in cs};seen.update(orb);sectors.append(orb)
 sums=[tuple(sum(profile(n,t)[k] for t in orb) for k in range(n)) for orb in sectors]
 assert len(set(sums))==1
 nc=cuts(n+1);birth_data=[];defects=[]
 for a in range(n):
  d=birth(n,(1<<m)-1,a)^birth(n,0,a)^((1<<len(edges(n+1)))-1)
  assert (d in nc)==(a==0)
  defects.append(d in nc)
  for orb in sectors:
   cnt=Counter(min(birth(n,t,a)^c for c in nc) for t in orb)
   assert len(cnt)==1<<(n-2) and set(cnt.values())=={2}
  birth_data.append(len(cnt))
 complement_same=all(min(t^c for c in cs)==min((t^((1<<m)-1))^c for c in cs) for t in range(1<<m))
 assert complement_same==(n==4)
 out[n]={'sectors':len(sectors),'orbit_size':len(cs),'common_profile_sum':sums[0],'birth_sectors_per_parent_sector_by_root':birth_data,'birth_multiplicity':2,'birth_complement_defect_is_cut_by_root':defects,'complement_in_same_sector':complement_same}
compiler={}
for n in range(1,7):
 N=n+1;e=edges(N);counts=[]
 for a in range(n):
  controls=[sum((((s>>i)^(s>>j))&1)<<k for k,(i,j) in enumerate(e)) for s in range(1<<N) if not (s>>a&1) and not (s>>(a+1)&1)]
  assert len(controls)==len(set(controls))==1<<(n-1)
  children=[birth(n,t,a)^c for t in range(1<<len(edges(n))) for c in controls]
  assert len(children)==len(set(children))==1<<len(e)
  counts.append(len(children))
 compiler[n]=counts
out['paired_root_compiler_children_by_root']=compiler
print(json.dumps(out,indent=2))
