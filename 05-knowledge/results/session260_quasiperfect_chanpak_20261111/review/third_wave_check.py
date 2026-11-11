#!/usr/bin/env python3
"""Focused independent controls, not a replay of the edge-shell censuses."""
from itertools import product,permutations
from pathlib import Path
import json

def check(q,label):
 if not q:raise RuntimeError(label)
def dom(w):
 n=len(w);a=[set() for _ in w]
 for i in range(n):
  for j in range(i+1,n):
   if w[i]=='A' and w[j]=='B':a[i].add(j);a[j].add(i)
 return sum(all(i in S or a[i]&S for i in range(n)) for bits in product((0,1),repeat=n) for S in [{i for i,b in enumerate(bits) if b}])
def count(adj):
 n=len(adj);d={}
 for i in range(n):d[1<<i,i]=1
 for mask in range(1,1<<n):
  for last in range(n):
   z=d.get((mask,last),0)
   if z:
    for v in adj[last]:
     if not mask>>v&1:d[mask|1<<v,v]=d.get((mask|1<<v,v),0)+z
 return sum(d.get(((1<<n)-1,i),0) for i in range(n))
core=[]
for n in range(2,7):
 for middle in product('AB',repeat=n-2):
  w='A'+''.join(middle)+'B';D=dom(w)
  F=1+sum(2 if j-i==1 else 2**(j-i-1) for i in range(n) for j in range(i+1,n) if w[i]=='A' and w[j]=='B')
  check(D==F,('adapter',w));core.append(D)
  if n>=3:
   P=sum(2**(j-1) for j in range(1,n-1) if w[j]=='B')
   Q=sum(2**(n-i-2) for i in range(1,n-1) if w[i]=='A')
   check(D-dom(w[1:-1])==2**(n-2)+P+Q+(w[1]=='B')+(w[-2]=='A'),('boundary',w))
check(29 not in core and 13 not in core,'core exclusions')
levels=[]
for bits in product((0,1),repeat=5):
 adj=[{j for j in range(i+1,5)} for i in range(5)]+[set()]
 for i,b in enumerate(bits):
  if b:adj[5].add(i)
  else:adj[i].add(5)
 # Direct permutations, not the OCF or bit recurrence.
 H=sum(all(p[i+1] in adj[p[i]] for i in range(5)) for p in permutations(range(6)))
 levels.append(H)
check(sorted(set(levels),reverse=True)[:2]==[23,19],'apex top gap a2')
abstract=[{1,2,3,4,5,6},{2,3,4,6},{4,6},{2,4},{5},{1,2,3,6},{3,4}]
cases=[]
for a in range(3,7):
 n=2*a+2;r=n-1;adj=[set() for _ in range(n)]
 for i in range(n):
  for j in range(i+1,n):
   if j==r and i<a:adj[j].add(i)
   else:adj[i].add(j)
 for i,j in [(0,r),(a-2,r),(a-1,2*a),(a-1,r)]:
  if j in adj[i]:adj[i].remove(j);adj[j].add(i)
  else:adj[j].remove(i);adj[i].add(j)
 groups=[[0],[a-2],[a-1],[2*a],[r],list(range(1,a-2)),list(range(a,2*a))]
 for i,group in enumerate(groups):
  for j,other in enumerate(groups):
   if i==j:continue
   check(all((y in adj[x])==(j in abstract[i]) for x in group for y in other),('rank interface',a,i,j))
 if a<=5:
  H=count(adj)
  F=9 if a==3 else 7-45*2**(a-3)-4*3**a+7*2**(2*a-3)+2**(a-1)*3**a
  check(H==F,('rank formula',a));cases.append([a,H])
out={'status':'PASS','chain_cores':len(core),'all_apex_cuts_literal_permutations':32,'rank_interfaces':[3,4,5,6],'rank_count_controls':cases,'core_values':sorted(set(core))}
Path(__file__).with_name('third_wave_receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
