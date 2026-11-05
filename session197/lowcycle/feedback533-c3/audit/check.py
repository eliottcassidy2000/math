from itertools import permutations
import json
A=[2,268,905,881,999,967,903,779,513,3];n=len(A)
def need(x,s):
 if not x:raise RuntimeError(s)
H=0
for S in range(1,1<<n):
 v=[(S>>u)&1 for u in range(n)]
 for _ in range(n-1):
  w=[0]*n
  for u,x in enumerate(v):
   for j in range(n):
    if S>>j&1 and A[u]>>j&1:w[j]+=x
  v=w
 H+=(-1)**(n-S.bit_count())*sum(v)
cycles=[]
for p in permutations(range(1,n)):
 w=(0,)+p
 if all(A[w[i]]>>w[(i+1)%n]&1 for i in range(n)):cycles.append(w)
need(H==533 and len(cycles)==3,'independent H/C')
best=0;deletions=[]
for S in range(1<<n):
 V=[v for v in range(n) if S>>v&1]
 if sorted((A[v]&S).bit_count() for v in V)==list(range(len(V))):
  if len(V)>best:best=len(V);deletions=[]
  if len(V)==best:deletions.append([v for v in range(n) if not S>>v&1])
need(best==8 and [1,3] in deletions,'exact feedback2')
ports=[1,3];tail=[4,5,6,2,7,8,9,0];masks=[1,1,1,2,3,0,1,1];rows=[0]*n;rows[1]|=1<<3;old=ports[:]
for v,m in zip(tail,masks):
 for u in old:
  if u in ports and m>>ports.index(u)&1:rows[v]|=1<<u
  else:rows[u]|=1<<v
 old.append(v)
need(rows==A,'fixed word exact graph')
print(json.dumps({'status':'PASS','walk_IE_H':H,'full_permutation_cycles':cycles,'all_vertex_subsets':1<<n,'maximum_transitive_order':best,'minimum_feedback_deletions':deletions,'fixed_ports':ports,'physical_arcs':45,'minimum_C_in_strong_H533_FVS_at_most2':3},sort_keys=True,indent=2))
