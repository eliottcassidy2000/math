#!/usr/bin/env python3
"""Exact all-order theorem controls; no producer code imports."""
from itertools import combinations,permutations
from pathlib import Path
import json
HERE=Path(__file__).resolve().parent

def need(ok,m):
 if not ok:raise RuntimeError(m)

def all_ends(rows):
 n=len(rows);full=(1<<n)-1;E=[[0]*n for _ in range(full+1)]
 incoming=[sum(1<<u for u in range(n) if rows[u]>>v&1) for v in range(n)]
 for S in range(1,full+1):
  vs=S
  while vs:
   b=vs&-vs;vs-=b;v=b.bit_length()-1;old=S^b
   if not old:E[S][v]=1;continue
   us=old&incoming[v];total=0
   while us:
    a=us&-us;us-=a;total+=E[old][a.bit_length()-1]
   E[S][v]=total
 return E

def strong(rows,S):
 for x in range(len(rows)):
  if not(S>>x&1):continue
  seen=1<<x
  while True:
   nxt=seen;vs=seen
   while vs:
    b=vs&-vs;vs-=b;nxt|=rows[b.bit_length()-1]&S
   if nxt==seen:break
   seen=nxt
  if seen!=S:return False
 return True

def tournaments(n):
 pairs=list(combinations(range(n),2))
 for mask in range(1<<len(pairs)):
  rows=[0]*n
  for k,(i,j) in enumerate(pairs):
   if mask>>k&1:rows[i]|=1<<j
   else:rows[j]|=1<<i
  yield rows

def literal_paths(rows,S):
 return [p for p in permutations([v for v in range(len(rows)) if S>>v&1]) if all(rows[u]>>v&1 for u,v in zip(p,p[1:]))]

def extend_fixed_endpoint(rows,S,p):
 full=(1<<len(rows))-1;p=list(p);endpoint=p[-1]
 while S!=full:
  choices=[v for v in range(len(rows)) if not(S>>v&1) and rows[v]&S]
  need(choices,'incoming boundary exists')
  z=choices[0]
  for j,v in enumerate(p):
   if rows[z]>>v&1 and (j==0 or rows[p[j-1]]>>z&1):
    p.insert(j,z);break
  else:raise RuntimeError('no endpoint-preserving insertion')
  S|=1<<z
 need(p[-1]==endpoint and all(rows[u]>>v&1 for u,v in zip(p,p[1:])),'valid extended path')
 return tuple(p)

def bridge_from_set(rows,S,x):
 # Multi-source BFS fixes the bridge independently of every Hamilton word.
 queue=[v for v in range(len(rows)) if S>>v&1]
 parent={v:None for v in queue}
 for u in queue:
  if x in parent:break
  for v in range(len(rows)):
   if rows[u]>>v&1 and v not in parent:
    parent[v]=u;queue.append(v)
 need(x in parent,'bridge exists')
 path=[x]
 while parent[path[-1]] is not None:path.append(parent[path[-1]])
 path.reverse()
 need(S>>path[0]&1 and all(not(S>>v&1) for v in path[1:]),'bridge avoids other S vertices')
 return path

counts={'exterior_bridge_injections':0,'labelled_tournaments' :0,'fixed_root_subset_inequalities':0,'strong_subset_endpoint_inequalities':0,'nontrivial_floor_inequalities':0,'literal_injections':0}
for n in range(1,7):
 full=(1<<n)-1
 for rows in tournaments(n):
  E=all_ends(rows);ends=E[full];counts['labelled_tournaments']+=1
  for S in range(1,full+1):
   for x in range(n):
    if S>>x&1 and ends[x]:
     need(E[S][x]<=ends[x],'fixed-root hereditary inequality')
     counts['fixed_root_subset_inequalities']+=1
   if strong(rows,S):
    q=min(E[S][x] for x in range(n) if S>>x&1)
    need(q>=1,'strong endpoint positivity')
    for e in ends:
     need(e==0 or e>=q,'strong-subtournament floor')
     counts['strong_subset_endpoint_inequalities']+=1
     if q>1 and S!=full:counts['nontrivial_floor_inequalities']+=1
   if n<=4:
    P=literal_paths(rows,S)
    for x in range(n):
     if not ends[x]:continue
     sources=[p for p in P if p[-1]==x]
     targets=[extend_fixed_endpoint(rows,S,p) for p in sources]
     need(len(targets)==len(set(targets)),'injection')
     need(all(tuple(v for v in p if S>>v&1) in sources for p in targets),'deletion inverse')
     counts['literal_injections']+=len(sources)
    if strong(rows,S):
     for x in range(n):
      if S>>x&1 or not ends[x]:continue
      bridge=bridge_from_set(rows,S,x);a=bridge[0]
      sources=[p for p in P if p[-1]==a]
      W=S|sum(1<<v for v in bridge[1:])
      targets=[extend_fixed_endpoint(rows,W,p+tuple(bridge[1:])) for p in sources]
      need(len(targets)==len(set(targets)),'bridge injection')
      need(all(tuple(v for v in p if S>>v&1) in sources for p in targets),'bridge deletion inverse')
      counts['exterior_bridge_injections']+=len(sources)
S=[10,12,9,0];T=[10,28,25,16,1]
need(all_ends(S)[-1]==[0,0,0,3],'nonstrong source')
need(all_ends(T)[-1][0]==1 and strong(T,31),'nonstrong counterexample')
need([p for p in literal_paths(T,31) if p[-1]==0]==[(1,2,3,4,0)],'unique witness')
B=[2,4,97,39,15,83,27]
def birth(B,mask):
 n=len(B);return [r|(0 if mask>>i&1 else 1<<n) for i,r in enumerate(B)]+[mask]
def reverse(B):
 n=len(B);return [((1<<n)-1)^(1<<i)^r for i,r in enumerate(B)]
banks={}
for d,C in [(0,B),(2,birth(B,59)),(4,birth(B,27))]:
 e=all_ends(C)[-1];s=all_ends(reverse(C))[-1]
 need(strong(C,(1<<len(C))-1),'strong special bank')
 need(min(e)==min(s)==(3 if d==0 else 7),'special bank floors')
 banks[d]={'rows':C,'ends':e,'starts':s,'endpoint_floor':min(e),'start_floor':min(s)}
# Exact join strongness criterion and retained five-core values.
def join(A,a,C,b):
 n=len(A);m=len(C)
 J=[r|(((1<<m)-1)<<n) for r in A]+[r<<n for r in reverse(C)]
 J[a]&=~(1<<(n+b));J[n+b]|=1<<a
 return J

def reaches_all(rows,a):
 seen=1<<a
 while True:
  nxt=seen
  for u in range(len(rows)):
   if seen>>u&1:nxt|=rows[u]
  if nxt==seen:return seen==(1<<len(rows))-1
  seen=nxt

counts['join_strongness_checks']=0
for n in range(1,4):
 for A in tournaments(n):
  for a in range(n):
   for C in [B,birth(B,59),birth(B,27)]:
    J=join(A,a,C,1)
    need(strong(J,(1<<len(J))-1)==reaches_all(A,a),'join strongness iff')
    counts['join_strongness_checks']+=1
counts['special_exterior_bridge_injections']=0
for C in [B,birth(B,59),birth(B,27)]:
 n=len(C);bank_mask=(1<<n)-1;x=n
 ambient=[r|(1<<n) for r in C]+[0,(1<<(n+1))-1]
 bridge=bridge_from_set(ambient,bank_mask,x);a=bridge[0]
 sources=[p for p in literal_paths(C,bank_mask) if p[-1]==a]
 targets=[extend_fixed_endpoint(ambient,bank_mask|(1<<x),p+(x,)) for p in sources]
 need(len(targets)==len(set(targets))>=min(all_ends(C)[-1]),'nontrivial-floor exterior bridge injection')
 need(all(tuple(v for v in p if v<n) in sources for p in targets),'special bridge deletion inverse')
 counts['special_exterior_bridge_injections']+=len(sources)
U=[24,9,3,20,6]
need(all_ends(U)[-1]==all_ends(reverse(U))[-1]==[3]*5,'regular five obstruction')
source_extension=[r<<1 for r in B]
source_extension=[((1<<(len(B)+1))-1)^1]+source_extension
need(all_ends(source_extension)[-1]==[0]+all_ends(B)[-1],'proper extension sharpness')
sink_extension=[r|(1<<len(B)) for r in B]+[0]
need(all_ends(sink_extension)[-1]==[0]*len(B)+[57],'zero endpoint caveat')
need(counts['labelled_tournaments']==33867,'exhaustive universe')
report={'status':'PASS','counts':counts,'universe':'All 33867 labelled tournaments of orders1 through6; all nonempty induced subsets; all final endpoint choices. Explicit path injections through order4.','nonstrong_counterexample':{'source_rows':S,'source_ends':all_ends(S)[-1],'target_rows':T,'target_ends':all_ends(T)[-1],'unique_path_ending_0':[1,2,3,4,0]},'banks':banks}
(HERE/'receipt.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,sort_keys=True))
