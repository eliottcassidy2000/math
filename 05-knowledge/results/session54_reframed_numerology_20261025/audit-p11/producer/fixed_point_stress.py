#!/usr/bin/env python3
import itertools as it,json,pathlib
D=pathlib.Path(__file__).resolve().parent
n=11;Q={x*x%n for x in range(1,n)};E=[(u,v) for u in range(n) for v in range(n) if (v-u)%n in Q];EI={e:i for i,e in enumerate(E)}
def ck(p,m):
 if not p:raise RuntimeError(m)
def act(B):return sum(1<<EI[4*u%n,4*v%n] for j,(u,v) in enumerate(E) if B>>j&1)
def path_dp(B):
 inc=[0]*n
 for i,(u,v) in enumerate(E):
  if not(B>>i&1):inc[v]|=1<<u
 end=[0]*(1<<n)
 for v in range(n):end[1<<v]=1<<v
 for S in range(1,1<<n):
  if S&(S-1)==0:continue
  for v in range(n):
   if S>>v&1 and end[S^(1<<v)]&inc[v]:end[S]|=1<<v
 if not end[-1]:return None
 S=(1<<n)-1;v=(end[-1]&-end[-1]).bit_length()-1;path=[v]
 while S&(S-1):
  S^=1<<v;p=end[S]&inc[v];v=(p&-p).bit_length()-1;path.append(v)
 return path[::-1]
def path_dfs(B):
 out={u:[v for v in range(n) if (u,v) in EI and not(B>>EI[u,v]&1)] for u in range(n)}
 def rec(P,S):
  if len(P)==n:return P
  for v in out[P[-1]]:
   if not(S>>v&1):
    ans=rec(P+[v],S|(1<<v))
    if ans:return ans
 for v in range(n):
  ans=rec([v],1<<v)
  if ans:return ans
 return None
left=set(range(55));orbs=[]
while left:
 j=min(left);O=set();u,v=E[j]
 for _ in range(5):O.add(EI[u,v]);u=4*u%n;v=4*v%n
 ck(len(O)==5,'orbit size');left-=O;orbs.append(sum(1<<i for i in O))
ck(len(orbs)==11,'arc orbit count')
five_rows=[]
for B in orbs:
 p=path_dp(B);q=path_dfs(B);ck(bool(p)==bool(q) and bool(p),'invariant five-arc hostile')
 five_rows.append({'mask':B,'hamiltonian_path':p})
rows=[];blocked=[]
for a,b in it.combinations(range(11),2):
 B=orbs[a]|orbs[b];ck(B.bit_count()==10 and act(B)==B,'invariant mask')
 p=path_dp(B);q=path_dfs(B);ck(bool(p)==bool(q),'DP versus DFS')
 if p:
  ck(len(set(p))==11 and all((u,v) in EI and not(B>>EI[u,v]&1) for u,v in zip(p,p[1:])),'path witness')
 else:blocked.append(B)
 rows.append({'arc_orbits':[a,b],'mask':B,'hamiltonian_path':p})
star=sum(1<<i for i,(u,v) in enumerate(E) if 0 in (u,v))
print(json.dumps({'status':'PASS','invariant_arc_orbits':orbs,'invariant_10_arc_sets':len(rows),'fixed_blockers':blocked,'star0':star,'field_F1024_nonzero_fixed_under_fourth_power':3,'rows':rows,'invariant_five_arc_survivors':five_rows},sort_keys=True,indent=2))
