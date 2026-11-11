"""Exact new-operation controls, Python standard library only."""
from itertools import combinations, permutations
from random import Random
import json
from pathlib import Path

def tournament(n,mask):
 a=[[0]*n for _ in range(n)]
 for k,(u,v) in enumerate(combinations(range(n),2)):
  a[u][v]=(mask>>k)&1;a[v][u]=1-a[u][v]
 return a

def data(a):
 n=len(a);N=1<<n;full=N-1
 E=[[0]*n for _ in range(N)];B=[[0]*n for _ in range(N)]
 for v in range(n):E[1<<v][v]=B[1<<v][v]=1
 for S in range(1,N):
  if S&(S-1)==0:continue
  for v in range(n):
   if S>>v&1:
    R=S^(1<<v)
    E[S][v]=sum(E[R][u]*a[u][v] for u in range(n))
    B[S][v]=sum(B[R][u]*a[v][u] for u in range(n))
 H=sum(E[full]);W=[[0]*n for _ in range(n)]
 for u,v in combinations(range(n),2):
  if not a[u][v]:u,v=v,u
  z=sum(E[S][u]*B[full^S][v] for S in range(N) if S>>u&1 and not S>>v&1)
  W[u][v]=W[v][u]=z
 return H,W,B[full],E[full],B,E

def components(a):
 n=len(a);r=[[bool(i==j or a[i][j]) for j in range(n)] for i in range(n)]
 for k in range(n):
  for i in range(n):
   for j in range(n):r[i][j]|=r[i][k] and r[k][j]
 rem=set(range(n));cs=[]
 while rem:
  u=min(rem);c={v for v in rem if r[u][v] and r[v][u]};rem-=c;cs.append(c)
 allcs=cs[:]
 cs.sort(key=lambda c: -sum(a[next(iter(c))][v] for d in allcs if d!=c for v in d))
 return cs

def apex(a,x):
 n=len(a);b=[r[:] + [int(i!=x)] for i,r in enumerate(a)];b.append([int(i==x) for i in range(n)]+[0]);return b

def check(a):
 n=len(a);H,W,s,e,B,E=data(a);full=(1<<n)-1;cs=components(a)
 assert sum(map(sum,W))==2*(n-1)*H
 assert all(sum(W[v])==2*H-s[v]-e[v] for v in range(n))
 assert data([[a[j][i] for j in range(n)] for i in range(n)])[1]==W
 tight=[]
 for S in range(1,1<<(n-1)):
  capacity=sum(W[u][v] for u in range(n) if S>>u&1 for v in range(n) if not S>>v&1)
  assert capacity>=H
  if capacity==H:tight.append(S)
 expected=[];p=0
 for c in cs[:-1]:
  p|=sum(1<<v for v in c);expected.append(p if not p>>(n-1)&1 else full^p)
 assert sorted(tight)==sorted(expected)
 if len(cs)==1 and n<=6:
  hc=sum(all(a[u][v] for u,v in zip((0,)+p, p+(0,))) for p in permutations(range(1,n)))
  capacities=[sum(W[u][v] for u in range(n) if S>>u&1 for v in range(n) if not S>>v&1) for S in range(1,1<<(n-1))]
  assert min(capacities)>=H+(n-2)*hc
 # Local factors and ordered rank-one interfaces.
 local=[]
 for c in cs:
  vs=sorted(c);local.append((vs,data([[a[u][v] for v in vs] for u in vs])))
 for i,(vs,(h,w,st,en,*_)) in enumerate(local):
  for u,x in enumerate(vs):
   for v,y in enumerate(vs):assert W[x][y]*h==w[u][v]*H
  if i+1<len(local):
   vt,(h2,w2,st2,en2,*_)=local[i+1]
   for u,x in enumerate(vs):
    for v,y in enumerate(vt):assert W[x][y]*h*h2==H*en[u]*st2[v]
 # New apex arithmetic consumer: even-size pointed prefix.
 for x in range(n):
  K=sum(B[S][x]*(sum(E[full^S]) if S!=full else 1) for S in range(1,full+1) if S>>x&1 and S.bit_count()%2==0)
  assert data(apex(a,x))[0]==H+2*K
 return len(tight)

def main():
 counts={};cuts=0
 for n in range(2,6):
  counts[n]=1<<(n*(n-1)//2)
  for m in range(counts[n]):cuts+=check(tournament(n,m))
 rng=Random(265);samples={}
 for n in range(6,10):
  samples[n]=12
  for _ in range(12):check(tournament(n,rng.getrandbits(n*(n-1)//2)))
 twins=[]
 for m in (1332,2744):
  a=tournament(6,m);check(a);H,W,s,e,B,E=data(a);full=63
  terms={k:sum(B[S][2]*(sum(E[full^S]) if S!=full else 1) for S in range(1,64) if S>>2&1 and S.bit_count()==k) for k in (2,4,6)}
  # Independent permutation implementation on the explicit witness and apex.
  for b in (a,apex(a,2)):
   paths=[p for p in permutations(range(len(b))) if all(b[u][v] for u,v in zip(p,p[1:]))]
   assert len(paths)==data(b)[0]
  twins.append(dict(mask=m,adjacency=a,H=H,W=[W[u][v] for u,v in combinations(range(6),2)],start=s,end=e,even_prefix_terms=terms,apex_H=data(apex(a,2))[0]))
 assert twins[0]['W']==twins[1]['W'] and twins[0]['apex_H']!=twins[1]['apex_H']
 result=dict(status='PASS',exhaustive_labelled_tournaments=counts,exhaustive_tight_cuts=cuts,seeded_samples=samples,seed=265,twins=twins)
 target=Path(__file__).with_name('RESULTS.json');target.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
