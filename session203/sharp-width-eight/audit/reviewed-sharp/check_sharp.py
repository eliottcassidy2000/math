#!/usr/bin/env python3
from pathlib import Path
from math import gcd
import json
P=Path(__file__).resolve().parent

def need(q,why):
 if not q:raise RuntimeError(why)

def profile(a):
 n=len(a);N=1<<n;F=N-1
 for i,x in enumerate(a):
  need(0<=x<N and not x>>i&1,'row')
  for j in range(i):need(((x>>j&1)+(a[j]>>i&1))==1,'orientation')
 end=[[0]*n for _ in range(N)];start=[[0]*n for _ in range(N)];H=[0]*N;H[0]=1
 for v in range(n):end[1<<v][v]=start[1<<v][v]=1
 for S in range(1,N):
  for v in range(n):
   if not S>>v&1:continue
   A=S^(1<<v)
   for w in range(n):
    if not A>>w&1:continue
    if a[w]>>v&1:end[S][v]+=end[A][w]
    else:start[S][v]+=start[A][w]
  H[S]=sum(end[S])
 L=sum(end[S][0]*H[F^S]for S in range(N))
 return H[F],start[F][0],L

def child(u):
 out=[2,8,3,5]+[0]*(len(u)-1)
 for j in range(1,len(u)):
  x=j+3
  for i in range(4):
   if i or u[0]>>j&1:out[i]|=1<<x
   else:out[x]|=1<<i
  for k in range(1,j):
   y=k+3
   if u[j]>>k&1:out[x]|=1<<y
   else:out[y]|=1<<x
 return out

primary={tuple(map(int,l.split()[:2]))for l in(P/'run153/STRONG_PROFILES.tsv').read_text().splitlines()}
secondary={tuple(map(int,l.split()))for l in(P/'rooted153/ROOTED_PROFILES.tsv').read_text().splitlines()}
need(primary==secondary and len(primary)==535,'two complete profile sets')
requests=[]
for j in range(8):
 B,s=153-2*j,1+3*j
 need(gcd(B,s)==1 and B>s,'weak reduction')
 need((B,s)not in primary,'absent strong source')
 requests.append({'j':j,'B':B,'s':s,'gcd':gcd(B,s),'status':'ABSENT_AT_ALL_SOURCE_ORDERS'})
u=[222,128,2,6,14,31,62,124];p=profile(u);need(p==(77,25,85),'repair profile');need((137,25)in primary,'repair occurrence');c=child(u);need(profile(c)[0]==461,'actual four-core child')
print(json.dumps({'status':'PASS','counterexample':{'condition':'F(7,137)','base_X':153,'residue_r':1,'absent_requests':requests},'repair':{'j':8,'B':137,'s':25,'root':0,'rows':u,'h_s_L':p,'child_rows':c,'child_H':461},'complete_strong_profile_count':535,'F_8_137':'OPEN'},sort_keys=True,indent=2))
