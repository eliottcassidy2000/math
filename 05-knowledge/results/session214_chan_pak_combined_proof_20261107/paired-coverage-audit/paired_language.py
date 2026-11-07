#!/usr/bin/env python3
"""Independent directed-word profiles plus integer matrix language checks."""
import json
from pathlib import Path

def need(x,why):
 if not x:raise RuntimeError(why)

# Coordinates h,s,e,L,R,d,C,b,T; matrices read from inherited nine-coordinate law.
O=((1,0,0,0,0,0,0,0,0),(0,1,0,0,0,0,0,0,0),(0,0,0,0,0,0,0,0,0),(0,0,0,1,0,0,0,0,0),(0,0,0,0,2,0,0,0,0),(0,0,0,0,0,1,0,0,0),(0,0,0,0,0,0,2,0,0),(0,0,0,0,0,0,0,2,0),(0,0,0,0,0,0,0,0,3))
Z=((1,0,-1,0,1,0,0,0,0),(0,1,0,0,0,0,0,0,0),(0,0,0,0,0,1,0,0,0),(0,0,0,1,0,0,0,1,0),(0,0,0,0,2,-1,0,0,0),(0,0,0,0,0,1,0,0,0),(0,0,0,-2,0,0,2,0,2),(0,0,0,0,0,0,0,2,0),(0,0,0,0,0,0,0,-1,3))
J=((1,0,-1,0,1,0,0,0,0),(0,1,0,0,0,0,0,0,0),(1,0,-1,0,0,0,0,0,0),(0,0,0,-1,0,0,1,0,0),(1,0,0,0,1,0,0,0,0),(1,0,0,0,0,0,0,0,0),(0,0,0,-2,0,0,2,0,2),(0,0,0,0,0,0,1,0,0),(0,0,0,0,0,0,1,0,1))
matrices=[O,Z,J]
sparse=[[[ (i,x) for i,x in enumerate(row) if x] for row in M] for M in matrices]
def mv(M,x):return tuple(sum(c*x[i] for i,c in row) for row in M)

def physical(rows,root):
 n=len(rows);F=(1<<n)-1;h=[0]*(F+1);e=[0]*(F+1);s=[0]*(F+1);h[0]=1
 def dfs(first,last,mask):
  h[mask]+=1;e[mask]+=last==root;s[mask]+=first==root
  todo=rows[last]&~mask
  while todo:
   bit=todo&-todo;todo-=bit
   dfs(first,bit.bit_length()-1,mask|bit)
 for i in range(n):dfs(i,i,1<<i)
 c=[0]*(F+1)
 for mask in range(F+1):
  a=mask
  while True:
   c[mask]+=h[a]*h[mask^a]
   if not a:break
   a=(a-1)&mask
 return (h[F],s[F],e[F],sum(e[a]*h[F^a] for a in range(F+1)),sum(s[a]*h[F^a] for a in range(F+1)),h[F^(1<<root)],c[F],c[F^(1<<root)],sum(s[a]*c[F^a] for a in range(F+1)))

def add(rows,root,m):
 n=len(rows);out=[x|(1<<n) for x in rows]+[0]
 if m:
  out[root]^=1<<n;out[n]=1<<root
 return out,n if m==2 else root

seeds=[(677,2,(3,3),[2,12,17,101,75,23,167,63],[2,12,17,165,43,71,31,119]),(759,1,(4,3),[2,44,49,133,11,89,31,119],[2,52,25,35,105,197,15,95]),(1547,0,(4,3),[2,148,97,7,45,75,27,125],[2,148,41,67,45,75,23,125]),(821,2,(5,5),[4,9,50,85,163,75,151,47],[4,9,50,85,99,139,167,31])]
report={'status':'PASS','literal_count_method':'enumerate all directed simple words','word_length_cap':11,'pairs':[]}
for X,root,js,*rows in seeds:
 ps=[physical(a,root) for a in rows];delta=tuple(a-b for a,b in zip(*ps));offset=-delta[1]+2*(js[0]-js[1]);literal=0
 for side in range(2):
  q=[(rows[side],root,ps[side],0)]
  while q:
   a,v,p,depth=q.pop();need(physical(a,v)==p,('literal',X,side,depth));literal+=1
   if depth<2:
    for m,M in enumerate(sparse):
     b,w=add(a,v,m);q.append((b,w,mv(M,p),depth+1))
 count=sync=0;q=[(delta,False,0)]
 while q:
  p,seen,depth=q.pop();D=p[0]+p[3]+offset
  need((D==0)==(not seen),('language',X,p,seen,depth))
  if seen:need(D>0 if X==821 else D<0,('sign',X,D))
  else:sync+=1
  count+=1
  if depth<11:
   for m,M in enumerate(sparse):q.append((mv(M,p),seen or m==2,depth+1))
 report['pairs'].append({'base':X,'seed_profiles':ps,'difference':delta,'literal_profiles':literal,'all_words':count,'accepted_words':sync})
print(json.dumps(report,indent=2,sort_keys=True))
