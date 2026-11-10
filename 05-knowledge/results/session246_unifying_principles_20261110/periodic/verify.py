#!/usr/bin/env python3
"""Exact independent fixed-apex tape, flag and periodic-Jordan checks."""
from itertools import permutations,product
from fractions import Fraction as Q
import json,sys

def ck(x,m):
 if not x: raise RuntimeError(m)
def mul(A,B):
 n=len(A);C=[[0]*n for _ in range(n)]
 for i in range(n):
  for z in range(n):
   if A[i][z]:
    for j in range(n): C[i][j]+=A[i][z]*B[z][j]
 return C
def eye(n):return [[int(i==j) for j in range(n)] for i in range(n)]
def shifted(A,c):return [[v-c*(i==j) for j,v in enumerate(row)] for i,row in enumerate(A)]
def rank(A):
 A=[[Q(x) for x in row] for row in A];r=0
 for j in range(len(A[0])):
  z=next((z for z in range(r,len(A)) if A[z][j]),None)
  if z is None:continue
  A[r],A[z]=A[z],A[r];v=A[r][j];A[r]=[x/v for x in A[r]]
  for z in range(r+1,len(A)):
   v=A[z][j]
   if v:A[z]=[x-v*y for x,y in zip(A[z],A[r])]
  r+=1
 return r

def letter(k,p,a):
 states=list(product(range(3),repeat=k));states=[s+(f,) for s in states for f in range(2)]
 ix={s:i for i,s in enumerate(states)};A=[[0]*len(states) for _ in states]
 for s in states:
  for g in range(k+1):
   b=(g==0 or bool(a>>p[g-1]&1))
   if g==k:
    terms=[(0,1)] if s[g]==0 or b else []
   elif s[g]==0:terms=[]
   else:
    good=not(a>>p[g]&1);factor=1 if s[g]==1 else int(b)
    terms=[(1,factor)]+([] if good else [(0,factor)])
   for out,v in terms:
    ss=list(s);ss[g]=out;A[ix[tuple(ss)]][ix[s]]+=v
 return states,A

def block(k,p,w):
 states,A=letter(k,p,w[0]);M=eye(len(A))
 for a in w:M=mul(letter(k,p,a)[1],M)
 return states,M

def H(k,w):
 total=0
 for p in permutations(range(k)):
  dp={(0,)*(k+1):1}
  for a in w:
   nd={}
   for s,v in dp.items():
    for g in range(k+1):
     if s[g]==0 and g>0 and not(a>>p[g-1]&1):continue
     ns=list(s);ns[g]=2 if g==k or not(a>>p[g]&1) else 1;ns=tuple(ns)
     nd[ns]=nd.get(ns,0)+v
   dp=nd
  for s,v in dp.items():
   if all(s[g]==2 or (s[g]==0 and (g==0 or p[g-1]<p[g])) for g in range(k)):total+=v
 return total

def graph(k,w):
 n=k+len(w);a=[0]*n
 for i in range(k):
  for j in range(i+1,k):a[i]|=1<<j
 for z,b in enumerate(w):
  v=k+z
  for i in range(k):
   if b>>i&1:a[i]|=1<<v
   else:a[v]|=1<<i
  for u in range(k,v):a[u]|=1<<v
 return a

def dpcount(a):
 n=len(a)
 if n==0:return 1
 dp=[[0]*n for _ in range(1<<n)]
 for i in range(n):dp[1<<i][i]=1
 for mask in range(1<<n):
  for i in range(n):
   if dp[mask][i]:
    for j in range(n):
     if not(mask>>j&1) and a[i]>>j&1:dp[mask|1<<j][j]+=dp[mask][i]
 return sum(dp[-1])

def solve(A,b):
 A=[[Q(v) for v in row]+[Q(v)] for row,v in zip(A,b)];n=len(b)
 for j in range(n):
  z=next(z for z in range(j,n) if A[z][j]);A[j],A[z]=A[z],A[j]
  t=A[j][j];A[j]=[v/t for v in A[j]]
  for z in range(n):
   if z!=j:
    t=A[z][j];A[z]=[v-t*u for v,u in zip(A[z],A[j])]
 return [row[-1] for row in A]

def main():
 if '--fail' in sys.argv:ck(False,'intentional hostile sentinel')
 nc=0;letters=0
 for k in range(3):
  for p in permutations(range(k)):
   for a in range(2**k):
    states,A=letter(k,p,a);letters+=1
    for i,s in enumerate(states):
     lab=sum(x==1 for x in s[:-1])+(s[-1]==0)
     ck(A[i][i]==lab,'diagonal label')
     for z,t in enumerate(states):
      if i!=z and A[i][z]:ck(sum(s)<sum(t),'strict rank decrease')
   for l in range(1,4):
    for w in product(range(2**k),repeat=l):
     states,M=block(k,p,w);R=eye(len(M))
     for _ in range(k+1):R=mul(M,R)
     for j in range(1,k+2):
      B=shifted(M,j**l)
      for _ in range(k-j+2):R=mul(B,R)
     ck(not any(any(row) for row in R),'universal annihilator');nc+=1
 ng=0
 for k in range(3):
  for l in range(5):
   for w in product(range(2**k),repeat=l):
    ck(H(k,w)==dpcount(graph(k,w)),'independent graph count');ng+=1
 examples=[]
 for k,w in [(1,(1,0)),(2,(3,2,1))]:
  cols=[(j,d) for j in range(1,k+2) for d in range(k-j+2)]
  # Account explicitly for the nilpotent zero transient by starting beyond k.
  start=k+1;n=len(cols);ms=range(start,start+n)
  cs=solve([[m**d*j**(len(w)*m) for j,d in cols] for m in ms],[H(k,w*m) for m in ms])
  for m in range(0,start+n+5):
   value=sum(c*m**d*j**(len(w)*m) for c,(j,d) in zip(cs,cols))
   ck(value==H(k,w*m),'exact scalar polynomial witness including m=0')
  ck(cs[cols.index((1,k))]!=0,'top degree survives scalar readout')
  if k==2:ck(cs[cols.index((2,1))]!=0,'middle degree survives')
  _,M=block(k,tuple(range(k)),w);B=shifted(M,1);P=eye(len(M));null=[]
  for r in range(1,k+3):P=mul(B,P);null.append(len(M)-rank(P))
  examples.append({'k':k,'mask_block':w,'coefficients':[[j,d,str(c)] for (j,d),c in zip(cols,cs)],'nullities_at_1_one_permutation':null,'counts':[H(k,w*m) for m in range(7)]})
 print(json.dumps({'status':'PASS','letter_blocks':letters,'annihilator_blocks':nc,'graph_controls':ng,'examples':examples},sort_keys=True,indent=2))
if __name__=='__main__':main()
