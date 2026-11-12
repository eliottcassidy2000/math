#!/usr/bin/env python3
"""Standalone exact receiver-carry and local-action controls. No project imports."""
from functools import lru_cache
from itertools import permutations,product
from math import gcd
import json,random

def check(ok,msg):
 if not ok:raise RuntimeError(msg)
def reverse(T):return tuple(sum(1<<u for u in range(len(T)) if T[u]>>v&1) for v in range(len(T)))
def hp_subsets(T):
 n=len(T);incoming=reverse(T)
 @lru_cache(None)
 def end(S,v):
  R=S^(1<<v)
  if not R:return 1
  return sum(end(R,u) for u in range(n) if R>>u&1 and incoming[v]>>u&1)
 H=[1]+[sum(end(S,v) for v in range(n) if S>>v&1) for S in range(1,1<<n)]
 return H,end
def packet(T):
 H,end=hp_subsets(T);RH,rend=hp_subsets(reverse(T));F=(1<<len(T))-1;D=F^1;vals=[]
 for HH,ee in ((RH,rend),(H,end)):
  mass=signed=0
  for S in range(1,F+1,2):
   x=ee(S,0)*HH[F^S];mass+=x;signed+=x*(1 if S.bit_count()%2 else -1)
  vals.extend((mass,signed))
 U,d,V,e=vals;Q=sum(H[S]*H[D^S] for S in range(F+1) if not S&1)
 return dict(H=H[F],h=H[D],Q=Q,U=U,V=V,d=d,e=e)
def permutation_packet(T):
 n=len(T);F=(1<<n)-1;paths={0:[()]}
 for S in range(1,F+1):paths[S]=[p for p in permutations([i for i in range(n) if S>>i&1]) if all(T[a]>>b&1 for a,b in zip(p,p[1:]))]
 U=V=d=e=0
 for S in range(1,F+1,2):
  sg=1 if S.bit_count()%2 else -1;o=len(paths[F^S]);a=sum(p[0]==0 for p in paths[S])*o;b=sum(p[-1]==0 for p in paths[S])*o
  U+=a;V+=b;d+=sg*a;e+=sg*b
 return dict(H=len(paths[F]),h=len(paths[F^1]),Q=sum(len(paths[S])*len(paths[(F^1)^S]) for S in range(F+1) if not S&1),U=U,V=V,d=d,e=e)
def mat(P):
 h=P['h'];return [[h,P['U']-h,P['d']-h,P['H']-h],[0,P['Q'],0,P['V']-h],[0,0,0,P['e']-h],[0,0,0,h]]
def mul(A,B):return [[sum(A[i][r]*B[r][j] for r in range(4)) for j in range(4)] for i in range(4)]
def decode(M):
 h=M[0][0];return dict(H=h+M[0][3],h=h,Q=M[1][1],U=h+M[0][1],V=h+M[1][3],d=h+M[0][2],e=h+M[2][3])
def word(Ps):
 M=[[int(i==j) for j in range(4)] for i in range(4)]
 for P in Ps:M=mul(M,mat(P))
 return decode(M)
def attach(D,cut):return (cut<<1,)+tuple((r<<1)|(1-((cut>>i)&1)) for i,r in enumerate(D))
def fuse(A,B):
 n=len(A);m=len(B);T=[0]*(n+m-1)
 for C,ids in ((A,list(range(n))),(B,[0]+list(range(n,n+m-1)))):
  for i,u in enumerate(ids):
   for j,v in enumerate(ids):
    if C[i]>>j&1:T[u]|=1<<v
 for i in range(1,n):
  for j in range(n,n+m-1):T[i]|=1<<j
 return tuple(T)
