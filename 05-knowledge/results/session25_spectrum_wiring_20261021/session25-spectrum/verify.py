#!/usr/bin/env python3
"""Independent exact controls. No repository writes or external packages."""
from itertools import product, permutations
from math import gcd
import json
from pathlib import Path

def require(b,msg):
 if not b: raise RuntimeError(msg)
def state(w):
 h,a,p=1,0,0
 for b in w:
  if b==0:h+=a+p
  a,p=2*a+b,b
 return h,a,p

def slots(w):
 if not w:return 1
 return w[0]+1-w[-1]+sum(a==0 and b==1 for a,b in zip(w,w[1:]))+sum((1<<(i-j-1)) for j in range(len(w)) for i in range(j+1,len(w)) if w[j]==1 and w[i]==0)
def adjacency(w):
 n=len(w); a=[[False]*(n+1) for _ in range(n+1)]
 for i in range(n):
  for j in range(i+1,n):a[i][j]=True
  a[n][i]=bool(w[i]); a[i][n]=not w[i]
 return a

def dp(w):
 a=adjacency(w); n=len(a); d=[[0]*n for _ in range(1<<n)]
 for i in range(n):d[1<<i][i]=1
 for mask in range(1,1<<n):
  for i in range(n):
   if d[mask][i]:
    for j in range(n):
     if not mask>>j&1 and a[i][j]:d[mask|1<<j][j]+=d[mask][i]
 return sum(d[-1])
def brute(w):
 a=adjacency(w)
 return sum(all(a[i][j] for i,j in zip(p,p[1:])) for p in permutations(range(len(a))))
def strong(w):
 a=adjacency(w);n=len(a)
 for k in range(n):
  for i in range(n):
   for j in range(n):a[i][j]|=a[i][k] and a[k][j]
 return all(a[i][j] for i in range(n) for j in range(n) if i!=j)
def closed(L,t):
 q=1<<L; num=q*q*((q**(2*t))-1)+2*t*(q+1)
 require(num%((q+1)**2)==0,'nonintegral closed expression')
 return 1+num//((q+1)**2)
def choose_L(M):
 m=M;k=0
 while m%2==0:k+=1;m//=2
 L=max(1,k)
 while pow(2,L,m)!=1%m:L+=1
 return L

def run():
 counts={'word_dp':0,'word_permutations':0,'slot_formula':0,'connectivity':0,'closed_form':0,'moduli':0,'residue_instances':0}
 for n in range(10):
  for w in product((0,1),repeat=n):
   h=state(w)[0]
   require(h==dp(w),f'dp {w}');counts['word_dp']+=1
   require(h==slots(w),f'slots {w}');counts['slot_formula']+=1
   if n<=6:
    require(h==brute(w),f'brute {w}');counts['word_permutations']+=1
   if n>=2 and w[0]==1 and w[-1]==0:
    require(strong(w),f'connectivity {w}');counts['connectivity']+=1
 for L in range(1,11):
  for t in range(1,21):
   require(closed(L,t)==state((1,)*L+(0,)*L if t==1 else ((1,)*L+(0,)*L)*t)[0],'closed');counts['closed_form']+=1
 for M in range(1,513):
  L=choose_L(M);q=pow(2,L,M);D=((q-1)**2+1)%M;g=gcd(2,M)
  require(gcd(D,M)==g,f'gcd {M}')
  reached={(1+t*D)%M for t in range(1,M//g+1)}
  expected={r for r in range(M) if g==1 or r%2==1}
  require(reached==expected,f'coverage {M}')
  a=h0=0;h=1
  for t in range(1,M//g+1):
   h=(h+q*(q-1)*a+(q-1)**2+1)%M
   a=(q*q*a+q*(q-1))%M
   require(a==0 and h==(1+t*D)%M,f'block {M} {t}')
   if M<=32:require(closed(L,t)%M==h,f'exact residue {M} {t}')
   counts['residue_instances']+=1
  counts['moduli']+=1
 # Hostiles: counts cannot predict continuation; arbitrary L fails.
 require(state((0,0))[0]==state((0,1))[0]==1,'scalar equal')
 require(state((0,0,0))[0]==1 and state((0,1,0))[0]==3,'scalar failure')
 require({closed(1,t)%3 for t in range(1,25)}=={0,1},'untuned L hostile')
 require(state((1,0))[0]==3 and slots((1,0))==3,'orphan triangle')
 return {'status':'PASS','counts':counts,'hostiles':{'equal_H_00_vs_01_append0':[1,3],'L1_mod3_residues':[0,1],'C3_H':3},'family_L1_first8':[closed(1,t) for t in range(1,9)]}
if __name__=='__main__':
 result=run(); print(json.dumps(result,indent=2,sort_keys=True))
