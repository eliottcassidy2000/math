#!/usr/bin/env python3
"""Standalone exact controls for central Apéry saturation, no third-party packages."""
from itertools import product
from math import gcd,prod
import json

def check(p,s):
 if not p:raise RuntimeError(s)
def residues(h,z):
 D=prod(h);w=[D//a for a in h]
 return [z*pow(b,-1,a)%a for a,b in zip(h,w)]
def normal(h,z):
 D=prod(h);w=[D//a for a in h];r=residues(h,z);R=sum(a*b for a,b in zip(w,r))
 return r,(z-R)//D

def decode(h,m,z):
 r,c=normal(h,z);caps=[(b-a)//d for a,b,d in zip(r,m,h)]
 if c<0 or min(caps)<0 or c>sum(caps):return None
 ans=[]
 for a,d,b in zip(r,h,caps):
  k=min(c,b);c-=k;ans.append(a+d*k)
 return ans

def cyclicrun(R,n):
 if len(R)==n:return None
 return max((next((k for k in range(n) if (a+k)%n not in R),n) for a in range(n)),default=0)
def maxrun(S):
 best=cur=0;last=None
 for x in sorted(S):
  cur=cur+1 if last is not None and x==last+1 else 1;best=max(best,cur);last=x
 return best

def verify():
 cases=0;tests=0
 for h in [(3,5),(3,7),(5,7),(3,5,7),(5,7,9),(3,5,11)]:
  D=prod(h);w=[D//a for a in h];F=(len(h)-1)*D-sum(w)
  G=[z for z in range(F+1) if normal(h,z)[1]<0]
  check(len(G)==(F+1)//2 and max(G)==F,'symmetric semigroup gaps')
  for z in range(F+1):check((z in G)!=(F-z in G),'symmetry')
  for q in range(1,4):
   m=[q*a-1 for a in h];S={sum(a*b for a,b in zip(w,d)) for d in product(*(range(x+1) for x in m))};M=sum(a*b for a,b in zip(w,m))
   for z in range(M+1):
    got=decode(h,m,z);pred=normal(h,z)[1]>=0 and normal(h,M-z)[1]>=0
    check((z in S)==pred and (got is not None)==pred,'two sided exact membership')
    if got is not None:check(sum(a*b for a,b in zip(w,got))==z,'physical digit decoder')
    tests+=1
   if M>=2*F+1:check(set(range(M+1))-S==set(G)|{M-g for g in G},'exact two boundary gaps')
   cases+=1
  for E in product(*([[0],[0,1],list(range(a-1)),list(range(a))] for a in h)):
   S={sum(a*b for a,b in zip(w,d)) for d in product(*E)}
   for a,b,e in zip(h,w,E):
    R={(b*x)%a for x in e};bound=cyclicrun(R,a)
    if bound is not None:check(maxrun(S)<=bound,'projected run bound')
   cases+=1
 # Exact known physical fusion packets from independently counted predecessor sources.
 # Determine output support using the known committed literal baseline, not arbitrary matrix realizability.
 support={474053+2*(169*a+15*b) for a,b in product(range(2),repeat=2)}
 check(support=={474053,474083,474391,474421},'published literal crosslane support')
 check(maxrun({(x-474053)//2 for x in support})==1,'known physical pair cannot be adjacent')
 # Small illustrative arithmetic box, expressly not a physically supplied bank.
 S={5*a+3*b for a in range(9) for b in range(9)}
 check(set(range(65))-S=={1,2,4,7,57,60,62,63},'finite boundary illustration')
 return dict(status='PASS',arithmetic_cases=cases,target_checks=tests,physical_support=sorted(support),illustrative_holes=[1,2,4,7,57,60,62,63],scope='Exact arithmetic and inherited physical support; no new physical complete-residue banks or Chan-Pak proof')
if __name__=='__main__':print(json.dumps(verify(),sort_keys=True,indent=2))
