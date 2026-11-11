from collections import Counter
from functools import lru_cache
from itertools import combinations
from math import comb

@lru_cache(None)
def wrap_partition(part):
 out=Counter({tuple(sorted(part+(1,1))):1});n=len(part)
 # components equal sized compressed by selected multiplicities
 counts=Counter(part);items=list(counts.items())
 def rec(k,selected,size,weight,rest):
  if k==len(items):
   if selected%2:out[tuple(sorted(rest+[size+2]))]+=2*weight
   return
  a,m=items[k]
  for t in range(m+1):rec(k+1,selected+t,size+t*a,weight*comb(m,t),rest+[a]*(m-t))
 rec(0,0,0,1,[])
 return out

def iterate(parts,r):
 for t in range(r):
  new=Counter()
  for p,c in parts.items():
   for q,d in wrap_partition(p).items():new[q]+=c*d
  parts=new
 return parts

def two_colour(parts):
 N=sum(next(iter(parts)));out=[0]*(N+1)
 for part,c in parts.items():
  p=[c]
  for j in part:
   v=[0]*(len(p)+j)
   for i,a in enumerate(p):v[i]+=a;v[i+j]+=a
   p=v
  for i,a in enumerate(p):out[i]+=a
 return out

def oddparts(n,minp=1):
 if n==0:yield ();return
 for a in range(minp,n+1,2):
  for tail in oddparts(n-a,a):yield(a,)+tail

def rem(p,m=16):return [sum(v*(-1)**(i//(m//2)) for i,v in enumerate(p) if i%(m//2)==j) for j in range(m//2)]

if __name__=='__main__':
 for r in range(7):
  p=two_colour(iterate({(1,):1},r));print(r,p[0],sum(p))
