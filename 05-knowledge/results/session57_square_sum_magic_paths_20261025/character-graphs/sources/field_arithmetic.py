"""Independent tuple-polynomial field arithmetic; no producer imports."""
import itertools,json,math
from pathlib import Path
class K:
 def __init__(self,p,f):
  self.p=p; self.f=f; self.n=len(f)-1; self.E=list(itertools.product(range(p),repeat=self.n)); self.z=(0,)*self.n;self.one=(1,)+(0,)*(self.n-1)
  self.sq={self.mul(x,x) for x in self.E}; self.ch={x:0 if x==self.z else 1 if x in self.sq else -1 for x in self.E}
  assert len(self.sq)==(len(self.E)+1)//2
  assert all(self.pow(x,len(self.E)-1)==self.one for x in self.E if x!=self.z)
 def add(self,a,b):return tuple((x+y)%self.p for x,y in zip(a,b))
 def scale(self,k,a):return tuple(k*x%self.p for x in a)
 def mul(self,a,b):
  c=[0]*(2*self.n-1)
  for i,x in enumerate(a):
   for j,y in enumerate(b):c[i+j]+=x*y
  for i in range(len(c)-1,self.n-1,-1):
   for j in range(self.n):c[i-self.n+j]-=c[i]*self.f[j]
  return tuple(x%self.p for x in c[:self.n])
 def pow(self,a,n):
  b=self.one
  while n:
   if n%2:b=self.mul(b,a)
   a=self.mul(a,a);n//=2
  return b
 def factors(self,u,v):return [self.add(self.one,self.add(self.scale(s,u),self.scale(t,v))) for s,t in [(1,0),(-1,0),(1,1),(1,-1),(-1,1),(-1,-1)]]
