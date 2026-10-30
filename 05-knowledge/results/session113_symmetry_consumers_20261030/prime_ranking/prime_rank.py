#!/usr/bin/env python3
"""Binary prime-bank ranking under the named O_k action, not graph isomorphism."""
from bank_rank import Bank, digits, number, period, least, necklace_ge_counts, need
from functools import lru_cache

class PrimeBank(Bank):
 def __init__(self,p,k):
  if p<5: raise ValueError('ramified prime-bank theorem requires p>=5')
  super().__init__(p,k,2)
  self.depths=[0]*self.h+[1]*self.h
  self.depths += [d for j in range(k-1) for d in (2*k-2*j-2,2*k-2*j-1) for _ in range(self.h)]
  self.full=(1<<(2*k-1))-1
  self._prime={}
  # Iterative table: no recursion depth depending on the input dimension.
  for H in self.states:
   for S in range(self.full+1):self._prime[len(self.blocks),H,S]=int(S==0)
  for i in range(len(self.blocks)-1,-1,-1):
   for H in self.states:
    kind,j,q,_,_=self.options(i,H)
    for S in range(self.full+1):
     self._prime[i,H,S]=sum(self.count_period(j,q,v)*self._prime[i+1,self.next_H(H,kind,q,v),self.next_S(i,S,v)] for v in range(j+1))
 def next_S(self,i,S,v):
  d=self.depths[i]
  return sum(1<<(r-1) for r in range(1,2*self.k) if S>>(r-1)&1 and (d>=r or v<=(r-d-1)//2))
 def next_H(self,H,kind,q,v):return H if v<=q else self.restrict(H,kind,v)
 def count_period(self,j,q,v):
  W=2**(self.p**v)-(2**(self.p**(v-1)) if v else 0)
  return W//self.p**max(0,v-q)
 def prime_total(self,i=0,H=('C',0),S=None):
  return self._prime[i,H,self.full if S is None else S]
 def ge_periods(self,j,q,t):
  L=self.p**j
  def dividing(v):
   base=2**(self.p**v); rep=(2**L-1)//(base-1)
   return base-max(0,min(base,(t+rep-1)//rep))
  result=[dividing(v)-(dividing(v-1) if v else 0) for v in range(q+1)]
  if q<j:
   ge=necklace_ge_counts(2**(self.p**q),self.p**(j-q),self.p,t)
   result.extend(ge[1:])
  return result
 def below_prime(self,i,H,S,t):
  kind,j,q,_,_=self.options(i,H)
  return sum((self.count_period(j,q,v)-g)*self._prime[i+1,self.next_H(H,kind,q,v),self.next_S(i,S,v)] for v,g in enumerate(self.ge_periods(j,q,t)))
 def step(self,i,H,S,w):
  kind,j,q,_,_=self.options(i,H);packet=self.p**q
  ws=[number(w[z:z+packet],2) for z in range(0,len(w),packet)]
  need(ws==least(ws),'block is not stabilizer canonical')
  e=period(w);v=0
  while self.p**v<e:v+=1
  return self.next_H(H,kind,q,v),self.next_S(i,S,v)
 def unrank(self,r):
  if not isinstance(r,int) or not 0<=r<self.prime_total():raise ValueError('rank outside prime-bank range')
  H=('C',0);S=self.full;out=[]
  for i,(_,j) in enumerate(self.blocks):
   lo=0;hi=2**(self.p**j)-1
   while lo<hi:
    mid=(lo+hi)//2
    if self.below_prime(i,H,S,mid+1)>r:hi=mid
    else:lo=mid+1
   r-=self.below_prime(i,H,S,lo)
   w=digits(lo,self.p**j,2);out.append(w);H,S=self.step(i,H,S,w)
  need(r==0 and S==0,'prime unrank terminal');return out
 def rank(self,words):
  if len(words)!=len(self.blocks) or any(len(w)!=self.p**j or any(a not in (0,1) for a in w) for w,(_,j) in zip(words,self.blocks)):raise ValueError('binary block shape')
  H=('C',0);S=self.full;out=0
  for i,w in enumerate(words):
   out+=self.below_prime(i,H,S,number(w,2));H,S=self.step(i,H,S,w)
  if S:raise ValueError('bank is composite')
  return out
