#!/usr/bin/env python3
"""Period-state canonical ranking of the shared-residue bank action."""
from functools import lru_cache
from itertools import product
from math import gcd, isqrt
import json, random
checks=0
def need(x,msg):
 global checks
 checks+=1
 if not x: raise RuntimeError(msg)

def mmul(A,B):
 n=len(A); C=[[0]*n for _ in range(n)]
 for i in range(n):
  for k,a in enumerate(A[i]):
   if a:
    for j,b in enumerate(B[k]):
     if b:C[i][j]+=a*b
 return C

def trace_power(A,e):
 n=len(A); R=[[int(i==j) for j in range(n)] for i in range(n)]
 while e:
  if e&1:R=mmul(R,A)
  e//=2
  if e:A=mmul(A,A)
 return sum(R[i][i] for i in range(n))

def digits(x,m,A):
 out=[0]*m
 for i in range(m-1,-1,-1):out[i]=x%A;x//=A
 return out

def number(w,A):
 out=0
 for x in w:out=out*A+x
 return out

def rotations(w):return [w[i:]+w[:i] for i in range(len(w))]
def least(w):return min(rotations(w))
def period(w):
 return next(r for r in range(1,len(w)+1) if len(w)%r==0 and all(w[i]==w[i%r] for i in range(len(w))))

def threshold_matrix(t,A):
 """Avoid every proper first-difference prefix smaller than t."""
 m=len(t); special=set(t); M=[[0]*m for _ in range(m)]
 for s in range(m):
  active=[j for j in range(s+1) if j==0 or t[s-j:s]==t[:j]]
  req=max(t[j] for j in active)
  # Allowed symbols not occurring in t always transition to empty prefix.
  M[s][0]=A-req-sum(a>=req for a in special)
  for a in special:
   if a<req:continue
   word=t[:s]+[a]
   dest=max(j for j in range(m) if j==0 or (j<=len(word) and word[-j:]==t[:j]))
   M[s][dest]+=1
 return M

@lru_cache(None)
def primitive(A,ell,p):
 return A if ell==1 else (A**ell-A**(ell//p))//ell

@lru_cache(None)
def necklace_ge_counts(A,m,p,t):
 """Return exact-period necklace counts with least length-m word >= t."""
 ells=[]; z=1
 while z<=m:ells.append(z);z*=p
 if t<=0:return tuple(primitive(A,e,p) for e in ells)
 if t>=A**m:return (0,)*len(ells)
 M=threshold_matrix(digits(t,m,A),A)
 traces=[trace_power(M,e) for e in ells]
 out=[]
 for i,e in enumerate(ells):
  v=traces[i]-(traces[i-1] if i else 0)
  need(v>=0 and v%e==0,'primitive trace nonnegative integral')
  out.append(v//e)
 return tuple(out)

class Bank:
 def __init__(self,p,k,Q=2):
  if not (isinstance(p,int) and p>=3 and p%2 and all(p%d for d in range(2,isqrt(p)+1))):
   raise ValueError('p must be an odd prime')
  if not (isinstance(k,int) and k>=2 and isinstance(Q,int) and Q>=2):
   raise ValueError('k>=2 and Q>=2 are required')
  self.p=p;self.k=k;self.Q=Q;self.n=p**(k-1);self.h=(p-1)//2
  # These named blocks are ordered. Equal-sized copies are not permuted.
  self.blocks=[('A',k-1)]*self.h+[('B',k-1)]*self.h
  self.blocks += [('A',j) for j in range(k-1) for _ in range(p-1)]
  self.states=[('C',r) for r in range(k-1)]+[('A',0),('B',0),('I',0)]
  # Bottom-up evaluation avoids recursion depth depending on block count.
  self._totals={(len(self.blocks),H):1 for H in self.states}
  for i in range(len(self.blocks)-1,-1,-1):
   for H in self.states:
    kind,j,q,m,A=self.options(i,H);out=0;e=1;r=0
    while e<=m:
     out+=primitive(A,e,p)*self._totals[i+1,self.restrict(H,kind,q+r)]
     e*=p;r+=1
    self._totals[i,H]=out
 def image_q(self,H,kind,j):
  typ,r=H
  if typ=='C':return min(r,j)
  if typ=='I' or typ==kind:return j
  return min(self.k-2,j)
 def restrict(self,H,kind,v):
  typ,r=H
  if typ=='C':return (kind,0) if v==self.k-1 else ('C',max(r,v))
  if typ=='I' or typ==kind or v<self.k-1:return H
  return ('I',0)
 def elements(self,H):
  typ,r=H
  return [(a,b) for a in range(self.n) for b in range(self.p) if
          (typ=='C' and a%self.p**r==0) or
          (typ=='A' and a==0) or
          (typ=='B' and (a+self.n//self.p*b)%self.n==0) or
          (typ=='I' and a==0 and b==0)]
 def options(self,i,H):
  kind,j=self.blocks[i];q=self.image_q(H,kind,j);m=self.p**(j-q);A=self.Q**(self.p**q)
  return kind,j,q,m,A
 def total(self,i,H):return self._totals[i,H]
 def below(self,i,H,t):
  kind,j,q,m,A=self.options(i,H);ge=necklace_ge_counts(A,m,self.p,t);out=0;e=1
  for r,g in enumerate(ge):
   out+=(primitive(A,e,self.p)-g)*self.total(i+1,self.restrict(H,kind,q+r));e*=self.p
  return out
 def word_state(self,i,H,w):
  kind,j,q,m,A=self.options(i,H);packet=self.p**q
  ws=[number(w[z:z+packet],self.Q) for z in range(0,len(w),packet)]
  need(ws==least(ws),'block not canonical in stabilizer')
  e=period(ws);r=0
  while self.p**r<e:r+=1
  return self.restrict(H,kind,q+r)
 def unrank(self,r):
  if not isinstance(r,int):raise ValueError('rank must be an integer')
  H=('C',0);out=[]
  need(0<=r<self.total(0,H),'rank range')
  for i,(_,j) in enumerate(self.blocks):
   lo=0;hi=self.Q**(self.p**j)-1
   while lo<hi:
    mid=(lo+hi)//2
    if self.below(i,H,mid+1)>r:hi=mid
    else:lo=mid+1
   r-=self.below(i,H,lo)
   w=digits(lo,self.p**j,self.Q);out.append(w);H=self.word_state(i,H,w)
  need(r==0,'unrank final remainder')
  return out
 def rank(self,words):
  if len(words)!=len(self.blocks) or any(len(w)!=self.p**j or any(not isinstance(a,int) or a<0 or a>=self.Q for a in w) for w,(_,j) in zip(words,self.blocks)):
   raise ValueError('wrong block shape or color digit')
  H=('C',0);out=0
  for i,w in enumerate(words):
   out+=self.below(i,H,number(w,self.Q));H=self.word_state(i,H,w)
  return out
 def act(self,words,a,b):
  out=[]
  for (kind,j),w in zip(self.blocks,words):
   shift=(a+(self.n//self.p*b if kind=='B' else 0))%len(w)
   out.append(w[shift:]+w[:shift])
  return out
 def canonical(self,words):
  return min(self.act(words,a,b) for a in range(self.n) for b in range(self.p))
 def flat(self,words):return tuple(x for w in words for x in w)
 def split(self,x):
  ws=[];idx=0
  for _,j in self.blocks:
   ell=self.p**j;ws.append(list(x[idx:idx+ell]));idx+=ell
  return ws

def test_automata():
 cases=0
 for p,m,A in [(2,1,3),(2,2,3),(2,4,2),(3,3,2),(5,5,2)]:
  reps={tuple(least(list(w))) for w in product(range(A),repeat=m)}
  periods=[p**r for r in range(m.bit_length()+1) if p**r<=m]
  for t in range(A**m+1):
   got=necklace_ge_counts(A,m,p,t)
   want=tuple(sum(period(list(w))==e and number(w,A)>=t for w in reps) for e in periods)
   need(got==want,'automaton all-threshold necklace enumeration');cases+=1
 return cases

def verify_bank(p,k):
 B=Bank(p,k);D=sum(p**j for _,j in B.blocks);rng=random.Random(p*100+k)
 for H in B.states:
  els=B.elements(H);need(bool(els),'nonempty subgroup')
  for kind,j in B.blocks:
   L=p**j;q=B.image_q(H,kind,j)
   shifts={(a+(B.n//p*b if kind=='B' else 0))%L for a,b in els}
   need(shifts==set(range(0,L,p**q)),'symbolic image subgroup')
   for v in range(q,j+1):
    new=B.restrict(H,kind,v)
    want={(a,b) for a,b in els if (a+(B.n//p*b if kind=='B' else 0))%p**v==0}
    need(set(B.elements(new))==want,'symbolic stabilizer transition')
 total=B.total(0,('C',0));num=0
 for a in range(B.n):
  for b in range(p):
   cycles=sum(gcd(a+(B.n//p*b if kind=='B' else 0),p**j) for kind,j in B.blocks)
   num+=2**cycles
 need(num==p*B.n*total,'period transfer vs literal Burnside')
 if D<=8:
  reps=sorted({B.flat(B.canonical(B.split(x))) for x in product(range(2),repeat=D)})
  need(len(reps)==total,'complete orbit cardinality')
  for r,x in enumerate(reps):
   need(B.flat(B.unrank(r))==x,'exhaustive unrank vs canonical enumeration')
   need(B.rank(B.split(x))==r,'exhaustive rank vs canonical enumeration')
  tested=len(reps)
 else:
  trials=8 if k==3 else 16
  ranks=[0,total-1]+[rng.randrange(total) for _ in range(trials)]
  for r in ranks:
   w=B.unrank(r);need(B.canonical(w)==w,'rank-selected word literal canonical')
   need(B.rank(w)==r,'rank-unrank inverse')
  for _ in range(trials):
   raw=B.split([rng.randrange(2) for _ in range(D)]);w=B.canonical(raw)
   need(B.unrank(B.rank(w))==w,'literal canonical random word selected back')
  tested=len(ranks)+trials
 return {'p':p,'k':k,'coordinates':D,'subgroup_states':len(B.states),'orbits':total,'rank_controls':tested}

if __name__=='__main__':
 auto=test_automata()
 data=[verify_bank(3,2),verify_bank(5,2),verify_bank(5,3)]
 print(json.dumps({'status':'PASS','checks':checks,'automaton_thresholds':auto,'banks':data,'note':'p=3 is a finite permutation-model control only, not a ramified-ring claim'},indent=2))
