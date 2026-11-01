"""Exact priority-aware acceleration of specified FRACTRAN instruction words."""
from fractions import Fraction
from math import prod

def factor(n):
 out={};p=2
 while p*p<=n:
  while n%p==0:out[p]=out.get(p,0)+1;n//=p
  p=3 if p==2 else p+2
 if n>1:out[n]=out.get(n,0)+1
 return out

class Program:
 def __init__(self,rows):
  self.rows=[Fraction(a,b) for a,b in rows]
  if any(f<=0 for f in self.rows):raise ValueError('positive fractions required')
  pairs=[(factor(f.numerator),factor(f.denominator)) for f in self.rows]
  self.primes=sorted(set().union(*(set(a)|set(b) for a,b in pairs))) if pairs else []
  self.need=[[b.get(p,0) for p in self.primes] for a,b in pairs]
  self.delta=[[a.get(p,0)-b.get(p,0) for p in self.primes] for a,b in pairs]
 def encode(self,N):
  if type(N)!=int or N<=0:raise ValueError('positive integer state required')
  x=[]
  for p in self.primes:
   e=0
   while N%p==0:N//=p;e+=1
   x.append(e)
  return x,N
 def decode(self,x,inert=1):return inert*prod(p**e for p,e in zip(self.primes,x))
 def first(self,x):
  return next((i for i,b in enumerate(self.need) if all(v>=a for v,a in zip(x,b))),None)
 def step(self,x):
  i=self.first(x)
  return None if i is None else (i,[v+d for v,d in zip(x,self.delta[i])])

def eligible_interval(base,direction,need):
 """Integer t>=0 such that base+t*direction satisfies one denominator guard."""
 lo=0;hi=None
 for x,d,b in zip(base,direction,need):
  if d>0:lo=max(lo,-((x-b)//d))
  elif d<0:
   bound=(x-b)//(-d);hi=bound if hi is None else min(hi,bound)
  elif x<b:return None
 if hi is not None and lo>hi:return None
 return lo,hi

def certify_word(p,x,word):
 if not word:raise ValueError('nonempty word required')
 if len(x)!=len(p.primes) or any(type(v)!=int or v<0 for v in x):raise ValueError('invalid valuation vector')
 if any(type(i)!=int or not 0<=i<len(p.rows) for i in word):raise ValueError('instruction index')
 D=[sum(p.delta[i][j] for i in word) for j in range(len(x))]
 offset=[0]*len(x);limit=None;limiter=None;guards=[]
 def restrict(t,reason):
  nonlocal limit,limiter
  if limit is None or t<limit:limit=t;limiter=reason
 for s,i in enumerate(word):
  b=[v+z for v,z in zip(x,offset)]
  if p.first(b)!=i:raise ValueError('word is not first-eligible at its first occurrence')
  own=eligible_interval(b,D,p.need[i])
  if own is None or own[0]!=0:raise RuntimeError('inconsistent own guard')
  if own[1] is not None:restrict(own[1]+1,{'phase':s,'kind':'own_denominator','instruction':i})
  earlier=[]
  for j in range(i):
   hit=eligible_interval(b,D,p.need[j])
   if hit is not None:
    if hit[0]==0:raise RuntimeError('earlier guard unexpectedly active')
    restrict(hit[0],{'phase':s,'kind':'earlier_instruction','instruction':j})
    earlier.append({'instruction':j,'interval':hit})
  guards.append({'phase':s,'instruction':i,'prefix':offset[:],'own_interval':own,'earlier_intervals':earlier})
  offset=[v+d for v,d in zip(offset,p.delta[i])]
 if limit is not None and limit<1:raise RuntimeError('zero valid repetitions')
 return {'word':list(word),'direction':D,'repetitions':limit,'limiter':limiter,'guards':guards,
  'primitive_steps':None if limit is None else len(word)*limit,
  'after':None if limit is None else [v+limit*d for v,d in zip(x,D)]}

def run(p,x,lookahead=12,max_blocks=10000):
 if type(lookahead)!=int or lookahead<1:raise ValueError('positive lookahead')
 x=list(x);steps=0;records=[]
 for _ in range(max_blocks):
  if p.first(x) is None:return {'status':'halted','state':x,'primitive_steps':steps,'blocks':records}
  word=[];probe=x[:];best=None
  for _ in range(lookahead):
   step=p.step(probe)
   if step is None:break
   i,probe=step;word.append(i);c=certify_word(p,x,word)
   if c['repetitions'] is None:
    return {'status':'certified_infinite_word','state':x,'primitive_steps_before_word':steps,'infinite_word':c,'blocks':records}
   if best is None or c['primitive_steps']>best['primitive_steps']:best=c
  records.append({'start':x,'word':best['word'],'repetitions':best['repetitions'],'primitive_steps':best['primitive_steps']})
  x=best['after'];steps+=best['primitive_steps']
 return {'status':'budget_exhausted','state':x,'primitive_steps':steps,'blocks':records}

PRIMEGAME=[(17,91),(78,85),(19,51),(23,38),(29,33),(77,29),(95,23),(77,19),(1,17),(11,13),(13,11),(15,2),(1,7),(55,1)]
