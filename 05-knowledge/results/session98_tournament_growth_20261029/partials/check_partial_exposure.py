#!/usr/bin/env python3
"""Stdlib exhaustive partial-star orbit audit, bases of orders 0..6."""
from itertools import permutations, product
from collections import Counter
from fractions import Fraction
from math import factorial, gcd
import json
from pathlib import Path
CHECKS=0
def need(b,m):
 global CHECKS
 CHECKS+=1
 if not b: raise RuntimeError(m)
def transform(a,p):
 n=len(a); b=[0]*n
 for i in range(n):
  for j in range(n):
   if a[i]>>j&1:b[p[i]]|=1<<p[j]
 return tuple(b)
def canon(a):
 n=len(a)
 # Only permutations taking degree-sorted vertices to the corresponding labels.
 blocks={}
 for i,row in enumerate(a):blocks.setdefault(row.bit_count(),[]).append(i)
 orders=product(*(permutations(blocks[d]) for d in sorted(blocks)))
 best=None
 for parts in orders:
  order=sum(parts,()); p=[0]*n
  for j,i in enumerate(order):p[i]=j
  b=transform(a,p)
  if best is None or b<best:best=b
 return best if best is not None else ()
def autos(a):
 return [p for p in permutations(range(len(a))) if all(a[i].bit_count()==a[p[i]].bit_count() for i in range(len(a))) and transform(a,p)==a]
def move(x,p):
 y=[0]*len(x)
 for i,v in enumerate(x):y[p[i]]=v
 return tuple(y)
def key(x,g):return min(move(x,p) for p in g)
def cycles(p):
 seen=set(); ans=[]
 for i in range(len(p)):
  if i in seen:continue
  c=[];j=i
  while j not in seen:seen.add(j);c.append(j);j=p[j]
  ans.append(c)
 return ans
def extend(a,s):
 n=len(a)
 return tuple(row|((1-(s>>i&1))<<n) for i,row in enumerate(a))+(s,)
def verify(a):
 n=len(a);g=autos(a)
 states={key(x,g) for x in product(range(3),repeat=n)} # 0 unknown, 1 -, 2 +
 sizes={x:len({move(x,p) for p in g}) for x in states}
 layers=Counter(sum(v!=0 for v in x) for x in states)
 burn=Counter()
 for p in g:
  poly=Counter({0:1})
  for c in cycles(p):
   new=poly.copy()
   for d,v in poly.items():new[d+len(c)]+=2*v
   poly=new
  burn.update(poly)
 need(all(Fraction(burn[r],len(g))==layers[r] for r in range(n+1)),'layer Burnside')
 up={};down={}
 for x in states:
  u=Counter();d=Counter()
  for i,v in enumerate(x):
   for w in ((1,2) if v==0 else (0,)):
    y=list(x);y[i]=w;y=key(y,g)
    (u if v==0 else d)[y]+=1
  up[x]=u;down[x]=d
  r=sum(v!=0 for v in x)
  need(sum(u.values())==2*(n-r) and sum(d.values())==r,'degree')
 for x in states:
  for y,v in up[x].items():need(sizes[x]*v==sizes[y]*down[y][x],'balance')
 paths={(0,)*n:1}
 for r in range(n):
  nxt=Counter()
  for x,m in paths.items():
   for y,v in up[x].items():nxt[y]+=m*v
  paths=nxt
 for x,m in paths.items():need(m==factorial(n)*sizes[x],'path mass')
 flagged_loss=0
 for x in states:
  h=[p for p in g if move(x,p)==x];unknown=[i for i,v in enumerate(x) if v==0]
  completions=[]
  for word in product((1,2),repeat=len(unknown)):
   y=list(x)
   for i,v in zip(unknown,word):y[i]=v
   completions.append(tuple(y))
  flagged={key(y,h) for y in completions};unflagged={key(y,g) for y in completions}
  b=sum(2**sum(all(i in unknown for i in c) for c in cycles(p)) for p in h)
  need(Fraction(b,len(h))==len(flagged),'completion Burnside')
  if n<=4:
   partial=tuple(row|((1 if x[i]==1 else 0)<<n) for i,row in enumerate(a))+(sum((1<<i) for i,v in enumerate(x) if v==2),)
   aut=autos(partial);d=len(unknown)
   if d>=2:need(len(aut)==len(h) and len(aut)%2==1,'star center is intrinsic')
   if d==1:
    j=unknown[0]
    twins=all(((partial[j]>>i)&1)==((partial[n]>>i)&1) for i in range(n) if i!=j)
    need((len(aut)%2==0)==twins,'one-defect twin iff even symmetry')
    need(len(aut)==len(h)*(2 if twins else 1),'one-defect exact group index')
  need(len(flagged)>=len(unflagged),'forgetting')
  flagged_loss+=len(flagged)>len(unflagged)
 return layers,flagged_loss

def oddparts(n,lo=1):
 if n==0:yield ()
 for j in range(lo,n+1,2):
  for tail in oddparts(n-j,j):yield (j,)+tail

def partition_layers(n):
 ans=Counter()
 for parts in oddparts(n):
  mult=Counter(parts);den=1
  for j,m in mult.items():den*=j**m*factorial(m)
  exponent=sum((j-1)//2 for j in parts)+sum(gcd(a,b) for i,a in enumerate(parts) for b in parts[i+1:])
  poly=Counter({0:1})
  for j in parts:
   nxt=poly.copy()
   for r,v in poly.items():nxt[r+j]+=2*v
   poly=nxt
  for r,v in poly.items():ans[r]+=Fraction((1<<exponent)*v,den)
 return ans

def main():
 banks=[{()}];results=[]
 for n in range(7):
  if n:banks.append({canon(extend(a,s)) for a in banks[-1] for s in range(1<<(n-1))})
  aggregate=Counter();loss=0
  for a in sorted(banks[-1]):
   l,k=verify(a);aggregate.update(l);loss+=k
  need(aggregate==partition_layers(n),'independent odd-partition partial-state GF')
  results.append({'base_order':n,'tournament_classes':len(banks[-1]),'partial_layers':[aggregate[r] for r in range(n+1)],'partial_states_with_forgetting_collision':loss})
 need([len(b) for b in banks]==[1,1,1,2,4,12,56],'independent small count control')
 # Triangle: H is trivial; +?? gives four flagged completions, three unflagged.
 a=(2,4,1);g=autos(a);x=(2,0,0);h=[p for p in g if move(x,p)==x]
 comp=[(2,u,v) for u,v in product((1,2),repeat=2)]
 hostile={'parent_rows':a,'partial':x,'group_order':len(g),'stabilizer_order':len(h),'flagged_completions':len({key(y,h) for y in comp}),'unflagged_completions':len({key(y,g) for y in comp})}
 need((hostile['flagged_completions'],hostile['unflagged_completions'])==(4,3),'C3 hostile')
 out={'status':'FINITE-EXACT','base_orders':'0..6, all isomorphism classes; all ternary star states; every completion','results':results,'hostile':hostile,'checks':CHECKS}
 print(json.dumps(out,indent=2));Path(__file__).with_name('partial_exposure.json').write_text(json.dumps(out,indent=2)+'\n')
if __name__=='__main__':main()
