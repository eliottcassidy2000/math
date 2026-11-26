"""Cycle-type Burnside and a Hamilton-path marked census through order 7."""
from collections import Counter
from fractions import Fraction
from itertools import combinations
from math import factorial,prod


def partitions(n,start=1):
 if not n:
  yield ()
 for q in range(start,n+1):
  for rest in partitions(n-q,q):
   yield (q,)+rest


def cycle_weight(parts,reverse=False):
 n=sum(parts)
 p=[]
 for s in parts:
  i=len(p);p.extend(range(i+1,i+s));p.append(i)
 pairs=set(combinations(range(n),2))
 b=0
 while pairs:
  i,j=min(pairs);u,v=i,j;t=0
  while True:
   pairs.remove(tuple(sorted((u,v))))
   u,v=p[u],p[v];t+=1
   if {u,v}=={i,j}:break
  swapped=u==j
  if swapped!=bool(reverse and t%2):return 0
  b+=1
 return 1<<b


def burn(n,reverse=False):
 s=0
 for parts in partitions(n):
  z=prod(k**v*factorial(v) for k,v in Counter(parts).items())
  s+=Fraction(cycle_weight(parts,reverse),z)
 assert int(s)==s
 return int(s)


def markeds(n):
 if n==0:return Counter({1:1}),Counter(),Counter({1:1}),Counter(),Counter()
 tiles=[(i,j) for i,j in combinations(range(n),2) if i+1<j]
 unseen=set(range(1<<len(tiles)))
 allc,strong,conv,cs,fiber=[Counter() for _ in range(5)]
 def code(a,order):
  return sum(((a[order[i]]>>order[j])&1)<<k for k,(i,j) in enumerate(tiles))
 while unseen:
  w=unseen.pop();a=[0]*n
  for i in range(n-1):a[i]|=1<<(i+1)
  for k,(i,j) in enumerate(tiles):
   if w>>k&1:a[i]|=1<<j
   else:a[j]|=1<<i
  b=set()
  def path(order,left):
   if not left:
    b.add(code(a,order))
    return 1
   possible=left&a[order[-1]];ct=0
   while possible:
    t=possible&-possible;possible^=t
    ct+=path(order+(t.bit_length()-1,),left^t)
   return ct
  U=(1<<n)-1
  h=sum(path((i,),U^(1<<i)) for i in range(n))
  assert w in b and h%len(b)==0
  unseen.difference_update(b)
  ar=[sum((a[j]>>i&1)<<j for j in range(n)) for i in range(n)]
  sc=code(ar,tuple(reversed(range(n)))) in b
  st=all(any(a[j]&((1<<i)-1) for j in range(i,n)) for i in range(1,n))
  allc[h]+=1;fiber[h,len(b),h//len(b)]+=1
  if st:strong[h]+=1
  if sc:conv[h]+=1
  if sc and st:cs[h]+=1
 return allc,strong,conv,cs,fiber


def mul(x,y):
 r=Counter()
 for a,c in x.items():
  for b,d in y.items():r[a*b]+=c*d
 return r


def replay():
 prev=[]
 print("n A000568 self-converse converse-merged strong H-covered; largest H fibre collision sizes")
 for n in range(8):
  A,S,C,T,F=markeds(n)
  a,c=burn(n),burn(n,True)
  assert sum(A.values())==a and sum(C.values())==c
  if n:
   x=Counter()
   for i in range(1,n+1):x+=mul(prev[i][1] if i<n else S,prev[n-i][0])
   assert x==A
   y=T.copy()
   for i in range(1,n//2+1):
    y+=mul(Counter({h*h:q for h,q in prev[i][1].items()}),prev[n-2*i][2])
   assert y==C
  prev.append((A,S,C,T))
  largest=sorted(((q,h) for h,q in A.items()),reverse=True)[:3]
  print(n,a,c,(a+c)//2,sum(S.values()),len(A),largest)
 print("Burnside at 10",burn(10),burn(10,True))
 f=prev[7][0]; print("n=7 distinct counts", sorted(f))
if __name__=='__main__':replay()
