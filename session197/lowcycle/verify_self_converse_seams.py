#!/usr/bin/env python3
"""Literal seam graphs, independent subset counters, and family minima."""
from pathlib import Path
from functools import lru_cache
import json,math
ROOT=Path(__file__).resolve().parent

def require(condition,message):
 if not condition:raise RuntimeError(message)

def factor(d):return 1+(1<<(d-1))
def predicted(gaps):return math.prod(factor(d) for d in gaps)
def graph(gaps):
 n=1+sum(gaps);rows=[((1<<n)-1)^((1<<(u+1))-1) for u in range(n)]
 a=0
 for d in gaps:
  b=a+d;rows[a]^=1<<b;rows[b]|=1<<a;a=b
 return rows

def validate(rows):
 n=len(rows)
 for u in range(n):
  require(not rows[u]>>u&1,'loop')
  for v in range(u+1,n):require(((rows[u]>>v)&1)+((rows[v]>>u)&1)==1,'orientation')

def reflect(rows):
 n=len(rows)
 return all(bool(rows[u]>>v&1)==bool(rows[n-1-v]>>(n-1-u)&1) for u in range(n) for v in range(n) if u!=v)

def physical(rows):
 n=len(rows);N=1<<n;root=n-1
 dp=[[0]*n for _ in range(N)];start=[[0]*n for _ in range(N)]
 for u in range(n):dp[1<<u][u]=1
 start[1<<root][root]=1
 for S in range(1,N):
  for v in range(n):
   if not S>>v&1:continue
   prev=S^(1<<v)
   if not prev:continue
   # Independent predecessor recurrence rather than extension of a path.
   dp[S][v]=sum(dp[prev][u] for u in range(n) if prev>>u&1 and rows[u]>>v&1)
   start[S][v]=sum(start[prev][u] for u in range(n) if prev>>u&1 and rows[u]>>v&1)
 H=[sum(row) for row in dp];H[0]=1
 R=sum(sum(start[S])*H[(N-1)^S] for S in range(1,N) if S>>root&1)
 C=sum(start[-1][u] for u in range(n) if rows[u]>>root&1)
 return H[-1],dp[-1][root],R,sum(start[-1]),C

@lru_cache(None)
def compositions(total):
 if total==0:return ((),)
 return tuple((d,)+rest for d in range(2,total+1) for rest in compositions(total-d))

def palindromes(total):
 for leftsum in range(total//2+1):
  center=total-2*leftsum
  if center==1:continue
  for left in compositions(leftsum):
   yield left+((center,) if center else ())+left[::-1]

def busch(n):
 if n%3==0:return 3*5**((n-3)//3)
 if n%3==1:return 5**((n-1)//3)
 return 9*5**((n-5)//3)

def sharp_gaps(n):
 if n%3==1:return (3,)*((n-1)//3)
 if n%3==2:return (2,)+(3,)*((n-5)//3)+(2,)
 if n%6==3:
  t=(n-3)//6;return (3,)*t+(2,)+(3,)*t
 if n==6:return (5,)
 t=n//6;return (3,)*(t-2)+(2,2,3,2,2)+(3,)*(t-2)

def main():
 physical_checks=palphysical=0
 for total in range(2,12):
  for gaps in compositions(total):
   rows=graph(gaps);validate(rows);h,e,R,s,C=physical(rows)
   require((h,R,s,C)==(predicted(gaps),predicted(gaps),1,1),('literal profile',gaps,(h,e,R,s,C)))
   # The final gap has E=1, so its endpoint readout equals the previous H.
   require(e==predicted(gaps[:-1]),('endpoint',gaps))
   physical_checks+=1
   if gaps==gaps[::-1]:require(reflect(rows),('reflection',gaps));palphysical+=1
 maxn=300;best=[None]*maxn;best[0]=1
 for total in range(2,maxn):
  candidates=[factor(d)*best[total-d] for d in range(2,total+1) if best[total-d] is not None]
  if candidates:best[total]=min(candidates)
 minima={};explicit=[];inventory=0;reflection_checks=0
 for n in range(3,maxn+1):
  total=n-1;candidates=[]
  for center in range(total+1):
   if center==1 or (total-center)%2:continue
   side=best[(total-center)//2]
   if side is not None:candidates.append(side*side*(factor(center) if center else 1))
  actual=min(candidates);gaps=sharp_gaps(n);want=predicted(gaps)
  require(sum(gaps)==total and gaps==gaps[::-1],('construction size',n))
  require(actual==want,('family minimum',n,actual,want))
  require(reflect(graph(gaps)),('large reflection',n));reflection_checks+=1
  if n%6:
   require(want==busch(n),('Busch attainment',n))
  elif n==6:require(want==17 and busch(n)+2==17,'six')
  else:
   t=n//6;require(want==405*25**(t-2),'multiple6 upper')
   require(25*want==27*busch(n),'ratio')
   require(want>=busch(n)+2,'rounded lower')
  if n<=32:
   words=list(palindromes(total));inventory+=len(words)
   require(min(map(predicted,words))==actual,('exhaustive minimum',n))
  if n<=18:
   explicit.append(dict(order=n,gaps=gaps,H=want,C=1,rows=graph(gaps),Busch=busch(n),self_converse_lower=busch(n)+(2 if n%6==0 else 0)))
  if n%6==0:minima[str(n)]=actual
 result=dict(status='PASS',literal_gap_compositions_checked=physical_checks,literal_max_order=12,literal_palindromic_reflections=palphysical,formula_and_family_DP_orders=[3,maxn],explicit_reflections=reflection_checks,exhaustive_palindromic_inventory_max_order=32,exhaustive_palindromic_lists=inventory,multiple6_family_minima=minima,small_witnesses=explicit)
 data=json.dumps(result,sort_keys=True,indent=2)+'\n';(ROOT/'self_converse_seams_certificate.json').write_text(data);print(data,end='')
if __name__=='__main__':main()
