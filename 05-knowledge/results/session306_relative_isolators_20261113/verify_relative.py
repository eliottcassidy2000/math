from itertools import combinations, permutations
from fractions import Fraction
from probe import decode,profile

def coeff(q,v=0):
 n=len(q); sides=[u for u in range(n) if u!=v]
 counts=[0]*(n//2+2)
 def rec(todo,w,j):
  if not todo and j:counts[j]+=1
  if not w or w[-1]!=v:
   if not w or q[w[-1]][v]:
    rec(todo,w+[v],j+1)
  for u in todo:
   if not w or q[w[-1]][u]:rec([x for x in todo if x!=u],w+[u],j)
 rec(sides,[],0)
 return counts
def cyc(q):
 ss=[(i,j) for i,n in enumerate(q) for j in range(n)]
 return [[int((k-i)%3==1) if i!=k else int(j<l) for k,l in ss] for i,j in ss]
def require(x):
 if not x: raise ValueError('identity failed')
for n in range(1,6):
 for c in range(1<<(n*(n-1)//2)):
  q=decode(n,c); k=coeff(q)
  for j,a in enumerate(k):
   if j>(n-1)//2+1:require(not a)
  require(k[1]==profile(q)[1])
  left=[i for i in range(1,n) if q[0][i]]
  right=[i for i in range(1,n) if q[i][0]]
  mu=max(j for j in range(1+min(len(left),len(right)))
         if any(all(q[x][y] for x,y in zip(u,w))
                for u in combinations(left,j) for w in permutations(right,j)))
  require([j for j,a in enumerate(k) if a]==list(range(1,mu+2)))
for r in range(2,6):
 k=coeff(cyc([1,r-1,r-1]))
 require(k[r]>0 and not any(k[r+1:]))
require(coeff(cyc([1,2,2]))==[0,11,17,4])
P=[profile(decode(7,c)) for c in (245600,168491,1104621)]
require([a[1:4] for a in P]==[[29,522,3402],[33,594,3714],[37,666,4026]])
require(len({Fraction(a[2],a[1]) for a in P})==1)
require(len({Fraction(a[3],a[1]) for a in P})==3)
print('PASS sharp probe, all 1099 quotient controls and fixed-order collision')
