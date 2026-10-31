"""Independent exhaustive layer-certificate checks, including zero coordinates."""
from itertools import combinations, product
from functools import lru_cache
import json

def require(t,m):
 if not t:raise RuntimeError(m)
@lru_cache(None)
def G(a,b):return sum(1<<(i&j).bit_count() for i in range(a) for j in range(b))
def patterns(n):
 out=set()
 for u in range(n//2+1):
  for v in range(u,n-u+1):
   p=(u,n-u,u,v,n-v)
   for k in range(5):out.add(p[k:]+p[:k])
 return sorted(out)
def layers(A,B,a,b):
 P=[{i,(i+2)%5} for i in range(5)]
 Q=[{i,(i+1)%5} for i in range(5)]
 S=[i for i in range(5) if a[i]%2 and b[i]%2]
 exceptional=False
 if A%2 and B%2:P=[{0,1}]*5;Q=[{0,2}]*5
 elif A%2:
  if len(S)==5:exceptional=True
  else:
   cover=next(set(c) for c in combinations(range(5),2) if all(set(c)&Q[i] for i in S))
   P=[cover]*5
 elif B%2:
  if len(S)==5:exceptional=True
  else:
   cover=next(set(c) for c in combinations(range(5),2) if all(set(c)&P[i] for i in S))
   Q=[cover]*5
 aa=[[a[i]//2+int(a[i]%2 and j in P[i]) for i in range(5)] for j in range(5)]
 bb=[[b[i]//2+int(b[i]%2 and j in Q[i]) for i in range(5)] for j in range(5)]
 budgets=[]
 for j in range(5):
  ac=max(aa[j][i]+aa[j][(i+1)%5] for i in range(5))
  bc=max(bb[j][i]+bb[j][(i+2)%5] for i in range(5))
  budgets.append((ac,bc))
 for i in range(5):require(sum(G(aa[j][i],bb[j][i]) for j in range(5))>=G(a[i],b[i]),('local',A,B,a,b,i))
 require(sum(G(x,y) for x,y in budgets)<=G(A,B),('global',A,B,a,b,budgets,exceptional))
 return exceptional
counts={'total':0,'all_double_odd_exception':0,'with_zero':0}
for A in range(2,11):
 for B in range(2,11):
  for a in patterns(A):
   for t in patterns(B):
    b=[0]*5
    for j in range(5):b[(2*j)%5]=t[j]
    exceptional=layers(A,B,a,b)
    counts['total']+=1;counts['all_double_odd_exception']+=exceptional;counts['with_zero']+=int(0 in a or 0 in b)
# Maximal odd-budget patterns cannot have all coordinates odd: explicitly test
# the nonmaximal exceptional family, rather than leaving its count at zero.
for A in range(2,9):
 for B in range(2,9):
  if A%2==B%2:continue
  ap=[p for p in product(range(1,A,2),repeat=5) if all(p[i]+p[(i+1)%5]<=A for i in range(5))]
  bp=[p for p in product(range(1,B,2),repeat=5) if all(p[i]+p[(i+2)%5]<=B for i in range(5))]
  for a in ap:
   for b in bp:
    require(layers(A,B,a,b),('exception branch',A,B,a,b))
    counts['total']+=1;counts['all_double_odd_exception']+=1
require(counts['all_double_odd_exception']>0,'missing exception controls')
print(json.dumps({'status':'PASS','budgets':'2..10','counts':counts},sort_keys=True))
