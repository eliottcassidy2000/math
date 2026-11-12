#!/usr/bin/env python3
"""Independent matrix controls, no producer imports."""
from itertools import combinations,permutations,product
from collections import Counter
from fractions import Fraction
from math import comb,factorial
import json
from functools import lru_cache

def check(x,msg):
 if not x:raise RuntimeError(msg)
def matrix(n,bits):
 a=[[0]*n for _ in range(n)]
 for z,(i,j) in enumerate(combinations(range(n),2)):a[i][j]=(bits>>z)&1;a[j][i]=1-a[i][j]
 return a
def homogeneous_pairs(a):
 return [(i,j) for i,j in combinations(range(len(a)),2) if all(a[i][k]==a[j][k] for k in range(len(a)) if k not in (i,j))]
def prime(a):
 n=len(a)
 for k in range(2,n):
  for s in combinations(range(n),k):
   if all(len({a[v][x] for x in s})==1 for v in range(n) if v not in s):return False
 return n>=3

def code(a,p):return tuple(a[p[i]][p[j]] for i,j in combinations(range(len(a)),2))
def canonical(a):return min(code(a,p) for p in permutations(range(len(a))))
def add(a,z):
 n=len(a);r=[row+[1-((z>>i)&1)] for i,row in enumerate(a)];r.append([(z>>i)&1 for i in range(n)]+[0]);return r

@lru_cache(None)
def inD(a):
 n=len(a)
 if n==1:return True
 def restrict(vs):return tuple(tuple(a[i][j] for j in vs) for i in vs)
 for mask in range(1,(1<<n)-1):
  left=[i for i in range(n) if mask>>i&1];right=[i for i in range(n) if not(mask>>i&1)]
  if all(a[i][j] for i in left for j in right) and inD(restrict(left)) and inD(restrict(right)):return True
 for tail in product(range(3),repeat=n-1):
  colors=(0,)+tail;parts=[[i for i in range(n) if colors[i]==c] for c in range(3)]
  if all(parts) and all(a[i][j] for c in range(3) for i in parts[c] for j in parts[(c+1)%3]) and all(inD(restrict(p)) for p in parts):return True
 return False

def main():
 # Every fixed labelled parent n<=4, k=2, not only a pooled ensemble.
 cases=0
 for n in range(1,5):
  k=2;N=n+k
  for c in range(1<<comb(n,2)):
   a=matrix(n,c);old=homogeneous_pairs(a);m=len(old);surv=Counter();sumrank=0
   for bits in range(1<<(n*k+1)):
    r=[row+[0]*k for row in a]+[[0]*N for _ in range(k)];p=0
    for x in range(n,N):
     for v in range(n):r[x][v]=(bits>>p)&1;r[v][x]=1-r[x][v];p+=1
    r[n][n+1]=(bits>>p)&1;r[n+1][n]=1-r[n][n+1]
    edges=homogeneous_pairs(r);surv[sum(e in edges for e in old)]+=1;sumrank+=N-len(edges);cases+=1
   den=1<<(n*k+1)
   for s in range(m+1):check(Fraction(surv[s],den)==Fraction(comb(m,s)*3**(m-s),4**m),'fixed-parent survival')
   check(Fraction(sumrank,den)==N-Fraction(m,4)-Fraction(n*k+1,2**(N-2)),'fixed-parent mean')
 # W5 exact prime attachment count and forgetting-root multiplicities.
 w='1010';q=len(w);a=[[int(i<j) for j in range(q)] for i in range(q)]
 a=[row+[1-int(w[i])] for i,row in enumerate(a)]+[[int(x) for x in w]+[0]]
 check(prime(a),'W5 prime');aut=[p for p in permutations(range(5)) if code(a,p)==code(a,tuple(range(5)))];check(len(aut)==1,'W5 rigid')
 fibres=Counter(canonical(add(a,z)) for z in range(32) if prime(add(a,z)))
 check(sum(fibres.values())==20,'prime augmentation');check(max(fibres.values())<=6,'root multiplicity')
 # Explicit 81-element permutation group of the nine-vertex ternary tower.
 a=[[0]*9 for _ in range(9)]
 for i,j in combinations(range(9),2):
  x,y=(i//3,j//3) if i//3!=j//3 else (i%3,j%3)
  a[i][j]=int((y-x)%3==1);a[j][i]=1-a[i][j]
 group=[]
 for root in range(3):
  for rs in product(range(3),repeat=3):
   p=tuple(3*((i//3+root)%3)+(i%3+rs[i//3])%3 for i in range(9))
   check(code(a,p)==code(a,tuple(range(9))),'tower automorphism');group.append(p)
 check(len(set(group))==81,'tower group order')
 floors={s:81 for s in range(10)}; colororders=Counter()
 for mask in range(512):
  pts=[i for i in range(9) if mask>>i&1]
  pw=sum(all(p[i]==i for i in pts) for p in group)
  floors[len(pts)]=min(floors[len(pts)],pw)
  co=sum(all(((mask>>p[i])&1)==((mask>>i)&1) for i in range(9)) for p in group)
  colororders[co]+=1
 for s in floors:check(floors[s]==3**(max(1-s,0)+max(3-s,0)),'pointwise group floor')
 # Weighted root-forgetting identity, computed from literal automorphisms.
 forgetting=[]
 for n,k in [(1,2),(2,2),(3,1),(3,2)]:
  reps={}
  for c in range(1<<comb(n,2)):
   a=matrix(n,c);reps.setdefault(canonical(a),a)
  total=Fraction(0);N=n+k
  for a in reps.values():
   ag=sum(code(a,p)==code(a,tuple(range(n))) for p in permutations(range(n)))
   for bits in range(1<<(n*k+comb(k,2))):
    r=[row+[0]*k for row in a]+[[0]*N for _ in range(k)];p=0
    for i,j in combinations(range(N),2):
     if j<n:continue
     r[i][j]=(bits>>p)&1;r[j][i]=1-r[i][j];p+=1
    ug=sum(code(r,q)==code(r,tuple(range(N))) for q in permutations(range(N)))
    total+=Fraction(ug,ag*factorial(k))
  answer=total/comb(N,k);check(answer=={3:2,4:4,5:12}[N],'weighted forgetting')
  forgetting.append({'parent_order':n,'batch':k,'classes':int(answer)})
 intersection=0
 for m in range(7):
  for word in product((0,1),repeat=m):
   a=[[int(i<j) for j in range(m)] for i in range(m)]
   a=[row+[1-word[i]] for i,row in enumerate(a)]+[list(word)+[0]]
   d=sum(word[i:i+2]==(1,0) for i in range(max(0,m-1)))
   check(inD(tuple(map(tuple,a)))==(d<=1),'independent actual D intersection');intersection+=1
 return {'status':'PASS','actual_intersection_words':intersection,'weighted_root_forgetting':forgetting,'nine_tower_pointwise_floors':floors,'nine_tower_coloring_group_orders':dict(colororders),'fixed_parent_batch_cases':cases,'W5_prime_attachments':20,'W5_distinct_prime_children':len(fibres),'W5_prime_fibre_sizes':sorted(fibres.values())}
if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
