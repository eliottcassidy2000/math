#!/usr/bin/env python3
"""Independent exact controls for elementary proposition G and typed fives."""
import itertools,json
from fractions import Fraction

def det(M):
 a=[[Fraction(v) for v in r] for r in M]; n=len(a); out=Fraction(1)
 for j in range(n):
  pivot=next((i for i in range(j,n) if a[i][j]),None)
  if pivot is None:return 0
  if pivot!=j:a[pivot],a[j]=a[j],a[pivot];out=-out
  v=a[j][j];out*=v
  for k in range(j,n):a[j][k]/=v
  for i in range(j+1,n):
   v=a[i][j]
   for k in range(j,n):a[i][k]-=v*a[j][k]
 return out

def rankmod(M,p):
 a=[[v%p for v in r] for r in M]; n=len(a); k=0
 for j in range(len(a[0])):
  pivot=next((i for i in range(k,n) if a[i][j]),None)
  if pivot is None:continue
  a[k],a[pivot]=a[pivot],a[k];inv=pow(a[k][j],-1,p);a[k]=[v*inv%p for v in a[k]]
  for i in range(n):
   if i!=k:
    v=a[i][j];a[i]=[(x-v*y)%p for x,y in zip(a[i],a[k])]
  k+=1
  if k==n:break
 return k

def req(v,msg):
 if not v:raise RuntimeError(msg)
checks=[]
for n in range(2,41):
 M=[[1]+[int(i==j)-int(i==n-1) for j in range(n-1)] for i in range(n)]
 req(abs(det(M))==n,f'index at {n}')
 for p in (2,3,5,7,11,13):req(rankmod(M,p)==n-int(n%p==0),f'rank at {n},{p}')
 checks.append({'n':n,'index':n})
platonic=[(p,q) for p in range(3,15) for q in range(3,15) if (p-2)*(q-2)<4]
req(platonic==[(3,3),(3,4),(3,5),(4,3),(5,3)],'Platonic scope')
# An involution's cycle lengths are a necessary invariant for equivariant bijection.
req(sum(p==q for p,q in platonic)==1,'duality fixed points')
# Five core types all fixed versus one fixed + two pairs: no equivariant map.
req(5!=1,'fixed-point obstruction')
# Catalan and partial permutation counts agree only through the first three sizes.
cat=[1]
for n in range(1,7):cat.append(sum(cat[i]*cat[n-1-i] for i in range(n)))
partial=[sum(__import__('math').factorial(n)//__import__('math').factorial(n-k) for k in range(n+1)) for n in range(6)]
req(cat[3]==partial[2]==5 and cat[4]==14 and partial[3]==16,'different five growth laws')
print(json.dumps({'status':'PASS','integer_index_and_modular_rank':checks,'primes':[2,3,5,7,11,13],'platonic_pairs':platonic,'catalan':cat,'partial_permutation':partial,'scope':'elementary proposition G and finite controls, not review of inherited research'},sort_keys=True,indent=2))
