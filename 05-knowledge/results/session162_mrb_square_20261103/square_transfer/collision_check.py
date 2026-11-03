from math import isqrt
A=[3,22,14,2,23,13,12,4,21,15,10,6,19,17,8,1,24,25,11,5,20,16,9,7,18]
B=[18,7,9,16,20,5,11,25,24,1,8,17,19,6,10,15,21,4,12,13,3,22,14,2,23]
def make(a):
 if set(a)!=set(range(1,26)): raise AssertionError("not a permutation")
 if not all(isqrt(x+y)**2==x+y for x,y in zip(a,a[1:])): raise AssertionError("non-square edge")
 es=[];prev=5
 for x,y in zip(a,a[1:]):
  d=abs(x-y);s=1
  while s*s%2!=d%2 or s*s<2*(prev+1)+d:s+=1
  u=(s*s-d)//2;v=u+d;es.append((u,v));prev=v
 return es
x,y=make(A),make(B)
if A == B or A == B[::-1]: raise AssertionError("same path")
if sorted(abs(a-b) for a,b in zip(A,A[1:])) != sorted(abs(a-b) for a,b in zip(B,B[1:])): raise AssertionError("different gaps")
if set(x) == set(y): raise AssertionError("identical interval edge sets")
if (x[-1][-1],y[-1][-1]) != (846,623): raise AssertionError("unexpected interval end")
print('different unoriented',A!=B[::-1]);print('same gap multiset',sorted(abs(a-b) for a,b in zip(A,A[1:]))==sorted(abs(a-b) for a,b in zip(B,B[1:])));print('ends',x[-1][-1],y[-1][-1]); print('edge sets equal',set(x)==set(y));print('first endpoints',x[:3],y[:3]);print('N',max(x[-1][-1],y[-1][-1]))
from collections import Counter
D=[abs(a-b) for a,b in zip(A,A[1:])]
g=846-5-sum(d+1 for d in D)
P=1; R=2**g
for d in D:
 P*=1 if d==1 else 1+2**(d-1)
 R*=2**(d+1)+2*(3**(d-1)-1)
print('gap multiset',sorted(Counter(D).items()), 'g', g)
print('P',P)
print('R=',f'2^{g}', '*', ' * '.join(f'{2**(d+1)+2*(3**(d-1)-1)}^{k}' for d,k in sorted(Counter(D).items())))
