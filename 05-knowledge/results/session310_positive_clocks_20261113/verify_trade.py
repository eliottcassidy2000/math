from collections import Counter
from itertools import permutations
def typ(rows):
 n=len(rows);cs=[]
 for mask in range(1,1<<n):
  s=mask.bit_count()
  if s<3 or s%2==0:continue
  ar=[i for i in range(n) if mask>>i&1]
  for p in permutations(ar[1:]):
   w=[ar[0],*p,ar[0]]
   if all(rows[u]>>v&1 for u,v in zip(w,w[1:])):cs.append((mask,s))
 d=Counter()
 def rec(i,mask,sizes):
  d[tuple(sorted(sizes))]+=1
  for j in range(i,len(cs)):
   x,l=cs[j]
   if not mask&x:rec(j+1,mask|x,sizes+[l])
 rec(0,0,[]);return d
A=[2,28,25,81,97,79,7]
B=[2,80,3,7,173,143,189,15]
if typ(A)!={():1,(3,):10,(5,):18,(7,):7,(3,3):6} or typ(B)!={():1,(3,):10,(5,):30,(7,):7}:raise ValueError("unexpected cycle trade")
print("PASS: literal directed cycle trade")
