#!/usr/bin/env python3
from itertools import permutations,combinations,product
from math import gcd,isqrt,prod
from functools import reduce
import json
checks=0
def check(x):
 global checks
 checks+=1
 if not x: raise RuntimeError('independent lattice audit failed')
def det(m):
 n=len(m)
 return sum((-1)**sum(p[i]>p[j] for i in range(n) for j in range(i+1,n))*prod(m[i][p[i]] for i in range(n)) for p in permutations(range(n)))
for nums,M,expect in [([101,35,175],[[1,0,0],[0,1,2],[0,1,1]],[1,1,1]),([91,175,325],[[0,2,2],[1,1,0],[1,0,1]],[1,1,4])]:
 ds=[1]
 for n in [1,2,3]:
  minors=[abs(det([[M[i][j] for j in J] for i in I])) for I in combinations(range(3),n) for J in combinations(range(3),n)]
  ds.append(reduce(gcd,minors))
 check([ds[i+1]//ds[i] for i in range(3)]==expect)
 square_masks=[]
 for mask in range(8):
  q=prod(x for i,x in enumerate(nums) if mask>>i&1)
  if isqrt(q)**2==q:square_masks.append(mask)
 check(square_masks==([0] if nums[0]==101 else [0,7]))
# The quotient map Z^3 -> Z/4 is surjective and kills each old column.
M=[[0,2,2],[1,1,0],[1,0,1]]
for j in range(3):check((M[0][j]-2*M[1][j]-2*M[2][j])%4==0)
check({(a-2*b-2*c)%4 for a,b,c in product(range(4),repeat=3)}==set(range(4)))
# Its kernel has index4, which matches |det M|, proving equality with image.
check(abs(det(M))==4)
print(json.dumps({'status':'PASS','checks':checks,'smith_focal':[1,1,1],'smith_inherited':[1,1,4],'independent_method':'permutation determinants; exact integer-square tests; quotient-map index'},indent=2,sort_keys=True))
