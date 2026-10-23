#!/usr/bin/env python3
"""Independent integer-matrix, tensor-basis and all ternary-valued 2D tables."""
from itertools import product
import json
N=0

def require(x,label):
 global N
 N+=1
 if not x:raise RuntimeError(label)

def mm(A,B):return [[sum(x*y for x,y in zip(row,col)) for col in zip(*B)] for row in A]
def sub(A,B):return [[x-y for x,y in zip(a,b)] for a,b in zip(A,B)]
def scale(a,A):return [[a*x for x in row] for row in A]
def power(A,n):
 R=[[int(i==j) for j in range(len(A))] for i in range(len(A))]
 for _ in range(n):R=mm(R,A)
 return R

def tensor(A,B):return [[a*b for a in rowA for b in rowB] for rowA in A for rowB in B]
def v3(x):
 if x==0:return 10**6
 n=0
 while x%3==0:x//=3;n+=1
 return n
E=[[1,0,0],[1,1,0],[1,2,1]];I=[[1,0,0],[0,1,0],[0,0,1]];F=[[1,0,0],[-1,1,0],[1,-2,1]]
U=[[0,1,0],[0,0,1],[1,0,0]];T=mm(mm(F,U),E);D=sub(T,I)
require(mm(E,F)==I and mm(F,E)==I,'integral inverse matrices')
require(T==[[1,1,0],[0,1,1],[0,-3,-2]],'cyclic wrap derived by conjugacy')
require(power(T,3)==I and mm(T,D)==mm(D,T),'order three and commuting difference')
require(power(D,3)==scale(-3,mm(T,D)),'integral cubic difference identity')
filtration=[]
for k in range(1,13):
 even=power(D,2*k);odd=power(D,2*k+1)
 require(even==scale((-3)**(k-1),mm(power(T,k-1),power(D,2))),'exact even exponent identity')
 require(odd==scale((-3)**k,mm(power(T,k),D)),'exact odd exponent identity')
 require(min(v3(x) for row in even for x in row)==k-1,'sharp lower divisibility in Newton basis')
 require(min(v3(x) for row in odd for x in row)==k,'sharp upper divisibility')
 m=3**k
 require(any(x%m for row in even for x in row) and not any(x%m for row in odd for x in row),'exact nilpotency index')
 filtration.append(dict(modulus=m,index=2*k+1))
A=[[1]];B=[[1]];tensor_sizes=[]
for r in range(5):
 size=3**r;want=[[int(i==j) for j in range(size)] for i in range(size)]
 require(mm(A,B)==want and mm(B,A)==want,'lossless tensor basis over integers')
 tensor_sizes.append(size);A=tensor(A,E);B=tensor(B,F)
Z=[[-1,1],[-1,-2]]
require(Z[0][0]*Z[1][1]-Z[0][1]*Z[1][0]==3,'zero-sum difference determinant')
for m in (2,4,5,8,10,16):
 image={tuple(sum(a*b for a,b in zip(row,v))%m for row in Z) for v in product(range(m),repeat=2)}
 require(len(image)==m*m,'zero-sum difference is bijective whenever 3 is a unit')
points=list(product(range(3),repeat=2));index={a:i for i,a in enumerate(points)}
shift_maps={v:[index[((a[0]+v[0])%3,(a[1]+v[1])%3)] for a in points] for v in points}
coefficient_matrix=tensor(F,F)
linear_maps={u:[(u[0]*a[0]+u[1]*a[1])%3 for a in points] for u in points if u!=(0,0)}
rank_hist={0:0,1:0,2:0}
for vals in product(range(3),repeat=9):
 coefficients=[sum(a*b for a,b in zip(row,vals))%9 for row in coefficient_matrix]
 for axis in (0,1):
  shift=(int(axis==0),int(axis==1));essential=any(vals[i]!=vals[j] for i,j in enumerate(shift_maps[shift]))
  require(essential==any(c!=0 and a[axis]>0 for a,c in zip(points,coefficients)),'essential axis from Newton coefficients')
 K=[v for v in points if all(vals[i]==vals[j] for i,j in enumerate(shift_maps[v]))]
 require(len(K) in (1,3,9),'literal translation kernel subgroup size')
 predicted={1:2,3:1,9:0}[len(K)]
 if len(set(vals))==1:actual=0
 else:
  sufficient=False
  for word in linear_maps.values():
   value_sets=[{v for v,b in zip(vals,word) if b==a} for a in range(3)]
   if all(len(s)==1 for s in value_sets):sufficient=True;break
  actual=1 if sufficient else 2
 require(actual==predicted,'minimum linear observation rank by independent exhaustive row maps')
 rank_hist[actual]+=1
hostile=tuple(sum(a)%3 for a in points)
coeff=[sum(a*b for a,b in zip(row,hostile))%9 for row in coefficient_matrix]
K=[v for v in points if all(hostile[i]==hostile[j] for i,j in enumerate(shift_maps[v]))]
require(K==[(0,0),(1,2),(2,1)],'hostile literal kernel')
require(coeff[index[1,0]] and coeff[index[0,1]],'hostile Newton indices rank two')
print(json.dumps(dict(status='PASS',checks=N,tensor_dimensions=tensor_sizes,sharp_filtration=filtration,
                     ternary_valued_2d_tables=3**9,minimum_rank_histogram=rank_hist,
                     hostile=dict(kernel=K,coefficients={str(a):c for a,c in zip(points,coeff)})),sort_keys=True,indent=2))
