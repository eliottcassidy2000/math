import json
from math import prod,gcd
from itertools import combinations,product
from fractions import Fraction

def need(c,m):
 if not c:raise RuntimeError(m)
def val(n,p):
 r=0
 while n%p==0:n//=p;r+=1
 return r
def det(A):
 if len(A)==1:return A[0][0]
 return sum((-1)**j*A[0][j]*det([r[:j]+r[j+1:] for r in A[1:]]) for j in range(len(A)))
def smith3(A):
 g1=0;g2=0
 for r in A:
  for x in r:g1=gcd(g1,x)
 for rs in combinations(range(3),2):
  for cs in combinations(range(3),2):g2=gcd(g2,abs(det([[A[i][j] for j in cs] for i in rs])))
 return [g1,g2//g1,abs(det(A))//g2]
def matrix(ns,ps):return [[val(n,p) for n in ns] for p in ps]
A=matrix([101,35,175],[101,5,7]);B=matrix([91,175,325],[5,7,13])
need(smith3(A)==[1,1,1],'unimodular');need(smith3(B)==[1,1,4],'index4')
need(Fraction(175,35)==5 and Fraction(35**2,175)==7,'basis inverse')
need(91*175*325==2275**2,'square relation')
need(Fraction(175*325,91)==5**4,'fourth power')
need(Fraction(91*175,325)==7**2 and Fraction(91*325,175)==13**2,'squares')
checks=6
for e in product(range(-5,6),repeat=3):
 s=[sum(B[i][j]*e[j] for j in range(3)) for i in range(3)]
 need((s[0]-2*s[1]-2*s[2])%4==0,'lattice condition');checks+=1
for s in product(range(-7,8),repeat=3):
 if (s[0]-2*s[1]-2*s[2])%4:continue
 y=Fraction(s[0]+2*s[1]-2*s[2],4);z=Fraction(s[0]-2*s[1]+2*s[2],4);x=Fraction(2*s[1]+2*s[2]-s[0],4)
 need(all(t.denominator==1 for t in (x,y,z)),'inverse integral')
 need([sum(B[i][j]*[x,y,z][j] for j in range(3)) for i in range(3)]==list(s),'inverse');checks+=2
for ns,ps,expected in [([101,35,175],[101,5,7],[0]),([91,175,325],[5,7,13],[0,7])]:
 masks=[]
 for mask in range(8):
  n=prod(x for i,x in enumerate(ns) if mask>>i&1)
  if all(val(n,p)%2==0 for p in ps):masks.append(mask)
 need(masks==expected,'binary kernel');checks+=1
print(json.dumps({'status':'FINITE-EXACT','checks':checks,'focal_matrix':A,'focal_smith':smith3(A),'orbit_triangle_matrix':B,'orbit_triangle_smith':smith3(B),'triangle_square_product':2275**2,'scope':'valuation lattices, not Collatz convergence or tournament isomorphism'},indent=2,sort_keys=True))
