from fractions import Fraction as F
from itertools import permutations
from math import factorial
import json

def sign(p):
 return (-1)**sum(p[i]>p[j] for i in range(len(p)) for j in range(i+1,len(p)))
def det(M):
 n=len(M); s=F(0)
 for p in permutations(range(n)):
  v=F(sign(p))
  for i,j in enumerate(p):v*=M[i][j]
  s+=v
 return s
def pf(M):
 n=len(M); s=F(0)
 for p in permutations(range(n)):
  v=F(sign(p))
  for i in range(0,n,2):v*=M[p[i]][p[i+1]]
  s+=v
 return s/(2**(n//2)*factorial(n//2))
def matrix(x):
 M=[[F(0) for j in range(6)] for i in range(6)]
 es={(0,1):F(x)}
 es.update(dict(zip([(2,3),(2,4),(2,5),(3,4),(3,5),(4,5)],[1,1,1,1,-1,1])))
 es.update(dict(zip([(0,2),(0,3),(0,4),(0,5),(1,2),(1,3),(1,4),(1,5)],[-1,-1,-1,-1,1,-1,1,1])))
 for (i,j),v in es.items():M[i][j]=F(v);M[j][i]=-F(v)
 return M
checks=0
for x in [F(-1),F(0),F(2,3),F(1)]:
 K=matrix(x)
 if pf(K)!=3*x-2:raise RuntimeError('Pfaffian polynomial')
 if det(K)!=pf(K)**2:raise RuntimeError('determinant square')
 checks+=2
K=matrix(F(2,3)); A=[[F(0) if i==j else (1+K[i][j])/2 for j in range(6)] for i in range(6)]
u=det([[A[i][j]+(i==j) for j in range(6)] for i in range(6)]);v=det(A)
if (u,v,u-abs(v))!=(F(32,9),F(-5,3),F(17,9)):raise RuntimeError('margin')
if pf([row[2:] for row in K[2:]])!=3:raise RuntimeError('complement')
if pf(K)!=0 or not F(0)<F(4,9):raise RuntimeError('hostile')
checks+=3
print(json.dumps({'status':'PASS','checks':checks,'method':'full antisymmetrization and Leibniz determinants','u':str(u),'v':str(v),'margin':str(u-abs(v))},sort_keys=True))
