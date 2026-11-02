from fractions import Fraction as F
from itertools import permutations,combinations
import sympy as S
n=5
w=S.Matrix([1,1,-2,0,0]); Q=S.eye(n)-2*w*w.T/(w.dot(w)); J=S.ones(n)
K=S.zeros(n)
K[0,3]=S.Rational(1,4);K[3,0]=-K[0,3]
K[1,4]=S.Rational(1,2);K[4,1]=-K[1,4]
A=(J-S.eye(n)+K)/2; B=Q*A*Q.T
z=S.symbols('z')
def inv(A):
 H=sum(S.prod(A[p[i],p[i+1]] for i in range(4)) for p in permutations(range(5)))
 V=sum(A[i,j]*A[j,i] for i,j in combinations(range(5),2))
 tau=sum(A[i,j]*A[j,k]*A[k,i]+A[i,k]*A[k,j]*A[j,i] for i,j,k in combinations(range(5),3))
 M23=0
 for i,j in combinations(range(5),2):
  a,b,c=[x for x in range(5) if x not in (i,j)]
  M23+=A[i,j]*A[j,i]*(A[a,b]*A[b,c]*A[c,a]+A[a,c]*A[c,b]*A[b,a])
 u=(S.eye(n)+A).det();v=-A.det();delta=H-u-abs(v)
 return {'H':H,'V':V,'tau':tau,'M23':M23,'u':u,'v':v,'delta':delta,'P':S.expand((S.eye(n)-z*A).det())}
print('Q',Q);print('A',A);print('B',B)
print('A invariants',inv(A));print('B invariants',inv(B))
