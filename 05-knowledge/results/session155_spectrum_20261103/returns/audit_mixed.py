from collections import deque,Counter
from sympy import Matrix
qs=[(240,339,90,127),(720,1101,270,413),(1020432,1560769,366879,561148),(211965,378780,76220,136205),(5,8,2,3),(11,18,4,7),(4773,8628,1789,3234)]
def prod(a,b,p):return ((a[0]*b[0]+a[1]*b[2])%p,(a[0]*b[1]+a[1]*b[3])%p,(a[2]*b[0]+a[3]*b[2])%p,(a[2]*b[1]+a[3]*b[3])%p)
def gate(A,p,k):
 a,b,c,d=(x%p for x in A)
 if (a*d-b*c)%p:return True
 return (a!=0 and b==0) or (c!=0 and d==0) if k==0 else ((b!=0 and a==0) or (d!=0 and c==0))
S={(1,0,0,1)};todo=deque(S)
while todo:
 a=todo.popleft()
 for q in qs:
  c=prod(a,q,30)
  if c not in S:S.add(c);todo.append(c)
assert all(prod(a,q,30) in S for a in S for q in qs)
z=[a for a in S if all(v%3==0 for v in a)]
profile=Counter((gate(a,2,0),gate(a,5,1)) for a in z)
print('closure',len(S),'zero3',len(z),'gates',sorted(profile.items()))
assert len(S)==6400 and len(z)==184 and profile=={(False,True):12,(True,False):111,(False,False):61}
r24=Matrix([[5,16,16],[18,28,30],[18,22,24]]);a=Matrix([[1,0,1],[2,0,1],[0,1,0]]);b=Matrix([[1,2,2],[4,3,4],[4,2,3]])
L2=r24*a*a*a;D=L2*b*L2;W=b*D
v=Matrix([[31,14,25]])*D;w=Matrix([[69,0,-21]]);t=Matrix([[2,-1,0]])
assert w*W==3*v+216*t
beta=Matrix([[0,1],[0,2],[1,-1]]); Qb=Matrix([[5,8],[2,3]]);QL=Matrix([[240,339],[90,127]])
R=Qb*QL*Qb*QL/3
assert all(x.q==1 for x in R) and R.det()==100
pred=v*beta*R.inv();assert pred==Matrix([[-21,90]])
print('P_D',list(D),'P_W',list(W));print('descendant',list(v),'unique rational predecessor bulk',list(pred),'carry',216)
print('PASS')
