from fractions import Fraction as F
from math import comb, prod
import json
import sympy as s
N=12
def mul(a,b):return [sum(a[j]*b[i-j] for j in range(i+1)) for i in range(N+1)]
def pw(a,n):
 r=[1]+[0]*N
 for _ in range(n):r=mul(r,a)
 return r
E=[1]+[240*sum(d**3 for d in range(1,k+1) if k%d==0) for k in range(1,N+1)]
D=[0,1]+[0]*(N-1)
for k in range(1,N+1):
 p=[0]*(N+1)
 for j in range(min(24,N//k)+1):p[k*j]=(-1)**j*comb(24,j)
 D=mul(D,p)
b=[mul(pw(E,12-3*j),pw(D,j)) for j in range(5)]
c=[1]
for i in range(1,5):c.append(-sum(c[j]*b[j][i] for j in range(i)))
theta=[sum(c[j]*b[j][i] for j in range(5)) for i in range(N+1)]
def moments(f,k):
 if f=='A':return [k*(k+1),-k*(k+1)*(2*k-1),F(4,3)*k*(k+1)*(5*k*k-7*k+3),-2*k*(k+1)*(14*k**3-37*k*k+39*k-15),F(4,5)*k*(k+1)*(42*k**4-175*k**3+301*k*k-245*k+85)]
 return [2*k*(k-1),-2*k*(k-1)*(4*k-7),F(8,3)*k*(k-1)*(20*k*k-74*k+69),-4*k*(k-1)*(112*k**3-672*k*k+1378*k-963),F(16,5)*k*(k-1)*(336*k**4-2800*k**3+8876*k*k-12705*k+6920)]
R1=[('A',1,2),('A',4,2),('A',5,1),('A',7,2),('A',8,1),('D',6,1)]
R2=[('A',2,2),('A',3,1),('A',6,2),('A',9,1),('D',5,3)]
val=lambda r:[sum(m*moments(f,k)[i] for f,k,m in r) for i in range(5)]
h=s.symbols('h');defect=s.factor(moments('A',h-1)[3]/(h-1)-moments('D',h/2+1)[3]/(h/2+1))
print('theta basis coefficients',c,'theta q^0..q^12',theta)
print('rank43',val(R1),val(R2),'delta5',val(R2)[4]-val(R1)[4]);print('Coxeter defect',defect)
print('symplectic24 vs split orthogonal24 maximal isotropic counts',prod(7**i+1 for i in range(1,25)),2*prod(7**i+1 for i in range(1,24)))
print('orthogonal fixed isotropic-line containment denominator',7**23+1)
print('Schottky normalization h30: factor should1; manuscript factor',F(15*896,4))
