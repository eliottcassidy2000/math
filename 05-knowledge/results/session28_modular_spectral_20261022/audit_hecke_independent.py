"""Independent q-product / E4,E6 basis verification. Run with Python3 + sympy."""
from fractions import Fraction as F
import json, math
from pathlib import Path
import sympy as s
N=36

def mul(a,b):
 c=[0]*(N+1)
 for i,x in enumerate(a):
  for j,y in enumerate(b[:N+1-i]):c[i+j]+=x*y
 return c

def powr(a,k):
 r=[1]+[0]*N
 for _ in range(k):r=mul(r,a)
 return r
E4=[1]+[240*sum(d**3 for d in range(1,n+1) if n%d==0) for n in range(1,N+1)]
E6=[1]+[-504*sum(d**5 for d in range(1,n+1) if n%d==0) for n in range(1,N+1)]
D=[0,1]+[0]*(N-1)
for m in range(1,N+1):
 factor=[0]*(N+1)
 for k in range(min(24,N//m)+1):factor[k*m]=(-1)**k*math.comb(24,k)
 D=mul(D,factor)
assert all(a-b==1728*d for a,b,d in zip(powr(E4,3),powr(E6,2),D))
B=[mul(powr(E4,11-3*j),powr(E6,2*j)) for j in range(4)]
A=s.Matrix([[B[j][i] for j in range(4)] for i in range(4)])
c=A.inv()*s.Matrix([1,0,0,0])
f=[sum(c[j]*B[j][i] for j in range(4)) for i in range(N+1)]
assert f[:5]==[1,0,0,0,168498000]
T=s.Matrix.hstack(*[A.inv()*s.Matrix([b[2*i]+(2**43*b[i//2] if i%2==0 else 0) for i in range(4)]) for b in B])
for j in range(4):
 for n in range(N//2+1):assert sum(T[i,j]*B[i][n] for i in range(4))==B[j][2*n]+(2**43*B[j][n//2] if n%2==0 else 0)
x=s.Symbol('x')
p=s.Poly(T.charpoly(x).as_expr(),x)
cp=p.exquo(s.Poly(x-(1+2**43),x))
assert cp.as_expr()==x**3+2209944*x**2-15663522502656*x-19976984434430705664
assert cp.is_irreducible and cp.count_roots(-s.oo,s.oo)==3
cyc=s.Matrix.hstack(c,T*c,T*T*c,T*T*T*c)
assert cyc.det()!=0
out={'basis':'E4^(11-3j) E6^(2j), j=0..3','method':'Delta independently from q product; Hecke verified through q18','theta_coefficients':list(map(str,f[:13])),'theta_basis_coordinates':list(map(str,c)),'cusp_characteristic_polynomial':str(cp.as_expr()),'cusp_irreducible':True,'cusp_real_roots':3,'cyclic_rank':cyc.rank(),'cyclic_determinant_in_independent_basis':str(cyc.det()),'T2_theta_q2':str(f[4]),'status':'Exact finite arithmetic audit; modular space dimension and Hecke theorem are cited separately; no rank88 minimum audit.'}
print(json.dumps(out,indent=2))
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
