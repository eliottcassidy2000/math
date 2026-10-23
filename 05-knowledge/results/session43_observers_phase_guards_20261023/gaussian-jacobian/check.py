import json, math
from fractions import Fraction

def need(x,msg):
 if not x: raise RuntimeError(msg)
def proj(poly,L):return {q:c for q,c in poly.items() if q%L==0 and c}
def mul(a,b):
 c={}
 for x,u in a.items():
  for y,v in b.items():c[x+y]=c.get(x+y,0)+u*v
 return {q:v for q,v in c.items() if v}
checks=0
for a in range(1,13):
 for b in range(1,13):
  f={q:Fraction(q*q+1) for q in range(-150,151)}
  need(proj(proj(f,a),b)==proj(f,math.lcm(a,b)),'lcm projection');checks+=1
for L in (4,8,9,72):
 for j in range(1,L):need(proj({j:1},L)=={},'premature alias');checks+=1
 need(proj({L:1},L)=={L:1},'missing alias');checks+=1
need(proj(mul({1:1},{-1:1}),8)=={0:1},'shared phase')
need(mul(proj({1:1},8),proj({-1:1},8))=={},'independent phase');checks+=2
# Newton identities exact for arbitrary explicit rational orbit lists.
for roots in ([1,2,-3],[-2,-1,0,1,2],[Fraction(1,2),Fraction(3,4),-1]):
 n=len(roots);e=[Fraction(1)]
 for k in range(1,n+1):e.append(sum((-1)**(j-1)*e[k-j]*sum(r**j for r in roots) for j in range(1,k+1))/k)
 for r in roots:need(sum((-1)**j*e[j]*r**(n-j) for j in range(n+1))==0,'Newton');checks+=1
# The formal cyclic list has power sums zero before N, p_N=N*z^N;
# Newton gives e_N=(-1)^(N-1)*z^N. Its coefficient certifies sharpness.
for n in range(1,73):
 e=[Fraction(1)]
 for k in range(1,n+1):e.append(sum((-1)**(j-1)*e[k-j]*(n if j==n else 0) for j in range(1,k+1))/k)
 need(e[n]==(-1)**(n-1),'Newton sharpness');checks+=1
# Symbolic independent determinant tests use SymPy, not floating sampling.
import sympy as s
x,y=s.symbols('x y')
for w in range(1,6):
 t=x**w*y; A=1+t+t*t;B=2-t
 F=x*A;G=y*B
 det=s.diff(F,x)*s.diff(G,y)-s.diff(F,y)*s.diff(G,x)
 u=s.symbols('u');aa=1+u+u*u;bb=2-u
 rhs=(aa*s.diff(u*bb,u)+w*u*s.diff(aa,u)*bb).subs(u,t)
 need(s.expand(det-rhs)==0,'quotient determinant');checks+=1
for L in (8,9,72):
 # integer charge of y^(L-1) is -(L-1), equal to 1 mod L but not 1.
 need((-(L-1)-1)%L==0 and -(L-1)!=1,'finite equivariance');checks+=1
print(json.dumps({'status':'PASS','checks':checks,'alias_72_radial_moment':str(math.factorial(36)),'scope':'finite regression checks; general proofs in THEOREMS.md'},sort_keys=True,indent=2))
