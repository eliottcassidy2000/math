"""Independent exact regressions, no asserts and no source mutations."""
import hashlib,json,math,pathlib,subprocess,sys
from fractions import Fraction
import sympy as s
HERE=pathlib.Path(__file__).resolve().parent
SOURCE=next((p for p in (HERE.parent/'session43-gaussian-jacobian', HERE.parent/'gaussian-jacobian') if (p/'MANIFEST.json').is_file()), None)
if SOURCE is None: raise RuntimeError('Gaussian source packet missing in workspace and publication layouts')
count=0
def require(ok,label):
 global count
 if not ok: raise RuntimeError(label)
 count+=1
manifest=json.loads((SOURCE/'MANIFEST.json').read_text())
for name,digest in manifest['files'].items():
 require(hashlib.sha256((SOURCE/name).read_bytes()).hexdigest()==digest,'source hash '+name)
require((SOURCE/'receipt.json').read_bytes()==(SOURCE/'receipt.optimized.json').read_bytes(),'producer modes')
# Actual cyclotomic ring summation, independent of charge-divisibility implementation.
t=s.symbols('t')
for N in range(1,17):
 cycl=s.Poly(s.cyclotomic_poly(N,t),t)
 for j in range(1,N+1):
  sm=s.Poly(sum(t**(g*j) for g in range(N)),t).rem(cycl)
  require(sm==s.Poly(N if j==N else 0,t),'cyclic sharpness')
# Formal Newton identities checked against a separately expanded product.
u=s.symbols('u')
for roots in ([1,2,-3],[0,0,2,-2],[Fraction(1,2),Fraction(-2,3),Fraction(7,5)]):
 E=[s.Integer(1)]
 for k in range(1,len(roots)+1):
  E.append(s.expand(sum((-1)**(j-1)*E[k-j]*sum(s.Rational(r)**j for r in roots) for j in range(1,k+1))/k))
 require(s.expand(s.prod(u-r for r in roots)-sum((-1)**k*e*u**(len(roots)-k) for k,e in enumerate(E)))==0,'Newton product')
# Composition counts each residue modulo lcm uniformly, not only a congruence test.
for a in range(1,13):
 for b in range(1,13):
  L=math.lcm(a,b); hist=[0]*L
  for i in range(a):
   for j in range(b):hist[(i*(L//a)+j*(L//b))%L]+=1
  require(hist==[math.gcd(a,b)]*L,'angle convolution')
# Gaussian radial normalization: |Z|^2=(X^2+Y^2)/2, standard real Gaussians.
def gm(n):return 0 if n%2 else math.prod(range(1,n,2))
for k in range(37):
 radial=sum(math.comb(k,j)*gm(2*j)*gm(2*(k-j)) for j in range(k+1))
 require(Fraction(radial,2**k)==math.factorial(k),'radial normalization')
# Shared versus independent phases: weights on Z1,W1,Z2,W2.
weights=[(1,0),(0,-1),(1,-1),(2,-2)]
for a,b in weights:
 require(((a+b)%72==0)==(a==-b),'shared low-charge selection')
require((1-1)%72==0 and not (1%72==0 and -1%72==0),'shared relative phase survives')
# Direct Jacobian and leading coefficient calculations.
x,y=s.symbols('x y')
for w in range(1,6):
 for da in range(4):
  for db in range(4):
   A=2+u**da; B=3+2*u**db; v=x**w*y
   F=x*A.subs(u,v);G=y*B.subs(u,v)
   J=s.expand(s.diff(F,x)*s.diff(G,y)-s.diff(F,y)*s.diff(G,x))
   D=s.expand(A*B+u*A*s.diff(B,u)+w*u*s.diff(A,u)*B)
   require(s.expand(J-D.subs(u,v))==0,'Jacobian identity')
   require(s.degree(D,u)==da+db,'no highest-degree cancellation')
for L in (2,3,8,9,72):
 F=x+y**(L-1);G=y
 require(s.det(s.Matrix([F,G]).jacobian([x,y]))==1,'Keller shear')
 require((-(L-1)-1)%L==0,'finite weight')
 require((s.Poly(F,x,y).total_degree()>1)==(L>=3),'nonlinear qualification')
# Long witness direct multinomial Wick contraction, not coefficient-identity shortcut.
z,b,T=s.symbols('z b T')
P=(1+z)*(b-s.Rational(1,2)*(2+z)*T**2)
def wick(poly):
 return sum(c*math.factorial(a)*gm(k) for (a,d,k),c in s.Poly(poly,z,b,T).terms() if a==d)
for m in range(1,9):
 pp=s.expand(P**m)
 require(wick(pp)==0,'Long pure moment')
 require(wick(z*pp)==math.factorial(m),'Long mixed moment')
print(json.dumps({'status':'PASS','checks':count,'producer_manifest_sha256':hashlib.sha256((SOURCE/'MANIFEST.json').read_bytes()).hexdigest(),'alias_72_radial_mean':str(math.factorial(36)),'scope':'Independent finite checks; general arguments separately reviewed; no Lean rebuild'},sort_keys=True,indent=2))
