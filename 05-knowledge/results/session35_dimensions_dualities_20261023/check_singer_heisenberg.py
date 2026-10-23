import sympy as s
from math import comb
from collections import Counter
import json

def require(v,msg):
 if not v: raise RuntimeError(msg)
x=s.symbols('x');N=57
f=sum(x**a-x**(N-a) for a in [1,7,49])
g=s.gcd(f,x**N-1);require(g==x-1,'only constant kernel')
# Exact Euclidean certificate over Q[x], more informative than a numeric rank.
u,v,d=s.gcdex(f,x**N-1,x)
require(s.expand(u*f+v*(x**N-1))==x-1,'Bezout certificate')
allowed={};triples=[]
for a in range(9):
 for b in range(9-a):
  c=8-a-b
  if (b+2*c)%3:continue
  fixed=2*(comb(a,2)+b*c)
  trace=(3*fixed-56)//2
  require((3*fixed-56)%2==0,'integrality')
  allowed.setdefault(fixed,set()).add(trace)
  triples.append([a,b,c,fixed,trace])
require(set(allowed)=={2,14,20,32,56},'allowed list')
fixed_aug=N//3-1
require(fixed_aug==18 and fixed_aug not in allowed,'restriction obstruction')
D=[8,15,19,38,48,50,51,56]
require(Counter((a-b)%N for a in D for b in D if a!=b)==Counter(range(1,N)),'Singer differences')
require(sorted(7*a%N for a in D)==D,'Frobenius line')
# Coefficient invariance under multiplication7 (Frobenius), skew reflection.
coeff={i:int(s.expand(f).coeff(x,i)) for i in range(N)}
require(all(coeff[i]==coeff[(7*i)%N] for i in range(N)),'Frobenius form')
require(all(coeff[i]==-coeff[(-i)%N] for i in range(N)),'alternating')
result={'N':N,'augmentation_dim':56,'heisenberg_dim':57,'form_exponents':[1,7,49],'form_gcd':str(g),'bezout_u':str(u),'bezout_v':str(v),'order3_augmentation_fixed_dim':fixed_aug,'allowed_E7_fixed_dims':sorted(allowed),'allowed_E7_traces':sorted({z for vs in allowed.values() for z in vs}),'sl8_multiplicity_triples':triples,'Singer_D':D,'status':'FINITE-EXACT controls; universal statements proved in companion note'}
print(json.dumps(result,indent=2))
