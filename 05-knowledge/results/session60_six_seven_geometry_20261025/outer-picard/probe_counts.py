"""Exploratory counts only. If b=0, the first two traces do NOT fix reciprocity sign."""
from itertools import product
import sympy as s

def count2(p):
 n=next(a for a in range(2,p) if pow(a,(p-1)//2,p)==p-1)
 es=list(product(range(p),repeat=2));one=(1,0)
 def add(x,y):return((x[0]+y[0])%p,(x[1]+y[1])%p)
 def neg(x):return((-x[0])%p,(-x[1])%p)
 def sub(x,y):return add(x,neg(y))
 def mul(x,y):return((x[0]*y[0]+n*x[1]*y[1])%p,(x[0]*y[1]+x[1]*y[0])%p)
 def chi(x):
  a=(x[0]*x[0]-n*x[1]*x[1])%p
  return 0 if a==0 else 1 if pow(a,(p-1)//2,p)==1 else -1
 ds={x:sub(one,mul(x,x)) for x in es}
 c=sum(1+chi(mul(mul(ds[x],ds[z]),add(add(one,x),z))) for x in es for z in es)+p*p+1
 return c+17*p*p,n

def count1(p):
 chi={0:0};chi.update({a: 1 if pow(a,(p-1)//2,p)==1 else -1 for a in range(1,p)})
 return sum(1+chi[((1-x*x)*(1-z*z)*(1+x+z))%p] for x in range(p) for z in range(p))+18*p+1

T=s.symbols('T')
for p in [5,7,11,13,17,19]:
 n1=count1(p);n2,nonresidue=count2(p)
 a=n1-1-p*p-18*p; a2=n2-1-p**4-18*p*p
 b=(a*a-a2)//2
 poly=T**4-a*T**3+b*T*T-p*p*a*T+p**4 if b else None
 print(dict(p=p,n1=n1,n2=n2,nonresidue=nonresidue,a=a,a2=a2,b=b,reciprocity_sign='plus' if b else 'UNDETERMINED FROM TWO TRACES',poly=str(s.factor(poly)) if poly is not None else None),flush=True)
