import sympy as s
from math import gcd
from collections import Counter
X,z=s.symbols('X z')
C=Counter()
for p in range(2,13):
 for r in range(1,p):
  if gcd(p,r)>1:continue
  F=2**r*X**p-3**p
  for t in range(p):
   A=(-t*r)%p; a=(-(t+1)*r)%p
   eta=lambda m:((t+m)*r+p-1)//p-(t*r+p-1)//p
   for eps in (-1,1):
    h=max(0,-eta(1)-eps)
    P=sum(3**(p-1-m)*2**(h+eta(m)+eps*(m==1))*X**m for m in range(p))
    Q=z**(a+1)-z**a+1 if eps==1 else z**(p-a)-z+1
    N=abs(s.resultant(z**p-2,Q,z))
    e=r*(p-1)-A if eps==1 else r*(p-2)
    assert abs(s.resultant(F,P,X))==2**e*3**(p*(p-1))*N,(p,r,t,eps)
    C['resultants']+=1
print(dict(C))
