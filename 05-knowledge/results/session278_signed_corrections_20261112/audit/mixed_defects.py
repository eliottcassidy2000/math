"""Independent mixed-weight formal-series check for cross equation (8)."""
import sympy as s
from functools import lru_cache
t,z,v=s.symbols('t z v');a=s.Rational(3,2)
R=z+2*v
@lru_cache(None)
def P(m,j,k):
 h=z**j*v**k
 if m==0:return h
 argument=s.expand((a+2*j+3*k+m-1)*h-R*h)
 ans=0
 for (u,w),c in s.Poly(argument,z,v).terms():ans+=c*P(m-1,u,w)
 if m>=2:
  for (u,w),c in s.Poly(R*h,z,v).terms():ans+=(m-1)*c*P(m-2,u,w)
 return s.expand(ans)
# Construct each coefficient from finite exponential products, avoiding symbolic infinite exp.
L=6
q=s.series(z*(1-(1-t)**-1)+v*(1-(1-t)**-2),t,0,L+1).removeO()
expq=sum(s.series(q**j,t,0,L+1).removeO()/s.factorial(j) for j in range(L+1))
series=s.series(s.series((1-t)**-a,t,0,L+1).removeO()*expq,t,0,L+1).removeO().expand()
for m in range(L+1):
 d=s.cancel(P(m,0,0)-s.factorial(m)*series.coeff(t,m))
 if d!=0:print(m,d)
 assert d==0
print('PASS: mixed quartic/sextic formal generating formula m=0..6, a=3/2')
