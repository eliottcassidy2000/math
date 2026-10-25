from itertools import product,combinations
from math import gcd
from functools import reduce
from collections import Counter
V=list(product((-1,0,1),repeat=2)); rows=[(1,i,j) for i,j in V]
def det(a,b,c):return a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0])
agg=Counter()
for mask in range(512):
 R=[r for i,r in enumerate(rows) if mask>>i&1]
 if not R:inv=()
 else:
  g2=reduce(gcd,(abs(a[i]*b[j]-a[j]*b[i]) for a,b in combinations(R,2) for i,j in combinations(range(3),2)),0)
  g3=reduce(gcd,(abs(det(*abc)) for abc in combinations(R,3)),0)
  inv=(1,) if not g2 else (1,g2) if not g3 else (1,g2,g3//g2)
 agg[inv]+=(-1)**len(R)
expected={():1,(1,):-9,(1,1):20,(1,2):8,(1,1,1):1,(1,1,2):-10,(1,1,3):-4,(1,1,4):-4,(1,2,2):-3}
if {k:v for k,v in agg.items() if v}!=expected:raise ValueError(agg)
print('signed invariant factors',dict(agg),'PASS')
for n in range(1,26):
 brute=sum(all((e+i*x+j*z)%n for i,j in V) for e,x,z in product(range(n),repeat=3))
 calc=sum(a*n**(3-len(inv))*__import__('math').prod(gcd(n,d) for d in inv) for inv,a in agg.items())
 if brute!=calc:raise ValueError((n,brute,calc))
 print('cyclic order',n,'complement',brute,'PASS')
for p in (2,3,5,7,11,13):
 brute=sum(all((e+i*x+j*z)%p for i,j in V) for e,x,z in product(range(p),repeat=3))
 calc=p**3-4*p*p+6*p-3 if p==2 else p**3-9*p*p+24*p-16 if p==3 else p**3-9*p*p+28*p-20
 if brute!=calc:raise ValueError(('additive',p,brute,calc))
 print('prime additive',p,brute,'PASS')
