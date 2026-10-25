from itertools import product,combinations
from collections import Counter
from math import gcd,prod
from functools import reduce
rows=[(2,0),(0,2),(2,2),(2,-2),(2,1),(2,-1),(1,2),(1,-2)]
agg=Counter()
for mask in range(256):
 R=[r for i,r in enumerate(rows) if mask>>i&1]
 g1=reduce(gcd,(abs(x) for r in R for x in r),0)
 g2=reduce(gcd,(abs(a[0]*b[1]-a[1]*b[0]) for a,b in combinations(R,2)),0)
 inv=() if not R else (g1,) if not g2 else(g1,g2//g1)
 agg[inv]+=(-1)**len(R)
print('signed Smith types',dict(agg))
for mods in [(n,) for n in range(1,31)]+[(2,2),(3,3),(2,4),(4,4),(2,2,2),(3,6),(5,5)]:
 elems=list(product(*(range(m) for m in mods)));N=len(elems)
 h=lambda d:prod(gcd(m,d) for m in mods)
 actual=sum(len({tuple((i*a+j*b)%m for a,b,m in zip(x,z,mods)) for i,j in product((-1,0,1),repeat=2)})==9 for x,z in product(elems,repeat=2))
 form=N*N-4*(1+h(2))*N-2-4*h(2)+4*h(4)+2*h(5)+4*h(6)+2*h(2)**2+h(2)*h(4)
 snf=sum(c*N**(2-len(inv))*prod(h(d) for d in inv) for inv,c in agg.items())
 if actual!=form or snf!=form:raise ValueError((mods,actual,form,snf))
 print(mods,N,actual,'distinct',N*actual,'distinct_nonzero',(N-9)*actual,'PASS')
