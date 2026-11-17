from slot_profiles import end_profiles,ordinal
from itertools import product
from generic import F
from math import isqrt,prod
C=(24,1,3,6,14)
bank={121:ordinal(C,C),311:(510,44,384,100,430,132,54,74,234),
      317:(482,372,481,439,389,272,312,354,0),
      383:(218,496,379,498,256,17,48,372,97)}
def profile(a):
    h,_=end_profiles(a);return h[-1],sum(h[s]*h[-1-s] for s in range(len(h)))
def cycle(a):
    n=len(a)
    return tuple(x|1<<n for x in a)+(1<<(n+1),(1<<n)-1)
if profile(C)!=(11,100):raise ValueError("decimal stock")
for h,a in bank.items():
    if profile(a)!=(h,10000) or end_profiles(cycle(a))[0][-1]!=10000+h:
        raise ValueError(h)
def inverse(N):
    j=0;t=1
    while t*10000<N:j+=1;t*=10000
    r=N-t;k=[]
    if j<1 or not 0<r<t:return None
    for q in bank:
        while r%q==0:r//=q;k.append(q)
    return tuple(k) if r==1 and len(k)==j else None
p=100098587
if inverse(p)!=(311,317) or any(p%d==0 for d in range(2,isqrt(p)+1)):
    raise ValueError("alternating stock")
for j in (1,2,3,5,10):
    k=(121,)*(j-1)+(317,)
    if inverse(10000**j+prod(k))!=k:raise ValueError(k)
print("prime selected by factor test",p,inverse(p))
def counts(b):
    n=len(b)
    return 11**n+sum(7921*11**(n-1-j+i)*100**(j-i-1)
                       +(121*11**(n-2) if j==i+1 else 0)
               for i in range(n) for j in range(i+1,n) if b[i] and not b[j])
def graph(b):
    n=len(b);a=()
    for _ in b:a=ordinal(a,C)
    rows=[r|((not b[i//5])<<len(a)) for i,r in enumerate(a)]
    return tuple(rows)+(sum(int(b[i//5])<<i for i in range(len(a))),)
for b in ((1,0),(1,1,0),(1,0,0)):
    if end_profiles(graph(b))[0][-1]!=counts(b):raise ValueError(b)
for n in range(2,12):
    high=max(counts((1,)+x+(0,)) for x in product(range(2),repeat=n-2))
    lo=n//2;hi=n-lo
    exact=100**n-100**lo*11**hi-100**hi*11**lo+3*11**n
    if high!=exact or exact>=100**n:raise ValueError((n,high))
print("decimal word formula and sorted maximum replay")
for k,p in ((1,5),(2,23),(10,8573437)):
    q=F(5**k-1)
    if q!=p or any(q%d==0 for d in range(2,isqrt(q)+1)):raise ValueError(k)
    print("decimal prime offset",k,q,10**k+q)
