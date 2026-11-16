from itertools import product,combinations
from math import gcd
from functools import cache
def rep(k):
    s=[1]
    for i in range(k):
        t=[0]*(len(s)+3*3**i)
        for n,u in enumerate(s):
            for d in range(4):t[n+d*3**i]+=u
        s=t
    return s[:3**k]
f,g=1,1
for k in range(1,10):
    t=rep(k)
    assert max(t)==g and t.count(1)==2**(k+1)-1
    f,g=g,f+g
for u in range(1,51,2):
    for v in range(1,51):
        if gcd(u,v)>1 or (u,v)==(1,1):continue
        a,b=(u,v-u) if v>u else (u-2*v,v) if u>2*v else (2*v-u,u-v)
        assert a%2 and gcd(a,b)==1 and a>0 and b>0 and a+b<u+v
assert [(a,b) for a in range(1,68) for b in range(a,68) if a*a+b*b==68*68]==[(32,60)]
def count(w,ab):
    n=len(w);m=n+2
    A=[((1<<n)-1)^((1<<(i+1))-1) for i in range(n)]+[0,0]
    A[n+(not ab)]|=1<<(n+int(ab))
    for i,z in enumerate(w):
        for d,j in zip(divmod(int(z),2),(n,n+1)):
            if d:A[j]|=1<<i
            else:A[i]|=1<<j
    @cache
    def f(S,i):
        U=S^(1<<i)
        return sum(f(U,j) for j in range(m) if U&A[i]>>0 & (1<<j)) if U else 1
    out=sum(f((1<<m)-1,i) for i in range(m));f.cache_clear()
    return out
for w,ab,h in [('1000',1,17),('10212',0,49),('100110',1,63),('1012000',0,343)]:
    assert count(w,ab)==h
print('Redundant Fibonacci law, primitive inverse chart, and literal 17/49/63/343 seeds passed.')

A=set(range(23,68,2))
for k in range(1,6):
    A={3*h+2*d for h in A for d in (0,1,2)}
    assert A==set(range(23*3**k,67*3**k+2*(3**k-1)+1,2))
print('Conditional absolute seed-band interval identity checked at five levels.')
