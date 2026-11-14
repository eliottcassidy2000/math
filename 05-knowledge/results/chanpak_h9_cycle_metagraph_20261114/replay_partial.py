from fractions import Fraction
from itertools import combinations, permutations
from math import factorial, gcd
from collections import Counter

def check(x):
    if not x: raise RuntimeError('partial-cycle control')

def parts(n, at=1):
    if n==0:
        yield ()
        return
    for i in range(at,n+1):
        for p in parts(n-i,i):
            yield (i,)+p

def z(p):
    r=1
    for i,j in Counter(p).items(): r*=i**j*factorial(j)
    return r

def formula(n,k,anti):
    m=n-k
    f=k%2
    if anti:
        u=(k-f)//2
        I=(k-3)*(k+1)//4 if f else k*k//4-k//2-1
        ret=0
        for p in parts(m):
            if Counter(p)[1]>1-f or any(i!=1 and i%4!=2 for i in p):continue
            E=sum(i//2 for i in p)+sum(gcd(i,j) for i,j in combinations(p,2))
            e=I+E+f*len(p)+u*sum(gcd(2,i) for i in p)
            ret+=Fraction(2**e,(2-f)*z(p))
        return ret
    ret=0
    for d in range(1,k+1,2):
        if k%d:continue
        phi=sum(gcd(d,i)==1 for i in range(1,d+1))
        for p in parts(m):
            if any(i%2==0 for i in p):continue
            E=sum(i//2 for i in p)+sum(gcd(i,j) for i,j in combinations(p,2))
            e=k*(k-3)//(2*d)+E+(k//d)*sum(gcd(d,i) for i in p)
            ret+=Fraction(phi*2**e,k*z(p))
    return ret

def brute_fix(n,k,perm,anti):
    forced={(i,i+1):0 for i in range(k-1)}
    forced[0,k-1]=1
    remaining=set(combinations(range(n),2)); dof=0
    while remaining:
        p=min(remaining); q=p; accumulated=0; anchor=None
        while True:
            remaining.remove(q)
            if q in forced:
                v=forced[q]^accumulated
                if anchor is not None and anchor!=v:return 0
                anchor=v
            i,j=perm[q[0]],perm[q[1]]
            accumulated^=int(i>j)^anti
            q=(min(i,j),max(i,j))
            if q==p:break
        if accumulated:return 0
        if anchor is None:dof+=1
    return 1<<dof

def brute(n,k,anti):
    m=n-k
    y=0
    for s in range(k):
        A=tuple((s-i if anti else i+s)%k for i in range(k))
        for B in permutations(range(k,n)):
            y+=brute_fix(n,k,A+B,anti)
    return Fraction(y,k*factorial(m))

for n in range(3,11):
    for k in range(3,n+1):
        for anti in (False,True):
            a,b=formula(n,k,anti),brute(n,k,anti)
            check(a==b and a.denominator==1)
print('Both partial-cycle Burnside sums match an independent signed edge-orbit solver through n=10.')
