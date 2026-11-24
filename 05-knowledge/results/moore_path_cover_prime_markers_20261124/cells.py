from math import gcd
from profile import mask, evaluate, profile, rows

def strong(T):
    for U in T,rows(len(T),lambda i,j:not(T[i]>>j&1)):
        r=1
        while True:
            a=r
            for i in range(len(T)):
                if r>>i&1:a|=U[i]
            if a==r:break
            r=a
        if r!=(1<<len(T))-1:return False
    return True

def scan(n):
    m=(n-1)*(n-2)//2
    P=[profile(T) if strong(T) else None for T in (mask(n,z) for z in range(1<<m))]
    return [(evaluate(a,1),i,j,a,P[j]) for i,a in enumerate(P) if a
            for b in range(m) if i<(j:=i^(1<<b)) and P[j] is not None
            and a!=P[j] and evaluate(a,1)==evaluate(P[j],1)]

def Z(n,k):
    return k**(n+2)+k*((k+1)**n-(k-1)**n)

def valuation(x,p):
    a=0
    while x%p==0:a+=1;x//=p
    return a

if __name__=='__main__':
    for n in range(3,7):
        s=scan(n)
        assert len(s)==(122 if n==6 else 0)
    P=profile(mask(6,40));Q=profile(mask(6,42))
    assert P==(0,0,24,0,12,0,1) and Q==(0,0,22,0,14,0,1)
    print('common-spine cells',mask(6,40),mask(6,42),'H 37, Q',evaluate(P,2),evaluate(Q,2),'R',evaluate(P,3),evaluate(Q,3))
    for k in range(1,500):
        a,b=evaluate(P,k)//k**2,evaluate(Q,k)//k**2
        assert gcd(a,b)==(1+(k%2==0))*(37 if k%37 in (1,36) else 1)
    A=(16,1,3,7,14); B=(20,1,2,7,14)
    assert profile(A)==profile(B)==(0,2,0,6,0,1)
    for n in range(1,9):
        T=rows(n+2,lambda i,j:(i!=n+1) if j<n+1 else i==n)
        for k in range(1,5):assert evaluate(profile(T),k)==Z(n,k)
    for d in (2,3,7,57):
        m=d*d+1;w=d+1
        print(d,'two marker exponents',m-w,w,'pair read',2*d-2,2*d)
    assert Z(9,1)==57*Z(3,1)
