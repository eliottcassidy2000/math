from math import factorial,comb
from functools import cache
from itertools import product
from random import Random

def require(ok,*data):
    if not ok: raise AssertionError(data)

def exts(a,k):
    # Return count of exterior skeletons with k visits to an independent slot
    # (an ordinary consecutive slot-slot transition would merge those runs).
    # a is a tournament; its last vertex is the slot.
    n=len(a)-1
    @cache
    def f(m,v,j):
        if m==(1<<n)-1 and j==k:return 1
        if v==-1: out=(1<<n)-1
        else: out=a[v]
        return sum(f(m|(1<<w),w,j) for w in range(n) if out>>w&1 and not m>>w&1) + (
          f(m,n,j+1) if j<k and (v==-1 or out>>n&1) else 0)
    return f(0,-1,0)

def reader(r):
    n=2*r-2
    return tuple(((1<<n)-(1<<(i+1)))|((i%2)<<n) for i in range(n)) + (
        sum(1<<i for i in range(0,n,2)),)

def insert(a,b):
    n=len(a)-1;m=len(b);ext=(1<<n)-1
    return tuple((x&ext)|(((1<<m)-1)<<n if x>>n&1 else 0) for x in a[:-1])+tuple(
                  (x<<n)|a[-1] for x in b)

def random_tournament(n,R):
    a=[0]*n
    for i in range(n):
        for j in range(i+1,n):
            k,l=(i,j) if R.randrange(2) else (j,i)
            a[k]|=1<<l
    return tuple(a)

@cache
def path_table(a):
    n=len(a)
    @cache
    def dp(mask,v):
        mask-=1<<v
        return sum(dp(mask,w) for w,x in enumerate(a) if mask>>w&1 and x>>v&1) if mask else 1
    return tuple(sum(dp(mask,v) for v in range(n) if mask>>v&1) if mask else 1 for mask in range(1<<n))

def covers(a):
    h=path_table(a)
    @cache
    def f(mask,r):
        if r==0:return int(mask==0)
        s=mask;x=0
        while s:
            x+=h[s]*f(mask^s,r-1)
            s=(s-1)&mask
        return x
    return tuple(f(len(h)-1,r) for r in range(1,len(a)+1))

def slots(C,q):
    return sum(comb(q,k)*c for k,c in enumerate(C,1) if k<=q)

def pext(q,m,M=None):
    if M:return (pow(q,m+2,M)+q*(pow(q+1,m,M)-pow(q-1,m,M)))%M
    return q**(m+2)+q*((q+1)**m-(q-1)**m)

def phi(n):
    from math import gcd
    return sum(gcd(n,j)==1 for j in range(1,n+1))

def depth(M):
    r=1
    while (factorial(r)*factorial(r-1))%M:r+=1
    return r-1

def main():
    for r in range(1,9):
        a=reader(r)
        coeff=[exts(a,i) for i in range(1,r+2)]
        require(coeff[-2:]==[factorial(r-1),0],r,coeff)
        require(all(c%factorial(i)==0 for i,c in enumerate(coeff)),r,coeff)
        print('reader',r,'coefficient rows',coeff[:-1])
    R=Random(726)
    for _ in range(80):
        a=random_tournament(R.randrange(1,9),R)
        b=random_tournament(R.randrange(1,7),R)
        C=covers(b)
        require(all(c%factorial(i)==0 for i,c in enumerate(C,1)),b,C)
        co=tuple(exts(a,i) for i in range(1,len(b)+1))
        require(all(c%factorial(i)==0 for i,c in enumerate(co)),a,co)
        require(sum(x*y for x,y in zip(C,co))==path_table(insert(a,b))[-1],a,b,C,co)
    A=(136,89,75,128,45,15,57,118)
    B=(42,32,163,150,39,200,31,83)
    a,b=covers(A),covers(B)
    require([slots(a,j) for j in (1,2,3)]==[345,5640,44505],a)
    require([slots(b,j) for j in (1,2,3)]==[345,5640,45225],b)
    Q=reader(3)
    require(path_table(insert(Q,B))[-1]-path_table(insert(Q,A))[-1]==1440)
    for M in range(2,101):
        e=M.bit_length()
        m=2*phi(M)*e
        t=phi(M)*e
        for q in range(1,depth(M)+1):
            require(pow(pext(q,m,M),t,M)==pow(q,t*(m+2),M),M,q)
    for m in range(1,7):
        n=m+2;a=list((1<<n)-(1<<(i+1)) for i in range(n))
        a[0]^=1<<(n-1);a[-1]=1
        C=covers(tuple(a))
        for q in range(1,7):require(slots(C,q)==pext(q,m),m,q)
    print('80 arbitrary literal substitutions; 99 modulus construction checks; 36 extreme-edge slot controls.')

if __name__=='__main__':main()
