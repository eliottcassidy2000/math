from functools import cache
from math import isqrt

def need(x):
    if not x:raise AssertionError

def table(rows):
    n=len(rows);full=(1<<n)-1
    @cache
    def end(S,j):
        Z=S^(1<<j)
        return 1 if not Z else sum(end(Z,i) for i,r in enumerate(rows) if Z>>i&1 and r>>j&1)
    @cache
    def H(S):
        return sum(end(S,j) for j in range(n) if S>>j&1) if S else 1
    return full,end,H

def menu(rows,v):
    n=len(rows)
    F,_,H=table(rows)
    rev=tuple(sum(1<<i for i,r in enumerate(rows) if r>>j&1) for j in range(n))
    _,start,_=table(rev)
    E=sum(start(S,v)*H(F^S) for S in range(1,F+1) if S>>v&1 and S.bit_count()%2==0)
    O=sum(start(S,v)*H(F^S) for S in range(1,F+1) if S>>v&1 and S.bit_count()%2==1)
    return H(F),E,O,H(F^(1<<v))

def strong(T):
    F=(1<<len(T))-1
    for s in range(len(T)):
        R=1<<s
        while True:
            U=R
            for i,r in enumerate(T):
                if R>>i&1: U|=r
            if U==R: break
            R=U
        if R!=F:return False
    return True

def adj(rows,v,N):
    n=len(rows)
    return tuple(r+((not N or i!=v)<<n) for i,r in enumerate(rows))+((1<<v) if N else 0,)

def add_root(D,w):
    n=len(D)
    return tuple(r+(int(not (w>>i&1))<<n) for i,r in enumerate(D))+(w,)

def hp(rows):
    F,_,H=table(rows);return H(F)

def extend(t,N):
    h,e,o,g=t
    return (h+2*e,e+o-g,e+o,g) if N else (h,e+o,e+o,g)

def prime(x):
    return x>1 and all(x%d for d in range(2,isqrt(x)+1))

def F(x):
    p=1
    while x:
        z=x//2
        if x&1==0:p+=z+(z&1)
        x=z
    return p

def program(B,M,res):
    a=(M&-M).bit_length()-1
    m=M>>a
    u,g=0,1
    def T(k):
        nonlocal u,g
        g+=(u+1)<<k
        u=((u+1)<<(k+1))-2
    if a:
        need(res&1 and B%(1<<a)==0)
        f=(res+482)*pow(25,-1,1<<a)%(1<<a)
        for k in range(1,a):
            if ((f^g)>>k)&1:T(k)
    if m>1:
        v,f=next((v,f) for v in range(m) for f in range(m) if (-482+B*v+25*f-res)%m==0)
        def exponent(e):
            j=max(1,a)
            while pow(2,j+1,m)!=e%m:j+=1
            return j
        k=exponent(1);l=exponent(2)
        for _ in range(u%m):T(k)
        steps=(1-v)%m
        drift=(pow(2,k,m)*sum(1-i for i in range(steps)))%m
        for _ in range((f-g-drift)%m):T(l)
        for _ in range(steps):T(k)
    need((-482+B*(u+1)+25*F(u)-res)%M==0 and F(u)==g)
    return u

for w in (549,549^(1<<15)):
    n=7;rows=[0]*n;bit=0
    for i in range(n):
        for j in range(i+1,n):
            if w>>bit&1:rows[i]|=1<<j
            else:rows[j]|=1<<i
            bit+=1
    print('one arc',w,rows,hp(rows))
    need(hp(rows)==(77 if w==549 else 91))
D=(262,268,280,49,35,7,319,383,56)
roots=[add_root(D,w) for w in (183,159)]
vals=[menu(A,9) for A in roots]
need(D[7]==((1<<9)-1)^(1<<7) and vals==[(1947,2698,2748,175),(1961,2698,2748,175)])
for word in ((1,),(1,1),(1,0,1),(1,0,0,1)):
    X=roots[:];V=vals[:]
    for N in word:
        X=[adj(x,9,N) for x in X];V=[extend(v,N) for v in V]
    need([menu(x,9) for x in X]==V and all(map(strong,X)))
    need(V[1][0]-V[0][0]==14 and V[0][1:]==V[1][1:])
    print('common word',word,V)
need(V[0][1]%7==0 and prime(1049) and prime(1051))
base=extend(vals[0],1)
need(base==(7343,5271,5446,175))
for x in range(1,128):
    V=base
    for j in range(x.bit_length()):V=extend(V,x>>j&1)
    need(V[0]//7==-482+1506*(x+1)+25*F(x))
for M in (3,5,7,8,9,14,30):
    a=(M&-M).bit_length()-1
    for L in (max(a,1),max(a,1)+3):
        B=1506<<L
        for res in range(M):
            if not res&1 and M%2==0:continue
            r=program(B,M,res)
            k=max(1,r.bit_length())
            for q in (0,2*M, 20*M):
                u=(q<<k)+r
                need((-482+B*(u+1)+25*F(u)-res-25*(F(q)-1))%M==0)
print('phase and low-bit splices verified')

# Clone vertex 0 of a single-arc template; two directed runs at most.
A=(32,13,25,113,99,198,135,31)
J=list(A);J[2]^=1<<6;J[6]^=1<<2
def clone(A,m):
    n=len(A)
    tags=[0]*m+list(range(1,n))
    R=[0]*len(tags)
    for i,x in enumerate(tags):
        for j,y in enumerate(tags):
            if (x==y and i<j) or (x!=y and A[x]>>y&1):R[i]|=1<<j
    return R
for m in (1,2,3):
    need((hp(clone(A,m)),hp(clone(J,m)))==(246*(1<<m)-95,246*(1<<m)-81))
print('single arc clone rails checked')
