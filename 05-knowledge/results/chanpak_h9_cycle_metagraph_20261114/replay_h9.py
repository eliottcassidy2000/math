from functools import cache
from itertools import permutations, product
from math import isqrt

def check(x):
    if not x: raise AssertionError('marked cut state')

def ordinal(A,B):
    return tuple(x | (((1<<len(B))-1)<<len(A)) for x in A)+tuple(
        x<<len(A) for x in B)

def apex(B,U):
    return tuple(x|((1<<len(B)) if not (U>>i)&1 else 0)
        for i,x in enumerate(B))+(U,)

def profile(B,U):
    n=len(B); full=(1<<n)-1
    reverse=[full^x^(1<<j) for j,x in enumerate(B)]
    @cache
    def p(s,j,end):
        z=s^(1<<j)
        if not z: return 1
        return sum(p(z,k,end) for k in range(n)
                   if (z&(reverse if end else B)[j])>>k&1)
    @cache
    def h(s):
        return sum(p(s,j,False) for j in range(n) if s>>j&1) if s else 1
    S=D=T=F=0
    for s in range(1,full+1):
        sign=1 if s.bit_count()%2 else -1
        q=h(full^s)
        for i in range(n):
            if s>>i&1:
                if U>>i&1:
                    x=q*p(s,i,False); S+=x; D+=sign*x
                else:
                    y=q*p(s,i,True); T+=y; F+=sign*y
    Q=sum(h(s)*h(full^s) for s in range(full+1))
    return h(full),Q,S,D,T,F

def total(B,U):
    A=apex(B,U)
    # Independent unmarked DP.
    @cache
    def path(S,v):
        T=S^(1<<v)
        return sum(path(T,w) for w in range(len(A)) if (T&A[v])>>w&1) if T else 1
    return sum(path((1<<len(A))-1,i) for i in range(len(A)))

def random_t(n,rng):
    B=[0]*n
    for i in range(n):
        for j in range(i+1,n):
            if rng.randrange(2): B[i]|=1<<j
            else: B[j]|=1<<i
    return tuple(B)


R=[(14,28,24,16,1),(10,28,25,16,1)]
tri=(2,4,1)

def orbit(T):
    n=len(T)
    for p in permutations(range(n)):
        yield tuple(sum(1<<j for j in range(n) if T[p[i]]>>p[j]&1) for i in range(n))

def candidates(n):
    edges=[(i,j) for i in range(n) for j in range(i+1,n)]
    for mask in range(1<<len(edges)):
        T=[0]*n
        for u,(i,j) in enumerate(edges):
            if mask>>u&1:T[i]|=1<<j
            else:T[j]|=1<<i
        yield tuple(T)

def strong(T):
    cov=1
    for _ in T:
        for i,b in enumerate(T):
            if cov>>i&1:cov|=b
    if cov!=(1<<len(T))-1:return False
    cov=1
    for _ in T:
        for i,b in enumerate(T):
            if b&cov:cov|=1<<i
    return cov==(1<<len(T))-1

def H(T):
    n=len(T); z=(1<<n)-1
    @cache
    def start(s,i):
        u=s^(1<<i)
        if not u:return 1
        q=u&T[i]; c=0
        while q:
            b=q&-q;q-=b;c+=start(u,b.bit_length()-1)
        return c
    ans=sum(start(z,i) for i in range(n))
    start.cache_clear()
    return ans

for n in (3,4,5):
    recorded={T for T in candidates(n) if strong(T) and H(T)==9}
    expected=set(orbit(R[0]))|set(orbit(R[1])) if n==5 else set()
    check(recorded==expected)
check(not set(orbit(R[0]))&set(orbit(R[1])))
check(all(profile(T,0)[:2]==(9,84) for T in R))

def g(bits):
    P=S=0; K=1
    for b in bits:
        if not b:K+=P+S
        P,S=2*P+b,b
    return K

def offset(W,i,j):
    vals=[]
    for z in (i,j):
        x=0
        for b in W[:z]:x=2*x+b
        r=2*W[z]+W[z+1]
        q=int(z+2<len(W) and not W[z+2])
        vals.append((r,(0,2+3*x-x%2-2*q,-4-3*x+x%2+2*q,0)[r]))
    (r,c),(s,d)=vals
    return 3*(c+d)+((0,-2,2,0)[r]*(0,-1,1,0)[s] if j==i+2 else 0)

def expand(W,i,j):
    pos=0; B=(); U=0
    while pos<len(W):
        if pos in (i,j):
            U|=((1<<(2*W[pos]+W[pos+1]))-1)<<len(B)
            B=ordinal(B,tri)
            pos+=2
        else:
            U|=W[pos]<<len(B);B=ordinal(B,(0,));pos+=1
    return B,U

for w in range(4,12):
    for W in product((0,1),repeat=w):
        vals=set()
        for i in range(w-3):
            for j in range(i+2,w-1):
                k=offset(W,i,j)
                vals.add(k)
                if w<=6:
                    B,U=expand(W,i,j)
                    check(H(apex(B,U))==9*g(W)+k)
        check(all(k+2 not in vals or k+4 not in vals for k in vals))
        check(all(k+2 not in vals or k%6==0 or (k+2)%6==0 for k in vals))
print('Source census and exact two-triangle interaction/no-triple controls passed.')

choices=[('0',27,'01'),('01',18,'1'),('',9,'101'),('',20,'101')]
def glue(left,R,M,right):
    B=()
    U=0
    for bit in left:
        U|=int(bit)<<len(B);B=ordinal(B,(0,))
    U|=M<<len(B);B=ordinal(B,R)
    for bit in right:
        U|=int(bit)<<len(B);B=ordinal(B,(0,))
    return B,U

for i,(l,m,r) in enumerate(choices):
    B,U=glue(l,R[0],m,r)
    check(profile(B,U)[:4]==(9,672,245,9))
    check(H(apex(B,U))==89+2*i)
    for extra,N in [('0',343+2*i),('000',245*8-147+2*i),
                    ('010010',11249+2*i)]:
        B,U=glue(l.lstrip('0'),R[0],m,r+extra)
        T=apex(B,U)
        check(strong(T) and H(T)==N)
def prime(n):
    return n>1 and all(n%d for d in range(2,isqrt(n)+1))
check(11249==7*1607 and prime(1607) and prime(11251))

table=[(i,u,*profile(r,u)[2:4]) for i,r in enumerate(R) for u in range(32)]
def inverse(f,P,S):
    check(f%84==0)
    b=(f//84).bit_length()-1
    for m in range(b+1):
        D=1<<m
        y=P*pow(9,-1,D)%D
        z=(P-9*y)//D
        x,r=divmod(z,84)
        if not 0<=x<(1<<(b-m)) or r>75:continue
        for i,U,p,s in table:
            if r==p and S==(9*(y%2) if m else s):
                yield i,U,x,m,y
check({(0,m,int(l or '0',2),len(r),int(r,2)) for l,m,r in choices}
      <=set(inverse(672,245,9)))
print('Four literal full windows, including the strong prime/7p pair; profile inverse controls passed.')
