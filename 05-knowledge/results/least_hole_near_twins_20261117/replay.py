from math import isqrt
from functools import cache

def need(x):
    if not x: raise AssertionError

def trans(n):
    return tuple(((1<<n)-1)^((2<<i)-1) for i in range(n))

def prepend(X):
    return ((1<<(len(X)+1))-2,)+tuple(x<<1 for x in X)

def lamp(X):
    n=len(X)
    return ((1<<(n+1))-2,)+tuple((x<<1)|(1<<(n+1)) for x in X)+(1,)

def table(rows):
    n=len(rows);U=(1<<n)-1
    @cache
    def end(S,j):
        R=S^(1<<j)
        return 1 if not R else sum(end(R,i) for i,r in enumerate(rows) if R>>i&1 and r>>j&1)
    @cache
    def H(S):
        return sum(end(S,j) for j in range(n) if S>>j&1) if S else 1
    return U,H,end

def H(T):
    U,h,_=table(T);return h(U)

def profile(T,v):
    U,h,d=table(T)
    rev=tuple(sum(1<<i for i,r in enumerate(T) if r>>j&1) for j in range(len(T)))
    _,_,start=table(rev)
    need(d(U,v)>0)
    return h(U),sum(d(S,v)*h(U^S) for S in range(1,U+1) if S>>v&1),start(U,v)

def p2(T):
    U,h,_=table(T)
    return sum(h(S)*h(U^S) for S in range(U+1))

def seam(A,v,D,u):
    n=len(A)
    W=list(A)+[(sum(1<<j for j,r in enumerate(D) if r>>i&1)<<n)|((1<<n)-1) for i in range(len(D))]
    W[v]|=1<<(n+u); W[n+u]^=1<<v
    return W

def prime(p):
    return p>1 and all(p%i for i in range(2,isqrt(p)+1))

for j in (2,3,6):
    A=lamp(trans(j));D=lamp(trans(j-1))
    x=2**(2*j)+3*2**j+1
    need((H(seam(A,0,D,0)),H(seam(A,1,D,0)))==(x,x+2))
    need(prime(x) and prime((x+2)//7 if (x+2)%7==0 else x+2))
    print('near square',j,x,x+2)
for k,m,A0,v0,A1,v1,K,d in (
    (5,1,(16,1,19,7,10),0,(16,17,3,7,12),1,51,-13),
    (3,7,(56,17,3,6,12,30),2,(48,49,3,7,12,28),1,161,-35)):
    t=2**k
    from_left=lamp(trans(k))
    root=k+1
    need(profile(from_left,root)==(t+1,2*t-1,1))
    for A,v,bias in ((A0,v0,d),(A1,v1,d+2*m)):
        a,b,c=profile(A,v)
        need((a+2*b,a-b-c)==(K,bias))
        need(H(seam(from_left,root,A,v))==K*t+bias)
        need((K*t+bias)%m==0 and prime((K*t+bias)//m))
    print('common coefficient',K,d,d+2*m,'axis',m,'quotients', (K*t+d)//m,(K*t+d)//m+2)
X=trans(3)
E=(14,4,8,1)
# Find a literal order-four 5-path core if the first candidate is not the match.
from itertools import combinations,product
for orient in product((0,1),repeat=6):
    A=[0]*4
    for (i,j),s in zip(combinations(range(4),2),orient):
        if s:A[i]|=1<<j
        else:A[j]|=1<<i
    if (H(A),p2(A))==(5,32):
        E=tuple(A);break
need((H(E),p2(E))==(5,32))
A=lamp(prepend(X));D=lamp(E)
u,R,v,S=H(X),p2(X),H(E),p2(E)
a=2*R+u
N=H(seam(A,0,D,0))
need((profile(A,0),profile(A,1),profile(D,0))==((a,a,u),(a,a+u,R),(S+v,S+v,v)))
need(u*(S//2)-v*(R//2-u)==1)
need(u*N==(4*R*R-2*R*u-3*u*u)*v+4*a and H(seam(A,1,D,0))==N+2)
print('Farey point',E,(u,R),(v,S),N,N+2)
for N,A,B in (
    (1000000007,(262,268,304,21,195,91,143,303,120),(2,272,291,23,357,331,143,319,73)),
    (1000000009,(268,65,42,162,15,339,157,55,222),(4,33,42,83,71,409,167,287,95))):
    u,R=H(A),2*p2(A)
    v,S=H(B),p2(B)
    need((R-u)*(S-v)+2*u*v==N and prime(N))
    print('literal 20-vertex cyclic witness',N,(u,R),(v,S))
Q=(24,1,3,6,14)
need((H(Q),p2(Q))==(11,100))
print('decimal module',H(Q),p2(Q))
