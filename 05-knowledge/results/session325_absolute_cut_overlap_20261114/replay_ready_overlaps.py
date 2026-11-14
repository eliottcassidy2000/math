from itertools import permutations

def req(x):
    if not x:raise RuntimeError('replay mismatch')

def end(A):
    n=len(A);D=[[0]*n for _ in range(1<<n)]
    for i in range(n):D[1<<i][i]=1
    for S in range(1,1<<n):
        for i,c in enumerate(D[S]):
            if c:
                a=A[i]&~S
                while a:
                    t=a&-a;a-=t;D[S|t][t.bit_length()-1]+=c
    return D

def converse(A):
    return tuple(((1<<len(A))-1)^r^(1<<i) for i,r in enumerate(A))

def relabel(A,p):
    return tuple(sum(1<<j for j in range(len(A)) if A[p[i]]>>p[j]&1) for i in range(len(A)))

def strong(A):
    for i in range(len(A)):
        S=1<<i
        for _ in A:
            for j,r in enumerate(A):
                if S>>j&1:S|=r
        if S!=(1<<len(A))-1:return False
    return True

def profile(A,v):
    E=end(A);S=end(converse(A));h=[1]+list(map(sum,E[1:]))
    C=((1<<len(A))-1)^(1<<v)
    L=h[C]
    m=C
    while m:
        L+=h[C^m]*sum(s for j,s in enumerate(S[m]) if A[v]>>j&1)
        m=(m-1)&C
    return h[-1],h[C],L

def attach(A,v,j):
    n=len(A)
    N=n+j+1
    out=[0]*N
    for i in range(n):
        out[i]=A[i]|(((1<<(j+1))-1)<<n)
    out[v]^=1<<(N-1)
    for i in range(n,N-1):
        out[i]=((1<<N)-1)^((1<<(i+1))-1)
    out[-1]=1<<v
    return tuple(out)

def check(A,v,want):
    req(strong(A) and profile(A,v)==want)

A=(2,48,35,7,13,25);B=(48,9,3,5,14,30)
req(all(relabel(A,p)!=B and relabel(converse(A),p)!=B for p in permutations(range(6))))
check(A,0,(33,15,36));check(B,4,(33,15,72))
check(attach(A,0,2),0,(177,15,273));check(attach(B,4,1),4,(177,15,273))
check(attach(attach(A,0,2),0,1),0,(723,15,1077))
check(attach(attach(B,4,1),4,1),4,(723,15,1077))

pair=[((8,17,3,6,13),3,13),((16,9,3,5,14),0,15)]
for T,v,H in pair:check(T,v,(H,3,15))
triple=[(72,17,67,38,13,23,58),(72,5,73,50,7,23,58),(66,24,67,37,13,23,58)]
for i,T in enumerate(triple):
    H=135+2*i
    check(T,5,(H,45,303))
    T=attach(attach(T,5,1),5,2)
    # L after two: 8*(4*303-45)-45
    check(T,5,(H+2*303+4*(4*303-45),45,8*(4*303-45)-45))

A=(8,1,3,6)
check(A,1,(5,3,6))
for s in (2,3):
    R=2*s+1
    X=attach(A,1,R-1)
    Y=attach(attach(A,1,s-1),1,s)
    h,x,L=profile(X,1);k,y,M=profile(Y,1)
    req(strong(X) and strong(Y) and h==k==5+6*4**s and x==y==3 and L!=M)
    p=sorted(r.bit_count() for r in X);q=sorted(r.bit_count() for r in Y)
    req(p==[1]+list(range(1,R+3))+[R+2] and q!=p and s+1 not in q)
print('nonconverse synchronization, ready pair/triple and asynchronous collisions verified')
