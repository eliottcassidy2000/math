from itertools import product, permutations, combinations
from math import comb

def req(v):
    if not v:raise RuntimeError('replay mismatch')

def greedy(W,r):
    a=[0]*(r+1)
    for j in range(r,-1,-1):a[j],W=divmod(W,4**j)
    return a

def actual(W,r):
    if not r:
        yield (W,);return
    for c in range(1+W//4**r):
        for v in actual(W-c*4**r,r-1):yield v+(c,)

def f(c,n,d,j):
    def Q(z):return z**n+2*z**d*sum(a*z**(2*i) for i,a in enumerate(c))
    return sum((-1)**(j-k)*comb(j,k)*Q(k) for k in range(j+1))

for r in range(1,4):
    n=2*r+3;d=1
    for W in range(200):
        a=greedy(W,r);G=sum(a);X=list(actual(W,r))
        req({sum(v) for v in X}==set(range(G,W+1,3)))
        choices=[]
        for j in range(1,r+1):
            if not a[j]:continue
            c=a.copy();c[j]-=1;c[j-1]+=4
            choices.append((j,tuple(c)))
        req({c for _,c in choices}=={v for v in X if sum(v)==G+3})
        for (j,A),(k,B) in combinations(choices,2):
            req(all(f(A,n,d,v)-f(B,n,d,v)>=0 for v in range(1,n+1)))
            req(f(A,n,d,3)>f(B,n,d,3))

def identities(A):
    n=len(A);N=1<<n;M=N-1
    good=lambda p:all(A[i]>>j&1 for i,j in zip(p,p[1:]))
    sets=[tuple(i for i in range(n) if S>>i&1) for S in range(N)]
    paths=[[p for p in permutations(X) if good(p)] for X in sets]
    H=len(paths[-1]);C=sum(A[p[-1]]>>p[0]&1 for p in paths[-1])//n
    Q=sum(len(paths[s])*len(paths[M^s]) for s in range(N));B0=0
    for S in range(1,M):
        if S>(M^S):continue
        for p in paths[S]:
            for q in paths[M^S]:
                B0+=not(A[p[-1]]>>q[0]&1 or A[q[-1]]>>p[0]&1)
    req(Q==2*n*H-n*(n-1)*C+2*B0)
for A in ((2,4,1),(8,1,3,6),(10,28,25,16,1),(36,109,8,81,103,76,5)):
    identities(A)
print('base-four exact fibres, first-carry total order and seam identity pass')
