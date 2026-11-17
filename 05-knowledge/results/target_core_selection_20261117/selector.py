from math import isqrt
from functools import cache

def need(ok):
    if not ok: raise AssertionError

@cache
def F(x):
    if not x: return 1
    z=x//2
    return F(z)+(0 if x&1 else z+(z&1))

def count(b,e):
    return b*e+F(b-1)*F(e-1)-(b&1)*(e&1)

def ceildiv(a,b):
    return -((-a)//b)

def selector(N):
    for b in range(2,isqrt(N)+1):
        a=F(b-1)
        L=2
        while L*L*b<N: L*=2
        k=L.bit_length()-1
        for r in range(L):
            A=L-1-int(f'{r:0{k}b}'[::-1],2)
            delta=int(r<L//2)
            B=b*L+a*A
            M=N-b*(r+1)-a*(F(r)-1)+(b&1)*int(r%2==0)
            for q in range(max(0,ceildiv(M-2*a,B+a)),(M-a)//B+1):
                e=L*q+r+1
                if e>=b and count(b,e)==N:
                    yield b,e

def core(b):
    x=b-1
    n=x.bit_length()
    rows=[((1<<n)-1)^((2<<i)-1) for i in range(n)]+[0]
    for i in range(n):
        if x>>i&1: rows[i]|=1<<n
        else: rows[-1]|=1<<i
    return tuple(rows),n

def join(A,v,D,u):
    n=len(A)
    W=list(A)+[(sum(1<<j for j,r in enumerate(D) if r>>i&1)<<n)|((1<<n)-1) for i in range(len(D))]
    W[v]|=1<<(n+u); W[n+u]^=1<<v
    return W

def seam(b,e):
    A,v=core(b);D,u=core(e)
    return join(A,v,D,u)

def F_direct(x):
    n=x.bit_length()
    return 1+sum((1<<(j-i-1))+int(j==i+1)
                 for i in range(n) for j in range(i+1,n)
                 if not (x>>i&1) and (x>>j&1))

def J(A,v):
    n=len(A)
    W=list(A)+[((1<<n)-1)^(1<<v)]
    W[v]|=1<<n
    return W,n

def slide_checks():
    for b in range(2,10):
        for e in range(2,10):
            A,v=core(b);D,u=core(e)
            left,V=J(A,v); right,U=J(D,u)
            # Shift the new vertex from either side onto the seam.
            W=join(left,V,D,u)
            X=join(A,v,right,U)
            old_to_new=list(range(len(A)))+[len(A)+len(D)]+list(range(len(A),len(A)+len(D)))
            need(all((W[i]>>j&1)==(X[old_to_new[i]]>>old_to_new[j]&1)
                     for i in range(len(W)) for j in range(len(W))))
    print('64 literal one-vertex seam slides')

def H(rows):
    n=len(rows)
    @cache
    def d(S,v):
        R=S^(1<<v)
        return 1 if not R else sum(d(R,u) for u,r in enumerate(rows) if R>>u&1 and r>>v&1)
    return sum(d((1<<n)-1,v) for v in range(n))

def hadamard_local():
    T=(102,252,88,49,225,68,136,45)
    U=(34,180,121,19,225,72,139,45)
    M=[[1 if i==j else (2*(T[i]>>j&1)-1 if i<j else 1-2*(U[j]>>i&1)) for j in range(8)] for i in range(8)]
    need(all(sum(x*y for x,y in zip(M[i],M[j]))==8*int(i==j) for i in range(8) for j in range(8)))
    need(all((R[i]>>(i+1)&1) for R in (T,U) for i in range(7)))
    X={(i,j) for i in range(8) for j in range(i+1,8) if (T[i]^U[i])>>j&1}
    need(X=={(0,2),(0,6),(1,3),(1,6),(2,5),(3,5)})
    need((H(T),H(U))==(345,469))
    print('Hadamard patch:',sorted(X),345,469)
    T=(3918,2732,360,3888,2023,3009,2442,3341,514,1220,2406,788)
    U=(4062,2732,1896,3504,1958,3009,410,3332,3586,1224,2146,596)
    n=12
    M=[[1 if i==j else (2*(T[i]>>j&1)-1 if i<j else 1-2*(U[j]>>i&1)) for j in range(n)] for i in range(n)]
    need(all(sum(x*y for x,y in zip(M[i],M[j]))==n*int(i==j) for i in range(n) for j in range(n)))
    need(all((R[i]>>(i+1)&1) for R in (T,U) for i in range(n-1)))
    X={(i,j) for i in range(n) for j in range(i+1,n) if (T[i]^U[i])>>j&1}
    V={x for pair in X for x in pair}
    need(len(V)==10 and all(sum(x in pair for pair in X)==2 for x in V))
    reached={min(V)}
    for _ in V:
        reached|={j if i in reached else i for i,j in X if i in reached or j in reached}
    need(reached==V)
    need((H(T),H(U))==(194243,181161) and 194243==7*27749)
    need(all(27749%d for d in range(2,isqrt(27749)+1)))
    print('Hadamard patch:',sorted(X),194243,181161)


if __name__=='__main__':
    cap=1500
    brute={}
    for b in range(2,isqrt(cap)+1):
        for e in range(b,cap//b+1):
            c=count(b,e)
            if c<=cap: brute.setdefault(c,set()).add((b,e))
    for N in range(1,cap+1,2):
        need(set(selector(N))==brute.get(N,set()))
    for p in [11,13,17,29,89,109,173,181,241,613,877,30011,30013,1000000007,1000000009]:
        reads=[]
        for N in (p,7*p):
            w=next(selector(N),None)
            if w:
                b,e=w
                need(b*e+F_direct(b-1)*F_direct(e-1)-(b&1)*(e&1)==N)
                if len(seam(*w))<=15: need(H(seam(*w))==N)
            reads.append(w)
        print(p,*reads)
    slide_checks()
    hadamard_local()
