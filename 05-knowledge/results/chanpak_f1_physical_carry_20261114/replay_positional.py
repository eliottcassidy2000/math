"""Exact 5-adic positional alphabets and optimal width in the fixed h=5 model."""
from itertools import product
from collections import defaultdict
from random import Random
from replay_f1 import Core,path_arrays,make,audit
from compile_absolute import binary_h as g

A=(8,1,3,6)
def masks(a,b):
    return tuple((1<<a)-(1<<(i+1)) for i in range(a))
def rev(x,k):
    y=0
    for _ in range(k): y=2*y+(x&1);x//=2
    return y
def puncture(a,b,x,z,U):
    R=list(masks(a,a))
    n=a+4+b
    for i in range(a):R[i]|=(1<<n)-(1<<a)
    R += [(r<<a)|((1<<n)-(1<<(a+4))) for r in A]
    R += [r<<(a+4) for r in masks(b,b)]
    return tuple(R),rev(x,a)|(U<<a)|(rev(z,b)<<(a+4))

DATA={
 0:(7,3232,0,0,[(1,6,1,32,12,4000,3169),(7,0,101,0,0,3676,2917),(6,1,50,0,12,3292,2613),
                  (2,5,3,0,1,4048,3205),(4,3,12,0,7,3724,2951)]),
 1:(7,3491,5,0,[(0,7,0,7,15,480,415),(4,3,13,7,12,276,245),(1,6,1,7,14,462,401),
                  (6,1,54,1,6,948,819),(2,5,3,7,9,444,385)]),
 2:(8,1438,0,5,[(2,6,0,6,14,5035,893),(6,2,11,2,1,5431,967),(7,1,22,0,6,6007,1067),
                  (5,3,5,6,12,4663,831),(3,5,1,6,9,4999,887)]),
 3:(7,600,0,5,[(4,3,2,0,3,3815,569),(0,7,0,120,0,3611,535),(2,5,0,24,6,3647,541),
                  (6,1,9,0,9,3143,471),(1,6,0,56,1,3629,539)]),
 4:(7,1438,0,5,[(1,6,0,6,14,2515,893),(6,1,22,0,6,3001,1067),(2,5,1,6,9,2497,887),
                  (5,2,11,2,1,2713,967),(4,3,5,6,12,2329,831)])}
def alphabet(r):
    D,P,S,T,rows=DATA[r];f=32<<D
    return [Core(*puncture(a,b,x,z,U),5,f,P,Q,S,T,K) for a,b,x,z,U,Q,K in rows]

def terms(x,k):
    C=alphabet(x%5)[0]
    p=[x]
    for j in range(k):p.append(C.f*p[-1]+C.P*5**j)
    return C,p

def add(C,j,x,p,d):
    return 5**j*(d.K-5)+p[j]*d.Q+(x%2)*C.T*(j==0)

def count(x,y,digits):
    C,p=terms(x,len(digits));H=g(x)
    for j,d in enumerate(digits): H=5*H+add(C,j,x,p,d)
    k=len(digits)
    return H+5**k*(g(y)-1)+p[k]*y+5**(k-1)*C.S*(y%2)

def solve(N,x,y,k):
    assert k>=1 and x>0 and y>0
    digs=alphabet(x%5)
    C,p=terms(x,k)
    z=N-5**k*(g(y)-1)-p[k]*y-5**(k-1)*C.S*(y%2)
    word=[]
    for j in reversed(range(k)):
        opts=[d for d in digs if (z-add(C,j,x,p,d))%5==0]
        assert len(opts)==1
        d=opts[0]
        z=(z-add(C,j,x,p,d))//5
        word.append(d)
    return word[::-1],z

def mins():
    h,starts=path_arrays(A)
    _,ends=path_arrays(tuple((15^r)^(1<<i) for i,r in enumerate(A)))
    P=[sum(row[i]*h[15^s] for s,row in enumerate(starts)) for i in range(4)]
    Q=[sum(row[i]*h[15^s] for s,row in enumerate(ends)) for i in range(4)]
    out=[]
    for U in range(16):
        L=[i for i in range(4) if U>>i&1];R=[i for i in range(4) if not U>>i&1]
        out.append((sum(P[i] for i in L),sum(Q[i] for i in R),
                   sum(ends[-1][i] for i in L),sum(starts[-1][i] for i in R)))
    for D in range(8):
        G=defaultdict(set)
        for a in range(D+1):
            b=D-a
            for x,z,U in product(range(1<<a),range(1<<b),range(16)):
                p,q,s,t=out[U]
                y=(1<<b)-1-rev(z,b)
                F=(1<<b)*(32*x+p)+5*z
                Gq=(1<<a)*(32*y+q)+5*((1<<a)-1-rev(x,a))
                S=5*(z%2) if b else s
                T=5*(1-(x>>(a-1))) if a else t
                G[F,S,T].add(Gq%5)
        yield len(G),sum(len(v)==5 for v in G.values())

if __name__=="__main__":
    rng=Random(25011)
    for r in range(5):
        digits=alphabet(r)
        for d in digits:audit(d)
        assert len({d.Q%5 for d in digits})==len({(d.K+r*d.Q)%5 for d in digits})==5
        assert (digits[0].f*r+digits[0].P)%5
        for _ in range(60):
            x=5*rng.randrange(1,50)+r;y=rng.randrange(1,60);k=rng.randrange(1,7)
            D=rng.choices(digits,k=k)
            N=count(x,y,D)
            assert solve(N,x,y,k)==(D,g(x))
            assert solve(N+5**k,x,y,k)==(D,g(x)+1)
        # Literal DP check of a representative of each bank.
        x=5+r;y=1
        for d in digits[:1]: assert path_arrays(make([d],x,y))[0][-1]==count(x,y,[d])
    N=[]
    for d in alphabet(0):
        v=path_arrays(make([d],10,1))[0][-1]
        assert v==count(10,1,[d]);N.append(v)
    assert N==[87401,83909,79765,87917,84423]
    assert list(mins())==[(16,0),(62,0),(176,0),(416,0),(896,0),(1856,0),(3784,0),(7624,40)]
    print("PASS: 25 literal profiles, 10 full literal counts, all digit permutations, 300 inverse checks, minimum width")
