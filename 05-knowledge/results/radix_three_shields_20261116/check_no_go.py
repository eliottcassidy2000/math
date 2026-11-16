"""Independent checks of the adjacent triangle shadow and no +2 law."""
from itertools import product
from random import Random

def combine(A,B):
    if A[1]==1:return B
    if B[1]==1:return A
    h,q,E,L,T,S,J=A
    H,Q,e,l,t,s,j=B
    return h*H,q*Q,q*e+H*E,h*l+Q*L,h*t,H*S,H*J+h*j+e*L+T*s

identity=(1,1,0,0,0,0,0)
def singleton(x):return 1,2,x,x,x,x,0
def cyclic(s):return 3,12,3*s,3*s,s,s,s*(s-1)
def tape(bits):
    p=identity
    for x in bits:p=combine(p,singleton(x))
    return p

def direct(rows,V):
    N=len(rows);F=(1<<N)-1
    def path(A):
        d=[[0]*N for _ in range(1<<N)]
        for i in range(N):d[1<<i][i]=1
        for m in range(1,F+1):
            for i in range(N):
                if not ((m>>i)&1):continue
                mask=m^(1<<i)
                if mask:d[m][i]=sum(d[mask][j] for j in range(N) if ((mask>>j)&1) and ((A[j]>>i)&1))
        return d
    d=path(rows); dual=path([sum(((r>>i)&1)<<j for j,r in enumerate(rows)) for i in range(N)])
    c=[1]+[sum(line) for line in d[1:]]
    e=[sum(v for j,v in enumerate(line) if j in V) for line in d]
    l=[sum(v for j,v in enumerate(line) if j in V) for line in dual]
    return c[F],sum(c[m]*c[F^m] for m in range(F+1)),sum(e[m]*c[F^m] for m in range(F+1)),sum(c[m]*l[F^m] for m in range(F+1)),e[F],l[F],sum(e[m]*l[F^m] for m in range(F+1))

def rows_with_triangle(bits,i):
    n=len(bits)+1
    rows=[sum(1<<j for j in range(t+1,n)) for t in range(n)]
    rows[i]^=1<<(i+2)
    rows[i+2]^=1<<i
    return rows

def mk_shadow(bits,i):
    u,v=bits[i:i+2]
    X,Y=tape(bits[:i]),tape(bits[i+2:])
    B=combine(combine(X,cyclic(u+2*v)),Y)
    mask={j for j in range(i) if bits[j]}|{i+j for j in range(u+2*v)}|{j+1 for j in range(i+2,len(bits)) if bits[j]}
    return B,mask

def check(bits,guard):
    N=len(bits)
    Z=tape(bits); E0,L0,T0,S0,J0=Z[2:]
    vals=[set(),set()]
    for i in range(N-1):
        u,v=bits[i:i+2];r=N-i-2
        X,Y=tape(bits[:i]),tape(bits[i+2:])
        B,mask=mk_shadow(bits,i)
        h,q,E,L,T,S,J=B
        if N<=5:
            assert B==direct(rows_with_triangle(bits,i),mask)
        assert E==3*E0 and L-3*L0==3*(1<<r)*(v-u)
        D=3*guard*(1<<r)+3*Y[2]+2*X[4]-Y[5]
        for b in (0,1):
            val=guard*L+J+b*S
            assert val-3*(guard*L0+J0+b*S0)==2*(u==0 and v==1)+(v-u)*(D+2*b*(i==0))
            vals[b].add(val)
    for vs in vals:
        assert not any(x+2 in vs for x in vs)

if __name__=='__main__':
    for n in range(2,12):
        for bits in product((0,1),repeat=n):check(bits,2)
    rng=Random(20261116)
    for n in range(12,31):
        for _ in range(10):check(tuple(rng.randrange(2) for _ in range(n)),1<<rng.randrange(1,8))
    print('Moving-triangle shadow and no +2: exhaustive bits to 11, random to 30; literal rows to 6 vertices.')
