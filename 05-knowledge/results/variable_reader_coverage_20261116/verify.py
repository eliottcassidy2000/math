from functools import lru_cache
from itertools import combinations,product
from fractions import Fraction
from math import isqrt

A0=((2,4,25,19,35,15),(6,4,24,19,35,15),(10,4,17,22,35,15))

def P(a,x=3):
    n=len(a);F=(1<<n)-1
    @lru_cache(None)
    def f(S,i):
        T=S^(1<<i)
        return sum(f(T,j) for j in range(n) if T>>j&1 and a[i]>>j&1) if T else 1
    @lru_cache(None)
    def h(S):
        return sum(f(S,j) for j in range(n) if S>>j&1) if S else 1
    R=S=0
    for M in range(1,F+1):
        if M>>x&1:
            z=f(M,x)*h(F^M)
            if M.bit_count()%2:S+=z
            else:R+=z
    return h(F),h(F^(1<<x)),R,S

def append(a,e,x=3):
    n=len(a)
    return tuple(row|((int(i!=x or e))<<n) for i,row in enumerate(a))+((0 if e else 1<<x),)

def vr(b,n):
    return int(f'{b:0{n}b}'[::-1],2)

def J(b):
    t=[i for i in range(b.bit_length()) if b>>i&1]
    return sum((1<<(j-i-1))+(j==i+1) for i,j in combinations(t,2))

def formula(b,n,z=82,g=9,I=-11):
    return I+z*(2**(n+1)-vr(b,n))+g*(2*b-J(b))

def anchor(N):
    D=1<<((N//88).bit_length()-1)
    n=D.bit_length()-1
    z=Fraction(N,D)
    x,y=((0,2-z/82) if z<=164 else ((z/2-82)/9,0))
    k=n//2;Q=1<<k
    u,v=(min(Q-1,int(Q*w)) for w in (x,y))
    b=(u<<(n-k))|vr(v,k)
    return formula(b,n)

if __name__=='__main__':
    assert [P(a) for a in A0]==[(23,9,24,26),(25,9,24,26),(29,9,23,27)]
    assert [sum((x^y).bit_count() for x,y in zip(a,b))//2
            for a,b in zip(A0,A0[1:])]==[1,3]
    for n in range(1,6):
        for w in product((0,1),repeat=n):
            b=int(''.join(map(str,w)),2)
            for j,a in enumerate(A0):
                a=append(a,0)
                for e in w+(0,):a=append(a,e)
                assert P(a)[0]==formula(b,n,I=-11+2*j)
    for m in range(1,8):
        L=1<<m;D=L*L
        rv=[vr(b,2*m) for b in range(D)]
        f=[18*b-9*J(b) for b in range(D)]
        mx=max(len({(164*L+18*p)*(2*D-rv[b])+f[b] for b in range(D)})
               for p in range(L))
        assert mx>=D*(164-18)/164
    for N in range(100000,1_000_000,111):
        assert abs(anchor(N)-N)<400*isqrt(N)
    print('All seed profiles, literal count triples, collision lower bound, and target bracketing pass.')

FIVE=((2,260,273,7,43,15,319,383,57),(0,257,267,35,271,87,415,319,41),
      (0,257,323,23,39,271,315,383,25),(272,5,9,259,14,31,319,383,54),
      (0,257,11,275,71,287,431,319,21))
def block(a,b,near,x=3):
    n=len(a);m=len(b); buf=((1<<m)-1)<<n
    orig=tuple(z|(0 if near and i==x else buf) for i,z in enumerate(a))
    return orig+tuple((z<<n)|(1<<x if near else 0) for z in b)

if __name__=='__main__':
    for j,a in enumerate(FIVE):
        h,g,r,s=P(a,8)
        assert (g,r+s-g,h-2*s+2*g)==(15,167,-75+2*j)
        assert P(append(a,0,8),8)[1:]==(15,167,182)
        c=append(a,0,8)
        for e in (1,0,1,0):c=append(c,e,8)
        assert P(c,8)[0]==formula(5,3,z=334,g=15,I=-75+2*j)
    tri=(2,4,1)
    a=append(A0[0],0)
    h,g,r,s=P(a);Z=r+s-g;I=h-2*s+2*g
    # H(C3)=3 and q_B=12 (two empty choices of 3, six proper splits)
    h1,g1,r1,s1=P(block(a,tri,True))
    assert (h1,g1,r1+s1-g1,h1-2*s1+2*g1)==(3*I+12*Z,3*g,12*Z,3*I)
    h2,g2,r2,s2=P(block(a,tri,False))
    assert (h2,g2,r2,s2)==(3*h,3*g,6*(r+s),6*(r+s))
    print('Five-reader bank, literal appended tapes and exact opaque-block transfer verified.')
