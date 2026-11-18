"""Exact replay of cyclic crossed readers and the square-norm observer."""
from fractions import Fraction
from itertools import product
from math import comb, gcd, isqrt, prod
from random import Random
from pathlib import Path
import sys
# Supports running next to the earlier packet or from this session folder.
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'moonshine_prime_blocks'))
from clone_replay import inflate, read, vectors
from frame_replay import cut
from gap_replay import rand_t
from prime_support_replay import LOWER, UPPER
from replay import profile, reads

A=(10,4,9,2)
B=cut(A,8)
MARKS=(0,2)
PAIRS=(
 ((254,188,248,112,96,192,130,24),(126,124,248,48,224,64,8,107)),
 ((222,188,120,112,224,193,2,76),(62,92,248,176,96,66,137,51)),
)
CASES=((65,1952,21537,65,1808,19953),(349,5104,39789,389,5984,46629))
SCORES=(11,Fraction(85,11))
ANSWERS=((6699827,6709187,2230556648),(55432627,55329667,9106026664))
CERTS={
 6699827: (2,7,373,1283), 6709187: (2,11,17,17939),
 55432627: (2,3,9238771), 55329667: (2,3,7,661,1993),
 4516644657971: (2,5,47,1093,8792207),
 4516644555011: (2,5,7,7,131,1103,63793),
}

def prime_cert(p):
    fs=CERTS[p]
    assert prod(fs)==p-1 and all(f>=2 and all(f%r for r in range(2,isqrt(f)+1)) for f in fs)
    assert pow(2,p-1,p)==1 and all(gcd(pow(2,(p-1)//f,p)-1,p)==1 for f in set(fs))

def evaluate(T,L,R,q):
    C=vectors(T,MARKS,q,q+2)
    def q_nonempty(U,n):
        return sum((-1)**(n-j)*comb(n,j)*read(U,j) for j in range(1,n+1))
    X=[0]+[q_nonempty(L,k) for k in range(1,q+3)]
    Y=[0]+[q_nonempty(R,k) for k in range(1,q+3)]
    return sum(w*X[r]*Y[s] for (r,s),w in C.items())

def norm_observer(T,q):
    """Direct principal-minor expansion and arithmetic in the squarefree ring."""
    import sympy as sp
    n=len(T)
    K=sp.Matrix([[(T[i]>>j&1)-(T[j]>>i&1) for j in range(n)] for i in range(n)])
    full=K.row_join(sp.ones(n,1)).col_join((-sp.ones(1,n)).row_join(sp.zeros(1,1)))
    D=[1]+[0]*((1<<n)-1);N=[0]*(1<<n)
    for m in range(1,1<<n):
        J=[i for i in range(n) if m>>i&1]
        if len(J)%2:
            N[m]=int(full.extract(J+[n],J+[n]).det())
            if len(J)==3:
                cyclic=all(sum((T[i]>>j)&1 for j in J)==1 for i in J)
                assert N[m]==(9 if cyclic else 1)
        else:D[m]=int(K.extract(J,J).det())
    def conv(x,y):
        z=[0]*len(x)
        for s in range(len(x)):
            t=s
            while True:
                z[s]+=x[t]*y[s^t]
                if t==0:break
                t=(t-1)&s
        return z
    def power(x,n):
        y=[1]+[0]*(len(x)-1)
        for _ in range(n):y=conv(x,y)
        return y
    P=power([d+v for d,v in zip(D,N)],q)
    Q=power([d-v for d,v in zip(D,N)],q)
    S=[0]*len(P)
    for m in range(len(P)):
        v=m;S[m]=P[m]
        while v:
            S[m]-=Q[v]*S[m^v]
            v=(v-1)&m
    assert S[-1]%(1<<n)==0
    return S[-1]>>n

def select(p, bank, neutral=False):
    """An exact target-first inverse over a supplied bank of physical blocks."""
    catalog={}
    by_h={}
    for U in bank:
        h,b,c=(read(U,j) for j in (1,2,3))
        catalog.setdefault((h,b),[]).append((c,U))
        by_h.setdefault(h,[]).append((b,c,U))
    for h, lefts in by_h.items():
        for x in by_h:
            residual=p-2*h*x
            if residual<=0:continue
            for b,c,L in lefts:
                v,r=divmod(residual,2*b-h)
                if r:continue
                y=v+x
                for z,R in catalog.get((x,y),[]):
                    if not neutral or (c-h)*y==(z-x)*b:
                        return inflate(A,{0:L,2:R}), inflate(B,{0:L,2:R})
    return None

def run():
    rng=Random(830911)
    for t in range(22):
        L,R=(rand_t(rng.randrange(1,4),rng) for _ in range(2))
        for q in (1,2,3):
            x,y=[read(inflate(T,{0:L,2:R}),q) for T in (A,B)]
            assert (x,y)==tuple(evaluate(T,L,R,q) for T in (A,B))
            F=lambda T,k:read(T,k) if k else 0
            assert x-y==q*(F(L,q)*(F(R,q+1)-F(R,q-1))-(F(L,q+1)-F(L,q-1))*F(R,q))
        V=rand_t(rng.randrange(1,7),rng)
        for q in (1,2):assert norm_observer(V,q)==read(V,q)
    print('66 cyclic transfer and 44 Pfaffian-norm checks against literal counting')
    for i,(L,R) in enumerate(PAIRS):
        h,b,c,d,e,f=(read(T,q) for T,q in product((L,R),(1,2,3)))
        assert (h,b,c,d,e,f)==CASES[i]
        assert Fraction(c-h,b)==Fraction(f-d,e)==SCORES[i]
        x,y,bb=ANSWERS[i]
        assert x==2*h*d+(2*b-h)*(e-d) and y==2*h*d+(2*e-d)*(b-h)
        assert y-x==b*d-h*e
        assert (evaluate(A,L,R,2),evaluate(B,L,R,2))==(bb,bb)
        assert tuple(profile(inflate(T,{0:L,2:R})) for T in (A,B))==((x,bb),(y,bb))
        selection=select(x,[L,R],neutral=True)
        assert selection is not None
        assert tuple(evaluate(T,L,R,1) for T in (A,B))==(x,y)
        assert selection==(inflate(A,{0:L,2:R}),inflate(B,{0:L,2:R}))
        for p in (x,y):prime_cert(p)
        print('score',SCORES[i],'prime placement',x,y,'height',bb,
              'factor certificate',x-2*h*d,(2*b-h,e-d))
    for scale,(D,I) in [(1,LOWER),(7,UPPER)]:
        d,g,v=reads(D,I)
        assert (g,d+v-g)==(62*scale,scale)
    for p, expect in zip((62*8*ANSWERS[1][2]+h for h in ANSWERS[1][:2]),
                         (4516644657971,4516644555011)):
        assert p==expect;prime_cert(p)
        print('same-input paired reader',p,7*p)

if __name__=='__main__':run()
