"""Variable reflected-arm exact three-reading inverse."""
from functools import lru_cache
from itertools import product
from fractions import Fraction
from math import comb, prod, gcd, isqrt
from random import Random
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'moonshine_prime_blocks'))
from clone_replay import vectors,TT,read,inflate
from frame_replay import cut,join
from replay import R as READERS,reads
from gap_replay import rand_t
A=(10,4,9,2);B=cut(A,8); MARK=(0,2,3)
L=(80,397,505,161,170,451,282,321,25)
R=(274,68,129,311,358,7,301,379,38)
W=TT(2)

@lru_cache(None)
def f(T,q):
    return read(T,q) if q else 0

def partition(T,q):
    return sum((-1)**(q-j)*comb(q,j)*f(T,j) for j in range(1,q+1))

def evaluate(T,cores,q):
    C=vectors(T,MARK,q,q+1)
    P={k:[0]+[partition(V,i) for i in range(1,q+2)] for k,V in enumerate(cores)}
    return sum(w*P[0][a]*P[1][b]*P[2][c] for (a,b,c),w in C.items())

def skew(L,R,W,q):
    x,y,z=(f(U,q) for U in (L,R,W))
    return q*sum((f(W,q+i)-z)*(x*f(R,q+i)-y*f(L,q+i)) for i in (1,-1))

def run():
    rng=Random(911306)
    for _ in range(43):
        C=tuple(rand_t(rng.randrange(1,4),rng) for i in MARK)
        for q in (1,2,3):
            x,y=(f(inflate(T,dict(zip(MARK,C))),q) for T in (A,B))
            assert x-y==skew(*C,q)
            assert tuple(evaluate(T,C,q) for T in (A,B))==(x,y)
    data=tuple(tuple(f(U,q) for q in (1,2,3)) for U in (L,R,W))
    assert data==((1223,23716,222309),(1073,19516,182979),(1,4,9))
    (h,b,c),(x,y,z),(j,k,m)=data
    assert Fraction(c,b)-Fraction(3*h,5*b)==Fraction(z,y)-Fraction(3*x,5*y)==Fraction(327,35)
    a0,a1=1729645521,1734383121
    assert tuple(evaluate(T,(L,R,W),1) for T in (A,B))==(a0,a1)
    assert tuple(evaluate(T,(L,R,W),2) for T in (A,B))==(630250594716,)*2
    for U in (L,R,W):assert len(U)>0
    for T in (A,B):
        S=inflate(T,{0:L,2:R,3:W})
        for k in (S,tuple(sum(1<<i for i,r in enumerate(S) if r>>j&1) for j in range(len(S)))):
            seen=1
            while True:
                temp=seen
                for i,row in enumerate(k):
                    if seen>>i&1:temp|=row
                if temp==seen:break
                seen=temp
            assert seen==(1<<len(S))-1
    print('129 literal and compressed comparisons; neutral two-vertex arm',a0,a1)
    P=join((92,105,42,224,302,129,52,343,111),(30,12,24,16,2))
    Q=join((126,220,24,112,224,198,132,13),(26,28,57,48,32,3))
    assert tuple(tuple(f(U,q) for q in (1,2,3)) for U in (P,Q))==(
        (8275,1831168,88274043),(7775,1512704,72922599))
    base=112463595873847652
    out=11010786697063
    assert tuple(evaluate(T,(P,Q,W),1) for T in (A,B))==(out,11015945813863)
    assert tuple(evaluate(T,(P,Q,W),2) for T in (A,B))==(base,base)
    for scale in (1,7):
        D,I,a,b=READERS[1,scale]
        h,g,q=reads(D,I)
        assert (g,h+q-g)==(6*scale,scale)==(a,b)
    p=12*base+out
    assert p==1349574161272868887 and 7*p==9447019128910082209
    fs=(2,3,2377,3433,5881,4686961)
    assert prod(fs)==p-1 and pow(3,p-1,p)==1
    assert all(all(q%d for d in range(2,isqrt(q)+1)) and gcd(pow(3,(p-1)//q,p)-1,p)==1 for q in fs)
    print('prime and paired exact count',p,7*p)
    P=join((140,241,186,226,9,465,149,272,95),(30,12,24,16,2))
    Q=join((126,188,168,112,100,192,134,25),(26,28,57,48,32,3))
    W2=(6,12,8,1)
    assert tuple(tuple(f(U,q) for q in (1,2,3)) for U in (P,Q,W2))==(
        (8635,1868032,89124867),(8375,1629440,77741775),(5,32,117))
    assert tuple(evaluate(T,(P,Q,W2),1) for T in (A,B))==(96833175110435,96875688057635)
    assert tuple(evaluate(T,(P,Q,W2),2) for T in (A,B))==(1578628339040509540,)*2
    print('nontransitive arm: 96833175110435, 96875688057635')
if __name__=='__main__':run()
