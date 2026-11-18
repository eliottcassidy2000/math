"""Simultaneously address the center and the number of repaired modules."""
import sys
from pathlib import Path
from math import isqrt
from itertools import product
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE.parent/'half_reflection'),str(HERE.parent/'ans_inverse')]
from reflection_replay import arm,rooted_last,glue,transitive_arm,apex_increment
from fixtures import state,forward
LOW=(62,124,24,112,96,68,5),2
HIGH=(16,65,3,71,14,31,53),1
SEED=(65,17,224,72,-6,135,-11)
G=17; B=224
def join(X,Y):
    h,g,b,u,v,a,c=X
    H,k,d,U,V,A,C=Y
    return (k*h+g*H-g*k+u*A+v*C,g*k,b*d,
            d*u+g*U,g*V,b*A+k*a,k*c)
def rows():
    m=1; S=SEED
    while True:
        h,g,b,u,v,_,_=S
        L=u//(4*g)
        yield m,S,L
        m+=1;S=join(S,SEED)
def dvalue(I):
    return -((-I)//2)  # ceil(I/2)
def bound(S):
    h,g,_,u,*_=S
    R=u*u
    return (R+16*g-1)//(16*g),R//g
def checksum(I):
    c=0
    while I:
        J,r=divmod(I,2)
        if r==0:c+=J+(J&1)
        I=J
    return c
def endpoints(m,S,I):
    h,g,_,u,v,_,_=S
    return h+u*I+g*checksum(I)-v*(I&1),4*(g//G)
def decode(N):
    """Restricted high tape annuli, m>=2. No isomorphism census."""
    for m,S,L in rows():
        if m<2:continue
        lower,upper=bound(S)
        if N<lower:return None
        if N>upper:continue
        h,_,_,u,_,_,_=S
        I=(N-h)//u
        if not dvalue(L)<=I<=L:return None
        C,D=endpoints(m,S,I)
        j,r=divmod(N-C,D)
        return (m,I,j) if r==0 and 0<=j<=m else None
def material(m,I,digits):
    arms=tuple(rooted_last(*X) for X in (LOW,HIGH))
    assert len(digits)==m
    T=arms[digits[0]]
    for bit in digits[1:]:T=glue(T,arms[bit])
    return glue(T,transitive_arm(I))
def main():
    assert arm(*LOW)==SEED and arm(*HIGH)==(69,*SEED[1:])
    for m,S,L in rows():
        if m>7:break
        if m<2:continue
        assert S[4]==-6*G**(m-1)
        assert S[3]==72*(B**m-G**m)//(B-G)
        lo,hi=bound(S)
        for I in {dvalue(L),L,*range(dvalue(L),min(L+1,dvalue(L)+15))}:
            C,D=endpoints(m,S,I)
            assert checksum(I)==apex_increment(I)
            for j in (0,1,m):
                N=C+D*j
                assert lo<=N<=hi and decode(N)==(m,I,j)
                assert decode(N+2) is None  # same address, wrong step
        print('layer',m,'address count',L-dvalue(L)+1,'palette',m+1)
    m,S,L=next(x for x in rows() if x[0]==2)
    for I in (dvalue(L), L):
        C,D=endpoints(m,S,I)
        A=[]
        for bits in product((0,1),repeat=m):
            T=material(m,I,bits);r=len(T)-1
            assert arm(T,r)[0]==C+D*sum(bits)
            Z=state(forward(T,r,(0,)),r)
            assert Z[0]==C+D*sum(bits)
            A.append(Z[1:])
        assert len(set(A))==1
        print('literal stack',m,I,'counts',C,C+D,C+2*D)
    for m,S,L in rows():
        if m<2:continue
        if m>3:break
        for I in range(dvalue(L),L+1):
            C,D=endpoints(m,S,I)
            for j in range(m+1):
                N=C+j*D
                if all(N%i for i in range(2,isqrt(N)+1)):
                    assert decode(N)==(m,I,j)
                    print('prime selected inside growing window',N,(m,I,j))
                    return
    raise AssertionError('expected prime')
if __name__=='__main__': main()
