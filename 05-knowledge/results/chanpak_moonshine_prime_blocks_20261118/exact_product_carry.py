"""Closed-form equal-height carry, with an exact continuation difference."""
from math import isqrt,gcd
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE.parent/x) for x in ('tape_carry','figurate_overlap','joint_growth','half_reflection','ans_inverse')]
from distinct_products import center
from fixtures import state,forward
from reflection_replay import glue,rev,transitive_arm,rooted_last,arm
from folded_centers import core,step,word,G
from joint_decode import checksum

def even_pair(S,r):
    _,g,B,*_=S
    assert r>0 and (r*g)%4==1 and B>=g
    a=8*r*B-(r*g+1)//2
    b=8*r*B+r*g
    return (4*a,4*b),(4*a+2,4*b-2)

def paired_check(S,r):
    _,g,B,U,V,Ub,Vb=S
    E,F=even_pair(S,r)
    x,y=E;X,Y=F;w=y.bit_length()
    assert x<=y<=2*x and X<=Y<=2*X and x+y==X+Y
    assert X*Y-x*y==12*r*g and checksum(X)+checksum(Y)-checksum(x)-checksum(y)==-12*r*B
    assert center(S,X,Y)-center(S,x,y)==2*(Ub-U)
    assert (x.bit_length(),w)==(X.bit_length(),Y.bit_length())
    inv=lambda n:int(format(n,'0%db'%w)[::-1],2)
    assert inv(Y)-inv(y)==1<<(w-3)
    return (16*B-g)<<(w-3)

def make(Q,x,y):
    return glue(glue(rev(transitive_arm(x)),Q),transitive_arm(y))

def rounded(m,L,eta):
    S=core(m);_,g,B,*_=S;P=step(m);M=P//2
    assert L>=2 and (2*L-1)%M==0
    eps=(1-L)&1
    K=2*B*((2*L-1)//M)+G*((L-1-eps-2*eta)//2)
    z=(K+(3*G)//2)//(3*G)
    b=4*z+2*eta+1;a=b-L
    E=(4*a,4*b);F=(4*a+2,4*b-2)
    assert b>=2*L and 1<=a
    delta=K-3*G*z
    assert abs(delta)<=1657 and center(S,*F)-center(S,*E)==delta*P
    return delta

def marked_suffix(T,T1):
    root=len(T)-1
    from ans_replay import run
    A,B=(state(forward(U,root,(0,)),root) for U in (T,T1))
    s=B[2]-A[2]
    assert A[0]==B[0] and A[1]==B[1]==0 and A[3]==B[3] and s==240
    for W in range(256):
        bits=tuple(W>>i&1 for i in range(W.bit_length()))
        a,b=run(A,bits),run(B,bits)
        assert b[0]-a[0]==s*W
        if W in (1,3):
            a1,b1=(state(forward(U,root,(0,*bits)),root) for U in (T,T1))
            assert (a,b)==(a1,b1)
            print('physical suffix address',W,'counts',a[0],b[0])

def main():
    S=(1,1,1,0,0,0,0);Q=(0,)
    E,F=even_pair(S,1)
    T,T1=[make(Q,*A) for A in (E,F)]
    A=state(T,len(T)-1);B=state(T1,len(T1)-1)
    assert (A,B)==((1063,1,1847,1),(1063,1,1967,1))
    assert B[2]-A[2]==paired_check(S,1)
    assert all(A[0]%d for d in range(2,isqrt(A[0])+1))
    R=make(Q,64,114)
    assert state(R,len(R)-1)[0]==7*A[0]
    print('literal masks',T,T1,R)
    print('literal prime and upper witnesses',A,B,7*A[0])
    marked_suffix(T,T1)
    # Nonsymmetric literal control, r=1 for g=17.
    raw=(62,124,24,112,96,68,5);root=2
    V=arm(raw,root)
    assert V[3]!=V[5] and paired_check(V,1)>0
    for m in (0,1,2,65,66,67,100):
        W=core(m);g=W[1]
        for r in (pow(g,-1,4),pow(g,-1,4)+4,pow(g,-1,4)+4*10**20):
            assert paired_check(W,r)>0
        if m>=65:
            M=step(m)//2
            slope=2*paired_check(W,pow(g,-1,4))
            assert gcd(slope,2*M)==2
            for R in (0,1,17,M//2,M-1):
                s=R*pow(slope//2,-1,M)%M
                assert slope*s%(2*M)==2*R
            for j in range(6):
                L=pow(2,-1,M)+M*j
                for eta in (0,1):rounded(m,L,eta)
            print('equal-height and all-gap controls at depth',m)
if __name__=='__main__':main()
