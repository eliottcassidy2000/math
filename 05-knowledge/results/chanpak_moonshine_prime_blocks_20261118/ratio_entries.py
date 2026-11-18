"""Physical ratio changes, target-factor entries, and a two-orientation return."""
import sys
from pathlib import Path
from itertools import product
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE.parent/q) for q in ('reflected_words','ans_inverse','exact_carry','figurate_overlap','joint_growth','half_reflection','tape_carry')]
from reflected_return import Return,SHORT
from reflection_search import find
from exact_product_carry import even_pair,make
from reflection_replay import transitive_arm
from fixtures import state,forward,dp,rev
from ans_replay import run

BASE=(((0,41,67,5,15,29,59),2),((0,65,83,7,11,31,57),6))
def prefix(S,T):
    n=len(S);m=len(T)
    return tuple(x|(((1<<m)-1)<<n) for x in S)+tuple(x<<n for x in T)
def cover(S):
    h=[1]+[sum(r) for r in dp(S)[1:]]
    f=len(h)-1
    return sum(h[j]*h[f^j] for j in range(f+1))
def deletion(T,r):
    I=[j for j in range(len(T)) if j!=r]
    return tuple(sum((T[j]>>i&1)<<k for k,i in enumerate(I)) for j in I)
def fixture(P,a,b,w):
    A,B=(state(*U) for U in P)
    h,e,t,d=A
    assert (B[0],B[1],B[3])==(h,e,d)
    assert t%a==0 and B[2]==b*(t//a)
    R=Return(a,b,w,d)
    H,u,j=R.macro(h,t//a,e,0)
    for ((T,r),bits,ratio,S) in zip(P,R.bitwords(0),(a,b),(A,B)):
        assert run(S,bits)==(H,j,ratio*u,d)
        if len(T)+len(bits)<=15:
            assert state(forward(T,r,bits),r)==(H,j,ratio*u,d)
    return H
def rotate(Z):
    H,u,v,B=Z
    return H,v,u,B
def transport(Z,j):
    H,u,v,B=Z
    h,a,_=SHORT.macro(H,u,3,j)
    return h,a,v+(B<<j), B<<(4+j)
def prev(Z):
    H,u,v,B=Z
    R=SHORT.inv(H,u,3)
    if R is None:return None
    h,a,_,j=R
    C,x=divmod(B,1<<(4+j))
    if x or C<1536 or v<=C<<j or h<863:return None
    return (h,a,v-(C<<j),C),j
def main():
    E,F=even_pair((1,1,1,0,0,0,0),1)
    A,B=(make((0,),*Q) for Q in (E,F))
    r=len(A)-1
    P=((A,r),(B,r));Q=((rev(B),r),(rev(A),r))
    assert tuple(state(*U) for U in P)==((1063,1,1847,1),(1063,1,1967,1))
    assert tuple(state(*U) for U in Q)==((1063,1,1105,1),(1063,1,1177,1))
    w,_=find(1105,1177)
    print('same arcs reversed; new 1063 continuation',fixture(Q,1105,1177,w), 'word length',w.bit_length())
    S=(2,4,1);s=3;C=cover(S)
    assert C==12
    for (T,r),(R,_) in zip(P,Q[::-1]):
        U=prefix(S,T);r+=len(S)
        h,e,t,d=state(T,r-len(S))
        assert state(U,r)==(s*h,s*e,s*t,s*d)
        assert state(rev(U),r)==(s*h,0,C*state(R,r-len(S))[2],s*d)
    X=transitive_arm(1154);Y=rev(X);r=len(X)-1
    assert state(Y,r)==(863,1,1155,1) and state(X,r)==(863,1,1527,1)
    assert find(385,509)[0]==260099
    print('targeted alternative 863 entry',fixture(((Y,r),(X,r)),385,509,260099))
    D=(110,124,120,48,97,64,8);q=4
    assert state(D,q)==(29,1,57,3)
    D1=forward(D,q,(0,))
    assert state(D1,q)==(29,0,114,3)
    print('targeted 29 switch keeping common branch',fixture((BASE[1],(D1,q)),13,19,57))
    P=tuple((forward(T,r,w),r) for (T,r),w in zip(BASE,SHORT.bitwords(0)))
    assert tuple(state(*U) for U in P)==((863,3,1023,3),(863,3,1209,3))
    Q=tuple((rev(T),r) for T,r in P[::-1])
    assert tuple(state(*U) for U in Q)==((863,3,1089,3),(863,3,1287,3))
    Z=(863,93,99,1536)
    for mode,P0 in enumerate((P,Q)):
        before=rotate(Z) if mode else Z
        h,u,v,C=transport(before,0)
        assert prev((h,u,v,C))==(before,0)
        A1=tuple((forward(T,r,w),r) for (T,r),w in zip(P0,SHORT.bitwords(0)))
        for i,((T,r),a,b) in enumerate(zip(A1,(11,13),(13,11))):
            assert state(T,r)==(h,3,a*u,3)
            assert state(rev(T),r)==(h,3,b*v,3)
            assert cover(deletion(T,r))==C
        print('both physical 863 dual-scale channels',mode,h,u,v,C)
    for m in (0,1,0,1):
        V=rotate(Z) if m else Z
        W=transport(V,m)
        assert prev(W)==(V,m)
        Z=W
if __name__=='__main__':main()
