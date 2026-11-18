"""Physical count-window overlap across a product-changing, fixed-sum carry."""
import sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE.parent/'figurate_overlap'),str(HERE.parent/'half_reflection'),
             str(HERE.parent/'joint_growth'),str(HERE.parent/'ans_inverse')]
from folded_centers import core,step,digits,word,CAP,G
from joint_decode import checksum
from reflection_replay import arm,rooted_last,glue,rev,transitive_arm
from math import gcd

def center(S,x,y):
    h,g,B,U,V,Ub,Vb=S
    return h+B*x*y+Ub*x+U*y+g*(checksum(x)+checksum(y))-Vb*(x&1)-V*(y&1)

def arithmetic(m,eta=1):
    h,g,B,U,V,Ub,Vb=core(m)
    P=step(m);M=P//2
    assert U==Ub and g==G*M and gcd(B,G)==1 and B%2==0
    assert -V-Vb==1040*P
    L0=pow(4,-1,M)
    k0=(4*L0-1)//M
    c=B//2*k0+G*(L0-3-4*eta)+1040
    A=2*B+g
    d=gcd(A,9*G)
    assert d in (1,3,9)
    delta=(c+d//2)%d-d//2
    assert delta==(-3-eta if m%3==0 else -eta)
    T=(9*G)//d
    s0=((delta-c)//d)*pow(A//d,-1,T)%T
    L=L0+M*s0
    z=(c+A*s0-delta)//(9*G)
    dL=M*T;dz=A*T//(9*G)
    db=4*dz
    slack=4*z+2*eta+1-2*L
    assert db>2*dL
    v=max(0,-(slack//(db-2*dL)))
    return delta,(L+v*dL,z+v*dz),(dL,dz)

def addresses(m,n,eta=1):
    delta,(L,z),(dL,dz)=arithmetic(m,eta)
    L+=n*dL;z+=n*dz
    b=4*z+2*eta+1;a=b-L
    assert 1<=a<=b<=2*a
    return (4*a,4*b),(4*a+1,4*b-1),delta

def witnesses(m,N):
    S=core(m);h,_,B,U,_,Ub,_=S
    def en(n):
        (x,y),_,_=addresses(m,n)
        return h+B*x*y+Ub*x+U*y
    if N<en(0):return ()
    a,b=0,1
    while en(b)<=N:a,b=b,2*b
    while b>a+1:
        c=(a+b)//2
        if en(c)<=N:a=c
        else:b=c
    E,F,_=addresses(m,a)
    out=[]
    for x,y in (E,F):
        q,r=divmod(N-center(S,x,y),step(m))
        d=digits(m,q)
        if r==0 and d is not None:out.append((a,x,y,*d))
    return tuple(out)

def main():
    A,r=CAP;C=rooted_last(A,r);S=arm(A,r)
    for x,y in ((4,4),(5,3)):
        T=glue(glue(rev(transitive_arm(x)),C),transitive_arm(y))
        assert arm(T,len(T)-1)[0]==center(S,x,y)
    for m in (1,2,3,9,65,66,67,128):
        P=step(m);S=core(m)
        for n in (0,1,2,27,10**30):
            (x,y),(a,b),d=addresses(m,n)
            assert x+y==a+b and x*y!=a*b
            assert center(S,a,b)-center(S,x,y)==d*P
            if m>=65:
                N=max(center(S,a,b),center(S,x,y))+2193*P
                Z=witnesses(m,N)
                assert len(Z)==2 and {r[0] for r in Z}=={n}
                assert not witnesses(m,N+2)
        E,F,d0=addresses(m,0,0)
        assert center(S,*F)-center(S,*E)==d0*P
        print('depth',m,'odd/even signed shifts',d,d0,'repair steps')
    # Two different literal arc tables at one jointly certified target.
    m=65
    (x,y),(a,b),d=addresses(m,0)
    N=max(center(core(m),a,b),center(core(m),x,y))+2193*step(m)
    for _,X,Y,j,k in witnesses(m,N):
        T=word(m,X,0,j,k,Y)
        assert all(not (A>>i&1) and A.bit_length()<=len(T) for i,A in enumerate(T))
        print('compiled order',len(T),'address widths',X.bit_length(),Y.bit_length(),
              'target decimal digits',len(str(N)))
if __name__=='__main__':main()
