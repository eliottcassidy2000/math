"""Two independent repairs -> exact gap-free growing window in one lattice."""
from joint_decode import join, checksum, material, dvalue
from reflection_replay import arm,rooted_last,glue,transitive_arm
from fixtures import state
from itertools import product,combinations
from math import gcd

ARMS=(((38,108,56,113,99,64,5),4),((96,1,99,7,15,26,58),4))
CORES=(((62,124,24,112,96,68,5),2),((16,65,3,71,14,31,53),1))
UNIT=(9689,493,71680,26406,-340,24747,-319)
_,G,B,U,V,*_=UNIT
SINKS=rooted_last((6,4,0),0)  # root then two ordinal sinks
PADDED=join(UNIT,(1,1,4,3,-1,0,0))
assert PADDED==(9689,493,286720,107103,-493,24747,-319)


def rows(base=UNIT):
    m=1;S=base
    while True:
        u,g=S[3],S[1]
        L=u//(4*g)
        yield m,S,L
        m+=1;S=join(S,base)

def bounds(S):
    u,g=S[3],S[1]
    return (u*u+16*g-1)//(16*g),u*u//g

def endpoints(m,S,I):
    h,g,_,u,v,_,_=S
    return h+u*I+g*checksum(I)-v*(I&1),4*(g//G)

def digits(m,r):
    j0=(12*r)%29
    k0=(r-17*j0)//29
    a=max(0,(k0-m+16)//17)
    j,k=j0+29*a,k0-17*a
    return (j,k) if 0<=j<=m and 0<=k<=m else None

def decode(N,base=UNIT):
    for m,S,L in rows(base):
        if m<2:continue
        lo,hi=bounds(S)
        if N<lo:return None
        if N>hi:continue
        h,g,_,u,*_=S
        I=(N-h)//u
        if not dvalue(L)<=I<=L:return None
        assert gcd(u,g)==1
        if (N-h-u*I)%(g//G):return None
        C,D=endpoints(m,S,I)
        r,rem=divmod(N-C,D)
        if rem:return None
        pair=digits(m,r)
        return (m,I,*pair) if pair is not None else None

def word(m,I,A,B,pad=False):
    assert 0<=A<=m and 0<=B<=m
    J=None
    for i in range(m):
        M=glue(rooted_last(*CORES[i<B]),rooted_last(*ARMS[i<A]))
        if pad:M=glue(M,SINKS)
        J=M if J is None else glue(J,M)
    return glue(J,transitive_arm(I))

def main():
    for m in range(1,121):
        D={17*j+29*k for j in range(m+1) for k in range(m+1)}
        assert all((digits(m,n) is not None)==(n in D) for n in range(-2,46*m+3))
        if m>=28:assert set(range(476,46*m-475))<=D
    m,S,L=next(rows()); I=L
    C,D=endpoints(m,S,I)
    for a,b in product((0,1),repeat=2):
        T=word(m,I,a,b);x=arm(T,len(T)-1)
        assert x==join((S[0]+68*a+116*b,*S[1:]),arm(transitive_arm(I),I.bit_length()))
        assert state(T,len(T)-1)[0]==C+D*(17*a+29*b)
    print('literal axes at shared tape',I,C,C+68,C+116,C+184)
    assert arm(SINKS,2)==(1,1,4,3,-1,0,0)
    C1=1870331984921
    assert decode(C1)==(2,981,1,2)
    assert all(C1%d for d in range(2,__import__('math').isqrt(C1)+1))
    print('prime from a prescribed joint window',C1,decode(C1))
    for m,S,L in rows():
        if m not in (2,3,28,30,31,64):
            if m>64:break
            continue
        I=dvalue(L); C,D=endpoints(m,S,I)
        assert S[3]==U*(B**m-G**m)//(B-G)
        assert S[4]==V*G**(m-1)
        r=476 if m>=28 else 17
        N=C+D*r
        m1,I1,a,b=decode(N)
        assert (m1,I1,17*a+29*b)==(m,I,r)
        assert decode(N+2) is None
        print('selected depth',m,'center address bits',I.bit_length(),
              'first central target digits',len(str(N)), 'repair',(a,b))
        if m==30:
            T=word(m,I,a,b)
            assert all((T[u]>>v&1)+(T[v]>>u&1)==1 for u,v in combinations(range(len(T)),2))
            print('physical constructed (count via transfer) order',len(T))
    m,S,L=next(x for x in rows(PADDED) if x[0]==30)
    h,g,_,u,*_=S
    P=g//G; assert L-dvalue(L)+1>=P
    low=dvalue(L)
    for r in (0,1,17,29,P//2,P-1):
        I=low+((r-h)*pow(u,-1,P)-low)%P
        assert I<=L
        C,D=endpoints(m,S,I)
        N=C+476*D
        assert N%P==r and decode(N,PADDED)==(m,I,28,0)
    print('padded profile',PADDED,'complete phase at depth',m,'modulus decimal digits',len(str(P)))

if __name__=='__main__':main()
