"""Independent quotient-word replay for the marked-core pad-score collision."""
from fractions import Fraction
from math import comb, gcd
from functools import cache
from padding_replay import slots, rooted
from frame_replay import join, add, cut
from replay import pad, profile, reads, wrap, prime

L=(226,236,73,81,103,76,128,60)
R=(174,68,136,210,167,78,149,34)
A=pad(L,2)
B=pad(R,2)
UPPER=((20,17,10,19,228,15,47,111),192)
LOWER=((16,1,3,7,78,31,47),96)
PADS=(A,B)
U=pad((),3)
C=5

def nonempty(f,m):
    return sum((-1)**(m-j)*comb(m,j)*f[j] for j in range(m+1))

def shells(q,flip=False):
    # Modules 0 and 4, ordinary exterior 1,2,3, contact vertex 5.
    M=add(pad((),5),C<<1)
    if flip:M=cut(M,1<<5)

    @cache
    def rec(counts,k,prev):
        if not any(counts):return 1
        s=rec(counts,k-1,-1) if k>1 else 0
        for x,z in enumerate(counts):
            if z and (prev<0 or M[prev]>>x&1):
                b=list(counts);b[x]-=1
                s+=rec(tuple(b),k,x)
        return s

    # One path sees either module at most twice, with the possible
    # extra visit using v. There is only one v across the q slots.
    result={(a,b):rec((a,1,1,1,b,1),q,-1)
            for a in range(1,q+2) for b in range(1,q+2)}
    return {ab:t for ab,t in result.items() if t}

def build(flip,scale):
    D=join(A,join(U,B));T=add(D,C<<len(A))
    if flip:T=cut(T,1<<len(D))
    K,I=LOWER if scale==1 else UPPER
    return wrap(K,I,pad(T,3))

def run():
    f,g=([v[-1] for v in slots(P,3)] for P in PADS)
    assert f==[0,359,22208,393903] and g==[0,355,22336,396171]
    assert rooted(U,C)[1:3]==[(5,5),(13,13)]
    score=lambda x: Fraction(13*x[3]-5*x[1],x[2])
    assert score(f)==score(g)
    u,v=Fraction(f[2],f[1]),Fraction(g[2],g[1])
    assert u.numerator==64*347 and v.numerator==64*349
    assert prime(347) and prime(349)
    print('common score',score(f),'new normalized height primes',347,349)
    profiles=[]
    for flip in (False,True):
        y=[]
        for q in (1,2):
            weights=shells(q,flip)
            y.append(sum(t*nonempty(f,a)*nonempty(g,b) for (a,b),t in weights.items()))
            print('templates',int(flip),q,sorted(weights.items()))
        profiles.append(tuple(y))
    assert profiles==[(39965675,232641769472),(39291755,232641769472)]
    for scale,(K,I) in [(1,LOWER),(7,UPPER)]:
        h,c,e=reads(K,I)
        assert (c,h+e-c)==(62*scale,scale)

    p,q=[62*y*2**3+x for x,y in profiles]
    # Deterministic Pocklington certificates, each with F>sqrt(target).
    for z,Fparts in [(p,(437501,6279733)),(q,(111441107,))]:
        F=1;s=[]
        for a in Fparts:
            assert prime(a)
            F*=a
            w=next(w for w in range(2,40)
                   if pow(w,z-1,z)==1 and gcd(pow(w,(z-1)//a,z)-1,z)==1)
            s.append((a,w))
        assert (z-1)%F==0 and F*F>z
        print('certified prime and partner',z,7*z,'Pocklington',s)

if __name__=='__main__':run()
