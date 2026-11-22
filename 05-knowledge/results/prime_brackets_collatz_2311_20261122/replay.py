from collections import Counter
from math import gcd, isqrt, prod

def v(x,p):
    assert x
    a=0
    while x%p==0:
        a+=1
        x//=p
    return a

def step(x):
    return (3*x-1)//2 if x&1 else x//2

def cycle(s):
    x=s; states=[]; A=[]
    while True:
        states.append(x)
        a=v(3*x-1,2)
        A.append(a)
        x=(3*x-1)>>a
        if x==s: return tuple(A),states
        assert x not in states

def check_orders(x,k,A,e,u):
    a,d=3**k,2**A
    assert v(a**e-d**e,11)==u
    for j in range(1,7):
        m=11**j
        z=a*pow(d,-1,m)%m
        r=e*11**max(0,j-u)
        assert pow(z,r,m)==1
        assert all(pow(z,r//p,m)!=1 for p in (2,5,11) if not r%p)
    for j in (1,2):
        r=e*11**max(0,j-u)
        for s in (1,3):
            d=1<<(A*r+s)
            x0=x+d; y=x; z=x0
            for _ in range(A*r):
                assert y%2==z%2
                y,z=step(y),step(z)
            assert y==x and z!=x0 and (z-x0)%(2**s*11**j)==0
            assert (z-x0)%(2**(s+1)*11**j)!=0
    print(x,"word",cycle(x)[0],"mod 11 order",e,"lift depth",u)

def spectrum(ps,a):
    remaining=set(__import__("itertools").product(*(range(p+1) for p in ps)))
    found=Counter()
    while remaining:
        x=next(iter(remaining))
        d=0
        while x in remaining:
            remaining.remove(x)
            x=tuple(a*u%p if u!=p else p for p,u in zip(ps,x))
            d+=1
        found[d]+=1
    return dict(sorted(found.items()))

if __name__=="__main__":
    for x,word,k,A,e,u in ((1,(1,),1,1,10,1),(5,(1,2),2,3,10,1),
                           (17,(1,1,1,2,1,1,4),7,11,2,2)):
        assert cycle(x)[0]==word and len(word)==k and sum(word)==A
        check_orders(x,k,A,e,u)
    assert (pow(3,70,11**3),pow(2,110,11**3))==(727,606)
    for p,e,z in ((5,4,1),(7,6,1),(11,5,2)):
        assert all(pow(3,d,p)!=1 for d in range(1,e))
        assert v(3**e-1,p)==z
    assert spectrum((3,11),2)=={1:4,2:2,10:4}
    assert spectrum((5,7),2)=={1:4,3:4,4:2,12:2}
    assert spectrum((2,3,11),35)=={1:12,2:6,10:12}
    assert spectrum((2,5,7),37)=={1:12,3:12,4:6,12:6}
    print("prime powers, prime-specific return, and projective spectra PASS")

    for n,t in ((4,712944),(8,57602544)):
        S=2**n+1
        assert 117*t-85680*S+986112==1024000*3**n
    assert 963*pow(1000,-1,121)%121==49
    assert pow(3,5,121)==1 and pow(3,3,121)==27
    for depth in range(1,5):
        M=11**depth
        D=11**(depth-1)
        U=set()
        for n0 in (4,8,12,16,20):
            items={}
            for j in range(D):
                n=n0+20*j
                S=(pow(2,n,M)+1)%M
                V=(1024000*pow(3,n,11*M))%(11*M)
                assert S not in items
                items[S]=V
                assert (1024000*pow(3,n+20*D,11*M)-V)%(11*M)==0
            assert len(items)==D and U.isdisjoint(items)
            U.update(items)
            assert all(S%11==(1+pow(2,n0,11))%11 for S in items)
            if n0==20:
                assert {(S-2)*(2*S-1)%M for S in items}==set(range(0,M,11))
        assert len(U)==5*D
    print("one-digit-delayed Farey coefficient and moving Euclid target controls PASS")
