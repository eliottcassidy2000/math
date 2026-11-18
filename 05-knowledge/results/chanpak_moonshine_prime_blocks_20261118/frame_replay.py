"""Literal checks and exact selectors for same-height switch frames."""
import argparse
from replay import R, dp, pad, prime, profile, transpose, wrap

D6 = (200,29,57,176,65,19,46,118)
D24 = (110,84,152,242,353,454,388,275,15)
TRIPLE=(18,100,89,51,162,197,153,15,255)
COUNT5=(2,4,9,3)
WITNESSES=((D6,101,38344,2393,2399),(D24,65,181928,9599,9623))

def add(D,I):
    return tuple(r|(((I>>i&1)^1)<<len(D)) for i,r in enumerate(D))+(I,)

def cut(A,b):
    f=(1<<len(A))-1
    return tuple(r^(f^b if b>>i&1 else b) for i,r in enumerate(A))

def join(A,B):
    n=len(A);z=((1<<len(B))-1)<<n
    return tuple(r|z for r in A)+tuple(r<<n for r in B)

def weights(D,q):
    E=dp(D);S=dp(transpose(D));f=len(E)-1
    h=[1]+[sum(x) for x in E[1:]]
    if q==1:
        Q=h
    else:
        assert q==2
        Q=[]
        for s in range(f+1):
            z=s;t=0
            while True:
                t+=h[z]*h[s^z]
                if not z:break
                z=(z-1)&s
            Q.append(t)
    a=[[0]*len(D) for _ in range(2)]
    for s in range(1,f+1):
        for x in range(len(D)):
            a[s.bit_count()%2][x]+=(E[s][x]-S[s][x])*Q[f^s]
    return a

def grow(d,positive,j,k,scale):
    U,I,_,_,_=WITNESSES[d]
    A=add(U,I if positive else I^((1<<len(U))-1))
    for _ in range(j):A=join(A,COUNT5)
    D,I,_,_=R[(1,scale)]
    return wrap(D,I,pad(A,k))

def select(N,scale):
    p,r=divmod(N,scale)
    if p<1 or r:return None
    for d,(_,_,b,h0,h1) in enumerate(WITNESSES):
        m=6*b;j=0;u,v=h0,h1
        while m<p:
            for h,positive in [(u,False),(v,True)]:
                a,r=divmod(p-h,m)
                if not r and a>0 and a&(a-1)==0:
                    return grow(d,positive,j,a.bit_length()-1,scale)
            j+=1;m*=32;u*=5;v*=5
    return None

def replay():
    assert profile(COUNT5)==(5,32)
    for i,(U,b,B,h0,h1) in enumerate(WITNESSES):
        T=add(U,b);Tp=cut(T,1<<len(U))
        for q,expect in [(1,(-3,1) if i==0 else (-12,-3)),
                         (2,(0,2) if i==0 else (0,15))]:
            a=weights(U,q)
            got=tuple(sum(z[x] for x in range(len(U)) if b>>x&1) for z in a)
            assert got==expect
        assert profile(T)==(h1,B) and profile(Tp)==(h0,B)
        print('literal elementary swap:',profile(T),'to',profile(Tp))
        for k in (1,2):
            P=join(pad((),k),join(U,pad((),k)))
            h,big=profile(add(P,b<<k))
            hp,bp=profile(add(P,((1<<len(P))-1)^(b<<k)))
            assert hp-h==(-(2 if i==0 else 15)<<k)
            assert bp-big==(2**(k+1)*(3**k-1)*(2 if i==0 else 15))
    assert [profile(cut(TRIPLE,s)) for s in (185,251,164)]==[
        (1127,22156),(1139,22156),(1151,22156)]
    for p in (1115147994227,1115147994239,116443543,116443519,7374013,7374043):
        assert prime(p)
        print('certified prime and partner:',p,7*p)
    assert all(select(p,1) and select(7*p,7) for p in (7374013,7374043))
    for p in (293015554074983,293015554075583):
        assert prime(p)
        assert select(p,1) and select(7*p,7)
        print('certified growing-width prime and partner:',p,7*p)

if __name__=='__main__':
    parse=argparse.ArgumentParser()
    parse.add_argument('N',nargs='?',type=int)
    parse.add_argument('--scale',choices=(1,7),type=int,default=1)
    a=parse.parse_args()
    if a.N is None:replay()
    else:print(select(a.N,a.scale))
