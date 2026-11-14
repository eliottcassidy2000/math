#!/usr/bin/env python3
"""Exact constructive congruence decoder in the actual Wave-20 eleven-block orbit."""
from itertools import product
from functools import cache
I=(1,0,0,0,1,0,0,0,1)
S0=(1,1,1)
WORDS=('','O','OO','OZ','OOO','OOZ','OZO','OZZ','OOOZ','OZOZ','OZZZ')
SEEDS=(('OO','OZ'),('OO','OOOZ'),('O','OZOZ'),('OOOZ','OZZ'),('O','OZZ'),('OZZ',),('','OZZZ'))
L=19344
P=(2,3,5,7)
BASE=(4,2,1,1)

def g(q):
    out=1
    while q:
        q,b=divmod(q,2)
        if not b:out+=q+(q&1)
    return out

def mat(w):
    p=int(''.join('1' if c=='O' else '0' for c in w) or '0',2)
    q=int(''.join('1' if c=='Z' else '0' for c in w[::-1]) or '0',2)
    return g(p),q+1,-int(q%2==0),p+1,1<<len(w),0,int(p%2==0),0,-int(not w)
J=mat('')

def mm(a,b,m=0):
    t=tuple(sum(a[3*i+k]*b[3*k+j] for k in range(3)) for i in range(3) for j in range(3))
    return tuple(x%m for x in t) if m else t

def mv(a,b,m=0):
    t=tuple(sum(a[3*i+k]*b[k] for k in range(3)) for i in range(3))
    return tuple(x%m for x in t) if m else t

def pw(a,e,m=0):
    q=I
    while e:
        if e%2:q=mm(a,q,m)
        a=mm(a,a,m); e//=2
    return q

def state(w):
    s=S0
    for v in w:
        assert v in WORDS
        s=mv(mat(v),s)
    return s

def crt(xs,ms):
    x=0;m=1
    for y,d in zip(xs,ms):
        if d==1:continue
        x+=m*(((y-x)*pow(m,-1,d))%d)
        m*=d
    return x

def dot(a,s):return sum(x*y for x,y in zip(a,s))

@cache
def table():
    """All 8*6*7*5 base residues and nonsingular choices."""
    K2=(14,17,-5);K5=(4,1,4);K7=(6,2,5)
    owners={}
    for w in SEEDS:
        s=state(w); v=s; m=16*9*7
        for j in range(624):
            h,_,_=v
            if dot(K2,v)%4==2 and h%3 and dot(K7,v)%7:
                key=(h%16,h%9,h%7)
                if key not in owners:owners[key]=(w,j)
            v=mv(J,v,m)
    assert set(owners)==set(product(range(1,16,2),(1,2,4,5,7,8),range(7)))
    both={}
    for key,(w,j) in owners.items():
        s=state(w);v=s
        for k in range(31):
            if dot(K5,v)%5:
                both.setdefault((*key,v[0]%5),(w,crt((j,k),(624,31))))
            v=mv(J,v,5)
    assert len(both)==8*6*7*5
    return both

def decode(N,a,b,c,d):
    if N%2==0 or N%3==0:raise ValueError('this decoder requires odd N prime to 3')
    power=tuple(max(i,j) for i,j in zip((a,b,c,d),BASE))
    w,j=table()[N%16,N%9,N%7,N%5]
    s=state(w); xs=[]; ms=[]
    for p,k,e in zip(P,power,BASE):
        x=0
        for t in range(k-e):
            m=p**(e+t+1);step=p**t
            options=[x+i*step for i in range(p)
                     if mv(pw(J,j+L*(x+i*step),m),s,m)[0]==N%m]
            assert len(options)==1
            x=options[0]
        xs.append(x);ms.append(p**(k-e))
    return w,j+L*crt(xs,ms)

def run():
    assert [state(w) for w in SEEDS]==[
        (27,35,1),(81,143,1),(53,75,1),(301,333,17),(33,37,1),(11,13,1),(39,41,1)]
    A=pw(J,8)
    K=tuple((x-y)//4 for x,y in zip(A,I))
    assert K==(14,17,-5,17,20,-6,5,6,-2)
    assert all(x%2==0 for x in mm(K,K)[:3])
    for p,l,k in [(3,13,(1,0,0)),(5,31,(4,1,4)),(7,48,(6,2,5))]:
        r=pw(J,l,p*p)
        assert tuple(x%p for x in r)==I
        assert tuple((x-y)//p%p for x,y in zip(r,I))[:3]==k
    import random
    rr=random.Random(5148)
    m=2**13*3**6*5**5*7**4
    for i in range(35):
        n=rr.randrange(1,m)|1
        if n%3==0:n+=2
        w,e=decode(n,13,6,5,4)
        assert mv(pw(J,e,m),state(w),m)[0]==n%m
    print('PASS: 1680 base phases, lift guards, 35 simultaneous 2/3/5/7 targets.')
    for a in (8,12,16,20):
        w,e=decode(67,a,2,1,1)
        print('target 67, binary exponent',a,'seed',w,'literal J letters',e)

if __name__=='__main__':run()
