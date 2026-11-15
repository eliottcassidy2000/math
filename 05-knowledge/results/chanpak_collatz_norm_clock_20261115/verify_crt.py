"""Small exact probes of identities in switch_residual_new.md (not its proof)."""
from math import gcd, lcm
from functools import cache

@cache
def p(x):
    if not x:
        return 0
    u,r=divmod(x,2)
    return p(u) if r else p(u)+u+(u&1)

def app(x,y,w,word,mod=None):
    for e in map(int,word):
        if e:
            w+=y+(y&1)
        else:
            w+=x+(x&1)
        x,y=2*x+e,2*y+1-e
        if mod is not None:
            x,y,w=x%mod,y%mod,w%mod
    return x,y,w

def F(x,y,w,j):
    return 432*x*y+129*x+258*y+15*(x&1)+30*(y&1)+45*w+135+2*j

def order(m):
    v=2%m
    n=1
    while v != 1:
        v=v*2%m;n+=1
    return n

def steers(m):
    L=order(m)
    A='1'; C='01'+'1'*(L-2)
    E='1'*(L-1)+C+'1'+C*(m-2)
    assert len(E)%L==len(C)%L==0
    return L,C,E

for m in (3,5,7,9,15,25,27,45):
    L,C,E=steers(m)
    inv=(m+1)//2
    for x in range(1,10,2):
        for y in range(2,10,2):
            a,b,w=app(x,y,0,E)
            assert ((a-x)%m,(b-y)%m,w%m)==(0,0,(-3*(x-y)*inv-inv)%m)
    for k in range(2,10):
        x,y=3<<k,(5<<k)-1
        for h in range(2,k):
            ell=((h+k+L)//L)*L
            u='0'*ell
            v='1'+'0'*(h-1)+'1'+'0'*(ell-h-1)
            a,b,w=app(x,y,0,u*m)
            A,B,W=app(x,y,0,v*m)
            assert ((A-a)%m,(B-b)%m,(W-w)%m)==(0,0,0)
            assert ((A-a)%(1<<k),(B-b)%(1<<k),(W-w)%(1<<k))==(0,0,(-m*(1<<h))%(1<<k))

for x in range(200):
    for y in range(1,100):
        if (x+y)%2!=1:continue
        w=p(x)+p(y)
        K=w-x+(x&1)
        for j in range(3):
            assert F(x,y,w,j)%4==(K+3+2*j)%4
            if x%2==0:
                assert F(2*x+1,2*y+1,w,j)%4==(K+3+3*x+2*j)%4
print('formal identities: ok')

def steering_check(m,k,x,y,N):
    # Same integer sum for every executable word of final length.
    j=2*N%3
    L,C,E=steers(m)
    prep='0'*k+'1'
    x,y,w=app(x,y,p(x)+p(y),prep)
    K=(w-x+1)%4
    finish='00' if (K+3+2*j-N)%4==0 else '10'
    slots=[]
    for h in range(2,k):
        ell=((h+k+L)//L)*L
        slots.append(('0'*(ell*m), ('1'+'0'*(h-1)+'1'+'0'*(ell-h-1))*m))
    ext='0'*k+''.join(U for U,V in slots)+finish
    B=(m-1)*(2*len(C)+len(E))
    R=(x+y+1) << (B+len(ext)+1)
    if m%5==0 and (3*R-N+2*j)%5:
        return False  # reject the wrong (externally chosen) sum phase
    for X in range(m):
        Xfinal=(pow(2,len(ext)+1,m)*X+2*int(ext,2)+1)%m
        # 45*W adds all multiples d and nothing else modulo m.
        a=432*Xfinal*(R-Xfinal)-129*Xfinal+258*R+180+2*j
        if (N-a)%gcd(45,m)==0:
            break
    else:
        raise RuntimeError('no x residue')
    r=(x-y)%m
    a,b,c=app(x,y,w,C*r)
    assert (a-b)%m==0
    t=2*(a-X)%m

    def full(i,tail):
        cmd=C*r+E*i+C*t
        cmd+='1'*(B-len(cmd))
        a,b,c=app(x,y,w,cmd+tail,m<<k)
        return F(2*a+1,2*b+1,c,j)

    H0=full(0,ext)
    i=next(i for i in range(m) if (H0-i*45*((m+1)//2)-N)%m==0)
    H=full(i,ext)
    assert (H-N)%m==(H-N)%4==0
    wanted=((H-N)//4*pow(45*m,-1,1<<(k-2)))%(1<<(k-2))
    tail='0'*k+''.join(V if wanted>>h&1 else U for h,(U,V) in enumerate(slots))+finish
    H1=full(i,tail)
    assert (H1-N)%(m<<k)==0
    return True

import random
rng=random.Random(10)
executed=0
for m in (3,9,15,21,27,45,63):
    for k in (2,3,5,9):
        for _ in range(10):
            x,y=rng.randrange(10,1000),rng.randrange(10,1000)
            N=rng.randrange(1,10000)|1
            executed+=steering_check(m,k,x,y,N)
print('complete even/odd CRT selector probes:',executed)
