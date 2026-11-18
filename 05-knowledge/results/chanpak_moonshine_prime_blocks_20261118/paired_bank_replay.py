"""Construct the lower witness and its genuine seven-rotation partner."""
import csv
from math import isqrt,gcd
from pathlib import Path
import sys
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'moonshine_prime_blocks'))
sys.path.insert(0,str(HERE.parent/'pole_prime_placement'))
from clone_replay import P,inflate,vectors
from dihedral_replay import C,count,polynomial
from frame_replay import join
from replay import pad,wrap,reads,profile
from crt_replay import phase,toggle,C11,TT

def bank():
    with (HERE/'rotational_bank.tsv').open() as f:
        for d in csv.reader(f,delimiter='\t'):
            if d[0]=='alpha':continue
            a,b,g,h,i,c,u=map(int,d[:7])
            D,Q=(tuple(map(int,w.split(','))) for w in d[7:])
            yield (a,b,g,h,i,c,u,D,Q)

def build(row,X):
    a,b,g,h,i,c,u,D,Q=row
    m=lambda a,b:(a//b).bit_length()-1
    return wrap(D,i,pad(X,m(a,g))),inflate(Q,{c:C,u:pad(X,m(a,h))})

def decode(p):
    """Complete target-first inverse for this bank over P_m and TT padding."""
    chosen=[]
    for row in bank():
        a,b,*_=row
        m=1
        while a*3**m<=p:
            h,B=1+2**m,4*2**m+2*3**m-2
            z,r=divmod(p-b*h,a*B)
            if not r and z>0 and not z&(z-1):
                chosen.append((row,m,z.bit_length()-1))
            m+=1
    return chosen

def prime_phase(t,L,i,j):
    """Physical same-source phase mod 2^L * 3^i * 5^j."""
    M=3**i*5**j<<L
    assert j>=1 and L>=1 and i>=1 and gcd(t,30)==1
    row=next(R for R in bank() if R[:2]==(5400,-877))
    r=(max(0,j-2)+1)//2
    u=t*pow((-877*11**r)%M,-1,M)%M
    T,h,b=phase(u,L)
    def add(a,k,c):
        nonlocal T,h,b
        U,s,z=toggle(a,k)
        for _ in range(c):T=join(T,U)
        h*=s**c;b*=z**c
    mod=5**j
    c=L+(5-L)%20
    gen=(1+3*2**c)%mod
    for n in range(4*5**(j-1)):
        if h*pow(gen,n,mod)%mod==u%mod:
            add(3,c,n);break
    else:assert False
    for _ in range(r):
        T=join(T,C11);h*=11;b*=100
    T=pad(T,L);b<<=L
    a=L+(1-L)%4
    if a<2:a+=4
    _,H,B=toggle(3,a)
    e=4*5**(j-1)
    genH,genB=pow(H,e,3**i),pow(B,e,3**i)
    def value(n,mod):return (5400*b*pow(genB,n,mod)-877*h*pow(genH,n,mod))%mod
    n=0
    assert value(n,3)==t%3
    for r in range(i-1):
        step=3**r;mod=3**(r+2)
        for d in range(3):
            if value(n+d*step,mod)==t%mod:
                n+=d*step;break
        else:assert False
    add(3,a,e*n)
    p=5400*b-877*h
    assert p%M==t%M
    return row,T,p

def run():
    rows=list(bank())
    assert len(rows)==26
    for row in rows:
        a,b,g,h,i,c,u,D,Q=row
        d,G,z=reads(D,i)
        assert (G,d+z-G)==(g,b)
        S=polynomial(Q,(c,u))
        assert S=={(1,):2*h+b,(2,):h}
        assert all(a//k==1<<((a//k).bit_length()-1) and a%k==0 for k in (g,h))
        T,W=build(row,(0,))
        assert count(T)==2*a+b and count(W)==7*(2*a+b)
    print('all 26 coefficient and literal singleton pairs verified')
    for p in (4441949209,26651515879,125864273,195414419,12152136127):
        assert all(p%i for i in range(2,isqrt(p)+1))
        sols=decode(p)
        assert sols
        R,m,k=sols[0]
        print('absolute target',p,7*p,'row',R[:2],'m,k',m,k)
    for L,i,j,t in ((3,1,1,17),(4,3,2,29),(2,4,2,221)):
        R,T,p=prime_phase(t,L,i,j)
        print('phase',L,i,j,t,'input order',len(T),'restored',p%(3**i*5**j<<L))
if __name__=='__main__':run()
