#!/usr/bin/env python3
"""Literal-guarded all-prime constructor; words encoded by signed 1-based macro labels."""
from itertools import product
from math import comb, lcm, prod
from argparse import ArgumentParser

TILES=("OI","OOI","OII")
O=[2,-1]; I=[-1,3]
def inv(w):return [-x for x in reversed(w)]
U=O+inv(I)
Z=[1]+inv(O+I)

def projection(p,t):
    c=[(int(i==0)-comb(p-1,i)*pow(-t,p-1-i,p))%p for i in range(p)]
    w=U*c[0]
    for a in c[1:]:w+=O+U*a
    return w+inv(O)*(p-1)

def positive_word(targets,even=False):
    primes=list(targets)
    E=lcm(*(p*(p-1) for p in primes))
    radical=prod(primes)
    word=O[:] if even else []
    for p,target in targets.items():
        d=(p-1)//2
        assert len(target)==d
        scale=2 if even else 1
        # Request delta in raw channel units, before multiplying by s^-2.
        delta=[(r*pow(s,-scale,p)-1)*s*s%p for s,r in enumerate(target,1)]
        isolate=radical//p*pow(radical//p,-1,p)%radical
        for j in range(d-1,0,-1):
            t=(1+pow(j,-1,p))%p;tinv=pow(t,-1,p)
            e=delta[j]*pow(-t,-1,p)%p
            if e:
                a,b=projection(p,t),projection(p,tinv)
                word+=(a+b+inv(a)+inv(b))*(e*isolate%radical)
                delta[j-1]=(delta[j-1]-e*tinv)%p
        word+=Z*(delta[0]*isolate%radical)
    if not word:return [(1,E)]
    # Nonnegative repetition list. Negative macro letters replaced by positive powers.
    out=[]
    for t in word:
        macro=abs(t); count=1 if t>0 else E-1
        if out and out[-1][0]==macro:out[-1]=(macro,out[-1][1]+count)
        else:out.append((macro,count))
    return out

def mul(u,v,p):
    a,b,d,g=u; A,B,D,G=v
    return ((A+a*D)%p,(b+d*B)%p,d*D%p,(g+G+a*B)%p)

def read(word,p):
    n=sum(len(TILES[t-1])*k for t,k in word)
    def channel(lam):
        if not lam:return sum(k for _,k in word)%p
        o=(1,0,lam,0);i=(0,1,lam,0)
        x=one=(0,0,1,0)
        for t,k in word:
            y=one
            for c in TILES[t-1]:y=mul(y,o if c=="O" else i,p)
            while k:
                if k&1:x=mul(x,y,p)
                y=mul(y,y,p);k>>=1
        return x[-1]
    result=[]
    for s in range(1,(p+1)//2):
        z=pow(s,-1,p)
        result.append(pow(s,n+1,p)*(1+z*z*(channel((1+z)%p)+channel((1-z)%p)))%p)
    return n+1,result

def literal(word,smax):
    w="".join(TILES[t-1]*k for t,k in word)
    n=len(w);f=(1<<n)-1
    out=[f^((1<<(i+1))-1) for i in range(n)]+[0]
    for j,t in enumerate(w):
        if t=="O":out[n]|=1<<j
        else:out[j]|=1<<n
    N=1<<(n+1);end=[[0]*(n+1) for _ in range(N)]
    for j in range(n+1):end[1<<j][j]=1
    for a in range(1,N):
        for j,v in enumerate(end[a]):
            if v:
                b=out[j]&~a
                while b:
                    c=b&-b;b-=c;end[a|c][c.bit_length()-1]+=v
    H=[1]+list(map(sum,end[1:]))
    Q=[1]+[0]*(N-1)
    answer=[]
    for _ in range(smax):
        nxt=[]
        for a in range(N):
            b=a;v=Q[a]
            while b:
                v+=H[b]*Q[a^b];b=(b-1)&a
            nxt.append(v)
        Q=nxt;answer.append(Q[-1])
    return answer

def check():
    count=0
    for size in range(1,4):
        for ts in product((1,2,3),repeat=size):
            word=[(t,1) for t in ts]
            L=literal(word,5)
            for p in (3,5,7,11):
                assert read(word,p)[1]==[x%p for x in L[:(p-1)//2]]
            count+=1
    print(count,'literal small path-cover tournament checks')
    for p in (3,5,7,11,13,17,19,23,31,41,53):
        d=(p-1)//2
        for even in (False,True):
            for target in ([0]*d,[int(i==d-1) for i in range(d)],[(i*i+2*i+3)%p for i in range(d)]):
                word=positive_word({p:target},even)
                n,actual=read(word,p)
                assert (n%2==0)==even and actual==target
    print('66 single-prime constructor and parity checks, including p=23')
    targets={p:[(i*i+4*p+i+2)%p for i in range((p-1)//2)] for p in (3,5,7,11)}
    for e in (False,True):
        w=positive_word(targets,e)
        for p,target in targets.items():
            n,actual=read(w,p)
            assert (n%2==0)==e and actual==target
    print('joint physical constructor passed at 3,5,7,11 for both parities')

if __name__=='__main__':
    ap=ArgumentParser(description=__doc__)
    ap.add_argument('--prime',type=int)
    ap.add_argument('--even',action='store_true')
    ap.add_argument('target',type=int,nargs='*')
    args=ap.parse_args()
    if args.prime:
        w=positive_word({args.prime:args.target},args.even)
        for t,k in w:print(TILES[t-1],k)
    else:check()
