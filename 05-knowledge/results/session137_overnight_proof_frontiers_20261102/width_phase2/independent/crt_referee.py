#!/usr/bin/env python3
"""Independent scalar controller audit; not an actual graph replay."""
from collections import deque
from functools import lru_cache
from math import gcd
from random import Random
from hashlib import sha256
from pathlib import Path

UNIT_MODULI=set()


def require(test,detail='check failed'):
    if not test:raise AssertionError(detail)


def parts(M):
    a=b=0;n=M
    while n%2==0:a+=1;n//=2
    while n%3==0:b+=1;n//=3
    return a,b,n


def factors(n):
    result=[];p=2
    while p*p<=n:
        q=1
        while n%p==0:q*=p;n//=p
        if q>1:result.append((p,q))
        p+=1
    if n>1:result.append((n,n))
    return result


def crt(equations):
    x=0;mod=1
    for y,m in equations:
        require(gcd(mod,m)==1)
        x+=mod*((y-x)*pow(mod,-1,m)%m);mod*=m;x%=mod
    return x


@lru_cache(None)
def unit_words(m):
    UNIT_MODULI.add(m)
    start=1%m;words={start:(0,0,0)};queue=deque([start])
    while queue:
        x=queue.popleft()
        for i,u in enumerate((5,17,23)):
            y=x*u%m
            if y not in words:
                ex=list(words[x]);ex[i]+=1;words[y]=tuple(ex);queue.append(y)
    require(len(words)==sum(gcd(x,m)==1 for x in range(m)),('unit bank',m))
    return words


def unit_order(u,m):
    t=u%m;k=1
    while t!=1%m:t=t*u%m;k+=1
    return k


WORDS={1:['R','RERER','E','ER','RE'],2:['R','ER','RERER','RE','E'],
       3:['R','E','RE','RERER','ER'],4:['R','RE','ER','E','RERER']}


def land(M,z):
    a,b,n=parts(M);two=2**a;three=3**b;r=two*three
    require(M==r*n and gcd(n,6)==1)
    if a:require(z%2==1)
    vv=0
    while vv<b and z%(3**(vv+1))==0:vv+=1
    mod=two*3**(b-vv)
    target=crt([(z*pow(3**vv,-1,two)%two,two),(z//3**vv,3**(b-vv))])
    k,l,e=unit_words(mod)[target]
    if k==0:k+=unit_order(5,mod)
    H=pow(3,vv,M)*pow(5,k,M)*pow(17,l,M)*pow(23,e,M)%M
    require(H%two==z%two and H%three==z%three,('anchors',M,z,H))
    h=H%n;g=pow(12,vv,n)*pow(32,k,n)*pow(224,l,n)*pow(296,e,n)%n
    require(gcd(gcd(h,g),n)==1,('seed gcd',M,z,h,g))
    tau=15*r*pow(4,-1,n)%n
    require(gcd(tau,n)==gcd(5,n))
    require(pow(-27,r,two)==1%two and pow(16,r,three)==1%three and pow(-2,r,three)==1%three)
    c=0;n0=n
    while n0%5==0:c+=1;n0//=5
    if c:require(h%5==0 and g%5!=0)
    equations=[]
    for p,q in factors(n0):
        kk=0 if g%p==0 else (1-h)*pow(tau*g,-1,q)%q
        equations.append((kk,q))
    kk=crt(equations);h=(h+tau*kk*g)%n
    require(gcd(h,n0)==1)
    h,g=((5*h-g)*pow(4,-1,n))%n,h
    require(gcd(g,n0)==1)
    steps=0
    if c:
        require(h%5!=0 and g%5==0)
        word=WORDS[h%5][z%5]
        for op in word:
            before=(h%5,g%5)
            desired=0 if op=='E' else g
            kk=(desired-h)*pow(tau*g,-1,n0)%n0
            h=(h+tau*kk*g)%n
            require((h%5,g%5)==before)
            require(h%n0==desired%n0)
            if op=='E':h,g=(11*h-3*g)*pow(8,-1,n)%n,(3*h-2*g)%n
            else:h,g=(5*h-g)*pow(4,-1,n)%n,h
            require(gcd(g,n0)==1,('safe step',M,z,op,h,g));steps+=1
        require(h%5==z%5 and gcd(g,n)==1)
    else:require(gcd(g,n)==1)
    d=gcd(tau*g,n)
    require((z-h)%d==0)
    kk=((z-h)//d)*pow((tau*g)//d,-1,n//d)%(n//d)
    h=(h+tau*kk*g)%n
    require(h==z%n,('landing',M,z,h))
    # Adding n to the final shear exponent preserves all controlled coordinates.
    require((h+tau*n*g)%n==h)
    return steps


def main():
    count=steps=0
    for M in range(1,257):
        for z in range(M):
            if M%2==0 and z%2==0:continue
            steps+=land(M,z);count+=1
    rng=Random(1371002)
    large=[2**6*3**5*5**3*7*13*37,5**6*13**3*7**2,
           2**7*3**4*5**2*7*13*17*23*37]
    for M in large:
        for _ in range(64):
            z=rng.randrange(M)
            if M%2==0:z|=1
            steps+=land(M,z);count+=1
    # Direct finite exceptional orbit and all 20 promised word-table endpoints.
    seen={(1,0)};todo=[(1,0)]
    while todo:
        h,g=todo.pop()
        for w in [(g,h),((2*h+4*g)%5,(3*h+3*g)%5)]:
            if w not in seen:seen.add(w);todo.append(w)
    require(seen=={(0,1),(0,4),(1,0),(1,2),(2,1),(2,3),(3,2),(3,4),(4,0),(4,3)})
    for t in range(1,5):
        for z,word in enumerate(WORDS[t]):
            h,g=t,0
            for op in word:
                h,g=(g,h) if op=='R' else ((2*h+4*g)%5,(3*h+3*g)%5)
            require(h==z and g!=0)
    require((2+4)%5==1 and (3+3)%5==1) # fixed diagonal hostile
    require((11*2-3*3)%13==0 and (3*2-2*3)%13==0) # singular E hostile
    print('PASS independent CRT scalar landings:',count)
    print('PASS safely preconditioned R/E steps:',steps)
    print('PASS all 20 mod-5 word endpoints, full ten-state orbit, fixed diagonal, singular E at 13')
    print('unit moduli exhaustively generated:',sorted(UNIT_MODULI))
    print('unit_word_cache',unit_words.cache_info())
    print('script_sha256',sha256(Path(__file__).read_bytes()).hexdigest())

if __name__=='__main__':main()
