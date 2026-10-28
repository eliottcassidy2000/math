#!/usr/bin/env python3
from math import gcd,isqrt

def need(b,m):
    if not b:raise RuntimeError(m)
def fac(n):
    f={};p=2
    while p*p<=n:
        while n%p==0:f[p]=f.get(p,0)+1;n//=p
        p+=1
    if n>1:f[n]=f.get(n,0)+1
    return f
def sig(n):
    v=1
    for p,a in fac(n).items():v*=(p**(a+1)-1)//(p-1)
    return v
def signature(n):return {p for p,a in fac(n).items() if a%2}
def main():
    F=[n for n in range(1,144) if sum(d for d in range(1,n+1) if n%d==0)==144]
    need(F==[66,70,94,115,119],'complete sigma144 fiber')
    first=[]
    for n in range(2,67):
        f=fac(n);tau=1
        for a in f.values():tau*=a+1
        if len(f)>1 and tau!=4 and isqrt(sig(n))**2==sig(n):first.append(n)
    need(first==[66],'first III square-sigma source')
    for a in range(1,101):
        for b in range(1,101):
            if gcd(a,b)!=1:continue
            need(sig(a*b)==sig(a)*sig(b),'multiplication')
            need(signature(sig(a*b))==signature(sig(a))^signature(sig(b)),'signature')
    N=100000;s=[0]*(N+1);tau=[0]*(N+1)
    for d in range(1,N+1):
        for n in range(d,N+1,d):s[n]+=d;tau[n]+=1
    got=[n for n in range(2,N+1) if tau[n]==8 and s[n]>2*n and isqrt(s[n])**2==s[n]]
    expected=[70];k=1
    while 6*(12*k*k-1)<=N:
        q=12*k*k-1
        if len(fac(q))==1 and fac(q).get(q)==1:expected.append(6*q)
        k+=1
    need(got==sorted(expected),'complete seven-term intersection finite')
    need(sig(22)==36 and sig(66)==sig(70)==144,'triples')
    print('PASS sigma144 fiber',F,'edges',[(n,144-n) for n in F],'seven-term intersection',got)
if __name__=='__main__':main()
