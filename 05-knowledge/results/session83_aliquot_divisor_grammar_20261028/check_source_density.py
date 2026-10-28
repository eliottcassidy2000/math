#!/usr/bin/env python3
from math import prod,isqrt
from collections import Counter

def need(b,m):
    if not b:raise RuntimeError(m)
def main():
    small=[2,3,5,7,11];M=prod(small);coeff=[1]
    for p in small:
        new=[0]*(len(coeff)+1)
        for j,v in enumerate(coeff):new[j]+=(p-1)*v;new[j+1]+=v
        coeff=new
    brute=Counter(sum(n%p==0 for p in small) for n in range(M))
    need(coeff==[brute[j] for j in range(6)],'CRT polynomial')
    N=100000;spf=list(range(N+1))
    for p in range(2,isqrt(N)+1):
        if spf[p]==p:
            for n in range(p*p,N+1,p):
                if spf[n]==n:spf[n]=p
    ps=[p for p in range(2,N+1) if spf[p]==p];I=set();II=set()
    for p in ps:
        n=p
        while n<=N:I.add(n);n*=p
    for i,p in enumerate(ps):
        if p*p>N:break
        for q in ps[i+1:]:
            if p*q>N:break
            II.add(p*q)
    need(not I&II,'partition overlap');counts=Counter()
    for n in range(2,N+1):
        x=n;omega=0;tau=1
        while x>1:
            p=spf[x];a=0
            while x%p==0:x//=p;a+=1
            omega+=1;tau*=a+1
        c=1 if omega==1 else 2 if tau==4 else 3
        expected=1 if n in I else 2 if n in II else 3
        need(c==expected,'independent classification');counts[c]+=1
    need(all(n not in I|II for n in [12,18,20,28,36]),'two-factor III')
    print('PASS CRT counts',coeff,'source cutoff',N,'class counts',dict(counts))
if __name__=='__main__':main()
