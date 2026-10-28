#!/usr/bin/env python3
from math import isqrt

def need(b,m):
    if not b:raise RuntimeError(m)
def prime(n):return n>=2 and all(n%d for d in range(2,isqrt(n)+1))
def factors(n):
    f={};p=2
    while p*p<=n:
        while n%p==0:f[p]=f.get(p,0)+1;n//=p
        p+=1
    if n>1:f[n]=f.get(n,0)+1
    return f
def stats(f):
    n=sig=tau=1
    for p,a in f.items():n*=p**a;sig*=(p**(a+1)-1)//(p-1);tau*=a+1
    return n,sig-n,tau-1
def main():
    records=[]
    for r in range(2,32):
        if not prime(r):continue
        M=2**r-1;mf=factors(M);f={2:r-1,**mf};n,s,k=stats(f)
        need((s==n)==prime(M),'Mersenne perfect guard')
        if prime(M):need(k==2*r-1,'term count')
        records.append((r,M,prime(M),k))
    need(factors(2047)=={23:1,89:1},'11 hostile')
    need(stats({2:10,23:1,89:1})[2]==43,'failed prime guard changes layer')
    N=100000;s=[0]*(N+1);tau=[0]*(N+1)
    for d in range(1,N+1):
        for n in range(d,N+1,d):tau[n]+=1
        for n in range(2*d,N+1,d):s[n]+=d
    perfect=[n for n in range(2,N+1) if s[n]==n]
    need(perfect==[6,28,496,8128],'finite perfect census')
    need(all(s[n]!=n for n in range(2,N+1) if tau[n]-1 in(7,21)),'excluded layers')
    print('PASS prime-parameter records',records,'finite perfect sources',perfect)
if __name__=='__main__':main()
