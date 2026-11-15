#!/usr/bin/env python3
"""Search actual (H,Q) pairs against all marked hosts through order six."""
import sys
from repair import coeffs

def run(profilefile,hostfile):
    cap=59999
    band=set(range(69,cap+1,2))
    central={n for n in band if n>32768}
    pq=set()
    direct=set()
    for line in open(profilefile):
        n,h,p,c,s=line.split()
        h,p=int(h),int(p)
        q=p+2*h
        if q<=cap:pq.add((q,h))
        if h in band:direct.add(h)
    pq=sorted(pq)
    readers={}
    for n,c in sorted(tuple(map(int,line.split())) for line in open(hostfile)):
        if n>6:continue
        for (r,M),mask in coeffs(n,c):
            assert M>0 and r>=1-M
            readers.setdefault((r,M),n)
    hit=set()
    for (r,M),size in sorted(readers.items(),key=lambda x:x[1]):
        for q,h in pq:
            if M*q>2*cap:break
            N=r*h+M*q
            if N in band:hit.add(N)
    miss=sorted(band-hit)
    print('readers',len(readers),'source pairs',len(pq))
    print('central bank coverage',len(hit&central),'/',len(central))
    print('bank coverage',len(hit),'/',len(band))
    print('without direct H:',miss[:200],len(miss))
    print('after direct H:',sorted(band-hit-direct)[:200],len(band-hit-direct))
    from math import isqrt
    missing_primes=[n for n in miss if all(n%d for d in range(3,isqrt(n)+1,2))]
    print('missing primes',missing_primes[:200],len(missing_primes))
if __name__=='__main__':run(*sys.argv[1:])
